import pandas as pd
import tensorflow as tf
from itertools import chain
import logging

from ms2rescore.feature_generators._base_classes import FeatureGeneratorBase
from psm_utils import PSMList
from ionmob.preprocess.data import to_tf_dataset_inference
from ionmob.utilities import tokenizer_from_json, reduced_mobility_to_ccs, get_ccs_shift
from ionmob.preprocess.data import calculate_mz
from ionmob.utilities.chemistry import VARIANT_DICT


logger = logging.getLogger(__name__)


class IonMobFeatureGenerator(FeatureGeneratorBase):
    """Ionmob Collision Cross Section (CCS)-based feature generator."""

    def __init__(
        self,
        *args,
        model: str = "pretrained_models/GRUPredictor",
        reference_dataset="example_data/reference.parquet",
        tokenizer_filepath="pretrained_models/tokenizers/tokenizer.json",
        processes: 1,
        **kwargs,
    ) -> None:
        """
        Ionmob Collision Cross Section (CCS)-based feature generator.

        Parameters
        ----------
        #TODO

        """
        super().__init__(*args, **kwargs)
        self.model = tf.keras.models.load_model(model)
        self.processes = processes
        self.reference_dataset = pd.read_parquet(reference_dataset)
        self.tokenizer = tokenizer_from_json(tokenizer_filepath)

    @property
    def feature_names(self):
        return [
            "ccs_predicted",
            "ccs_observed",
            "ccs_error",
            "abs_ccs_error",
            "perc_ccs_error",
        ]

    @property
    def allowed_mods(self):
        return [token for aa_tokens in VARIANT_DICT.values() for token in aa_tokens]

    def add_features(self, psm_list: PSMList) -> None:
        """
        Add Ionmob-derived features to PSMs.

        Parameters
        ----------
        psm_list
            PSMs to add features to.

        """
        logger.info("Adding Ionmob-derived features to PSMs.")
        psm_dict = psm_list.get_psm_dict()
        current_run = 0
        total_runs = len(psm_dict.values())

        for runs in psm_dict.values():
            for run, psms in runs.items():
                logger.info(
                    f"Running Ionmob for PSMs from run ({current_run}/{total_runs}): `{run}`..."
                )
                psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
                psm_list_run_df = psm_list_run.to_dataframe()

                # prepare dataframes for CCS prediction
                psm_list_run_df["charge"] = [
                    peptidoform.charge for peptidoform in psm_list_run_df["peptidoform"]
                ]
                psm_list_run_df = psm_list_run_df[
                    psm_list_run_df["charge"] < 5
                ]  # predictions do not go higher for ionmob

                psm_list_run_df["sequence-tokenized"] = psm_list_run_df["peptidoform"].apply(
                    self.proforma_tokenizer, axis=1
                )
                psm_list_run_df = psm_list_run_df[~(self.invalid_mods(psm_list_run_df))]
                psm_list_run_df["mz"] = psm_list_run_df.apply(
                    lambda x: calculate_mz(x["sequence-tokenized"], x["charge"]), axis=1
                )  # use precursor m/z from PSMs?

                psm_list_run_df["ccs_observed"] = psm_list_run_df.apply(
                    lambda x: reduced_mobility_to_ccs(x["ion_mobility"], x["mz"], x["charge"]),
                    axis=1,
                )

                # calibrate CCS values
                shift_factor = self.calculate_ccs_shift(self, psm_list_run_df)
                psm_list_run_df["ccs_observed"] = psm_list_run_df.apply(
                    lambda r: r["ccs_observed"] + shift_factor, axis=1
                )

                # predict CCS values
                tf_ds = to_tf_dataset_inference(
                    psm_list_run_df["mz"],
                    psm_list_run_df["charge"],
                    psm_list_run_df["sequence-tokenized"],
                    self.tokenizer,
                )

                psm_list_run_df["ccs_predicted"], _ = self.model.predict(tf_ds)

                # calculate CCS features
                ccs_features = self._calculate_features(psm_list_run_df)

                # add CCS features to PSMs
                for psm in psms.values():
                    try:
                        psm["rescoring_features"].update(ccs_features[psm.spectrum_id])
                    except KeyError:
                        psm["rescoring_features"].update({})
                current_run += 1

    def _calculate_features(self, feature_df):
        """Get ccs features for PSMs."""

        ccs_features = {}
        for row in feature_df.iterrows():
            ccs_features[row["spectrum_id"]] = {
                "ccs_predicted": row["ccs_predicted"],
                "ccs_observed": row["ccs_observed"],
                "ccs_error": row["ccs_observed"] - row["ccs_predicted"],
                "abs_ccs_error": abs(row["ccs_observed"] - row["ccs_predicted"]),
                "perc_ccs_error": (row["abs_ccs_error"] / row["ccs_observed"]) * 100,
            }

        return ccs_features

    @staticmethod
    def proforma_tokenizer(peptidoform):
        """
        Tokenize proforma sequence and add modifications.

        Args:
            seq (str): Peptide sequence.
            peprec_mod (str): Peptide modifications in the format "loc1|mod1|loc2|mod2|...".

        Returns:
            list: A list of tokenized and modified peptide sequence.
        """
        tokenized_seq = []

        if peptidoform.properties["n_term"]:
            tokenized_seq.append(
                f"<START>[UNIMOD:{peptidoform.properties['n_term'].definition['id']}]"
            )
        else:
            tokenized_seq.append("<START>")

        if peptidoform.properties["c_term"]:
            pass  # provide if c-term mods are supported

        for amino_acid, modification in peptidoform.parsed_sequence:
            tokenized_seq.append(amino_acid)
            if modification:
                tokenized_seq[-1] = tokenized_seq[-1] + tokenized_seq.append(
                    f"[UNIMOD:{modification[0].definition['id']}]"
                )

        return tokenized_seq

    def calculate_ccs_shift(self, psm_dataframe):
        """
        Apply CCS shift to CCS values.

        Args:
            peprec (pandas.DataFrame): Modified and parsed Peprec data.
            reference (str): Path to the reference data used for CCS shift calculation.

        Returns:
            pandas.DataFrame: Peprec data with CCS values after applying the shift.
        """
        df = psm_dataframe.copy()
        df["charge"] = [peptidoform.charge for peptidoform in df["peptidoform"]]
        high_conf_hits = list(
            psm_dataframe["spectrum_id"][psm_dataframe["score"].rank(pct=True) > 0.95]
        )
        logger.debug(
            f"Number of high confidence hits for calculating shift: {len(high_conf_hits)}"
        )

        shift_factor = get_ccs_shift(
            psm_dataframe[psm_dataframe["spec_id"].isin(high_conf_hits)][
                ["charge", "sequence-tokenized", "ccs"]
            ],
            self.reference_data,
        )

        logger.debug(f"CCS shift factor: {shift_factor}")

        return shift_factor

    def invalid_mods(self, tokenized_seq):
        """
        Check if peptide sequence contains invalid modifications.

        Args:
            tokenized_seq (list): Tokenized peptide sequence.

        Returns:
            bool: True if invalid modifications are present, False otherwise.
        """
        for token in tokenized_seq:
            if token not in self.allowed_mods:
                logger.debug(f"Invalid modification found: {token}")
                return True
        return False


def add_ccs_predictions(peprec, model_path, tokenizer_filepath):
    """
    Add CCS predictions to peprec file.

    Args:
        peprec (pandas.DataFrame): Modified and parsed Peprec data.
        model_path (str): Path to the CCS prediction model.
        tokenizer_filepath (str): Path to the tokenizer file.

    Returns:
        pandas.DataFrame: Peprec data with added CCS predictions and error values.
    """
    logger.info(f"Adding CCS predictions to peprec file")

    tokenizer = tokenizer_from_json(tokenizer_filepath)
    tf_ds = to_tf_dataset_inference(
        peprec["mz"], peprec["charge"], peprec["sequence-tokenized"], tokenizer
    )
    gruModel = tf.keras.models.load_model(model_path)
    peprec["ccs_predicted"], _ = gruModel.predict(tf_ds)

    peprec["ccs_error"] = peprec["ccs_observed"] - peprec["ccs_predicted"]
    peprec["abs_ccs_error"] = abs(peprec["ccs_observed"] - peprec["ccs_predicted"])
    peprec["perc_ccs_error"] = (peprec["abs_ccs_error"] / peprec["ccs_observed"]) * 100

    peprec["sequence-tokenized"] = peprec.apply(lambda x: "".join(x["sequence-tokenized"]), axis=1)

    return peprec


def write_pin_files(peprec, pin, pin_filepath):
    """
    Write pin files.

    Args:
        peprec (pandas.DataFrame): Peprec data.
        pin (pandas.DataFrame): Percolator In data.
        pin_filepath (str): Path to the Percolator In file.

    Returns:
        None
    """
    ccs_filename, non_ccs_filename = create_filenames(pin_filepath)
    ccs_features = [
        "ccs_predicted",
        "ccs_observed",
        "ccs_error",
        "abs_ccs_error",
        "perc_ccs_error",
    ]
    final_pin = pd.merge(
        pin, peprec[["spec_id"] + ccs_features], left_on="SpecId", right_on="spec_id"
    ).drop("spec_id", axis=1)

    final_pin = final_pin[
        [c for c in final_pin.columns if c not in ["Peptide", "Proteins"]]
        + ["Peptide", "Proteins"]
    ]

    logger.info(f"Writing pin files to {ccs_filename} and {non_ccs_filename}")
    final_pin.to_csv(ccs_filename, sep="\t", index=False, header=True)
    redo_pin_tabs(str(ccs_filename))
    final_pin.drop(ccs_features, axis=1).to_csv(
        non_ccs_filename, sep="\t", index=False, header=True
    )
    redo_pin_tabs(str(non_ccs_filename))
