import pandas as pd
import tensorflow as tf
from itertools import chain
import logging
from pathlib import Path
import contextlib
import os

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from psm_utils import PSMList
from ionmob.preprocess.data import to_tf_dataset_inference
from ionmob.utilities.utility import get_ccs_shift
from ionmob.utilities.tokenization import tokenizer_from_json
from ionmob.utilities.chemistry import VARIANT_DICT, reduced_mobility_to_ccs, calculate_mz
import ionmob


logger = logging.getLogger(__name__)

ionmob_dir = Path(ionmob.__file__).parent

DEFAULT_MODELS_IONMOB = {
    Path("pretrained_models/DeepTwoMerModel"),
    Path("pretrained_models/GRUPredictor"),
    Path("pretrained_models/SqrtModel"),
}
DEFAULT_MODELS_DICT = {
    str(mod_path.stem): ionmob_dir.joinpath(mod_path) for mod_path in DEFAULT_MODELS_IONMOB
}
DEFAULT_TOKENIZER = ionmob_dir.joinpath("pretrained_models/tokenizers/tokenizer.json")
DEFAULT_REFERENCE_DATASET = ionmob_dir.joinpath("example_data/Tenzer_unimod.parquet")


class IonMobFeatureGenerator(FeatureGeneratorBase):
    """Ionmob Collision Cross Section (CCS)-based feature generator."""

    def __init__(
        self,
        *args,
        ionmob_model: str = "GRUPredictor",
        reference_dataset: str = DEFAULT_REFERENCE_DATASET,
        tokenizer: str = DEFAULT_TOKENIZER,
        **kwargs,
    ) -> None:
        """
        Ionmob Collision Cross Section (CCS)-based feature generator.

        Parameters
        ----------
        #TODO

        """
        super().__init__(*args, **kwargs)
        try:
            self.ionmob_model = tf.keras.models.load_model(
                DEFAULT_MODELS_DICT[ionmob_model].__str__()
            )
        except KeyError:
            self.ionmob_model = tf.keras.models.load_model(Path(ionmob_model).absolute().__str__())

        self.reference_dataset = pd.read_parquet(reference_dataset)
        self.tokenizer = tokenizer_from_json(tokenizer)

        self._verbose = logger.getEffectiveLevel() <= logging.DEBUG

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
        current_run = 1
        total_runs = len(list(chain.from_iterable([runs.keys() for runs in psm_dict.values()])))

        for runs in psm_dict.values():
            for run, psms in runs.items():
                logger.info(
                    f"Running Ionmob for PSMs from run ({current_run}/{total_runs}): `{run}`..."
                )
                with contextlib.redirect_stdout(
                    open(os.devnull, "w")
                ) if not self._verbose else contextlib.nullcontext():
                    psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
                    psm_list_run_df = psm_list_run.to_dataframe()

                    # prepare dataframes for CCS prediction
                    psm_list_run_df["charge"] = [
                        peptidoform.precursor_charge
                        for peptidoform in psm_list_run_df["peptidoform"]
                    ]
                    psm_list_run_df = psm_list_run_df[
                        psm_list_run_df["charge"] < 5
                    ]  # predictions do not go higher for ionmob

                    psm_list_run_df["sequence-tokenized"] = psm_list_run_df.apply(
                        lambda x: self.tokenize_peptidoform(x["peptidoform"]), axis=1
                    )
                    psm_list_run_df = psm_list_run_df[
                        psm_list_run_df.apply(
                            lambda x: self._is_valid_tokenized_sequence(x["sequence-tokenized"]),
                            axis=1,
                        )
                    ]

                    psm_list_run_df["mz"] = psm_list_run_df.apply(
                        lambda x: calculate_mz(x["sequence-tokenized"], x["charge"]), axis=1
                    )  # use precursor m/z from PSMs?

                    psm_list_run_df["ccs_observed"] = psm_list_run_df.apply(
                        lambda x: reduced_mobility_to_ccs(x["ion_mobility"], x["mz"], x["charge"]),
                        axis=1,
                    )
                    # calibrate CCS values
                    shift_factor = self.calculate_ccs_shift(psm_list_run_df)
                    psm_list_run_df["ccs_observed"] = psm_list_run_df.apply(
                        lambda x: x["ccs_observed"] + shift_factor, axis=1
                    )
                    # predict CCS values
                    tf_ds = to_tf_dataset_inference(
                        psm_list_run_df["mz"],
                        psm_list_run_df["charge"],
                        psm_list_run_df["sequence-tokenized"],
                        self.tokenizer,
                    )

                    psm_list_run_df["ccs_predicted"], _ = self.ionmob_model.predict(tf_ds)

                    # calculate CCS features
                    ccs_features = self._calculate_features(psm_list_run_df)

                    # add CCS features to PSMs
                    for psm in psm_list_run:
                        try:
                            psm["rescoring_features"].update(ccs_features[psm.spectrum_id])
                        except KeyError:
                            psm["rescoring_features"].update({})
                    current_run += 1

    def _calculate_features(self, feature_df):
        """Get ccs features for PSMs."""

        ccs_features = {}
        for row in feature_df.itertuples():
            ccs_features[row.spectrum_id] = {
                "ccs_predicted": row.ccs_predicted,
                "ccs_observed": row.ccs_observed,
                "ccs_error": row.ccs_observed - row.ccs_predicted,
                "abs_ccs_error": abs(row.ccs_observed - row.ccs_predicted),
                "perc_ccs_error": ((abs(row.ccs_observed - row.ccs_predicted)) / row.ccs_observed)
                * 100,
            }

        return ccs_features

    @staticmethod
    def tokenize_peptidoform(peptidoform):
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
                f"<START>[UNIMOD:{peptidoform.properties['n_term'][0].definition['id']}]"
            )
        else:
            tokenized_seq.append("<START>")

        for amino_acid, modification in peptidoform.parsed_sequence:
            tokenized_seq.append(amino_acid)
            if modification:
                tokenized_seq[-1] = (
                    tokenized_seq[-1] + f"[UNIMOD:{modification[0].definition['id']}]"
                )

        if peptidoform.properties["c_term"]:
            pass  # provide if c-term mods are supported
        else:
            tokenized_seq.append("<END>")

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
        df.rename({"ccs_observed": "ccs"}, axis=1, inplace=True)
        high_conf_hits = list(df["spectrum_id"][df["score"].rank(pct=True) > 0.95])
        logger.debug(
            f"Number of high confidence hits for calculating shift: {len(high_conf_hits)}"
        )

        shift_factor = get_ccs_shift(
            df[["charge", "sequence-tokenized", "ccs"]][df["spectrum_id"].isin(high_conf_hits)],
            self.reference_dataset,
        )

        logger.debug(f"CCS shift factor: {shift_factor}")

        return shift_factor

    def _is_valid_tokenized_sequence(self, tokenized_seq):
        """
        Check if peptide sequence contains invalid tokens.

        Args:
            tokenized_seq (list): Tokenized peptide sequence.

        Returns:
            bool: False if invalid tokens are present, True otherwise.
        """
        for token in tokenized_seq:
            if token not in self.allowed_mods:
                logger.debug(f"Invalid modification found: {token}")
                return False
        return True
