"""
Placeholder
"""

import contextlib
import logging
import os
from collections import defaultdict
from inspect import getfullargspec
from itertools import chain
from typing import List, Optional, Union

# Hardcoded for now because only one model is available and not yet in a separate package
im2deep_dir = os.path.dirname(os.path.realpath(__file__))
# TODO: put in im2deep package
DEFAULT_MODELS = [
    "../im2deep_data/mods/full_hc_peprec_CCS_1fd8363d9af9dcad3be7553c39396960.hdf5",
    "../im2deep_data/mods/full_hc_peprec_CCS_8c22d89667368f2f02ad996469ba157e.hdf5",
    "../im2deep_data/mods/full_hc_peprec_CCS_cb975cfdd4105f97efa0b3afffe075cc.hdf5",
]
DEFAULT_MODELS = [os.path.join(im2deep_dir, mod) for mod in DEFAULT_MODELS]

import numpy as np
import pandas as pd
from psm_utils.peptidoform import Peptidoform
from psm_utils import PSMList
from psm_utils.io import peptide_record
from numpy import ndarray

# TODO: put in im2deep package
DEFAULT_REFERENCE_DATASET = pd.read_csv(
    os.path.join(im2deep_dir, "../im2deep_data/reference_ccs.csv")
)

from ms2rescore.feature_generators.base import FeatureGeneratorBase

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
logger = logging.getLogger(__name__)


class IM2DeepFeatureGenerator(FeatureGeneratorBase):
    """IM2Deep collision cross section feature generator."""

    def __init__(
        self,
        *args,
        lower_score_is_better: bool = False,
        reference_dataset: Optional[str] = None,
        spectrum_path: Optional[str] = None,
        processes: int = 1,
        **kwargs,
    ):
        """Placeholder"""  # TODO
        super().__init__(*args, **kwargs)
        self.lower_score_is_better = lower_score_is_better
        self.spectrum_path = spectrum_path
        self.processes = processes
        self.deeplc_kwargs = kwargs or {}

        self._verbose = logger.getEffectiveLevel() <= logging.DEBUG

        # Lazy-load DeepLC
        from deeplc import DeepLC

        self.im2deep = DeepLC

        # Remove any kwargs that are not DeepLC arguments
        self.im2deep_kwargs = {
            k: v for k, v in self.deeplc_kwargs.items() if k in getfullargspec(DeepLC).args
        }
        self.im2deep_kwargs.update({"config_file": None})

        # TODO: Implement deeplc_retrain?

        self.im2deep_predictor = None
        self.im2deep_model = DEFAULT_MODELS
        self.reference_dataset = (
            reference_dataset if reference_dataset else DEFAULT_REFERENCE_DATASET
        )

    # TODO name differently than ionmob so it doesnt get overwritten
    @property
    def feature_names(self) -> List[str]:
        return [
            "ccs_observed_im2deep",
            "ccs_predicted_im2deep",
            "ccs_error_im2deep",
            "abs_ccs_error_im2deep",
            "perc_ccs_error_im2deep",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        """Add IM2Deep-derived features to PSMs"""

        logger.info("Adding IM2Deep-derived features to PSMs")

        # Get easy-access nested version of PSMlist
        psm_dict = psm_list.get_psm_dict()

        # Run IM2Deep for each spectrum file
        current_run = 1
        total_runs = sum(len(runs) for runs in psm_dict.values())

        for runs in psm_dict.values():
            # Reset IM2Deep predictor for each collection of runs
            self.im2deep_predictor = None
            self.selected_model = None
            for run, psms in runs.items():
                logger.info(
                    f"Running IM2Deep for PSMs from run ({current_run}/{total_runs}): `{run}`..."
                )

                # Disable wild logging to stdout by TensorFlow, unless in debug mode
                with contextlib.redirect_stdout(
                    open(os.devnull, "w")
                ) if not self._verbose else contextlib.nullcontext():
                    # Make new PSM list for this run (chain PSMs per spectrum to flat list)
                    psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))

                    logger.debug("Calibrating IM2Deep...")

                    self.im2deep_predictor = self.im2deep(
                        n_jobs=self.processes,
                        verbose=self._verbose,
                        path_model=self.im2deep_model,
                        # TODO: only one model as of now, so these models should be available somewhere else, I guess DeepLC for now until IM2Deep becomes its own package?
                        predict_ccs=True,
                        **self.im2deep_kwargs,
                    )

                    # Convert ion mobility to CCS and calibrate CCS values
                    psm_list_run_df = psm_list_run.to_dataframe()
                    psm_list_run_df["charge"] = [
                        peptidoform.precursor_charge
                        for peptidoform in psm_list_run_df["peptidoform"]
                    ]
                    psm_list_run_df["ccs_observed"] = psm_list_run_df.apply(
                        lambda x: self.im2ccs(
                            x["ion_mobility"],
                            x["precursor_mz"],  # TODO: Why does ionmob use calculated mz?
                            x["charge"],
                        ),
                        axis=1,
                    )

                    shift_factor = self.calculate_ccs_shift(psm_list_run_df)
                    psm_list_run_df["ccs_observed"] = psm_list_run_df.apply(
                        lambda x: x["ccs_observed"] + shift_factor, axis=1
                    )

                    # TODO: probably not required since only one model as of now
                    if not self.selected_model:
                        # self.selected_model = list(self.im2deep_predictor.model.keys())
                        self.selected_model = self.im2deep_predictor.model
                        self.im2deep_kwargs["deeplc_retrain"] = False
                        logger.debug(
                            f"Selected model: {self.selected_model}. Using this model (after new "
                            "calibrations) for all subsequent runs."
                        )

                    logger.debug("Predicting CCS values...")
                    predictions = np.array(
                        self.im2deep_predictor.make_preds(
                            seq_df=self._psm_list_df_to_deeplc_peprec(psm_list_run_df),
                            calibrate=False,
                        )
                    )

                    observations = psm_list_run_df["ccs_observed"]
                    ccs_diffs_run = np.abs(predictions - observations)

                    logger.debug("Adding features to PSMs...")

                    for i, psm in enumerate(psm_list_run):
                        psm["rescoring_features"].update(
                            {
                                "ccs_observed_im2deep": observations[i],
                                "ccs_predicted_im2deep": predictions[i],
                                "ccs_error_im2deep": ccs_diffs_run[i],
                                "abs_ccs_error_im2deep": np.abs(ccs_diffs_run[i]),
                                "perc_ccs_error_im2deep": np.abs(ccs_diffs_run[i])
                                / observations[i]
                                * 100,
                            }
                        )

                current_run += 1

    # TODO: Remove when either DeepLC or psm_list supports 'observed_ccs' directly
    @staticmethod
    def _psm_list_df_to_deeplc_peprec(psm_list_df: pd.DataFrame) -> pd.DataFrame:
        # TODO: As for now, DeepLC uses the column 'tr' as observed CCS. This should be changed
        psm_list_df["seq"] = psm_list_df["peptidoform"].apply(lambda x: x.sequence)
        psm_list_df["charge"] = psm_list_df["peptidoform"].apply(lambda x: x.precursor_charge)
        psm_list_df["modifications"] = psm_list_df["peptidoform"].apply(
            lambda x: peptide_record.proforma_to_peprec(x)[1]
        )
        peprec = psm_list_df.rename(
            columns={
                "ccs_observed": "tr",
            }
        )[["tr", "charge", "seq", "modifications"]]
        return peprec

    # TODO: What to do when no overlapping peptide-charge pairs are found?
    def get_ccs_shift(
        self, df: pd.DataFrame, reference_dataset: pd.DataFrame, use_charge_state: int = 2
    ) -> ndarray:
        """
        Calculate CCS shift factor, i.e. a constant offset based on identical precursors as in reference.

        Parameters
        ----------
        df : pd.DataFrame
            PSMs with CCS values.
        reference_dataset : pd.DataFrame
            Reference dataset with CCS values.
        use_charge_state : int, optional
            Charge state to use for CCS shift calculation, needs to be [2,4], by default 2.
        return_shift_factor : float
            CCS shift factor.

        """
        if use_charge_state > 4 or use_charge_state < 2:
            raise ValueError("Charge state needs to be higher than 1 and lower than 5.")
        else:
            logger.debug(f"Using charge state {use_charge_state} for CCS shift calculation.")

        tmp_df = df.copy(deep=True)
        tmp_ref_df = reference_dataset.copy(deep=True)

        tmp_df["sequence"] = tmp_df["peptidoform"].apply(lambda x: x.proforma.split("\\")[0])
        tmp_df["charge"] = tmp_df["peptidoform"].apply(lambda x: x.precursor_charge)
        tmp_ref_df["sequence"] = tmp_ref_df["peptidoform"].apply(
            lambda x: Peptidoform(x).proforma.split("\\")[0]
        )
        tmp_ref_df["charge"] = tmp_ref_df["peptidoform"].apply(
            lambda x: Peptidoform(x).precursor_charge
        )

        reference_tmp = tmp_ref_df[tmp_ref_df["charge"] == use_charge_state]
        df_tmp = tmp_df[tmp_df["charge"] == use_charge_state]
        both = pd.merge(
            left=reference_tmp,
            right=df_tmp,
            right_on=["sequence", "charge"],
            left_on=["sequence", "charge"],
            how="inner",
            suffixes=("_ref", "_data"),
        )
        logger.debug(
            "Calculating CCS shift based on {} overlapping peptide-charge pairs found between PSMs and reference dataset".format(
                both.shape[0]
            )
        )

        return 0 if both.shape[0] == 0 else np.mean(both["CCS"] - both["ccs_observed"])

    def calculate_ccs_shift(self, psm_df: pd.DataFrame) -> float:
        """
        Apply CCS shift to CCS values.

        Parameters
        ----------
        psm_list : PSMList

        """
        psm_df = psm_df[psm_df["charge"] < 5]  # predictions do not go higher for IM2Deep
        high_conf_hits = list(psm_df["spectrum_id"][psm_df["score"].rank(pct=True) > 0.95])
        logger.debug(
            f"Number of high confidence hits for calculating shift: {len(high_conf_hits)}"
        )

        # Filter df for high_conf_hits
        psm_df = psm_df[psm_df["spectrum_id"].isin(high_conf_hits)]
        shift_factor = self.get_ccs_shift(
            psm_df,
            self.reference_dataset,
        )

        logger.debug(f"CCS shift factor: {shift_factor}")

        return shift_factor

    # TODO ion mobility to CCS conversion? Is this necessary?
    def im2ccs(self, reverse_im, mz, charge, mass_gas=28.013, temp=31.85, t_diff=273.15):
        """
        Convert ion mobility to CCS.  #TODO: Took this from ionmob. how to reference?

        Parameters
        ----------
        reverse_im : float
            Reduced ion mobility.
        mz : float
            Precursor m/z.
        charge : int
            Precursor charge.
        mass_gas : float, optional
            Mass of gas, by default 28.013
        temp : float, optional
            Temperature in Celsius, by default 31.85
        t_diff : float, optional
            Factor to convert Celsius to Kelvin, by default 273.15"""

        SUMMARY_CONSTANT = 18509.8632163405
        reduced_mass = (mz * charge * mass_gas) / (mz * charge + mass_gas)
        return (SUMMARY_CONSTANT * charge) / (
            np.sqrt(reduced_mass * (temp + t_diff)) * 1 / reverse_im
        )
