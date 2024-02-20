"""
IM2Deep ion mobility-based feature generator.

IM2Deep is a fully modification-aware peptide ion mobility predictor. It uses a deep convolutional
neural network to predict retention times based on the atomic composition of the (modified) amino
acid residues in the peptide. See
`github.com/compomics/deeplc <https://github.com/compomics/IM2Deep>`_ for more information.

"""

import contextlib
import logging
import os
from inspect import getfullargspec
from itertools import chain
from typing import List, Optional

import numpy as np
from im2deep.im2deep import predict_ccs
from psm_utils import PSMList

from ms2rescore.feature_generators.base import FeatureGeneratorBase

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
logger = logging.getLogger(__name__)


class IM2DeepFeatureGenerator(FeatureGeneratorBase):
    """IM2Deep collision cross section feature generator."""

    def __init__(
        self,
        *args,
        lower_score_is_better: bool = False,
        spectrum_path: Optional[str] = None,
        processes: int = 1,
        calibrate_per_charge: bool = True,
        **kwargs,
    ):
        """
        Initialize the IM2DeepFeatureGenerator.

        Parameters
        ----------
        lower_score_is_better : bool, optional
            A boolean indicating whether lower scores are better for the generated features.
        spectrum_path : str or None, optional
            Optional path to the spectrum file used for IM2Deep predictions.
        processes : int, optional
            Number of parallel processes to use for IM2Deep predictions.
        calibrate_per_charge : bool, optional
            A boolean indicating whether to calibrate CCS values per charge state.
        **kwargs : dict, optional
            Additional keyword arguments.

        Returns
        -------
        None
        """
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

        # TODO: Implement im2deep_retrain?

        self.im2deep_predictor = None
        self.calibrate_per_charge = calibrate_per_charge

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
                with (
                    contextlib.redirect_stdout(open(os.devnull, "w"))
                    if not self._verbose
                    else contextlib.nullcontext()
                ):
                    # Make new PSM list for this run (chain PSMs per spectrum to flat list)
                    psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))

                    logger.debug("Calibrating IM2Deep...")

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

                    # Create dataframe with high confidence hits for calibration
                    cal_psm_df = self.make_cal_df(psm_list_run_df)

                    # Make predictions with IM2Deep
                    logger.debug("Predicting CCS values...")
                    calibrated_predictions = predict_ccs(
                        psm_list_run, cal_psm_df, write_output=False
                    )

                    # Add features to PSMs
                    logger.debug("Adding features to PSMs...")
                    predictions = calibrated_predictions
                    observations = psm_list_run_df["ccs_observed"]
                    ccs_diffs_run = np.abs(predictions - observations)
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

    def im2ccs(self, reverse_im, mz, charge, mass_gas=28.013, temp=31.85, t_diff=273.15):
        """
        Convert ion mobility to CCS.

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
            Factor to convert Celsius to Kelvin, by default 273.15

        Notes
        -----
        Adapted from theGreatHerrLebert/ionmob (https://doi.org/10.1093/bioinformatics/btad486)

        """

        SUMMARY_CONSTANT = 18509.8632163405
        reduced_mass = (mz * charge * mass_gas) / (mz * charge + mass_gas)
        return (SUMMARY_CONSTANT * charge) / (
            np.sqrt(reduced_mass * (temp + t_diff)) * 1 / reverse_im
        )

    # TODO: replace threshold by identified psms?
    def make_cal_df(self, psm_list_df, threshold=0.95):
        """Make dataframe for calibration of IM2Deep predictions.

        Parameters
        ----------
        psm_list_df : pd.DataFrame
            DataFrame with PSMs.
        threshold : float, optional
            Threshold for high confidence hits, by default 0.95.

        Returns
        -------
        pd.DataFrame
            DataFrame with high confidence hits for calibration."""

        psm_list_df = psm_list_df[
            psm_list_df["charge"] < 5
        ]  # predictions do not go higher for IM2Deep
        high_conf_hits = list(
            psm_list_df["spectrum_id"][psm_list_df["score"].rank(pct=True) > threshold]
        )
        logger.debug(
            f"Number of high confidence hits for calculating shift: {len(high_conf_hits)}"
        )
        # Filter df for high_conf_hits
        cal_psm_df = psm_list_df[psm_list_df["spectrum_id"].isin(high_conf_hits)]
        return cal_psm_df
