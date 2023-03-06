"""DeepLC retention time-based feature generator."""

import logging
import os
from typing import Union, Optional
from pathlib import Path
from itertools import chain

import numpy as np
import pandas as pd
from collections import defaultdict
from psm_utils import PSMList
from psm_utils.io import peptide_record

from ms2rescore.feature_generators import FeatureGenerator
from ms2rescore.utils import infer_spectrum_path
from ms2rescore.parse_mgf import parse_mgf_title_rt
from ms2rescore.exceptions import MS2RescoreError

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
logger = logging.getLogger(__name__)


class DeepLCFeatureGenerator(FeatureGenerator):
    def __init__(
        self,
        config,
        higher_psm_score_better: bool = True,
        calibration_set_size: Optional[Union[int, float]] = 0.20,
        num_cpu: int = 1,
        *args, **kwargs
    ) -> None:
        """
        Generate DeepLC-based features for rescoring.

        Parameters
        ----------
        psm_list: psm_utils.PSMList
            List of PSMs to which to add features.
        higher_psm_score_better: bool
            Whether a higher PSM score (`psm_score` column in PEPREC) denotes a better
            score. (default: True)
        calibration_set_size: int or float
            Amount of best PSMs to use for DeepLC calibration. If this value is lower
            than the number of available PSMs, all PSMs will be used. (default: 0.20)
        num_cpu: {int, None}
            Number of processes to use in DeepLC

        """
        super().__init__(*args, **kwargs)
        self.config = config
        self.tmp_file_root = str(
            Path(self.config["ms2rescore"]["tmp_path"])
            / Path(self.config["ms2rescore"]["psm_file"]).stem
        )
        self.feature_names = []
        self.higher_psm_score_better = higher_psm_score_better
        self.calibration_set_size = calibration_set_size
        self.num_cpu = num_cpu

        self.deeplc_predictor = None
        self.selected_model = None

        self.feature_names = [
            "observed_retention_time",
            "predicted_retention_time",
            "rt_diff",
            "observed_retention_time_best",
            "predicted_retention_time_best",
            "rt_diff_best",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        """Add DeepLC-derived features to PSMs."""
        from deeplc import DeepLC
        logger.info("Adding DeepLC-derived features to PSMs.")

        # Get easy-access nested version of PSMList
        psm_dict = psm_list.get_psm_dict()
        peptide_rt_diff_dict = defaultdict(lambda: {"observed_retention_time_best": np.Inf, "predicted_retention_time_best": np.Inf, "rt_diff_best": np.Inf}) 
        
        # Run MS²PIP for each spectrum file
        for collection, runs in psm_dict.items():
            # Reset DeepLC predictor for each collection of runs
            self.deeplc_predictor = None
            self.selected_model = None
            for run, psms in runs.items():
                logger.info(f"Processing {run}")
                # Prepare PSM file

                psm_list_run = PSMList(
                    psm_list=list(chain.from_iterable(psms.values()))
                )

                if not all(psm_list["retention_time"]):
                    retention_times = psm_list_run["retention_time"]
                    spec_ids = psm_list_run["retention_time"]
                    retention_time_dict = self.read_rt_from_spectrum_file(run)

                    try:
                        # psm_list["retention_time"] = [rt if rt else retention_time_dict[spec_ids[i]] for i, rt in enumerate(retention_times)] 
                        psm_list_run["retention_time"] = [retention_time_dict[psm_id] for psm_id in psm_list_run["spectrum_id"]] # Probably faster to replace all
                    except KeyError:
                        raise MS2RescoreError("Could not find all map spectrum ids to retention times")

                psm_list_calibration = self.get_calibration_psms(psm_list_run)

                logger.debug("Calibrating DeepLC")
                self.deeplc_predictor = DeepLC(
                    split_cal=10,
                    n_jobs=self.num_cpu,
                    cnn_model=True,
                    verbose=False,
                    path_model=self.selected_model,
                    pygam_calibration=True,
                    deeplc_retrain=True,
                )
                self.deeplc_predictor.calibrate_preds(
                    seq_df=self._psm_list_to_deeplc_peprec(psm_list_calibration)
                )
                # Still calibrate for each run, but do not try out all model options.
                # Just use model that was selected based on first run
                if not self.selected_model:
                    self.selected_model = list(self.deeplc_predictor.model.keys())
                    logger.debug(
                        f"Selected DeepLC model {self.selected_model} based on "
                        "calibration of first run. Using this model (after new "
                        "calibrations) for the remaining runs."
                    )

                predictions = np.array(self.deeplc_predictor.make_preds(
                    seq_df=self._psm_list_to_deeplc_peprec(psm_list_run)
                ))
                observations = psm_list_run["retention_time"]
                rt_diffs_run = np.abs(predictions - observations)
                
                for i, psm in enumerate(psm_list_run):
                    psm["rescoring_features"].update(
                        {"observed_retention_time": observations[i],
                        "predicted_retention_time": predictions[i],
                        "rt_diff": rt_diffs_run[i]
                    })
                    peptide = psm.peptidoform.proforma.split("\\")[0] # remove charge
                    if peptide_rt_diff_dict[peptide]["rt_diff_best"] > rt_diffs_run[i]:
                        peptide_rt_diff_dict[peptide] = {"observed_retention_time_best": observations[i], "predicted_retention_time_best": predictions[i], "rt_diff_best": rt_diffs_run[i]}

        for psm in psm_list:
            psm["rescoring_features"].update(peptide_rt_diff_dict[psm.peptidoform.proforma.split("\\")[0]])


    def read_rt_from_spectrum_file(self, run_name):

        # Prepare spectrum filenames
        spectrum_filename = infer_spectrum_path(
            self.config["ms2rescore"]["spectrum_path"], run_name
        )

        return parse_mgf_title_rt(spectrum_filename)




    # TODO: Remove when DeepLC supports PSMList directly
    @staticmethod
    def _psm_list_to_deeplc_peprec(psm_list: PSMList) -> pd.DataFrame:
        peprec = peptide_record.to_dataframe(psm_list)
        peprec = peprec.rename(
            columns={
                "observed_retention_time": "tr",
                "peptide": "seq",
            }
        )[["tr", "seq", "modifications"]]
        return peprec

    def get_calibration_psms(self, psm_list: PSMList):
        """Get N best scoring target PSMs for calibration."""
        psm_list_targets = psm_list[~psm_list["is_decoy"]]
        n_psms = self.get_number_of_calibration_psms(psm_list_targets)
        indices = np.argsort(psm_list_targets["score"])
        indices = indices[-n_psms:] if self.higher_psm_score_better else indices[:n_psms]
        return psm_list[indices]

    def get_number_of_calibration_psms(self, psm_list):
        """Get number of calibration PSMs given `calibration_set_size` and total number of PSMs."""
        if isinstance(self.calibration_set_size, float):
            if not 0 < self.calibration_set_size <= 1:
                raise ValueError(
                    "If `calibration_set_size` is a float, it cannot be smaller than "
                    "or equal to 0 or larger than 1."
                )
            else:
                num_calibration_psms = round(
                    len(psm_list) * self.calibration_set_size
                )
        elif isinstance(self.calibration_set_size, int):
            if self.calibration_set_size > len(psm_list):
                logger.warning(
                    f"Requested number of calibration PSMs ({self.calibration_set_size}"
                    f") is larger than total number of PSMs ({len(psm_list)}). Using "
                    "all PSMs for calibration."
                )
                num_calibration_psms = len(psm_list)
            else:
                num_calibration_psms = self.calibration_set_size
        else:
            raise TypeError(
                "Expected float or int for `calibration_set_size`. Got "
                f"{type(self.calibration_set_size)} instead. "
            )
        logger.debug(f"Using {num_calibration_psms} PSMs for calibration")
        return num_calibration_psms
