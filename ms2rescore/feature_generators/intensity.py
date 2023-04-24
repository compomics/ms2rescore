"""MS²PIP fragmentation intensity-based feature generator."""

import logging
import multiprocessing
import warnings
from itertools import chain
from pathlib import Path

import numpy as np
from ms2pip.core import correlate
from psm_utils import PSMList
from psm_utils.io import peptide_record
from rich.progress import track
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error as mse
from collections import defaultdict

from ms2rescore.exceptions import MS2RescoreError
from ms2rescore.feature_generators import FeatureGenerator
from ms2rescore.utils import infer_spectrum_path

logger = logging.getLogger(__name__)


class MS2PIPFeatureGenerator(FeatureGenerator):
    def __init__(self, config, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.config = config
        self.tmp_file_root = str(
            Path(self.config["ms2rescore"]["tmp_path"])
            / Path(self.config["ms2rescore"]["psm_file"]).stem
        )
        self.feature_names = [
            "spec_pearson_norm",
            "ionb_pearson_norm",
            "iony_pearson_norm",
            "spec_mse_norm",
            "ionb_mse_norm",
            "iony_mse_norm",
            "min_abs_diff_norm",
            "max_abs_diff_norm",
            "abs_diff_Q1_norm",
            "abs_diff_Q2_norm",
            "abs_diff_Q3_norm",
            "mean_abs_diff_norm",
            "std_abs_diff_norm",
            "ionb_min_abs_diff_norm",
            "ionb_max_abs_diff_norm",
            "ionb_abs_diff_Q1_norm",
            "ionb_abs_diff_Q2_norm",
            "ionb_abs_diff_Q3_norm",
            "ionb_mean_abs_diff_norm",
            "ionb_std_abs_diff_norm",
            "iony_min_abs_diff_norm",
            "iony_max_abs_diff_norm",
            "iony_abs_diff_Q1_norm",
            "iony_abs_diff_Q2_norm",
            "iony_abs_diff_Q3_norm",
            "iony_mean_abs_diff_norm",
            "iony_std_abs_diff_norm",
            "dotprod_norm",
            "dotprod_ionb_norm",
            "dotprod_iony_norm",
            "cos_norm",
            "cos_ionb_norm",
            "cos_iony_norm",
            "spec_pearson",
            "ionb_pearson",
            "iony_pearson",
            "spec_spearman",
            "ionb_spearman",
            "iony_spearman",
            "spec_mse",
            "ionb_mse",
            "iony_mse",
            "min_abs_diff_iontype",
            "max_abs_diff_iontype",
            "min_abs_diff",
            "max_abs_diff",
            "abs_diff_Q1",
            "abs_diff_Q2",
            "abs_diff_Q3",
            "mean_abs_diff",
            "std_abs_diff",
            "ionb_min_abs_diff",
            "ionb_max_abs_diff",
            "ionb_abs_diff_Q1",
            "ionb_abs_diff_Q2",
            "ionb_abs_diff_Q3",
            "ionb_mean_abs_diff",
            "ionb_std_abs_diff",
            "iony_min_abs_diff",
            "iony_max_abs_diff",
            "iony_abs_diff_Q1",
            "iony_abs_diff_Q2",
            "iony_abs_diff_Q3",
            "iony_mean_abs_diff",
            "iony_std_abs_diff",
            "dotprod",
            "dotprod_ionb",
            "dotprod_iony",
            "cos",
            "cos_ionb",
            "cos_iony",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        """Add MS²PIP-derived features to PSMs."""
        logger.info("Adding MS²PIP-derived features to PSMs.")

        """Generate rescoring features and add to psm_list."""
        # Prepare MS²PIP configuration
        self.config["ms2pip"].update(
            {
                "ptm": self._get_modification_config(psm_list),
                "sptm": [],
                "gptm": [],
            }
        )

        # Get easy-access nested version of PSMList
        psm_dict = psm_list.get_psm_dict()
        # Run MS²PIP for each spectrum file
        for collection, runs in psm_dict.items():
            for run, psms in runs.items():
                # Prepare PSM file
                logger.info(f"Preparing MS²PIP PSM file for {run}.")
                psm_list_run = PSMList(
                    psm_list=list(chain.from_iterable(psms.values()))
                )
                psm_id_mapper = {i: psm for i, psm in enumerate(psm_list_run)}
                # Prepare spectrum filenames
                spectrum_filename = infer_spectrum_path(
                    self.config["ms2rescore"]["spectrum_path"], run
                )
                output_filename = "-".join(
                    [self.tmp_file_root, str(collection), str(run)]
                )  # TODO use this if files have to be written

                # Run MS²PIP
                ms2pip_processing_results = correlate(
                    psms=psm_list_run,
                    spectrum_file=spectrum_filename,
                    spectrum_id_pattern=self.config["ms2rescore"][
                        "spectrum_id_pattern"
                    ],
                    model=self.config["ms2pip"]["model"],
                    compute_correlations=False,
                    processes=self.config["ms2rescore"]["num_cpu"],
                )

                # Calculate features
                logger.debug("Calculating features from predicted spectra")
                features = self._calculate_features(
                    ms2pip_processing_results,
                    self.config["ms2rescore"]["num_cpu"],
                    show_progress_bar=True,
                )
                # Add features to PSMs
                for psm_id, psm in psm_id_mapper.items():
                    try:
                        psm["rescoring_features"].update(features[psm_id])
                    except TypeError:  # if none we get type error
                        psm["rescoring_features"] = None

    @staticmethod
    def _get_modification_config(psm_list):
        """
        Get MS²PIP-style modification configuration from PSMList.

        Notes
        -----
        Fixed, labile, and unlocalized modifications are ignored. Fixed modifications
        should therefore already have been applied.

        """
        unique_modifications = set()
        for psm in psm_list:
            for aa, mods in psm.peptidoform.parsed_sequence:
                if mods:
                    unique_modifications.update([(aa, mod) for mod in mods])
            if psm.peptidoform.properties["n_term"]:
                unique_modifications.update(
                    [("N-term", mod) for mod in psm.peptidoform.properties["n_term"]]
                )
            if psm.peptidoform.properties["c_term"]:
                unique_modifications.update(
                    [("C-term", mod) for mod in psm.peptidoform.properties["c_term"]]
                )
        return [
            ",".join([mod.value, str(mod.mass), "opt", target])
            for target, mod in unique_modifications
        ]

    def _compute_features(self, processing_result):
        """Compute features from observed and predicted intensities."""
        # Suppress RuntimeWarnings about invalid values
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            if (
                processing_result.observed_intensity is None
                or processing_result.predicted_intensity is None
            ):
                return (processing_result.psm_index, None)

            # Convert intensities to arrays
            target_b = processing_result.predicted_intensity["b"].clip(np.log2(0.001))
            target_y = processing_result.predicted_intensity["y"].clip(np.log2(0.001))
            target_all = np.concatenate([target_b, target_y])
            prediction_b = processing_result.observed_intensity["b"].clip(
                np.log2(0.001)
            )
            prediction_y = processing_result.observed_intensity["y"].clip(
                np.log2(0.001)
            )
            prediction_all = np.concatenate([prediction_b, prediction_y])

            # Prepare 'unlogged' intensity arrays
            target_b_unlog = 2**target_b - 0.001
            target_y_unlog = 2**target_y - 0.001
            target_all_unlog = 2**target_all - 0.001
            prediction_b_unlog = 2**prediction_b - 0.001
            prediction_y_unlog = 2**prediction_y - 0.001
            prediction_all_unlog = 2**prediction_all - 0.001

            # Calculate absolute differences
            abs_diff_b = np.abs(target_b - prediction_b)
            abs_diff_y = np.abs(target_y - prediction_y)
            abs_diff_all = np.abs(target_all - prediction_all)
            abs_diff_b_unlog = np.abs(target_b_unlog - prediction_b_unlog)
            abs_diff_y_unlog = np.abs(target_y_unlog - prediction_y_unlog)
            abs_diff_all_unlog = np.abs(target_all_unlog - prediction_all_unlog)

            # Compute features
            feature_values = [
                # Features between spectra in log space
                pearsonr(target_all, prediction_all)[0],  # Pearson all ions
                pearsonr(target_b, prediction_b)[0],  # Pearson b ions
                pearsonr(target_y, prediction_y)[0],  # Pearson y ions
                mse(target_all, prediction_all),  # MSE all ions
                mse(target_b, prediction_b),  # MSE b ions
                mse(target_y, prediction_y),  # MSE y ions
                np.min(abs_diff_all),  # min_abs_diff_norm
                np.max(abs_diff_all),  # max_abs_diff_norm
                np.quantile(abs_diff_all, 0.25),  # abs_diff_Q1_norm
                np.quantile(abs_diff_all, 0.5),  # abs_diff_Q2_norm
                np.quantile(abs_diff_all, 0.75),  # abs_diff_Q3_norm
                np.mean(abs_diff_all),  # mean_abs_diff_norm
                np.std(abs_diff_all),  # std_abs_diff_norm
                np.min(abs_diff_b),  # ionb_min_abs_diff_norm
                np.max(abs_diff_b),  # ionb_max_abs_diff_norm
                np.quantile(abs_diff_b, 0.25),  # ionb_abs_diff_Q1_norm
                np.quantile(abs_diff_b, 0.5),  # ionb_abs_diff_Q2_norm
                np.quantile(abs_diff_b, 0.75),  # ionb_abs_diff_Q3_norm
                np.mean(abs_diff_b),  # ionb_mean_abs_diff_norm
                np.std(abs_diff_b),  # ionb_std_abs_diff_norm
                np.min(abs_diff_y),  # iony_min_abs_diff_norm
                np.max(abs_diff_y),  # iony_max_abs_diff_norm
                np.quantile(abs_diff_y, 0.25),  # iony_abs_diff_Q1_norm
                np.quantile(abs_diff_y, 0.5),  # iony_abs_diff_Q2_norm
                np.quantile(abs_diff_y, 0.75),  # iony_abs_diff_Q3_norm
                np.mean(abs_diff_y),  # iony_mean_abs_diff_norm
                np.std(abs_diff_y),  # iony_std_abs_diff_norm
                np.dot(target_all, prediction_all),  # Dot product all ions
                np.dot(target_b, prediction_b),  # Dot product b ions
                np.dot(target_y, prediction_y),  # Dot product y ions
                np.dot(target_all, prediction_all)
                / (
                    np.linalg.norm(target_all, 2) * np.linalg.norm(prediction_all, 2)
                ),  # Cos similarity all ions
                np.dot(target_b, prediction_b)
                / (
                    np.linalg.norm(target_b, 2) * np.linalg.norm(prediction_b, 2)
                ),  # Cos similarity b ions
                np.dot(target_y, prediction_y)
                / (
                    np.linalg.norm(target_y, 2) * np.linalg.norm(prediction_y, 2)
                ),  # Cos similarity y ions
                # Same features in normal space
                pearsonr(target_all_unlog, prediction_all_unlog)[0],  # Pearson all ions
                pearsonr(target_b_unlog, prediction_b_unlog)[0],  # Pearson b ions
                pearsonr(target_y_unlog, prediction_y_unlog)[0],  # Pearson y ions
                spearmanr(target_all_unlog, prediction_all_unlog)[
                    0
                ],  # Spearman all ions
                spearmanr(target_b_unlog, prediction_b_unlog)[0],  # Spearman b ions
                spearmanr(target_y_unlog, prediction_y_unlog)[0],  # Spearman y ions
                mse(target_all_unlog, prediction_all_unlog),  # MSE all ions
                mse(target_b_unlog, prediction_b_unlog),  # MSE b ions
                mse(target_y_unlog, prediction_y_unlog),  # MSE y ions,
                0
                if np.min(abs_diff_b_unlog) <= np.min(abs_diff_y_unlog)
                else 1,  # Ion type with min absolute difference
                0
                if np.max(abs_diff_b_unlog) >= np.max(abs_diff_y_unlog)
                else 1,  # Ion type with max absolute difference
                np.min(abs_diff_all_unlog),  # min_abs_diff
                np.max(abs_diff_all_unlog),  # max_abs_diff
                np.quantile(abs_diff_all_unlog, 0.25),  # abs_diff_Q1
                np.quantile(abs_diff_all_unlog, 0.5),  # abs_diff_Q2
                np.quantile(abs_diff_all_unlog, 0.75),  # abs_diff_Q3
                np.mean(abs_diff_all_unlog),  # mean_abs_diff
                np.std(abs_diff_all_unlog),  # std_abs_diff
                np.min(abs_diff_b_unlog),  # ionb_min_abs_diff
                np.max(abs_diff_b_unlog),  # ionb_max_abs_diff_norm
                np.quantile(abs_diff_b_unlog, 0.25),  # ionb_abs_diff_Q1
                np.quantile(abs_diff_b_unlog, 0.5),  # ionb_abs_diff_Q2
                np.quantile(abs_diff_b_unlog, 0.75),  # ionb_abs_diff_Q3
                np.mean(abs_diff_b_unlog),  # ionb_mean_abs_diff
                np.std(abs_diff_b_unlog),  # ionb_std_abs_diff
                np.min(abs_diff_y_unlog),  # iony_min_abs_diff
                np.max(abs_diff_y_unlog),  # iony_max_abs_diff
                np.quantile(abs_diff_y_unlog, 0.25),  # iony_abs_diff_Q1
                np.quantile(abs_diff_y_unlog, 0.5),  # iony_abs_diff_Q2
                np.quantile(abs_diff_y_unlog, 0.75),  # iony_abs_diff_Q3
                np.mean(abs_diff_y_unlog),  # iony_mean_abs_diff
                np.std(abs_diff_y_unlog),  # iony_std_abs_diff
                np.dot(target_all_unlog, prediction_all_unlog),  # Dot product all ions
                np.dot(target_b_unlog, prediction_b_unlog),  # Dot product b ions
                np.dot(target_y_unlog, prediction_y_unlog),  # Dot product y ions
                np.dot(target_all_unlog, prediction_all_unlog)
                / (
                    np.linalg.norm(target_all_unlog, 2)
                    * np.linalg.norm(prediction_all_unlog, 2)
                ),  # Cos similarity all ions
                np.dot(target_b_unlog, prediction_b_unlog)
                / (
                    np.linalg.norm(target_b_unlog, 2)
                    * np.linalg.norm(prediction_b_unlog, 2)
                ),  # Cos similarity b ions
                np.dot(target_y_unlog, prediction_y_unlog)
                / (
                    np.linalg.norm(target_y_unlog, 2)
                    * np.linalg.norm(prediction_y_unlog, 2)
                ),  # Cos similarity y ions
            ]

        features = dict(
            zip(
                self.feature_names,
                [0.0 if ft is np.nan else ft for ft in feature_values],
            )
        )
        return (processing_result.psm_index, features)

    def _calculate_features(
        self, processing_results, num_cpu=1, show_progress_bar=True
    ):
        """Calculate MS²PIP-based features in parallelized fashion."""
        logger.debug("Computing features")

        # Do not use multiprocessing for small amount of features
        if len(processing_results) < 10000:
            feature_result_list = [
                self._compute_features(result) for result in processing_results
            ]
            all_features = {
                psm_id: features for psm_id, features in feature_result_list
            }
        else:
            # Split up df into list of chunk_size df's (will be divided over num_cpu)
            # Use imap, so we can use a progress bar
            with multiprocessing.Pool(int(num_cpu)) as pool:
                feature_result_list = track(
                    pool.imap(
                        self._compute_features, processing_results, chunksize=100
                    ),
                    total=len(processing_results),
                    disable=True if not show_progress_bar else False,
                    description="Calculating features...",
                    transient=True,
                )
                all_features = {
                    psm_id: features for psm_id, features in feature_result_list
                }
        return all_features
