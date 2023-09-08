"""
Feature generator for PSMs from the MaxQuant search engine.

MaxQuant msms.txt files contain various metrics from peptide-spectrum matching that can be used
to generate rescoring features. These include features related to the mass errors of the seven
fragment ions with the highest intensities, and features related to the ion current of the
identified fragment ions.

"""

import logging
from typing import List, Tuple

import numpy as np
from psm_utils import PSMList

from ms2rescore.exceptions import MS2RescoreError
from ms2rescore.feature_generators.base import FeatureGeneratorBase

logger = logging.getLogger(__name__)


class MaxQuantFeatureGenerator(FeatureGeneratorBase):
    """Generate MaxQuant-derived features."""

    def __init__(self, *args, **kwargs) -> None:
        """
        Generate MaxQuant-derived features.

        Attributes
        ----------
        feature_names: list[str]
            Names of the features that will be added to the PSMs.

        Raises
        ------
        MissingMetadataError
            If the required metadata entries are not present in the PSMs.

        """
        super().__init__(*args, **kwargs)

    @property
    def feature_names(self) -> List[str]:
        return [
            "mean_error_top7",
            "sq_mean_error_top7",
            "stdev_error_top7",
            "ln_explained_ion_current",
            "ln_nterm_ion_current_ratio",
            "ln_cterm_ion_current_ratio",
            "ln_ms2_ion_current",
        ]

    def add_features(self, psm_list: PSMList):
        """
        Add MS²PIP-derived features to PSMs.

        Parameters
        ----------
        psm_list
            PSMs to add features to.

        """
        logger.info("Adding MaxQuant-derived features to PSMs.")

        # Infer mass deviations column name
        for column_name in [
            "Mass deviations [Da]",
            "Mass Deviations [Da]",
            "Mass deviations [ppm]",
            "Mass Deviations [ppm]",
        ]:
            if column_name in psm_list[0]["metadata"].keys():
                self._mass_deviations_key = column_name
                break
        else:
            raise MissingMetadataError(
                "No mass deviations entry in PSM metadata. Cannot compute MaxQuant features."
            )

        # Check other columns
        for column_name in ["Intensities", "Matches", "Intensity coverage"]:
            if column_name not in psm_list[0]["metadata"].keys():
                raise MissingMetadataError(
                    f"Missing {column_name} entry in PSM metadata. Cannot compute MaxQuant features."
                )

        # Add features to PSMs
        for psm in psm_list:
            psm["rescoring_features"].update(self._compute_features(psm["metadata"]))

    def _compute_features(self, psm_metadata):
        """Compute features from derived from intensities and mass errors."""
        features = {}
        if all(k in psm_metadata.keys() for k in ["Intensities", self._mass_deviations_key]):
            (
                features["mean_error_top7"],
                features["sq_mean_error_top7"],
                features["stdev_error_top7"],
            ) = self._calculate_top7_peak_features(
                psm_metadata["Intensities"], psm_metadata[self._mass_deviations_key]
            )

        if all(k in psm_metadata.keys() for k in ["Intensities", "Matches", "Intensity coverage"]):
            (
                features["ln_explained_ion_current"],
                features["ln_nterm_ion_current_ratio"],
                features["ln_cterm_ion_current_ratio"],
                features["ln_ms2_ion_current"],
            ) = self._calculate_ion_current_features(
                psm_metadata["Matches"],
                psm_metadata["Intensities"],
                psm_metadata["Intensity coverage"],
            )

        return features

    @staticmethod
    def _calculate_top7_peak_features(intensities: str, mass_errors: str) -> Tuple[np.ndarray]:
        """
        Calculate "top 7 peak"-related search engine features.
        The following features are calculated:
        - mean_error_top7: Mean of mass errors of the seven fragment ion peaks with the
          highest intensities
        - sq_mean_error_top7: Squared MeanErrorTop7
        - stdev_error_top7: Standard deviation of mass errors of the seven fragment ion
          peaks with the highest intensities
        """
        try:
            intensities = [float(i) for i in intensities.split(";")]
            mass_errors = [float(i) for i in mass_errors.split(";")]
        except ValueError:
            return 0.0, 0.0, 0.0

        indices_most_intens = np.array(intensities).argsort()[-1:-8:-1]
        mass_errors_top7 = [(mass_errors[i]) for i in indices_most_intens]
        mean_error_top7 = np.mean(mass_errors_top7)
        sq_mean_error_top7 = mean_error_top7**2
        stdev_error_top7 = np.std(mass_errors_top7)

        return mean_error_top7, sq_mean_error_top7, stdev_error_top7

    @staticmethod
    def _calculate_ion_current_features(
        matches: str, intensities: str, intensity_coverage: str
    ) -> Tuple[np.ndarray]:
        """
        Calculate ion current related search engine features.
        The following features are calculated:
        - ln_explained_ion_current: Summed intensity of identified fragment ions,
          divided by that of all fragment ions, logged
        - ln_nterm_ion_current_ratio: Summed intensity of identified N-terminal
          fragments, divided by that of all identified fragments, logged
        - ln_cterm_ion_current_ratio: Summed intensity of identified N-terminal
          fragments, divided by that of all identified fragments, logged
        - ln_ms2_ion_current: Summed intensity of all observed fragment ions, logged
        """
        pseudo_count = 0.00001
        try:
            ln_explained_ion_current = float(intensity_coverage) + pseudo_count
            summed_intensities = sum([float(i) for i in intensities.split(";")])
        except ValueError:
            return 0.0, 0.0, 0.0, 0.0

        # Calculate ratio between matched b- and y-ion intensities
        y_ion_int = sum(
            [
                float(intensities.split(";")[i])
                for i, m in enumerate(matches.split(";"))
                if m.startswith("y")
            ]
        )
        y_int_ratio = y_ion_int / summed_intensities

        ln_nterm_ion_current_ratio = (y_int_ratio + pseudo_count) * ln_explained_ion_current
        ln_cterm_ion_current_ratio = (1 - y_int_ratio + pseudo_count) * ln_explained_ion_current
        ln_ms2_ion_current = summed_intensities / ln_explained_ion_current

        out = [
            ln_explained_ion_current,
            ln_nterm_ion_current_ratio,
            ln_cterm_ion_current_ratio,
            ln_ms2_ion_current,
        ]

        return tuple([np.log(x) for x in out])


class MissingMetadataError(MS2RescoreError):
    """Exception raised when a required metadata entry is missing."""

    pass
