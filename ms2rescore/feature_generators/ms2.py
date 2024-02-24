"""
MS2-based feature generator.

"""

import logging
import re
from typing import List, Optional, Union
from itertools import chain
from collections import defaultdict


import numpy as np
from psm_utils import PSMList
from pyteomics.mass import calculate_mass

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from ms2rescore.utils import infer_spectrum_path

logger = logging.getLogger(__name__)


class MS2FeatureGenerator(FeatureGeneratorBase):
    """DeepLC retention time-based feature generator."""

    def __init__(
        self,
        *args,
        spectrum_path: Optional[str] = None,
        spectrum_id_pattern: str = "(.*)",
        **kwargs,
    ) -> None:
        """
        Generate MS2-based features for rescoring.

        DeepLC retraining is on by default. Add ``deeplc_retrain: False`` as a keyword argument to
        disable retraining.

        Parameters
        ----------

        Attributes
        ----------
        feature_names: list[str]
            Names of the features that will be added to the PSMs.

        """
        super().__init__(*args, **kwargs)
        self.spectrum_path = spectrum_path
        self.spectrum_id_pattern = spectrum_id_pattern

    @property
    def feature_names(self) -> List[str]:
        return [
            "ln_explained_intensity",
            "ln_total_intensity",
            "ln_explained_intensity_ratio",
            "ln_explained_b_ion_ratio",
            "ln_explained_y_ion_ratio",
            "longest_b_ion_sequence",
            "longest_y_ion_sequence",
            "matched_b_ions",
            "matched_b_ions_pct",
            "matched_y_ions",
            "matched_y_ions_pct",
            "matched_ions_pct",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        """Add DeepLC-derived features to PSMs."""

        logger.info("Adding MS2-derived features to PSMs.")
        psm_dict = psm_list.get_psm_dict()
        current_run = 1
        total_runs = sum(len(runs) for runs in psm_dict.values())

        for runs in psm_dict.values():
            for run, psms in runs.items():
                logger.info(
                    f"Running MS2 for PSMs from run ({current_run}/{total_runs}) `{run}`..."
                )
                psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
                spectrum_filename = infer_spectrum_path(self.spectrum_path, run)
                logger.debug(f"Using spectrum file `{spectrum_filename}`")
                self._add_features_run(psm_list_run, spectrum_filename)
                current_run += 1

    def _add_features_run(self, psm_list: PSMList, spectrum_filename: str) -> None:
        pass


def longest_sequence_of_trues(lst):
    max_sequence = 0
    current_sequence = 0

    for value in lst:
        current_sequence = current_sequence + 1 if value else 0
        max_sequence = max(max_sequence, current_sequence)

    return max_sequence


def extract_spectrum_features(annotated_spectrum, peptidoform):
    features = defaultdict(list)
    b_ions_matched = [False] * (len(peptidoform.sequence))
    y_ions_matched = [False] * (len(peptidoform.sequence))

    pseudo_count = 0.00001

    for annotated_peak in annotated_spectrum.spectrum:
        features["total_intensity"].append(annotated_peak.intensity)

        if annotated_peak.annotation:
            features["matched_intensity"].append(annotated_peak.intensity)
            for matched_ion in annotated_peak.annotation:
                if "y" in matched_ion.ion:
                    features["y_ion_matched"].append(annotated_peak.intensity)
                    y_ions_matched[int(re.search(r"\d+", matched_ion.ion).group())] = True
                elif "b" in matched_ion.ion:
                    features["b_ion_matched"].append(annotated_peak.intensity)
                    b_ions_matched[int(re.search(r"\d+", matched_ion.ion).group())] = True

    return {
        "ln_explained_intensity": np.log(np.sum(features["matched_intensity"]) + pseudo_count),
        "ln_total_intensity": np.log(np.sum(features["total_intensity"]) + pseudo_count),
        "ln_explained_intensity_ratio": np.log(
            np.sum(features["matched_intensity"]) / np.sum(features["total_intensity"])
            + pseudo_count
        ),
        "ln_explained_b_ion_ratio": np.log(
            np.sum(features["b_ion_matched"]) / np.sum(features["matched_intensity"])
            + pseudo_count
        ),
        "ln_explained_y_ion_ratio": np.log(
            np.sum(features["y_ion_matched"]) / np.sum(features["matched_intensity"])
            + pseudo_count
        ),
        "longest_b_ion_sequence": longest_sequence_of_trues(b_ions_matched),
        "longest_y_ion_sequence": longest_sequence_of_trues(y_ions_matched),
        "matched_b_ions": sum(b_ions_matched),
        "matched_b_ions_pct": sum(b_ions_matched) / len(b_ions_matched),
        "matched_y_ions": sum(y_ions_matched),
        "matched_y_ions_pct": sum(y_ions_matched) / len(y_ions_matched),
        "matched_ions_pct": (sum(b_ions_matched) + sum(y_ions_matched))
        / (len(b_ions_matched) + len(y_ions_matched)),
    }


# TODO: keep this here?
def modification_evidence():
    return


# TODO: move to basic feature generator
def extract_psm_features(psm):
    return {
        "peptide_len": len(psm.peptidoform.sequence),
        "expmass": (psm.precursor_mz * psm.get_precursor_charge())
        - (psm.get_precursor_charge() * 1.007276),  # TODO: replace by constant
        "calcmass": calculate_mass(psm.peptidoform.composition),
    }
