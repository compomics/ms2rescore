from ms2rescore.feature_generators.base import FeatureGeneratorBase
from typing import List
from psm_utils import PSM, PSMList

import numpy as np
from tqdm import tqdm

from .evidence import PeptideEvidence

from pyteomics import mgf
import polars as pl
import os
from glob import glob
import numpy as np
from rustyms import RawSpectrum, LinearPeptide, FragmentationModel

from .utils import (
    get_annotated_spectrum,
    annot_peaks_to_fragments,
    fragments_to_polars,
    parse_to_iontype_dict,
    ion_dict_to_matrix,
    matrix_to_ion_dict,
    calculate_ppm
)

class SpectrumVector():
    def __init__(
            self,
            ion_types,
            neutral_losses
    ):
        self.annotated = False
        self.ion_types = ion_types
        self.neutral_losses = neutral_losses

    def parse(self, psm, annot_spec, theo_frags):

        self.peptidoform = psm.peptidoform
        self.spectrum_id = psm.spectrum_id

        # Convert the list of fragments to a polars object
        annot_frags, mz_array, intensity_array = annot_peaks_to_fragments(annot_spec)
        spec_polars_theo = fragments_to_polars(
            fragment_list=theo_frags,
            ion_types=self.ion_types,
            neutral_losses=self.neutral_losses
        )
        spec_polars_annot = fragments_to_polars(
            fragment_list=annot_frags,
            ion_types=self.ion_types,
            neutral_losses=self.neutral_losses,
            mz_array=mz_array,
            intensity_array=intensity_array
        )
        
        # Parse the polars spectrum to a dictionary
        result_dict_theo_mz = parse_to_iontype_dict(
            pl_df=spec_polars_theo,
            len_pep=len(psm.peptidoform),
            ion_types=self.ion_types,
            value="mz"
        )
        result_dict_annot_mz = parse_to_iontype_dict(
            pl_df=spec_polars_annot,
            len_pep=len(psm.peptidoform),
            ion_types=self.ion_types,
            value="mz"
        )
        result_dict_annot_intensity = parse_to_iontype_dict(
            pl_df=spec_polars_annot,
            len_pep=len(psm.peptidoform),
            ion_types=self.ion_types,
            value="intensity"
        )

        self.theoretical_mz = result_dict_theo_mz
        self.experimental_mz = self._complete_ion_dict(result_dict_annot_mz)
        self.experimental_intensity = self._complete_ion_dict(result_dict_annot_intensity)
        self.annotated=True

    def _complete_ion_dict(self, ion_dict):
        for nl in self.neutral_losses:
            if nl not in ion_dict.keys():
                ion_dict[nl] = {i: np.zeros(self.n) for i in self.ion_types}
            else:
                for ion_type in self.ion_types:
                    if ion_type not in ion_dict[nl].keys():
                        ion_dict[nl][ion_type] = np.zeros(self.n)
        return ion_dict

    def load(
            self,
            theoretical_mz,
            experimental_mz,
            experimental_intensity,
            tic,
            peptidoform,
            spectrum_id,
            precursor_mz,
            precursor_ppm
        ):
        """
        Load a parsed spectrum.
        """
        self.annotated = True
        self.theoretical_mz = theoretical_mz
        self.experimental_mz = experimental_mz
        self.experimental_intensity = experimental_intensity
        self.peptidoform = peptidoform
        self.spectrum_id = spectrum_id

    def delete_shared_peaks(self, method):
        self.experimental_mz, self.experimental_intensity = method(self.experimental_mz, self.experimental_intensity)

    @property
    def ppm_diff(self):
        matrix_theoretical = ion_dict_to_matrix(
            self.theoretical_mz[""],
            ion_types=self.ion_types,
            n=len(self.theoretical_mz[""]["y1"])
        )
        matrix_theoretical[matrix_theoretical==0] = 1
        matrix_experimental = ion_dict_to_matrix(
            self.experimental_mz[""],
            ion_types=self.ion_types,
            n=len(self.theoretical_mz[""]["y1"])
        )

        ppm_diff = np.where(
            matrix_experimental!=0, 
            calculate_ppm(
                matrix_theoretical,
                matrix_experimental,
            ),
            np.nan
        )
        return matrix_to_ion_dict(ppm_diff, ion_types=self.ion_types)

    @property
    def n(self):
        return len(self.peptidoform)-1
    
    def to_ion_matrix(self):
        return np.concatenate(
            [
                ion_dict_to_matrix(
                    ion_dict=m,
                    ion_types=self.ion_types,
                    n=self.n
                ) 
                for m in self.experimental_intensity.values()
            ],
            axis=0
        )

class MissingFragmentationFeatures(FeatureGeneratorBase):
    """MS2Rescore type feature generator for peak-type features."""

    def __init__(self, *args, **kwargs):
        """Initialize feature generator class."""
        super().__init__(*args, **kwargs)

    @property
    def feature_names(self) -> List[str]:
        """Names of features added to rescoring_features dict."""
        return [
            "missing_frag_count",
            "missing_frag_pct",
            "missing_frag_sites",
            "missing_frag_longest_sequence",
            "missing_frag_from_n",
            "missing_frag_from_c",
        ]

    @property
    def name(self) -> str:
        return "missing_frag"
    @property
    def input_type(self) -> str:
        return "peptide_evidence"

    def add_features(self, psm_list: PSMList, peptide_evidence_list: List[PeptideEvidence]) -> None:
        for psm, pe in tqdm(zip(psm_list, peptide_evidence_list)):

            psm.rescoring_features.update(
                self.calculate_features(pe.evidence)
            )

    def calculate_features(self, evidence):
        # Ensure it's a boolean array
        arr = np.asarray(evidence, dtype=bool)
        
        # Number of consecutive False entries starting from the start
        consecutive_false_start = np.argmax(arr) if not arr.all() else 0

        # Number of consecutive False entries starting from the end
        consecutive_false_end = np.argmax(arr[::-1]) if not arr.all() else 0

        # Number of False entries in the array
        num_false_entries = np.size(arr) - np.sum(arr)

        # Percentage of False entries in the array
        percentage_false = num_false_entries / len(arr) *100

        # Finding the longest sequence of consecutive False values
        # Identify the regions where values change (False <-> True)
        # Create an array to mark transitions
        if not arr.all():
            # Extend the array with True at both ends to handle edge cases
            padded_arr = np.r_[True, arr, True]
            diff_arr = np.diff(padded_arr.astype(int)) # Transitions
            false_starts = np.flatnonzero(diff_arr == -1)
            false_ends = np.flatnonzero(diff_arr == 1)
            longest_false_seq = (false_ends - false_starts).max() if len(false_starts) > 0 else 0
        else:
            longest_false_seq = 0

        # Number of islands of False entries (bounded by True or at boundaries)
        # This means we count all transitions from True to False.
        # Adjusting the island definition to include boundary cases:
        if not arr.all():
            # A False island starts at -1 or at the beginning of the array
            # and ends at 1 or the end of the array
            num_false_islands = len(false_starts)
        else:
            num_false_islands = 0

        return {
            'missing_frag_from_n': consecutive_false_start,
            'missing_frag_from_c': consecutive_false_end,
            'missing_frag_sites': num_false_entries,
            'missing_frag_pct': percentage_false,
            'missing_frag_longest_sequence': longest_false_seq,
            'missing_frag_count': num_false_islands
        }