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
from pyteomics import mass, mzml, mgf
from rustyms import RawSpectrum, LinearPeptide, FragmentationModel, MassMode

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from ms2rescore.utils import infer_spectrum_path
from ms2rescore.exceptions import ParseSpectrumError

logger = logging.getLogger(__name__)

FRAGMENTATION_MODELS = {
    "cidhcd": FragmentationModel.CidHcd,
    "etcid": FragmentationModel.Etcid,
    "etd": FragmentationModel.Etd,
    "ethcd": FragmentationModel.Ethcd,
    "all": FragmentationModel.All,
}
MASS_MODES = {
    "average": MassMode.Average,
    "monoisotopic": MassMode.Monoisotopic,
}


class MS2FeatureGenerator(FeatureGeneratorBase):
    """DeepLC retention time-based feature generator."""

    def __init__(
        self,
        *args,
        spectrum_path: Optional[str] = None,
        spectrum_id_pattern: str = "(.*)",
        fragmentationmodel: str = "All",
        massmode: str = "Monoisotopic",
        **kwargs,
    ) -> None:
        """
        Generate MS2-based features for rescoring.

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
        self.fragmentation_model = FRAGMENTATION_MODELS[fragmentationmodel.lower()]
        self.mass_mode = MASS_MODES[massmode.lower()]

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
        """Add MS2-derived features to PSMs."""

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

                annotated_spectra = self._annotate_spectra(psm_list_run, spectrum_filename)
                self._calculate_features(psm_list_run, annotated_spectra)
                current_run += 1

    def _calculate_features(self, psm_list: PSMList, annotated_spectra: List) -> None:
        """Calculate features from annotated spectra and add them to the PSMs."""
        for psm, spectrum in zip(psm_list, annotated_spectra):
            psm.features.update(self.extract_spectrum_features(spectrum, psm.peptidoform))

    def _annotate_spectra(self, psm_list: PSMList, spectrum_file: str) -> List:
        """Retrieve annotated spectra for all psms."""

        annotated_spectra = []
        spectrum_reader = self._read_spectrum_file(spectrum_file)

        spectrum_id_pattern = re.compile(
            self.spectrum_id_pattern if self.spectrum_id_pattern else r"(.*)"
        )
        logger.info(spectrum_id_pattern)
        logger.info(spectrum_reader._offset_index.mapping["spectrum"].keys())

        try:
            mapper = {
                spectrum_id_pattern.search(spectrum_id).group(1): spectrum_id
                for spectrum_id in spectrum_reader._offset_index.mapping["spectrum"].keys()
            }
        except AttributeError:
            raise ParseSpectrumError(
                "Could not parse spectrum IDs using ´spectrum_id_pattern´. Please make sure that there is a capturing in the pattern."
            )

        for psm in psm_list:
            pyteomics_spectrum = spectrum_reader.get_by_id(mapper[psm.spectrum_id])
            annotated_spectrum = self._annotate_spectrum(psm, pyteomics_spectrum)
            annotated_spectra.append(annotated_spectrum)

        return annotated_spectra

    @staticmethod
    def _read_spectrum_file(spectrum_filepath: str) -> Union[mzml.PreIndexedMzML, mgf.IndexedMGF]:

        if spectrum_filepath.suffix == ".mzML":
            return mzml.PreIndexedMzML(str(spectrum_filepath))
        elif spectrum_filepath.suffix == ".mgf":
            return mgf.IndexedMGF(str(spectrum_filepath))

    @staticmethod
    def _longest_ion_sequence(lst):
        max_sequence = 0
        current_sequence = 0

        for value in lst:
            current_sequence = current_sequence + 1 if value else 0
            max_sequence = max(max_sequence, current_sequence)

        return max_sequence

    def extract_spectrum_features(self, annotated_spectrum, peptidoform):
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
            "longest_b_ion_sequence": self._longest_ion_sequence(b_ions_matched),
            "longest_y_ion_sequence": self._longest_ion_sequence(y_ions_matched),
            "matched_b_ions": sum(b_ions_matched),
            "matched_b_ions_pct": sum(b_ions_matched) / len(b_ions_matched),
            "matched_y_ions": sum(y_ions_matched),
            "matched_y_ions_pct": sum(y_ions_matched) / len(y_ions_matched),
            "matched_ions_pct": (sum(b_ions_matched) + sum(y_ions_matched))
            / (len(b_ions_matched) + len(y_ions_matched)),
        }

    def _annotate_spectrum(self, psm, spectrum):

        spectrum = RawSpectrum(
            title=psm.spectrum_id,
            num_scans=1,
            rt=psm.retention_time,
            precursor_charge=psm.get_precursor_charge(),
            precursor_mass=mass.calculate_mass(psm.peptidoform.compsoition),
            mz_array=spectrum["m/z array"],
            intensity_array=spectrum["intensity array"],
        )

        annotated_spectrum = spectrum.annotate(
            peptide=LinearPeptide(psm.peptidoform.proforma),
            model=self.fragmentationmodel,
            mode=self.massmode,
        )

        return annotated_spectrum


# TODO: keep this here?
def modification_evidence():
    return


# TODO: move to basic feature generator
def extract_psm_features(psm):

    expmass = (psm.precursor_mz * psm.get_precursor_charge()) - (
        psm.get_precursor_charge() * 1.007276
    )  # TODO: replace by constant
    calcmass = mass.calculate_mass(psm.peptidoform.composition)
    delta_mass = expmass - calcmass

    return {
        "peptide_len": len(psm.peptidoform.sequence),
        "expmass": expmass,
        "calcmass": calcmass,
        "delta_mass": delta_mass,
    }
