"""
MS2-based feature generator.

"""

import logging
import re
from collections import defaultdict
from itertools import chain
from typing import List, Optional, Union

import numpy as np
import pyopenms as oms
from psm_utils import PSM, Peptidoform, PSMList
from psm_utils.peptidoform import PeptidoformException
from pyteomics.proforma import ProFormaError
from pyteomics import mass, mgf, mzml
from rustyms import FragmentationModel, LinearPeptide, MassMode, RawSpectrum
from pathlib import Path

from ms2rescore.exceptions import ParseSpectrumError
from ms2rescore.feature_generators.base import FeatureGeneratorBase
from ms2rescore.utils import infer_spectrum_path

logger = logging.getLogger(__name__)

FRAGMENTATION_MODELS = {
    "cidhcd": FragmentationModel.CidHcd,
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
        fragmentation_model: str = "All",
        mass_mode: str = "Monoisotopic",
        processes: int = 1,
        calculate_hyperscore: bool = True,  # Allow optional ?
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
        self.fragmentation_model = FRAGMENTATION_MODELS[fragmentation_model.lower()]
        self.mass_mode = MASS_MODES[mass_mode.lower()]
        self.processes = processes
        self.calculate_hyperscore = calculate_hyperscore


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
            "hyperscore",
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

                self._calculate_features(psm_list_run, spectrum_filename)
                current_run += 1

    def _calculate_features(self, psm_list: PSMList, spectrum_file: str) -> List:
        """Retrieve annotated spectra for all psms."""

        spectrum_reader = self._read_spectrum_file(spectrum_file)

        spectrum_id_pattern = re.compile(
            self.spectrum_id_pattern if self.spectrum_id_pattern else r"(.*)"
        )

        try:
            mapper = {
                spectrum_id_pattern.search(spectrum_id).group(1): spectrum_id
                for spectrum_id in spectrum_reader._offset_index.mapping["spectrum"].keys()
            }
        except AttributeError:
            raise ParseSpectrumError(
                "Could not parse spectrum IDs using ´spectrum_id_pattern´. Please make sure that there is a capturing in the pattern."
            )
        annotated_spectra = [
            self._annotate_spectrum(psm, spectrum_reader.get_by_id(mapper[psm.spectrum_id]))
            for psm in psm_list
        ]
        for psm, annotated_spectrum in zip(psm_list, annotated_spectra):
            psm.rescoring_features.update(
                self._calculate_spectrum_features(psm, annotated_spectrum)
            )

            if self.calculate_hyperscore:
                # Filters out peaks which are unnannotated (can be specified, but keep at b-y ions of any charge ?)
                mz_list, intensity_list = _annotated_spectrum_to_mzint(
                    annotated_spectrum=annotated_spectrum
                )
                hyperscore = calculate_hyperscore(
                    psm=psm, observed_mz=mz_list, observed_intensity=intensity_list
                )
                if hyperscore is None:
                    continue
                psm.rescoring_features.update(
                    {
                        "hyperscore": hyperscore
                    }
                )


    @staticmethod
    def _read_spectrum_file(spectrum_filepath: Path) -> Union[mzml.PreIndexedMzML, mgf.IndexedMGF]:

        if spectrum_filepath.suffix.lower() == ".mzml":
            return mzml.PreIndexedMzML(str(spectrum_filepath))
        elif spectrum_filepath.suffix.lower() == ".mgf":
            return mgf.IndexedMGF(str(spectrum_filepath))

    @staticmethod
    def _longest_ion_sequence(lst):
        max_sequence = 0
        current_sequence = 0

        for value in lst:
            current_sequence = current_sequence + 1 if value else 0
            max_sequence = max(max_sequence, current_sequence)

        return max_sequence

    def _calculate_spectrum_features(self, psm, annotated_spectrum):

        if not annotated_spectrum:
            return {}

        features = defaultdict(list)
        b_ions_matched = [False] * (len(psm.peptidoform.sequence))
        y_ions_matched = [False] * (len(psm.peptidoform.sequence))

        pseudo_count = 0.00001
        ion_fragment_regex = re.compile(r"\d+")

        for peak in annotated_spectrum:
            features["total_intensity"].append(peak.intensity)

            if peak.annotation:
                features["matched_intensity"].append(peak.intensity)
                for matched_ion in peak.annotation:
                    if "y" in matched_ion.ion:
                        features["y_ion_matched"].append(peak.intensity)
                        y_ions_matched[int(ion_fragment_regex.search(matched_ion.ion).group())] = (
                            True
                        )
                    elif "b" in matched_ion.ion:
                        features["b_ion_matched"].append(peak.intensity)
                        b_ions_matched[int(ion_fragment_regex.search(matched_ion.ion).group())] = (
                            True
                        )

        total_intensity_sum = np.sum(features["total_intensity"])
        matched_intensity_sum = np.sum(features["matched_intensity"])
        b_ion_matched_sum = np.sum(features["b_ion_matched"])
        y_ion_matched_sum = np.sum(features["y_ion_matched"])

        return {
            "ln_explained_intensity": np.log(matched_intensity_sum + pseudo_count),
            "ln_total_intensity": np.log(total_intensity_sum + pseudo_count),
            "ln_explained_intensity_ratio": np.log(
                matched_intensity_sum / total_intensity_sum + pseudo_count
            ),
            "ln_explained_b_ion_ratio": np.log(
                b_ion_matched_sum / matched_intensity_sum + pseudo_count
            ),
            "ln_explained_y_ion_ratio": np.log(
                y_ion_matched_sum / matched_intensity_sum + pseudo_count
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

    def _annotate_spectrum(self, psm, pyteomics_spectrum):

        spectrum = RawSpectrum(
            title=psm.spectrum_id,
            num_scans=1,
            rt=psm.retention_time,
            precursor_charge=psm.get_precursor_charge(),
            precursor_mass=mass.calculate_mass(psm.peptidoform.composition),
            mz_array=pyteomics_spectrum["m/z array"],
            intensity_array=pyteomics_spectrum["intensity array"],
        )
        try:
            linear_peptide = LinearPeptide(psm.peptidoform.proforma.split("/")[0])
            annotated_spectrum = spectrum.annotate(
                peptide=linear_peptide,
                model=self.fragmentation_model,
                mode=self.mass_mode,
            )
        except:  # noqa E722
            return []

        return annotated_spectrum.spectrum


############################
### HYPERSCORE FUNCTIONS ###
############################
def _annotated_spectrum_to_mzint(annotated_spectrum, ion_types=["b", "y"]):

    mz_list = []
    intensity_list = []

    for peak in annotated_spectrum:

        annotations = peak.annotation
        for fragment in annotations:
            ion = fragment.ion
            if ion[0] not in ion_types:
                continue

            mz_list.append(peak.experimental_mz)
            intensity_list.append(peak.intensity)
            break

    return mz_list, intensity_list


def _peptidoform_to_oms(peptidoform: Peptidoform) -> tuple[oms.AASequence, Optional[int]]:
    """
    Parse a peptide sequence in proforma format to pyOpenMS compatible format.

    Only supports UNIMOD format.

    Parameter
    ---------
    peptide: str
        Peptide string in proforma format

    Returns
    -------
    AASequence (pyOpenMS):
        A peptide sequence in pyOpenMS format
    """
    peptide = peptidoform.proforma

    # Reformat unimod modifications
    pattern_unimod = r"\[UNIMOD:(\d+)\]"

    def replace_unimod(match):
        return f"(UniMod:{match.group(1)})"

    peptide_oms_str = re.sub(
        pattern=pattern_unimod, repl=replace_unimod, string=peptide
    )

    # Parse N-terminal modifications
    if ")-" in peptide_oms_str:
        peptide_oms_list = peptide_oms_str.split(")-")
        nterm_modification, peptide_oms_str = peptide_oms_list[-2], peptide_oms_list[-1]
        nterm_modification += ")"
        peptide_oms_str = "." + nterm_modification + peptide_oms_str + "."
    elif "]-" in peptide_oms_str:
        peptide_oms_list = peptide_oms_str.split("]-")
        nterm_modification, peptide_oms_str = peptide_oms_list[-2], peptide_oms_list[-1]
        nterm_modification += "]"
        peptide_oms_str = "." + nterm_modification + peptide_oms_str + "."

    # Split the charge from the peptide string
    if "/" in peptide_oms_str:
        peptide_oms_str, _ = peptide_oms_str.split("/")

    try:
        peptide_oms = oms.AASequence.fromString(peptide_oms_str)
        return peptide_oms

    except:
        return


# def _peptidoform_to_oms(peptidoform: Peptidoform) -> tuple[oms.AASequence, Optional[int]]:
#     """
#     Parse a peptidoform object to pyOpenMS compatible format.

#     Parameter
#     ---------
#     Peptidoform: Peptidoform
#         Peptide string in Peptidoform format

#     Returns
#     -------
#     AASequence (pyOpenMS):
#         A peptide sequence in pyOpenMS format
#     int:
#         charge of the peptide
#     """

#     n_term = peptidoform.properties["n_term"]
#     peptide_oms_str = f"[{sum([mod.mass for mod in n_term])}]" if n_term else ""

#     for aa, mods in peptidoform.parsed_sequence:
#         peptide_oms_str += aa
#         if isinstance(mods, list):
#             peptide_oms_str += f"[{sum([mod.mass for mod in mods])}]"

#     c_term = peptidoform.properties["c_term"]
#     peptide_oms_str += f"[{sum([mod.mass for mod in c_term])}]" if c_term else ""

#     peptide_oms = oms.AASequence.fromString(peptide_oms_str)

#     return peptide_oms


def _peptidoform_to_theoretical_spectrum(peptidoform: str) -> oms.MSSpectrum:
    """
    Create a theoretical spectrum from a peptide sequence.

    Parameter
    ---------
    peptide: str
        Peptide sequence in proforma format
    engine: str
        The engine to use to create theoretical spectrum.
        Can only be 'pyopenms' or 'spectrum-utils' (default)

    Return
    ------
    MSSpectrum
        Spectrum object in pyOpenMS format
    """
    # Reformat peptide sequence in pyOpenMS format
    peptide_oms = _peptidoform_to_oms(peptidoform=peptidoform)
    if peptide_oms is None:
        return

    # Initialize the required objects to create the spectrum
    spectrum = oms.MSSpectrum()
    tsg = oms.TheoreticalSpectrumGenerator()
    p = oms.Param()

    p.setValue("add_b_ions", "true")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(param=p)

    # Create the theoretical spectrum
    tsg.getSpectrum(spec=spectrum, peptide=peptide_oms, min_charge=1, max_charge=2)
    return spectrum


def calculate_hyperscore(
    psm: PSM,
    observed_mz: List[float],
    observed_intensity: List[float],
    fragment_tol_mass=20,
    fragment_tol_mode="ppm",
):
    """
    Calculate the hyperscore as defined in the X!Tandem search engine.

    It is a metric of how good two spectra match with each other (matching peaks).

    Parameters
    ----------
    psm: psm_utils.PSM
        The PSM used to extract 'spectrum_id' (for MGF spectrum extraction)
        and 'Peptidoform' (the peptide sequence)
    observed_mz: List[float]
        List of observed mz values with matching order as observed intensity
    observed_intensity: List[float]
        List of observed intensity values
    fragment_tol_mass: int
        The allowed tolerance to match peaks
    fragment_tol_mode: str
        'ppm' for parts-per-million mode. 'Da' for fragment_tol_mass in Dalton.
    Return
    ------
    int
        The hyperscore
    """
    if fragment_tol_mode == "ppm":
        fragment_mass_tolerance_unit_ppm = True
    elif fragment_tol_mode == "Da":
        fragment_mass_tolerance_unit_ppm = False
    else:
        raise Exception(
            "fragment_tol_mode can only take 'Da' or 'ppm'. {} was provided.".format(
                fragment_tol_mode
            )
        )
    if len(observed_intensity) == 0:
        logging.warning(f"PSM ({psm.spectrum_id}) has no annotated peaks.")
        return

    theoretical_spectrum = _peptidoform_to_theoretical_spectrum(peptidoform=psm.peptidoform)
    # This is mainly the cause of the modification not being allowed according to pyopenms
    # pyOpenMS sets stringent criteria on modification being placed on annotated amino acids
    # according to the unimod database
    if theoretical_spectrum is None:
        logging.warning(f'Peptidoform has either unsupported modifications or is being placed on non-allowed residue: {psm.peptidoform.proforma}')
        return

    observed_spectrum_oms = oms.MSSpectrum()
    observed_spectrum_oms.set_peaks([observed_mz, observed_intensity])
    hyperscore = oms.HyperScore()
    result = hyperscore.compute(
        fragment_mass_tolerance=fragment_tol_mass,
        fragment_mass_tolerance_unit_ppm=fragment_mass_tolerance_unit_ppm,
        exp_spectrum=observed_spectrum_oms,
        theo_spectrum=theoretical_spectrum,
    )
    return result


###########################
### FRAG SITE FUNCTIONS ###
###########################
