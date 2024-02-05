"""Parse MGF files."""

import logging
import re
from itertools import chain
from typing import Dict, Tuple

from psm_utils import PSMList
from pyteomics.mgf import IndexedMGF
from pyteomics.mzml import PreIndexedMzML
from rich.progress import track

from ms2rescore.exceptions import MS2RescoreError
from ms2rescore.utils import infer_spectrum_path

logger = logging.getLogger(__name__)


def get_missing_values(config, psm_list, missing_rt=False, missing_im=False):
    """Get missing RT/IM features from spectrum file."""
    logger.debug("Extracting missing RT/IM values from spectrum file(s).")

    psm_dict = psm_list.get_psm_dict()
    for runs in psm_dict.values():
        for run, psms in track(runs.items(), description="Extracting RT/IM values..."):
            psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
            spectrum_file = infer_spectrum_path(config["spectrum_path"], run)

            if missing_im or missing_rt:
                if spectrum_file.suffix.lower() == ".mzml":
                    _parse_values_from_mzml(
                        psm_list_run, spectrum_file, config, missing_rt, missing_im
                    )
                elif spectrum_file.suffix.lower() == ".mgf":
                    _parse_values_from_mgf(
                        psm_list_run, spectrum_file, config, missing_rt, missing_im
                    )
                else:
                    raise MS2RescoreError(
                        f"Could not parse retention time and/or ion mobility from spectrum file for run {run}. "
                        "Please make sure that the spectrum file is either in mzML or MGF format."
                    )


def _parse_values_from_mgf(
    psm_list_run, spectrum_file, config, missing_rt, missing_im
) -> Tuple[Dict, Dict]:
    """Parse retention time and/or ion mobility from an mzML file."""

    mgf = IndexedMGF(str(spectrum_file))
    spectrum_id_pattern = re.compile(
        config["spectrum_id_pattern"] if config["spectrum_id_pattern"] else r"(.*)"
    )

    try:
        mapper = {
            spectrum_id_pattern.search(spectrum_id).group(1): spectrum_id
            for spectrum_id in mgf._offset_index.mapping["spectrum"].keys()
        }
    except AttributeError:
        raise ParseMGFError(
            "Could not parse spectrum IDs using ´spectrum_id_pattern´. Please make sure that there is a capturing in the pattern."
        )

    spectra = {spectrum_id: mgf.get_by_id(spectrum_id) for spectrum_id in mapper.values()}

    for psm in psm_list_run:
        spectrum = spectra.get(mapper[psm.spectrum_id])
        if spectrum is None:
            raise ParsingError(f"Could not find spectrum with ID {psm.spectrum_id} in MGF file.")

        if missing_rt and "params" in spectrum and "rtinseconds" in spectrum["params"]:
            psm.retention_time = float(spectrum["params"]["rtinseconds"])

        if missing_im and "params" in spectrum and "ion_mobility" in spectrum["params"]:
            psm.ion_mobility = float(spectrum["params"]["ion_mobility"])


def _parse_values_from_mzml(
    psm_list_run, spectrum_file, config, missing_rt, missing_im
) -> Tuple[Dict, Dict]:
    """Parse retention time and/or ion mobility from an mzML file."""

    mzml = PreIndexedMzML(str(spectrum_file))
    spectrum_id_pattern = re.compile(
        config["spectrum_id_pattern"] if config["spectrum_id_pattern"] else r"(.*)"
    )

    try:
        mapper = {
            spectrum_id_pattern.search(spectrum_id).group(1): spectrum_id
            for spectrum_id in mzml._offset_index.mapping["spectrum"].keys()
        }
    except AttributeError as e:
        raise ParseMGFError(
            "Could not parse spectrum IDs using ´spectrum_id_pattern´. Please make sure that there is a capturing in the pattern."
        ) from e

    spectra = {spectrum_id: mzml.get_by_id(spectrum_id) for spectrum_id in mapper.values()}

    for psm in psm_list_run:
        spectrum = spectra.get(mapper[psm.spectrum_id])
        if spectrum is None:
            raise ParsingError(f"Could not find spectrum with ID {psm.spectrum_id} in mzML file.")

        if (
            missing_rt
            and "scanList" in spectrum
            and "scan" in spectrum["scanList"]
            and spectrum["scanList"]["scan"]
        ):
            psm.retention_time = float(spectrum["scanList"]["scan"][0].get("scan start time", 0))

        if missing_im:
            if (
                "precursorList" in spectrum
                and "precursor" in spectrum["precursorList"]
                and spectrum["precursorList"]["precursor"]
            ):
                psm.ion_mobility = float(
                    spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][
                        0
                    ].get("inverse reduced ion mobility", 0)
                )
            elif (
                "scanList" in spectrum
                and "scan" in spectrum["scanList"]
                and spectrum["scanList"]["scan"]
            ):
                psm.ion_mobility = float(
                    spectrum["scanList"]["scan"][0].get("reverse ion mobility", 0)
                )


class ParseMGFError(MS2RescoreError):
    """Error parsing MGF file."""

    pass


class ParsingError(MS2RescoreError):
    """Error parsing retention time from spectrum file."""

    pass
