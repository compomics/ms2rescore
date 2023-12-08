"""Parse MGF files."""

import logging
import re
from itertools import chain
from typing import Dict, Tuple

from psm_utils import PSMList
from pyteomics.mgf import MGF
from pyteomics.mzml import MzML
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

            rt_dict = None
            im_dict = None
            rt_dict_mzml = None
            im_dict_mzml = None
            rt_dict_mgf = None
            im_dict_mgf = None

            rt_dict, im_dict, missing_rt, missing_im = _parse_values_spectrum_id(
                config, psm_list, missing_rt, missing_im
            )

            if missing_im or missing_rt:
                if spectrum_file.suffix.lower() == ".mzml":
                    rt_dict_mzml, im_dict_mzml = _parse_values_from_mzml(
                        spectrum_file, config, run, missing_rt, missing_im
                    )
                elif spectrum_file.suffix.lower() == ".mgf":
                    rt_dict_mgf, im_dict_mgf = _parse_values_from_mgf(
                        spectrum_file, config, run, missing_rt, missing_im
                    )
            rt_dict = rt_dict or rt_dict_mzml or rt_dict_mgf
            im_dict = im_dict or im_dict_mzml or im_dict_mgf

            for value_dict, value in zip([rt_dict, im_dict], ["retention_time", "ion_mobility"]):
                if value_dict:
                    try:
                        psm_list_run[value] = [value_dict[psm.spectrum_id] for psm in psm_list_run]
                    except KeyError:
                        raise ParsingError(
                            f"Could not parse {value} values from spectrum file for run {run}."
                        )


def _parse_values_from_mgf(
    spectrum_file, config, run, missing_rt, missing_im
) -> Tuple[Dict, Dict]:
    """
    Parse retention time and/or ion mobility from an MGF file.

    Notes
    -----
    - Extracting values (e.g., ion mobility) according to the Matrix documentation:
      http://www.matrixscience.com/help/data_file_help.html

    """
    rt_dict = {}
    im_dict = {}

    spectrum_id_pattern = re.compile(
        config["spectrum_id_pattern"] if config["spectrum_id_pattern"] else r"(.*)"
    )

    for spectrum in MGF(str(spectrum_file)):
        matched_id = spectrum_id_pattern.match(spectrum["params"]["title"]).group()
        if missing_rt:
            try:
                rt_dict[matched_id] = float(spectrum["params"]["rtinseconds"])
            except KeyError:
                raise ParsingError(
                    "Could not parse retention time (`rtinseconds`) from spectrum file for "
                    f"run {run}. Please make sure that the retention time key is present in the "
                    "spectrum file or disable the relevant feature generator."
                )
        if missing_im:
            try:
                im_dict[matched_id] = float(spectrum["params"]["ion_mobility"])
            except KeyError:
                raise ParsingError(
                    "Could not parse ion mobility (`ion_mobility`) from spectrum file "
                    f"for run {run}. Please make sure that the ion mobility key is present in the "
                    "spectrum file or disable the relevant feature generator."
                )

    return rt_dict, im_dict


def _parse_values_from_mzml(
    spectrum_file, config, run, missing_rt, missing_im
) -> Tuple[Dict, Dict]:
    """Parse retention time and/or ion mobility from an mzML file."""
    rt_dict = {}
    im_dict = {}

    spectrum_id_pattern = re.compile(
        config["spectrum_id_pattern"] if config["spectrum_id_pattern"] else r"(.*)"
    )

    for spectrum in MzML(str(spectrum_file)):
        matched_id = spectrum_id_pattern.match(spectrum["id"]).group()
        if missing_rt:
            try:
                rt_dict[matched_id] = float(spectrum["scanList"]["scan"][0]["scan start time"])
            except KeyError:
                raise ParsingError(
                    "Could not parse retention time (`scan start time`) from spectrum file for "
                    f"run {run}. Please make sure that the retention time key is present in the "
                    "spectrum file or disable the relevant feature generator."
                )
        if missing_im:
            try:
                im_dict[matched_id] = float(
                    spectrum["scanList"]["scan"][0]["reverse ion mobility"]
                )
            except KeyError:
                raise ParsingError(
                    "Could not parse ion mobility (`reverse ion mobility`) from spectrum file "
                    f"for run {run}. Please make sure that the ion mobility key is present in the "
                    "spectrum file or disable the relevant feature generator."
                )

    return rt_dict, im_dict


def _parse_values_spectrum_id(config, psm_list, missing_rt, missing_im):
    """Parse the missing values from the spectrum_id."""
    rt_pattern = re.compile(r'RetentionTime:"([\d\.]+)"')
    im_pattern = re.compile(r'IonMobility:"([\d\.]+)"')

    if missing_rt and "retention_time" in psm_list["spectrum_id"][0]:
        rt_dict = {
            psm.spectrum_id: float(rt_pattern.search(psm.spectrum_id).group(1)) for psm in psm_list
        }
        missing_rt = False
    else:
        rt_dict = {}
    if missing_im and "IonMobility" in psm_list["spectrum_id"][0]:
        im_dict = {
            psm.spectrum_id: float(im_pattern.search(psm.spectrum_id).group(1)) for psm in psm_list
        }
        missing_im = False
    else:
        im_dict = {}

    return rt_dict, im_dict, missing_rt, missing_im


class ParseMGFError(MS2RescoreError):
    """Error parsing MGF file."""

    pass


class ParsingError(MS2RescoreError):
    """Error parsing retention time from spectrum file."""

    pass
