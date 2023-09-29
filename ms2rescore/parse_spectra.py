"""Parse MGF files."""

import logging
import re
from typing import Union, Tuple, Dict

from rich.progress import track
from pyteomics.mgf import MGF
from pyteomics.mzml import MzML
from itertools import chain

from ms2rescore.exceptions import MS2RescoreError
from ms2rescore.utils import infer_spectrum_path
from psm_utils import PSMList

logger = logging.getLogger(__name__)


class ParseMGFError(MS2RescoreError):
    """Error parsing MGF file."""

    pass


class ParsingError(MS2RescoreError):
    """Error parsing retention time from spectrum file."""

    pass


def get_missing_values(config, psm_list, missing_rt_values=False, missing_im_values=False):
    """Get missing features from spectrum file."""

    rt_dict = {}
    im_dict = {}

    psm_dict = psm_list.get_psm_dict()
    for runs in psm_dict.values():
        for run, psms in track(runs.items()):
            psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
            spectrum_file = infer_spectrum_path(config["spectrum_path"], run)

            if spectrum_file.suffix.lower() == ".mzml":
                rt_dict, im_dict = _parse_values_from_mzml(
                    spectrum_file, config, run, missing_rt_values, missing_im_values
                )
            elif spectrum_file.suffix.lower() == ".mgf":
                rt_dict, im_dict = _parse_values_from_mgf(
                    spectrum_file, config, run, missing_rt_values, missing_im_values
                )

            for value_dict, value in zip([rt_dict, im_dict], ["retention_time", "ion_mobility"]):
                if value_dict:
                    try:
                        psm_list_run[value] = [value_dict[psm.spectrum_id] for psm in psm_list_run]
                    except KeyError:
                        raise ParsingError(
                            f"Could not parse {value} values from spectrum file for run {run}."
                        )


def _parse_values_from_mgf(
    spectrum_file, config, run, missing_rt_values, missing_im_values
) -> Tuple[Dict, Dict]:
    """Parse retention time and/or ion mobility from MGF file."""
    rt_dict = None
    im_dict = None

    for spectrum in MGF(str(spectrum_file)):
        if missing_rt_values:
            rt_dict = {}
            try:
                rt_dict[
                    re.match(config["spectrum_id_pattern"], spectrum["params"]["title"]).group()
                ] = float(spectrum["params"]["rtinseconds"])
            except KeyError:
                raise ParsingError(
                    f"Could not parse retention time key `rtinsecondes` from spectrum file for run {run}."
                    "Please make sure that the retention time key is present in the spectrum file."
                )
        if missing_im_values:
            im_dict = {}
            try:
                im_dict[
                    re.match(config["spectrum_id_pattern"], spectrum["params"]["title"]).group()
                ] = float(spectrum["params"]["ionmobility"])
            except KeyError:
                raise ParsingError(
                    f"Could not parse ion mobility key `ionmobility` from spectrum file for run {run}."
                    "Please make sure that the ion mobility key is present in the spectrum file."
                )

    return rt_dict, im_dict


def _parse_values_from_mzml(
    spectrum_file, config, run, missing_rt_values, missing_im_values
) -> Tuple[Dict, Dict]:
    """Parse retention time and/or ion mobility from MGF file."""
    rt_dict = None
    im_dict = None

    for spectrum in MGF(str(spectrum_file)):
        if missing_rt_values:
            rt_dict = {}
            try:
                rt_dict[re.match(config["spectrum_id_pattern"], spectrum["id"]).group()] = float(
                    spectrum["scanList"]["scan"][0][
                        "scan start time"
                    ]  # is rt in minutes by default?
                )
            except KeyError:
                raise ParsingError(
                    f"Could not parse retention time key `scan start time` from spectrum file for run {run}."
                    "Please make sure that the retention time key is present in the spectrum file."
                )
        if missing_im_values:
            im_dict = {}
            try:
                im_dict[re.match(config["spectrum_id_pattern"], spectrum["id"]).group()] = float(
                    spectrum["scanList"]["scan"][0]["reverse ion mobility"]
                )
            except KeyError:
                raise ParsingError(
                    f"Could not parse ion mobility key `reverse ion mobility` from spectrum file for run {run}."
                    "Please make sure that the ion mobility key is present in the spectrum file."
                )

    return rt_dict, im_dict
