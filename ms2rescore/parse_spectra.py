"""Parse MGF files."""

import logging
import re
from itertools import chain

from ms2rescore_rs import get_precursor_info
from psm_utils import PSMList

from ms2rescore.exceptions import MS2RescoreError
from ms2rescore.utils import infer_spectrum_path

logger = logging.getLogger(__name__)


def get_missing_values(
    psm_list: PSMList, config: dict, rt_required: bool = False, im_required: bool = False
):
    """Get missing RT/IM features from spectrum file."""
    psm_dict = psm_list.get_psm_dict()
    for runs in psm_dict.values():
        for run, psms in runs.items():
            psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
            spectrum_file = infer_spectrum_path(config["spectrum_path"], run)
            logger.debug("Reading spectrum file: '%s'", spectrum_file)
            precursors = get_precursor_info(str(spectrum_file))

            if config["spectrum_id_pattern"]:
                spectrum_id_pattern = re.compile(config["spectrum_id_pattern"])
                precursors = {
                    spectrum_id_pattern.search(spectrum_id).group(1): precursor
                    for spectrum_id, precursor in precursors.items()
                }

            for psm in psm_list_run:
                try:
                    if rt_required:
                        psm.retention_time = precursors[psm.spectrum_id].rt
                    if im_required:
                        psm.ion_mobility = precursors[psm.spectrum_id].im
                    if not psm.precursor_mz:
                        psm.precursor_mz = precursors[psm.spectrum_id].mz
                except KeyError as e:
                    raise SpectrumParsingError(
                        f"Could not extract missing RT/IM values from spectrum file for run {run}."
                    ) from e


class SpectrumParsingError(MS2RescoreError):
    """Error parsing retention time from spectrum file."""

    pass
