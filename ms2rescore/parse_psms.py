import logging
import re
from itertools import chain
from typing import Dict, Union

import psm_utils.io
from psm_utils import PSMList

from ms2rescore.exceptions import MS2RescoreConfigurationError

logger = logging.getLogger(__name__)


def parse_psms(config: Dict, psm_list: Union[PSMList, None]) -> PSMList:
    """
    Parse PSMs and prepare for rescoring.

    Parameters
    ----------
    config
        Dictionary containing general ms2rescore configuration (everything under ``ms2rescore``
        top-level key).
    psm_list
        PSMList object containing PSMs. If None, PSMs will be read from ``psm_file``.

    """
    # Read PSMs, find decoys, calculate q-values
    psm_list = _read_psms(config, psm_list)
    _find_decoys(config, psm_list)
    _calculate_qvalues(config, psm_list)

    # Store scoring values for comparison later
    for psm in psm_list:
        psm.provenance_data.update(
            {
                "before_rescoring_score": psm.score,
                "before_rescoring_qvalue": psm.qvalue,
                "before_rescoring_pep": psm.pep,
                "before_rescoring_rank": psm.rank,
            }
        )

    logger.debug("Parsing modifications...")
    psm_list.rename_modifications(config["modification_mapping"])
    psm_list.add_fixed_modifications(config["fixed_modifications"])
    psm_list.apply_fixed_modifications()

    logger.debug("Applying `psm_id_pattern`...")
    if config["psm_id_pattern"]:
        pattern = re.compile(config["psm_id_pattern"])
        new_ids = [_match_psm_ids(old_id, pattern) for old_id in psm_list["spectrum_id"]]
        psm_list["spectrum_id"] = new_ids

    # TODO: Temporary fix until implemented in psm_utils
    # Ensure that spectrum IDs are strings (Pydantic 2.0 does not coerce int to str)
    psm_list["spectrum_id"] = [str(spec_id) for spec_id in psm_list["spectrum_id"]]

    return psm_list


def _read_psms(config, psm_list):
    if isinstance(psm_list, PSMList):
        return psm_list
    else:
        logger.info("Reading PSMs from file...")
        current_file = 1
        total_files = len(config["psm_file"])
        psm_list_list = []
        for psm_file in config["psm_file"]:
            logger.info(
                f"Reading PSMs from PSM file ({current_file}/{total_files}): `{psm_file}`..."
            )
            try:
                id_file_psm_list = psm_utils.io.read_file(
                    psm_file,
                    filetype=config["psm_file_type"],
                    show_progressbar=True,
                    **config["psm_reader_kwargs"],
                )
            except psm_utils.io.PSMUtilsIOException:
                raise MS2RescoreConfigurationError(
                    "Error occurred while reading PSMs. Please check the `psm_file` and "
                    "`psm_file_type` settings. See "
                    "https://ms2rescore.readthedocs.io/en/latest/userguide/input-files/"
                    " for more information."
                )

            psm_list_list.append(id_file_psm_list)
            current_file += 1

        return PSMList(psm_list=list(chain.from_iterable(p.psm_list for p in psm_list_list)))


def _find_decoys(config, psm_list):
    """Find decoys in PSMs, log amount, and raise error if none found."""
    logger.debug("Finding decoys...")
    if config["id_decoy_pattern"]:
        psm_list.find_decoys(config["id_decoy_pattern"])

    n_psms = len(psm_list)
    percent_decoys = sum(psm_list["is_decoy"]) / n_psms * 100
    logger.info(f"Found {n_psms} PSMs, of which {percent_decoys:.2f}% are decoys.")

    if not any(psm_list["is_decoy"]):
        raise MS2RescoreConfigurationError(
            "No decoy PSMs found. Please check if decoys are present in the PSM file and that "
            "the `id_decoy_pattern` option is correct. See "
            "https://ms2rescore.readthedocs.io/en/latest/userguide/configuration/#selecting-decoy-psms"
            " for more information."
        )


def _calculate_qvalues(config, psm_list):
    """Calculate q-values for PSMs if not present."""
    # Calculate q-values if not present
    if None in psm_list["qvalue"]:
        logger.debug("Recalculating q-values...")
        psm_list.calculate_qvalues(reverse=not config["lower_score_is_better"])


def _match_psm_ids(old_id, regex_pattern):
    """Match PSM IDs to regex pattern or raise Exception if no match present."""
    match = re.search(regex_pattern, str(old_id))
    try:
        return match[1]
    except (TypeError, IndexError):
        raise MS2RescoreConfigurationError(
            "`psm_id_pattern` could not be matched to all PSM spectrum IDs."
            " Ensure that the regex contains a capturing group?"
        )
