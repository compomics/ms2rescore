import logging
import re
from typing import Dict, Union, Optional
from pathlib import Path
from glob import glob

import numpy as np
import psm_utils.io
from psm_utils import PSMList

from ms2rescore.exceptions import MS2RescoreConfigurationError, MS2RescoreError

logger = logging.getLogger(__name__)


def parse_psms(config: Dict, psm_list: Union[PSMList, None], output_file_root: str) -> PSMList:
    """
    Parse PSMs and prepare for rescoring.

    Parameters
    ----------
    config
        Dictionary containing general ms2rescore configuration (everything under ``ms2rescore``
        top-level key).
    psm_list
        PSMList object containing PSMs. If None, PSMs will be read from ``psm_file``.
    output_file_root
        Path to output file root (without file extension).

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

    # Add filename if all values are none
    # if (psm_list["run"] == None).all():  # noqa: E711
    #     # Map inferred spectrum paths
    spectrum_path_mapping = {
        run: infer_spectrum_path(configured_path=config["spectrum_path"], run_name=run)
        for run in set(psm_list["run"])
    }
    logger.debug(f"Mapped PSM list runs to spectrum file paths: {spectrum_path_mapping}")
    psm_list["run"] = np.vectorize(spectrum_path_mapping.get)(psm_list["run"])
    exit()

    # TODO: Temporary fix until implemented in psm_utils
    # Ensure that spectrum IDs are strings (Pydantic 2.0 does not coerce int to str)
    psm_list["spectrum_id"] = [str(spec_id) for spec_id in psm_list["spectrum_id"]]

    return psm_list


def _read_psms(config, psm_list):
    logger.info("Reading PSMs...")
    if isinstance(psm_list, PSMList):
        return psm_list
    else:
        try:
            return psm_utils.io.read_file(
                config["psm_file"],
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
        raise MS2RescoreError(
            "`psm_id_pattern` could not be matched to all PSM spectrum IDs."
            " Ensure that the regex contains a capturing group?"
        )


def infer_spectrum_path(
    configured_path: Union[str, Path, None],
    run_name: Optional[str] = None,
) -> Union[str, Path]:
    """
    Infer spectrum path from passed path and expected filename (e.g. from PSM file).

    Parameters
    ----------
    configured_path: str, Path, None
        User-defined path to spectrum file or directory containing spectrum file
    run_name : str, optional
        MS run name (stem of spectrum filename), e.g., as expected from PSM file.

    """
    # If no spectrum path configured, use expected run_name in default dir
    if not configured_path:
        if run_name:
            resolved_path = Path(".").joinpath(run_name)
        else:
            raise MS2RescoreConfigurationError(
                "Could not resolve spectrum file name: No spectrum path configured "
                "and no run name in PSM file found."
            )

    else:
        configured_path = Path(configured_path)
        # If passed path is directory, join with run name
        if configured_path.is_dir():
            if run_name:
                resolved_path = configured_path.joinpath(run_name)
            else:
                raise MS2RescoreConfigurationError(
                    "Could not resolve spectrum file name: Spectrum path is directory "
                    "but no run name in PSM file found."
                )

        # If passed path is file, use that, but warn if basename doesn't match expected
        elif configured_path.is_file():
            if run_name and configured_path.stem != Path(run_name).stem:
                logger.warning(
                    "Passed spectrum path (`%s`) does not match run name found in PSM "
                    "file (`%s`). Continuing with passed spectrum path.",
                    configured_path,
                    run_name,
                )
            resolved_path = configured_path
        else:
            raise MS2RescoreConfigurationError(
                "Configured `spectrum_path` must be `None` or a path to an existing file "
                "or directory. If `None` or path to directory, spectrum run information "
                "should be present in the PSM file."
            )

    # Match with file extension if not in resolved_path yet
    if not re.match(".mgf$|.mzml$", resolved_path, flags=re.IGNORECASE):
        for filename in glob(resolved_path + "*"):
            if re.match(r".*(\.mgf$|\.mzml$)", filename, flags=re.IGNORECASE):
                resolved_path = filename
                break
        else:
            raise MS2RescoreConfigurationError(
                f"Resolved spectrum filename '{resolved_path}' does not contain a supported file "
                "extension (mgf or mzml) and could not find any matching existing "
                "files."
            )

    return Path(resolved_path).as_posix()
