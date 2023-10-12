import logging
import os
import re
from glob import glob
from pathlib import Path
from typing import Optional, Union

from ms2rescore.exceptions import MS2RescoreConfigurationError

logger = logging.getLogger(__name__)


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
            resolved_path = os.path.join(".", run_name)
        else:
            raise MS2RescoreConfigurationError(
                "Could not resolve spectrum file name: No spectrum path configured "
                "and no run name in PSM file found."
            )

    # If passed path is directory, join with run name
    elif os.path.isdir(configured_path):
        if run_name:
            resolved_path = os.path.join(configured_path, run_name)
        else:
            raise MS2RescoreConfigurationError(
                "Could not resolve spectrum file name: Spectrum path is directory "
                "but no run name in PSM file found."
            )

    # If passed path is file, use that, but warn if basename doesn't match expected
    elif os.path.isfile(configured_path):
        if run_name and Path(configured_path).stem != Path(run_name).stem:
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
                "Resolved spectrum filename does not contain a supported file "
                "extension (mgf or mzml) and could not find any matching existing "
                "files."
            )

    return Path(resolved_path)
