"""Parse MGF files."""

import logging
import mmap
import os.path
from typing import Union, Tuple, Dict

from rich.progress import track

from ms2rescore.exceptions import MS2RescoreError

logger = logging.getLogger(__name__)


class ParseMGFError(MS2RescoreError):
    """Error parsing MGF file."""

    pass

def parse_mgf_title_rt(
    path_to_mgf: Union[str, os.PathLike]
) -> Dict[str, float]:
    """Parse MGF file to extract title and retention time fields, by spectrum index."""
    logger.debug("Parsing MGF file to extract retention times.")

    title = None
    retention_time = None
    retention_times = dict()
    spectrum_header = False
    with open(path_to_mgf, "rt") as mgf_in:
        for line in mgf_in:
            if not line[0].isdigit():
                spectrum_header = True
                if line[0] == "T":
                    if line.startswith("TITLE="):
                        title = line[6:].strip()
                elif line[0] == "R":
                    if line.startswith("RTINSECONDS="):
                        retention_time = float(line[12:].strip())
            elif line[0].isdigit() and spectrum_header:
                if title is None:
                    raise ParseMGFError("Missing `TITLE` for `RTINSECONDS` entry.")
                if title is None:
                    raise ParseMGFError("Missing `RTINSECONDS` for `TITLE` entry.")
                retention_times[title] = retention_time
                title = None
                spectrum_header = False
            else:
                raise ParseMGFError("Expected spectrum header for spectrum")

    return retention_times


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines
