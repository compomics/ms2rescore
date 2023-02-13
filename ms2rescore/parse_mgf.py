"""Parse MGF files."""

import logging
import mmap
import os.path
import random
import re

from rich.progress import track

from ms2rescore.exceptions import MS2RescoreError

logger = logging.getLogger(__name__)


class ParseMGFError(MS2RescoreError):
    """Error parsing MGF file."""

    pass

def parse_mgf_title_rt(
    path_to_mgf: Union[str, os.PathLike]
) -> Tuple[Dict[int, str], Dict[int, float]]:
    """Parse MGF file to extract title and retention time fields, by spectrum index."""
    title = None
    retention_times = dict()
    with open(path_to_mgf, "rt") as mgf_in:
        for line in mgf_in:
            if line[0] == "T":
                if line.startswith("TITLE="):
                    title = line[6:].strip()
            if line[0] == "R":
                if line.startswith("RTINSECONDS="):
                    if not title:
                        raise ParseMGFError("Missing `TITLE` for `RTINSECONDS` entry.")
                    retention_times[title] = float(line[12:].strip())
                    title = None  # Reset to detect potential missing titles
    return retention_times


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines
