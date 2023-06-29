"""Parse MGF files."""

import logging
import mmap
import os.path
from typing import Union, Tuple, Dict

from rich.progress import track
from pyteomics.mgf import MGF

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
    mgf_reader = MGF(path_to_mgf, read_charges=False, read_ions=False)
    retention_times = {}
    for spectrum in mgf_reader:
        try:
            title = spectrum["params"]["title"]
        except KeyError:
            raise ParseMGFError("MGF file missing title field.")
        try:
            rt = float(spectrum["params"]["rtinseconds"])
        except KeyError:
            rt = None
        retention_times[title] = rt

    print(retention_times)
    if any(list(retention_times.values())):
        return retention_times
    else:
        raise ParseMGFError("MGF file missing rtinseconds field.")

def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines
