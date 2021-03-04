"""Parse MGF files."""

import logging
import mmap
import os.path
import re

from tqdm import tqdm

from ms2rescore._exceptions import MS2ReScoreError

logger = logging.getLogger(__name__)


class ParseMGFError(MS2ReScoreError):
    """Error parsing MGF file."""

    pass


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


def title_parser(line, method='full', run=None):
    """
    Take an MGF TITLE line and return the spectrum title.

    Depending on the software tool that was used to write the MGF files and the
    search engine, we need to extract different parts of the TITLE field. E.g,
    for MaxQuant, everything up until the first space is required

    line: string, line from MGF file containing 'TITLE='
    method: string, one of the following:
    - 'full': take everything after 'TITLE='
    - 'first_space': take everything between 'TITLE=' and first space.
    - 'first_space_no_charge': take everything between 'TITLE=' and first space,
      but leave out everything after last dot. (required for MaxQuant pipeline).
    - 'run.scan.scan': Extract scan number and merge with run name (for MaxQuant IDs).
    """

    if method == 'full':
        title = line[6:].strip()
    elif method == 'first_space':
        title = line[6:].split(' ')[0].strip()
    elif method == 'first_space_no_charge':
        title = '.'.join(line[6:].split(' ')[0].split('.')[:-1]).strip()
    elif method == 'run.scan.scan':
        if not run:
            raise TypeError("If `method` is `run.scan.scan`, `run` cannot be None.")
        scan_m = re.match(r"TITLE=.*scan=([0-9]+).*$", line)
        if scan_m:
            scan = scan_m.group(1)
        else:
            raise ParseMGFError(
                f"Could not extract scan number from TITLE field: `{line.strip()}`"
            )
        title = '.'.join([run, scan, scan])
    else:
        raise ValueError("method '{}' is not a valid title parsing method".format(
            method
        ))
    return title


def parse_mgf(df_in, mgf_folder, outname='scan_mgf_result.mgf',
              filename_col='mgf_filename', spec_title_col='spec_id',
              title_parsing_method='full',
              show_progress_bar=True):

    if not os.path.isdir(mgf_folder):
        raise NotADirectoryError(mgf_folder)

    if df_in[spec_title_col].duplicated().any():
        logger.warning("Duplicate spec_id's found in PeptideRecord.")

    if df_in[filename_col].iloc[0][-4:].lower() == '.mgf':
        file_suffix = ''
    else:
        file_suffix = '.mgf'

    runs = df_in[filename_col].unique()
    logger.info("Parsing %i MGF files to single MGF containing all PSMs.", len(runs))

    with open(outname, 'w') as out:
        count = 0
        for run in runs:
            current_mgf_file = os.path.join(mgf_folder, run + file_suffix)
            spec_set = set(df_in[(df_in[filename_col] == run)][spec_title_col].values)

            # Temporary fix: replace charges in MGF with ID'ed charges
            # Until MS2PIP uses ID'ed charge instead of MGF charge
            id_charges = df_in[(df_in[filename_col] == run)].set_index('spec_id')['charge'].to_dict()

            found = False
            with open(current_mgf_file, 'r') as f:
                iterator = tqdm(f, total=get_num_lines(current_mgf_file)) if show_progress_bar else f
                for line in iterator:
                    if 'TITLE=' in line:
                        title = title_parser(line, method=title_parsing_method, run=run)
                        if title in spec_set:
                            found = True
                            line = "TITLE=" + title + "\n"
                            out.write("BEGIN IONS\n")
                            out.write(line)
                            count += 1
                            continue
                    if 'END IONS' in line:
                        if found:
                            out.write(line + '\n')
                            found = False
                            continue
                    # Temporary fix (see above)
                    if 'CHARGE=' in line:
                        if found:
                            charge = id_charges[title]
                            out.write("CHARGE=" + str(charge) + "+\n")
                            continue
                    # Only print lines when spectrum is found and intensity != 0
                    if found and line[-4:] != '0.0\n':
                        out.write(line)

    num_expected = len(df_in[spec_title_col].unique())
    logger.debug(
        "%i/%i spectra found and written to new MGF file.", count, num_expected
    )
    if not count == num_expected:
        raise ParseMGFError("Not all PSMs could be found in the provided MGF files.")
