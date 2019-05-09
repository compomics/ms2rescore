# Standard library
import os.path
import logging
import mmap

# Third party
from tqdm import tqdm


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


def title_parser(line, method='full'):
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
    - 'TRFP_MQ': For MGF parsed with ThermoRawFileParser and spec_ids from
      MaxQuant msms.txt.
    """

    if method == 'full':
        title = line[6:].strip()
    elif method == 'first_space':
        title = line[6:].split(' ')[0].strip()
    elif method == 'first_space_no_charge':
        title = '.'.join(line[6:].split(' ')[0].split('.')[:-1]).strip()
    elif method == 'TRFP_MQ':
        line = line.strip().split('mzspec=')[1].split(' ')
        filename = line[0].replace('.raw:', '')
        scan = line[3].replace('scan=', '')
        title = '.'.join([filename, scan, scan])
    else:
        raise ValueError("method '{}' is not a valid title parsing method".format(
            method
        ))
    return title


def parse_mgf(df_in, mgf_folder, outname='scan_mgf_result.mgf',
              filename_col='mgf_filename', spec_title_col='spec_id',
              title_parsing_method='full',
              show_progress_bar=True):

    if df_in[filename_col].iloc[0][-4:] in ['.mgf', '.MGF']:
        file_suffix = ''
    else:
        file_suffix = '.mgf'

    runs = df_in[filename_col].unique()
    logging.info("Parsing %i MGF files to single MGF containing all PSMs.", len(runs))

    with open(outname, 'w') as out:
        count = 0
        for run in runs:
            found = False
            current_mgf_file = os.path.join(mgf_folder, run + file_suffix)
            assert os.path.isfile(current_mgf_file), "MGF file {} could not be found.".format(current_mgf_file)

            spec_set = set(df_in[(df_in[filename_col] == run)][spec_title_col].values)

            # Temporary fix: replace charges in MGF with ID'ed charges
            # Until MS2PIP uses ID'ed charge instead of MGF charge
            id_charges = df_in[(df_in[filename_col] == run)].set_index('spec_id')['charge'].to_dict()

            with open(current_mgf_file, 'r') as f:
                iterator = tqdm(f, total=get_num_lines(current_mgf_file)) if show_progress_bar else f
                for line in iterator:
                    if 'TITLE=' in line:
                        title = title_parser(line, method=title_parsing_method)
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

    logging.info("%i/%i spectra found and written to new MGF file.", count, len(df_in))
    assert count == len(df_in), "Not all PSMs could be found in the provided MGF files"
