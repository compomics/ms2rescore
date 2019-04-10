__author__ = "Ralf Gabriels"
__credits__ = ["Ralf Gabriels", "Sven Degroeve", "Lennart Martens"]
__license__ = "Apache License, Version 2.0"
__version__ = "0.1.1"
__email__ = "Ralf.Gabriels@UGent.be"


# Native
import argparse
import logging
import mmap
import re

# Third party
import numpy as np
import pandas as pd

# tqdm
try:
    from tqdm import tqdm
except ImportError:
    USE_TQDM = False
else:
    USE_TQDM = True


def argument_parser():
    parser = argparse.ArgumentParser(description='Extract features and MS2PIP input data from MaxQuant msms.txt file for ReScore.')
    parser.add_argument('msms_filename', metavar='<msms.txt file>',
                        help='Path to msms.txt file.')
    parser.add_argument('-m', metavar='<mgf folder>', dest='mgf_folder', action='store', default=None,
                        help='Folder containing MGF files to parse.')
    parser.add_argument('-o', metavar='<output filename>', dest='outname', action='store', default='maxquant_to_rescore_out',
                        help='Output filename for PEPREC and parsed MGF.')
    args = parser.parse_args()
    return args


def calc_top7_peak_features(intens, mass_errors):
    """
    Calculate top7-related features for Percolator.

    MeanErrorTop7: Mean of mass errors of the seven fragment ion peaks with the highest intensities
    sqMeanErrorTop7: Squared MeanErrorTop7
    StdevErrorTop7: Standard deviation of mass errors of the seven fragment ion peaks with the highest intensities
    """
    if (type(intens) != list) or (type(mass_errors) != list):
        return np.nan, np.nan, np.nan
    else:
        intens = [float(i) for i in intens]
        mass_errors = [float(i) for i in mass_errors]

        indices_most_intens = np.array(intens).argsort()[-1:-8:-1]
        mass_errors_top7 = [(mass_errors[i]) for i in indices_most_intens]
        mean_error_top7 = np.mean(mass_errors_top7)
        sq_mean_error_top7 = mean_error_top7 ** 2
        stdev_error_top7 = np.std(mass_errors_top7)

        return mean_error_top7, sq_mean_error_top7, stdev_error_top7


def calc_ion_current_features(matches, intensities, intensity_coverage):
    """
    Calculate ion current related features for Percolator.

    lnExplainedIonCurrent: Summed intensity of identified fragment ions, divided by that of all fragment ions, logged
    lnNTermIonCurrentRatio: Summed intensity of identified N-terminal fragments, divided by that of all identified fragments, logged
    lnCTermIonCurrectRatio: Summed intensity of identified N-terminal fragments, divided by that of all identified fragments, logged
    lnMS2IonCurrent: Summed intensity of all observed fragment ions, logged
    """
    if type(intensities) != list:
        return np.nan, np.nan, np.nan, np.nan
    else:
        lnExplainedIonCurrent = intensity_coverage
        summed_intensities = sum([float(i) for i in intensities])

        # Calculate ratio between matched b- and y-ion intensities
        y_ion_int = sum([float(intensities[i]) for i, m in enumerate(matches) if m.startswith('y')])
        y_int_ratio = y_ion_int / summed_intensities

        lnNTermIonCurrentRatio = y_int_ratio * lnExplainedIonCurrent
        lnCTermIonCurrentRatio = (1 - y_int_ratio) * lnExplainedIonCurrent
        lnMS2IonCurrent = summed_intensities / lnExplainedIonCurrent

        out = [lnExplainedIonCurrent, lnNTermIonCurrentRatio, lnCTermIonCurrentRatio, lnMS2IonCurrent]

    return tuple([np.log2(x) for x in out])


def msms_to_peprec(msms_filename, fixed_modifications=None, ptm_mapping=None,
                   validate_amino_acids=True):
    """
    Extract features and MS2PIP input data from MaxQuant msms.txt file for ReScore.'

    Percolator features are derived from the MSGF2PIN script. See table 1 of
    Percolator-MSGF+ article (doi.org/10.1021/pr400937n).

    Positional arguments:
    `msms_file`: str with the file location of the MSMS.txt file

    Keyword arguments:
    `ptm_mapping` (dict) is used to convert the MaxQuant PTM labels to PSI-MS
    modification names. For correct parsing, the key should always include the
    two brackets.
    `fixed_modifications` (list of tuples, [(aa, ptm)]) can contain fixed
    modifications to be added to the peprec. E.g. `[('C', 'cm')]`. The first
    tuple element contains the one-letter amino acid code. The second tuple
    element contains a two-character label for the PTM. This PTM also needs
    to be present in the `ptm_mapping` dictionary.
    `validate_amino_acids`: Remove PSMs where the sequence includes an invalid
    amino acid (B, J, O, U, X, Z); required for MS2PIP compatibility.
    """

    if not fixed_modifications:
        fixed_modifications = [
            ('C', 'cm'),
        ]
    if not ptm_mapping:
        ptm_mapping = {
            '(ox)': 'Oxidation',
            '(ac)': 'Acetyl',
            '(cm)': 'Carbamidomethyl'
        }

    msms_cols = [
        'Raw file', 'Scan number', 'Charge', 'Length', 'Sequence', 'Modified sequence',
        'Proteins', 'Missed cleavages', 'Mass', 'Mass error [Da]',
        'Reverse', 'PEP', 'Score', 'Delta score', 'Localization prob', 'Matches',
        'Intensities', 'Mass Deviations [Da]', 'Intensity coverage', 'id',
    ]

    logging.debug("Reading msms file")
    msms = pd.read_csv(msms_filename, sep='\t', usecols=msms_cols)

    # Filter for rank 1 PSMs
    msms = msms.sort_values('Score', ascending=False)
    msms = msms[~msms.duplicated(['Raw file', 'Scan number'], keep='first')]
    msms = msms.sort_index().reset_index(drop=True)

    # Remove PSMs with invalid amino acids
    if validate_amino_acids:
        msms = msms[~(msms['Sequence'].str.contains('[BJOUXZ]', regex=True))]\
            .reset_index(drop=True)

    logging.info("Found {} PSMs of which {:.0%} are decoy hits.".format(
        len(msms),
        len(msms[msms['Reverse'] == '+']) / len(msms)
    ))

    logging.debug("Calculating Search Engine features")
    # Calculate peak intensity related features
    top7_features = pd.DataFrame(
        [calc_top7_peak_features(i, md)
        for i, md in zip(
            msms['Intensities'].str.split(';'), msms['Mass Deviations [Da]'].str.split(';')
        )],
        columns=['MeanErrorTop7', 'sqMeanErrorTop7', 'StdevErrorTop7']
    )

    # Ion current related features
    ion_current_features = pd.DataFrame(
        [calc_ion_current_features(m, i, ic) for m, i, ic in zip(
            msms['Matches'].str.split(';'),
            msms['Intensities'].str.split(';'),
            msms['Intensity coverage'],
        )], columns=['lnExplainedIonCurrent', 'lnNTermIonCurrentRatio', 'lnCTermIonCurrentRatio', 'lnMS2IonCurrent']
    )

    # Other features
    msms['absdM'] = msms['Mass error [Da]'].abs()

    # Other MS2PIP / ReScore columns
    msms['Label'] = msms['Reverse'].isna().apply(lambda x: 1 if x else -1)
    msms['charge_ms2pip'] = msms['Charge']
    msms['spec_id'] = msms['Raw file'] + '.' + msms['Scan number'].astype(str) + '.' + msms['Scan number'].astype(str)
    msms['Proteins'] = msms['Proteins'].str.split(';')
    msms['Peptide'] = msms['Sequence']

    # Parse modifications for MS²PIP
    logging.debug("Parsing modifications to MS2PIP format")
    for aa, mod in fixed_modifications:
        msms['Modified sequence'] = msms['Modified sequence'].str.replace(aa, '{}({})'.format(aa, mod))
    pattern = r'\([a-z].\)'
    msms['Parsed modifications'] = ['|'.join(['{}|{}'.format(m.start(0) - 1 - i*4, ptm_mapping[m.group()]) for i, m in enumerate(re.finditer(pattern, s))]) for s in msms['Modified sequence']]
    msms['Parsed modifications'] = ['-' if mods == '' else mods for mods in msms['Parsed modifications']]

    # Bringing it all together
    msms = pd.concat([msms.reset_index(drop=True), top7_features, ion_current_features], axis=1)
    peprec_columns = ['spec_id', 'Parsed modifications', 'Sequence', 'charge_ms2pip']
    percolator_columns = [
        'Label', 'Peptide', 'Proteins', 'Score', 'Delta score',
        'Localization prob', 'PEP', 'lnExplainedIonCurrent',
        'lnNTermIonCurrentRatio', 'lnCTermIonCurrentRatio', 'lnMS2IonCurrent',
        'Mass', 'Length', 'Mass error [Da]', 'absdM', 'MeanErrorTop7',
        'sqMeanErrorTop7', 'StdevErrorTop7', 'Charge', 'Missed cleavages'
    ]
    col_mapping = {
        'Parsed modifications': 'modifications',
        'Sequence': 'peptide',
        'charge_ms2pip': 'charge',
        'Score': 'RawScore',
        'Delta score': 'RawDeltaScore',
        'Localization prob': 'RawModLocProb',
        'PEP': 'MaxQuantPEP',
        'Mass': 'Mass',
        'Length': 'PepLen',
        'Mass error [Da]': 'dM',
        'Charge': 'ChargeN',
        'Missed cleavages': 'enzInt',
    }

    peprec_percolator = msms[peprec_columns + percolator_columns + ['Raw file']].rename(columns=col_mapping)
    logging.debug("Finished parsing msms.txt file")

    return peprec_percolator


def scan_mgf(df_in, mgf_folder, outname='scan_mgf_result.mgf',
             filename_col='mgf_filename', spec_title_col='spec_id',
             use_tqdm=False):


    def get_num_lines(file_path):
        fp = open(file_path, "r+")
        buf = mmap.mmap(fp.fileno(), 0)
        lines = 0
        while buf.readline():
            lines += 1
        return lines


    if df_in[filename_col].iloc[0][-4:] in ['.mgf', '.MGF']:
        file_suffix = ''
    else:
        file_suffix = '.mgf'

    with open(outname, 'w') as out:
        count_runs = 0
        count = 0
        runs = df_in[filename_col].unique()
        print("Scanning MGF files: {} runs to do. Now working on run: ".format(len(runs)), end='')
        for run in runs:
            count_runs += 1
            if count_runs % 10 == 0:
                print(str(count_runs), end='')
            else:
                print('.', end='')

            spec_set = set(df_in[(df_in[filename_col] == run)][spec_title_col].values)

            # Temporary fix: replace charges in MGF with ID'ed charges
            # Until MS2PIP uses ID'ed charge instead of MGF charge
            id_charges = df_in[(df_in[filename_col] == run)].set_index('spec_id')['charge'].to_dict()

            found = False
            current_mgf_file = '{}/{}{}'.format(mgf_folder, str(run), file_suffix)
            with open(current_mgf_file, 'r') as f:
                if use_tqdm:
                    mgf_iterator = tqdm(f, total=get_num_lines(current_mgf_file))
                else:
                    mgf_iterator = f
                for line in mgf_iterator:
                    if 'TITLE=' in line:
                        # Take everything between `TITLE=` and the first space
                        # as the title, but exclude the charge, as the measured
                        # charge is not always the identified charge...
                        title = '.'.join(line[6:].split(' ')[0].split('.')[:-1])
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

    print("\n{}/{} spectra found and written to new MGF file.".format(count, len(df_in)))


def main():
    args = argument_parser()
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.DEBUG
    )
    logging.info("Parsing msms.txt file")
    peprec_percolator = msms_to_peprec(args.msms_filename)

    # If MGF folder is provide, scan MGF files for spectra to include in the one MGF
    if args.mgf_folder:
        logging.info("Parsing MGF files")
        scan_mgf(
            peprec_percolator,
            args.mgf_folder, outname='{}.mgf'.format(args.outname),
            filename_col='Raw file', spec_title_col='spec_id', use_tqdm=USE_TQDM
        )

    logging.info("Writing PEPREC file")
    peprec_percolator.drop('Raw file', axis=1, inplace=True)
    peprec_percolator.to_csv('{}.peprec'.format(args.outname), sep=' ', index=False)

if __name__ == '__main__':
    main()
