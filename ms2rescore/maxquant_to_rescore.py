# Native
import logging
import re

# Third party
import numpy as np
import pandas as pd

# Project
import ms2rescore


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


def msms_to_peprec(msms_filename, modifications_mapping=None,
                   fixed_modifications=None, validate_amino_acids=True):
    """
    Extract features and MS2PIP input data from MaxQuant msms.txt file for ReScore.'

    Percolator features are derived from the MSGF2PIN script. See table 1 of
    Percolator-MSGF+ article (doi.org/10.1021/pr400937n).

    Positional arguments:
    `msms_file`: str with the file location of the MSMS.txt file

    Keyword arguments:
    `modifications_mapping` (dict) is used to convert the two-letter MaxQuant
    modification labels to PSI-MS modification names.
    `fixed_modifications` (dict, {aa: mod}) can contain fixed
    modifications to be added to the peprec. E.g. `{'C': 'Carbamidomethyl'}`,
    as the MaxQuant output does not include modifications that were set as fixed
    during the search. The first tuple element contains the one-letter amino
    acid code. The second tuple element contains the full modification name, as
    listed in the values of `modifications_mapping`.
    `validate_amino_acids`: Remove PSMs where the sequence includes an invalid
    amino acid (B, J, O, U, X, Z); required for MS2PIP compatibility.
    """

    # Parse modification input
    if not modifications_mapping:
        modifications_mapping = {}
    else:
        modifications_mapping = {"({})".format(k): v for k, v in modifications_mapping.items()}

    if not fixed_modifications:
        fixed_modifications = []
    else:
        modifications_mapping_rev = {v: k for k, v in modifications_mapping.items()}
        fixed_modifications = [(k, modifications_mapping_rev[v]) for k, v in fixed_modifications.items()]

    msms_cols = {
        'Raw file', 'Scan number', 'Charge', 'Length', 'Sequence', 'Modified sequence',
        'Proteins', 'Missed cleavages', 'Mass', 'Mass error [Da]',
        'Reverse', 'PEP', 'Score', 'Delta score', 'Localization prob', 'Matches',
        'Intensities', 'Mass Deviations [Da]', 'Intensity coverage', 'id',
    }

    # Read first line of msms, to check is Mass Error is present in Da or ppm
    with open(msms_filename, 'rt') as f:
        header = f.readline().split('\t')
        if 'Mass error [Da]' in header:
            mass_error_unit = 'Da'
        elif 'Mass Error [ppm]' in header:
            mass_error_unit = 'ppm'
            msms_cols.remove('Mass error [Da]')
            msms_cols.add('Mass Error [ppm]')
        f.seek(0)

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

    logging.info("Found {} rank 1 PSMs of which {:.0%} are decoy hits.".format(
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
    msms['absdM'] = msms['Mass error [{}]'.format(mass_error_unit)].abs()

    # Other MS2PIP / ReScore columns
    msms['Label'] = msms['Reverse'].isna().apply(lambda x: 1 if x else -1)
    msms['charge_ms2pip'] = msms['Charge']
    msms['spec_id'] = msms['Raw file'] + '.' + msms['Scan number'].astype(str) + '.' + msms['Scan number'].astype(str)
    msms['Proteins'] = msms['Proteins'].str.split(';')
    msms['Peptide'] = msms['Sequence']

    # Fill NaN values in Proteins column for decoy PSMs
    # But first check that NaN Proteins only occur for decoy PSMs, if so:
    # fill these without the "REV_"
    if (msms['Proteins'].isna() & msms['Reverse'].isna()).any():
        req_cols = zip(msms['Proteins'], msms['Reverse'], msms['Modified sequence'])
        msms['Proteins'] = [[modseq] if (type(rev) == float) & (type(prot) == float) else prot for prot, rev, modseq in req_cols]

    msms['Proteins'] = msms['Proteins'].fillna('REV_' + msms['Modified sequence'])

    # Parse modifications for MS²PIP
    logging.debug("Parsing modifications to MS2PIP format")
    for aa, mod in fixed_modifications:
        msms['Modified sequence'] = msms['Modified sequence'].str.replace(aa, '{}{}'.format(aa, mod))
    pattern = r'\([a-z].\)'
    msms['Parsed modifications'] = ['|'.join(['{}|{}'.format(m.start(0) - 1 - i*4, modifications_mapping[m.group()]) for i, m in enumerate(re.finditer(pattern, s))]) for s in msms['Modified sequence']]
    msms['Parsed modifications'] = ['-' if mods == '' else mods for mods in msms['Parsed modifications']]

    # Bringing it all together
    msms = pd.concat([msms.reset_index(drop=True), top7_features, ion_current_features], axis=1)
    peprec_columns = ['spec_id', 'Parsed modifications', 'Sequence', 'charge_ms2pip']
    percolator_columns = [
        'Label', 'Modified sequence', 'Proteins', 'Score', 'Delta score',
        'Localization prob', 'PEP', 'lnExplainedIonCurrent',
        'lnNTermIonCurrentRatio', 'lnCTermIonCurrentRatio', 'lnMS2IonCurrent',
        'Mass', 'Length', 'Mass error [Da]', 'absdM', 'MeanErrorTop7',
        'sqMeanErrorTop7', 'StdevErrorTop7', 'Charge', 'Missed cleavages'
    ]
    col_mapping = {
        'Parsed modifications': 'modifications',
        'Sequence': 'peptide',
        'Modified sequence': 'ModPeptide',
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


def maxquant_pipeline(config):
    """
    Interface to maxquant_to_rescore: Prepare PEPREC and single MGF for
    MS2ReScore.
    """
    logging.info("Parsing msms.txt file")
    outname = config['general']['output_filename']
    peprec = msms_to_peprec(
        config['general']['identification_file'],
        modifications_mapping=config['maxquant_to_rescore']['modifications_mapping'],
        fixed_modifications=config['maxquant_to_rescore']['fixed_modifications'],
        validate_amino_acids=True
    )
    peprec.to_csv(outname + '.peprec', sep=' ', index=False)

    logging.info("Parsing MGF files")
    ms2rescore.parse_mgf.parse_mgf(
        peprec, config['maxquant_to_rescore']['mgf_dir'],
        outname=outname + '.mgf',
        filename_col='Raw file', spec_title_col='spec_id',
        title_parsing_method='TRFP_MQ',
        show_progress_bar=config['general']['show_progress_bar']
    )

    peprec.drop('Raw file', axis=1, inplace=True)
    peprec.to_csv(outname + '.peprec', sep=' ', index=False)

    return outname + '.peprec', outname + '.mgf'
