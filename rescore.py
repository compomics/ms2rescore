"""
Functions necessary to run the rescore algorithm. Currently supports MSGF+ with
concatenated searches.
"""

import subprocess
import sys
import re
import pandas as pd

from mapper import mapper  # TODO shouldn't have to do this, check __init__.py


def run_msgfplus(msgf_dir, outfile, mgffile, fastafile, modsfile, frag='HCD'):
    """
    Runs MSGF+ with some fixed settings: 10ppm precursor mass tolerance,
    concatenated search, minimum peptide length of 8aa, minimum charge is 2 and
    maximum is 4, 1 match per spectrum, include additional features, no protocol
    and use 23 threads. MSGF+ is called with subprocess.

    :param msgfdir: string, path to MSGFPlus.jar
    :param outfile: string, path where to store the search results (.mzid)
    :param mgffile: string, spectrum file (MGF,PKL,DTA,mzXML,mzDATA or mzML)
    :param fastafile: string, the database to search against in .fasta format
    :param modsfile: string, path to the file with modifications
    :param frag: fragmentation method (HCD or CID) on which some settings depend
    """
    # TODO pipe the output somewhere to make the terminal look cleaner
    if frag == 'HCD':
        m = 3
        inst = 1
    elif frag == 'CID':
        m = 1
        inst = 0
    if modsfile != '':
        mods = '-mod {} '.format(modsfile)
    else:
        mods = ''
    msgf_command = "java -Xmx8000M -jar {}/MSGFPlus.jar {}-s {} -d {} -o {} -t \
        10ppm -tda 1 -m {} -inst {} -minLength 8 -minCharge 2 -maxCharge 4 -n \
        1 -addFeatures 1 -protocol 0 -thread 23 \n".format(msgf_dir, mods,
        mgffile, fastafile, outfile + ".mzid", m, inst)

    sys.stdout.write("Running search with MSGF+: {}".format(msgf_command))
    sys.stdout.flush()
    subprocess.run(msgf_command, shell=True)

    return None


def make_pepfile(path_to_pin, modsfile=None):
    """
    Read a pin file and create the corresponding pepfile dataframe, which will
    be saved as a MS2PIP PEPREC file.

    :param path_to_pin: path to pin file

    Returns
    :pd.DataFrame pepfile, columns ['TITLE', 'Peptide', 'Charge', 'Label',
        'modifications']
    """
    pin = pd.read_csv(path_to_pin, sep='\t', header=0, skiprows=[1])
    # TODO these columns depend on the file, so should not be hard coded
    pin.loc[:, 'Charge'] = [2 if r[1].Charge2 == 1
                            else 3 if r[1].Charge3 == 1
                            else 4 if r[1].Charge4 == 1
                            else 5 if r[1].Charge5 == 1
                            else 6 if r[1].Charge6 == 1
                            else 'X' for r in pin.iterrows()]

    pepfile = pin[['TITLE', 'Peptide', 'Charge', 'Label']]

    # Add modifications column to PEPREC file
    # TODO: read modifications from MSGF+ modifications file
    # the keys correspond to the UNIMOD keys for each modification
    modifications = {'4': 'Cam', '35': 'Oxidation'}

    modlist = []
    for _, row in pepfile.iterrows():
        if 'UNIMOD' in row['Peptide']:
            pep = row['Peptide'].split('.')[1]
            mods = re.findall(r'\[([^]]*)\]', pep)
            modstring = ''
            for mod in mods:
                mod = '[' + mod + ']'
                modstring += str(pep.find(mod)) + '|' + \
                    modifications[mod.split(':')[1].rstrip(']')] + '|'
                pep = pep.replace(mod, '', 1)
            modlist.append(modstring.rstrip('|'))
        else:
            modlist.append('')
    pepfile.loc[:, 'modifications'] = modlist

    # Rewrite peptide sequences without the UNIMOD modification ids
    pep_list = []
    for _, row in pepfile.iterrows():
        pep = row['Peptide']
        pep = pep.split('.')[1]
        if 'UNIMOD' in pep:
            mods = re.findall(r'\[([^]]*)\]', pep)
            for mod in mods:
                pep = pep.replace('[' + mod + ']', '', 1)
        pep_list.append(pep)
    pepfile.loc[:, 'peptide'] = pep_list

    return pepfile


def write_PEPREC(pepfile, path_to_pep, concat=True):
    """
    Write the PEPREC file, which will be the input to MS2PIP.

    :param pepfile: pd.DataFrame created by rescore.make_pepfile()
    :param path_to_pep: string, path where to save the PEPREC file(s)
    :param concat: boolean, True if the search was concatenated
    """
    # TODO SettingWithoutCopyWarning Try using .iloc instead bla bla
    if concat:
        pepfile_tosave = pepfile[['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
        pepfile_tosave.to_csv(path_to_pep + '.PEPREC', sep=' ', index=False)

    else:
        pepfile_tosave = pepfile[pepfile.Label == 1][['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
        pepfile_tosave.to_csv(path_to_pep + '.targets.PEPREC', sep=' ', index=False)

        pepfile_tosave = pepfile[pepfile.Label == -1][['TITLE', 'modifications','peptide', 'Charge', 'Label']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
        pepfile_tosave.to_csv(path_to_pep + '.decoys.PEPREC', sep=' ', index=False)

    return None


def join_features(path_to_target_features, path_to_pin, path_to_decoy_features=None):
    """
    Combine the features table obtained from Percolator (i.e. in the pin file)
    and the ones obtained from MS2PIP_rescore in one DataFrame

    :param path_to_target_features: string, path to MS2PIP features
    :param path_to_pin: string, path to pin file
    :param path_to_decoy_features: string, if non-concatenated search then there
         will be a separate feature file with the decoy features

    Returns
    :pd.DataFrame all_features, includes all the Percolator and MS2PIP features
    """
    # read pin file - should not need the lazy pin parser as this pin already
    # has the TITLE which means it was processed by mapper
    pin = pd.read_csv(path_to_pin, sep='\t')

    # Read rescore_features.csv file and fillna
    rescore_targets = pd.read_csv(path_to_target_features)
    rescore_targets = rescore_targets.fillna(0)
    # If not concat searches, do that for target and decoy files
    if path_to_decoy_features != None:
        rescore_decoys = pd.read_csv(path_to_decoy_features)
        rescore_decoys = rescore_decoys.fillna(0)

        # join target and decoy tables
        all_features = pd.concat([rescore_decoys.merge(pin[pin.Label == -1], left_on='spec_id', right_on='TITLE'),
                                  rescore_targets.merge(pin[pin.Label == 1], left_on='spec_id', right_on='TITLE')])
    else:
        all_features = rescore_targets.merge(pin, left_on='spec_id', right_on='TITLE')

    return all_features


def write_pin_files(all_features, savepath):
    """
    Given a dataframe with all the features, writes five different pin files:
    _only_rescore.pin with only the rescore features
    _all_percolator.pin with all the percolator features
    _percolator_default.pin with only the default percolator features
    _all_features.pin with all the rescore and percolator features
    _default_and_rescore.pin with the rescore and the default percolator features

    :param all_features: pd.DataFrame obtained from rescore.join_features()
    """
    # columns to save
    rescore_features = list(all_features.columns[4:42])
    percolator_features = list(all_features.columns[45:-5])
    percolator_default = percolator_features[:27]

    # Writing files with appropriate columns
    all_features[['SpecId', 'Label', 'ScanNr'] + rescore_features + ['Peptide',
                'Proteins']].to_csv('{}_only_rescore.pin'.format(savepath), sep='\t', index=False)
    all_features[['SpecId', 'Label', 'ScanNr'] + percolator_features + ['Peptide',
                'Proteins']].to_csv('{}_all_percolator.pin'.format(savepath), sep='\t', index=False)
    all_features[['SpecId', 'Label', 'ScanNr'] + percolator_default + ['Peptide', 'Proteins']
                 ].to_csv('{}_percolator_default.pin'.format(savepath), sep='\t', index=False)
    all_features[['SpecId', 'Label', 'ScanNr'] + percolator_features + rescore_features +
                 ['Peptide', 'Proteins']].to_csv('{}_all_features.pin'.format(savepath), sep='\t', index=False)
    all_features[['SpecId', 'Label', 'ScanNr'] + percolator_default + rescore_features + ['Peptide',
                 'Proteins']].to_csv('{}_default_and_rescore.pin'.format(savepath), sep='\t', index=False)

    return None
