# Standard library
import logging
import subprocess
import os
import re

# Third party
import pandas as pd

# Package
import ms2rescore.mapper as mapper

def run_msgfplus(msgfplus_config, outname, log=False, num_cpu=23):
    """
    Runs MSGFPlus with some fixed settings: concatenated search, include
    additional features, no protocol. MSGF+ is called with subprocess.

    :param msgfplus_config: dictionary, contains search engine options
    :param outname: string, name for output files
    :param log: boolean, write MSGFPlus log to file or not
    """

    if msgfplus_config['search_params']['frag'] == 'HCD':
        m = 3
        inst = 1
    elif msgfplus_config['search_params']['frag'] == 'CID':
        m = 1
        inst = 0
    if msgfplus_config['search_params']['path_to_modsfile'] != '':
        mods = '-mod {} '.format(msgfplus_config['search_params']['path_to_modsfile'])
    else:
        mods = ''
    if log:
        log_cmd = " > {}.log".format(outname)
    else:
        log_cmd = ''

    jar_file = msgfplus_config['search_params']['jar_file']
    mgf_file = msgfplus_config['mgf_file']
    fasta_file = msgfplus_config['search_params']['fasta_file']
    min_length = msgfplus_config['search_params']['min_length']
    min_charge = msgfplus_config['search_params']['min_charge']
    max_charge = msgfplus_config['search_params']['max_charge']
    ms1_tol = msgfplus_config['search_params']['ms1_tolerance']

    msgf_command = f"java -Xmx28000M -jar {jar_file} {mods}-s {mgf_file} \
        -d {fasta_file} -o {outname}.mzid -t {ms1_tol} -tda 1 -m {m} \
        -inst {inst} -minLength {min_length} -minCharge {min_charge} \
        -maxCharge {max_charge} -n 1 -e 1 -addFeatures 1 -protocol 0 \
        -thread {num_cpu}{log_cmd}"

    logging.info("MSGFPlus command: %s", msgf_command)
    subprocess.run(msgf_command, shell=True, check=True)


def make_pepfile(path_to_pin, options):
    """
    Read a pin file and create the corresponding pepfile dataframe, which will
    be saved as a MS2PIP PEPREC file.

    :param path_to_pin: path to pin file

    Returns
    :pd.DataFrame pepfile, columns ['TITLE', 'Peptide', 'Charge', 'Label',
        'modifications', 'Proteins']
    """
    pin = pd.read_csv(path_to_pin, sep='\t', header=0, skiprows=[1])

    charge_states = []
    for a in pin.columns:
        if a.startswith('Charge'):
            charge_states.append(a)
        else: continue

    pin.loc[:, 'Charge'] = [None] * len(pin)

    for ch in charge_states:
        value = int(ch.lstrip('Charge'))
        pin.loc[pin[ch] == 1, 'Charge'] = value

    pepfile = pin.loc[:, ['TITLE', 'Peptide', 'Charge', 'Label', 'Proteins']]

    # Add modifications column to PEPREC file
    # the keys correspond to the UNIMOD keys for each modification
    modifications = {}
    mods = options["ms2pip"]["modifications"]
    for mod in mods:
        modifications[str(mod["unimod_accession"])] = mod["name"]

    modlist = []
    # TODO get rid of iterrows!
    for _, row in pepfile.iterrows():
        if 'UNIMOD' in row['Peptide']:
            pep = row['Peptide'].split('.')[1]
            mods = re.findall(r'\[([^]]*)\]', pep)
            modstring = ''
            for mod in mods:
                mod = '[' + mod + ']'
                key = mod.split(':')[1].rstrip(']')
                try:
                    if key == '21':
                        phospholoc = pep[pep.find(mod)-1]
                        modstring += str(pep.find(mod)) + '|' + modifications[key] + phospholoc + '|'
                        pep = pep.replace(mod, '', 1)
                    else:
                        modstring += str(pep.find(mod)) + '|' + modifications[key] + '|'
                        pep = pep.replace(mod, '', 1)
                except:
                    print('Modification not expected: {}'.format(mod))
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

    write_PEPREC(pepfile, path_to_pin)


def write_PEPREC(pepfile, path_to_pep, concat=True):
    """
    Write the PEPREC file, which will be the input to MS2PIP.

    :param pepfile: pd.DataFrame created by rescore.make_pepfile()
    :param path_to_pep: string, path where to save the PEPREC file(s)
    :param concat: boolean, True if the search was concatenated
    """

    if concat:
        pepfile_tosave = pepfile.loc[:, ['TITLE', 'modifications', 'peptide', 'Charge']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge']
        pepfile_tosave.to_csv(path_to_pep + '.peprec', sep=' ', index=False)

    else:
        pepfile_tosave = pepfile.loc[pepfile.Label == 1, ['TITLE', 'modifications', 'peptide', 'Charge']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge']
        pepfile_tosave.to_csv(path_to_pep + '.targets.peprec', sep=' ', index=False)

        pepfile_tosave = pepfile.loc[pepfile.Label == -1, ['TITLE', 'modifications', 'peptide', 'Charge']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge']
        pepfile_tosave.to_csv(path_to_pep + '.decoys.peprec', sep=' ', index=False)

    return None


def join_features(path_to_pin, path_to_pep):
    pin = pd.read_csv(path_to_pin, sep='\t')
    pep = pd.read_csv(path_to_pep, sep=' ')

    # Reorganize PIN file to have all features after ID cols (as it has to be
    # in PEPREC)
    pin.drop('SpecId', axis=1, inplace=True)
    pin.rename(columns={'TITLE': 'spec_id', 'Peptide': 'ModPeptide'}, inplace=True)
    id_cols = ['spec_id', 'Label', 'ScanNr', 'ModPeptide', 'Proteins']
    feature_cols = [col for col in pin.columns if col not in id_cols]
    pin = pin[id_cols + feature_cols]

    # Merge PIN to PEPREC
    pep = pep.merge(pin, on='spec_id')

    pep.to_csv(path_to_pep, sep=' ', index=False)


def msgf_pipeline(config):
    outname = config['general']['output_filename']
    mzid_file = config['general']['identification_file']
    if not os.path.isfile(config['general']['mgf_file']):
        logging.critical(
            "MGF file %s not found. Please specify the correct path to the MGF file.",
            config['general']['mgf_file']
        )
        exit(1)

    logging.info("Running msgf2pin")
    # Convert .mzid to pin. XXX is the decoy pattern from MSGF+
    convert_command = "msgf2pin -P XXX {} > {}_original.pin".format(mzid_file, outname)
    subprocess.run(convert_command, shell=True)

    # PIN FILE: "Proteins" column has tab-separated values which makes the file
    # cumbersome to read. mapper.fix_pin_tabs replaces those tabs with ";"
    logging.info("Fixing tabs on pin file")
    mapper.fix_pin_tabs(outname + "_original.pin", prot_sep='|||')
    os.rename(outname + "_original_fixed.pin", outname + "_edited.pin")

    # Percolator generates its own spectrum ID, so we must match it to the mgf
    # file's TITLE field.
    logging.info("Adding mgf TITLE to pin file")
    mapper.map_mgf_title(outname + "_edited.pin", mzid_file, msgs=False)
    os.rename(outname + "_edited.pin_title", outname + "_edited.pin")

    logging.info("Writing PEPREC file")
    make_pepfile(outname + "_edited.pin", config)
    os.rename(outname + "_edited.pin.peprec", outname + ".peprec")

    # Add features from the pin file to the PEPREC file
    logging.info("Adding pin features to PEPREC file")
    join_features(outname + "_edited.pin", outname + ".peprec")
    os.remove(outname + '_edited.pin')

    return outname + ".peprec", config['general']['mgf_file']
