# Standard library
import logging
import subprocess
import os
import re

# Third party
import pandas as pd

# Package
from ms2rescore import mapper

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

    pepfile = pin.loc[:, ['SpecId', 'Peptide', 'Charge', 'Label', 'Proteins']]

    pepfile.to_csv(open(path_to_pin[:-4] + ".peprecpnomod", "w"), index=False)


    # Add modifications column to PEPREC file
    # the keys correspond to the mass shift for each modification
    modifications = {}
    mods = options["ms2pip"]["modifications"]
    for mod in mods:
        modifications[str(mod["mass_shift"])] = mod["name"]

    pep_list = list()
    mod_list = list()

    mass_regex = re.compile(r'(?:\[-?\d*\.?\d+]){2,}')
    dict_regex = re.compile(r'\[([^][]*)]')
    for _, row in pepfile.iterrows():
        pep = row['Peptide'].split('.', maxsplit=1)[1].rsplit('.', maxsplit=1)[0]
        mod = []

        # X!Tandem sequences can contain multiple mods per AA
        if "][" in pep:
            pep = mass_regex.sub(lambda m: "[{:g}]".format(sum([float(n) for n in m.group()[1:-1].split('][')])), pep)

        # Replace by name in dict
        if "[" in pep:
            pep_nomod = dict_regex.sub('', pep)
            m = dict_regex.search(pep)
            tmp = pep
            while m:
                mod.append("{}|{}".format(m.start(), modifications[str(float(m.group(1)))]))
                tmp = "".join([tmp[:m.start()], tmp[m.end():]])
                m = dict_regex.search(tmp)

        pep_list.append(pep_nomod)
        mod_list.append("|".join(mod))

    pepfile.loc[:, 'modifications'] = mod_list
    pepfile.loc[:, 'peptide'] = pep_list

    pepfile.to_csv("test-pyro.pepfile")

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


def tandem_pipeline(config, outname):

    xml_file = config['xtandem']['xml_file']

    logging.info("Running tandem2pin")
    # Convert .x.xml to pin. DECOY_ is the decoy pattern from MSGF+
    convert_command = "tandem2pin -P DECOY {} > {}_original.pin".format(xml_file, outname)
    subprocess.run(convert_command, shell=True)

    # PIN FILE: "Proteins" column has tab-separated values which makes the file
    # cumbersome to read. mapper.fix_pin_tabs replaces those tabs with ";"
    logging.info("Fixing tabs on pin file")
    mapper.fix_pin_tabs(outname + "_original.pin", prot_sep='|||')
    os.rename(outname + "_original_fixed.pin", outname + "_edited.pin")

    # Percolator generates its own spectrum ID, so we must match it to the mgf
    # file's TITLE field.
    logging.info("Adding mgf TITLE to pin file")
    mapper.map_mgf_title(outname + "_edited.pin", xml_file, msgs=False)
    os.rename(outname + "_edited.pin_title", outname + "_edited.pin")

    logging.info("Writing PEPREC file")
    make_pepfile(outname + "_edited.pin", config)
    os.rename(outname + "_edited.pin.peprec", outname + ".peprec")

    # Add features from the pin file to the PEPREC file
    logging.info("Adding pin features to PEPREC file")
    join_features(outname + "_edited.pin", outname + ".peprec")
    os.remove(outname + '_edited.pin')

    return outname + ".peprec", config['msgfplus']['mgf_file']


make_pepfile("/home/compomics/extra_disk/rescore-pyro-tandem/tandem/pyro-pyro.pin", options="/home/compomics/PycharmProjects/scripts/ms2rescore/config_tandem.json")