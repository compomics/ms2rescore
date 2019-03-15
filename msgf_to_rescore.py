__author__ = "Ana Silva"
__credits__ = ["Ana Silva", "Ralf Gabriels", "Sven Degroeve", "Lennart Martens"]
__license__ = "Apache License, Version 2.0"
__version__ = "0.1"
__email__ = "Ralf.Gabriels@UGent.be"

import argparse
import sys
import logging
import subprocess
import os
import re
import json
import pandas as pd

import mapper

def argument_parser():
    parser = argparse.ArgumentParser(
        description="Run MSGF+ to get PSMs and get PEPREC and pin files")
    parser.add_argument("spec_file", metavar="spectrum-file",
                        help="file containing MS2 spectra (MGF)")
    parser.add_argument("fasta_file", metavar="FASTA-file",
                        help="file containing protein sequences")
    parser.add_argument("config_file", metavar="config-file",
                        help="json file containing configurable variables")

    args = parser.parse_args()
    return(args)

def run_msgfplus(msgf_dir, mgffile, fastafile, options, log=False):
    """
    Runs MSGF+ with some fixed settings: 10ppm precursor mass tolerance,
    concatenated search, minimum peptide length of 8aa, minimum charge is 2 and
    maximum is 4, 1 match per spectrum, include additional features, no protocol
    and use 23 threads. MSGF+ is called with subprocess.

    :param msgfdir: string, path to MSGFPlus.jar
    :param mgffile: string, spectrum file (MGF,PKL,DTA,mzXML,mzDATA or mzML)
    :param fastafile: string, the database to search against in .fasta format
    :param options: dictionary, contains search engine options
    """

    if options["frag"] == 'HCD':
        m = 3
        inst = 1
    elif options["frag"] == 'CID':
        m = 1
        inst = 0
    if options["path_to_modsfile"] != '':
        mods = '-mod {} '.format(options["path_to_modsfile"])
    else:
        mods = ''
    if log:
        l = " > {}.log".format(mgffile.rstrip(".mgf"))
    else:
        l = ''


    outfile = mgffile.rstrip(".mgf") + ".mzid"

    msgf_command = "java -Xmx28000M -jar {}/MSGFPlus.jar {}-s {} -d {} -o {} -t\
         10ppm -tda 1 -m {} -inst {} -minLength {} -minCharge {} -maxCharge {}\
         -n 1 -e 1 -addFeatures 1 -protocol 0 -thread 23{}".format(msgf_dir,
         mods, mgffile, fastafile, outfile, m, inst, options["min_length"],
         options["min_charge"], options["max_charge"], l)

    logging.info("MSGF+ command: {}".format(msgf_command))
    subprocess.run(msgf_command, shell=True)

    return None

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
        if a.startswith('Charge'): charge_states.append(a)
        else: continue

    pin.loc[:, 'Charge'] = [None] * len(pin)

    for ch in charge_states:
        value = int(ch.lstrip('Charge'))
        pin.loc[pin[ch]==1, 'Charge'] = value

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
        pepfile_tosave.to_csv(path_to_pep + '.PEPREC', sep=' ', index=False)

    else:
        pepfile_tosave = pepfile.loc[pepfile.Label == 1, ['TITLE', 'modifications', 'peptide', 'Charge']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge']
        pepfile_tosave.to_csv(path_to_pep + '.targets.PEPREC', sep=' ', index=False)

        pepfile_tosave = pepfile.loc[pepfile.Label == -1, ['TITLE', 'modifications', 'peptide', 'Charge']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge']
        pepfile_tosave.to_csv(path_to_pep + '.decoys.PEPREC', sep=' ', index=False)

    return None

def join_features(path_to_pin, path_to_pep):

    pin = pd.read_csv(path_to_pin, sep='\t')
    pep = pd.read_csv(path_to_pep, sep=' ')

    pep = pep.join(pin)

    pep.to_csv(path_to_pep, sep=' ', index=False)

    return None

def main():
    # Parse arguments
    args = argument_parser()
    fname = re.sub('.mgf', '', args.spec_file, flags=re.IGNORECASE)

    # Parse config.json
    with open(args.config_file) as f:
        config = json.load(f)

    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.DEBUG
    )

    logging.info("Running MS-GF+")
    MSGF_DIR = config["search_engine_options"]["dir"]
    run_msgfplus(MSGF_DIR, args.spec_file, args.fasta_file, config["search_engine_options"])

    logging.info("Running msgf2pin")
    # Convert .mzid to pin. XXX is the decoy pattern from MSGF+
    convert_command = "msgf2pin -P XXX {}.mzid > {}.pin".format(fname, fname)
    subprocess.run(convert_command, shell=True)

    # PIN FILE: "Proteins" column has tab-separated values which makes the file
    # cumbersome to read. mapper.fix_pin_tabs replaces those tabs with ";"
    logging.info("Fixing tabs on pin file")
    mapper.fix_pin_tabs(fname + ".pin")
    # os.remove(fname + ".pin")
    os.rename(fname + "_fixed.pin", fname + ".pin")

    # Percolator generates its own spectrum ID, so we must match it to the mgf
    # file's TITLE field.
    logging.info("Adding mgf TITLE to pin file")
    mapper.map_mgf_title(fname + ".pin", fname + ".mzid", msgs=False)
    os.rename(fname + ".pin_title", fname + ".pin")

    logging.info("Writing PEPREC file")
    make_pepfile(fname + ".pin", config)
    os.rename(fname + ".pin.PEPREC", fname + ".PEPREC")

    # Add features from the pin file to the PEPREC file
    logging.info("Addin pin features to PEPREC file")
    join_features(fname + ".pin", fname + ".PEPREC")

if __name__ == '__main__':
    main()
