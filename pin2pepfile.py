"""
Create a PEPREC file, which is part of the input for ms2pip.
Takes either one (concatenated) or two <mzid> files and a <pin> file. Outputs
a <PEPREC> file, which is a list of all the identified peptides, their spectrum
ID, charge and modifications. This file is part of the input necessary to run
ms2pip.
"""

import argparse
import re
from mapper import mapper

def add_modifications(pepfile):
    """
    Take the "Peptide" column from a pepfile and remove modifications from the
    sequence adding them to a separate column in the ms2pip format.
    """
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
                modstring += str(pep.find(mod)) + '|' + modifications[mod.split(':')[1].rstrip(']')] + '|'
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

def make_pepfile(path_to_pin):
    """
    Read a pin file and create the corresponding pepfile.
    """
    pin = mapper.lazy_pin_parser(path_to_pin)

    pin.loc[:, 'Charge'] = [2 if r[1].Charge2 == '1' else 3 if r[1].Charge3 == '1' else 4 if r[1].Charge4 == '1' else 5 if r[1].Charge5 == '1' else 6 if r[1].Charge6 == '1' else 0 for r in pin.iterrows()]

    pepfile = pin[['TITLE', 'Peptide', 'Charge', 'Label']]

    # Add modifications column to PEPREC file
    pepfile = add_modifications(pepfile)

    return pepfile

def write_PEPREC(pepfile, path_to_pep, concat=True):
    """
    Write the PEPREC file
    """
    if concat:
        pepfile_tosave = pepfile[['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
        pepfile_tosave.to_csv(path_to_pep + '.PEPREC', sep=' ', index=False)

    else:
        pepfile_tosave = pepfile[pepfile.Label == 1][['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
        pepfile_tosave.to_csv(path_to_pep + '.targets.PEPREC', sep=' ', index=False)

        pepfile_tosave = pepfile[pepfile.Label == -1][['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
        pepfile_tosave.to_csv(path_to_pep + '.decoys.PEPREC', sep=' ', index=False)

    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get an ms2pip-ready PEPREC file from a pin file")
    parser.add_argument('-p', dest='pin', help='Path to pin file')
    # Don't really need the mzid files, only to know if the search was concatenated
    parser.add_argument('-m', dest='mzid', help='Path to single mzid file (concatenated search)')
    parser.add_argument('-t', dest='targets', help='Path to target search mzid file')
    parser.add_argument('-d', dest='decoys', help='Path to decoy search mzid file')


    args = parser.parse_args()

    peprec = make_pepfile(args.pin)
    write_PEPREC(peprec, args.pin)
