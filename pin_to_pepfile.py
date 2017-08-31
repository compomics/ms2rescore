import argparse
import xmltodict
import pandas as pd
import sys
import re

# Arguments: percolator features, .mzid
parser = argparse.ArgumentParser(description="Get an ms2pip-ready PEPREC file from a pin file")
parser.add_argument('-m', dest='mzid', help='Path to single mzid file (concatenated search)')
parser.add_argument('-t', dest='targets', help='Path to target search mzid file')
parser.add_argument('-d', dest='decoys', help='Path to decoy search mzid file')
parser.add_argument('-p', dest='pin', help='Path to pin file')

args = parser.parse_args()

def get_indices(doc):
    """
    Given a dictionary with the .mzid file, go through every
    SpectrumIdentificationResult and match the mgf TITLE to the initial part of
    the Percolator index. Save these correspondences in a dictionary
    """
    mapper = {}
    for i in range(len(doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'])):
        if type(doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'][i]['SpectrumIdentificationItem']) is list:
            for j in range(len(doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'][i]['SpectrumIdentificationItem'])):
                spectrum = doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'][i]
                hit = spectrum['SpectrumIdentificationItem'][j]
                # perc_id = hit['@id'] + '_' + hit['@chargeState'] + '_' + hit['@rank']
                perc_id = hit['@id']
                title = spectrum['cvParam']['@value']
                mapper[perc_id] = title
        else:
            spectrum = doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'][i]
            hit = spectrum['SpectrumIdentificationItem']
            # perc_id = hit['@id'] + '_' + hit['@chargeState'] + '_' + hit['@rank']
            perc_id = hit['@id']
            title = spectrum['cvParam']['@value']
            mapper[perc_id] = title

    return mapper

def lazy_pin_parser(path):
    """
    To parse the pin file. In some rows, the "Proteins" column contains
    tab-separated values. Because of this normal parsers don't work too well.
    This is a lazily built parser that addresses this.
    """
    f = open(path)
    rows = f.readlines()
    for i, row in enumerate(rows):
        if i == 0:
            data = pd.DataFrame(columns=[r.replace('"', '').rstrip('\n') for r in row.split('\t')])
            n_rows = len(row.split('\t'))
        elif i == 1: continue # row 1 is initial direction
        else:
            r = row.split('\t')
            tmp = []
            for j in range(n_rows-1):
                tmp.append(r[j])
            tmp.append('|'.join(r[n_rows-1:]).rstrip('\n'))
            data.loc[i-1] = tmp
    return data

# open percolator features; add column with mgf title
sys.stdout.write('Parsing pin file... '); sys.stdout.flush()
pin = lazy_pin_parser(args.pin)
pin['TITLE'] = [None] * len(pin)
sys.stdout.write('Done! \n')
# parse mzid file: xmltodict imports it as a dictionary
# concatenated searches yield one mzid
if args.mzid:
    sys.stdout.write('Concatenated search results; parsing .mzid... '); sys.stdout.flush()
    with open(args.mzid) as fd:
         doc = xmltodict.parse(fd.read())

    # Adding mgf "TITLE" column
    mapper = get_indices(doc)
    sys.stdout.write('Done! \n')

    sys.stdout.write('Adding "TITLE" to pin file... '); sys.stdout.flush()
    for i in range(1, len(pin)+1): # because index starts at 1
        k = '_'.join(pin.loc[i, 'SpecId'].split('_')[-6:-3])
        if k in mapper.keys():
            pin.loc[i, 'TITLE'] = mapper[k]
        else:
            continue
    sys.stdout.write('Done! \n')

# for separate target-decoy there are two mzid
elif args.targets and args.decoys:
    sys.stdout.write('Separate target and decoy search; \n')

    sys.stdout.write('parsing targets .mzid... '); sys.stdout.flush()
    with open(args.targets) as fd:
         doc = xmltodict.parse(fd.read())

    mapper_targets = get_indices(doc)
    sys.stdout.write('Done! \n')

    sys.stdout.write('parsing decoys .mzid... '); sys.stdout.flush()
    with open(args.decoys) as fd:
         doc = xmltodict.parse(fd.read())

    mapper_decoys = get_indices(doc)
    sys.stdout.write('Done! \n')

    sys.stdout.write('Adding "TITLE" to pin file... '); sys.stdout.flush()
    for i in range(1, len(pin)+1): # because index starts at 1
        k = '_'.join(pin.loc[i, 'SpecId'].split('_')[-6:-3])
        if pin.loc[i, 'Label'] == "-1":
            if k in mapper_decoys.keys():
                pin.loc[i, 'TITLE'] = mapper_decoys[k]
            else:
                sys.stdout.write('oops\n')
                continue
        elif pin.loc[i, 'Label'] == "1":
            if k in mapper_targets.keys():
                pin.loc[i, 'TITLE'] = mapper_targets[k]
            else:
                sys.stdout.write('oops\n')
                continue

    sys.stdout.write('Done! \n')


sys.stdout.write('Creating PEPREC file... '); sys.stdout.flush()
pin.loc[:,'Charge'] = [2 if r[1].Charge2=='1' else 3 if r[1].Charge3=='1' else 4 if r[1].Charge4=='1'
                 else 5 if r[1].Charge5=='1' else 6 if r[1].Charge6=='1' else 0 for r in pin.iterrows()]

pepfile = pin[['TITLE', 'Peptide', 'Charge', 'Label']]

# Add modifications column to PEPREC file
# TODO: read modifications from MSGF+ modifications file
# the keys correspond to the UNIMOD keys for each modification
modifications = {'4': 'Cam', '35': 'Oxidation'}

modlist = []
for _, r in pepfile.iterrows():
# for _, r in pepfile[pepfile.Peptide == 'R.EYWAPGHAAC[UNIMOD:4]AGC[UNIMOD:4]GC[UNIMOD:4]ATALR.L'].iterrows():
    if 'UNIMOD' in r['Peptide']:
        pep = r['Peptide'].split('.')[1]
        mods = re.findall(r'\[([^]]*)\]', pep)
        modstring = ''
        for mod in mods:
            mod ='[' + mod + ']'
            modstring += str(pep.find(mod)) + '|' + modifications[mod.split(':')[1].rstrip(']')] + '|'
            # print(str(pep.find(mod)-1) + '|' + modifications[mod.split(':')[1].rstrip(']')])
            pep = pep.replace(mod, '', 1)
        modlist.append(modstring.rstrip('|'))
    else:
        modlist.append('')

pepfile.loc[:,'modifications'] = modlist

# Rewrite peptide sequences without the UNIMOD modification ids
pep_list = []
for _, r in pepfile.iterrows():
    pep = r['Peptide']
    pep = pep.split('.')[1]
    if 'UNIMOD' in pep:
        mods = re.findall(r'\[([^]]*)\]', pep)
        for m in mods:
            pep = pep.replace('[' + m + ']', '', 1)
    pep_list.append(pep)

pepfile.loc[:,'peptide'] = pep_list
sys.stdout.write('Done! \n')

sys.stdout.write('Writing .PEPREC file... '); sys.stdout.flush()
if args.mzid:
    pepfile_tosave = pepfile[['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
    pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
    pepfile_tosave.to_csv(args.mzid + '.PEPREC', sep=' ', index=False)

elif args.targets and args.decoys:
    pepfile_tosave = pepfile[pepfile.Label == 1][['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
    pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
    pepfile_tosave.to_csv(args.targets + '.PEPREC', sep=' ', index=False)

    pepfile_tosave = pepfile[pepfile.Label == -1][['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
    pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
    pepfile_tosave.to_csv(args.decoys + '.PEPREC', sep=' ', index=False)
sys.stdout.write('Done! \n')
