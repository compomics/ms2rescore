import argparse
import xmltodict
import pandas as pd
import sys
import os
import gc
from pyteomics import mzid

def get_indices(path_to_mzid):
    """
    Given a dictionary with the .mzid file, go through every
    SpectrumIdentificationResult and match the mgf TITLE to the initial part of
    the Percolator index. Save these correspondences in a dictionary.
    https://github.com/percolator/percolator/issues/147
    """
    index_map = {}
    tid = ""

    with open(path_to_mzid) as f:
    	for row in f:
    		if "<SpectrumIdentificationItem" in row:
    			l=row.rstrip().split('id=')
    			tid = l[1][1:-2] + row.split('rank=')[1][1]
    		if 'name="spectrum title"' in row:
    			l=row.rstrip().split('value=')[1].split(' ')[0]
    			#id_map[tid] = l[1:-1]
    			index_map[tid] = l[1:].rstrip('"/>')

    return index_map

def get_indices_old(path_to_mzid):
    """
    Given a dictionary with the .mzid file, go through every
    SpectrumIdentificationResult and match the mgf TITLE to the initial part of
    the Percolator index. Save these correspondences in a dictionary.
    https://github.com/percolator/percolator/issues/147
    """
    with open(path_to_mzid) as fd:
        doc = xmltodict.parse(fd.read())

    index_map = {}
    for i in tqdm(range(len(doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult']))):
        if type(doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'][i]['SpectrumIdentificationItem']) is list:
            for j in range(len(doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'][i]['SpectrumIdentificationItem'])):
                spectrum = doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'][i]
                hit = spectrum['SpectrumIdentificationItem'][j]
                # perc_id = hit['@id'] + '_' + hit['@chargeState'] + '_' + hit['@rank']
                perc_id = hit['@id']
                if type(spectrum['cvParam']) == list:
                    for k,d in enumerate(spectrum['cvParam']):
                        if 'spectrum title' in d.values(): title = spectrum['cvParam'][k]['@value']
                        else: continue
                else:
                    title = spectrum['cvParam']['@value']
                index_map[perc_id] = title
        else:
            spectrum = doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'][i]
            hit = spectrum['SpectrumIdentificationItem']
            perc_id = hit['@id']
            if type(spectrum['cvParam']) == list:
                for k,d in enumerate(spectrum['cvParam']):
                    if 'spectrum title' in d.values(): title = spectrum['cvParam'][k]['@value']
                    else: continue
            else:
                title = spectrum['cvParam']['@value']
            index_map[perc_id] = title
        del doc['MzIdentML']['DataCollection']['AnalysisData']['SpectrumIdentificationList']['SpectrumIdentificationResult'][i]
        gc.collect()

    return index_map

def get_indices_pyteomics(path_to_mzid):
    """
    Given a dictionary with the .mzid file, go through every
    SpectrumIdentificationResult and match the mgf TITLE to the initial part of
    the Percolator index. Save these correspondences in a dictionary.
    https://github.com/percolator/percolator/issues/147
    """
    index_map = {}
    for a in mzid.read(path_to_mzid):
         index_map[a['spectrumID'].split('=')[1] + '_' + str(a['SpectrumIdentificationItem'][0]['rank'])] = a['spectrum title']
    return index_map


def fix_pin_tabs(path):
    """
    Takes a pin file and re-writes it, replacing the tabs that separate the
    Proteins column with pipes
    """
    f = open(path)
    rows = f.readlines()
    outfile = path.rstrip('.pin') + '_fixed.pin'
    out = open(outfile, 'w+')

    for i, row in enumerate(rows):
        if i == 0:
            numcol = len(row.split('\t'))
            out.write(row)
        elif i == 1:
            out.write(row)
        else:
            r = row.rstrip('\n').split('\t')
            tmp = []
            for j in range(numcol):
                tmp.append(r[j])
            tmp.append(';'.join(r[numcol:]))
            out.write('\t'.join(tmp[:numcol]))
            out.write('\n')
    f.close()
    out.close()
    return None

def map_mgf_title(path_to_pin, path_to_mzid, path_to_decoy_mzid=None, msgs=False):
    """
    Add the TITLE column to the pin file. Processes the MzIdentML file (one if
    the search was concatenated, two if the target and decoy searches were ran
    separately).
    """

    # parse mzid file: xmltodict imports it as a dictionary
    # concatenated searches yield one mzid
    if not path_to_decoy_mzid:
        # Use get_indices() to get a dictionary that corresponds each percolator
        #  SpecId to its mgf TITLE
        if msgs: sys.stdout.write('\n parsing title map from mzid...\n')
        title_map = get_indices(path_to_mzid)
        gc.collect()
        # Adding mgf "TITLE" column. Avoiding using pandas due to memory issues
        if msgs: sys.stdout.write(' reading pin/writing pin_title...\n')
        pin_title = path_to_pin + '_title'
        pin_out = open(pin_title, 'w+')
        pin_in = open(path_to_pin, 'r')
        for line in pin_in:
            if line.startswith('SpecId'):
                pin_out.write(line.rstrip('\n')+'\t'+'TITLE'+'\n')
                continue
            elif line.startswith('Default'): continue
            k = line.split('\t')[0]
            k = k.split('_')[-6:-3]
            k = '_'.join(k)
            if k in title_map.keys():
                pin_out.write(line.rstrip('\n')+'\t'+title_map[k]+'\n')
            else:
                sys.stderr.write("No match found in MGF file for SpecId {}\n".format(line.split('\t')[0]))
        pin_in.close()
        pin_out.close()
        os.rename(path_to_pin + '_title', path_to_pin)

    # for separate target-decoy there are two mzid (haven't removed pandas here)
    else:
        pin = pd.read_csv(path_to_pin, header=0, skiprows=[1], sep='\t')
        pin['TITLE'] = [None] * len(pin)

        title_map_targets = get_indices(path_to_mzid)
        title_map_decoys = get_indices(path_to_decoy_mzid)

        for i in range(0, len(pin)):
            k = '_'.join(pin.loc[i, 'SpecId'].split('_')[-6:-3])
            if pin.loc[i, 'Label'] == "-1":
                if k in title_map_decoys.keys():
                    pin.loc[i, 'TITLE'] = title_map_decoys[k]
                else:
                    sys.stdout.write('oops\n')
                    continue
            elif pin.loc[i, 'Label'] == "1":
                if k in title_map_targets.keys():
                    pin.loc[i, 'TITLE'] = title_map_targets[k]
                else:
                    sys.stdout.write('oops\n')
                    continue
        pin.to_csv(path_to_pin, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Take pin file built with Percolator's msgf2pin and add the 'TITLE' from the mgf file")
    parser.add_argument('-m', dest='mzid', help='Path to single mzid file (concatenated search)')
    parser.add_argument('-t', dest='targets', help='Path to target search mzid file')
    parser.add_argument('-d', dest='decoys', help='Path to decoy search mzid file')
    parser.add_argument('-p', dest='pin', help='Path to pin file')

    args = parser.parse_args()

    # open percolator features; add column with mgf title
    sys.stdout.write('Fixing tabs on pin file... ')
    sys.stdout.flush()
    fix_pin_tabs(args.pin)
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    sys.stdout.write('Parsing pin file... ')
    sys.stdout.flush()
    pin = pd.read_csv(args.pin, header=0, sep='\t')
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    sys.stdout.write('Mapping spectrum TITLE on to pin file... \n')
    sys.stdout.flush()
    if args.decoys:
        pin = map_mgf_title(pin, args.targets, args.decoys)
    else:
        pin = map_mgf_title(pin, args.mzid)
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    sys.stdout.write('Saving pin file... ')
    sys.stdout.flush()
    pin.to_csv(args.pin, sep='\t', index=False)
    sys.stdout.write('Done! \n')
    sys.stdout.flush()
