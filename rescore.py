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
    msgf_command = "java -Xmx28000M -jar {}/MSGFPlus.jar {}-s {} -d {} -o {} -t \
        10ppm -tda 1 -m {} -inst {} -minLength 8 -minCharge 2 -maxCharge 4 -n \
        1 -addFeatures 1 -protocol 0 -thread 23 > {}.log\n".format(msgf_dir, mods,
        mgffile, fastafile, outfile + ".mzid", m, inst, outfile)

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

    charge_states = []
    for a in pin.columns:
        if a.startswith('Charge'): charge_states.append(a)
        else: continue

    pin.loc[:, 'Charge'] = [None] * len(pin)

    for ch in charge_states:
        value = int(ch.lstrip('Charge'))
        pin.loc[pin[ch]==1, 'Charge'] = value

    pepfile = pin.loc[:, ['TITLE', 'Peptide', 'Charge', 'Label']]

    # Add modifications column to PEPREC file
    # TODO: read modifications from MSGF+ modifications file OR write dict with
    # all UNIMOD modifications
    # the keys correspond to the UNIMOD keys for each modification
    modifications = {'1': 'Acetylation', '4': 'Cam', '35': 'Oxidation'}

    modlist = []
    for _, row in pepfile.iterrows():
        if 'UNIMOD' in row['Peptide']:
            pep = row['Peptide'].split('.')[1]
            mods = re.findall(r'\[([^]]*)\]', pep)
            modstring = ''
            for mod in mods:
                mod = '[' + mod + ']'
                key = mod.split(':')[1].rstrip(']')
                try:
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

    return pepfile


def write_PEPREC(pepfile, path_to_pep, concat=True):
    """
    Write the PEPREC file, which will be the input to MS2PIP.

    :param pepfile: pd.DataFrame created by rescore.make_pepfile()
    :param path_to_pep: string, path where to save the PEPREC file(s)
    :param concat: boolean, True if the search was concatenated
    """

    if concat:
        pepfile_tosave = pepfile.loc[:, ['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
        pepfile_tosave.to_csv(path_to_pep + '.PEPREC', sep=' ', index=False)

    else:
        pepfile_tosave = pepfile.loc[pepfile.Label == 1, ['TITLE', 'modifications', 'peptide', 'Charge', 'Label']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
        pepfile_tosave.to_csv(path_to_pep + '.targets.PEPREC', sep=' ', index=False)

        pepfile_tosave = pepfile.loc[pepfile.Label == -1, ['TITLE', 'modifications','peptide', 'Charge', 'Label']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label']
        pepfile_tosave.to_csv(path_to_pep + '.decoys.PEPREC', sep=' ', index=False)

    return None

def calculate_features(path_to_pred_and_emp):
    rescore_features = pd.DataFrame(columns=['spec_id', 'peplen', 'charge',
        'spec_pearson_norm', 'ionb_pearson_norm', 'iony_pearson_norm',
        'spec_spearman_norm', 'ionb_spearman_norm', 'iony_spearman_norm',
        'spec_mse_norm', 'ionb_mse_norm', 'iony_mse_norm', 'min_abs_diff_iontype_norm',
        'max_abs_diff_iontype_norm', 'min_abs_diff_norm', 'max_abs_diff_norm',
        'abs_diff_Q1_norm', 'abs_diff_Q2_norm', 'abs_diff_Q3_norm', 'mean_abs_diff_norm',
        'std_abs_diff_norm', 'ionb_min_abs_diff_norm', 'ionb_max_abs_diff_norm',
        'ionb_abs_diff_Q1_norm', 'ionb_abs_diff_Q2_norm', 'ionb_abs_diff_Q3_norm',
        'ionb_mean_abs_diff_norm', 'ionb_std_abs_diff_norm', 'iony_min_abs_diff_norm',
        'iony_max_abs_diff_norm', 'iony_abs_diff_Q1_norm', 'iony_abs_diff_Q2_norm',
        'iony_abs_diff_Q3_norm', 'iony_mean_abs_diff_norm', 'iony_std_abs_diff_norm',
        'dotprod_norm', 'dotprod_ionb_norm', 'dotprod_iony_norm', 'cos_norm',
        'cos_ionb_norm', 'cos_iony_norm', 'spec_pearson', 'ionb_pearson', 'iony_pearson',
        'spec_spearman', 'ionb_spearman', 'iony_spearman', 'spec_mse', 'ionb_mse',
        'iony_mse', 'min_abs_diff_iontype', 'max_abs_diff_iontype', 'min_abs_diff',
        'max_abs_diff', 'abs_diff_Q1', 'abs_diff_Q2', 'abs_diff_Q3', 'mean_abs_diff',
        'std_abs_diff', 'ionb_min_abs_diff', 'ionb_max_abs_diff', 'ionb_abs_diff_Q1',
        'ionb_abs_diff_Q2', 'ionb_abs_diff_Q3', 'ionb_mean_abs_diff', 'ionb_std_abs_diff',
        'iony_min_abs_diff', 'iony_max_abs_diff', 'iony_abs_diff_Q1', 'iony_abs_diff_Q2',
        'iony_abs_diff_Q3', 'iony_mean_abs_diff', 'iony_std_abs_diff', 'dotprod', 'dotprod_ionb',
        'dotprod_iony', 'cos', 'cos_ionb', 'cos_iony'])

    df = pd.read_csv(path_to_pred_and_emp)

    for peptide in df.spec_id.unique():
        tmp = df[df.spec_id == peptide]
        tmp.loc[tmp.prediction < np.log2(0.001), 'prediction'] = np.log2(0.001)
        tmp["abs_diff"] = np.abs(tmp["target"] - tmp["prediction"])

        feats = {}

        feats["spec_id"] = tmp["spec_id"][0]
        feats["peplen"] = tmp["peplen"][0]
        feats["charge"] = tmp["charge"][0]

        # calculation of features between normalized spectra
        feats["spec_pearson_norm"] = tmp["target"].corr(tmp["prediction"])
        feats["ionb_pearson_norm"] = tmp[tmp.ion == "b"]["target"].corr(tmp[tmp.ion == "b"]["prediction"])
        feats["iony_pearson_norm"] = tmp[tmp.ion == "y"]["target"].corr(tmp[tmp.ion == "y"]["prediction"])

        feats["spec_spearman_norm"] = tmp["target"].corr(tmp["prediction"], "spearman")
        feats["ionb_spearman_norm"] = tmp[tmp.ion == "b"]["target"].corr(tmp[tmp.ion == "b"]["prediction"], "spearman")
        feats["iony_spearman_norm"] = tmp[tmp.ion == "y"]["target"].corr(tmp[tmp.ion == "y"]["prediction"], "spearman")

        feats["spec_mse_norm"] = mean_squared_error(tmp["target"], tmp["prediction"])
        feats["ionb_mse_norm"] = mean_squared_error(tmp[tmp.ion == "b"]["target"], tmp[tmp.ion == "b"]["prediction"])
        feats["iony_mse_norm"] = mean_squared_error(tmp[tmp.ion == "y"]["target"], tmp[tmp.ion == "y"]["prediction"])

        feats["min_abs_diff_iontype_norm"] = tmp[tmp.abs_diff == np.min(tmp["abs_diff"])]["ion"].values[0]
        feats["max_abs_diff_iontype_norm"] = tmp[tmp.abs_diff == np.max(tmp["abs_diff"])]["ion"].values[0]

        feats["min_abs_diff_norm"] = np.min(tmp["abs_diff"])
        feats["max_abs_diff_norm"] = np.max(tmp["abs_diff"])
        feats["abs_diff_Q1_norm"] = tmp.quantile(q=0.25)["abs_diff"]
        feats["abs_diff_Q2_norm"] = tmp.quantile(q=0.5)["abs_diff"]
        feats["abs_diff_Q3_norm"] = tmp.quantile(q=0.75)["abs_diff"]
        feats["mean_abs_diff_norm"] = np.mean(tmp["abs_diff"])
        feats["std_abs_diff_norm"] = np.std(tmp["abs_diff"])

        feats["ionb_min_abs_diff_norm"] = np.min(tmp[tmp.ion == "b"]["abs_diff"])
        feats["ionb_max_abs_diff_norm"] = np.max(tmp[tmp.ion == "b"]["abs_diff"])
        feats["ionb_abs_diff_Q1_norm"] = tmp[tmp.ion == "b"].quantile(q=0.25)["abs_diff"]
        feats["ionb_abs_diff_Q2_norm"] = tmp[tmp.ion == "b"].quantile(q=0.5)["abs_diff"]
        feats["ionb_abs_diff_Q3_norm"] = tmp[tmp.ion == "b"].quantile(q=0.75)["abs_diff"]
        feats["ionb_mean_abs_diff_norm"] = np.mean(tmp[tmp.ion == "b"]["abs_diff"])
        feats["ionb_std_abs_diff_norm"] = np.std(tmp[tmp.ion == "b"]["abs_diff"])

        feats["iony_min_abs_diff_norm"] = np.min(tmp[tmp.ion == "y"]["abs_diff"])
        feats["iony_max_abs_diff_norm"] = np.max(tmp[tmp.ion == "y"]["abs_diff"])
        feats["iony_abs_diff_Q1_norm"] = tmp[tmp.ion == "y"].quantile(q=0.25)["abs_diff"]
        feats["iony_abs_diff_Q2_norm"] = tmp[tmp.ion == "y"].quantile(q=0.5)["abs_diff"]
        feats["iony_abs_diff_Q3_norm"] = tmp[tmp.ion == "y"].quantile(q=0.75)["abs_diff"]
        feats["iony_mean_abs_diff_norm"] = np.mean(tmp[tmp.ion == "y"]["abs_diff"])
        feats["iony_std_abs_diff_norm"] = np.std(tmp[tmp.ion == "y"]["abs_diff"])

        feats["dotprod_norm"] = np.dot(tmp["target"], tmp["prediction"])
        feats["dotprod_ionb_norm"] = np.dot(tmp[tmp.ion == "b"]["target"], tmp[tmp.ion == "b"]["prediction"])
        feats["dotprod_iony_norm"] = np.dot(tmp[tmp.ion == "y"]["target"], tmp[tmp.ion == "y"]["prediction"])

        feats["cos_norm"] = feats["dotprod_norm"] / (np.linalg.norm(tmp["target"], 2) * np.linalg.norm(tmp["prediction"], 2))
        feats["cos_ionb_norm"] = feats["dotprod_ionb_norm"] / (np.linalg.norm(tmp[tmp.ion == "b"]["target"], 2) * np.linalg.norm(tmp[tmp.ion == "b"]["prediction"], 2))
        feats["cos_iony_norm"] = feats["dotprod_iony_norm"] / (np.linalg.norm(tmp[tmp.ion == "y"]["target"], 2) * np.linalg.norm(tmp[tmp.ion == "y"]["prediction"], 2))

        # same features but between un-normalized spectral
        tmp['target'] = 2**tmp['target']-0.001
        tmp['prediction'] = 2**tmp['prediction']-0.001
        tmp["abs_diff"] = np.abs(tmp["target"] - tmp["prediction"])

        feats["spec_pearson"] = tmp["target"].corr(tmp["prediction"])
        feats["ionb_pearson"] = tmp[tmp.ion == "b"]["target"].corr(tmp[tmp.ion == "b"]["prediction"])
        feats["iony_pearson"] = tmp[tmp.ion == "y"]["target"].corr(tmp[tmp.ion == "y"]["prediction"])

        feats["spec_spearman"] = tmp["target"].corr(tmp["prediction"], "spearman")
        feats["ionb_spearman"] = tmp[tmp.ion == "b"]["target"].corr(tmp[tmp.ion == "b"]["prediction"], "spearman")
        feats["iony_spearman"] = tmp[tmp.ion == "y"]["target"].corr(tmp[tmp.ion == "y"]["prediction"], "spearman")

        feats["spec_mse"] = mean_squared_error(tmp["target"], tmp["prediction"])
        feats["ionb_mse"] = mean_squared_error(tmp[tmp.ion == "b"]["target"], tmp[tmp.ion == "b"]["prediction"])
        feats["iony_mse"] = mean_squared_error(tmp[tmp.ion == "y"]["target"], tmp[tmp.ion == "y"]["prediction"])

        feats["min_abs_diff_iontype"] = tmp[tmp.abs_diff == np.min(tmp["abs_diff"])]["ion"].values[0]
        feats["max_abs_diff_iontype"] = tmp[tmp.abs_diff == np.max(tmp["abs_diff"])]["ion"].values[0]

        feats["min_abs_diff"] = np.min(tmp["abs_diff"])
        feats["max_abs_diff"] = np.max(tmp["abs_diff"])
        feats["abs_diff_Q1"] = tmp.quantile(q=0.25)["abs_diff"]
        feats["abs_diff_Q2"] = tmp.quantile(q=0.5)["abs_diff"]
        feats["abs_diff_Q3"] = tmp.quantile(q=0.75)["abs_diff"]
        feats["mean_abs_diff"] = np.mean(tmp["abs_diff"])
        feats["std_abs_diff"] = np.std(tmp["abs_diff"])

        feats["ionb_min_abs_diff"] = np.min(tmp[tmp.ion == "b"]["abs_diff"])
        feats["ionb_max_abs_diff"] = np.max(tmp[tmp.ion == "b"]["abs_diff"])
        feats["ionb_abs_diff_Q1"] = tmp[tmp.ion == "b"].quantile(q=0.25)["abs_diff"]
        feats["ionb_abs_diff_Q2"] = tmp[tmp.ion == "b"].quantile(q=0.5)["abs_diff"]
        feats["ionb_abs_diff_Q3"] = tmp[tmp.ion == "b"].quantile(q=0.75)["abs_diff"]
        feats["ionb_mean_abs_diff"] = np.mean(tmp[tmp.ion == "b"]["abs_diff"])
        feats["ionb_std_abs_diff"] = np.std(tmp[tmp.ion == "b"]["abs_diff"])

        feats["iony_min_abs_diff"] = np.min(tmp[tmp.ion == "y"]["abs_diff"])
        feats["iony_max_abs_diff"] = np.max(tmp[tmp.ion == "y"]["abs_diff"])
        feats["iony_abs_diff_Q1"] = tmp[tmp.ion == "y"].quantile(q=0.25)["abs_diff"]
        feats["iony_abs_diff_Q2"] = tmp[tmp.ion == "y"].quantile(q=0.5)["abs_diff"]
        feats["iony_abs_diff_Q3"] = tmp[tmp.ion == "y"].quantile(q=0.75)["abs_diff"]
        feats["iony_mean_abs_diff"] = np.mean(tmp[tmp.ion == "y"]["abs_diff"])
        feats["iony_std_abs_diff"] = np.std(tmp[tmp.ion == "y"]["abs_diff"])

        feats["dotprod"] = np.dot(tmp["target"], tmp["prediction"])
        feats["dotprod_ionb"] = np.dot(tmp[tmp.ion == "b"]["target"], tmp[tmp.ion == "b"]["prediction"])
        feats["dotprod_iony"] = np.dot(tmp[tmp.ion == "y"]["target"], tmp[tmp.ion == "y"]["prediction"])

        feats["cos"] = feats["dotprod"] / (np.linalg.norm(tmp["target"], 2) * np.linalg.norm(tmp["prediction"], 2))
        feats["cos_ionb"] = feats["dotprod_ionb"] / (np.linalg.norm(tmp[tmp.ion == "b"]["target"], 2) * np.linalg.norm(tmp[tmp.ion == "b"]["prediction"], 2))
        feats["cos_iony"] = feats["dotprod_iony"] / (np.linalg.norm(tmp[tmp.ion == "y"]["target"], 2) * np.linalg.norm(tmp[tmp.ion == "y"]["prediction"], 2))

        rescore_features = rescore_features.append(feats, ignore_index=True)


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
