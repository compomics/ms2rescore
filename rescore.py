"""
Functions necessary to run the rescore algorithm. Currently supports MSGF+ with
concatenated searches.
"""

import subprocess
import multiprocessing
import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler

def run_msgfplus(msgf_dir, mgffile, fastafile, options):
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

    outfile = mgffile.rstrip(".mgf") + ".mzid"

    msgf_command = "java -Xmx28000M -jar {}/MSGFPlus.jar {}-s {} -d {} -o {} -t\
         10ppm -tda 1 -m {} -inst {} -minLength {} -minCharge {} -maxCharge {}\
         -n 1 -e 1 -addFeatures 1 -protocol 0 -thread 23 > {}.log".format(msgf_dir,
         mods, mgffile, fastafile, outfile, m, inst, options["min_length"],
         options["min_charge"], options["max_charge"], mgffile.rstrip(".mgf"))

    sys.stdout.write("Running search with MSGF+: {}".format(msgf_command))
    sys.stdout.flush()
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
        pepfile_tosave = pepfile.loc[:, ['TITLE', 'modifications', 'peptide', 'Charge', 'Label', 'Proteins']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label', 'Proteins']
        pepfile_tosave.to_csv(path_to_pep + '.PEPREC', sep=' ', index=False)

    else:
        pepfile_tosave = pepfile.loc[pepfile.Label == 1, ['TITLE', 'modifications', 'peptide', 'Charge', 'Label', 'Proteins']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label', 'Proteins']
        pepfile_tosave.to_csv(path_to_pep + '.targets.PEPREC', sep=' ', index=False)

        pepfile_tosave = pepfile.loc[pepfile.Label == -1, ['TITLE', 'modifications', 'peptide', 'Charge', 'Label', 'Proteins']]
        pepfile_tosave.columns = ['spec_id', 'modifications', 'peptide', 'charge', 'Label', 'Proteins']
        pepfile_tosave.to_csv(path_to_pep + '.decoys.PEPREC', sep=' ', index=False)

    return None

def make_ms2pip_config(options):
    """
    write configuration file for ms2pip based on what's on the rescore config
    file.
    """
    cwd = os.getcwd()
    ms2pip_config = open(cwd + "/rescore_config.txt", 'wt')

    if options["ms2pip"]["frag"] == "CID":
        ms2pip_config.write("frag_method=CID\n")
        ms2pip_config.write("frag_error=0.8\n")
    elif options["ms2pip"]["frag"] == "phospho":
        ms2pip_config.write("frag_method=phospho\n")
        ms2pip_config.write("frag_error=0.02\n")
    else:
        ms2pip_config.write("frag_method=HCD\n")
        ms2pip_config.write("frag_error=0.02\n")

    ms2pip_config.write("\n")

    modifications = options["ms2pip"]["modifications"]
    for mod in modifications:
        if mod["amino_acid"] == None and mod["n_term"] == True:
            aa = "N-term"
        else:
            aa = mod["amino_acid"]
        tmp = ','.join([mod["name"], str(mod["mass_shift"]), "opt", aa])
        ms2pip_config.write("ptm=" + tmp + "\n")

    ms2pip_config.close()

def compute_features(df):
    conv = {}
    conv['b'] = 0
    conv['y'] = 1
    rescore_features = pd.DataFrame(columns=['spec_id', 'charge',
        'spec_pearson_norm', 'ionb_pearson_norm', 'iony_pearson_norm',
        'spec_spearman_norm', 'ionb_spearman_norm', 'iony_spearman_norm',
        'spec_mse_norm', 'ionb_mse_norm', 'iony_mse_norm',
        'min_abs_diff_iontype_norm', 'max_abs_diff_iontype_norm',
        'min_abs_diff_norm', 'max_abs_diff_norm', 'abs_diff_Q1_norm',
        'abs_diff_Q2_norm', 'abs_diff_Q3_norm', 'mean_abs_diff_norm',
        'std_abs_diff_norm', 'ionb_min_abs_diff_norm', 'ionb_max_abs_diff_norm',
        'ionb_abs_diff_Q1_norm', 'ionb_abs_diff_Q2_norm',
        'ionb_abs_diff_Q3_norm', 'ionb_mean_abs_diff_norm',
        'ionb_std_abs_diff_norm', 'iony_min_abs_diff_norm',
        'iony_max_abs_diff_norm', 'iony_abs_diff_Q1_norm',
        'iony_abs_diff_Q2_norm', 'iony_abs_diff_Q3_norm',
        'iony_mean_abs_diff_norm', 'iony_std_abs_diff_norm', 'dotprod_norm',
        'dotprod_ionb_norm', 'dotprod_iony_norm', 'cos_norm', 'cos_ionb_norm',
        'cos_iony_norm', 'spec_pearson', 'ionb_pearson', 'iony_pearson',
        'spec_spearman', 'ionb_spearman', 'iony_spearman', 'spec_mse',
        'ionb_mse', 'iony_mse', 'min_abs_diff_iontype', 'max_abs_diff_iontype',
        'min_abs_diff', 'max_abs_diff', 'abs_diff_Q1', 'abs_diff_Q2',
        'abs_diff_Q3', 'mean_abs_diff', 'std_abs_diff', 'ionb_min_abs_diff',
        'ionb_max_abs_diff', 'ionb_abs_diff_Q1', 'ionb_abs_diff_Q2',
        'ionb_abs_diff_Q3', 'ionb_mean_abs_diff', 'ionb_std_abs_diff',
        'iony_min_abs_diff', 'iony_max_abs_diff', 'iony_abs_diff_Q1',
        'iony_abs_diff_Q2', 'iony_abs_diff_Q3', 'iony_mean_abs_diff',
        'iony_std_abs_diff', 'dotprod', 'dotprod_ionb', 'dotprod_iony', 'cos',
        'cos_ionb', 'cos_iony'])

    for peptide in df.spec_id.unique():
        tmp = df[df.spec_id == peptide].copy()
        tmp.loc[tmp.prediction < np.log2(0.001), 'prediction'] = np.log2(0.001)
        tmp.loc[:, "abs_diff"] = np.abs(tmp["target"] - tmp["prediction"])

        feats = {}

        feats["spec_id"] = tmp["spec_id"].unique()[0]
        feats["charge"] = tmp["charge"].unique()[0]

        # calculation of features between normalized spectra
        feats["spec_pearson_norm"] = tmp["target"].corr(tmp["prediction"])
        feats["ionb_pearson_norm"] = tmp.loc[tmp.ion == "b", "target"].corr(tmp.loc[tmp.ion == "b", "prediction"])
        feats["iony_pearson_norm"] = tmp.loc[tmp.ion == "y", "target"].corr(tmp.loc[tmp.ion == "y", "prediction"])

        feats["spec_spearman_norm"] = tmp["target"].corr(tmp["prediction"], "spearman")
        feats["ionb_spearman_norm"] = tmp.loc[tmp.ion == "b", "target"].corr(tmp.loc[tmp.ion == "b", "prediction"], "spearman")
        feats["iony_spearman_norm"] = tmp.loc[tmp.ion == "y", "target"].corr(tmp.loc[tmp.ion == "y", "prediction"], "spearman")

        feats["spec_mse_norm"] = mean_squared_error(tmp["target"], tmp["prediction"])
        feats["ionb_mse_norm"] = mean_squared_error(tmp.loc[tmp.ion == "b", "target"], tmp.loc[tmp.ion == "b", "prediction"])
        feats["iony_mse_norm"] = mean_squared_error(tmp.loc[tmp.ion == "y", "target"], tmp.loc[tmp.ion == "y", "prediction"])

        feats["min_abs_diff_iontype_norm"] = conv[tmp[tmp.abs_diff == np.min(tmp["abs_diff"])]["ion"].values[0]]
        feats["max_abs_diff_iontype_norm"] = conv[tmp[tmp.abs_diff == np.max(tmp["abs_diff"])]["ion"].values[0]]

        feats["min_abs_diff_norm"] = np.min(tmp["abs_diff"])
        feats["max_abs_diff_norm"] = np.max(tmp["abs_diff"])
        feats["abs_diff_Q1_norm"] = tmp.quantile(q=0.25)["abs_diff"]
        feats["abs_diff_Q2_norm"] = tmp.quantile(q=0.5)["abs_diff"]
        feats["abs_diff_Q3_norm"] = tmp.quantile(q=0.75)["abs_diff"]
        feats["mean_abs_diff_norm"] = np.mean(tmp["abs_diff"])
        feats["std_abs_diff_norm"] = np.std(tmp["abs_diff"])

        feats["ionb_min_abs_diff_norm"] = np.min(tmp.loc[tmp.ion == "b", "abs_diff"])
        feats["ionb_max_abs_diff_norm"] = np.max(tmp.loc[tmp.ion == "b", "abs_diff"])
        feats["ionb_abs_diff_Q1_norm"] = tmp[tmp.ion == "b"].quantile(q=0.25)["abs_diff"]
        feats["ionb_abs_diff_Q2_norm"] = tmp[tmp.ion == "b"].quantile(q=0.5)["abs_diff"]
        feats["ionb_abs_diff_Q3_norm"] = tmp[tmp.ion == "b"].quantile(q=0.75)["abs_diff"]
        feats["ionb_mean_abs_diff_norm"] = np.mean(tmp.loc[tmp.ion == "b", "abs_diff"])
        feats["ionb_std_abs_diff_norm"] = np.std(tmp.loc[tmp.ion == "b", "abs_diff"])

        feats["iony_min_abs_diff_norm"] = np.min(tmp.loc[tmp.ion == "y", "abs_diff"])
        feats["iony_max_abs_diff_norm"] = np.max(tmp.loc[tmp.ion == "y", "abs_diff"])
        feats["iony_abs_diff_Q1_norm"] = tmp[tmp.ion == "y"].quantile(q=0.25)["abs_diff"]
        feats["iony_abs_diff_Q2_norm"] = tmp[tmp.ion == "y"].quantile(q=0.5)["abs_diff"]
        feats["iony_abs_diff_Q3_norm"] = tmp[tmp.ion == "y"].quantile(q=0.75)["abs_diff"]
        feats["iony_mean_abs_diff_norm"] = np.mean(tmp.loc[tmp.ion == "y", "abs_diff"])
        feats["iony_std_abs_diff_norm"] = np.std(tmp.loc[tmp.ion == "y", "abs_diff"])

        feats["dotprod_norm"] = np.dot(tmp["target"], tmp["prediction"])
        feats["dotprod_ionb_norm"] = np.dot(tmp.loc[tmp.ion == "b"]["target"], tmp[tmp.ion == "b"]["prediction"])
        feats["dotprod_iony_norm"] = np.dot(tmp[tmp.ion == "y"]["target"], tmp[tmp.ion == "y"]["prediction"])

        feats["cos_norm"] = feats["dotprod_norm"] / (np.linalg.norm(tmp["target"], 2) * np.linalg.norm(tmp["prediction"], 2))
        feats["cos_ionb_norm"] = feats["dotprod_ionb_norm"] / (np.linalg.norm(tmp[tmp.ion == "b"]["target"], 2) * np.linalg.norm(tmp[tmp.ion == "b"]["prediction"], 2))
        feats["cos_iony_norm"] = feats["dotprod_iony_norm"] / (np.linalg.norm(tmp[tmp.ion == "y"]["target"], 2) * np.linalg.norm(tmp[tmp.ion == "y"]["prediction"], 2))

        # same features but between un-normalized spectral
        tmp.loc[:, 'target'] = 2**tmp['target']-0.001
        tmp.loc[:, 'prediction'] = 2**tmp['prediction']-0.001
        tmp.loc[:, "abs_diff"] = np.abs(tmp["target"] - tmp["prediction"])

        feats["spec_pearson"] = tmp["target"].corr(tmp["prediction"])
        feats["ionb_pearson"] = tmp.loc[tmp.ion == "b", "target"].corr(tmp.loc[tmp.ion == "b", "prediction"])
        feats["iony_pearson"] = tmp.loc[tmp.ion == "y", "target"].corr(tmp.loc[tmp.ion == "y", "prediction"])

        feats["spec_spearman"] = tmp["target"].corr(tmp["prediction"], "spearman")
        feats["ionb_spearman"] = tmp.loc[tmp.ion == "b", "target"].corr(tmp.loc[tmp.ion == "b", "prediction"], "spearman")
        feats["iony_spearman"] = tmp.loc[tmp.ion == "y", "target"].corr(tmp.loc[tmp.ion == "y", "prediction"], "spearman")

        feats["spec_mse"] = mean_squared_error(tmp["target"], tmp["prediction"])
        feats["ionb_mse"] = mean_squared_error(tmp.loc[tmp.ion == "b", "target"], tmp.loc[tmp.ion == "b", "prediction"])
        feats["iony_mse"] = mean_squared_error(tmp.loc[tmp.ion == "y", "target"], tmp.loc[tmp.ion == "y", "prediction"])

        feats["min_abs_diff_iontype"] = conv[tmp[tmp.abs_diff == np.min(tmp["abs_diff"])]["ion"].values[0]]
        feats["max_abs_diff_iontype"] = conv[tmp[tmp.abs_diff == np.max(tmp["abs_diff"])]["ion"].values[0]]

        feats["min_abs_diff"] = np.min(tmp["abs_diff"])
        feats["max_abs_diff"] = np.max(tmp["abs_diff"])
        feats["abs_diff_Q1"] = tmp.quantile(q=0.25)["abs_diff"]
        feats["abs_diff_Q2"] = tmp.quantile(q=0.5)["abs_diff"]
        feats["abs_diff_Q3"] = tmp.quantile(q=0.75)["abs_diff"]
        feats["mean_abs_diff"] = np.mean(tmp["abs_diff"])
        feats["std_abs_diff"] = np.std(tmp["abs_diff"])

        feats["ionb_min_abs_diff"] = np.min(tmp.loc[tmp.ion == "b", "abs_diff"])
        feats["ionb_max_abs_diff"] = np.max(tmp.loc[tmp.ion == "b", "abs_diff"])
        feats["ionb_abs_diff_Q1"] = tmp[tmp.ion == "b"].quantile(q=0.25)["abs_diff"]
        feats["ionb_abs_diff_Q2"] = tmp[tmp.ion == "b"].quantile(q=0.5)["abs_diff"]
        feats["ionb_abs_diff_Q3"] = tmp[tmp.ion == "b"].quantile(q=0.75)["abs_diff"]
        feats["ionb_mean_abs_diff"] = np.mean(tmp.loc[tmp.ion == "b", "abs_diff"])
        feats["ionb_std_abs_diff"] = np.std(tmp.loc[tmp.ion == "b", "abs_diff"])

        feats["iony_min_abs_diff"] = np.min(tmp.loc[tmp.ion == "y", "abs_diff"])
        feats["iony_max_abs_diff"] = np.max(tmp.loc[tmp.ion == "y", "abs_diff"])
        feats["iony_abs_diff_Q1"] = tmp[tmp.ion == "y"].quantile(q=0.25)["abs_diff"]
        feats["iony_abs_diff_Q2"] = tmp[tmp.ion == "y"].quantile(q=0.5)["abs_diff"]
        feats["iony_abs_diff_Q3"] = tmp[tmp.ion == "y"].quantile(q=0.75)["abs_diff"]
        feats["iony_mean_abs_diff"] = np.mean(tmp.loc[tmp.ion == "y", "abs_diff"])
        feats["iony_std_abs_diff"] = np.std(tmp.loc[tmp.ion == "y", "abs_diff"])

        feats["dotprod"] = np.dot(tmp["target"], tmp["prediction"])
        feats["dotprod_ionb"] = np.dot(tmp.loc[tmp.ion == "b", "target"], tmp.loc[tmp.ion == "b", "prediction"])
        feats["dotprod_iony"] = np.dot(tmp.loc[tmp.ion == "y", "target"], tmp.loc[tmp.ion == "y", "prediction"])

        feats["cos"] = feats["dotprod"] / (np.linalg.norm(tmp["target"], 2) * np.linalg.norm(tmp["prediction"], 2))
        feats["cos_ionb"] = feats["dotprod_ionb"] / (np.linalg.norm(tmp.loc[tmp.ion == "b", "target"], 2) * np.linalg.norm(tmp.loc[tmp.ion == "b", "prediction"], 2))
        feats["cos_iony"] = feats["dotprod_iony"] / (np.linalg.norm(tmp.loc[tmp.ion == "y", "target"], 2) * np.linalg.norm(tmp.loc[tmp.ion == "y", "prediction"], 2))

        rescore_features = rescore_features.append(feats, ignore_index=True)
    return rescore_features

def calculate_features(path_to_pred_and_emp, path_to_out, num_cpu):
    """
    parallelize calculation of features and write them into a csv file
    """

    myPool = multiprocessing.Pool(num_cpu)

    df = pd.read_csv(path_to_pred_and_emp)

    peptides = list(df.spec_id.unique())
    split_peptides = [peptides[i * len(peptides) // num_cpu: (i + 1) * len(peptides) // num_cpu] for i in range(num_cpu)]

    results = []

    for i in range(num_cpu):
        tmp = df[df["spec_id"].isin(split_peptides[i])]
        results.append(myPool.apply_async(compute_features, args=(tmp, )))

    myPool.close()
    myPool.join()

    all_results = []
    for r in results:
        all_results.append(r.get())
    all_results = pd.concat(all_results)

    all_results.to_csv(path_to_out, index=False)

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
    all_features = all_features.drop(all_features.loc[:,all_features.columns.str.endswith('.1')], axis=1)

    norm_cols = ['spec_pearson_norm', 'ionb_pearson_norm', 'iony_pearson_norm',
        'spec_spearman_norm', 'ionb_spearman_norm', 'iony_spearman_norm',
        'spec_mse_norm', 'ionb_mse_norm', 'iony_mse_norm', 'min_abs_diff_norm',
        'max_abs_diff_norm', 'abs_diff_Q1_norm', 'abs_diff_Q2_norm',
        'abs_diff_Q3_norm', 'mean_abs_diff_norm', 'std_abs_diff_norm',
        'ionb_min_abs_diff_norm', 'ionb_max_abs_diff_norm',
        'ionb_abs_diff_Q1_norm', 'ionb_abs_diff_Q2_norm',
        'ionb_abs_diff_Q3_norm', 'ionb_mean_abs_diff_norm',
        'ionb_std_abs_diff_norm', 'iony_min_abs_diff_norm',
        'iony_max_abs_diff_norm', 'iony_abs_diff_Q1_norm',
        'iony_abs_diff_Q2_norm', 'iony_abs_diff_Q3_norm',
        'iony_mean_abs_diff_norm', 'iony_std_abs_diff_norm', 'dotprod_norm',
        'dotprod_ionb_norm', 'dotprod_iony_norm', 'cos_norm', 'cos_ionb_norm',
        'cos_iony_norm', 'spec_pearson', 'ionb_pearson', 'iony_pearson',
        'spec_spearman', 'ionb_spearman', 'iony_spearman', 'spec_mse',
        'ionb_mse', 'iony_mse', 'min_abs_diff', 'max_abs_diff', 'abs_diff_Q1',
        'abs_diff_Q2', 'abs_diff_Q3', 'mean_abs_diff', 'std_abs_diff',
        'ionb_min_abs_diff', 'ionb_max_abs_diff', 'ionb_abs_diff_Q1',
        'ionb_abs_diff_Q2', 'ionb_abs_diff_Q3', 'ionb_mean_abs_diff',
        'ionb_std_abs_diff', 'iony_min_abs_diff', 'iony_max_abs_diff',
        'iony_abs_diff_Q1', 'iony_abs_diff_Q2', 'iony_abs_diff_Q3',
        'iony_mean_abs_diff', 'iony_std_abs_diff', 'dotprod', 'dotprod_ionb',
        'dotprod_iony', 'cos', 'cos_ionb', 'cos_iony']

    all_features.to_csv(path_to_pin.rstrip('.pin')+'_all_features.csv', sep=',', index=False)
    norm_features = pd.DataFrame(StandardScaler().fit_transform(X=all_features.loc[:, norm_cols]), columns=norm_cols)
    all_features.loc[:, norm_features.columns] = norm_features
    all_features.to_csv(path_to_pin.rstrip('.pin')+'_all_features_norm.csv', sep=',', index=False)


def write_pin_files(path_to_features, savepath):
    """
    Given a dataframe with all the features, writes three different pin files:
    _only_rescore.pin with only the rescore features
    _all_percolator.pin with all the percolator features
    _all_features.pin with all the rescore and percolator features
    :param path_to_features: pd.DataFrame obtained from rescore.join_features()
    :param savepath: path to save the pin files
    """
    # columns to save
    percolator_features = ['ExpMass', 'CalcMass', 'RawScore', 'DeNovoScore',
        'ScoreRatio', 'Energy', 'lnEValue', 'IsotopeError',
        'lnExplainedIonCurrentRatio', 'lnNTermIonCurrentRatio',
        'lnCTermIonCurrentRatio', 'lnMS2IonCurrent', 'Mass', 'PepLen', 'dM',
        'absdM', 'MeanErrorTop7', 'sqMeanErrorTop7', 'StdevErrorTop7',
        'Charge2', 'Charge3', 'Charge4', 'Charge5', 'Charge6', 'enzN', 'enzC',
        'enzInt', 'ptm', 'A-Freq', 'C-Freq', 'D-Freq', 'E-Freq', 'F-Freq',
        'G-Freq', 'H-Freq', 'I-Freq', 'K-Freq', 'L-Freq', 'M-Freq', 'N-Freq',
        'P-Freq', 'Q-Freq', 'R-Freq', 'S-Freq', 'T-Freq', 'V-Freq', 'W-Freq',
        'Y-Freq', 'B-Freq', 'Z-Freq', 'J-Freq', 'X-Freq', 'U-Freq', 'O-Freq']

    rescore_features = ['spec_pearson_norm', 'ionb_pearson_norm',
        'iony_pearson_norm', 'spec_spearman_norm', 'ionb_spearman_norm',
        'iony_spearman_norm', 'spec_mse_norm', 'ionb_mse_norm', 'iony_mse_norm',
         'min_abs_diff_iontype_norm', 'max_abs_diff_iontype_norm',
        'min_abs_diff_norm', 'max_abs_diff_norm', 'abs_diff_Q1_norm',
        'abs_diff_Q2_norm', 'abs_diff_Q3_norm', 'mean_abs_diff_norm',
        'std_abs_diff_norm', 'ionb_min_abs_diff_norm', 'ionb_max_abs_diff_norm',
        'ionb_abs_diff_Q1_norm', 'ionb_abs_diff_Q2_norm',
        'ionb_abs_diff_Q3_norm', 'ionb_mean_abs_diff_norm',
        'ionb_std_abs_diff_norm', 'iony_min_abs_diff_norm',
        'iony_max_abs_diff_norm', 'iony_abs_diff_Q1_norm',
        'iony_abs_diff_Q2_norm', 'iony_abs_diff_Q3_norm',
        'iony_mean_abs_diff_norm', 'iony_std_abs_diff_norm', 'dotprod_norm',
        'dotprod_ionb_norm', 'dotprod_iony_norm', 'cos_norm', 'cos_ionb_norm',
        'cos_iony_norm', 'spec_pearson', 'ionb_pearson', 'iony_pearson',
        'spec_spearman', 'ionb_spearman', 'iony_spearman', 'spec_mse',
        'ionb_mse', 'iony_mse', 'min_abs_diff_iontype', 'max_abs_diff_iontype',
        'min_abs_diff', 'max_abs_diff', 'abs_diff_Q1', 'abs_diff_Q2',
        'abs_diff_Q3', 'mean_abs_diff', 'std_abs_diff', 'ionb_min_abs_diff',
        'ionb_max_abs_diff', 'ionb_abs_diff_Q1', 'ionb_abs_diff_Q2',
        'ionb_abs_diff_Q3', 'ionb_mean_abs_diff', 'ionb_std_abs_diff',
        'iony_min_abs_diff', 'iony_max_abs_diff', 'iony_abs_diff_Q1',
        'iony_abs_diff_Q2', 'iony_abs_diff_Q3', 'iony_mean_abs_diff',
        'iony_std_abs_diff', 'dotprod', 'dotprod_ionb', 'dotprod_iony', 'cos',
        'cos_ionb', 'cos_iony']

    all_features = pd.read_csv(path_to_features, sep=',')
    # Writing files with appropriate columns
    all_features.loc[:, ['SpecId', 'Label', 'ScanNr'] + rescore_features + ['Peptide', 'Proteins']].fillna(value=0).to_csv('{}_rescore.pin'.format(savepath), sep='\t', index=False)
    all_features.loc[:, ['SpecId', 'Label', 'ScanNr'] + percolator_features + ['Peptide', 'Proteins']].fillna(value=0).to_csv('{}_percolator.pin'.format(savepath), sep='\t', index=False)
    all_features.loc[:, ['SpecId', 'Label', 'ScanNr'] + percolator_features + rescore_features + ['Peptide', 'Proteins']].fillna(value=0).to_csv('{}_all_features.pin'.format(savepath), sep='\t', index=False)
    return None

def format_output(path_to_pout, search_engine, savepath, fname, fig=True):
    out = pd.concat([pd.read_csv(path_to_pout, sep='\t'), pd.read_csv(path_to_pout+"_dec", sep='\t')])
    inp = pd.read_csv(fname + '_all_features.csv')

    if search_engine == "MSGFPlus":
        out['Label'] = [-1 if p.startswith('XXX') else 1 for p in out.proteinIds]
        score = 'lnEValue'

    if fig:
        f, axes = plt.subplots(4, 3, figsize=(16, 21))

        sns.boxplot(data=inp, y=score, x='Label', ax=axes[0][0])

        sns.distplot(inp.loc[inp.Label == -1, score], kde=False, ax=axes[0][1])
        sns.distplot(inp.loc[inp.Label == 1, score], kde=False, ax=axes[0][1])

        sns.distplot(inp.loc[inp.Label == -1, score], hist=False, ax=axes[0][2])
        sns.distplot(inp.loc[inp.Label == 1, score], hist=False, ax=axes[0][2])

        sns.boxplot(data=out, y='score', x='Label', ax=axes[1][0])

        sns.distplot(out.loc[out.Label == -1, 'score'], kde=False, ax=axes[1][1])
        sns.distplot(out.loc[out.Label == 1, 'score'], kde=False, ax=axes[1][1])

        sns.distplot(out.loc[out.Label == -1, 'score'], hist=False, ax=axes[1][2])
        sns.distplot(out.loc[out.Label == 1, 'score'], hist=False, ax=axes[1][2])

        sns.boxplot(data=out, y='q-value', x='Label', ax=axes[2][0])

        sns.distplot(out.loc[out.Label == -1, 'q-value'], kde=False, ax=axes[2][1])
        sns.distplot(out.loc[out.Label == 1, 'q-value'], kde=False, ax=axes[2][1])

        sns.distplot(out.loc[out.Label == -1, 'q-value'], hist=False, ax=axes[2][2])
        sns.distplot(out.loc[out.Label == 1, 'q-value'], hist=False, ax=axes[2][2])

        sns.boxplot(data=out, y='posterior_error_prob', x='Label', ax=axes[3][0])

        sns.distplot(out.loc[out.Label == -1, 'posterior_error_prob'], kde=False, ax=axes[3][1])
        sns.distplot(out.loc[out.Label == 1, 'posterior_error_prob'], kde=False, ax=axes[3][1])

        sns.distplot(out.loc[out.Label == -1, 'posterior_error_prob'], hist=False, ax=axes[3][2])
        sns.distplot(out.loc[out.Label == 1, 'posterior_error_prob'], hist=False, ax=axes[3][2])

        axes[0][1].set_xlim([np.min(inp[score]), np.max(inp[score])])
        axes[0][2].set_xlim([np.min(inp[score]), np.max(inp[score])])

        axes[1][1].set_xlim([np.min(out['score']), np.max(out['score'])])
        axes[1][2].set_xlim([np.min(out['score']), np.max(out['score'])])

        axes[2][1].set_xlim([np.min(out['q-value']), np.max(out['q-value'])])
        axes[2][2].set_xlim([np.min(out['q-value']), np.max(out['q-value'])])

        axes[3][1].set_xlim([np.min(out['posterior_error_prob']), np.max(out['posterior_error_prob'])])
        axes[3][2].set_xlim([np.min(out['posterior_error_prob']), np.max(out['posterior_error_prob'])])

        f.savefig(savepath)
