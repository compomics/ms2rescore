"""
Functions necessary to run the rescore algorithm. Currently supports MSGF+ with
concatenated searches.
"""

import multiprocessing
import os
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler

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
        print([mod["name"], str(mod["mass_shift"]), "opt", aa])
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


def norm_features(path_to_features):
    """
    Normalize the features obtained from MS2PIP

    :param path_to_features: string, path to MS2PIP features

    Returns
    :pd.DataFrame features, includes normalized MS2PIP features
    """

    # Read features csv file and fillna
    features = pd.read_csv(path_to_features)
    features = features.fillna(0)

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

    norm_features = pd.DataFrame(StandardScaler().fit_transform(X=features.loc[:, norm_cols]), columns=norm_cols)
    features.loc[:, norm_features.columns] = norm_features
    features.to_csv(path_to_features, sep=',', index=False)

def write_pin_files(path_to_features, path_to_pep, savepath):
    """
    Given a dataframe with all the features, writes a pin file
    """
    # feature columns
    features = ['spec_pearson_norm', 'ionb_pearson_norm',
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

    pep = pd.read_csv(path_to_pep, sep=' ')
    complete_df = pd.merge(all_features, pep, on=['spec_id', 'charge'])
    complete_df.rename(mapper={'spec_id': 'SpecId', 'peptide': 'Peptide'}, axis='columns', inplace=True)

    # Add artificial protein column, scan numbers and flanking aminoacids to
    # peptide sequences
    complete_df['ScanNr'] = list(range(len(complete_df)))
    complete_df['Proteins'] = complete_df['Peptide']
    complete_df['Peptide'] = 'X.' + complete_df.Peptide + '.X'

    # Writing files with ordered columns
    complete_df.loc[:, ['SpecId', 'Label', 'ScanNr'] + features + ['Peptide', 'Proteins']].fillna(value=0).to_csv('{}.pin'.format(savepath), sep='\t', index=False)
    return None
