"""
Functions necessary to run the rescore algorithm. Currently supports MSGF+ with
concatenated searches.
"""

# Standard library
import logging
import warnings
import sys
import multiprocessing
import os

# Third party
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error as mse
from tqdm import tqdm


def make_ms2pip_config(options, filename='ms2pip_config.txt'):
    """
    write configuration file for ms2pip based on what's on the rescore config
    file.
    """
    cwd = os.getcwd()
    ms2pip_config = open(os.path.join(cwd, filename), 'wt')

    if "frag" in options["ms2pip"]:
        ms2pip_config.write("frag_method={}\n".format(options["ms2pip"]["frag"]))
    else:
        # Assume HCD
        ms2pip_config.write("frag_method=HCD\n")

    if "frag_error" in options["ms2pip"]:
        ms2pip_config.write("frag_error={}\n".format(options["ms2pip"]["frag_error"]))
    else:
        if options["ms2pip"]["frag"] == "CID":
            ms2pip_config.write("frag_error=0.8\n")
        elif options["ms2pip"]["frag"] == "phospho":
            ms2pip_config.write("frag_error=0.02\n")
        else:
            # Assume HCD
            ms2pip_config.write("frag_error=0.02\n")

    ms2pip_config.write("\n")

    modifications = options["ms2pip"]["modifications"]
    for mod in modifications:
        if mod["amino_acid"] is None and mod["n_term"]:
            aa = "N-term"
        else:
            aa = mod["amino_acid"]
        tmp = ','.join([mod["name"], str(mod["mass_shift"]), "opt", aa])
        ms2pip_config.write("ptm=" + tmp + "\n")

    ms2pip_config.close()


def df_to_dict(df):
    """
    Create easy to access dict from pred_and_emp
    """
    preds_dict = {}
    preds_list = df[['spec_id', 'charge', 'ion', 'target', 'prediction']].values.tolist()

    for row in preds_list:
        spec_id = row[0]
        if spec_id in preds_dict.keys():
            if row[2] in preds_dict[spec_id]['target']:
                preds_dict[spec_id]['target'][row[2]].append(row[3])
                preds_dict[spec_id]['prediction'][row[2]].append(row[4])
            else:
                preds_dict[spec_id]['target'][row[2]] = [row[3]]
                preds_dict[spec_id]['prediction'][row[2]] = [row[4]]
        else:
            preds_dict[spec_id] = {
                'charge': row[1],
                'target': {row[2]: [row[3]]},
                'prediction': {row[2]: [row[4]]}
            }
    return preds_dict


def compute_features(df):
    """
    Compute ReScore features
    """

    preds_dict = df_to_dict(df)

    rescore_features = []
    spec_ids = []
    charges = []

    feature_names = [
        'spec_pearson_norm', 'ionb_pearson_norm', 'iony_pearson_norm',
        'spec_mse_norm', 'ionb_mse_norm', 'iony_mse_norm',
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
        'cos_ionb', 'cos_iony'
    ]

    # Suppress RuntimeWarnings about invalid values
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        for spec_id, preds in preds_dict.items():
            spec_ids.append(spec_id)
            charges.append(preds['charge'])

            # Create numpy arrays
            target_b = np.array(preds['target']['B'])
            target_y = np.array(preds['target']['Y'])
            target_all = np.concatenate([target_b, target_y])
            prediction_b = np.array(preds['prediction']['B'])
            prediction_y = np.array(preds['prediction']['Y'])
            prediction_all = np.concatenate([prediction_b, prediction_y])

            target_b_unlog = 2 ** target_b - 0.001
            target_y_unlog = 2 ** target_y - 0.001
            target_all_unlog = 2 ** target_all - 0.001
            prediction_b_unlog = 2 ** prediction_b - 0.001
            prediction_y_unlog = 2 ** prediction_y - 0.001
            prediction_all_unlog = 2 ** prediction_all - 0.001

            # Calculate absolute differences
            abs_diff_b = np.abs(target_b - prediction_b)
            abs_diff_y = np.abs(target_y - prediction_y)
            abs_diff_all = np.abs(target_all - prediction_all)

            abs_diff_b_unlog = np.abs(target_b_unlog - prediction_b_unlog)
            abs_diff_y_unlog = np.abs(target_y_unlog - prediction_y_unlog)
            abs_diff_all_unlog = np.abs(target_all_unlog - prediction_all_unlog)

            # Add features
            feats = np.array([
                #spec_id,
                #preds['charge'],

                # Features between spectra in log space
                pearsonr(target_all, prediction_all)[0],  # Pearson all ions
                pearsonr(target_b, prediction_b)[0],  # Pearson b ions
                pearsonr(target_y, prediction_y)[0],  # Pearson y ions

                mse(target_all, prediction_all),  # MSE all ions
                mse(target_b, prediction_b),  # MSE b ions
                mse(target_y, prediction_y),  # MSE y ions

                np.min(abs_diff_all),  # min_abs_diff_norm
                np.max(abs_diff_all),  # max_abs_diff_norm
                np.quantile(abs_diff_all, 0.25),  # abs_diff_Q1_norm
                np.quantile(abs_diff_all, 0.5),  # abs_diff_Q2_norm
                np.quantile(abs_diff_all, 0.75),  # abs_diff_Q3_norm
                np.mean(abs_diff_all),  # mean_abs_diff_norm
                np.std(abs_diff_all),  # std_abs_diff_norm

                np.min(abs_diff_b),  # ionb_min_abs_diff_norm
                np.max(abs_diff_b),  # ionb_max_abs_diff_norm
                np.quantile(abs_diff_b, 0.25),  # ionb_abs_diff_Q1_norm
                np.quantile(abs_diff_b, 0.5),  # ionb_abs_diff_Q2_norm
                np.quantile(abs_diff_b, 0.75),  # ionb_abs_diff_Q3_norm
                np.mean(abs_diff_b),  # ionb_mean_abs_diff_norm
                np.std(abs_diff_b),  # ionb_std_abs_diff_norm

                np.min(abs_diff_y),  # iony_min_abs_diff_norm
                np.max(abs_diff_y),  # iony_max_abs_diff_norm
                np.quantile(abs_diff_y, 0.25),  # iony_abs_diff_Q1_norm
                np.quantile(abs_diff_y, 0.5),  # iony_abs_diff_Q2_norm
                np.quantile(abs_diff_y, 0.75),  # iony_abs_diff_Q3_norm
                np.mean(abs_diff_y),  # iony_mean_abs_diff_norm
                np.std(abs_diff_y),  # iony_std_abs_diff_norm

                np.dot(target_all, prediction_all),  # Dot product all ions
                np.dot(target_b, prediction_b),  # Dot product b ions
                np.dot(target_y, prediction_y),  # Dot product y ions

                np.dot(target_all, prediction_all) / (np.linalg.norm(target_all, 2) * np.linalg.norm(prediction_all, 2)),  # Cos similarity all ions
                np.dot(target_b, prediction_b) / (np.linalg.norm(target_b, 2) * np.linalg.norm(prediction_b, 2)),  # Cos similarity b ions
                np.dot(target_y, prediction_y) / (np.linalg.norm(target_y, 2) * np.linalg.norm(prediction_y, 2)),  # Cos similarity y ions

                # Same features in normal space
                pearsonr(target_all_unlog, prediction_all_unlog)[0],  # Pearson all ions
                pearsonr(target_b_unlog, prediction_b_unlog)[0],  # Pearson b ions
                pearsonr(target_y_unlog, prediction_y_unlog)[0],  # Pearson y ions

                spearmanr(target_all_unlog, prediction_all_unlog)[0],  # Spearman all ions
                spearmanr(target_b_unlog, prediction_b_unlog)[0],  # Spearman b ions
                spearmanr(target_y_unlog, prediction_y_unlog)[0],  # Spearman y ions

                mse(target_all_unlog, prediction_all_unlog),  # MSE all ions
                mse(target_b_unlog, prediction_b_unlog),  # MSE b ions
                mse(target_y_unlog, prediction_y_unlog),  # MSE y ions,

                0 if np.min(abs_diff_b_unlog) <= np.min(abs_diff_y_unlog) else 1,  # Ion type with min absolute difference
                0 if np.max(abs_diff_b_unlog) >= np.max(abs_diff_y_unlog) else 1,  # Ion type with max absolute difference

                np.min(abs_diff_all_unlog),  # min_abs_diff
                np.max(abs_diff_all_unlog),  # max_abs_diff
                np.quantile(abs_diff_all_unlog, 0.25),  # abs_diff_Q1
                np.quantile(abs_diff_all_unlog, 0.5),  # abs_diff_Q2
                np.quantile(abs_diff_all_unlog, 0.75),  # abs_diff_Q3
                np.mean(abs_diff_all_unlog),  # mean_abs_diff
                np.std(abs_diff_all_unlog),  # std_abs_diff

                np.min(abs_diff_b_unlog),  # ionb_min_abs_diff
                np.max(abs_diff_b_unlog),  # ionb_max_abs_diff_norm
                np.quantile(abs_diff_b_unlog, 0.25),  # ionb_abs_diff_Q1
                np.quantile(abs_diff_b_unlog, 0.5),  # ionb_abs_diff_Q2
                np.quantile(abs_diff_b_unlog, 0.75),  # ionb_abs_diff_Q3
                np.mean(abs_diff_b_unlog),  # ionb_mean_abs_diff
                np.std(abs_diff_b_unlog),  # ionb_std_abs_diff

                np.min(abs_diff_y_unlog),  # iony_min_abs_diff
                np.max(abs_diff_y_unlog),  # iony_max_abs_diff
                np.quantile(abs_diff_y_unlog, 0.25),  # iony_abs_diff_Q1
                np.quantile(abs_diff_y_unlog, 0.5),  # iony_abs_diff_Q2
                np.quantile(abs_diff_y_unlog, 0.75),  # iony_abs_diff_Q3
                np.mean(abs_diff_y_unlog),  # iony_mean_abs_diff
                np.std(abs_diff_y_unlog),  # iony_std_abs_diff

                np.dot(target_all_unlog, prediction_all_unlog),  # Dot product all ions
                np.dot(target_b_unlog, prediction_b_unlog),  # Dot product b ions
                np.dot(target_y_unlog, prediction_y_unlog),  # Dot product y ions

                np.dot(target_all_unlog, prediction_all_unlog) / (np.linalg.norm(target_all_unlog, 2) * np.linalg.norm(prediction_all_unlog, 2)),  # Cos similarity all ions
                np.dot(target_b_unlog, prediction_b_unlog) / (np.linalg.norm(target_b_unlog, 2) * np.linalg.norm(prediction_b_unlog, 2)),  # Cos similarity b ions
                np.dot(target_y_unlog, prediction_y_unlog) / (np.linalg.norm(target_y_unlog, 2) * np.linalg.norm(prediction_y_unlog, 2)),  # Cos similarity y ions
            ], dtype=np.float64)

            rescore_features.append(feats)

    rescore_features = np.vstack(rescore_features)
    rescore_features = pd.DataFrame(rescore_features, columns=feature_names)
    rescore_features['spec_id'] = spec_ids
    rescore_features['charge'] = charges

    return rescore_features


def calculate_features(path_to_pred_and_emp, path_to_out, num_cpu):
    """
    parallelize calculation of features and write them into a csv file
    """

    logging.debug("Reading MS2PIP predictions file")
    df = pd.read_csv(path_to_pred_and_emp)
    df[['prediction', 'target']] = df[['prediction', 'target']].clip(lower=np.log2(0.001))

    logging.debug("Computing features")
    if len(df['spec_id'].unique()) < 10000:
        all_results = compute_features(df)
    else:
        # Split up df into list of chunk_size df's (will be divided over num_cpu)
        chunk_size = 100
        spec_ids = list(df['spec_id'].unique())
        split_df = [df[df['spec_id'].isin(spec_ids[i * len(spec_ids) // chunk_size: (i + 1) * len(spec_ids) // chunk_size])] for i in range(chunk_size)]

        # Use imap, so we can use a tqdm progress bar
        with multiprocessing.Pool(int(num_cpu)) as p:
            all_results = list(tqdm(p.imap(compute_features, split_df), total=chunk_size))

        all_results = pd.concat(all_results)

    logging.debug("Writing to file")
    all_results.to_csv(path_to_out, index=False)


def redo_pin_tabs(pin_filename):
    """
    Replaces triple pipe (`|||`) with tab in given file.
    """
    with open(pin_filename, 'rt') as pin_in:
        with open(pin_filename + '_tmp', 'wt') as pin_out:
            for line in pin_in:
                pin_out.write(line.replace('|||', '\t'))
    os.remove(pin_filename)
    os.rename(pin_filename + '_tmp', pin_filename)


def write_pin_files(path_to_features, path_to_pep, savepath, feature_sets=None):
    """
    Given a dataframe with all the features, writes three PIN files.
    Writes three PIN files: search_engine_features, ms2pip_features and
    all_features
    Arguments:
    path_to_features: str of path to features or pd.DataFrame with features
    path_to_pep: str of path to peprec or pd.DataFrame with peprec
    save_path: path to save PIN files into
    feature_sets: list with feature sets for which to write PIN files. Can
    contain `all`, `search_engine` and `ms2pip`.
    """

    if not feature_sets:
        feature_sets = ['all', 'ms2pip', 'searchengine']

    all_features = pd.read_csv(path_to_features, sep=',')

    if type(path_to_pep) == str:
        with open(path_to_pep, 'rt') as f:
            line = f.readline()
            if line[:7] != 'spec_id':
                logging.critical('PEPREC file should start with header column')
                exit(1)
            sep = line[7]

        # Convert Protein literal list to pseudo-PIN notation:
        # not yet tab-separated, but pipe-separated
        # to avoid isues in pd.DataFrame.to_csv()
        protein_converter = lambda prot_list: '|||'.join(prot_list.strip("[]'").split("', '"))
        pep = pd.read_csv(path_to_pep,
                           sep=sep,
                           index_col=False,
                           dtype={"spec_id": str, "modifications": str},
                           converters={
                               'Proteins': protein_converter,
                               'protein_list': protein_converter
                            })
    else:
        pep = path_to_pep

    pep['modifications'].fillna("-", inplace=True)
    pep = pep.replace(-np.inf, 0.0).replace(np.inf, 0.0)

    # Prepare Proteins column for PIN
    # If no Proteins column in PEPREC, fill with peptide
    if not 'Proteins' in pep.columns:
        if 'protein_list' in pep.columns:
            pep['Proteins'] = pep['protein_list']
        elif 'ModPeptide' in pep.columns:
            pep['Proteins'] = pep['ModPeptide']
        else:
            pep['Proteins'] = pep['peptide']

    peprec_cols = ['spec_id', 'peptide', 'modifications', 'charge']
    pin_columns = ['SpecId', 'ScanNr', 'Label', 'ModPeptide', 'Proteins', 'TITLE']

    # Get list with feature names split by type of feature
    search_engine_feature_names = [col for col in pep.columns if (col not in peprec_cols) and (col not in pin_columns)]
    ms2pip_feature_names = [col for col in all_features.columns if col not in peprec_cols]

    # Merge ms2pip_features and peprec DataFrames
    complete_df = pd.merge(all_features, pep, on=['spec_id', 'charge'])
    complete_df = complete_df.fillna(value=0).reset_index(drop=True)

    # Add missing columns if necessary
    if not 'ScanNr' in complete_df.columns:
        complete_df['ScanNr'] = complete_df.index
    if not 'SpecId' in complete_df.columns:
        complete_df['SpecId'] = complete_df['spec_id']

    # Write PIN files with ordered columns
    if 'all' in feature_sets:
        complete_df[['SpecId', 'Label', 'ScanNr'] + ms2pip_feature_names + search_engine_feature_names + ['ModPeptide', 'Proteins']]\
            .to_csv('{}_allfeatures.pin'.format(savepath), sep='\t', header=True, index=False)
    if 'ms2pip' in feature_sets:
        complete_df[['SpecId', 'Label', 'ScanNr'] + ms2pip_feature_names + ['ModPeptide', 'Proteins']]\
            .to_csv('{}_ms2pipfeatures.pin'.format(savepath), sep='\t', header=True, index=False)
    if 'searchengine' in feature_sets:
        complete_df[['SpecId', 'Label', 'ScanNr'] + search_engine_feature_names + ['ModPeptide', 'Proteins']]\
            .to_csv('{}_searchenginefeatures.pin'.format(savepath), sep='\t', header=True, index=False)

    for features_version in feature_sets:
        redo_pin_tabs('{}_{}features.pin'.format(savepath, features_version))
