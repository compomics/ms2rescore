"""
Functions necessary to run the rescore algorithm. Currently supports MSGF+ with
concatenated searches.
"""

import logging
import multiprocessing
import os
import warnings
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error as mse
from tqdm import tqdm

from ms2rescore._exceptions import MS2ReScoreError

logger = logging.getLogger(__name__)


def make_ms2pip_config(
    options: Dict,
    filename: Optional[Union[str, os.PathLike]] = "ms2pip_config.txt"
):
    """Write ms2pip configuration file based on ms2pip config dict."""
    with open(filename, "wt") as ms2pip_config:
        if "frag" in options:
            ms2pip_config.write("frag_method={}\n".format(options["frag"]))
        if "model" in options:
            ms2pip_config.write("model={}\n".format(options["model"]))
        else:
            # Assume HCD
            ms2pip_config.write("frag_method=HCD\n")

        if "frag_error" in options:
            ms2pip_config.write("frag_error={}\n".format(options["frag_error"]))
        else:
            if options["frag"] == "CID":
                ms2pip_config.write("frag_error=0.8\n")
            elif options["frag"] == "phospho":
                ms2pip_config.write("frag_error=0.02\n")
            else:
                # Assume HCD
                ms2pip_config.write("frag_error=0.02\n")

        ms2pip_config.write("\n")

        modifications = options["modifications"]
        for mod in modifications:
            if mod["amino_acid"] is None and mod["n_term"]:
                aa = "N-term"
            elif mod["amino_acid"] is None and mod["c_term"]:
                aa = "C-term"
            else:
                aa = mod["amino_acid"]
            tmp = ",".join([mod["name"], str(mod["mass_shift"]), "opt", aa])
            ms2pip_config.write("ptm=" + tmp + "\n")


def df_to_dict(df):
    """Create easy to access dict from pred_and_emp."""
    preds_dict = {}
    preds_list = df[
        ["spec_id", "charge", "ion", "target", "prediction"]
    ].values.tolist()

    for row in preds_list:
        spec_id = row[0]
        if spec_id in preds_dict.keys():
            if row[2] in preds_dict[spec_id]["target"]:
                preds_dict[spec_id]["target"][row[2]].append(row[3])
                preds_dict[spec_id]["prediction"][row[2]].append(row[4])
            else:
                preds_dict[spec_id]["target"][row[2]] = [row[3]]
                preds_dict[spec_id]["prediction"][row[2]] = [row[4]]
        else:
            preds_dict[spec_id] = {
                "charge": row[1],
                "target": {row[2]: [row[3]]},
                "prediction": {row[2]: [row[4]]},
            }
    return preds_dict


def compute_features(df):
    """Compute ReScore features."""

    preds_dict = df_to_dict(df)

    rescore_features = []
    spec_ids = []
    charges = []

    feature_names = [
        "spec_pearson_norm",
        "ionb_pearson_norm",
        "iony_pearson_norm",
        "spec_mse_norm",
        "ionb_mse_norm",
        "iony_mse_norm",
        "min_abs_diff_norm",
        "max_abs_diff_norm",
        "abs_diff_Q1_norm",
        "abs_diff_Q2_norm",
        "abs_diff_Q3_norm",
        "mean_abs_diff_norm",
        "std_abs_diff_norm",
        "ionb_min_abs_diff_norm",
        "ionb_max_abs_diff_norm",
        "ionb_abs_diff_Q1_norm",
        "ionb_abs_diff_Q2_norm",
        "ionb_abs_diff_Q3_norm",
        "ionb_mean_abs_diff_norm",
        "ionb_std_abs_diff_norm",
        "iony_min_abs_diff_norm",
        "iony_max_abs_diff_norm",
        "iony_abs_diff_Q1_norm",
        "iony_abs_diff_Q2_norm",
        "iony_abs_diff_Q3_norm",
        "iony_mean_abs_diff_norm",
        "iony_std_abs_diff_norm",
        "dotprod_norm",
        "dotprod_ionb_norm",
        "dotprod_iony_norm",
        "cos_norm",
        "cos_ionb_norm",
        "cos_iony_norm",
        "spec_pearson",
        "ionb_pearson",
        "iony_pearson",
        "spec_spearman",
        "ionb_spearman",
        "iony_spearman",
        "spec_mse",
        "ionb_mse",
        "iony_mse",
        "min_abs_diff_iontype",
        "max_abs_diff_iontype",
        "min_abs_diff",
        "max_abs_diff",
        "abs_diff_Q1",
        "abs_diff_Q2",
        "abs_diff_Q3",
        "mean_abs_diff",
        "std_abs_diff",
        "ionb_min_abs_diff",
        "ionb_max_abs_diff",
        "ionb_abs_diff_Q1",
        "ionb_abs_diff_Q2",
        "ionb_abs_diff_Q3",
        "ionb_mean_abs_diff",
        "ionb_std_abs_diff",
        "iony_min_abs_diff",
        "iony_max_abs_diff",
        "iony_abs_diff_Q1",
        "iony_abs_diff_Q2",
        "iony_abs_diff_Q3",
        "iony_mean_abs_diff",
        "iony_std_abs_diff",
        "dotprod",
        "dotprod_ionb",
        "dotprod_iony",
        "cos",
        "cos_ionb",
        "cos_iony",
    ]

    # Suppress RuntimeWarnings about invalid values
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        for spec_id, preds in preds_dict.items():
            spec_ids.append(spec_id)
            charges.append(preds["charge"])

            # Create numpy arrays
            target_b = np.array(preds["target"]["B"])
            target_y = np.array(preds["target"]["Y"])
            target_all = np.concatenate([target_b, target_y])
            prediction_b = np.array(preds["prediction"]["B"])
            prediction_y = np.array(preds["prediction"]["Y"])
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
            feats = np.array(
                [
                    # spec_id,
                    # preds['charge'],
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
                    np.dot(target_all, prediction_all)
                    / (
                        np.linalg.norm(target_all, 2)
                        * np.linalg.norm(prediction_all, 2)
                    ),  # Cos similarity all ions
                    np.dot(target_b, prediction_b)
                    / (
                        np.linalg.norm(target_b, 2) * np.linalg.norm(prediction_b, 2)
                    ),  # Cos similarity b ions
                    np.dot(target_y, prediction_y)
                    / (
                        np.linalg.norm(target_y, 2) * np.linalg.norm(prediction_y, 2)
                    ),  # Cos similarity y ions
                    # Same features in normal space
                    pearsonr(target_all_unlog, prediction_all_unlog)[
                        0
                    ],  # Pearson all ions
                    pearsonr(target_b_unlog, prediction_b_unlog)[0],  # Pearson b ions
                    pearsonr(target_y_unlog, prediction_y_unlog)[0],  # Pearson y ions
                    spearmanr(target_all_unlog, prediction_all_unlog)[
                        0
                    ],  # Spearman all ions
                    spearmanr(target_b_unlog, prediction_b_unlog)[0],  # Spearman b ions
                    spearmanr(target_y_unlog, prediction_y_unlog)[0],  # Spearman y ions
                    mse(target_all_unlog, prediction_all_unlog),  # MSE all ions
                    mse(target_b_unlog, prediction_b_unlog),  # MSE b ions
                    mse(target_y_unlog, prediction_y_unlog),  # MSE y ions,
                    0
                    if np.min(abs_diff_b_unlog) <= np.min(abs_diff_y_unlog)
                    else 1,  # Ion type with min absolute difference
                    0
                    if np.max(abs_diff_b_unlog) >= np.max(abs_diff_y_unlog)
                    else 1,  # Ion type with max absolute difference
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
                    np.dot(
                        target_all_unlog, prediction_all_unlog
                    ),  # Dot product all ions
                    np.dot(target_b_unlog, prediction_b_unlog),  # Dot product b ions
                    np.dot(target_y_unlog, prediction_y_unlog),  # Dot product y ions
                    np.dot(target_all_unlog, prediction_all_unlog)
                    / (
                        np.linalg.norm(target_all_unlog, 2)
                        * np.linalg.norm(prediction_all_unlog, 2)
                    ),  # Cos similarity all ions
                    np.dot(target_b_unlog, prediction_b_unlog)
                    / (
                        np.linalg.norm(target_b_unlog, 2)
                        * np.linalg.norm(prediction_b_unlog, 2)
                    ),  # Cos similarity b ions
                    np.dot(target_y_unlog, prediction_y_unlog)
                    / (
                        np.linalg.norm(target_y_unlog, 2)
                        * np.linalg.norm(prediction_y_unlog, 2)
                    ),  # Cos similarity y ions
                ],
                dtype=np.float64,
            )

            rescore_features.append(feats)

    rescore_features = np.vstack(rescore_features)
    rescore_features = pd.DataFrame(rescore_features, columns=feature_names)
    rescore_features["spec_id"] = spec_ids
    rescore_features["charge"] = charges

    return rescore_features


def calculate_features(
    path_to_pred_and_emp, path_to_out, num_cpu, show_progress_bar=True
):
    """
    parallelize calculation of features and write them into a csv file
    """

    logger.debug("Reading MS2PIP predictions file")
    df = pd.read_csv(path_to_pred_and_emp)
    df[["prediction", "target"]] = df[["prediction", "target"]].clip(
        lower=np.log2(0.001)
    )

    logger.debug("Computing features")
    if len(df["spec_id"].unique()) < 10000:
        all_results = compute_features(df)
    else:
        # Split up df into list of chunk_size df's (will be divided over num_cpu)
        chunk_size = 100
        spec_ids = list(df["spec_id"].unique())
        split_df = [
            df[
                df["spec_id"].isin(
                    spec_ids[
                        i
                        * len(spec_ids)
                        // chunk_size: (i + 1)
                        * len(spec_ids)
                        // chunk_size
                    ]
                )
            ]
            for i in range(chunk_size)
        ]

        # Use imap, so we can use a tqdm progress bar
        with multiprocessing.Pool(int(num_cpu)) as p:
            if show_progress_bar:
                all_results = list(
                    tqdm(p.imap(compute_features, split_df), total=chunk_size)
                )
            else:
                all_results = list(p.imap(compute_features, split_df))
        all_results = pd.concat(all_results)

    logger.debug("Writing to file")
    all_results.to_csv(path_to_out, index=False)


def redo_pin_tabs(pin_filename):
    """
    Replaces triple pipe (`|||`) with tab in given file.
    """
    with open(pin_filename, "rt") as pin_in:
        with open(pin_filename + "_tmp", "wt") as pin_out:
            for line in pin_in:
                pin_out.write(line.replace("|||", "\t"))
    os.remove(pin_filename)
    os.rename(pin_filename + "_tmp", pin_filename)


def write_pin_files(
    peprec_path: str,
    savepath: str,
    searchengine_features_path: Optional[str] = None,
    ms2pip_features_path: Optional[str] = None,
    rt_features_path: Optional[str] = None,
    feature_sets: Optional[list] = None,
):
    """
    Write Percolator IN file.

    Write PIN files for each requested feature_set: search_engine_features,
    ms2pip_features, rt_features and all_features.

    Parameters
    ----------
    peprec_path: str
        Path to PEPREC file.
    save_path: str
        Directory and basename for PIN files.
    searchengine_features_path: {str, None}
        Path to CSV with search engine features.
    ms2pip_features_path: {str, None}
        Path to CSV with ms2pip features.
    rt_features_path: {str, None}
        Path to CSV with rt features.
    feature_sets: {list, None}
        Feature sets for which to write PIN files. Can contain `all`, `search_engine`,
        `ms2pip`, and `rt`.

    """
    if not feature_sets:
        feature_sets = ["all", "searchengine", "ms2pip", "rt"]

    # Read search engine features
    if ("searchengine" in feature_sets) or ("all" in feature_sets):
        searchengine_features = pd.read_csv(
            searchengine_features_path,
            sep=",",
            index_col=None
        )
    else:
        searchengine_features = None

    # Read MS²PIP features
    if ("ms2pip" in feature_sets) or ("all" in feature_sets):
        ms2pip_features = pd.read_csv(ms2pip_features_path, sep=",", index_col=None)
    else:
        ms2pip_features = None

    # Read RT features
    if ("rt" in feature_sets) or ("all" in feature_sets):
        rt_features = pd.read_csv(rt_features_path, sep=",", index_col=None)
    else:
        rt_features = None

    # Read PEPREC
    with open(peprec_path, "rt") as f:
        line = f.readline()
        if line[:7] != "spec_id":
            logger.critical("PEPREC file should start with header column")
            exit(1)
        sep = line[7]

    # Convert Protein literal list to pseudo-PIN notation:
    # not yet tab-separated, but pipe-separated
    # to avoid isues in pd.DataFrame.to_csv()
    protein_converter = lambda prot_list: "|||".join(
        prot_list.strip("[]'").split("', '")
    )
    pep = pd.read_csv(
        peprec_path,
        sep=sep,
        index_col=False,
        dtype={"spec_id": str, "modifications": str},
        converters={"Proteins": protein_converter, "protein_list": protein_converter},
    )

    pep["modifications"].fillna("-", inplace=True)
    pep = pep.replace(-np.inf, 0.0).replace(np.inf, 0.0)

    # Prepare Proteins column for PIN
    # If no Proteins column in PEPREC, fill with peptide
    if "Proteins" not in pep.columns:
        if "protein_list" in pep.columns:
            pep["Proteins"] = pep["protein_list"]
        elif "ModPeptide" in pep.columns:
            pep["Proteins"] = pep["ModPeptide"]
        else:
            pep["Proteins"] = pep["peptide"]

    # ModPeptide contains sequence with modifications (e.g. from MaxQuant)
    if "ModPeptide" not in pep.columns:
        pep["ModPeptide"] = pep["peptide"]
    pep.rename(
        columns={
            "peptide": "peptide_peprec",
            "ModPeptide": "Peptide",
            "label": "Label",
        },
        inplace=True,
    )

    # TODO: Fix duality in `observed_retention_time`: peprec column, also rt feature
    peprec_cols = [
        "spec_id",
        "peptide",
        "peptide_peprec",
        "modifications",
        "charge",
        "label",
        "psm_score",
        "observed_retention_time",
    ]
    pin_columns = [
        "SpecId",
        "ScanNr",
        "Label",
        "Peptide",
        "ModPeptide",
        "Proteins",
        "TITLE",
    ]

    # Get list with feature names split by type of feature
    search_engine_feature_names = [
        col
        for col in searchengine_features.columns
        if (col not in peprec_cols) and (col not in pin_columns)
    ]

    if ("ms2pip" in feature_sets) or ("all" in feature_sets):
        ms2pip_feature_names = [
            col
            for col in ms2pip_features.columns
            if (col not in peprec_cols) and (col not in pin_columns)
        ]
    else:
        ms2pip_feature_names = []

    if ("rt" in feature_sets) or ("all" in feature_sets):
        rt_feature_names = ["observed_retention_time"] + [
            col
            for col in rt_features.columns
            if (col not in peprec_cols) and (col not in pin_columns)
        ]
    else:
        rt_feature_names = []

    # Merge features and peprec DataFrames
    complete_df = pep
    for feature_set in [searchengine_features, ms2pip_features, rt_features]:
        if isinstance(feature_set, pd.DataFrame):
            on_cols = ["spec_id", "charge"]
            cols_to_use = (
                set(feature_set.columns)
                .difference(set(complete_df.columns))
                .union(set(on_cols))
            )
            # If spec_id and charge in feature_set columns, use merge
            if all(c in feature_set.columns for c in on_cols):
                complete_df = pd.merge(
                    complete_df, feature_set[cols_to_use], on=on_cols
                )
            # Else, just use concat
            else:
                if not len(complete_df) == len(feature_set):
                    raise MS2ReScoreError("Feature sets do not match.")
                complete_df = pd.concat(
                    [
                        complete_df.reset_index(drop=True),
                        feature_set.reset_index(drop=True)
                    ], axis=1
                )
    complete_df = complete_df.fillna(value=0).reset_index(drop=True)
    complete_df = complete_df.loc[:, ~complete_df.columns.duplicated()]

    # Add missing columns if necessary
    if "ScanNr" not in complete_df.columns:
        complete_df["ScanNr"] = complete_df.index
    if "SpecId" not in complete_df.columns:
        complete_df["SpecId"] = complete_df["spec_id"]
    if "Label" not in complete_df.columns:
        complete_df["Label"] = complete_df["label"]

    # Write PIN files with ordered columns
    # From Percolator documentation:
    # SpecId <tab> Label <tab> ScanNr <tab> feature1name <tab> ... <tab> featureNname <tab> Peptide <tab> Proteins

    if "all" in feature_sets:
        complete_df[
            ["SpecId", "Label", "ScanNr"]
            + ms2pip_feature_names
            + rt_feature_names
            + search_engine_feature_names
            + ["Peptide", "Proteins"]
        ].to_csv(
            "{}_allfeatures.pin".format(savepath), sep="\t", header=True, index=False
        )
    if "ms2pip_rt" in feature_sets:
        complete_df[
            ["SpecId", "Label", "ScanNr"]
            + ms2pip_feature_names
            + rt_feature_names
            + ["Peptide", "Proteins"]
        ].to_csv(
            "{}_ms2pip_rtfeatures.pin".format(savepath),
            sep="\t",
            header=True,
            index=False,
        )
    if "ms2pip" in feature_sets:
        complete_df[
            ["SpecId", "Label", "ScanNr"]
            + ms2pip_feature_names
            + ["Peptide", "Proteins"]
        ].to_csv(
            "{}_ms2pipfeatures.pin".format(savepath), sep="\t", header=True, index=False
        )
    if "rt" in feature_sets:
        complete_df[
            ["SpecId", "Label", "ScanNr"] + rt_feature_names + ["Peptide", "Proteins"]
        ].to_csv(
            "{}_rtfeatures.pin".format(savepath), sep="\t", header=True, index=False
        )
    if "searchengine" in feature_sets:
        complete_df[
            ["SpecId", "Label", "ScanNr"]
            + search_engine_feature_names
            + ["Peptide", "Proteins"]
        ].to_csv(
            "{}_searchenginefeatures.pin".format(savepath),
            sep="\t",
            header=True,
            index=False,
        )

    for features_version in feature_sets:
        redo_pin_tabs("{}_{}features.pin".format(savepath, features_version))
