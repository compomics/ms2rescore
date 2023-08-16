"""Mokapot integration for MS²Rescore."""

import logging
from typing import Any, List, Optional, Tuple

import mokapot
import numpy as np
import pandas as pd
import psm_utils
from mokapot.brew import brew
from mokapot.dataset import LinearPsmDataset
from pyteomics.mass import nist_mass

logger = logging.getLogger(__name__)

logging.getLogger("numba").setLevel(logging.WARNING)
logging.getLogger("mokapot.model").setLevel(logging.WARNING)
logging.getLogger("mokapot.dataset").setLevel(logging.WARNING)


def rescore(
    psm_list: psm_utils.PSMList,
    output_file_root: str = "ms2rescore",
    fasta_file: Optional[str] = None,
    write_weights: bool = False,
    write_txt: bool = False,
    write_flashlfq: bool = False,
    **kwargs: Any,
):
    """
    Rescore PSMs with Mokapot.

    Parameters
    ----------
    psm_list
        PSMs to be rescored.
    output_file_root
        Root of output file names. Defaults to ``"ms2rescore"``.
    fasta_file
        Path to FASTA file with protein sequences to use for protein inference. Defaults to
        ``None``.
    write_weights
        Write model weights to a text file. Defaults to ``False``.
    write_txt
        Write Mokapot results to a text file. Defaults to ``False``.
    write_flashlfq
        Write Mokapot results to a FlashLFQ-compatible file. Defaults to ``False``.
    **kwargs
        Additional keyword arguments are passed to the Mokapot ~:py:function:`mokapot.brew.brew`
        function.

    """
    # Convert PSMList to Mokapot dataset
    feature_names = list(psm_list[0].rescoring_features.keys())
    lin_psm_data = convert_psm_list(psm_list, feature_names)

    # Add proteins
    if fasta_file:
        proteins = mokapot.read_fasta(fasta_file)
        lin_psm_data.add_proteins(proteins)

    # Rescore
    confidence_results, models = brew(lin_psm_data, **kwargs)

    # Reshape confidence estimates to match PSMList
    mokapot_values_targets = (
        confidence_results.confidence_estimates["psms"]
        .set_index("index")
        .sort_index()[["mokapot score", "mokapot q-value", "mokapot PEP"]]
    )
    mokapot_values_decoys = (
        confidence_results.decoy_confidence_estimates["psms"]
        .set_index("index")
        .sort_index()[["mokapot score", "mokapot q-value", "mokapot PEP"]]
    )
    q = np.full((len(psm_list), 3), np.nan)
    q[mokapot_values_targets.index] = mokapot_values_targets.values
    q[mokapot_values_decoys.index] = mokapot_values_decoys.values

    # Add Mokapot results to PSMList
    psm_list["score"] = q[:, 0]
    psm_list["qvalue"] = q[:, 1]
    psm_list["pep"] = q[:, 2]

    # Write results
    if write_weights:
        save_model_weights(models, feature_names, output_file_root)
    if write_txt:
        confidence_results.to_txt(file_root=output_file_root, decoys=True)
    if write_flashlfq:
        confidence_results.to_flashlfq(output_file_root + ".mokapot.flashlfq.txt")


def convert_psm_list(
    psm_list: psm_utils.PSMList,
    feature_names: List[str],
) -> LinearPsmDataset:
    """
    Convert a PSM list to a Mokapot dataset.

    Parameters
    ----------
    psm_list
        PSMList to rescore.
    feature_names
        List of feature names to use. Items must be keys in the PSM `rescoring_features` dict.

    Returns
    -------
    mokapot.dataset.LinearPsmDataset

    """

    psm_df = psm_list.to_dataframe()
    psm_df = psm_df.reset_index(drop=True).reset_index()

    psm_df["peptide"] = (
        psm_df["peptidoform"].astype(str).str.replace(r"(/\d+$)", "", n=1, regex=True)
    )
    psm_df["is_target"] = ~psm_df["is_decoy"]
    psm_df["charge"] = psm_df["peptidoform"].apply(lambda x: x.precursor_charge)
    psm_df["calcmass"] = psm_df["peptidoform"].apply(lambda x: x.theoretical_mass)
    psm_df["expmass"] = _mz_to_mass(psm_df["precursor_mz"], psm_df["charge"])
    psm_df = pd.concat(
        [
            psm_df.drop(columns=["rescoring_features"]),
            pd.DataFrame(list(psm_df["rescoring_features"])).fillna(0.0),
        ],
        axis=1,
    )

    required_columns = [
        "index",
        "peptide",
        "is_target",
        "protein_list",
        "run",
        "calcmass",
        "expmass",
        "retention_time",
        "charge",
    ]
    required_columns.extend(feature_names)

    lin_psm_data = LinearPsmDataset(
        psms=psm_df[required_columns],
        target_column="is_target",
        spectrum_columns="index",  # Use artificial index to allow multi-rank rescoring
        peptide_column="peptide",
        protein_column="protein_list",
        feature_columns=feature_names,
        filename_column="run",
        scan_column="index",  # Keep as spectrum_id?
        calcmass_column="calcmass",
        expmass_column="expmass",
        rt_column="retention_time",
        charge_column="charge",
    )

    return lin_psm_data


def save_model_weights(
    models: Tuple[mokapot.model.Model], feature_names: List[str], output_file_root: str
):
    """
    Save model weights to a file.

    Parameters
    ----------
    models
        Tuple of Mokapot models (one for each fold) to save.
    feature_names
        List of feature names that were used to train the models.
    output_file_root
        Root of output file names.

    """
    pd.DataFrame(
        np.stack([m.estimator.coef_[0] for m in models]),
        columns=list(feature_names),
    ).to_csv(output_file_root + ".mokapot.weights.tsv", sep="\t", index=False)


def _mz_to_mass(mz: float, charge: int) -> float:
    """Convert m/z to mass."""
    return mz * charge - charge * nist_mass["H"][1][0]
