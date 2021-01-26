"""Plot MS²ReScore results."""

from typing import List, Optional

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pyteomics.auxiliary import qvalues
from statsmodels.distributions.empirical_distribution import ECDF

from ms2rescore.percolator import PercolatorIn


def target_decoy_distribution(
    df: pd.DataFrame,
    q_label: str = "q",
    score_label: str = "score",
    decoy_label: str = "is decoy",
    score_name: str = "Score",
    fdr_threshold: Optional[float] = None,
    plot_title: str = ""
):
    """
    Plot target-decoy distributions.

    Plot for a given search engine output the target-decoy score distributions, the
    relation between the q-values and the PSM scores and a PP plot between the target
    and the decoy distribution.

    Parameters
    ----------
    df : pd.DataFrame
        Pandas DataFrame containing all rank 1 PSMs
    q_label : str, optional
        Label of column containing each PSM's q-value, default `q`
    score_label : str, optional
        Label of column containing the search engine scores, default `score`
    decoy_label : str, optional
        Label of column marking decoy PSMs as True and target PSMs as False, default
        `is decoy`
    score_name : str, optional
        Score name used in axis labels
    fdr_threshold : float, optional
        FDR threshold to plot as vertical dotted line
    plot_title : str, optional
        plot suptitle

    Returns
    -------
    fig : Matplotlib figure
    axes : Matplotlib Axes

    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    if fdr_threshold:
        score_cutoff = (
            df[(df[q_label] <= 0.01) & (~df[decoy_label])]
            .sort_values(q_label)
            .iloc[-1][score_label]
        )
    else:
        score_cutoff = None

    # Score distribution plot
    plot_list = [
        list(x)
        for x in [df[df[decoy_label]][score_label], df[~df[decoy_label]][score_label]]
    ]
    axes[0].hist(
        plot_list,
        bins=24,
        label=["Decoy", "Target"],
        color=["r", "blue"],
        lw=1,
        rwidth=1,
    )
    if fdr_threshold:
        axes[0].vlines(
            x=score_cutoff, ymin=0, ymax=axes[0].get_ylim()[1], linestyles="dashed"
        )
    axes[0].legend()
    axes[0].set_ylabel("Number of matches")
    axes[0].set_xlabel(score_name)

    # Q value plot
    axes[1].plot(
        df.sort_values(score_label)[score_label], df.sort_values(score_label)[q_label]
    )
    if fdr_threshold:
        axes[1].vlines(
            x=score_cutoff, ymin=0, ymax=axes[1].get_ylim()[1], linestyles="dashed"
        )
    axes[1].set_ylabel("q-value")
    axes[1].set_xlabel(score_name)

    # PP plot
    ratio = df[decoy_label].value_counts()[True] / df[decoy_label].value_counts()[False]
    Ft = ECDF(df[~df[decoy_label]][score_label])
    Fd = ECDF(df[df[decoy_label]][score_label])
    x = df[~df[decoy_label]][score_label]
    Fdp = Fd(x)
    Ftp = Ft(x)
    axes[2].scatter(Fdp, Ftp, s=4)
    axes[2].plot((0, 1), (0, ratio), color="r")
    axes[2].set_xlabel("Decoy percentile")
    axes[2].set_ylabel("Target percentile")

    plt.suptitle(plot_title)
    sns.despine()

    return fig, axes


def qvalue_comparison(
    datasets: List[pd.DataFrame],
    dataset_labels: Optional[List[str]] = None,
    q_label: str = "q",
    decoy_label: str = "is decoy",
    fdr_thresholds: Optional[List[float]] = None,
    log_scale: bool = True,
    title: str = "",
    ax: Optional[matplotlib.axes.Axes] = None,
):
    """
    Plot identification count in function of q-value threshold for multiple datasets.

    Parameters
    ----------
    datasets : List[pd.DataFrame]
        list of datasets in the form of `pandas.DataFrame`s
    dataset_labels : List[str]
        list of dataset labels to use in figure legend
    q_label : str, optional
        label of column containing each PSM's q-value, default `q`
    decoy_label : str, optional
        label of column marking decoy PSMs as True and target PSMs as False, default
        `is decoy`
    fdr_thresholds : List[float], optional
        list of FDR thresholds to plot as vertical, dotted lines
    log_scale : bool
        plot x-axis (q-values) in log scale or not.
    ax : matplotlib Axes, optional
        axes object to draw the plot onto, otherwise uses the current axes.

    Returns
    -------
    ax : matplotlib Axes


    """
    if not isinstance(datasets, list):
        raise TypeError("`datasets` should be of type `list`.")
    if not datasets:
        raise ValueError("`datasets` cannot be empty.")

    if not fdr_thresholds:
        fdr_thresholds = [0.01, 0.001]

    if ax is None:
        ax = plt.gca()

    max_count = 0

    for i, df_in in enumerate(datasets):
        # Cumulatively count target IDs at each FDR threshold
        df = (
            df_in[~df_in[decoy_label]]
            .reset_index(drop=True)
            .sort_values(q_label, ascending=True)
            .copy()
        )
        df["count"] = (~df[decoy_label]).cumsum()

        # Plot counts
        label = dataset_labels[i] if dataset_labels else None
        ax.plot(df[q_label], df["count"], label=label, alpha=0.5)

        # Get maximum count, required for vertical line height
        tmp_max = np.max(df["count"])
        if tmp_max > max_count:
            max_count = tmp_max

    # Plot FDR thresholds as dotted lines
    if fdr_thresholds:
        for fdr in fdr_thresholds:
            ax.plot(
                [fdr] * 2, np.linspace(0, max_count, 2), linestyle="--", color="black",
            )

    # Figure labels and legend
    ax.set_xlim(0.00001, 1)
    ax.set_ylabel("Number of identified spectra")
    ax.set_title(title)
    if log_scale:
        ax.set_xlabel("FDR threshold (log scale)")
        ax.set_xscale("log")
    else:
        ax.set_xlabel("FDR threshold")
    if dataset_labels:
        ax.legend()

    return ax


def plot_rescoring_results(
    path_to_pin,
    path_to_target_pout,
    path_to_decoy_pout
):

    def _read_pin_file(path_to_pin):
        pin = PercolatorIn("Linfeng_010411_HapMap31.pin")
        pin_qvalues = pd.DataFrame(qvalues(
            pin.df,
            key=pin.df["lnEValue"],
            is_decoy=pin.df["Label"] == -1,
            reverse=True,
            remove_decoy=False,
            formula=1
        ))
        return pin_qvalues

    def _read_pout_file(path_to_target_pout, path_to_decoy_pout):
        """Read target and decoy pout files and combine into single pandas DataFrame."""
        def _read_pout(path):
            pout = pd.read_csv(
                PercolatorIn.fix_tabs(path, id_column="PSMId"),
                sep="\t"
            )
            return pout

        target_pout = _read_pout(path_to_target_pout)
        decoy_pout = _read_pout(path_to_decoy_pout)
        target_pout["Label"] = 1
        decoy_pout["Label"] = -1
        pout = pd.concat([target_pout, decoy_pout])

        pout_qvalues = pout[["score", "q-value", "Label"]].rename(
            columns={"q-value": "q", "Label": "is decoy"}
        )
        pout_qvalues["is decoy"] = pout["Label"] == -1

        return pout_qvalues

    pin_qvalues = _read_pin_file(path_to_pin)
    pout_qvalues = _read_pout_file(path_to_target_pout, path_to_decoy_pout)

    fig, axes = target_decoy_distribution(
        pin_qvalues, plot_title="Before rescoring"
    )
    plt.savefig("target_decoy_dist_before_rescoring.svg")

    fig, axes = target_decoy_distribution(
        pout_qvalues, plot_title="After rescoring"
    )
    plt.savefig("target_decoy_dist_after_rescoring.svg")

    ax = qvalue_comparison(
        datasets=[pin_qvalues, pout_qvalues],
        dataset_labels=["Before rescoring", "After rescoring"],
    )
    plt.savefig("qvalue_comparison.svg")
