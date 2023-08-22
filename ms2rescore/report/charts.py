import warnings
from collections import defaultdict
from typing import Dict, Optional, Tuple

import mokapot
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.subplots
import pyteomics.auxiliary
from psm_utils.psm_list import PSMList


def feature_ecdf_auc_bar(
    feature_ecdf_auc: pd.DataFrame, color_discrete_map: Optional[Dict[str, str]] = None
) -> go.Figure:
    """Plot bar chart of feature q-value ECDF AUCs."""
    return px.bar(
        data_frame=feature_ecdf_auc.sort_values("ecdf_auc", ascending=True),
        x="ecdf_auc",
        y="feature",
        color="feature_generator",
        orientation="h",
        hover_name="feature",
        labels={
            "ecdf_auc": "Q-value ECDF AUC",
            "feature_generator": "Feature generator",
            "feature": "Feature",
        },
        color_discrete_map=color_discrete_map,
    )


def feature_weights(
    feature_weights: pd.DataFrame, color_discrete_map: Optional[Dict[str, str]] = None
) -> go.Figure:
    """Plot bar chart of feature weights."""
    bar_data = (
        feature_weights.groupby(["feature", "feature_generator"])
        .median(numeric_only=True)
        .abs()
        .sort_values("weight")
        .reset_index()
    )

    return px.bar(
        data_frame=bar_data,
        x="weight",
        y="feature",
        color="feature_generator",
        orientation="h",
        hover_name="feature",
        title="Absolute median weights by feature",
        labels={
            "weight": "Absolute median weight",
            "feature_generator": "Feature generator",
            "feature": "Feature",
        },
        color_discrete_map=color_discrete_map,
    )


def feature_weights_by_generator(
    feature_weights: pd.DataFrame, color_discrete_map: Optional[Dict[str, str]] = None
) -> go.Figure:
    bar_data = (
        feature_weights.groupby(["feature", "feature_generator"])
        .median()
        .abs()
        .reset_index()
        .groupby("feature_generator")
        .sum(numeric_only=True)
        .reset_index()
        .sort_values("weight")
    )

    return px.bar(
        data_frame=bar_data,
        x="weight",
        y="feature_generator",
        color="feature_generator",
        orientation="h",
        hover_name="feature_generator",
        title="Absolute median weights, summed by feature generator",
        labels={
            "weight": "Absolute median weight",
            "feature_generator": "Feature generator",
            "feature": "Feature",
        },
        color_discrete_map=color_discrete_map,
    )


def calculate_feature_qvalues(
    features: pd.DataFrame,
    is_decoy: pd.Series,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate q-values and ECDF AUC for all rescoring features.

    Q-values are calculated for each feature as if it was directly used PSM score. For each q-value
    distribution, the ECDF AUC is calculated as a measure of overall individual performance of the
    feature.

    As it is not known whether higher or lower values are better for each feature, q-values are
    calculated for both the original and reversed scores. The q-values and ECDF AUC are returned
    for the calculation with the highest ECDF AUC.
    """
    feature_qvalues = dict()
    feature_ecdf_auc = dict()
    for fname in features:
        # Calculate q-values for reversed and non-reversed scores
        q_values = []
        for reverse in [False, True]:
            with warnings.catch_warnings():  # Ignore divide by zero warning
                warnings.simplefilter("ignore", category=RuntimeWarning)
                q_values.append(
                    pyteomics.auxiliary.qvalues(
                        features,
                        key=fname,
                        is_decoy=is_decoy,
                        remove_decoy=True,
                        reverse=reverse,
                        full_output=False,
                    )["q"]
                )

        # Calculate ECDF AUC as measure of overall individual performance of feature
        ecdf_aucs = [np.trapz(y=np.max(q) - np.sort(q)) for q in q_values]

        # Select and save q-value calculation with best AUC (score reversed or not)
        idx_best = np.argmax(ecdf_aucs)
        feature_qvalues[fname] = q_values[idx_best]
        feature_ecdf_auc[fname] = ecdf_aucs[idx_best]

    # Restructure as data frames
    feature_qvalues = pd.DataFrame(feature_qvalues)
    feature_ecdf_auc = (
        pd.DataFrame([feature_ecdf_auc])
        .transpose()
        .reset_index()
        .rename(columns={"index": "feature", 0: "ecdf_auc"})
    )

    return feature_qvalues, feature_ecdf_auc


def get_feature_values(psm_list: PSMList, feature_names: Optional[list] = None) -> pd.DataFrame:
    """Get feature values for all PSMs in a PSM list."""
    if not feature_names:
        feature_names = list(psm_list[0].rescoring_features.keys())
    features = pd.DataFrame(
        {fname: psm.rescoring_features[fname] for fname in feature_names} for psm in psm_list
    ).astype("float32")
    return features


class ECDF:
    """
    Return the Empirical CDF of an array as a step function.

    Parameters
    ----------
    x : array_like
        Observations
    """

    def __init__(self, x):
        # Get ECDF
        x = np.array(x, copy=True)
        x.sort()
        nobs = len(x)
        y = np.linspace(1.0 / nobs, 1, nobs)

        # Make into step function
        _x = np.asarray(x)
        _y = np.asarray(y)

        if _x.shape != _y.shape:
            msg = "x and y do not have the same shape"
            raise ValueError(msg)
        if len(_x.shape) != 1:
            msg = "x and y must be 1-dimensional"
            raise ValueError(msg)

        self.x = np.r_[-np.inf, _x]
        self.y = np.r_[0.0, _y]
        self.n = self.x.shape[0]

    def __call__(self, time):
        tind = np.searchsorted(self.x, time, side="right") - 1
        return self.y[tind]


def score_histogram(psm_df):
    psm_df = psm_df.copy()
    psm_df["is_decoy"] = psm_df["is_decoy"].map({True: "decoy", False: "target"})
    fig = px.histogram(
        psm_df,
        x="score",
        color="is_decoy",
        barmode="overlay",
        histnorm="",
        labels={"is_decoy": "PSM type", "False": "target", "True": "decoy"},
        opacity=0.5,
    )
    return fig


def pp_plot(psm_df):
    """Generate PP plot for given PSM dataframe."""
    n_decoys = np.count_nonzero(psm_df["is_decoy"])
    n_targets = len(psm_df) - n_decoys
    pi_zero = n_decoys / n_targets
    if n_decoys == 0:
        raise ValueError("No decoy PSMs found in PSM file.")
    target_scores = psm_df["score"][~psm_df["is_decoy"]]
    decoy_scores = psm_df["score"][psm_df["is_decoy"]]
    if len(psm_df) > 1000:
        target_scores_quantiles = psm_df["score"][~psm_df["is_decoy"]].quantile(
            np.linspace(0, 1, 1000)
        )
    else:
        target_scores_quantiles = target_scores
    target_ecdf = ECDF(target_scores)(target_scores_quantiles)
    decoy_ecdf = ECDF(decoy_scores)(target_scores_quantiles)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=decoy_ecdf,
            y=target_ecdf,
            mode="markers",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, pi_zero],
            mode="lines",
            line=go.scatter.Line(color="red"),
            showlegend=True,
            name="pi0",
        )
    )
    fig.update_layout(
        xaxis_title="Fdp",
        yaxis_title="Ftp",
        showlegend=False,
    )
    return fig


def fdr_plot(psm_df, fdr_threshold):
    """Plot number of identifications in function of FDR threshold."""
    df = (
        psm_df[~psm_df["is_decoy"]]
        .reset_index(drop=True)
        .sort_values("qvalue", ascending=True)
        .copy()
    )
    df["count"] = (~df["is_decoy"]).cumsum()
    fig = px.line(
        df,
        x="qvalue",
        y="count",
        log_x=True,
        labels={"count": "Number of identified target PSMs", "qvalue": "FDR threshold"},
    )
    fig.add_vline(x=fdr_threshold, line=go.layout.shape.Line(color="red"))
    return fig


def score_scatter_plot(
    before: mokapot.LinearConfidence,
    after: mokapot.LinearConfidence,
    level: str = "psms",
    indexer: str = "index",
):
    """Plot PSM scores before and after rescoring."""
    # Restructure data
    merge_columns = [indexer, "mokapot score", "mokapot q-value", "mokapot PEP"]
    ce_psms_targets = pd.merge(
        left=before.confidence_estimates[level],
        right=after.confidence_estimates[level][merge_columns],
        on=indexer,
        suffixes=(" before", " after"),
    )
    ce_psms_decoys = pd.merge(
        left=before.decoy_confidence_estimates[level],
        right=after.decoy_confidence_estimates[level][merge_columns],
        on=indexer,
        suffixes=(" before", " after"),
    )
    ce_psms_targets["PSM type"] = "target"
    ce_psms_decoys["PSM type"] = "decoy"
    ce_psms = pd.concat([ce_psms_targets, ce_psms_decoys], axis=0)

    # Get score thresholds
    score_threshold_before = (
        ce_psms[ce_psms["mokapot q-value before"] <= 0.01]
        .sort_values("mokapot q-value before", ascending=False)["mokapot score before"]
        .iloc[0]
    )
    score_threshold_after = (
        ce_psms[ce_psms["mokapot q-value after"] <= 0.01]
        .sort_values("mokapot q-value after", ascending=False)["mokapot score after"]
        .iloc[0]
    )

    # Plot
    fig = px.scatter(
        data_frame=ce_psms,
        x="mokapot score before",
        y="mokapot score after",
        color="PSM type",
        marginal_x="histogram",
        marginal_y="histogram",
        opacity=0.1,
        labels={
            "mokapot score before": "PSM score (before rescoring)",
            "mokapot score after": "PSM score (after rescoring)",
        },
    )
    # draw FDR thresholds
    fig.add_vline(x=score_threshold_before, line_dash="dash", row=1, col=1)
    fig.add_hline(y=score_threshold_after, line_dash="dash", row=1, col=1)
    fig.add_vline(x=score_threshold_before, line_dash="dash", row=2, col=1)
    fig.add_hline(y=score_threshold_after, line_dash="dash", row=1, col=2)

    return fig


def fdr_plot_comparison(
    before: mokapot.LinearConfidence,
    after: mokapot.LinearConfidence,
    level: str = "psms",
    indexer: str = "index",
) -> go.Figure:
    # Prepare data
    ce_psms_targets_melted = (
        pd.merge(
            left=before.confidence_estimates[level],
            right=after.confidence_estimates[level][
                [indexer, "mokapot score", "mokapot q-value", "mokapot PEP"]
            ],
            on=indexer,
            suffixes=(" before", " after"),
        )
        .rename(
            columns={
                "mokapot q-value before": "before rescoring",
                "mokapot q-value after": "after rescoring",
            }
        )
        .melt(
            id_vars=["index", "peptide", "is_target"],
            value_vars=["before rescoring", "after rescoring"],
            var_name="before/after",
            value_name="q-value",
        )
    )

    # Plot
    fig = px.ecdf(
        data_frame=ce_psms_targets_melted,
        x="q-value",
        color="before/after",
        log_x=True,
        ecdfnorm=None,
        labels={
            "q-value": "FDR threshold",
            "before/after": "",
        },
        color_discrete_map={
            "before rescoring": "#316395",
            "after rescoring": "#319545",
        },
    )
    fig.add_vline(x=0.01, line_dash="dash", line_color="black")
    fig.update_layout(yaxis_title="Identified PSMs")
    return fig


def identification_overlap(
    before: mokapot.LinearConfidence,
    after: mokapot.LinearConfidence,
) -> go.Figure():
    levels = before.levels  # ["psms", "peptides", "proteins"] if all available
    indexers = ["index", "index", "mokapot protein group"]

    overlap_data = defaultdict(dict)
    for level, indexer in zip(levels, indexers):
        df_before = before.confidence_estimates[level]
        df_after = after.confidence_estimates[level]
        set_before = set(df_before[df_before["mokapot q-value"] <= 0.01][indexer])
        set_after = set(df_after[df_after["mokapot q-value"] <= 0.01][indexer])

        overlap_data["removed"][level] = -len(set_before - set_after)
        overlap_data["shared"][level] = len(set_before | set_after)
        overlap_data["added"][level] = len(set_after - set_before)

    colors = ["#953331", "#316395", "#319545"]
    fig = plotly.subplots.make_subplots(rows=3, cols=1)

    for i, level in enumerate(levels):
        for (item, data), color in zip(overlap_data.items(), colors):
            fig.add_trace(
                go.Bar(
                    y=["protein groups" if level == "proteins" else level],
                    x=[data[level]],
                    marker={"color": color},
                    orientation="h",
                    name=item,
                    showlegend=True if i == 0 else False,
                ),
                row=i + 1,
                col=1,
            )
    fig.update_layout(barmode="relative")

    return fig
