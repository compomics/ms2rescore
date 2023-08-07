"""Plot MS²Rescore results."""

import logging
import os
from abc import ABC
from collections import defaultdict
from typing import List, Optional

import click
import matplotlib.axes
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from pyteomics.auxiliary import qvalues

from ms2rescore.exceptions import MS2RescoreError
from ms2rescore.rescoring_engines.percolator import PercolatorIn

sns.set_style("whitegrid")
logger = logging.getLogger(__name__)

MS2PIP_FEATURES = [
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

DEEPLC_FEATURES = [
    "observed_retention_time",
    "predicted_retention_time",
    "rt_diff",
    "rt_diff_best",
    "observed_retention_time_best",
    "predicted_retention_time_best",
]


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


class RescoreRecord(ABC):
    # TODO remove support for different samples
    rerecs = []
    weights = None
    loss_gain_df = pd.DataFrame()
    unique_df = pd.DataFrame()
    count_df = pd.DataFrame()

    def __init__(self) -> None:
        pass

    @classmethod
    def empty_rerecs(cls):
        cls.rerecs = []
        cls.weights = None

    @classmethod
    def show_rerec_items(cls):
        for rerec in cls.rerecs:
            print(rerec)
        if cls.weights:
            print(cls.weights)

    def target_decoy_distribution(
        self,
        q_label: str = "q",
        score_label: str = "score",
        decoy_label: str = "is decoy",
        score_name: str = "Score",
        fdr_threshold: Optional[float] = None,
        plot_title: str = "",
    ):
        """
        Plot target-decoy distributions.

        Plot for a given search engine output the target-decoy score distributions, the
        relation between the q-values and the PSM scores and a PP plot between the
        target and the decoy distribution.

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
                self.df[(self.df[q_label] <= 0.01) & (~self.df[decoy_label])]
                .sort_values(q_label)
                .iloc[-1][score_label]
            )
        else:
            score_cutoff = None

        # Score distribution plot
        plot_list = [
            list(self.df[self.df[decoy_label]][score_label]),
            list(self.df[~self.df[decoy_label]][score_label]),
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
            axes[0].vlines(x=score_cutoff, ymin=0, ymax=axes[0].get_ylim()[1], linestyles="dashed")
        axes[0].legend()
        axes[0].set_ylabel("Number of matches")
        axes[0].set_xlabel(score_name)

        # Q value plot
        axes[1].plot(
            self.df.sort_values(score_label)[score_label],
            self.df.sort_values(score_label)[q_label],
        )
        if fdr_threshold:
            axes[1].vlines(x=score_cutoff, ymin=0, ymax=axes[1].get_ylim()[1], linestyles="dashed")
        axes[1].set_ylabel("q-value")
        axes[1].set_xlabel(score_name)

        # PP plot
        ratio = (
            self.df[decoy_label].value_counts()[True] / self.df[decoy_label].value_counts()[False]
        )
        Ft = ECDF(self.df[~self.df[decoy_label]][score_label])
        Fd = ECDF(self.df[self.df[decoy_label]][score_label])
        x = self.df[~self.df[decoy_label]][score_label]
        Fdp = Fd(x)
        Ftp = Ft(x)
        axes[2].scatter(Fdp, Ftp, s=4)
        axes[2].plot((0, 1), (0, ratio), color="r")
        axes[2].set_xlabel("Decoy percentile")
        axes[2].set_ylabel("Target percentile")

        plt.suptitle(plot_title)
        sns.despine()

        return fig, axes

    @classmethod
    def qvalue_comparison(
        cls,
        dataset_labels: Optional[List[str]] = None,
        q_label: str = "q",
        decoy_label: str = "is decoy",
        fdr_thresholds: Optional[List[float]] = None,
        log_scale: bool = True,
        title: str = "",
        ax: Optional[matplotlib.axes.Axes] = None,
    ):
        """
        Plot identification count in function of q-value threshold for multiple
        datasets.

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
        fig = plt.figure()
        if not isinstance(cls.rerecs, list):
            raise TypeError("`datasets` should be of type `list`.")
        if not cls.rerecs:
            raise ValueError("`datasets` cannot be empty.")

        if not fdr_thresholds:
            fdr_thresholds = [0.01, 0.001]
        else:
            fdr_thresholds = list(set(fdr_thresholds + [0.01, 0.001]))
            fdr_thresholds.sort(reverse=True)

        if ax is None:
            ax = plt.gca()

        max_count = 0

        for i, rerec in enumerate(cls.rerecs):
            # Cumulatively count target IDs at each FDR threshold
            df = (
                rerec.df[~rerec.df[decoy_label]]
                .reset_index(drop=True)
                .sort_values(q_label, ascending=True)
                .copy()
            )
            df["count"] = (~df[decoy_label]).cumsum()

            # Plot counts
            if dataset_labels:
                label = dataset_labels[i]
            else:
                label = rerec.rescore_features
            ax.plot(df[q_label], df["count"], label=label, alpha=0.5)

            # Get maximum count, required for vertical line height
            tmp_max = np.max(df["count"])
            if tmp_max > max_count:
                max_count = tmp_max

        # Plot FDR thresholds as dotted lines
        if fdr_thresholds:
            for fdr in fdr_thresholds:
                ax.plot(
                    [fdr] * 2,
                    np.linspace(0, max_count, 2),
                    linestyle="--",
                    color="black",
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

        ax.legend(frameon=True, ncol=4, loc=9)
        fig.set_size_inches(12, 10)
        fig.tight_layout()
        return ax

    @classmethod
    def _count_identifications(cls, FDR_threshold=[0.01]):
        """
        Count the number of identifications for each PIN, POUT that is listed
        in the rescore record for certain FDR thresholds.

        Parameters
        ----------
        FDR_threshold : List[floats]
            list of FDR threshold for which to count number of identifications

        """
        identification_counts = []
        for rerec in cls.rerecs:
            for threshold in FDR_threshold:
                tmp = []
                tmp.append(rerec.name)
                tmp.append(rerec.rescore_features)
                tmp.append(sum(rerec.df.q[~rerec.df["is decoy"]] < threshold))
                tmp.append(threshold)
                identification_counts.append(tmp)
        cls.count_df = pd.DataFrame(
            identification_counts, columns=["sample", "rescoring", "count", "FDR"]
        )

    @classmethod
    def _separate_unique_peptides(cls, FDR_threshold=[0.01]):
        """
        Count the number of unique identifications for each PIN, POUT that is listed in
        the rescore record for certain FDR thresholds.

        Parameters
        ----------
        FDR_threshold : List[floats]
            list of FDR threshold for which to count number of unique identifications

        """
        unique_samples = []
        for rerec in cls.rerecs:
            for threshold in FDR_threshold:
                tmp = []
                tmp.append(rerec.name)
                tmp.append(rerec.rescore_features)
                tmp.append(
                    rerec.df["peptide"][
                        ~(rerec.df["is decoy"]) & (rerec.df["q"] < threshold)
                    ].unique()
                )
                tmp.append(threshold)
                unique_samples.append(tmp)
        cls.unique_df = pd.DataFrame(
            unique_samples, columns=["sample", "rescoring", "upeps", "FDR"]
        )
        cls.unique_df["count"] = cls.unique_df["upeps"].apply(len)

    @classmethod
    def count_plot(cls, unique=False):
        """
        Barplot of the number of (unique) identifications for each FDR threshold.

        Parameters
        ----------
        unique : Boolean
            If true only plot unique identifications
            If false plot total amount of identifications
        """
        if unique:
            if cls.unique_df.empty:
                cls._separate_unique_peptides()
            y_label = "number of unique identified peptides"
            count_df = cls.unique_df
        else:
            if cls.count_df.empty:
                cls._count_identifications()
            y_label = "number of identified PSMs"
            count_df = cls.count_df

        g = sns.catplot(
            x="FDR",
            y="count",
            data=count_df,
            hue="rescoring",
            kind="bar",
            # col="sample",
            legend=False,
        )

        g.set_ylabels(y_label)
        g.add_legend(loc=9, ncol=4)
        g.fig.set_size_inches(12, 10)
        g.fig.tight_layout()

        return g

    @classmethod
    def calculate_loss_gain_df(cls, reference=None, FDR_threshold=[0.01]):
        """
        Calculate relative amount of gained, lossed and shared unique peptides for each
        rescoring method relative to the reference.

        Parameters
        ----------
        reference : string
            The rescoring method be used as a reference

        FDR_threshold : List[floats]
            list of FDR threshold for which to count relative gain, loss

        """
        if cls.unique_df.empty:
            cls._separate_unique_peptides()
        loss_gain = defaultdict(list)
        for sample in cls.unique_df["sample"].unique():
            for threshold in FDR_threshold:
                ft_dict = {}
                for feature in cls.unique_df["rescoring"][
                    cls.unique_df["sample"] == sample
                ].unique():
                    ft_dict[feature] = set(
                        list(
                            cls.unique_df["upeps"][
                                (cls.unique_df["sample"] == sample)
                                & (cls.unique_df["rescoring"] == feature)
                                & (cls.unique_df["FDR"] == threshold)
                            ].item()
                        )
                    )
                if reference:
                    pass
                elif "After rescoring: Searchengine" in ft_dict.keys():
                    reference = "After rescoring: Searchengine"
                else:
                    reference = "Before rescoring"

                total = len(ft_dict[reference])
                if total == 0:
                    continue
                for feature in cls.unique_df["rescoring"][
                    cls.unique_df["sample"] == sample
                ].unique():
                    loss_gain["sample"].append(sample)
                    loss_gain["feature"].append(feature)
                    loss_gain["FDR"].append(threshold)
                    loss_gain["shared"].append(
                        (len(ft_dict[reference] & ft_dict[feature]) / total) * 100
                    )
                    loss_gain["gain"].append(
                        (len(ft_dict[feature] - ft_dict[reference]) / total) * 100
                    )
                    loss_gain["loss"].append(
                        (len(ft_dict[reference] - ft_dict[feature]) / total) * -100
                    )

        cls.loss_gain_df = pd.DataFrame(loss_gain)
        cls.loss_gain_df["gain"] = cls.loss_gain_df["gain"] + cls.loss_gain_df["shared"]

    @classmethod
    def loss_gain_plot(cls, FDR):
        """
        Barplot of the relative count of unique peptides for a fixed FDR
        for each rescoring feature set.

        Parameters
        ----------
        FDR : float
            single FDR threshold used for the plot
        """
        fig = plt.figure()

        if cls.loss_gain_df.empty:
            cls.calculate_loss_gain_df(FDR_threshold=[FDR])
        if FDR not in cls.loss_gain_df["FDR"].unique():
            cls._separate_unique_peptides(FDR_threshold=[FDR])
            cls.calculate_loss_gain_df(FDR_threshold=[FDR])
        tmp = cls.loss_gain_df[cls.loss_gain_df["FDR"] == FDR]

        samples = tmp["sample"].unique()
        number_samples = len(samples)

        for i, sample in enumerate(samples, start=1):
            if i == number_samples:
                ax = fig.add_subplot(int(f"{number_samples}1{i}"), frameon=False)
                ax.title.set_text(f"FDR={FDR}")
                sns.barplot(
                    y="feature",
                    x="gain",
                    data=tmp[tmp["sample"] == sample],
                    color="#2FA92D",
                    order=tmp[tmp["sample"] == sample].sort_values("gain").feature,
                    ax=ax,
                )
                sns.barplot(
                    y="feature",
                    x="shared",
                    data=tmp[tmp["sample"] == sample],
                    color="#1AA3FF",
                    order=tmp[tmp["sample"] == sample].sort_values("gain").feature,
                    ax=ax,
                )
                sns.barplot(
                    y="feature",
                    x="loss",
                    data=tmp[tmp["sample"] == sample],
                    color="#FF0000",
                    order=tmp[tmp["sample"] == sample].sort_values("gain").feature,
                    ax=ax,
                )
                ax.axes.set_xlabel("unique identified peptides (%)")
                ax.axes.set_ylabel("")
            else:
                ax = fig.add_subplot(int(f"{len(samples)}1{i}"), frameon=False)
                sns.barplot(
                    y="feature",
                    x="gain",
                    data=tmp[tmp["sample"] == sample],
                    color="#2FA92D",
                    order=tmp[tmp["sample"] == sample].sort_values("gain").feature,
                    ax=ax,
                )
                sns.barplot(
                    y="feature",
                    x="shared",
                    data=tmp[tmp["sample"] == sample],
                    color="#1AA3FF",
                    order=tmp[tmp["sample"] == sample].sort_values("gain").feature,
                    ax=ax,
                )
                sns.barplot(
                    y="feature",
                    x="loss",
                    data=tmp[tmp["sample"] == sample],
                    color="#FF0000",
                    order=tmp[tmp["sample"] == sample].sort_values("gain").feature,
                    ax=ax,
                )
                ax.axes.set_xlabel("")
                ax.axes.set_ylabel("")
        fig.set_size_inches(12, 6)
        return ax

    @classmethod
    def save_plots_to_pdf(cls, filename, FDR_thresholds=[0.01]):
        """
        Create pdf file with rescore plots.

        Parameters
        ----------
        filename : string
            filename for the pdf file

        FDR_thresholds : list[float]
            list of FDR thresholds for wich to create plots

        """
        if not cls.rerecs:
            raise MS2RescoreError("no pin/pout files listed")

        pdf = matplotlib.backends.backend_pdf.PdfPages(filename)
        cls._separate_unique_peptides(FDR_thresholds)
        cls._count_identifications(FDR_thresholds)

        cls.count_plot(unique=False)
        pdf.savefig()

        cls.qvalue_comparison(fdr_thresholds=FDR_thresholds)
        plt.tight_layout()
        pdf.savefig()

        for fdr in FDR_thresholds:
            cls.loss_gain_plot(FDR=fdr)
            plt.tight_layout()
            pdf.savefig()

        if cls.weights:
            cls.plot_rescore_weights()
            pdf.savefig()

        pdf.close()

    @classmethod
    def plot_rescore_weights(cls):
        """Plot weights of rescoring run"""

        fig = plt.figure(figsize=(16, 6))

        ax = plt.subplot2grid((1, 15), (0, 0), colspan=13, fig=fig)
        cls.weights.plot_individual_feature_weights(ax)

        ax1 = plt.subplot2grid((1, 15), (0, 13), colspan=1, fig=fig)
        cls.weights.plot_feature_set_weights(ax1)

        fig.tight_layout()


class PIN(RescoreRecord):
    """PIN file record."""

    def __init__(self, path_to_file, sample_name, score_metric=None) -> None:
        super().__init__()
        self.type = "pin"
        self.rescore_features = "Before rescoring"
        self.name = sample_name
        if self._infer_filetype(path_to_file) == "pin":
            self.score_metric = score_metric
            self.df = self._read_pin_file(path_to_file)
        elif self._infer_filetype(path_to_file) == "peprec":
            self.score_metric = "psm_score"
            self.df = self._read_pin_from_peprec(path_to_file)
        else:
            raise MS2RescoreError("Not a peprec or pin file")

        self.rerecs.append(self)  # add PIN record to RescoreRecord

        """
        Parameters
        ----------
        path_to_file : string
            path to input file, either peprec or pin file

        sample_name : string
            name of the sample

        score_metric: string
            score column label (only required when providing pin file)
        """

    def __str__(self) -> str:
        return (
            f"Sample name: {self.name}\nType: {self.type}\nRescore status: {self.rescore_features}"
        )

    def _read_pin_file(self, path_to_pin):
        """Read pin file, calculate qvalues and write into single pandas DataFrame."""
        pin = PercolatorIn(path_to_pin)
        pin_qvalues = pd.DataFrame(
            qvalues(
                pin.df,
                key=pin.df[self.score_metric],
                is_decoy=pin.df["Label"] == -1,
                reverse=True,
                remove_decoy=False,
                formula=1,
                full_output=True,
            )
        )
        return pin_qvalues[["SpecId", "is decoy", "score", "q", "Peptide"]].rename(
            columns={"SpecId": "PSMId", "Peptide": "peptide"}
        )

    def _read_pin_from_peprec(self, path_to_peprec):
        peprec = pd.read_table(path_to_peprec, sep=" ")
        pin_qvalues = pd.DataFrame(
            qvalues(
                peprec,
                key=peprec["psm_score"],
                is_decoy=peprec["Label"] == -1,
                reverse=True,
                remove_decoy=False,
                formula=1,
                full_output=True,
            )
        )
        return pin_qvalues[["spec_id", "is decoy", "score", "q", "peptide"]].rename(
            columns={"spec_id": "PSMId"}
        )

    def _infer_filetype(self, filepath):
        """
        Extract filetype from file

        Parameters
        ----------
        filepath : string
            filepath from file

        """
        filetype = filepath.rsplit(".", 1)[1]
        return filetype


class POUT(RescoreRecord):
    """POUT file record."""

    def __init__(
        self, path_to_target_pout, path_to_decoy_pout, sample_name, rescoring_features
    ) -> None:
        super().__init__()
        self.target_pout = path_to_target_pout
        self.decoy_pout = path_to_decoy_pout
        self.df = self._read_pout_file
        self.type = "pout"
        self.rescore_features = "After rescoring: " + rescoring_features
        self.name = sample_name
        self.df = self._read_pout_file()
        self.rerecs.append(self)  # add POUT record to RescoreRecord

        """
        Parameters
        ----------
        path_to_target_pout : string
            path to target pout file

        path_to_decoy_pout : string
            path to decoy pout file

        sample_name : string
            name of the sample

        rescoring_features: list[strings]
            list of the used rescoring features
        """

    def __str__(self) -> str:
        return (
            f"Sample name: {self.name}\nType: {self.type}\nRescore status: {self.rescore_features}"
        )

    @staticmethod
    def _read_pout(path):
        pout = pd.read_csv(PercolatorIn.fix_tabs(path, id_column="PSMId"), sep="\t")
        return pout

    def _read_pout_file(self):
        """Read target and decoy pout files and combine into single pandas DataFrame."""

        target_pout = self._read_pout(self.target_pout)
        decoy_pout = self._read_pout(self.decoy_pout)
        target_pout["Label"] = 1
        decoy_pout["Label"] = -1
        pout = pd.concat([target_pout, decoy_pout])

        pout_qvalues = pout[["PSMId", "score", "q-value", "Label", "peptide"]].rename(
            columns={"q-value": "q", "Label": "is decoy"}
        )
        pout_qvalues["is decoy"] = pout["Label"] == -1

        return pout_qvalues


class PERCWEIGHT(RescoreRecord):
    """Percolator weights file record"""

    def __init__(
        self,
        path_to_weights_file,
        rescoring_features,
        sample_name,
        normalized_weights=True,
    ) -> None:
        super().__init__()
        self.type = "weights"
        self.rescore_features = rescoring_features
        self.name = sample_name
        self.weights_df = self.read_weights_file(path_to_weights_file, normalized_weights)
        RescoreRecord.weights = self

        """
        Parameters
        ----------
        path_to_file : string
            path to input percolator weights file

        sample_name : string
            name of the sample

        rescoring_features: list[strings]
            list of the used rescoring features
        """

    def __str__(self) -> str:
        return (
            f"Sample name: {self.name}\nType: {self.type}\nRescore status: {self.rescore_features}"
        )

    @staticmethod
    def _check_file_type(filepath):
        """Check if file exists and is a percolator weights file"""

        if os.path.exists(filepath) & filepath.endswith(".weights"):
            return True

        return False

    @property
    def _get_se_features(self):
        return [
            ft for ft in self.weights_df.columns if ft not in MS2PIP_FEATURES + DEEPLC_FEATURES
        ]

    def read_weights_file(self, filename, use_norm_weights: bool):
        """Read weights file to dataframe"""

        if not self._check_file_type(filename):
            raise FileNotFoundError(
                "Not a valid filepath or not a percolator weights file (.weights)"
            )

        if use_norm_weights:
            remove_rows = [1, 2, 4, 5, 7]
        else:
            remove_rows = [0, 2, 3, 5, 6]

        weights_df = pd.read_table(filename, sep="\t")
        weights_df.drop(remove_rows, axis=0, inplace=True)
        weights_df = weights_df.drop("m0", axis=1).astype(float)
        weights_df.loc["mean"] = weights_df.mean()

        return weights_df

    def calculate_precentage_ft_weights(self):
        """Return a dict with the percentage of weight for each feature set"""
        feature_weights = {}
        total_weight = sum(self.weights_df.loc["mean", :].abs())
        for name, feature_set in zip(
            ["MS²PIP", "DeepLC", "Search engine"],
            [MS2PIP_FEATURES, DEEPLC_FEATURES, self._get_se_features],
        ):
            try:
                feature_weights[name] = (
                    sum(np.abs([self.weights_df.loc["mean", feature] for feature in feature_set]))
                    / total_weight
                    * 100
                )
            except KeyError:
                feature_weights[name] = 0
                continue

        return feature_weights

    def plot_feature_set_weights(self, ax=None):
        """Plot the total weigths as percentages for the different features sets"""

        if not ax:
            ax = plt.gca()

        ft_weights = pd.DataFrame(self.calculate_precentage_ft_weights(), index=[0])
        ft_weights["Search engine"] = ft_weights["Search engine"] + ft_weights["MS²PIP"]
        ft_weights["DeepLC"] = ft_weights["Search engine"] + ft_weights["DeepLC"]
        sns.barplot(
            y="DeepLC",
            data=ft_weights,
            palette=sns.color_palette(["#28ea22"]),
            ax=ax,
            label="DeepLC",
        )
        sns.barplot(
            y="Search engine",
            data=ft_weights,
            palette=sns.color_palette(["#FFCD27"]),
            ax=ax,
            label="Search engine",
        )
        sns.barplot(
            y="MS²PIP",
            data=ft_weights,
            palette=sns.color_palette(["#1AA3FF"]),
            ax=ax,
            label="ms2pip",
        )
        sns.despine(left=True, right=True, top=True, ax=ax)
        ax.set_xlabel("Percolator weights (%)")
        ax.set_ylabel("")
        return ax

    def plot_individual_feature_weights(self, ax=None, absolute=True):
        """Plot the individual feature weights"""

        if not ax:
            ax = plt.gca()

        reindex_list = []
        for feature_set in [MS2PIP_FEATURES, DEEPLC_FEATURES, self._get_se_features]:
            try:
                feature_reindex = list(
                    self.weights_df.loc["mean", feature_set]
                    .abs()
                    .sort_values(ascending=True)
                    .index
                )
                reindex_list.extend(feature_reindex)
            except KeyError:
                continue

        if absolute:
            mean_row = self.weights_df.loc["mean", :].abs().reindex(reindex_list)
        else:
            mean_row = self.weights_df.loc["mean", :].reindex(reindex_list)

        color_map = self._get_color_map()
        feature_cmap = [color_map[ft] for ft in mean_row.index]

        mean_row.plot(kind="bar", color=feature_cmap, ax=ax)

        ms2pip_l = Patch(facecolor="#1AA3FF", edgecolor="#1AA3FF")
        deeplc_l = Patch(facecolor="#28ea22", edgecolor="#28ea22")
        searchengine_l = Patch(facecolor="#FFCD27", edgecolor="#FFCD27")

        ax.legend(
            [ms2pip_l, deeplc_l, searchengine_l],
            ["MS²PIP", "DeepLC", "Search engine"],
            frameon=True,
            ncol=3,
            loc=2,
        )

        return ax

    def _get_color_map(self, feature_colors: dict = None):
        """Get color for each feature"""

        if not feature_colors:
            feature_colors = {
                "MS²PIP": "#1AA3FF",
                "DeepLC": "#28ea22",
                "Search engine": "#FFCD27",
            }
        color_mapping = {}
        for feature in self.weights_df.columns:
            if feature in MS2PIP_FEATURES:
                color_mapping[feature] = feature_colors["MS²PIP"]
            elif feature in DEEPLC_FEATURES:
                color_mapping[feature] = feature_colors["DeepLC"]
            else:
                color_mapping[feature] = feature_colors["Search engine"]

        return color_mapping


@click.command()
@click.argument("pin_file", required=True)
@click.option(
    "-p",
    "--pout",
    multiple=True,
    required=True,
    help=".pout MS²Rescore file, multiple .pout files possible with multiple flags",
)
@click.option(
    "-d",
    "--pout_dec",
    multiple=True,
    required=True,
    help=".pout_dec MS²Rescore file, multiple .pout_dec files possible with multiple flags",
)
@click.option(
    "-f",
    "--feature_sets",
    multiple=True,
    required=True,
    help="Features sets used for rescoring, if multiple pout files than multiple feature set names are required",
)
@click.option("-s", "--score_metric", required=True, help="Score metric used in the pin file")
@click.option("-o", "--output_filename", default="MS²Rescore_plots", help="output_name")
@click.option(
    "-w", "--weights_file", default=None, help="Percolator weight file to plot feature importances"
)
@click.option(
    "-n", "--sample_name", default="MS²Rescore run", help="Sample name used for generating plots"
)
@click.option("--fdr", default="0.01", help="Comma separated FDR values to plot PSMs")
def main(**kwargs):
    """
    Plot different analysis plots for the PIN_FILE, POUT_FILE and POUT_DEC_FILE from MS²Rescore
    """

    if not (len(kwargs["pout"]) == len(kwargs["pout_dec"])) & (
        len(kwargs["pout"]) == len(kwargs["feature_sets"])
    ):
        raise MS2RescoreError("Pout, pout_dec and feature_sets should be of equal length")

    kwargs["fdr"] = [float(fdr) for fdr in kwargs["fdr"].split(",")]
    logger.info(f"Create plots with these FDR vales: {kwargs['fdr']}")
    RescoreRecord.empty_rerecs()
    PIN(kwargs["pin_file"], kwargs["sample_name"], kwargs["score_metric"])
    for pout, pout_dec, feature_sets in zip(
        kwargs["pout"], kwargs["pout_dec"], kwargs["feature_sets"]
    ):
        POUT(pout, pout_dec, kwargs["sample_name"], feature_sets)
    if kwargs["weights_file"]:
        weights = PERCWEIGHT(
            kwargs["weights_file"],
            "MS²Rescore",
            kwargs["sample_name"],
        )
    logger.info(f"Saving plots to {kwargs['output_filename']}.pdf")
    RescoreRecord.save_plots_to_pdf(kwargs["output_filename"] + ".pdf", list(kwargs["fdr"]))


if __name__ == "__main__":
    main()
