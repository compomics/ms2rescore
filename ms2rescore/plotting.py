"""Plot MS²ReScore results."""

from collections import defaultdict
from typing import List, Optional

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pyteomics.auxiliary import qvalues
from statsmodels.distributions.empirical_distribution import ECDF

from ms2rescore.percolator import PercolatorIn


class _REREC:
    rerecs = []

    def __init__(self) -> None:
        self.loss_gain_df = None
        self.unique_df = None

    @classmethod
    def empty_rerecs(cls):
        cls.rerecs = []

    @classmethod
    def show_rerec_items(cls):
        for rerec in cls.rerecs:
            print(rerec)

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
            list(x)
            for x in [
                self.df[self.df[decoy_label]][score_label],
                self.df[~self.df[decoy_label]][score_label],
            ]
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
            self.df.sort_values(score_label)[score_label],
            self.df.sort_values(score_label)[q_label],
        )
        if fdr_threshold:
            axes[1].vlines(
                x=score_cutoff, ymin=0, ymax=axes[1].get_ylim()[1], linestyles="dashed"
            )
        axes[1].set_ylabel("q-value")
        axes[1].set_xlabel(score_name)

        # PP plot
        ratio = (
            self.df[decoy_label].value_counts()[True]
            / self.df[decoy_label].value_counts()[False]
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

        ax.legend()
        fig.set_size_inches(15, 12)

        return ax

    @classmethod
    def _separate_unique_peptides(cls, FDR_threshold=[0.01]):
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

    @classmethod
    def unique_count_plot(cls):
        if not any(cls.unique_df):
            cls._separate_unique_peptides()
        cls.unique_df["count"] = cls.unique_df["upeps"].apply(len)

        g = sns.catplot(
            x="sample",
            y="count",
            data=cls.unique_df,
            hue="rescoring",
            kind="bar",
            col="FDR",
        )
        g.set_ylabels("number of unique peptides identified")
        g.set_xlabels("")

        return g

    @classmethod
    def calculate_loss_gain_df(
        cls, reference="Before rescoring", FDR_threshold=[0.01]
    ):
        if not any(cls.unique_df):
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
                total = len(ft_dict[reference])

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
        if not any(cls.loss_gain_df):
            cls.calculate_loss_gain_df(FDR_threshold=[FDR])
        if FDR not in cls.loss_gain_df["FDR"].unique():
            cls.calculate_loss_gain_df(FDR_threshold=[FDR])

        fig = plt.figure()
        tmp = cls.loss_gain_df[cls.loss_gain_df["FDR"] == FDR]

        for sample in zip(
            tmp["sample"].unique(), list(range(1, len(tmp["sample"].unique())))
        ):
            if sample[1] == len(tmp["sample"].unique()):
                ax = fig.add_subplot(
                    int(f"{tmp['sample'].unique()}1{sample[1]}"), frameon=False
                )
                ax.title.set_text(sample[0])
                sns.barplot(
                    y="feature",
                    x="gain",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#2FA92D",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                )
                sns.barplot(
                    y="feature",
                    x="shared",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#1AA3FF",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                )
                sns.barplot(
                    y="feature",
                    x="loss",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#FF0000",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                )
                ax.axes.set_xlabel("unique identified peptides (%)")
                ax.axes.set_ylabel("")
            else:
                ax = fig.add_subplot(
                    int(f"{tmp['sample'].unique()}1{sample[1]}"), frameon=False
                )
                ax.title.set_text(sample[0])
                sns.barplot(
                    y="feature",
                    x="gain",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#2FA92D",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                )
                sns.barplot(
                    y="feature",
                    x="shared",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#1AA3FF",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                )
                sns.barplot(
                    y="feature",
                    x="loss",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#FF0000",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                )
                ax.axes.set_xlabel("")
                ax.axes.set_ylabel("")
        fig.set_size_inches(15, 12)

        return fig


class PIN(_REREC):
    def __init__(self, path_to_pin, score_metric, sample_name) -> None:
        super().__init__()
        self.score_metric = score_metric
        self.type = "pin"
        self.rescore_features = "Before rescoring"
        self.name = sample_name
        self.df = self._read_pin_file(path_to_pin)
        self.rerecs.append(self)

    def __str__(self) -> str:
        return f"Sample name: {self.name}\nType: {self.type}\nRescore status: {self.rescore_features}"

    def _read_pin_file(self, path_to_pin):
        pin = PercolatorIn(path_to_pin)
        pin_qvalues = pd.DataFrame(
            qvalues(
                pin.df,
                key=pin.df[self.score_metric],
                is_decoy=pin.df["Label"] == -1,
                reverse=True,
                remove_decoy=False,
                formula=1,
                full_output=True
            )
        )
        return pin_qvalues[["SpecId", "is decoy", "score", "q", "Peptide"]].rename(columns={"SpecId" : "PSMId", "Peptide": "peptide"})


class POUT(_REREC):
    def __init__(
        self, path_to_target_pout, path_to_decoy_pout_, rescoring_features, sample_name
    ) -> None:
        super().__init__()
        self.target_pout = path_to_target_pout
        self.decoy_pout = path_to_decoy_pout_
        self.df = self._read_pout_file
        self.type = "pout"
        self.rescore_features = "After rescoring: " + rescoring_features
        self.name = sample_name
        self.df = self._read_pout_file()
        self.rerecs.append(self)

    def __str__(self) -> str:
        return f"Sample name: {self.name}\nType: {self.type}\nRescore status: {self.rescore_features}"

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

        pout_qvalues = pout[["PSMId","score", "q-value", "Label", "peptide"]].rename(
            columns={"q-value": "q", "Label": "is decoy"}
        )
        pout_qvalues["is decoy"] = pout["Label"] == -1

        return pout_qvalues
