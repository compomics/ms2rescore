"""Plot MS²ReScore results."""

from abc import ABC
from collections import defaultdict
from typing import List, Optional

import matplotlib.axes
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np
import pandas as pd
import seaborn as sns
from pyteomics.auxiliary import qvalues
from statsmodels.distributions.empirical_distribution import ECDF


from ms2rescore.percolator import PercolatorIn
from ms2rescore._exceptions import MS2RescoreError


class RescoreRecord(ABC):
    rerecs = []
    loss_gain_df = pd.DataFrame()
    unique_df = pd.DataFrame()
    count_df = pd.DataFrame()

    def __init__(self) -> None:
        pass

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
        fig.set_size_inches(12, 10)

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
            y_label = "number of identified peptides"
            count_df = cls.count_df

        g = sns.catplot(
            x="sample",
            y="count",
            data=count_df,
            hue="rescoring",
            kind="bar",
            col="FDR",
            legend=False,
        )

        g.set_ylabels(y_label)
        g.set_xlabels("")
        g.add_legend(loc=7)
        g.fig.set_size_inches(12, 10)

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

        for sample in zip(
            tmp["sample"].unique(), list(range(1, len(tmp["sample"].unique()) + 1))
        ):
            if sample[1] == len(tmp["sample"].unique()):
                ax = fig.add_subplot(
                    int(f"{len(tmp['sample'].unique())}1{sample[1]}"), frameon=False
                )
                ax.title.set_text(f"{sample[1]}\nFDR={FDR}")
                sns.barplot(
                    y="feature",
                    x="gain",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#2FA92D",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                    ax=ax,
                )
                sns.barplot(
                    y="feature",
                    x="shared",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#1AA3FF",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                    ax=ax,
                )
                sns.barplot(
                    y="feature",
                    x="loss",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#FF0000",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                    ax=ax,
                )
                ax.axes.set_xlabel("unique identified peptides (%)")
                ax.axes.set_ylabel("")
            else:
                ax = fig.add_subplot(
                    int(f"{len(tmp['sample'].unique())}1{sample[1]}"), frameon=False
                )
                ax.title.set_text(sample[0])
                sns.barplot(
                    y="feature",
                    x="gain",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#2FA92D",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                    ax=ax,
                )
                sns.barplot(
                    y="feature",
                    x="shared",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#1AA3FF",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                    ax=ax,
                )
                sns.barplot(
                    y="feature",
                    x="loss",
                    data=tmp[tmp["sample"] == sample[0]],
                    color="#FF0000",
                    order=tmp[tmp["sample"] == sample[0]].sort_values("gain").feature,
                    ax=ax,
                )
                ax.axes.set_xlabel("")
                ax.axes.set_ylabel("")
        plt.grid(axis="x")
        ax.set_axisbelow(True)
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

        cls.qvalue_comparison()
        plt.tight_layout()
        pdf.savefig()

        cls.loss_gain_plot(FDR=0.01)
        plt.tight_layout()
        pdf.savefig()
        pdf.close()


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
        return f"Sample name: {self.name}\nType: {self.type}\nRescore status: {self.rescore_features}"

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

        pout_qvalues = pout[["PSMId", "score", "q-value", "Label", "peptide"]].rename(
            columns={"q-value": "q", "Label": "is decoy"}
        )
        pout_qvalues["is decoy"] = pout["Label"] == -1

        return pout_qvalues
