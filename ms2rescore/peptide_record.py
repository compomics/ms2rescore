"""
Peptide record (PEPREC).

TODO: Move module to ms2pip
"""

from typing import List, Optional

import pandas as pd


class InvalidPeprecError(Exception):
    """Invalid PEPREC file."""

    pass


class PeptideRecord:
    """Peptide record (PEPREC)."""

    def __init__(
        self,
        path: Optional[str] = None,
        context: str = "default",
        extra_required_columns: Optional[List[str]] = None,
    ):
        """
        Peptide record (PEPREC).

        Parameters
        ----------
        path: str, None
            Path to PEPREC file. If not None, file is read on instance creation.
        context: str
            Context of PEPREC. Is used for determining the required columns. Can be any
            of the following: {default, ms2rescore}.
        extra_required_columns: List[str]
            Extra required columns to be validated.

        Attributes
        ----------
        df: pandas.DataFrame
            DataFrame containing peptide record content.

        Methods
        -------
        from_csv(path, **kwargs)
            Read PEPREC from CSV.
        to_csv(path, **kwargs)
            Save PEPREC to CSV, if Path is None, overwrite existing PEPREC file.

        """
        self.path = path
        self.context = context
        self.extra_required_columns = extra_required_columns
        self.df = None

        if self.path:
            self.read_csv(self.path)

    def __repr__(self):
        """Get string representation of PeptideRecord."""
        return self.df.__repr__()

    @property
    def df(self) -> Optional[pd.DataFrame]:
        """Get DataFrame with PeptideRecord."""
        return self._df

    @df.setter
    def df(self, value: Optional[pd.DataFrame]):
        """Set DataFrame with PeptideRecord."""
        if isinstance(value, pd.DataFrame):
            self._validate_column_names(value)
        self._df = value

    def read_csv(self, path: str, **kwargs):
        """Read PEPREC from CSV."""
        self._validate_header(path)
        sep = self._infer_separator(path)
        self.df = pd.read_csv(path, sep=sep, index_col=None, **kwargs)
        self.df["modifications"] = self.df["modifications"].fillna("-")

    def to_csv(self, path: Optional[str] = None, **kwargs):
        """Save PEPREC to CSV, if Path is None, overwrite existing PEPREC file."""
        if not path and self.path:
            path = self.path
        elif not path and not self.path:
            raise ValueError("No known path to write PEPREC file.")
        self.df.to_csv(path, sep=" ", index=False, float_format="%.3f", **kwargs)

    def reorder_columns(self):
        """Reorder columns with canonical PeptideRecord columns first."""
        original_columns = self.df.columns.to_list()
        for col in self._required_columns:
            if col in original_columns:
                original_columns.remove(col)
            else:
                self._required_columns.remove(col)
        reordered_columns = self._required_columns + original_columns
        self.df = self.df[reordered_columns]

    @property
    def _required_columns(self):
        """Get required columns depending on PeptideRecord context."""
        required_columns = [
            "spec_id",
            "peptide",
            "modifications",
            "charge",
        ]
        if self.context == "ms2rescore":
            required_columns.extend(
                ["psm_score", "observed_retention_time",]
            )
        elif self.context == "retention_time":
            required_columns.remove("modifications")
        return required_columns

    def _validate_header(self, path: str):
        """Validate PEPREC header."""
        with open(path, "rt") as f:
            line = f.readline()
            if line[:7] != "spec_id":
                raise InvalidPeprecError("PEPREC header should start with `spec_id`.")

    def _validate_column_names(self, df: pd.DataFrame):
        """Validate header of PEPREC file."""
        for col in self._required_columns:
            if col not in df.columns:
                raise InvalidPeprecError(
                    f"Required column `{col}` missing from header."
                )

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame):
        """Create PeptideRecord instance with data from a pandas DataFrame."""
        peprec = cls()
        peprec.df = dataframe
        peprec.reorder_columns()
        return peprec

    @classmethod
    def from_csv(cls, *args, **kwargs):
        """Create PeptideRecord instance with data from CSV file."""
        peprec = cls()
        peprec.read_csv(*args, **kwargs)
        return peprec

    @staticmethod
    def _infer_separator(path: str) -> str:
        """Infer separator in PEPREC file."""
        with open(path, "rt") as f:
            line = f.readline()
            separator = line[7]
        return separator
