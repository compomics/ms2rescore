"""
Peptide record (PEPREC).

TODO: Move module to ms2pip
"""

import pandas as pd


class InvalidPeprecError(Exception):
    """Invalid PEPREC file."""

    pass


class PeptideRecord:
    """Peptide record (PEPREC)."""

    def __init__(self, path=None, extra_required_columns=None):
        """Peptide record (PEPREC)."""
        self.path = path
        self.extra_required_columns = extra_required_columns
        self.df = None

        if path:
            self.from_csv(path)

    def __repr__(self):
        """Get string representation of PeptideRecord."""
        return self.df.__repr__()

    def _validate_header(self, path: str):
        """Validate PEPREC header."""
        with open(path, "rt") as f:
            line = f.readline()
            if line[:7] != "spec_id":
                raise InvalidPeprecError("PEPREC header should start with `spec_id`")

    def _infer_separater(self, path: str) -> str:
        """Infer separator in PEPREC file."""
        with open(path, "rt") as f:
            line = f.readline()
            separater = line[7]
        return separater

    def _validate_columns(self, path: str):
        """Validate header of PEPREC file."""
        required_columns_default = [
            "spec_id",
            "peptide",
            "modifications",
            "charge",
        ]

        required_columns_rt = [
            "psm_score",
            "observed_retention_time",
        ]

        all_columns = required_columns_default + required_columns_rt
        separator = self._infer_separater(path)

        with open(path, "rt") as f:
            line = f.readline()
            line = line.split(separator)
            for col in all_columns:
                if col not in line:
                    raise InvalidPeprecError("Required column missing from header", col)

    def from_csv(self, path, **kwargs):
        """Read PEPREC from CSV."""
        self._validate_header(path)
        self._validate_columns(path)
        sep = self._infer_separater(path)
        self.df = pd.read_csv(path, sep=sep, index_col=None, **kwargs)

        self.df["modifications"] = self.df["modifications"].fillna("-")

    def to_csv(self, path=None, **kwargs):
        """Save PEPREC to CSV, if Path is None, overwrite existing PEPREC file."""
        if not path and self.path:
            path = self.path
        elif not path and not self.path:
            raise ValueError("No known path to write PEPREC file.")
        self.df.to_csv(path, sep=" ", index=False, float_format="%.3f", **kwargs)
