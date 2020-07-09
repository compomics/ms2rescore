"""Percolator integration."""

import re
import logging
from io import StringIO
from typing import Dict, List, Optional

import pandas as pd

from ms2rescore.peptide_record import PeptideRecord


class PercolatorIn:
    """Percolator In (PIN)."""

    def __init__(
        self,
        path: Optional[str] = None,
        modification_mapping: Optional[Dict[float, str]] = None,
    ):
        """
        Percolator In (PIN).

        Parameters
        ----------
        path: str, optional
            Path to PIN file.
        modification_mapping: dict, optional
            Mapping of mass shifts -> modification name, e.g.
            {-17.02655: "Gln->pyro-Glu"}. Keys can either be strings or floats
            (representing mass shifts). If the keys are floats, they are rounded to
            three decimals to avoid rounding issues while matching. If None, the
            original modification labels from the PIN file will be used.

        """
        # Attributes
        self.modification_pattern = r"\[([^\[^\]]*)\]"

        # Parameters
        self.path = path
        self.modification_mapping = modification_mapping

        if self.path:
            self.read()

    @property
    def modification_mapping(self):
        """Get modification_mapping."""
        return self._modification_mapping

    @modification_mapping.setter
    def modification_mapping(self, value):
        """Set modification_mapping."""
        if value:
            if all([isinstance(label, str) for label in value.keys()]):
                self.modification_pattern = r"\[([^\[^\]]*)\]"
                self._modification_label_type = "str"
                self._modification_mapping = value
            elif all([isinstance(label, float) for label in value.keys()]):
                self.modification_pattern = r"\[([0-9\-\.]*)\]"
                self._modification_label_type = "float"
                self._modification_mapping = {
                    round(shift, 3): name for shift, name in value.items()
                }
            else:
                raise TypeError(
                    "modification_mapping keys cannot have mixed types and should "
                    "either be string or float."
                )
        else:
            self.modification_pattern = r"\[([^\[^\]]*)\]"
            self._modification_mapping = None
            self._modification_label_type = None

    def _find_mods_recursively(
        self, mod_seq: str, mod_list: Optional[List[str]] = None
    ) -> List[str]:
        """
        Find modifications in modified sequence string recursively.

        TODO: Fix handling of modifications on different residues with identical mass
        shifts, while also handing residue-aspecific modifications (e.g. N-terminal).
        """
        if not mod_list:
            mod_list = []
        match = re.search(self.modification_pattern, mod_seq)
        if match:
            mod_location = str(match.start())
            if self.modification_mapping:
                if self._modification_label_type == "float":
                    mod_label = round(float(match.group(1)), 3)
                elif self._modification_label_type == "str":
                    mod_label = match.group(1)
                else:
                    raise ValueError(
                        "Unsupported modification_label_type {}".format(
                            self._modification_label_type
                        )
                    )
                mod_name = self.modification_mapping[mod_label]
            else:
                mod_name = match.group(1)
            mod_list.extend([mod_location, mod_name])
            mod_seq = re.sub(self.modification_pattern, "", mod_seq, count=1)
            mod_list = self._find_mods_recursively(mod_seq, mod_list)
        return mod_list

    def _get_peprec_modifications(self, modified_sequence: str) -> str:
        """Parse modified sequence to get PEPREC-style modifications."""
        mod_seq = ".".join(modified_sequence.split(".")[1:-1])
        mod_list = self._find_mods_recursively(mod_seq)

        # In PIN, N-terminal mods are placed on position 1, instead of 0. This would
        # lead to issues in ms2pip if another modification is also present on the first
        # amino acid. Fix by setting N-term modification to 0.
        if len(mod_list) > 3:
            if mod_list[0] == "1" and mod_list[2] == "1":
                mod_list[0] = "0"

        mods = "|".join(mod_list)

        if not mods:
            mods = "-"
        return mods

    def _get_unmodified_sequence(self, modified_sequence: str) -> pd.Series:
        """Parse modified sequence to get unmodified sequence."""
        unmod_seq = ".".join(modified_sequence.split(".")[1:-1])
        unmod_seq = re.sub(self.modification_pattern, "", unmod_seq, count=0)
        return unmod_seq

    def _get_peprec_modifications_column(self) -> pd.Series:
        """Get PEPREC-style modifications column to PIN DataFrame."""
        return self.df["Peptide"].apply(self._get_peprec_modifications)

    def _get_sequence_column(self) -> pd.Series:
        """Get unmodified sequence column to PIN DataFrame."""
        return self.df["Peptide"].apply(self._get_unmodified_sequence)

    def _get_charge_column(self) -> pd.Series:
        """Get charge column from one-hot encoded `ChargeX` columns."""
        charge_cols = [col for col in self.df.columns if col.startswith("Charge")]
        assert (self.df[charge_cols] == 1).any(axis=1).all(), (
            "Not all PSMs have" " an assigned charge state."
        )
        return (
            self.df[charge_cols]
            .rename(
                columns={col: int(col.replace("Charge", "")) for col in charge_cols}
            )
            .idxmax(1)
        )

    def _get_spectrum_index_column(self, pattern: Optional[str] = None) -> pd.Series:
        """Get spectrum index column, parsed from Percolator SpecId."""
        if not pattern:
            pattern = r".+_([0-9]+)_[0-9]+_[0-9]+"
        id_col = self.df["SpecId"].str.extract(pattern, expand=False).astype(int)
        assert (
            ~id_col.duplicated().any()
        ), "Issue in matching spectrum IDs, duplicates found."
        assert ~id_col.isna().any(), "Issue in matching spectrum IDs, NaN found."
        return id_col

    def add_peprec_modifications_column(self):
        """Add PEPREC-style modifications column to PIN DataFrame."""
        self.df["modifications"] = self._get_peprec_modifications_column()

    def add_sequence_column(self):
        """Add unmodified sequence column to PIN DataFrame."""
        self.df["sequence"] = self._get_sequence_column()

    def add_charge_column(self):
        """Add charge column from one-hot encoded `ChargeX` columns."""
        self.df["charge"] = self._get_charge_column()

    def add_spectrum_index_column(
        self, pattern: Optional[str] = None, label: Optional[str] = "spectrum_index"
    ):
        """
        Add spectrum index column, parsed from Percolator SpecId.

        See https://github.com/percolator/percolator/issues/147 for pattern explanation.
        """
        self.df[label] = self._get_spectrum_index_column(pattern=pattern)

    def get_spectrum_filename(
        self, pattern: Optional[str] = r"(.+)_[0-9]+_[0-9]+_[0-9]+"
    ) -> str:
        """
        Get spectrum filename, parsed from Percolator SpecId.

        See https://github.com/percolator/percolator/issues/147 for pattern explanation.
        """
        spectrum_filenames = (
            self.df["SpecId"].str.extract(pattern, expand=False).unique()
        )
        if len(spectrum_filenames) == 1:
            return spectrum_filenames[0]
        else:
            raise ValueError("Multiple spectrum filenames found in single PIN file.")

    @staticmethod
    def fix_tabs(path: str, id_column: str = "SpecId", prot_sep: Optional[str] = "|||"):
        """
        Return StringIO instance of PIN/POUT file with fixed Proteins column separator.

        In a PIN/POUT file, multiple proteins in the Proteins column are separated by a
        tab, which is the same separator used to separate different columns in the file.
        This makes it impossible to read with pandas. This function makes a temporary
        copy of the PIN file with the Proteins tab-separations replaced with another
        string.

        Parameters
        ----------
        path: str
            Path to input file
        id_column: str
            Label of the ID column: `SpecId` for PIN files, `PSMId` for POUT files.
        prot_sep: str
            Separator to use in proteins column.

        """
        fixed_file = StringIO()
        with open(path, "rt") as file_in:
            numcol = None
            for i, line in enumerate(file_in):
                if i == 0 & line.startswith(id_column):
                    numcol = len(line.split("\t"))
                    fixed_file.write(line)
                elif i == 1 & line.startswith("DefaultDirection"):
                    pass
                    # fixed_file.write(line)
                else:
                    r = line.strip().split("\t")
                    r_cols = r[: numcol - 1]
                    r_proteins = r[numcol - 1 :]
                    r_cols.append(prot_sep.join(r_proteins))
                    fixed_file.write("\t".join(r_cols) + "\n")

        fixed_file.seek(0)
        return fixed_file

    @staticmethod
    def write_with_tabs(
        file_object: StringIO, path: str, prot_sep: Optional[str] = "|||"
    ):
        """Write PIN io.StringIO object to file with tab-separated Proteins column."""
        raise NotImplementedError

    def _get_default_direction_row(self):
        """Get default direction row from PIN file."""
        if not self.path:
            raise ValueError("PIN path is None. First set path.")
        default_direction = None
        with open(self.path, "rt") as pin_in:
            for i, line in enumerate(pin_in):
                if i == 1 & line.startswith("DefaultDirection"):
                    default_direction = line.strip().split("\t")
                if i > 1:
                    break
        return default_direction

    def read(self, path: Optional[str] = None):
        """Read PIN file into pandas.DataFrame."""
        if path:
            self.path = path
        if not self.path:
            raise ValueError("No path for PIN file defined.")
        self.df = pd.read_csv(self.fix_tabs(self.path), sep="\t")

    def write(self, path: Optional[str] = None):
        """Write PIN to file."""
        raise NotImplementedError

    def to_peptide_record(
        self,
        extract_spectrum_index: Optional[bool] = True,
        spectrum_index_pattern: Optional[str] = None,
        score_column_label: Optional[str] = None,
    ) -> PeptideRecord:
        """
        Convert PercolatorIn to PeptideRecord.

        Parameters
        ----------
        extract_spectrum_index: bool (optional)
            If True, extract spectrum index from spec_id `SpecId` with the
            `spectrum_index_pattern` regex pattern, or the default pattern:
            `.+_([0-9]+)_[0-9]+_[0-9]+`. If False, copy the original values in `SpecId`.
        spectrum_index_pattern: str (optional)
            Regex pattern with one group matching the MGF spectrum index in the PIN
            `SpecId` column.
        score_column_label: str (optional)
            Label of column to copy to PeptideRecord as `psm_score`. If None, common
            column labels will be tested and used if available.

        Returns
        -------
        ms2rescore.peptide_record.PeptideRecord

        """
        # Assign one of the default score column labels, if available
        if not score_column_label:
            for col in ["lnEValue", "Score", "hyperscore"]:
                if col in self.df.columns:
                    score_column_label = col
                    logging.debug(
                        "Found score column in PIN with label %s", score_column_label
                    )
            else:
                logging.debug("No known score column found.")

        # Create new peprec
        peprec_df = pd.DataFrame()
        if extract_spectrum_index:
            peprec_df["spec_id"] = self._get_spectrum_index_column(
                spectrum_index_pattern
            )
        else:
            peprec_df["spec_id"] = self.df["SpecId"]
        peprec_df["peptide"] = self._get_sequence_column()
        peprec_df["modifications"] = self._get_peprec_modifications_column()
        peprec_df["charge"] = self._get_charge_column()
        if score_column_label:
            peprec_df["psm_score"] = self.df[score_column_label]
        peprec_df["is_decoy"] = self.df["Label"] == -1

        return PeptideRecord.from_dataframe(peprec_df)
