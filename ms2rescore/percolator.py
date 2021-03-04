"""Percolator integration."""


import logging
import os
import re
import subprocess
from io import StringIO
from math import floor
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd

from ms2rescore._exceptions import MS2ReScoreError
from ms2rescore.peptide_record import PeptideRecord


logger = logging.getLogger(__name__)


class PercolatorInError(MS2ReScoreError):
    """Error while processing Percolator In file."""


class UnknownModificationLabelStyleError(PercolatorInError):
    """Could not infer modification label style."""

    pass


class UnknownModificationError(PercolatorInError):
    """Modification not found in `modification_mapping`."""

    pass


class PercolatorIn:
    """Percolator In (PIN)."""

    def __init__(
        self,
        path: Optional[str] = None,
        modification_mapping: Optional[
            Dict[Tuple[Union[str, None], Union[float, str]], str]
        ] = None,
        invalid_amino_acids: Optional[str] = r"[BJOUXZ]"
    ):
        """
        Percolator In (PIN).

        Parameters
        ----------
        path: str, optional
            Path to PIN file.
        modification_mapping: dict, optional
            Mapping of modification residue and PIN label to the modification name.
            Modification residue is the one-letter amino acid code, or None.
            Modification labels can either be strings (e.g.
            `{("Q", "UNIMOD:28"): "Gln->pyro-Glu"}`) or floats representing mass shifts
            (e.g. `{("Q", -17.02655): "Gln->pyro-Glu"}). If the keys are floats, they
            are rounded up to three decimals to avoid rounding issues while matching.
            If None, the original modification labels from the PIN file will be used.
        invalid_amino_acids: str, optional
            regex pattern of invalid amino acids. PSMs containing these amino acids will
            be dropped. (default: `r"[BJOUXZ]"`)

        """
        # Attributes
        self.modification_pattern = r"\[([^\[^\]]*)\]"
        self.invalid_amino_acids = invalid_amino_acids

        # Parameters
        self.path = path
        self.modification_mapping = modification_mapping

        if self.path:
            self.read()

    @staticmethod
    def _round(n, decimals=0):
        multiplier = 10 ** decimals
        return floor(n * multiplier + 0.5) / multiplier

    @property
    def modification_mapping(self):
        """Get modification_mapping."""
        return self._modification_mapping

    @modification_mapping.setter
    def modification_mapping(
        self, value: Optional[Dict[Tuple[Union[str, None], Union[float, str]], str]]
    ):
        """Set modification_mapping."""
        if value:
            mod_labels = [key[1] for key in value.keys()]
            if all([isinstance(label, str) for label in mod_labels]):
                self.modification_pattern = r"\[([^\[^\]]*)\]"
                self._modification_label_type = "str"
                self._modification_mapping = value
            elif all([isinstance(label, float) for label in mod_labels]):
                self.modification_pattern = r"\[([0-9\-\.]*)\]"
                self._modification_label_type = "float"
                self._modification_mapping = {
                    (aa, self._round(shift)): name
                    for (aa, shift), name in value.items()
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

    def modification_mapping_from_list(
        self, modification_config, label_style="infer"
    ):
        """
        Generate modification_mapping from ms2rescore configuration.

        Parameters
        ----------
        config: Dict
            ms2rescore configuration
        label_style: str
            Modification label style, any of {"infer", "unimod_accession", "mass_shift"}
            If `infer`, the style will be attempted to be infered from the PIN file.

        """
        if label_style == "infer":
            label_style = self._infer_modification_label_style()
            self.modification_mapping_from_list(
                modification_config, label_style=label_style
            )
        elif label_style == "unimod_accession":
            self.modification_mapping = {
                (mod["amino_acid"], f"UNIMOD:{mod['unimod_accession']}"): mod["name"]
                for mod in modification_config
            }
        elif label_style == "mass_shift":
            self.modification_mapping = {
                (mod["amino_acid"], mod["mass_shift"]): mod["name"]
                for mod in modification_config
            }
        else:
            raise ValueError(label_style)

    def _infer_modification_label_style(self) -> str:
        """Infer modification label style (e.g. unimod_accession or mass shift)."""
        # TODO: What happens if there are no mods? Crashes?
        modified_peptides = self.df[
            self.df["Peptide"].str.contains(r"\[([^\[^\]]*)\]", regex=True)
        ]["Peptide"]
        if (
            modified_peptides.str.extract(r"\[([^\[^\]]*)\]", expand=False)
            .str.startswith("UNIMOD:")
            .all()
        ):
            mod_notation = "unimod_accession"
        elif (
            pd.to_numeric(
                modified_peptides.str.extract(r"\[([^\[^\]]*)\]", expand=False),
                errors="coerce",
            )
            .notnull()
            .all()
        ):
            mod_notation = "mass_shift"
        else:
            raise UnknownModificationLabelStyleError()
        return mod_notation

    def _find_mods_recursively(
        self, mod_seq: str, mod_list: Optional[List[str]] = None
    ) -> List[str]:
        """
        Find modifications in modified sequence string recursively.
        """
        if not mod_list:
            mod_list = []
        match = re.search(self.modification_pattern, mod_seq)
        if match:
            mod_location = match.start()
            mod_aa = mod_seq[mod_location - 1] if mod_location > 0 else mod_seq[match.end()]
            if self.modification_mapping:
                if self._modification_label_type == "float":
                    mod_label = self._round(float(match.group(1)))
                elif self._modification_label_type == "str":
                    mod_label = match.group(1)
                else:
                    raise ValueError(
                        "Unsupported modification_label_type {}".format(
                            self._modification_label_type
                        )
                    )
                # First try with specific amino acid, then try with None (e.g. N-term)
                if (mod_aa, mod_label) in self.modification_mapping:
                    mod_name = self.modification_mapping[mod_aa, mod_label]
                elif (None, mod_label) in self.modification_mapping:
                    mod_name = self.modification_mapping[None, mod_label]
                else:
                    raise UnknownModificationError(
                        f"No modification with label `{mod_label}` for amino acid "
                        f"`{mod_aa}` found in modification_mapping "
                        f"`{self.modification_mapping}`."
                    )
                # Set C-terminal modifications location to "-1"
                if match.end() == len(mod_seq):
                    mod_location = -1
            else:
                mod_name = match.group(1)
            mod_list.extend([str(mod_location), mod_name])
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
        if not (self.df[charge_cols] == 1).any(axis=1).all():
            raise PercolatorInError("Not all PSMs have an assigned charge state.")
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
        if id_col.duplicated().any():
            raise PercolatorInError("Issue in matching spectrum IDs, duplicates found.")
        if id_col.isna().any():
            raise PercolatorInError("Issue in matching spectrum IDs, NaN found.")
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

    def drop_invalid_amino_acids(self):
        """Drop all PSMs (rows) with peptides containing invalid amino acids."""
        if "sequence" in self.df.columns:
            sequences = self.df['sequence']
        else:
            sequences = self._get_sequence_column()
        to_drop = sequences[
            sequences.str.contains(self.invalid_amino_acids, regex=True)
        ].index
        if len(to_drop) > 0:
            logger.warning(
                "Dropping %i PSMs from PIN due to invalid amino acids (%s)",
                len(to_drop),
                self.invalid_amino_acids
            )
            self.df = self.df.drop(index=to_drop)

    @staticmethod
    def fix_tabs(
        path: str, id_column: str = "SpecId", prot_sep: Optional[str] = "|||"
    ) -> StringIO:
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

        Returns
        -------
        io.StringIO
            io.StringIO object with fixed PIN/POUT file.

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
        if self.invalid_amino_acids:
            self.drop_invalid_amino_acids()

    def write(self, path: Optional[str] = None):
        """Write PIN to file."""
        raise NotImplementedError

    def get_feature_table(self) -> pd.DataFrame:
        """Return pandas.DataFrame with PIN features."""
        pin_cols = ["SpecId", "Label", "ScanNr", "Peptide", "Proteins"]
        peprec_cols = [
            "spec_id",
            "peptide",
            "modifications",
            "charge",
            "observed_retention_time"
        ]
        misc_cols = ["tandem_id", "ModPeptide"]
        non_feature_cols = pin_cols + peprec_cols + misc_cols
        feature_cols = [col for col in self.df.columns if col not in non_feature_cols]
        return self.df[feature_cols]

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
            for col in ["lnEValue", "hyperscore", "Score"]:
                if col in self.df.columns:
                    score_column_label = col
                    logger.debug(
                        "Found score column in PIN with label %s", score_column_label
                    )
                    break
            if not score_column_label:
                logger.debug("No known score column found.")

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
        peprec_df["label"] = self.df["Label"]
        peprec_df["Proteins"] = self.df["Proteins"]

        peprec = PeptideRecord.from_dataframe(peprec_df)

        return peprec


def run_percolator_converter(
    converter_name: str,
    path_to_id_file: str,
    path_to_pin_file: str,
    id_decoy_pattern: Optional[str] = None,
    log_level: str = "info",
):
    """
    Run Percolator PIN Converter.

    Parameters
    ----------
    converter_name: str
        name of the converter to use; currently supported are {"msgf2pin", "tandem2pin"}
    path_to_id_file: str
        path to identification file (input to converter)
    path_to_pin_file: str'
        path to pin file (output of converter)
    id_decoy_pattern: str, optional
        identifications file decoy pattern (e.g. REVERSED), passed to `-P` option of
        Percolator converter; default depends on input filetype
    log_level: str, optional
        log_level, default is `info`

    """
    default_decoy_patterns = {"msgf2pin": "XXX", "tandem2pin": "REVERSED"}
    if converter_name not in default_decoy_patterns.keys():
        raise NotImplementedError(converter_name)
    if not id_decoy_pattern:
        id_decoy_pattern = default_decoy_patterns[converter_name]

    command = [
        converter_name,
        "-P",
        id_decoy_pattern,
        "-o",
        os.path.abspath(path_to_pin_file),
        os.path.abspath(path_to_id_file),
    ]

    logger.info("Running Percolator PIN converter")
    logger.debug(' '.join(command))

    subprocess.run(command, capture_output=log_level == "debug", check=True)
