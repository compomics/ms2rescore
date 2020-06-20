"""Percolator integration."""

import re
from io import StringIO
from typing import Dict, List, Optional

import pandas as pd

class PercolatorIn:
    """Percolator In (PIN)."""

    def __init__(
        self,
        path: Optional[str] = None,
        modification_mapping: Optional[Dict[float, str]] = None
    ):
        """
        Percolator In (PIN).

        Parameters
        ----------
        modification_mapping: dict, optional
            Mapping of mass shifts -> modification name, e.g.
            {-17.02655: "Gln->pyro-Glu"}. For matching, mass shifts are rounded to three
            decimals to avoid rounding issues.
        """
        # Parameters
        self.path = path
        self.modification_mapping = modification_mapping

        # Attributes
        self.modification_pattern = r"\[([0-9\-\.]*)\]"

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
            self._modification_mapping = {
                round(shift, 3): name for shift, name in value.items()
            }
        else:
            self._modification_mapping = value

    def _find_mods_recursively(
        self, mod_seq: str, mod_list: Optional[List[str]] = None
    ) -> List[str]:
        """Find modifications in modified sequence string recursively."""
        if not mod_list:
            mod_list = []
        match = re.search(self.modification_pattern, mod_seq)
        if match:
            mod_location = str(match.start())
            mod_shift = round(float(match.group(1)), 3)
            mod_name = self.modification_mapping[mod_shift]
            mod_list.extend([mod_location, mod_name])
            mod_seq = re.sub(self.modification_pattern, "", mod_seq, count=1)
            mod_list = self._find_mods_recursively(mod_seq, mod_list)
        return mod_list

    def _get_peprec_modifications(self, modified_sequence: str):
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

    def _get_unmodified_sequence(self, modified_sequence: str):
        """Parse modified sequence to get unmodified sequence."""
        unmod_seq = ".".join(modified_sequence.split(".")[1:-1])
        unmod_seq = re.sub(self.modification_pattern, "", unmod_seq, count=0)
        return unmod_seq

    def add_peprec_modifications_column(self):
        """Add PEPREC-style modifications column to PIN DataFrame."""
        self.df["modifications"] = self.df["Peptide"].apply(
            self._get_peprec_modifications
        )

    @staticmethod
    def fix_tabs(path: str, prot_sep: Optional[str] = '|||'):
        """
        Make tmp version of PIN file with fixed Proteins column separator.

        In a PIN file, multiple proteins in the Proteins column are separated by a tab,
        which is the same separator used to separate different columns in the file. This
        makes it impossible to read with pandas. This function makes a temporary copy of
        the PIN file with the Proteins tab-separations replaced with another string.

        """
        fixed_pin = StringIO()
        with open(path, 'rt') as pin_in:
            numcol = None
            for i, line in enumerate(pin_in):
                if i == 0 & line.startswith('SpecId'):
                    numcol = len(line.split('\t'))
                    fixed_pin.write(line)
                elif i == 1 & line.startswith('DefaultDirection'):
                    pass
                    #fixed_pin.write(line)
                else:
                    r = line.strip().split('\t')
                    r_cols = r[:numcol-1]
                    r_proteins = r[numcol-1:]
                    r_cols.append(prot_sep.join(r_proteins))
                    fixed_pin.write('\t'.join(r_cols) + '\n')

        fixed_pin.seek(0)
        return fixed_pin

    @staticmethod
    def write_with_tabs(file_object: StringIO, path: str, prot_sep: Optional[str] = "|||"):
        """Write PIN io.StringIO object to file with tab-separated Proteins column."""
        raise NotImplementedError

    def _get_default_direction_row(self):
        """Get default direction row from PIN file."""
        if not self.path:
            raise ValueError("PIN path is None. First set path.")
        default_direction = None
        with open(self.path, 'rt') as pin_in:
            for i, line in enumerate(pin_in):
                if i == 1 & line.startswith('DefaultDirection'):
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
        self.df = pd.read_csv(self.fix_tabs(self.path), sep='\t')

    def write(self, path: Optional[str] = None):
        """Write PIN to file."""
        raise NotImplementedError
