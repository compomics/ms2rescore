"""Extract PeptideRecord and search engine features from identification file."""

import re
import os
import logging
from abc import ABC, abstractmethod
from typing import Optional, Union, Dict, Tuple, List

import numpy as np
import pandas as pd
from pyteomics import tandem

from ms2rescore.percolator import PercolatorIn, run_percolator_converter
from ms2rescore.peptide_record import PeptideRecord
from ms2rescore.parse_mgf import parse_mgf
from ms2rescore.maxquant import MSMS
from ms2rescore.peptideshaker import ExtendedPsmReport


def parse_mgf_title_rt(
    path_to_mgf: Union[str, os.PathLike]
) -> Tuple[Dict[int, str], Dict[int, float]]:
    """Parse MGF file to extract title and retention time fields, by spectrum index."""
    titles = dict()
    retention_times = dict()
    with open(path_to_mgf, "rt") as mgf_in:
        index = 0
        for line in mgf_in:
            if line[0] == "B":
                if line.strip() == "BEGIN IONS":
                    index += 1
            if line[0] == "T":
                if line.startswith("TITLE="):
                    titles[index] = line[6:].strip()
            if line[0] == "R":
                if line.startswith("RTINSECONDS="):
                    retention_times[index] = float(line[12:].strip())
    return titles, retention_times


class _Pipeline(ABC):
    """ABC for pipeline from ID file to PeptideRecord and search engine features."""

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        # Attributes from config and output_basename
        self.path_to_id_file = config["general"]["identification_file"]
        self.passed_mgf_path = config["general"]["mgf_path"]
        self.output_basename = output_basename
        self.modification_list = config["ms2pip"]["modifications"]
        self.id_decoy_pattern = config["general"]["id_decoy_pattern"]
        self.log_level = config["general"]["log_level"]

        # General private attributes
        self._pin_spec_id_patterns = {
            "generic": r".+_([0-9]+)_[0-9]+_[0-9]+",
            "tandem": r".+_([0-9]+)_[0-9]+_[0-9]+",
            "msgfplus": r".+_SII_([0-9]+)_[0-9]+_[0-9]+_[0-9]+",
        }

        # Private attributes specific to pipeline, override these in each subclass
        self._path_to_original_pin = None
        self._pin_converter_name = None
        self._pin_spec_id_style = "generic"
        self._pin_modification_style = "infer"

        # Private attributes used by methods
        self._original_pin = None

    @staticmethod
    def _validate_mgf_path(
        passed_path: Union[None, str, os.PathLike],
        default_dir: Union[str, os.PathLike],
        expected_rootname: str,
        expected_ext: str = ".mgf"
    ) -> Union[str, os.PathLike]:
        """
        Infer MGF path from passed path and expected filename (e.g. from ID file).

        Parameters
        ----------
        passed_path : None, string, os.PathLike
            user-defined path to MGF file or directory containing MGF file
        default_dir : str, os.PathLike
            default directory for MGF file, used in combination with `expected_rootname`
            and `expected_ext` if `passed_path` is None
        expected_rootname : str, os.PathLike
            rootname of MGF file, as expected from, e.g., identification file
        expected_ext : str, optional
            expected filename extension, including period, default: ".mgf"

        """
        # Make sure that expected_rootname is in fact rootname without expected ext
        expected_rootname = os.path.basename(expected_rootname)
        if expected_rootname.endswith(expected_ext):
            expected_rootname = re.sub(".mgf$", "", expected_rootname)

        # If no mgf path configured, return expected rootname + ext in default directory
        if not passed_path:
            path_to_mgf_file = os.path.join(
                default_dir, expected_rootname + expected_ext
            )

        # If passed path is directory, join with expected rootname + ext
        elif os.path.isdir(passed_path):
            path_to_mgf_file = os.path.join(
                passed_path,
                expected_rootname + expected_ext
            )

        # If passed path is file, use that, but warn if basename doesn't match expected
        elif os.path.isfile(passed_path):
            passed_rootname = os.path.splitext(os.path.basename(passed_path))[0]
            if passed_rootname != expected_rootname:
                logging.warning(
                    "Passed MGF name root `%s` does not match MGF name root `%s` from "
                    "identifications file. Continuing with passed MGF name.",
                    passed_rootname,
                    expected_rootname
                )
            path_to_mgf_file = passed_path

        else:
            raise ValueError(
                "Configured `mgf_path` must be None or a path to an existing file or "
                "directory."
            )

        return path_to_mgf_file

    def _run_percolator_converter(self):
        """Run Percolator converter with settings specific to pipeline."""
        run_percolator_converter(
            self._pin_converter_name,
            self.path_to_id_file,
            self._path_to_original_pin,
            id_decoy_pattern=self.id_decoy_pattern,
            log_level=self.log_level
        )

    @property
    def path_to_mgf_file(self) -> Union[str, os.PathLike]:
        """Get path to MGF file, inferred from mgf_path (which can also be to a dir)."""
        path_to_mgf_file = self._validate_mgf_path(
            self.passed_mgf_path,
            os.path.dirname(self.path_to_id_file),
            os.path.basename(os.path.splitext(self.path_to_id_file)[0])
        )
        return path_to_mgf_file

    @property
    def original_pin(self) -> PercolatorIn:
        """Get PercolatorIn object from identification file."""
        if not self._original_pin:
            if not self._path_to_original_pin:
                raise TypeError("To load PIN, path_to_pin cannot be None.")
            elif not os.path.isfile(self._path_to_original_pin):
                self._run_percolator_converter()
            self._original_pin = PercolatorIn(self._path_to_original_pin)
            self._original_pin.modification_mapping_from_list(
                self.modification_list,
                label_style=self._pin_modification_style
            )
        return self._original_pin

    def peprec_from_pin(self) -> PeptideRecord:
        """Get PeptideRecord from PIN file and MGF file."""
        # Get peprec
        peprec = self.original_pin.to_peptide_record(
            spectrum_index_pattern=self._pin_spec_id_patterns[self._pin_spec_id_style]
        )

        # Map MGF titles and observed retention times
        titles, retention_times = parse_mgf_title_rt(self.path_to_mgf_file)
        peprec.df["observed_retention_time"] = peprec.df["spec_id"].map(retention_times)
        peprec.df["spec_id"] = peprec.df["spec_id"].map(titles)
        assert (
            ~peprec.df["observed_retention_time"].isna().any()
        ), "Could not map all MGF retention times to spectrum indices."
        assert (
            ~peprec.df["spec_id"].isna().any()
        ), "Could not map all MGF titles to spectrum indices."

        return peprec

    @abstractmethod
    def get_peprec(self) -> PeptideRecord:
        """Get PeptideRecord."""
        raise NotImplementedError

    @abstractmethod
    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        raise NotImplementedError


class PinPipeline(_Pipeline):
    """Percolator IN file to PeptideRecord and search engine features."""

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        super().__init__(config, output_basename)

        # Private attributes specific to pipeline
        self._path_to_original_pin = self.path_to_id_file
        self._pin_converter_name = None
        self._pin_spec_id_style = "generic"
        self._pin_modification_style = "infer"

    def get_peprec(self) -> PeptideRecord:
        """Get PeptideRecord from original PIN file."""
        return self.peprec_from_pin()

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        return self.original_pin.get_feature_table()


class MSGFPipeline(_Pipeline):
    """MSGFPlus mzid to PeptideRecord and search engine features."""

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        super().__init__(config, output_basename)

        # Private attributes specific to pipeline
        self._path_to_original_pin = output_basename + "_original.pin"
        self._pin_converter_name = "msgf2pin"
        self._pin_spec_id_style = "msgfplus"
        self._pin_modification_style = "infer"

    def get_peprec(self) -> PeptideRecord:
        """Get PeptideRecord from original PIN file."""
        return self.peprec_from_pin()

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        return self.original_pin.get_feature_table()


class TandemPipeline(_Pipeline):
    """X!Tandem XML to PeptideRecord and search engine features."""

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        super().__init__(config, output_basename)

        # Private attributes specific to pipeline
        self._path_to_original_pin = output_basename + "_original.pin"
        self._pin_converter_name = "tandem2pin"
        self._pin_spec_id_style = "tandem"
        self._pin_modification_style = "mass_shift"

    # Override method to use spectrum filenames found in PIN spec_ids as expected_root
    @property
    def path_to_mgf_file(self) -> Union[str, os.PathLike]:
        """Get path to MGF file, inferred from mgf_path (which can also be to a dir)."""
        path_to_mgf_file = self._validate_mgf_path(
            self.passed_mgf_path,
            os.path.dirname(self.path_to_id_file),
            self.original_pin.get_spectrum_filename()
        )
        return path_to_mgf_file

    def get_peprec(self) -> PeptideRecord:
        """Convert X!Tandem XML file and PIN to PEPREC."""
        # Load tandem dataframe with Pyteomics
        logging.debug("Converting X!Tandem XML to PEPREC...")
        tandem_df = tandem.DataFrame(self.path_to_id_file)
        tandem_df["id"] = tandem_df["id"].astype(int)
        if "RTINSECONDS" in tandem_df["scan"].loc[0]:
            # Sometimes "RTINSECONDS" is in scan number column...
            tandem_df["scan"] = tandem_df["scan"].str.replace(" RTINSECONDS.*", "")

        tandem_peprec_mapping = {
            "scan": "spec_id",
            "seq": "peptide",
            "z": "charge",
            "rt": "observed_retention_time",
            "expect": "psm_score",
            "hyperscore": "hyperscore_tandem",
            "id": "tandem_id"
        }

        peprec_df = tandem_df[tandem_peprec_mapping.keys()].rename(
            columns=tandem_peprec_mapping
        )
        # Set PSM score as -log(e-value)
        peprec_df["psm_score"] = - np.log(peprec_df["psm_score"])

        logging.debug("Adding modifications from original PIN to PEPREC...")
        pin = self.original_pin
        pin.add_peprec_modifications_column()
        pin.add_spectrum_index_column(label="tandem_id")
        peprec_df = peprec_df.merge(
            pin.df["modifications", "tandem_id", "hyperscore"],
            on="tandem_id"
        )
        # Validate merge by comparing the hyperscore columns
        assert (peprec_df["hyperscore_tandem"] == peprec_df["hyperscore"]).all()
        peprec_df.drop(
            columns=["tandem_id", "hyperscore_tandem"],
            axis="columns",
            inplace=True
        )

        peprec = PeptideRecord.from_dataframe(peprec_df)
        peprec.reorder_columns()
        return peprec

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        return self.original_pin.get_feature_table()


class MaxQuantPipeline(_Pipeline):
    """MaxQuant msms.txt to PeptideRecord and search engine features."""

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        super().__init__(config, output_basename)

        # Private attributes, specific to this pipeline
        mq_conf = config["maxquant_to_rescore"]
        self._modification_mapping = mq_conf["modification_mapping"]
        self._fixed_modifications = mq_conf["fixed_modifications"]
        self._path_to_new_mgf = None
        self._msms = None

    @property
    def path_to_mgf_file(self):
        """Get path to single unified MGF file."""
        if not self._path_to_new_mgf:
            logging.warning(
                "`_path_to_new_mgf` is not set yet; first run the `parse_mgf_files` "
                "method"
            )
        else:
            return self._path_to_new_mgf

    @property
    def original_pin(self):
        """Get PercolatorIn object from identification file."""
        raise NotImplementedError(
            "Property `original_pin` is not implemented in class `MaxQuantPipeline`."
        )

    def peprec_from_pin(self):
        """Get PeptideRecord from PIN file and MGF file."""
        raise NotImplementedError(
            "Method `peprec_from_pin` is not implemented in class `MaxQuantPipeline`."
        )

    @property
    def msms(self):
        """Get msms.txt identification results."""
        if self._msms is None:
            self._msms = MSMS.from_file(self.path_to_id_file)
        return self._msms

    def parse_mgf_files(self, peprec):
        """Parse multiple MGF files into one for MS²PIP."""
        logging.debug("Parsing MGF files into one for MS²PIP")
        path_to_new_mgf = self.output_basename + "unified.mgf"
        parse_mgf(
            peprec.df,
            self.passed_mgf_path,
            outname=path_to_new_mgf,
            filename_col='Raw file', spec_title_col='spec_id',
            title_parsing_method='TRFP_MQ',
        )
        self._path_to_new_mgf = path_to_new_mgf

    def get_peprec(self, parse_mgf: bool = True) -> PeptideRecord:
        """Get PeptideRecord from msms.txt file, optionally parse MGF files into one."""
        logging.debug("Converting MaxQuant msms.txt to PEPREC...")
        peprec = self.msms.to_peprec(
            modification_mapping=self._modification_mapping,
            fixed_modifications=self._fixed_modifications
        )
        if parse_mgf:
            self.parse_mgf_files(peprec)
        peprec.df.drop("Raw file", axis=1, inplace=True)
        return peprec

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        return self.msms.get_search_engine_features()


class PeptideShakerPipeline(_Pipeline):
    """PeptideShaker extended PSM report to PeptideRecord and search engine features."""

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        super().__init__(config, output_basename)

        # Private attributes, specific to this pipeline
        self._extended_psm_report = None

    @property
    def original_pin(self):
        """Get PercolatorIn object from identification file."""
        raise NotImplementedError(
            "Property `original_pin` is not implemented in class "
            "`PeptideShakerPipeline`."
        )

    def peprec_from_pin(self):
        """Get PeptideRecord from PIN file and MGF file."""
        raise NotImplementedError(
            "Method `peprec_from_pin` is not implemented in class "
            "`PeptideShakerPipeline`."
        )

    @property
    def extended_psm_report(self):
        """Get Extended PSM Report with identification results."""
        if self._extended_psm_report is None:
            self._extended_psm_report = ExtendedPsmReport.from_tsv(self.path_to_id_file)
        return self._extended_psm_report

    def get_peprec(self) -> PeptideRecord:
        """Get PeptideRecord."""
        return self.extended_psm_report.to_peprec()

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        return self.extended_psm_report.get_search_engine_features()


class CometPipeline(_Pipeline):
    """Comet txt file to peprec and search engine features."""

    default_columns = {
        "Spectrum",
        "ScanNr",
        'Charge',
        'RT',
        'HLA.Class',
        'Peptide',
        'Sequence',
        'Proteins',
        'IsVariant',
        'IsDecoy',
        'Comet.Rank',
        'Comet.XCorr',
        'Comet.DeltaCn',
        'Comet.SpScore',
        'Comet.NegLogPv',
        'Comet.massdiff',
        'Comet.tot_num_ions',
        'Comet.num_matched_ions',
        'Comet.lFDR',
    }

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        super().__init__(config, output_basename)

        # Private attributes, specific to this pipeline

    @property
    def original_pin(self):
        """Get PercolatorIn object from identification file."""
        raise NotImplementedError(
            "Property `original_pin` is not implemented in class "
            "`CometPipeline`."
        )

    def peprec_from_pin(self):
        """Get PeptideRecord from PIN file and MGF file."""
        raise NotImplementedError(
            "Method `peprec_from_pin` is not implemented in class "
            "`CometPipeline`."
        )

    @staticmethod
    def _get_peprec_modifications(sequences, mods_requiring_suffix=None) -> List:
        """Get peprec-formatted modifications."""
        if not mods_requiring_suffix:
            mods_requiring_suffix = []
        parsed_modifications = []
        for sequence in sequences:
            mod = []
            if "(" not in sequence:
                parsed_modifications.append("-")
            else:
                while re.match(r".*\(([^)]*)\).*", sequence):
                    x = re.search(r"\(([^)]*)\)", sequence)
                    loc = int(x.start())
                    name = x.group().strip("()")
                    if name in mods_requiring_suffix:
                        name = name + sequence[loc-1]
                    mod.extend([str(loc), name])
                    sequence = re.sub(r"\(([^)]*)\)", "", sequence, 1)
                mod = "|".join(mod)
                parsed_modifications.append(mod)
        return parsed_modifications

    def get_peprec(self):

        peprec = pd.DataFrame(
           columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "protein_list",
                "psm_score",
                "observed_retention_time",
                "Label",
                "Raw file"
            ]
        )
        comet_df = pd.read_table(self.path_to_id_file, sep="\t")
        peprec["spec_id"] = "controllerType=0 controllerNumber=1 scan=" + comet_df[
            "ScanNr"
        ].astype(str)
        peprec["peptide"] = comet_df["Sequence"]
        peprec["modifications"] = comet_df["Peptide"].apply(
            CometPipeline._get_peprec_modifications
            )
        peprec["charge"] = comet_df["charge"]
        peprec["protein_list"] = comet_df["Proteins"]
        peprec["psm_score"] = comet_df["comet.SpScore"]
        peprec["observed_retention_time"] = comet_df["RT"]
        peprec["Label"] = comet_df["IsDecoy"].apply(lambda x: 1 if x else -1)
        raw_files = []
        for i in comet_df.Spectrum:
            raw_file, _, _ = i.partition(".")
            raw_files.append(raw_file)
        peprec["Raw File"] = raw_files

        return peprec

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        raise NotImplementedError(
            "Method `get_search_engine_features` is not implemented in class "
            "`CometPipeline`."
        )


class MzidPipeline(_Pipeline):
    """Mzid identification output file from PEAKS to PeptideRecord
    and search engine features."""

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        super().__init__(config, output_basename)

        # Private attributes, specific to this pipeline

    @property
    def original_pin(self):
        """Get PercolatorIn object from identification file."""
        raise NotImplementedError(
            "Property `original_pin` is not implemented in class "
            "`MzidPipeline`."
        )

    def peprec_from_pin(self):
        """Get PeptideRecord from PIN file and MGF file."""
        raise NotImplementedError(
            "Method `peprec_from_pin` is not implemented in class "
            "`MzidPipeline`."
        )

    @staticmethod
    def _get_peprec_modifications(modifications: List):
        if isinstance(modifications, float):
            return "-"
        elif isinstance(modifications, List):
            suffix_list = ["Phospho"]
            return "|".join(f"{b['location']}|{b['name'] + b['residues'][0]}"
                            if b['name'] in suffix_list
                            else f"{b['location']}|{b['name']}" for b in modifications)
        else:
            raise TypeError

    def get_peprec(self) -> PeptideRecord:

        peprec = pd.DataFrame(
           columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "protein_list",
                "psm_score",
                "observed_retention_time",
                "Label",
                "Raw file"
            ]
        )
        mzid_df = pd.read_table(self.path_to_id_file, sep="\t")
        peprec["spec_id"] = mzid_df.spectrum_ID
        [lambda x: ''.join(re.findall(r"\d+", x))]
        peprec["peptide"] = mzid_df["PeptideSequence"]
        peprec["modifications"] = mzid_df["Modifications"].apply(
            MzidPipeline._get_peprec_modifications)
        peprec["charge"] = mzid_df["chargeState"]
        peprec["protein_list"] = mzid_df["accession"]
        peprec["psm_score"] = mzid_df["PEAKS:peptideScore"]

        return peprec

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        raise NotImplementedError(
            "Method `get_search_engine_features` is not implemented in class "
            "`MzidPipeline`."
        )


class SpectrumMillPipeline(_Pipeline):
    """Pipeline for Spectrum Mill identification output files(.ssv)."""

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        super().__init__(config, output_basename)

        # Private attributes, specific to this pipeline

    @property
    def original_pin(self):
        """Get PercolatorIn object from identification file."""
        raise NotImplementedError(
            "Property `original_pin` is not implemented in class "
            "`SpectrumMillPipeline`."
        )

    def peprec_from_pin(self):
        """Get PeptideRecord from PIN file and MGF file."""
        raise NotImplementedError(
            "Method `peprec_from_pin` is not implemented in class "
            "`SpectrumMillPipeline`."
        )

    @staticmethod
    def _get_peprec_modifications(sequence: str):

        modification_map = {
                        'm': ('M', 'Oxidation'),
                        'q': ('Q', 'Glu->pyro-Glu'),
                        'n': ('N', 'Deamidated'),
                        'y': ('Y', 'PhosphoY'),
                        't': ('T', 'PhosphoT'),
                        's': ('S', 'PhosphoS')
                    }
        modification = []
        if sequence.isupper():
            return sequence, "-"
        else:
            while re.match(".*[a-z].*", sequence):
                lc = re.search("[a-z]", sequence)
                modification.append(
                    [str(lc.start()+1), modification_map[lc.group()][1]]
                )
                sequence = re.sub('[a-z]', modification_map[lc.group()][0], sequence, 1)
            modification = ["|".join(tups) for tups in modification]
            modification = "|".join(modification)
        return sequence, modification

    @staticmethod
    def _get_filename_and_scannumber(filename: str, separator: str):
        rawfilename, _, scannumber = filename.partition(".")
        scannumber = scannumber.split('.')
        scannumber = scannumber[0]
        return rawfilename, scannumber

    def get_peprec(self) -> PeptideRecord:
        peprec = pd.DataFrame(
           columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "protein_list",
                "psm_score",
                "observed_retention_time",
                "Label",
                "Raw file"
            ]
        )
        ssv_df = pd.read_table(self.path_to_id_file, sep="\t")
        peprec[['raw file', 'spec_id']] = pd.DataFrame(
            ssv_df['filename']
            .apply(SpectrumMillPipeline._get_filename_and_scannumber)
            .tolist(), index=ssv_df.index)
        peprec[["peptide", "modifications"]] = pd.DataFrame(
            ssv_df.sequence.apply(SpectrumMillPipeline._get_peprec_modifications)
            .tolist(), index=peprec.index)
        peprec['charge'] = ssv_df["parent_charge"]
        peprec["protein_list"] = ssv_df["accession_numbers"]
        peprec["psm_score"] = ssv_df["score"]
        peprec["observed_retention_time"] = ssv_df["retentionTimeMin"]

        return peprec

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        raise NotImplementedError(
            "Method `get_search_engine_features` is not implemented in class "
            "`SpectrumMillPipeline`."
        )
