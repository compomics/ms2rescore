"""Extract PeptideRecord and search engine features from identification file."""

import logging
import os
import re
from abc import ABC, abstractmethod
from typing import Dict, Tuple, Union, List

import numpy as np
import pandas as pd
from pyteomics import tandem, mzid
from tqdm import tqdm

from ms2rescore._exceptions import MS2RescoreError
from ms2rescore.maxquant import MSMSAccessor
from ms2rescore.parse_mgf import parse_mgf
from ms2rescore.peptide_record import PeptideRecord
from ms2rescore.peptideshaker import ExtendedPsmReportAccessor
from ms2rescore.percolator import PercolatorIn, run_percolator_converter


logger = logging.getLogger(__name__)


class IDFileParserError(MS2RescoreError):
    """Error parsing ID file."""

    pass


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
        self.mgf_title_pattern = config["maxquant_to_rescore"]["mgf_title_pattern"]

        # General private attributes
        self._pin_spec_id_patterns = {
            "generic": r".+_([0-9]+)_[0-9]+_[0-9]+",
            "tandem": r".+_([0-9]+)_[0-9]+_[0-9]+",
            "msgfplus": r".+_SII_([0-9]+)_[0-9]+_[0-9]+_[0-9]+",
            # "USI": r"mzspec:PXD[0-9]{6}:[^\s\:]*:scan:([0-9]+)"
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
        expected_ext: str = ".mgf",
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
                passed_path, expected_rootname + expected_ext
            )

        # If passed path is file, use that, but warn if basename doesn't match expected
        elif os.path.isfile(passed_path):
            passed_rootname = os.path.splitext(os.path.basename(passed_path))[0]
            if passed_rootname != expected_rootname:
                logger.debug(
                    "Passed MGF name root `%s` does not match MGF name root `%s` from "
                    "identifications file. Continuing with passed MGF name.",
                    passed_rootname,
                    expected_rootname,
                )
            path_to_mgf_file = passed_path

        else:
            raise IDFileParserError(
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
            log_level=self.log_level,
        )

    @property
    def path_to_mgf_file(self) -> Union[str, os.PathLike]:
        """Get path to MGF file, inferred from mgf_path (which can also be to a dir)."""
        path_to_mgf_file = self._validate_mgf_path(
            self.passed_mgf_path,
            os.path.dirname(self.path_to_id_file),
            os.path.basename(os.path.splitext(self.path_to_id_file)[0]),
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
                self.modification_list, label_style=self._pin_modification_style
            )
        return self._original_pin

    def peprec_from_pin(self) -> PeptideRecord:
        """Get PeptideRecord from PIN file and MGF file."""
        # Get peprec
        peprec = self.original_pin.to_peptide_record(
            spectrum_index_pattern=self._pin_spec_id_patterns[self._pin_spec_id_style]
        )
        titles, retention_times = parse_mgf_title_rt(self.path_to_mgf_file)
        peprec.df["spec_id"] = peprec.df["spec_id"].map(titles)
        if "observed_retention_time" not in peprec.df.columns:
            # Map MGF titles and observed retention times
            peprec.df["observed_retention_time"] = peprec.df["spec_id"].map(
                retention_times
            )

        if not ~peprec.df["observed_retention_time"].isna().any():
            raise IDFileParserError(
                "Could not map all MGF retention times to spectrum indices."
            )
        if not ~peprec.df["spec_id"].isna().any():
            raise IDFileParserError("Could not map all MGF titles to spectrum indices.")

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
            self.original_pin.get_spectrum_filename(),
        )
        return path_to_mgf_file

    def get_peprec(self) -> PeptideRecord:
        """Convert X!Tandem XML file and PIN to PEPREC."""
        # Load tandem dataframe with Pyteomics
        logger.debug("Converting X!Tandem XML to PEPREC...")
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
            "id": "tandem_id",
        }

        peprec_df = tandem_df[tandem_peprec_mapping.keys()].rename(
            columns=tandem_peprec_mapping
        )
        # Set PSM score as -log(e-value)
        peprec_df["psm_score"] = -np.log(peprec_df["psm_score"])

        logger.debug("Adding modifications from original PIN to PEPREC...")
        pin = self.original_pin
        pin.add_peprec_modifications_column()
        pin.add_spectrum_index_column(label="tandem_id")
        pin.df["ModPeptide"] = pin.df["Peptide"]
        peprec_df = peprec_df.merge(
            pin.df[
                [
                    "modifications",
                    "tandem_id",
                    "hyperscore",
                    "Label",
                    "ModPeptide",
                    "Proteins",
                ]
            ],
            on="tandem_id",
        )
        # Validate merge by comparing the hyperscore columns
        if not (peprec_df["hyperscore_tandem"] == peprec_df["hyperscore"]).all():
            raise IDFileParserError(
                "Could not merge tandem xml and generated pin files."
            )
        peprec_df.drop(columns=["tandem_id", "hyperscore_tandem"], inplace=True)

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
        self._msms_df = None

    @property
    def path_to_mgf_file(self):
        """Get path to single unified MGF file."""
        if not self._path_to_new_mgf:
            logger.warning(
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
    def msms_df(self):
        """Get msms.txt identification results."""
        if self._msms_df is None:
            self._msms_df = pd.DataFrame.msms.from_file(self.path_to_id_file)
        return self._msms_df

    def parse_mgf_files(self, peprec):
        """Parse multiple MGF files into one for MS²PIP."""
        logger.debug("Parsing MGF files into one for MS²PIP")
        path_to_new_mgf = self.output_basename + "_unified.mgf"
        parse_mgf(
            peprec.df,
            self.passed_mgf_path,
            outname=path_to_new_mgf,
            filename_col="Raw file",
            spec_title_col="spec_id",
            title_parsing_method="run.scan.scan",
            mgf_title_pattern=self.mgf_title_pattern
        )
        self._path_to_new_mgf = path_to_new_mgf

    def get_peprec(self, parse_mgf: bool = True) -> PeptideRecord:
        """Get PeptideRecord from msms.txt file, optionally parse MGF files into one."""
        logger.debug("Converting MaxQuant msms.txt to PEPREC...")
        peprec = self.msms_df.msms.to_peprec(
            modification_mapping=self._modification_mapping,
            fixed_modifications=self._fixed_modifications,
        )
        if parse_mgf:
            self.parse_mgf_files(peprec)
        return peprec

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        return self.msms_df.msms.get_search_engine_features()


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
            self._extended_psm_report = pd.DataFrame.ext_psm_report.from_file(
                self.path_to_id_file
            )
        return self._extended_psm_report

    def get_peprec(self) -> PeptideRecord:
        """Get PeptideRecord."""
        return self.extended_psm_report.ext_psm_report.to_peprec()

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""
        return self.extended_psm_report.ext_psm_report.get_search_engine_features()


class PeaksPipeline(_Pipeline):
    # TODO: move code to separate file and create mzid class
    """Peaks mzid report to PeptideRecord and search engine features."""

    def __init__(self, config: Dict, output_basename: Union[str, os.PathLike]) -> None:
        super().__init__(config, output_basename)

        # Private attributes, specific to this pipeline
        self.df = None

    @property
    def original_pin(self):
        """Get PercolatorIn object from identification file."""
        raise NotImplementedError(
            "Property `original_pin` is not implemented in class `PeaksPipeline`."
        )

    def peprec_from_pin(self):
        """Get PeptideRecord from PIN file and MGF file."""
        raise NotImplementedError(
            "Method `peprec_from_pin` is not implemented in class " "`PeaksPipeline`."
        )

    @staticmethod
    def _get_peprec_modifications(modifications: List):
        # TODO: Add unit tests for this function
        """get peprec modifications out of the peaks id file."""
        suffix_list = ["Phospho", "Dioxidation"]
        if isinstance(modifications, List):
            mods = []
            for m in modifications:
                if m["name"] in suffix_list:
                    mods.append(f"{m['location']}|{m['name'] + m['residues'][0]}")
                elif "residues" not in m.keys() and m["location"] != 0:
                    mods.append(f"-1|{m['name']}")
                else:
                    mods.append(f"{m['location']}|{m['name']}")
            return "|".join(mods)
        else:
            raise TypeError("`modifications` should be of type list.")

    def _convert_to_flat_dict(self, nested_dict, parent_key="", sep="_"):
        """Convert nested dict to flat dict with concatenated keys."""

        items = []
        for k, v in nested_dict.items():
            new_key = parent_key + sep + k if parent_key else k

            if isinstance(v, Dict):
                items.extend(self._convert_to_flat_dict(v, new_key, sep=sep))
            elif isinstance(v, List):
                if isinstance(v[0], Dict) and k != "Modification":
                    items.extend(self._convert_to_flat_dict(v[0], new_key, sep=sep))
                elif isinstance(v[0], Dict) and k == "Modification":
                    items.append((new_key, v))
                else:
                    items.append((new_key, v[0]))
            else:
                items.append((new_key, v))
        return items

    def read_df_from_mzid(self) -> pd.DataFrame:
        """Read mzid to Dataframe."""
        logger.info("Processing mzid file")
        psm_list = []
        with mzid.read(self.path_to_id_file) as reader:
            for spectrum_identification_result in tqdm(reader):
                psm = {}
                flat_dict = dict(
                    self._convert_to_flat_dict(spectrum_identification_result)
                )
                spec_id = (
                    flat_dict["location"]
                    .rsplit("/", 1)[1]
                    .split(".", 1)[0]
                    .replace(",", "_", 1)
                    + ":"
                    + flat_dict["spectrumID"]
                )
                psm["spec_id"] = spec_id
                psm["peptide"] = flat_dict[
                    "SpectrumIdentificationItem_PeptideSequence"
                ]
                psm["peptide_length"] = len(psm["peptide"])
                try:
                    psm["modifications"] = self._get_peprec_modifications(
                        flat_dict["SpectrumIdentificationItem_Modification"]
                    )
                except KeyError:
                    psm["modifications"] = "-"
                psm["charge"] = flat_dict[
                    "SpectrumIdentificationItem_chargeState"
                ]
                psm["protein_list"] = [
                    d["accession"]
                    for d in spectrum_identification_result[
                        "SpectrumIdentificationItem"
                    ][0]["PeptideEvidenceRef"]
                    if "accession" in d.keys()
                ]
                psm["PEAKS:peptideScore"] = flat_dict[
                    "SpectrumIdentificationItem_PEAKS:peptideScore"
                ]
                psm["Label"] = flat_dict[
                    "SpectrumIdentificationItem_PeptideEvidenceRef_isDecoy"
                ]
                psm["Raw file"] = (
                    flat_dict["location"]
                    .rsplit("/", 1)[1]
                    .split(".", 1)[0]
                    .replace(",", "_", 1)
                )
                psm["calculatedMassToCharge"] = flat_dict["SpectrumIdentificationItem_calculatedMassToCharge"]
                psm["experimentalMassToCharge"] = flat_dict["SpectrumIdentificationItem_experimentalMassToCharge"]

                psm_list.append(psm)
                
            df = pd.DataFrame(psm_list)

            df["Label"] = df["Label"].apply(lambda x: -1 if x else 1)
            df["dM"] = df["experimentalMassToCharge"].astype(float) - df["calculatedMassToCharge"].astype(float)
            df["absdM"] = abs(df["dM"])

            return df

    def parse_mgf_files(self, peprec):
        """Parse multiple MGF files into one for MS²PIP."""
        logger.debug("Parsing MGF files into one for MS²PIP")
        path_to_new_mgf = self.output_basename + "_unified.mgf"
        parse_mgf(
            peprec,
            self.passed_mgf_path,
            outname=path_to_new_mgf,
            filename_col="Raw file",
            spec_title_col="spec_id",
            title_parsing_method="run_index",
        )
        self.passed_mgf_path = path_to_new_mgf

    def get_peprec(self) -> PeptideRecord:
        """Get PeptideRecord."""

        if not self.df:
            self.df = self.read_df_from_mzid()
        peprec_df = self.df[
            [
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "protein_list",
                "PEAKS:peptideScore",
                "Label",
                "Raw file",
                "calculatedMassToCharge",
                "experimentalMassToCharge"
            ]
        ].rename({"PEAKS:peptideScore": "psm_score"}, axis=1)
        self.parse_mgf_files(peprec_df)
        titles, rt = parse_mgf_title_rt(self.passed_mgf_path)
        id_rt_dict = {
            "spec_id": list(titles.values()),
            "observed_retention_time": list(rt.values()),
        }
        id_rt_df = pd.DataFrame.from_dict(id_rt_dict)
        peprec_df = pd.merge(peprec_df, id_rt_df, on="spec_id", how="inner")
        
        return PeptideRecord.from_dataframe(peprec_df)

    def get_search_engine_features(self) -> pd.DataFrame:
        """Get pandas.DataFrame with search engine features."""

        peprec_cols = [
            "peptide",
            "modifications",
            "observed_retention_time",
            "Raw file",
            "protein_list",
            "Label",
        ]

        feature_cols = [col for col in self.df.columns if col not in peprec_cols]
        return self.df[feature_cols]
