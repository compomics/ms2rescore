import logging
import re
import subprocess
from multiprocessing import cpu_count
from pathlib import Path
from typing import Dict, Optional

import psm_utils.io
from rich.console import Console

from ms2rescore import setup_logging
from ms2rescore.config_parser import parse_config
from ms2rescore.exceptions import MS2RescoreConfigurationError, MS2RescoreError
from ms2rescore.feature_generators import FEATURE_GENERATORS
from ms2rescore.rescoring_engines.percolator import PercolatorRescoring

logger = logging.getLogger(__name__)

id_file_parser = None


class MS2Rescore:
    """
    MS²Rescore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs.

    Parameters
    ----------
    parse_cli_args : bool, optional
        parse command line arguments, default True
    configuration : dict, optional
        dict containing general ms2rescore configuration; should at least contain
        `identification_file`; required if `parse_cli_args` is False
    set_logger : bool, optional
        set custom logger or not, default False
    """

    def __init__(
        self,
        parse_cli_args: bool = True,
        configuration: Optional[Dict] = None,
        set_logger: bool = False,
        rich_console: Optional[Console] = None,
    ) -> None:
        """Initialize MS2ReScore object."""
        self.config = parse_config(parse_cli_args=parse_cli_args, config_class=configuration)

        # Set output and temporary paths
        self.output_path = self.config["ms2rescore"]["output_path"]
        self.output_file_root = str(
            Path(self.output_path) / Path(self.config["ms2rescore"]["psm_file"]).stem
        )
        self.tmp_path = self.config["ms2rescore"]["tmp_path"]
        self.tmp_file_root = str(
            Path(self.tmp_path) / Path(self.config["ms2rescore"]["psm_file"]).stem
        )

        # Set logger
        self._rich_console = rich_console or Console(record=True)
        self.log_level = self.config["ms2rescore"]["log_level"]
        if set_logger:
            setup_logging.setup_logging(
                self.log_level,
                log_file=self.output_file_root + "-ms2rescore-log.txt",
                rich_console=self._rich_console,
            )

        logger.debug(
            "Using %i of %i available CPUs.",
            self.config["ms2rescore"]["processes"],
            cpu_count(),
        )

    def run(self):
        # Read PSMs
        logger.info("Reading PSMs...")
        psm_list = psm_utils.io.read_file(
            self.config["ms2rescore"]["psm_file"],
            filetype=self.config["ms2rescore"]["psm_file_type"],
            show_progressbar=True,
        )

        logger.debug("Finding decoys...")
        if self.config["ms2rescore"]["id_decoy_pattern"]:
            psm_list.find_decoys(self.config["ms2rescore"]["id_decoy_pattern"])
        n_psms = len(psm_list)
        percent_decoys = sum(psm_list["is_decoy"]) / n_psms * 100
        logger.info(f"Found {n_psms} PSMs, of which {percent_decoys:.2f}% are decoys.")
        if not any(psm_list["is_decoy"]):
            raise MS2RescoreConfigurationError(
                "No decoy PSMs found. Please check if decoys are present in the PSM file and that "
                "the `id_decoy_pattern` option is correct."
            )

        logger.debug("Parsing modifications...")
        psm_list.rename_modifications(self.config["ms2rescore"]["modification_mapping"])
        psm_list.add_fixed_modifications(self.config["ms2rescore"]["fixed_modifications"])
        psm_list.apply_fixed_modifications()

        logger.debug("Applying `psm_id_pattern`...")
        if self.config["ms2rescore"]["psm_id_pattern"]:
            pattern = re.compile(self.config["ms2rescore"]["psm_id_pattern"])
            new_ids = [_match_psm_ids(old_id, pattern) for old_id in psm_list["spectrum_id"]]
            psm_list["spectrum_id"] = new_ids

        psm_list["spectrum_id"] = [str(spec_id) for spec_id in psm_list["spectrum_id"]]

        # Add rescoring features
        feature_names = dict()
        for fgen_name, fgen_config in self.config["ms2rescore"]["feature_generators"].items():
            # TODO: Handle this somewhere else, more generally? Warning required?
            if fgen_name == "maxquant" and not (psm_list["source"] == "msms").all():
                continue
            conf = self.config["ms2rescore"].copy()
            conf.update(fgen_config)
            fgen = FEATURE_GENERATORS[fgen_name](**conf)
            fgen.add_features(psm_list)
            feature_names[fgen_name] = fgen.feature_names

        # Filter out psms that do not have all added features
        all_feature_names = set([f for fgen in feature_names.values() for f in fgen])
        psm_list = psm_list[
            [(set(psm.rescoring_features.keys()) == all_feature_names) for psm in psm_list]
        ]

        if self.config["ms2rescore"]["USI"]:
            logging.debug(f"Creating USIs for {len(psm_list)} PSMs")
            psm_list["spectrum_id"] = [psm.get_usi(as_url=False) for psm in psm_list]

        if "percolator" in self.config["ms2rescore"]["rescoring_engine"]:
            logging.debug(f"Writing {self.output_file_root}.pin file")
            percolator = PercolatorRescoring(
                self.output_file_root,
                log_level=self.config["ms2rescore"]["log_level"],
                processes=self.config["ms2rescore"]["processes"],
                percolator_kwargs=self.config["ms2rescore"]["rescoring_engine"]["percolator"],
            )
            percolator.rescore(psm_list)

        elif "mokapot" in self.config["ms2rescore"]["rescoring_engine"]:
            raise NotImplementedError()

    @staticmethod
    def _validate_cli_dependency(command):
        """Validate that command returns zero exit status."""
        if subprocess.getstatusoutput(command)[0] != 0:
            logger.critical(
                "`%s` returned non-zero exit status. Please verify installation.",
                command,
            )
            exit(1)

    def save_log(self) -> None:
        """Save full rich-text log to HTML."""
        if self._rich_console:
            self._rich_console.save_html(
                self.output_file_root + "-ms2rescore-log.html",
            )
        else:
            logger.warning("Could not write logs to HTML: rich console is not defined.")


def _match_psm_ids(old_id, regex_pattern):
    """Match PSM IDs to regex pattern or raise Exception if no match present."""
    match = re.search(regex_pattern, str(old_id))
    try:
        return match[1]
    except (TypeError, IndexError):
        raise MS2RescoreError(
            "`psm_id_pattern` could not be matched to all PSM spectrum IDs."
            " Ensure that the regex contains a capturing group?"
        )
