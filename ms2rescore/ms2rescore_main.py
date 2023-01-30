import json
import logging
import os
import re
import subprocess
import tempfile
from glob import glob
from multiprocessing import cpu_count
from pathlib import Path
from typing import Dict, Optional, Union

import psm_utils.io
from pandas.errors import EmptyDataError
from rich.console import Console

from ms2rescore import setup_logging
from ms2rescore.config_parser import parse_config
from ms2rescore.exceptions import MS2RescoreConfigurationError, MS2RescoreError
from ms2rescore.feature_generators.ms2pip import MS2PIPFeatureGenerator

logger = logging.getLogger(__name__)

id_file_parser = None


class MS2Rescore:
    """
    MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs.

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
        self.config = parse_config(
            parse_cli_args=parse_cli_args, config_class=configuration
        )
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
            self.config["ms2rescore"]["num_cpu"],
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

        logger.debug("Parsing modifications...")
        psm_list.set_ranks(lower_score_better=False)  # TODO make config option
        psm_list = psm_list.get_rank1_psms()
        psm_list.rename_modifications(self.config["ms2rescore"]["modification_mapping"])
        psm_list.add_fixed_modifications(
            self.config["ms2rescore"]["fixed_modifications"]
        )
        psm_list.apply_fixed_modifications()

        logger.debug("Applying `psm_id_pattern`...")
        if self.config["ms2rescore"]["psm_id_pattern"]:
            pattern = re.compile(self.config["ms2rescore"]["psm_id_pattern"])

            def _match_ids(old_id):
                match = re.search(pattern, str(old_id))
                try:
                    return match[1]
                except (TypeError, IndexError):
                    raise MS2RescoreError(
                        "`psm_id_pattern` could not be matched to all PSM spectrum IDs."
                        " Are you sure that the regex contains a capturing group?"
                    )

            new_ids = [_match_ids(old_id) for old_id in psm_list["spectrum_id"]]
            psm_list["spectrum_id"] = new_ids

        # fgen = MS2PIPFeatureGenerator(config=self.config)
        # fgen.add_features(psm_list)

        psm_utils.io.write_file(
            psm_list,
            filename="test.pin",
            filetype="percolator",
            style="pin",
            feature_names=psm_list[0].rescoring_features.keys(),
        )

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
