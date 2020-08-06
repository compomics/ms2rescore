#! python
"""MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs."""

# Standard library
import logging
from os import pipe
import subprocess
import os
import tempfile
from multiprocessing import Value, cpu_count
from typing import Optional, Union, Dict

# From package
from ms2rescore._version import __version__
from ms2rescore._exceptions import MS2ReScoreError
from ms2rescore.config_parser import parse_config
from ms2rescore.retention_time import RetentionTimeIntegration

from ms2rescore import (
    setup_logging,
    rescore_core,
    pin_to_rescore,
    maxquant_to_rescore,
    msgf_to_rescore,
    tandem_to_rescore,
    peptideshaker_to_rescore,
)


class MS2ReScore:
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
    ) -> None:
        """Initialize MS2ReScore object."""
        self.config = parse_config(
            parse_cli_args=parse_cli_args, config_class=configuration
        )
        if set_logger:
            setup_logging.setup_logging(self.config["general"]["log_level"])

        self._validate_cli_dependency("percolator -h")
        self._validate_cli_dependency("ms2pip -h")

        logging.debug(
            "Using %i of %i available CPUs.",
            self.config["general"]["num_cpu"],
            cpu_count(),
        )

        if not self.config["general"]["tmp_path"]:
            self.tmp_path = tempfile.TemporaryDirectory()
            self.config["general"]["tmp_path"] = self.tmp_path.name
        else:
            self.tmp_path = None

        self.tmpfile_basepath = os.path.join(
            self.config["general"]["tmp_path"],
            os.path.basename(
                os.path.splitext(self.config["general"]["identification_file"])[0]
            ),
        )

        self.pipeline = self._select_pipeline()
        logging.info("Using %s pipeline", self.config["general"]["pipeline"])

    @staticmethod
    def _validate_cli_dependency(command):
        """Validate that command returns zero exit status."""
        if subprocess.getstatusoutput(command)[0] != 0:
            logging.critical(
                "`%s` returned non-zero exit status. Please verify installation.",
                command,
            )
            exit(1)

    def _select_pipeline(self):
        """Select specific rescoring pipeline."""
        if self.config["general"]["pipeline"] == "infer":
            pipeline = self._infer_pipeline(
                self.config["general"]["identification_file"]
            )
        elif self.config["general"]["pipeline"] == "pin":
            pipeline = pin_to_rescore.pipeline
        elif self.config["general"]["pipeline"] == "maxquant":
            pipeline = maxquant_to_rescore.maxquant_pipeline
        elif self.config["general"]["pipeline"] == "msgfplus":
            pipeline = msgf_to_rescore.msgf_pipeline
        elif self.config["general"]["pipeline"] == "tandem":
            pipeline = tandem_to_rescore.tandem_pipeline
        elif self.config["general"]["pipeline"] == "peptideshaker":
            pipeline = peptideshaker_to_rescore.pipeline
        else:
            raise NotImplementedError(self.config["general"]["pipeline"])
        return pipeline

    @staticmethod
    def _infer_pipeline(identification_file: Union[str, os.PathLike]):
        """Infer pipeline from identification file."""
        if identification_file.lower().endswith(".pin"):
            pipeline = pin_to_rescore.pipeline
        elif identification_file.lower().endswith(".t.xml"):
            pipeline = tandem_to_rescore.tandem_pipeline
        elif identification_file == "msms.txt":
            pipeline = maxquant_to_rescore.maxquant_pipeline
        elif identification_file.lower().endswith(".mzid"):
            pipeline = msgf_to_rescore.msgf_pipeline
        else:
            raise MS2ReScoreError(
                "Could not infer pipeline from identification filename. Please specify "
                "`general` > `pipeline` in your configuration file."
            )
        return pipeline

    @staticmethod
    def get_ms2pip_features(
        ms2pip_config: Dict,
        peprec_filename: Union[str, os.PathLike],
        mgf_filename: Union[str, os.PathLike],
        output_filename: Union[str, os.PathLike],
        num_cpu: int,
    ):
        """Get predicted MS² peak intensities from MS2PIP."""
        logging.info("Adding MS2 peak intensity features with MS²PIP.")
        ms2pip_config_filename = output_filename + "_ms2pip_config.txt"
        rescore_core.make_ms2pip_config(ms2pip_config, filename=ms2pip_config_filename)
        ms2pip_command = "ms2pip {} -c {} -s {} -n {}".format(
            peprec_filename, ms2pip_config_filename, mgf_filename, num_cpu,
        )

        logging.debug("Running MS2PIP: %s", ms2pip_command)
        subprocess.run(ms2pip_command, shell=True, check=True)

        logging.info("Calculating features from predicted spectra")
        preds_filename = (
            peprec_filename.replace(".peprec", "")
            + "_"
            + ms2pip_config["model"]
            + "_pred_and_emp.csv"
        )
        rescore_core.calculate_features(
            preds_filename, output_filename + "_ms2pipfeatures.csv", num_cpu,
        )

    @staticmethod
    def get_rt_features(
        peprec_filename: Union[str, os.PathLike],
        output_filename: Union[str, os.PathLike],
        num_cpu: int,
    ):
        """Get retention time features with DeepLC."""
        logging.info("Adding retention time features with DeepLC.")
        rt_int = RetentionTimeIntegration(
            peprec_filename, output_filename + "_rtfeatures.csv", num_cpu=num_cpu,
        )
        rt_int.run()

    def _run_percolator(self):
        """Run Percolator with different feature subsets."""
        for subset in self.config["general"]["feature_sets"]:
            subname = (
                self.config["general"]["output_filename"] + "_" + subset + "features"
            )
            percolator_cmd = "percolator "
            for op in self.config["percolator"].keys():
                percolator_cmd = percolator_cmd + "--{} {} ".format(
                    op, self.config["percolator"][op]
                )
            percolator_cmd = (
                percolator_cmd
                + "{} -m {} -M {} -w {} -v 0 -U --post-processing-tdc\n".format(
                    subname + ".pin",
                    subname + ".pout",
                    subname + ".pout_dec",
                    subname + ".weights",
                )
            )

            logging.info("Running Percolator: %s", percolator_cmd)
            subprocess.run(percolator_cmd, shell=True)

            if not os.path.isfile(subname + ".pout"):
                logging.error("Error running Percolator")

    def run(self):
        """Run MS²ReScore."""
        peprec_filename, mgf_filename = self.pipeline(self.config)

        if any(
            fset in self.config["general"]["feature_sets"]
            for fset in ["ms2pip", "all", "ms2pip_rt"]
        ):
            self.get_ms2pip_features(
                self.config["ms2pip"],
                peprec_filename,
                mgf_filename,
                self.tmpfile_basepath,
                self.config["general"]["num_cpu"],
            )

        if any(
            fset in self.config["general"]["feature_sets"]
            for fset in ["rt", "all", "ms2pip_rt"]
        ):
            self.get_rt_features(
                peprec_filename,
                self.tmpfile_basepath,
                self.config["general"]["num_cpu"],
            )

        logging.info("Generating PIN files")
        rescore_core.write_pin_files(
            peprec_filename,
            self.config["general"]["output_filename"],
            ms2pip_features_path=self.tmpfile_basepath + "_ms2pipfeatures.csv",
            rt_features_path=self.tmpfile_basepath + "_rtfeatures.csv",
            feature_sets=self.config["general"]["feature_sets"],
        )

        if self.config["general"]["run_percolator"]:
            self._run_percolator()

        logging.info("MS²ReScore finished!")
