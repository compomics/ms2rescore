#! python
"""MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and retention times."""

# Standard library
import logging
import subprocess
import os
from multiprocessing import cpu_count
from typing import Optional, Union, Dict

# From package
from ms2rescore.config_parser import parse_config
from ms2rescore.retention_time import RetentionTimeIntegration
from ms2rescore._version import __version__

from ms2rescore import (
    setup_logging,
    rescore_core,
    pin_to_rescore,
    maxquant_to_rescore,
    msgf_to_rescore,
    tandem_to_rescore,
    peptideshaker_to_rescore,
)


class MS2ReScore():
    """
    MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and retention times.

    Parameters
    ----------
    parse_cli_args : bool, optional
        parse command line arguments, default True
    configuration : dict
        dict containing general ms2rescore configuration. Should at least contain
        `identification_file`. Required if `parse_cli_args` is False.
    """

    def __init__(
        self,
        parse_cli_args: bool = True,
        configuration: Optional[Dict] = None,
        set_logger: bool = True
    ) -> None:
        """Initialize MS2ReScore object."""
        self.config = parse_config(
            parse_cli_args=parse_cli_args,
            config_class=configuration
        )
        if set_logger:
            setup_logging.setup_logging(self.config["general"]["log_level"])

        self._validate_cli_dependency("percolator -h")
        self._validate_cli_dependency("ms2pip -h")

        logging.debug(
            "Using %i of %i available CPUs.",
            self.config["general"]["num_cpu"],
            cpu_count()
        )

        self.pipeline = self._select_pipeline(self.config)
        logging.info("Using %s pipeline", self.config["general"]["pipeline"])

    @staticmethod
    def _validate_cli_dependency(command):
        """Validate that command returns zero exit status."""
        if subprocess.getstatusoutput(command)[0] != 0:
            logging.critical(
                "`%s` returned non-zero exit status. Please verify installation.",
                command
            )
            exit(1)

    @staticmethod
    def _select_pipeline(config):
        """Select specific rescoring pipeline."""
        if config["general"]["pipeline"] == "pin":
            pipeline = pin_to_rescore.pipeline
        elif config["general"]["pipeline"] == "maxquant":
            pipeline = maxquant_to_rescore.maxquant_pipeline
        elif config["general"]["pipeline"] == "msgfplus":
            pipeline = msgf_to_rescore.msgf_pipeline
        elif config["general"]["pipeline"] == "tandem":
            pipeline = tandem_to_rescore.tandem_pipeline
        elif config["general"]["pipeline"] == "peptideshaker":
            pipeline = peptideshaker_to_rescore.pipeline
        else:
            raise NotImplementedError(config["general"]["pipeline"])
        return pipeline

    @staticmethod
    def get_ms2pip_features(
        ms2pip_config: Dict,
        peprec_filename: Union[str, os.PathLike],
        mgf_filename: Union[str, os.PathLike],
        ms2pip_config_filename: Union[str, os.PathLike],
        output_filename: Union[str, os.PathLike],
        preds_filename: Union[str, os.PathLike],
        num_cpu: int
    ):
        """Get predicted MS² peak intensities from MS2PIP."""
        logging.info("Adding MS2 peak intensity features with MS²PIP.")
        rescore_core.make_ms2pip_config(ms2pip_config, filename=ms2pip_config_filename)
        ms2pip_command = "ms2pip {} -c {} -s {} -n {}".format(
            peprec_filename,
            ms2pip_config_filename,
            mgf_filename,
            num_cpu,
        )

        logging.debug("Running MS2PIP: %s", ms2pip_command)
        subprocess.run(ms2pip_command, shell=True, check=True)

        logging.info("Calculating features from predicted spectra")
        rescore_core.calculate_features(
            preds_filename,
            output_filename + "_ms2pipfeatures.csv",
            num_cpu,
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
            peprec_filename,
            output_filename + "_rtfeatures.csv",
            num_cpu=num_cpu,
        )
        rt_int.run()

    @staticmethod
    def _remove_tmp_files(tmp_files):
        logging.debug("Removing temporary files")
        for filename in tmp_files:
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

    def _run_percolator(self, output_filename: Union[str, os.PathLike]):
        """Run Percolator with different feature subsets."""
        for subset in self.config["general"]["feature_sets"]:
            subname = output_filename + "_" + subset + "features"
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

        outname = self.config["general"]["output_filename"]
        ms2pip_config_filename = outname + "_ms2pip_config.txt"
        preds_filename = (
            peprec_filename.replace(".peprec", "")
            + "_"
            + self.config["ms2pip"]["model"]
            + "_pred_and_emp.csv"
        )

        if any(
            fset in self.config["general"]["feature_sets"]
            for fset in ["ms2pip", "all", "ms2pip_rt"]
        ):
            self.get_ms2pip_features(
                self.config["ms2pip"],
                peprec_filename,
                mgf_filename,
                ms2pip_config_filename,
                outname,
                preds_filename,
                self.config["general"]["num_cpu"]
            )

        if any(
            fset in self.config["general"]["feature_sets"] for fset in ["rt", "all", "ms2pip_rt"]
        ):
            self.get_rt_features(
                peprec_filename, outname, self.config["general"]["num_cpu"]
            )

        logging.info("Generating PIN files")
        rescore_core.write_pin_files(
            peprec_filename,
            outname,
            ms2pip_features_path=outname + "_ms2pipfeatures.csv",
            rt_features_path=outname + "_rtfeatures.csv",
            feature_sets=self.config["general"]["feature_sets"],
        )

        if not self.config["general"]["keep_tmp_files"]:
            tmp_files = [
                ms2pip_config_filename,
                preds_filename,
                outname + "_ms2pipfeatures.csv",
                outname + "_" + self.config["ms2pip"]["model"] + "_correlations.csv",
                outname + ".mgf",
                outname + ".peprec",
                outname + "_rtfeatures.csv",
            ]
            self._remove_tmp_files(tmp_files)

        if self.config["general"]["run_percolator"]:
            self._run_percolator(outname)

        logging.info("MS²ReScore finished!")
