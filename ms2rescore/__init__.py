#! python
"""MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs."""

import logging
import os
import subprocess
import tempfile
from multiprocessing import cpu_count
from typing import Dict, Optional, Union

from pandas.errors import EmptyDataError

from ms2pip.ms2pipC import MS2PIP

from ms2rescore import id_file_parser, plotting, rescore_core, setup_logging
from ms2rescore._exceptions import MS2RescoreError
from ms2rescore._version import __version__
from ms2rescore.config_parser import parse_config
from ms2rescore.retention_time import RetentionTimeIntegration

logger = logging.getLogger(__name__)


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

        if self.config["general"]["run_percolator"]:
            self._validate_cli_dependency("percolator -h")

        logger.debug(
            "Using %i of %i available CPUs.",
            self.config["general"]["num_cpu"],
            cpu_count(),
        )

        if not self.config["general"]["tmp_path"]:
            self.tmp_path = tempfile.mkdtemp()
            self.config["general"]["tmp_path"] = self.tmp_path
        else:
            self.tmp_path = self.config["general"]["tmp_path"]
            os.makedirs(self.tmp_path, exist_ok=True)

        self.tmpfile_basepath = os.path.join(
            self.tmp_path,
            os.path.basename(
                os.path.splitext(self.config["general"]["identification_file"])[0]
            ),
        )

        selected_pipeline = self._select_pipeline()
        self.pipeline = selected_pipeline(self.config, self.tmpfile_basepath)
        logger.info("Using %s.", selected_pipeline.__name__)

    @staticmethod
    def _validate_cli_dependency(command):
        """Validate that command returns zero exit status."""
        if subprocess.getstatusoutput(command)[0] != 0:
            logger.critical(
                "`%s` returned non-zero exit status. Please verify installation.",
                command,
            )
            exit(1)

    @staticmethod
    def _infer_pipeline(identification_file: str):
        """Infer pipeline from identification file."""
        logger.debug("Inferring pipeline from identification filename...")
        if identification_file.lower().endswith(".pin"):
            pipeline = id_file_parser.PinPipeline
        elif identification_file.lower().endswith(".t.xml"):
            pipeline = id_file_parser.TandemPipeline
        elif identification_file.endswith("msms.txt"):
            pipeline = id_file_parser.MaxQuantPipeline
        elif identification_file.lower().endswith(".mzid"):
            pipeline = id_file_parser.MSGFPipeline
        else:
            raise MS2RescoreError(
                "Could not infer pipeline from identification filename. Please specify "
                "`general` > `pipeline` in your configuration file."
            )
        return pipeline

    def _select_pipeline(self):
        """Select specific rescoring pipeline."""
        if self.config["general"]["pipeline"] == "infer":
            pipeline = self._infer_pipeline(
                self.config["general"]["identification_file"]
            )
        elif self.config["general"]["pipeline"] == "pin":
            pipeline = id_file_parser.PinPipeline
        elif self.config["general"]["pipeline"] == "maxquant":
            pipeline = id_file_parser.MaxQuantPipeline
        elif self.config["general"]["pipeline"] == "msgfplus":
            pipeline = id_file_parser.MSGFPipeline
        elif self.config["general"]["pipeline"] == "tandem":
            pipeline = id_file_parser.TandemPipeline
        elif self.config["general"]["pipeline"] == "peptideshaker":
            pipeline = id_file_parser.PeptideShakerPipeline
        elif self.config["general"]["pipeline"] == "peaks":
            pipeline = id_file_parser.PeaksPipeline
        else:
            raise NotImplementedError(self.config["general"]["pipeline"])
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
        logger.info("Adding MS2 peak intensity features with MS²PIP.")
        ms2pip_config = rescore_core.make_ms2pip_config_dict(ms2pip_config)

        # Check if input files exist
        for f in [peprec_filename, mgf_filename]:
            if not os.path.isfile(f):
                raise FileNotFoundError(f)

        logger.debug("Running MS²PIP...")
        ms2pip = MS2PIP(
            peprec_filename,
            spec_file=mgf_filename,
            params=ms2pip_config,
            num_cpu=num_cpu,
            add_retention_time=False,
            compute_correlations=False,
        )
        try:
            ms2pip.run()
        finally:
            ms2pip.cleanup()
        logger.debug("MS²PIP finished.")

        logger.info("Calculating features from predicted spectra")
        preds_filename = (
            peprec_filename.replace(".peprec", "")
            + "_"
            + ms2pip_config["ms2pip"]["model"]
            + "_pred_and_emp.csv"
        )
        rescore_core.calculate_features(
            preds_filename,
            output_filename + "_ms2pipfeatures.csv",
            num_cpu,
            show_progress_bar=False,
        )

    @staticmethod
    def get_rt_features(
        peprec_filename: Union[str, os.PathLike],
        output_filename: Union[str, os.PathLike],
        num_cpu: int,
    ):
        """Get retention time features with DeepLC."""
        logger.info("Adding retention time features with DeepLC.")
        rt_int = RetentionTimeIntegration(
            peprec_filename,
            output_filename + "_rtfeatures.csv",
            num_cpu=num_cpu,
        )
        rt_int.run()

    def _run_percolator(self):
        """Run Percolator with different feature subsets."""


        for subset in self.config["general"]["feature_sets"]:
            subname = (
                self.config["general"]["output_filename"]
                + "_"
                + "_".join(subset)
                + "_features"
            )

            kwargs = {
                    "results-psms": subname + ".pout",
                    "decoy-results-psms": subname + ".pout_dec",
                    "weights": subname + ".weights",
                    "verbose": 0,
                    "post-processing-tdc": True
                    }
            kwargs.update(self.config["percolator"])

            percolator_cmd = ["percolator"]
            for key, value in kwargs.items():

                if not isinstance(value,bool):
                    percolator_cmd.append(f"--{key}")
                    percolator_cmd.append(str(value))
                elif isinstance(value, bool) & value==False:
                    continue
                else:
                    percolator_cmd.append(f"--{key}")
            percolator_cmd.append(subname+".pin")
            logger.info("Running Percolator: %s", " ".join(percolator_cmd))
            output = subprocess.run(percolator_cmd, capture_output=True, check=True)

            logger.debug(output.stdout)
            logger.info(output.stderr)


    def run(self):
        """Run MS²ReScore."""
        peprec = self.pipeline.get_peprec()
        peprec.validate_decoy_presence()
        peprec_filename = self.tmpfile_basepath + ".peprec"
        peprec.to_csv(peprec_filename)

        search_engine_features = self.pipeline.get_search_engine_features()
        search_engine_features_filename = (
            self.tmpfile_basepath + "_search_engine_features.csv"
        )
        search_engine_features.to_csv(search_engine_features_filename, index=False)

        if any("ms2pip" in fst for fst in self.config["general"]["feature_sets"]):
            self.get_ms2pip_features(
                self.config["ms2pip"],
                peprec_filename,
                self.pipeline.path_to_mgf_file,
                self.tmpfile_basepath,
                self.config["general"]["num_cpu"],
            )

        if any("rt" in fst for fst in self.config["general"]["feature_sets"]):
            self.get_rt_features(
                peprec_filename,
                self.tmpfile_basepath,
                self.config["general"]["num_cpu"],
            )

        logger.info("Generating PIN files")
        rescore_core.write_pin_files(
            peprec_filename,
            self.config["general"]["output_filename"],
            searchengine_features_path=search_engine_features_filename,
            ms2pip_features_path=self.tmpfile_basepath + "_ms2pipfeatures.csv",
            rt_features_path=self.tmpfile_basepath + "_rtfeatures.csv",
            feature_sets=self.config["general"]["feature_sets"],
        )

        if self.config["general"]["run_percolator"]:
            self._run_percolator()

            # Only use plotting module when run_percolator is true
            if self.config["general"]["plotting"]:
                logger.info("Generating Rescore plots")

                plotting.PIN(
                    peprec_filename, self.config["general"]["output_filename"]
                )

                for fset in self.config["general"]["feature_sets"]:
                    pout_file = (
                        self.config["general"]["output_filename"]
                        + "_"
                        + "_".join(fset)
                        + "_features.pout"
                    )
                    pout_decoy_file = (
                        self.config["general"]["output_filename"]
                        + "_"
                        + "_".join(fset)
                        + "_features.pout_dec"
                    )
                    try:
                        plotting.POUT(
                            pout_file,
                            pout_decoy_file,
                            self.config["general"]["output_filename"],
                            " ".join(fset)
                        )
                    except EmptyDataError:
                        logger.warn(f"Feature set: {'_'.join(fset)} returned empty pout file")
                        continue

                plotting.RescoreRecord.save_plots_to_pdf(
                    self.config["general"]["output_filename"] + "_plots.pdf",
                    FDR_thresholds=[0.01, 0.001],
                )
        if not self.config["general"]["run_percolator"] and self.config["general"]["plotting"]:
            logger.warn("To plot rescore results run_percolator should be enabled")
        logger.info("MS²ReScore finished!")
