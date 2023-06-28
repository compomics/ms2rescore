import logging
import subprocess
from pathlib import Path

import psm_utils

from ms2rescore.rescoring_engines import Rescoringengine
from ms2rescore.exceptions import MS2RescoreError

logger = logging.getLogger(__name__)


LOG_LEVEL_MAP = {
    "critical": 0,
    "error": 0,
    "warning": 0,
    "info": 1,
    "debug": 2,
}


class PercolatorRescoring(Rescoringengine):
    def __init__(self, psm_list, config, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        self.psm_list = psm_list
        self.config = config
        self.output_file_root = str(
            Path(self.config["ms2rescore"]["output_path"])
            / Path(self.config["ms2rescore"]["psm_file"]).stem
        )  # TODO add feature generators to the name?
        self.kwargs = {
            "results-psms": self.output_file_root + "_target_psms.pout",
            "decoy-results-psms": self.output_file_root + "_decoy_psms.pout",
            "results-peptides": self.output_file_root + "_target_peptides.pout",
            "decoy-results-peptides": self.output_file_root + "_decoy_peptides.pout",
            "weights": self.output_file_root + ".weights",
            "verbose": LOG_LEVEL_MAP[self.config["ms2rescore"]["log_level"]],
            "num-threads": self.config["ms2rescore"]["processes"],
            "post-processing-tdc": True,
        }

    def rescore(self):
        """Run Percolator"""
        self.write_pin_file()
        try:
            output = subprocess.run(
                self.parse_percolator_command(),
                capture_output=True,
            )
        except subprocess.CalledProcessError:
            logger.warn(f"Percolator was not run properly:\n {output.stdout}")
            raise MS2RescoreError("Percolator error")
        stderr = self._decode_string(output.stderr)

        logger.info(
            "Percolator output: \n" + stderr,
            extra={"highlighter": None},
        )

    def write_pin_file(self):
        """Write pin file for rescoring"""
        logger.debug(f"Writing pin file to {self.output_file_root}.pin")
        psm_utils.io.write_file(
            self.psm_list,
            filename=f"{self.output_file_root}.pin",
            filetype="percolator",
            style="pin",
            feature_names=self.psm_list[0].rescoring_features.keys(),
        )

    def parse_percolator_command(self):
        """Create Percolator command"""

        self.kwargs.update(self.config["percolator"])

        percolator_cmd = ["percolator"]
        for key, value in self.kwargs.items():
            if not isinstance(value, bool):
                percolator_cmd.append(f"--{key}")
                percolator_cmd.append(str(value))
                if key == "init-weights":
                    percolator_cmd.append("--static")
            elif isinstance(value, bool) & value == False:
                continue
            else:
                percolator_cmd.append(f"--{key}")

        percolator_cmd.append(f"{self.output_file_root}.pin")
        logger.debug(f"Running percolator command {' '.join(percolator_cmd)}")

        return percolator_cmd

    @staticmethod
    def _decode_string(encoded_string):
        encodings = [
            "utf-8",
            "latin-1",
            "ascii",
            "iso-8859-15",
        ]

        for encoding in encodings:
            try:
                decoded_string = encoded_string.decode(encoding)
                logger.debug(f"Decoded stderr with {encoding}")
                return decoded_string
            except UnicodeDecodeError:
                pass

        raise MS2RescoreError("Could not infer encoding of Percolator output")
