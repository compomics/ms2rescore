import logging
import subprocess
from typing import Any, Dict, Optional

import psm_utils

from ms2rescore.exceptions import MS2RescoreError
from ms2rescore.rescoring_engines import RescoringEngineBase

logger = logging.getLogger(__name__)


LOG_LEVEL_MAP = {
    "critical": 0,
    "error": 0,
    "warning": 0,
    "info": 1,
    "debug": 2,
}


class PercolatorRescoring(RescoringEngineBase):
    def __init__(
        self,
        output_file_root: str,
        *args,
        log_level: str = "info",
        processes: int = 1,
        percolator_kwargs: Optional[Dict[str, Any]] = None,
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)

        self.output_file_root = output_file_root
        self.log_level = log_level
        self.processes = processes
        self.percolator_kwargs = {
            "results-psms": self.output_file_root + "_target_psms.pout",
            "decoy-results-psms": self.output_file_root + "_decoy_psms.pout",
            "results-peptides": self.output_file_root + "_target_peptides.pout",
            "decoy-results-peptides": self.output_file_root + "_decoy_peptides.pout",
            "weights": self.output_file_root + ".weights",
            "verbose": LOG_LEVEL_MAP[self.log_level],
            "num-threads": self.processes,
            "post-processing-tdc": True,
        }
        if percolator_kwargs:
            self.percolator_kwargs.update(percolator_kwargs)

    def rescore(self, psm_list: psm_utils.PSMList):
        """Run Percolator on a PSMList."""
        self.write_pin_file(psm_list, self.output_file_root)
        try:
            output = subprocess.run(
                self._parse_percolator_command(),
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

    @staticmethod
    def write_pin_file(psm_list: psm_utils.PSMList, output_file_root: str):
        """Write PIN file for rescoring"""
        logger.debug(f"Writing pin file to {output_file_root}.pin")
        psm_utils.io.write_file(
            psm_list,
            filename=output_file_root + ".pin",
            filetype="percolator",
            style="pin",
            feature_names=psm_list[0].rescoring_features.keys(),
        )

    def _parse_percolator_command(self):
        """Create Percolator command"""
        percolator_cmd = ["percolator"]
        for key, value in self.percolator_kwargs.items():
            if not isinstance(value, bool):
                percolator_cmd.append(f"--{key}")
                percolator_cmd.append(str(value))
                if key == "init-weights":
                    percolator_cmd.append("--static")
            elif isinstance(value, bool) & value is False:
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
        else:
            raise MS2RescoreError("Could not infer encoding of Percolator logs.")
