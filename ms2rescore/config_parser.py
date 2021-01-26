"""Parse configuration from command line arguments and configuration files."""

import argparse
import json
import multiprocessing as mp
import os
from typing import Dict, Optional

from cascade_config import CascadeConfig

from ms2rescore import package_data
from ms2rescore._exceptions import MS2ReScoreConfigurationError
from ms2rescore._version import __version__

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources


def _parse_arguments() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="MS²ReScore: Sensitive PSM rescoring with predicted MS²\
            peak intensities."
    )
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument(
        "identification_file",
        type=str,
        help="path to identification file (pin, mzid, msms.txt, tandem xml...)",
    )
    parser.add_argument(
        "-m",
        metavar="FILE",
        action="store",
        type=str,
        dest="mgf_path",
        help="path to MGF file or directory with MGF files (default: derived from\
            identification file)",
    )
    parser.add_argument(
        "-c",
        metavar="FILE",
        action="store",
        type=str,
        dest="config_file",
        help="path to MS²ReScore configuration file (see README.md)",
    )
    parser.add_argument(
        "-t",
        metavar="PATH",
        action="store",
        type=str,
        dest="tmp_path",
        help="path to directory to place temporary files"
    )
    parser.add_argument(
        "-o",
        metavar="FILE",
        action="store",
        type=str,
        dest="output_filename",
        help="name for output files (default: derive from identification file)",
    )
    parser.add_argument(
        "-l",
        metavar="LEVEL",
        action="store",
        type=str,
        dest="log_level",
        default="info",
        help="logging level (default: `info`)",
    )
    return parser.parse_args()


def _validate_filenames(config: Dict) -> Dict:
    """Validate and infer input/output filenames."""
    # identification_file should exist
    id_file = config["general"]["identification_file"]
    if not os.path.isfile(id_file):
        raise FileNotFoundError(id_file)

    # MGF path should either be None, or existing path to file or dir
    mgf_path = config["general"]["mgf_path"]
    if mgf_path:
        if not os.path.exists(mgf_path):
            raise FileNotFoundError(mgf_path)

    # Output filename should be None or its path should exist. If not, make path.
    if config["general"]["output_filename"]:
        output_path = os.path.abspath(config["general"]["output_filename"])
        if not os.path.isdir(output_path):
            os.makedirs(output_path, exist_ok=True)
    else:
        config["general"]["output_filename"] = os.path.splitext(id_file)[0]

    return config


def _validate_num_cpu(config: Dict) -> Dict:
    """Validate requested num_cpu with available cpu count."""
    n_available = mp.cpu_count()
    if (config["general"]["num_cpu"] == -1) or (
        config["general"]["num_cpu"] > n_available
    ):
        config["general"]["num_cpu"] = n_available
    return config


def parse_config(parse_cli_args: bool = True, config_class: Optional[Dict] = None) -> Dict:
    """
    Parse and validate MS²ReScore configuration files and arguments.

    Default configuration, user configuration files, and CLI/class arguments are parsed
    in cascading order.

    Parameters
    ----------
    parse_cli_args : bool
        parse command line arguments or not, default True
    config_class : Dict
        dictionary with arguments from the Python class; required if `parse_cli_args`
        is False
    """
    config_schema = pkg_resources.open_text(package_data, "config_schema.json")
    config_default = pkg_resources.open_text(package_data, "config_default.json")

    # MS²ReScore can be run from the CLI, or as a Python module
    if parse_cli_args:
        args = _parse_arguments()
        config_user = args.config_file
        if config_class:
            raise MS2ReScoreConfigurationError(
                "If `parse_cli_args` is True, `config_class` must be None."
            )
    elif config_class:
        args = None
        config_user = config_class["config_file"]
    else:
        raise MS2ReScoreConfigurationError(
            "If `parse_cli_args` is False, `config_class` arguments are required."
        )

    cascade_conf = CascadeConfig(validation_schema=json.load(config_schema))
    cascade_conf.add_dict(json.load(config_default))
    if config_user:
        cascade_conf.add_json(config_user)
    if parse_cli_args:
        cascade_conf.add_namespace(args, subkey="general")
    elif config_class:
        cascade_conf.add_dict(config_class, subkey="general")
    config = cascade_conf.parse()

    config = _validate_filenames(config)
    config = _validate_num_cpu(config)

    return config
