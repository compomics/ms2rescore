"""Parse configuration from command line arguments and configuration files."""

import argparse
import importlib.resources as pkg_resources
import json
import multiprocessing as mp
import tempfile
from pathlib import Path
from typing import Dict, Optional

import tomlkit
from cascade_config import CascadeConfig

from ms2rescore import __version__, package_data
from ms2rescore.exceptions import MS2RescoreConfigurationError


def _parse_arguments() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="MS²Rescore: Sensitive PSM rescoring with predicted MS²\
            peak intensities."
    )
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument(
        "-p",
        metavar="FILE",
        action="store",
        type=str,
        dest="psm_file",
        help="path to PSM file (pin, mzid, msms.txt, tandem xml...)",
    )
    parser.add_argument(
        "-m",
        metavar="FILE",
        action="store",
        type=str,
        dest="spectrum_path",
        help="path to MGF file or directory with MGF files (default: derived from\
            identification file)",
    )
    parser.add_argument(
        "-c",
        metavar="FILE",
        action="store",
        type=str,
        dest="config_file",
        help="path to MS²Rescore configuration file (see README.md)",
    )
    parser.add_argument(
        "-t",
        metavar="PATH",
        action="store",
        type=str,
        dest="tmp_path",
        help="path to directory to place temporary files",
    )
    parser.add_argument(
        "-o",
        metavar="FILE",
        action="store",
        type=str,
        dest="output_path",
        help="name for output files (default: derive from identification file)",
    )
    parser.add_argument(
        "-l",
        metavar="LEVEL",
        action="store",
        type=str,
        dest="log_level",
        help="logging level (default: `info`)",
    )
    parser.add_argument(
        "-n",
        metavar="VALUE",
        action="store",
        type=int,
        dest="processes",
        default=None,
        help="number of parallel processes available to MS²Rescore",
    )
    parser.add_argument(
        "--psm_file_type",
        metavar="FILE",
        action="store",
        type=str,
        dest="psm_file_type",
        default=None,
        help="determines psm parser to use from PSM_utils (default: 'infer')",
    )

    return parser.parse_args()


def _validate_filenames(config: Dict) -> Dict:
    """Validate and infer input/output filenames."""
    # psm_file should be provided
    if not config["ms2rescore"]["psm_file"]:
        raise MS2RescoreConfigurationError("`psm_file` should be provided.")

    # psm_file should exist
    id_file = Path(config["ms2rescore"]["psm_file"])
    if not id_file.is_file():
        raise FileNotFoundError(id_file)
    config["ms2rescore"]["psm_file"] = str(id_file)

    # spectrum_path should either be None, or existing path to file or dir
    if config["ms2rescore"]["spectrum_path"]:
        spectrum_path = Path(config["ms2rescore"]["spectrum_path"])
        if not spectrum_path.exists():
            raise FileNotFoundError(spectrum_path)
        config["ms2rescore"]["spectrum_path"] = str(spectrum_path)

    # Output filename should be None or its path should exist. If not, make path
    if config["ms2rescore"]["output_path"]:
        output_path = Path(config["ms2rescore"]["output_path"])
        if not output_path.is_dir():
            output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = Path(id_file).parent
    config["ms2rescore"]["output_path"] = str(output_path)

    # tmp_path should either be None, or path to dir. If not, make path.
    if config["ms2rescore"]["tmp_path"]:
        tmp_path = Path(config["ms2rescore"]["tmp_path"])
        tmp_path.mkdir(parents=True, exist_ok=True)
    else:
        tmp_path = tempfile.mkdtemp()
    config["ms2rescore"]["tmp_path"] = str(tmp_path)

    return config


def _validate_processes(config: Dict) -> Dict:
    """Validate requested processes with available cpu count."""
    n_available = mp.cpu_count()
    if (config["ms2rescore"]["processes"] == -1) or (
        config["ms2rescore"]["processes"] > n_available
    ):
        config["ms2rescore"]["processes"] = n_available
    return config


def parse_config(parse_cli_args: bool = True, config_class: Optional[Dict] = None) -> Dict:
    """
    Parse and validate MS²Rescore configuration files and arguments.

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

    # MS²Rescore can be run from the CLI, or as a Python module
    if parse_cli_args:
        args = _parse_arguments()
        config_user = args.config_file
        if config_class:
            raise MS2RescoreConfigurationError(
                "If `parse_cli_args` is True, `config_class` must be None."
            )
    elif config_class:
        args = None
        config_user = config_class["ms2rescore"]["config_file"]
    else:
        raise MS2RescoreConfigurationError(
            "If `parse_cli_args` is False, `config_class` arguments are required."
        )

    cascade_conf = CascadeConfig(
        validation_schema=json.load(config_schema),
        none_overrides_value=False,
        max_recursion_depth=1,
    )
    cascade_conf.add_dict(json.load(config_default))
    if config_user:
        if Path(config_user).suffix.lower() == ".json":
            cascade_conf.add_json(config_user)
        elif Path(config_user).suffix.lower() == ".toml":
            with open(config_user, "rt") as toml_file:
                cascade_conf.add_dict(tomlkit.load(toml_file))
        else:
            raise MS2RescoreConfigurationError(
                "Unknown file extension for configuration file. Should be `json` or " "`toml`."
            )
    if parse_cli_args:
        cascade_conf.add_namespace(args, subkey="ms2rescore")
    elif config_class:
        cascade_conf.add_dict(config_class)
    config = cascade_conf.parse()

    config = _validate_filenames(config)
    config = _validate_processes(config)

    config["feature_generators"] = {k.lower(): v for k, v in config["feature_generators"].items()}
    config["rescoring_engine"] = {k.lower(): v for k, v in config["rescoring_engine"]}


    return config
