"""MS²Rescore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs."""

import argparse
import logging
import sys
from pathlib import Path
from typing import Union

from rich.console import Console
from rich.logging import RichHandler
from rich.text import Text

from ms2rescore import __version__
from ms2rescore.config_parser import parse_configurations
from ms2rescore.core import rescore
from ms2rescore.exceptions import MS2RescoreConfigurationError

try:
    import matplotlib.pyplot as plt

    plt.set_loglevel("warning")
except ImportError:
    pass

LOG_MAPPING = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}
LOGGER = logging.getLogger(__name__)
CONSOLE = Console(record=True)


def _print_credits():
    """Print software credits to terminal."""
    text = Text()
    text.append("\n")
    text.append("MS²Rescore", style="bold link https://github.com/compomics/ms2rescore")
    text.append(f" (v{__version__})\n", style="bold")
    text.append("Developed at CompOmics, VIB / Ghent University, Belgium.\n")
    text.append("Please cite: ")
    text.append(
        "Declercq et al. MCP (2022)", style="link https://doi.org/10.1016/j.mcpro.2022.100266"
    )
    text.append("\n")
    text.stylize("cyan")
    CONSOLE.print(text)


def _argument_parser() -> argparse.ArgumentParser:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="MS²Rescore: Sensitive PSM rescoring with predicted features.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=42),
    )
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument(
        "-p",
        "--psm-file",
        metavar="FILE",
        action="store",
        type=str,
        nargs="*",
        dest="psm_file",
        help="path to PSM file (PIN, mzIdentML, MaxQuant msms, X!Tandem XML...)",
    )
    parser.add_argument(
        "-t",
        "--psm-file-type",
        metavar="STR",
        action="store",
        type=str,
        dest="psm_file_type",
        help="PSM file type (default: 'infer')",
    )
    parser.add_argument(
        "-s",
        "--spectrum-path",
        metavar="FILE/DIR",
        action="store",
        type=str,
        dest="spectrum_path",
        help="path to MGF/mzML spectrum file or directory with spectrum files (default: derived\
            from identification file)",
    )
    parser.add_argument(
        "-c",
        "--config-file",
        metavar="FILE",
        action="store",
        type=str,
        dest="config_file",
        help="path to MS²Rescore configuration file (see README.md)",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        metavar="FILE",
        action="store",
        type=str,
        dest="output_path",
        help="Path and stem for output file names (default: derive from identification file)",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        metavar="STR",
        action="store",
        type=str,
        dest="log_level",
        help="logging level (default: `info`)",
    )
    parser.add_argument(
        "-n",
        "--processes",
        metavar="INT",
        action="store",
        type=int,
        dest="processes",
        help="number of parallel processes available to MS²Rescore",
    )
    parser.add_argument(
        "-f",
        "--fasta-file",
        metavar="FILE",
        action="store",
        type=str,
        dest="fasta_file",
        help="path to FASTA file",
    )

    return parser


def _setup_logging(passed_level: str, log_file: Union[str, Path]):
    """Setup logging for writing to log file and Rich Console."""
    if passed_level not in LOG_MAPPING:
        raise MS2RescoreConfigurationError(
            f"Invalid log level '{passed_level}'. "
            f"Valid levels are: {', '.join(LOG_MAPPING.keys())}"
        )
    logging.basicConfig(
        format="%(name)s // %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=LOG_MAPPING[passed_level],
        handlers=[
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
            RichHandler(rich_tracebacks=True, console=CONSOLE, show_path=False),
        ],
    )


def main():
    """Run MS²Rescore command-line interface."""
    _print_credits()

    # Parse CLI arguments and configuration file
    parser = _argument_parser()
    cli_args = parser.parse_args()
    try:
        if cli_args.config_file:
            config = parse_configurations([cli_args.config_file, cli_args])
        else:
            config = parse_configurations(cli_args)
    except MS2RescoreConfigurationError as e:
        LOGGER.critical(e)
        sys.exit(1)

    # Setup logging
    _setup_logging(
        config["ms2rescore"]["log_level"], config["ms2rescore"]["output_path"] + ".log.txt"
    )

    # Run MS²Rescore
    try:
        rescore(configuration=config)
    except Exception as e:
        LOGGER.exception(e)
        sys.exit(1)
    finally:
        CONSOLE.save_html(config["ms2rescore"]["output_path"] + ".log.html")


if __name__ == "__main__":
    main()
