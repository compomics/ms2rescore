"""Graphical user interface using Gooey."""

import argparse
import logging

from gooey import Gooey, GooeyParser, local_resource_path

from ms2rescore import MS2ReScore, package_data
from ms2rescore._exceptions import MS2ReScoreConfigurationError
from ms2rescore._version import __version__

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

logger = logging.getLogger(__name__)

with pkg_resources.path(package_data, "img") as img_dir:
    img_dir = local_resource_path(img_dir)


@Gooey(
    program_name="MS²Rescore",
    program_description="Sensitive PSM rescoring with predicted MS² peak intensities.",
    image_dir=img_dir,
)
def main():
    """Run MS²Rescore."""
    conf = _parse_arguments().__dict__
    rescore = MS2ReScore(parse_cli_args=False, configuration=conf, set_logger=True)
    rescore.run()


def _parse_arguments() -> argparse.Namespace:
    """Parse GUI arguments."""
    parser = GooeyParser()
    parser.add_argument(
        "identification_file",
        metavar="Identification file",
        type=str,
        help="Path to identification file (pin, mzid, msms.txt, tandem xml...)",
        widget="FileChooser",
    )
    parser.add_argument(
        "-m",
        metavar="Spectrum file directory",
        action="store",
        type=str,
        dest="mgf_path",
        help=(
            "Path to MGF file or directory with MGF files (default: derived from "
            "identification file)"
        ),
        widget="DirChooser",
    )
    parser.add_argument(
        "-c",
        metavar="Configuration file",
        action="store",
        type=str,
        dest="config_file",
        help="Path to MS²Rescore configuration file (see README.md)",
        widget="FileChooser",
    )
    parser.add_argument(
        "-t",
        metavar="Temporary file directory",
        action="store",
        type=str,
        dest="tmp_path",
        help="Path to directory to place temporary files",
        widget="DirChooser",
    )
    parser.add_argument(
        "-o",
        metavar="Output filename prefix",
        action="store",
        type=str,
        dest="output_filename",
        help="Name for output files (default: derive from identification file)",
        widget="FileSaver",
    )
    parser.add_argument(
        "-l",
        metavar="Logging level",
        action="store",
        type=str,
        dest="log_level",
        default="info",
        # help="Logging level",
        widget="Dropdown",
        choices=["debug", "info", "warning", "error", "critical"],
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
