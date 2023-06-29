"""MS²Rescore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs."""

import logging

from rich.console import Console
from rich.text import Text

from ms2rescore import MS2Rescore, __version__

try:
    import matplotlib.pyplot as plt

    plt.set_loglevel("warning")
except ImportError:
    pass

logger = logging.getLogger(__name__)
console = Console(record=True)


def _build_credits():
    """Build credits."""
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
    return text


def main():
    """Run MS²Rescore."""
    console.print(_build_credits())

    rescore = None
    try:
        rescore = MS2Rescore(
            parse_cli_args=True, configuration=None, set_logger=True, rich_console=console
        )
        rescore.run()
    except Exception as e:
        logger.exception(e)
        raise e
    finally:
        if rescore:
            rescore.save_log()


if __name__ == "__main__":
    main()
