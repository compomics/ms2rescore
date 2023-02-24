import logging
from logging import FileHandler
from pathlib import Path
from typing import Union

from rich.console import Console
from rich.logging import RichHandler

LOG_MAPPING = {
        "critical": logging.CRITICAL,
        "error": logging.ERROR,
        "warning": logging.WARNING,
        "info": logging.INFO,
        "debug": logging.DEBUG,
    }

def setup_logging(
    passed_level: str, log_file: Union[str, Path], rich_console: Console = None
):
    

    if passed_level not in LOG_MAPPING:
        print(
            "Invalid log level. Should be one of the following: ",
            ", ".join(LOG_MAPPING.keys()),
        )
        exit(1)

    if not rich_console:
        rich_console = Console(record=True)

    logging.basicConfig(
        format="%(name)s // %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=LOG_MAPPING[passed_level],
        handlers=[
            FileHandler(log_file, mode="w"),
            RichHandler(rich_tracebacks=True, console=rich_console, show_path=False),
        ],
    )
