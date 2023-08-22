import logging

import click
from rich.logging import RichHandler

from ms2rescore.report.generate import generate_report

logger = logging.getLogger(__name__)


@click.command()
@click.argument("output_prefix", type=str)
def main(output_prefix: str):
    logging.getLogger("mokapot").setLevel(logging.WARNING)
    logging.basicConfig(
        level=logging.INFO,
        handlers=[RichHandler(rich_tracebacks=True)],
        format="%(message)s",
    )

    output_prefix = "examples/id/msms"

    try:
        generate_report(output_prefix)
    except Exception as e:
        logger.exception(e)
        exit(1)


if __name__ == "__main__":
    main()
