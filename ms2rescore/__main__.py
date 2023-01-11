"""MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs."""

import logging

from ms2rescore import MS2Rescore

logger = logging.getLogger(__name__)


def main():
    """Run MS²ReScore."""
    rescore = None
    try:
        rescore = MS2Rescore(parse_cli_args=True, configuration=None, set_logger=True)
        rescore.run()
    except Exception:
        logger.exception("Critical error occurred in MS²Rescore")
    finally:
        if rescore:
            rescore.save_log()


if __name__ == "__main__":
    main()
