"""MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs."""

import logging

from ms2rescore import MS2ReScore

logger = logging.getLogger(__name__)


def main():
    """Run MS²ReScore."""
    rescore = None
    try:
        rescore = MS2ReScore(parse_cli_args=True, configuration=None, set_logger=True)
        rescore.run()
    except Exception:
        logger.exception("Critical error occured in MS²Rescore")
    finally:
        if rescore:
            rescore.save_log()


if __name__ == "__main__":
    main()
