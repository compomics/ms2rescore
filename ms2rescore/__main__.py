"""MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities and retention times."""

from ms2rescore import MS2ReScore


def main():
    """Run MS²ReScore."""
    rescore = MS2ReScore(parse_cli_args=True, configuration=None)
    rescore.run()


if __name__ == "__main__":
    main()
