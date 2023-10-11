"""Entrypoint for MS²Rescore GUI."""

import multiprocessing
import os
import contextlib

from ms2rescore.gui.app import app


def main():
    """Entrypoint for MS²Rescore GUI."""
    multiprocessing.freeze_support()
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        app()


if __name__ == "__main__":
    main()
