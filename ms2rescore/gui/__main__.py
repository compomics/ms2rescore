"""Entrypoint for MS²Rescore GUI."""

from __future__ import annotations

import contextlib
import multiprocessing
import os

from ms2rescore.gui.app import app


def main():
    """Entrypoint for MS²Rescore GUI."""
    multiprocessing.freeze_support()
    # Redirect stdout/stderr when running GUI (packaged app might not have console attached)
    if os.devnull is not None:
        with (
            contextlib.redirect_stdout(open(os.devnull, "w", encoding="utf-8")),
            contextlib.redirect_stderr(open(os.devnull, "w", encoding="utf-8")),
        ):
            app()
    else:
        app()


if __name__ == "__main__":
    main()
