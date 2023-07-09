"""Entrypoint for MS²Rescore GUI."""

import multiprocessing

from ms2rescore.gui.app import app

if __name__ == "__main__":
    multiprocessing.freeze_support()
    app()
