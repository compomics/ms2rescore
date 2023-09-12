"""Entrypoint for MS²Rescore GUI."""

import multiprocessing

from ms2rescore.gui.app import app


def main():
    """Entrypoint for MS²Rescore GUI."""
    multiprocessing.freeze_support()
    app()


if __name__ == "__main__":
    main()
    main()
