"""MS²Rescore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs."""

__version__ = "3.0.0-dev4"

from warnings import filterwarnings

# mzmlb is not used, so hdf5plugin is not needed
filterwarnings(
    "ignore",
    message="hdf5plugin is missing",
    category=UserWarning,
    module="psims.mzmlb",
)

from ms2rescore.core import rescore
