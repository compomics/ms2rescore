"""MS²ReScore exceptions."""


class MS2ReScoreError(Exception):
    """Generic MS2ReScore error."""

    pass


class MS2ReScoreConfigurationError(MS2ReScoreError):
    """Invalid MS2ReScore configuration."""

    pass
