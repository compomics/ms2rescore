"""MS²Rescore exceptions."""


class MS2RescoreError(Exception):
    """Generic MS2Rescore error."""

    pass


class MS2RescoreConfigurationError(MS2RescoreError):
    """Invalid MS2Rescore configuration."""

    pass
