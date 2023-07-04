"""Rescoring engines for MS²Rescore."""

import logging
from abc import ABC, abstractmethod

from psm_utils import PSMList

logger = logging.getLogger(__name__)


class RescoringEngineBase(ABC):
    """Base class for rescoring engines."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__()

    @abstractmethod
    def rescore(psm_list: PSMList):
        raise NotImplementedError()
