"""Feature generation for MS²Rescore."""

import logging
from abc import ABC, abstractmethod

from psm_utils import PSMList

logger = logging.getLogger(__name__)


class FeatureGenerator(ABC):
    """Base class for feature generators."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__()

    @abstractmethod
    def add_features(psm_list: PSMList):
        raise NotImplementedError()

