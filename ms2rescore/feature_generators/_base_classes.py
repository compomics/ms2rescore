from abc import ABC, abstractmethod

from psm_utils import PSMList


class FeatureGeneratorBase(ABC):
    """Base class for feature generators."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__()

    @property
    @abstractmethod
    def feature_names(self):
        pass

    @abstractmethod
    def add_features(psm_list: PSMList):
        pass
