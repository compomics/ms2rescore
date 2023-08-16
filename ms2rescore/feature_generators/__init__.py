"""Feature generation for MS²Rescore."""

from ms2rescore.feature_generators.basic import BasicFeatureGenerator
from ms2rescore.feature_generators.deeplc import DeepLCFeatureGenerator
from ms2rescore.feature_generators.maxquant import MaxQuantFeatureGenerator
from ms2rescore.feature_generators.ms2pip import MS2PIPFeatureGenerator

FEATURE_GENERATORS = {
    "basic": BasicFeatureGenerator,
    "ms2pip": MS2PIPFeatureGenerator,
    "deeplc": DeepLCFeatureGenerator,
    "maxquant": MaxQuantFeatureGenerator,
}
