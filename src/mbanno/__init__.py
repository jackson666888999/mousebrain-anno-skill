"""
mbanno — Mouse Brain Consensus Annotator

A reproducible consensus annotation framework for mouse whole-brain
single-cell and spatial transcriptomics using Allen/BICCN reference taxonomies.
"""

__version__ = "0.1.0"
__author__ = "huangyanling"

from .qc import QCProcessor
from .references import ReferenceManager
from .consensus import ConsensusAnnotator
from .confidence import ConfidenceEvaluator

__all__ = [
    "QCProcessor",
    "ReferenceManager",
    "ConsensusAnnotator",
    "ConfidenceEvaluator",
    "__version__",
]
