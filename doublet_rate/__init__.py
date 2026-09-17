"""RNA-only doublet Score, Call, and Predicted Doublet Rate."""

from .annotate import annotate_anndata, detect_doublets
from .io_counts import InputError, SIZE_GATE

__all__ = [
    "detect_doublets",
    "annotate_anndata",
    "InputError",
    "SIZE_GATE",
]
