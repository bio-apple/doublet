# Report scores, calls, and predicted rate; do not remove cells

The skill measures doublets, it does not clean the dataset. For one RNA input it must emit a Doublet Score per cell, a data-driven Doublet Call, and a Predicted Doublet Rate derived from those calls. Removal is out of scope because a rate analysis that silently drops cells mixes measurement with QC, and OSCA warns that deleting high-scoring cells often leaves residual doublet structure.
