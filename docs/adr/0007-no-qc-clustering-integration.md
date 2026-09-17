# This skill does not QC, cluster, or integrate

The input Sample is already cell-called raw counts. QC, clustering, annotation, and integration are outside this context. Doing them here would turn rate measurement into a preprocessing pipeline, and cluster labels would leak into Calls that are supposed to come from Scores.
