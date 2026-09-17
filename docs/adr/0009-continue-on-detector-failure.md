# One Detector failing does not abort the Sample

Join on barcode and keep going. An erroring Detector is recorded as `failed` with a reason; expected absences stay `skipped` (size gate, missing package) or `skipped:no_gpu` (Solo without CUDA or Apple MPS). Score and Call stay empty. Aborting the whole run would throw away successful Detectors because Solo or DoubletFinder crashed, which is a common event in the mixed R/Python roster.
