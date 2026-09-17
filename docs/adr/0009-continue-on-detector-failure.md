# One Detector failing does not abort the Sample

Join on barcode and keep going. An erroring Detector becomes Skipped with a reason; its Score and Call are empty. Aborting the whole run would throw away successful Detectors because Solo or DoubletFinder crashed, which is a common event in the mixed R/Python roster.