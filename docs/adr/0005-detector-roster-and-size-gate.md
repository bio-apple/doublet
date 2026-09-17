# Detector roster: score-capable methods, size gate, no DoubletDecon

Run scDblFinder, Scrublet, cxds, bcds, and hybrid on every Sample. Also run DoubletDetection, DoubletFinder (pANN score only), and Solo when the Sample has at most 20,000 cells. Drop DoubletDecon because it does not emit a continuous Score. Above 20,000 cells, skip the slow Detectors and record them as Skipped rather than pretending they ran. DoubletFinder may not consume `nExp`; scDblFinder and Solo may not consume an Expected Doublet Rate to place the Call.
