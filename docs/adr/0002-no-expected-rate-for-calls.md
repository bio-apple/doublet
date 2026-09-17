# Do not use Expected Doublet Rate to make Calls

Calls must come from a data-driven threshold on each method's Score. Feeding `nExp`, `dbr`, 10x 0.8%/1000 cells, or Solo's expected-number calibration as a cutoff would make Predicted Doublet Rate echo the loading formula instead of measuring the data. Methods that can only emit a binary label by consuming that formula are out of scope for calling.
