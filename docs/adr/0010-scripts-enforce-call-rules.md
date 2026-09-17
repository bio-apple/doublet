# Ship scripts; the skill executes them rather than regenerating Detector calls

The forbidden Expected Doublet Rate cutoff and the native-then-MAD Call rule are easy to violate with package defaults (`nExp`, `dbr`, top-x%). Agent-written one-off code will reintroduce those defaults. Root `scripts/` is the contract; `SKILL.md` tells the agent to run those scripts.