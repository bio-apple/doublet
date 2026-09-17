# Do not reimplement Detector calls; packages are the contract

Status: superseded in part by [ADR 0013](0013-dual-r-python-packages.md).

The forbidden Expected Doublet Rate cutoff and the native-then-MAD Call rule are easy to violate with package defaults (`nExp`, `dbr`, top-x%). Agent-written one-off code will reintroduce those defaults. The contract is now the installable packages (`python -m doublet_rate`, `doubletRate::annotate_doublets`); `SKILL.md` tells the agent to run those. Root `scripts/` remain thin shims so older commands keep working, not a second implementation.
