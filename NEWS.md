# irtident 0.1.0

First version (MVP).

* `id_local()`: item-level variance-inflation decomposition (intrinsic, within-item, between-item) from the information matrix, boundary flags and weakly identified directions.
* `id_profile()`: profile likelihood intervals with closed / bound / open sides and a flatness status, traced by continuation.
* `id_multistart()`: multi-start estimation, clustering of solutions, separation of equivalent-fit modes from inferior local optima.
* `id_samplesize()`: simulation of the sample size needed to recover selected parameters.
* `id_diagnose()`: runs the above and gives one verdict and advice per item.
* Supports unidimensional single-group models from `mirt::mirt()` with itemtypes `Rasch`, `2PL`, `3PL`, `3PLu`, `4PL` and `graded`.
