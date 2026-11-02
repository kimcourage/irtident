# irtident 0.1.0

Initial release.

* `id_local()`: item-level variance-inflation decomposition (conditional, within-item, between-item) from the information matrix, boundary flags and weakly identified directions from the scaled information matrix.
* `id_profile()`: profile likelihood intervals with closed / bound / open sides, traced by continuation. Bounds of `g` and `u` are judged with the 50:50 mixture of chi-squared(0) and chi-squared(1).
* `id_multistart()`: multi-start estimation and clustering of solutions into equivalent-fit clusters and inferior optima.
* `id_samplesize()`: repeated-estimate mean and SD, mirt SE, fitting success, convergence and SE failure across simulated sample sizes.
* `id_diagnose()`: combines local, profile and multi-start diagnostics by item.
* Supports unidimensional single-group models from `mirt::mirt()` with itemtypes `Rasch`, `2PL`, `3PL`, `3PLu`, `4PL` and `graded`.
