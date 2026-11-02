# irtident 0.1.0

* `id_samplesize()` now reports recovery without an RMSE cutoff. Removed
  `rmse_target`, `min_success`, the `required` table, and the `target`/`meets`
  columns. Printing no longer recommends a sample size, and plots no longer
  include target lines. Calls using the removed arguments receive an error.

Initial release.

* `id_local()`: item-level variance-inflation decomposition (intrinsic, within-item, between-item) from the information matrix, boundary flags and weakly identified directions.
* `id_profile()`: profile likelihood intervals with closed / bound / open sides and a flatness status, traced by continuation.
* `id_multistart()`: multi-start estimation, clustering of solutions, separation of equivalent-fit modes from inferior local optima.
* `id_samplesize()`: bias, RMSE and fitting success across simulated sample sizes.
* `id_diagnose()`: combines local, profile and multi-start diagnostics by item.
* Supports unidimensional single-group models from `mirt::mirt()` with itemtypes `Rasch`, `2PL`, `3PL`, `3PLu`, `4PL` and `graded`.
