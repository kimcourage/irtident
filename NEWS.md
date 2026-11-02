# irtident 0.1.0

* `id_samplesize()` replaces bias and RMSE with the current estimate,
  repeated-estimate mean and SD, and mirt SE (`sqrt(mean(SE^2))`) across N.
  Summaries use the same converged fits with usable SEs on the reporting
  scale. The output includes the number used, convergence and SE failure;
  plots compare repeated-estimate SD with mirt SE. Raw estimates and SEs
  remain available. The `truth` argument specifies the generating scenario.
* README examples now include actual diagnostic and sample-size output.

* `id_samplesize()` now reports recovery without an RMSE cutoff. Removed
  `rmse_target`, `min_success`, the `required` table, and the `target`/`meets`
  columns. Printing no longer recommends a sample size, and plots no longer
  include target lines. Calls using the removed arguments receive an error.

Initial release.

* `id_local()`: item-level variance-inflation decomposition (intrinsic, within-item, between-item) from the information matrix, boundary flags and weakly identified directions.
* `id_profile()`: profile likelihood intervals with closed / bound / open sides and a flatness status, traced by continuation.
* `id_multistart()`: multi-start estimation, clustering of solutions, separation of equivalent-fit modes from inferior local optima.
* `id_samplesize()`: parameter precision and fitting success across simulated sample sizes.
* `id_diagnose()`: combines local, profile and multi-start diagnostics by item.
* Supports unidimensional single-group models from `mirt::mirt()` with itemtypes `Rasch`, `2PL`, `3PL`, `3PLu`, `4PL` and `graded`.
