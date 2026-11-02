#' irtident: practical identification diagnostics for IRT models
#'
#' Examines parameter uncertainty using variance decomposition, profile
#' likelihoods, multiple starting values and sample-size simulations.
#' Supports unidimensional, single-group models fitted with [mirt::mirt()]
#' (itemtypes `Rasch`, `2PL`, `3PL`, `3PLu`, `4PL` and `graded`).
#'
#' Main functions:
#' * [id_diagnose()]: combines local, profile and multi-start results by item;
#' * [id_local()]: item-level variance-inflation decomposition from the
#'   information matrix;
#' * [id_profile()]: profile likelihood intervals and flatness;
#' * [id_multistart()]: multi-start estimation with clustering of solutions;
#' * [id_samplesize()]: parameter precision across simulated sample sizes.
#'
#' @keywords internal
"_PACKAGE"
