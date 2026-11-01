#' irtident: practical identification diagnostics for IRT models
#'
#' A fitted model can converge, pass mirt's second-order test and still have
#' parameters that the data barely determine. irtident tells you *which*
#' item parameters are weakly identified and *why*, for unidimensional,
#' single-group models fitted with [mirt::mirt()] (itemtypes `Rasch`, `2PL`,
#' `3PL`, `3PLu`, `4PL` and `graded`).
#'
#' Main functions:
#' * [id_diagnose()]: runs the diagnostics below and gives one verdict per
#'   item, with advice;
#' * [id_local()]: item-level variance-inflation decomposition from the
#'   information matrix;
#' * [id_profile()]: profile likelihood intervals and flatness;
#' * [id_multistart()]: multi-start estimation with clustering of solutions;
#' * [id_samplesize()]: simulation of the sample size needed to recover a
#'   parameter.
#'
#' @keywords internal
"_PACKAGE"
