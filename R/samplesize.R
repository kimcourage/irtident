#' Parameter recovery across sample sizes
#'
#' Generates data sets at several sample sizes from an assumed true model
#' (by default the fitted model), refits the model and summarises parameter
#' recovery. Results are conditional on that model and the
#' fitting settings.
#'
#' For each parameter and sample size the function reports the share of
#' replications in which the model could be fitted (`success`), the share
#' that converged within the iteration limit (`conv`), and the bias and root
#' mean squared error (RMSE) over all fitted replications, plus the share of
#' estimates that end up at a bound. Non-converged fits that return estimates
#' are included in the RMSE; failed fits are excluded. Fits use
#' `technical = list(NCYCLES = 2000)` unless you pass your own `technical`.
#'
#' Results are descriptive: the function does not apply an RMSE cutoff or
#' select a required sample size. RMSE is an average loss, not a bound on an
#' individual estimate's error. Interpret it alongside the fitting success
#' and convergence rates and allow for simulation variability.
#'
#' @param mod A unidimensional single-group model fitted with [mirt::mirt()].
#' @param params Parameters to study; `NULL` (default) uses those flagged by
#'   [id_local()]. See [id_profile()] for the accepted formats.
#' @param N Sample sizes to simulate.
#' @param reps Replications per sample size.
#' @param truth A fitted mirt model with the same items and itemtypes from
#'   which data are generated. Defaults to `mod`. To ask "what if the true
#'   guessing parameter were 0.2?", edit `mirt::mod2values(mod)` and refit
#'   with `pars = ..., TOL = NaN`.
#' @param local Optional `irtident_local` object (to reuse flagged
#'   parameters).
#' @param parallel If `TRUE`, replications are run in parallel with
#'   `future.apply`.
#' @param ... Further arguments passed to [mirt::mirt()] for every fit.
#'
#' @return An object of class `irtident_samplesize` with `summary` (one row
#'   per parameter and sample size), the raw `estimates` and `reps`.
#' @examples
#' \donttest{
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' fit <- mirt(dat, 1, itemtype = "3PL", verbose = FALSE)
#' set.seed(1)
#' ss <- id_samplesize(fit, params = "Item.2:g", N = c(500, 2000), reps = 10)
#' ss
#' }
#' @export
id_samplesize <- function(mod, params = NULL,
                          N = c(250, 500, 1000, 2000, 4000), reps = 30,
                          truth = mod, local = NULL,
                          parallel = FALSE, ...) {
  extra <- list(...)
  removed <- intersect(names(extra), c("rmse_target", "min_success"))
  if (length(removed)) {
    stop("Removed argument(s): ", paste(removed, collapse = ", "),
         ". id_samplesize() reports recovery without a cutoff or required N.",
         call. = FALSE)
  }
  .check_mod(mod)
  .check_mod(truth)
  if (!identical(mirt::extract.mirt(mod, "itemtype"),
                 mirt::extract.mirt(truth, "itemtype"))) {
    stop("`truth` must have the same items and itemtypes as `mod`.", call. = FALSE)
  }
  pt <- .par_table(mod)
  if (is.null(params)) {
    if (is.null(local)) local <- id_local(mod)
    params <- local$params$param[local$params$status %in% c("weak", "boundary")]
  }
  params <- .resolve_params(params, pt, default = character())
  if (!length(params)) {
    stop("No parameters to study: none were flagged by id_local(). ",
         "Pass `params` explicitly.", call. = FALSE)
  }
  N <- sort(unique(as.integer(N)))
  tt <- .par_table(truth)
  true_val <- tt$value[match(params, tt$param)]
  itemtype <- mirt::extract.mirt(mod, "itemtype")
  K <- mirt::extract.mirt(mod, "K")

  # simulate all data up front so results do not depend on `parallel`
  jobs <- expand.grid(rep = seq_len(reps), N = N)
  datasets <- lapply(jobs$N, function(n) mirt::simdata(model = truth, N = n))
  if (is.null(extra$technical)) extra$technical <- list(NCYCLES = 2000)

  res <- .lapply(seq_along(datasets), function(j) {
    dat <- datasets[[j]]
    args <- utils::modifyList(list(data = dat, model = 1, itemtype = itemtype,
                                   verbose = FALSE), extra)
    fit <- tryCatch(suppressWarnings(suppressMessages(do.call(mirt::mirt, args))),
                    error = function(e) NULL)
    na <- rep(NA_real_, length(params))
    if (is.null(fit) || !identical(as.integer(mirt::extract.mirt(fit, "K")),
                                   as.integer(K))) {
      return(list(ok = FALSE, conv = FALSE, values = na))
    }
    fp <- .par_table(fit)
    vals <- .align_reflection(fp$value, fp$name)
    list(ok = TRUE, conv = isTRUE(mirt::extract.mirt(fit, "converged")),
         values = vals[match(params, fp$param)])
  }, parallel = parallel)

  est <- t(vapply(res, `[[`, numeric(length(params)), "values"))
  if (length(params) == 1L) est <- matrix(est, ncol = 1L)
  colnames(est) <- params
  ok <- vapply(res, `[[`, logical(1), "ok")
  conv <- vapply(res, `[[`, logical(1), "conv")
  estimates <- cbind(jobs, ok = ok, converged = conv, as.data.frame(est,
                                                       check.names = FALSE))

  rows <- list()
  for (k in seq_along(params)) {
    p <- params[k]
    row <- pt[pt$param == p, ]
    for (n in N) {
      sel <- jobs$N == n
      use <- sel & ok & !is.na(est[, k])
      e <- est[use, k]
      at_b <- (is.finite(row$lbound) & e - row$lbound < 0.01) |
              (is.finite(row$ubound) & row$ubound - e < 0.01)
      rows[[length(rows) + 1L]] <- data.frame(
        param = p, N = n, true = true_val[k], success = mean(ok[sel]),
        conv = mean(conv[sel]),
        mean = if (length(e)) mean(e) else NA_real_,
        bias = if (length(e)) mean(e) - true_val[k] else NA_real_,
        rmse = if (length(e)) sqrt(mean((e - true_val[k])^2)) else NA_real_,
        p_bound = if (length(e)) mean(at_b) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  summary <- do.call(rbind, rows)
  structure(list(summary = summary, estimates = estimates, reps = reps),
            class = "irtident_samplesize")
}
