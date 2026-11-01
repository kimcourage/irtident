#' Sample size needed to recover selected parameters
#'
#' A parametric-bootstrap style simulation: data sets of several sizes are
#' generated from a true model (by default the fitted model itself), the same
#' model is refitted, and the recovery of the selected parameters is
#' summarised. The result answers the practical question "how many
#' respondents would I need before this item's guessing parameter can be
#' estimated to within X?".
#'
#' For each parameter and sample size the function reports the share of
#' replications in which the model could be fitted (`success`), the share
#' that converged within the iteration limit (`conv`), and the bias and root
#' mean squared error (RMSE) over all fitted replications, plus the share of
#' estimates that end up at a bound. Non-converged fits are kept in the RMSE
#' because they are what an analyst would actually obtain; for 3PL/4PL
#' models EM often needs many cycles along a flat direction, so fits use
#' `technical = list(NCYCLES = 2000)` unless you pass your own `technical`.
#'
#' The required sample size is the smallest `N` in the grid from which on
#' the RMSE stays below `rmse_target` and the success rate stays above
#' `min_success`.
#'
#' @param mod A unidimensional single-group model fitted with [mirt::mirt()].
#' @param params Parameters to study; `NULL` (default) uses those flagged by
#'   [id_local()]. See [id_profile()] for the accepted formats.
#' @param N Sample sizes to simulate.
#' @param reps Replications per sample size.
#' @param rmse_target Named vector of target RMSEs by parameter type
#'   (`a1`, `d`, `g`, `u`, `COV_11`); graded intercepts `d1`, `d2`, ... use
#'   the `d` target. Values are on the reporting scale.
#' @param min_success Minimum share of replications in which the model must
#'   be fitted successfully (no error, no empty response category).
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
#'   per parameter and sample size), `required` (one row per parameter) and
#'   the raw `estimates`.
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
                          rmse_target = c(a1 = 0.3, d = 0.3, g = 0.05,
                                          u = 0.05, COV_11 = 0.2),
                          min_success = 0.9, truth = mod, local = NULL,
                          parallel = FALSE, ...) {
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
  extra <- list(...)
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
  summary$target <- .target_for(sub("^.*:", "", summary$param), rmse_target)
  summary$meets <- !is.na(summary$rmse) & summary$rmse <= summary$target &
                   summary$success >= min_success

  required <- do.call(rbind, lapply(params, function(p) {
    s <- summary[summary$param == p, ]
    s <- s[order(s$N), ]
    # smallest N from which on every larger N also meets the targets
    ok_from <- rev(cumprod(rev(s$meets))) == 1
    n_req <- if (any(ok_from)) s$N[which(ok_from)[1]] else NA_integer_
    data.frame(param = p, true = s$true[1], target_rmse = s$target[1],
               N_required = n_req,
               note = if (is.na(n_req)) paste0("> ", max(N), " (not reached)") else "",
               stringsAsFactors = FALSE)
  }))
  structure(list(summary = summary, required = required, estimates = estimates,
                 reps = reps, min_success = min_success),
            class = "irtident_samplesize")
}

# look up a per-parameter-type value; graded intercepts d1, d2, ... use "d"
.target_for <- function(name, values) {
  key <- ifelse(grepl("^d[0-9]+$", name), "d", name)
  out <- unname(values[key])
  out[is.na(out)] <- Inf
  out
}
