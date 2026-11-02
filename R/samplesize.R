#' Parameter precision across sample sizes
#'
#' Generates data sets at several sample sizes from a fitted model, refits
#' each data set and compares the variability of estimates with mirt's
#' reported standard errors. Results describe the generating model and
#' fitting settings used in the simulation.
#'
#' For each parameter and N, `current` is the estimate in `mod`. `mean` and
#' `sd` summarise the repeated estimates; `mirt_se` is
#' `sqrt(mean(SE^2))` across the same replications. These three summaries use
#' only converged fits with a finite estimate and a finite, positive SE for
#' that parameter. `n_used` reports their number. The SD is `NA` with fewer
#' than two usable replications.
#'
#' `success` and `conv` are proportions of all attempted replications.
#' `se_failure` is the proportion of converged, successful fits without a
#' usable estimate/SE pair, or `NA` if none converged. `p_bound` is the share
#' of usable estimates within 0.01 of a parameter bound. Inspect these rates
#' alongside the precision summaries: exclusions can affect the comparison.
#'
#' Estimates and SEs use the reporting scale. For `g` and `u`, mirt's logit
#' SE is converted by the delta method: `p * (1 - p) * SE_logit`. Exact
#' probability boundaries have no usable delta-method SE; near a boundary
#' the approximation may be unreliable. No bias, RMSE, cutoff or required
#' sample size is calculated.
#'
#' Fits use `technical = list(NCYCLES = 2000)` unless supplied in `...`.
#' SEs are evaluated at each converged fit without further optimisation,
#' using the same fitting arguments and mirt's default SE method unless
#' `SE.type` is supplied. Failed SE calculations do not discard the fit.
#'
#' @param mod A unidimensional single-group model fitted with [mirt::mirt()].
#' @param params Parameters to study; `NULL` uses those flagged by
#'   [id_local()]. See [id_profile()] for the accepted formats.
#' @param N Sample sizes to simulate.
#' @param reps Replications per sample size.
#' @param truth Generating model with the same items and itemtypes as `mod`.
#'   Defaults to `mod`. The argument name is retained for compatibility;
#'   these values specify a simulation scenario, not known population values.
#' @param local Optional `irtident_local` object to reuse flagged parameters.
#' @param parallel If `TRUE`, run replications with `future.apply`.
#' @param ... Further arguments to [mirt::mirt()], including `SE.type`.
#'   `SE` is managed by this function and must not be supplied here.
#'
#' @return An `irtident_samplesize` object with `summary`, raw `estimates`,
#'   `standard_errors` on the reporting scale, `reps`, `current_N` and
#'   `generating` (the parameter values used to generate data).
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
         ". id_samplesize() reports precision without a cutoff or required N.",
         call. = FALSE)
  }
  if ("SE" %in% names(extra)) {
    stop("`SE` is managed by id_samplesize(); use `SE.type` to select its method.",
         call. = FALSE)
  }
  .check_mod(mod)
  .check_mod(truth)
  if (!identical(mirt::extract.mirt(mod, "itemtype"),
                 mirt::extract.mirt(truth, "itemtype")) ||
      !identical(colnames(mirt::extract.mirt(mod, "data")),
                 colnames(mirt::extract.mirt(truth, "data")))) {
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
  # Include fixed generating parameters when a supplied scenario fixes them.
  tt <- mirt::mod2values(truth)
  generating <- data.frame(param = params, value = tt$value[
    match(params, paste0(tt$item, ":", tt$name))])
  itemtype <- mirt::extract.mirt(mod, "itemtype")
  K <- mirt::extract.mirt(mod, "K")

  # Generate data before dispatch so sequential and parallel inputs agree.
  jobs <- expand.grid(rep = seq_len(reps), N = N)
  datasets <- lapply(jobs$N, function(n) mirt::simdata(model = truth, N = n))
  if (is.null(extra$technical)) extra$technical <- list(NCYCLES = 2000)

  res <- .lapply(seq_along(datasets), function(j) {
    args <- utils::modifyList(list(data = datasets[[j]], model = 1,
                                   itemtype = itemtype, verbose = FALSE), extra)
    args$SE <- FALSE
    fit <- .samplesize_fit(args)
    na <- rep(NA_real_, length(params))
    if (is.null(fit) || !identical(as.integer(mirt::extract.mirt(fit, "K")),
                                   as.integer(K))) {
      return(list(ok = FALSE, conv = FALSE, values = na, se = na))
    }
    fp <- .par_table(fit)
    vals <- .align_reflection(fp$value, fp$name)
    conv <- isTRUE(mirt::extract.mirt(fit, "converged"))
    se <- na
    if (conv) {
      # Keep estimates even if the subsequent information calculation fails.
      args$SE <- TRUE
      args$TOL <- NaN
      args$pars <- mirt::mod2values(fit)
      info <- .samplesize_fit(args)
      V <- if (is.null(info)) NULL else mirt::extract.mirt(info, "vcov")
      se <- .samplesize_reporting_se(fp, V)[match(params, fp$param)]
    }
    list(ok = TRUE, conv = conv, values = vals[match(params, fp$param)], se = se)
  }, parallel = parallel)

  matrix_for <- function(key) {
    out <- t(vapply(res, `[[`, numeric(length(params)), key))
    if (length(params) == 1L) out <- matrix(out, ncol = 1L)
    colnames(out) <- params
    out
  }
  est <- matrix_for("values")
  se <- matrix_for("se")
  ok <- vapply(res, `[[`, logical(1), "ok")
  conv <- vapply(res, `[[`, logical(1), "conv")
  estimates <- cbind(jobs, ok = ok, converged = conv,
                     as.data.frame(est, check.names = FALSE))
  standard_errors <- cbind(jobs, as.data.frame(se, check.names = FALSE))

  rows <- list()
  for (k in seq_along(params)) {
    p <- params[k]
    row <- pt[pt$param == p, ]
    for (n in N) {
      sel <- jobs$N == n
      z <- .samplesize_precision(est[sel, k], se[sel, k], ok[sel], conv[sel])
      e <- est[sel, k][z$used]
      at_b <- (is.finite(row$lbound) & e - row$lbound < 0.01) |
              (is.finite(row$ubound) & row$ubound - e < 0.01)
      rows[[length(rows) + 1L]] <- data.frame(
        param = p, N = n, current = row$value,
        mean = z$mean, sd = z$sd, mirt_se = z$mirt_se, n_used = sum(z$used),
        success = mean(ok[sel]), conv = mean(conv[sel]),
        se_failure = z$se_failure,
        p_bound = if (length(e)) mean(at_b) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  structure(list(summary = do.call(rbind, rows), estimates = estimates,
                 standard_errors = standard_errors, reps = reps,
                 current_N = nrow(mirt::extract.mirt(mod, "data")),
                 generating = generating), class = "irtident_samplesize")
}

.samplesize_fit <- function(args) {
  tryCatch(suppressWarnings(suppressMessages(do.call(mirt::mirt, args))),
           error = function(e) NULL)
}

.samplesize_reporting_se <- function(pt, V) {
  out <- rep(NA_real_, nrow(pt))
  if (!is.matrix(V) || is.null(colnames(V))) return(out)
  pn <- suppressWarnings(as.integer(sub(".*\\.", "", colnames(V))))
  variance <- diag(V)[match(pt$parnum, pn)]
  valid <- is.finite(variance) & variance > 0
  out[valid] <- sqrt(variance[valid])
  logit <- .is_logit_par(pt$name)
  out[logit] <- out[logit] * pt$value[logit] * (1 - pt$value[logit])
  out[logit & (pt$value <= 0 | pt$value >= 1)] <- NA_real_
  out[!is.finite(out) | out <= 0] <- NA_real_
  out
}

.samplesize_precision <- function(est, se, ok, conv) {
  eligible <- ok & conv
  used <- eligible & is.finite(est) & is.finite(se) & se > 0
  list(used = used,
       mean = if (any(used)) mean(est[used]) else NA_real_,
       sd = if (sum(used) >= 2L) stats::sd(est[used]) else NA_real_,
       mirt_se = if (any(used)) sqrt(mean(se[used]^2)) else NA_real_,
       se_failure = if (any(eligible)) mean(!used[eligible]) else NA_real_)
}
