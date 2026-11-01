#' Profile likelihood diagnostics for selected parameters
#'
#' For each selected parameter, fixes it at a grid of values, re-estimates
#' all other parameters, and records the profile log-likelihood. The profile
#' is traced outward from the estimate, each refit starting from the previous
#' solution (continuation), which keeps refits fast and stable.
#'
#' Following Raue et al. (2009), a parameter is *practically identified* when
#' its profile likelihood interval is closed on both sides. Each side of the
#' interval is classified as
#' * `"closed"`: the deviance crosses the \eqn{\chi^2_1} cutoff;
#' * `"bound"`: the grid reaches a natural bound of the parameter (`g = 0`,
#'   `u = 1`, `a1 = 0`) without crossing the cutoff;
#' * `"open"`: the grid ends without crossing the cutoff.
#'
#' The resulting `status` is `"identified"` (both sides closed), `"bounded"`
#' (one side at a natural bound, the other closed: the data rule out large
#' values but not the bound itself), `"partial"` (one side open) or `"flat"`
#' (neither side closed: the data barely distinguish values across the
#' whole grid).
#'
#' Deviances are measured from the highest log-likelihood seen (the original
#' fit or any profile point). If a profile reaches a log-likelihood clearly
#' above the original fit (more than `tol_ll`), the original estimate was
#' not the maximum; this is reported in `improved`.
#'
#' @param mod A unidimensional single-group model fitted with [mirt::mirt()].
#' @param params Parameters to profile. `NULL` (default) profiles the
#'   parameters flagged by [id_local()]; `"all"` profiles every item
#'   parameter; otherwise a character vector of labels such as
#'   `"Item.2:g"`, or a data frame with columns `item` and `name`.
#' @param n_grid Number of grid points per parameter.
#' @param level Confidence level of the profile likelihood interval.
#' @param grid Optional named list of numeric vectors giving the grid for
#'   specific parameters (names are labels such as `"Item.2:g"`). Values are
#'   on the reporting scale (`g` and `u` as probabilities).
#' @param local An `irtident_local` object from [id_local()], reused to pick
#'   flagged parameters, size the grids and compute Wald intervals.
#' @param stop_dev Tracing in one direction stops once the deviance exceeds
#'   `stop_dev` times the cutoff; the remaining grid points carry no
#'   information about the interval and are left as `NA`.
#' @param tol_ll Log-likelihood gain above which a profile is reported as
#'   having found a better optimum than the original fit (`improved`).
#'   Smaller gains are treated as EM stopping slightly early on a flat ridge.
#' @param parallel If `TRUE`, parameters are profiled in parallel with
#'   `future.apply` under the current `future::plan()`.
#' @param ... Further arguments passed to [mirt::mirt()] for every refit
#'   (e.g. `technical = list(NCYCLES = 1000)`).
#'
#' @return An object of class `irtident_profile`: a list with `summary` (one
#'   row per profiled parameter), `curves` (the profile for each parameter)
#'   and the `cutoff` used.
#' @references Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M.,
#'   Klingmuller, U., & Timmer, J. (2009). Structural and practical
#'   identifiability analysis of partially observed dynamical models by
#'   exploiting the profile likelihood. *Bioinformatics, 25*(15), 1923-1929.
#' @seealso [id_local()], [plot.irtident_profile()]
#' @examples
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' fit <- mirt(dat, 1, itemtype = "3PL", verbose = FALSE)
#' prof <- id_profile(fit, params = "Item.2:g", n_grid = 7)
#' prof
#' plot(prof)
#' @export
id_profile <- function(mod, params = NULL, n_grid = 11, level = 0.95,
                       grid = NULL, local = NULL, stop_dev = 2, tol_ll = 0.1,
                       parallel = FALSE, ...) {
  .check_mod(mod)
  pt <- .par_table(mod)
  if (is.null(local)) local <- id_local(mod)
  flagged <- local$params$param[local$params$status %in% c("weak", "boundary")]
  params <- .resolve_params(params, pt, default = flagged)
  cutoff <- stats::qchisq(level, df = 1)
  ll_hat <- mirt::extract.mirt(mod, "logLik")

  if (!length(params)) {
    message("No parameter was flagged by id_local(); nothing to profile. ",
            "Use `params = \"all\"` to profile every item parameter.")
  }

  extra <- list(...)
  curves <- .lapply(params, function(p) {
    row <- pt[pt$param == p, ]
    lrow <- local$params[local$params$param == p, ]
    g <- if (!is.null(grid[[p]])) sort(unique(c(grid[[p]], row$value)))
         else .default_grid(row, lrow$se, n_grid)
    do.call(.trace_profile,
            c(list(mod = mod, row = row, grid = g,
                   stop_ll = ll_hat - stop_dev * cutoff / 2), extra))
  }, parallel = parallel)
  names(curves) <- params

  # EM stops early on flat ridges, so a profile can end slightly above the
  # original fit. Deviances are measured from the best log-likelihood seen.
  ll_seen <- unlist(lapply(curves, function(cv) cv$loglik))
  ll_ref <- max(c(ll_hat, ll_seen), na.rm = TRUE)
  rows <- lapply(params, function(p) {
    curves[[p]]$deviance <<- 2 * (ll_ref - curves[[p]]$loglik)
    s <- .summarise_profile(pt[pt$param == p, ],
                            local$params[local$params$param == p, ],
                            curves[[p]], cutoff, level)
    s$improved <- max(curves[[p]]$loglik, na.rm = TRUE) - ll_hat > tol_ll
    s
  })
  summary <- do.call(rbind, rows)
  if (is.null(summary)) summary <- .empty_profile_summary()
  structure(list(summary = summary, curves = curves, cutoff = cutoff,
                 level = level, loglik = ll_hat, loglik_ref = ll_ref,
                 improvement = ll_ref - ll_hat, tol_ll = tol_ll),
            class = "irtident_profile")
}

.natural_bounds <- function(row) {
  lo <- row$lbound
  hi <- row$ubound
  if (row$name == "a1") lo <- 0
  c(lo, hi)
}

.default_grid <- function(row, se, n_grid) {
  v0 <- row$value
  nm <- row$name
  if (nm == "g") {
    hi <- max(0.5, min(0.95, v0 + 0.2))
    g <- c(0, seq(0.01, hi, length.out = n_grid - 1))
  } else if (nm == "u") {
    lo <- min(0.5, max(0.05, v0 - 0.2))
    g <- c(seq(lo, 0.99, length.out = n_grid - 1), 1)
  } else if (grepl("^COV", nm)) {
    g <- exp(seq(log(v0 / 4), log(v0 * 4), length.out = n_grid))
  } else {
    hw <- if (is.na(se)) 1.5 else min(max(4 * se, 0.75), 4)
    lo <- v0 - hw
    if (nm == "a1") lo <- max(0, lo)
    g <- seq(lo, v0 + hw, length.out = n_grid)
  }
  sort(unique(c(g, v0)))
}

# keep graded intercepts strictly descending around a fixed one
.repair_graded <- function(pv, item, name, value, gap = 0.05) {
  if (!grepl("^d[0-9]+$", name)) return(pv)
  rows <- which(pv$item == item & grepl("^d[0-9]+$", pv$name))
  if (length(rows) < 2L) return(pv)
  k <- as.integer(sub("^d", "", pv$name[rows]))
  rows <- rows[order(k)]
  k <- sort(k)
  kf <- as.integer(sub("^d", "", name))
  for (j in seq_along(rows)) {
    if (k[j] < kf) {
      pv$value[rows[j]] <- max(pv$value[rows[j]], value + gap * (kf - k[j]))
    } else if (k[j] > kf) {
      pv$value[rows[j]] <- min(pv$value[rows[j]], value - gap * (k[j] - kf))
    }
  }
  pv
}

.trace_profile <- function(mod, row, grid, stop_ll = -Inf, ...) {
  v0 <- row$value
  base <- mirt::mod2values(mod)
  ri <- which(base$item == row$item & base$name == row$name)
  out <- data.frame(value = grid, loglik = NA_real_, converged = NA,
                    point = "skipped", stringsAsFactors = FALSE)
  at0 <- which(grid == v0)
  out$loglik[at0] <- mirt::extract.mirt(mod, "logLik")
  out$converged[at0] <- isTRUE(mirt::extract.mirt(mod, "converged"))
  out$point[at0] <- "estimate"
  walk <- function(idx) {
    cur <- base
    for (i in idx) {
      pv <- cur
      pv$est[ri] <- FALSE
      pv$value[ri] <- grid[i]
      pv <- .repair_graded(pv, row$item, row$name, grid[i])
      fit <- .refit(mod, pv, ...)
      if (is.null(fit)) {
        out$point[i] <<- "failed"
        next
      }
      ll <- mirt::extract.mirt(fit, "logLik")
      out$loglik[i] <<- ll
      out$converged[i] <<- isTRUE(mirt::extract.mirt(fit, "converged"))
      out$point[i] <<- "fitted"
      cur <- mirt::mod2values(fit)
      if (ll < stop_ll) break  # far past the cutoff: nothing more to learn
    }
  }
  walk(which(grid > v0))
  walk(rev(which(grid < v0)))
  out
}

.crossing <- function(x, d, cutoff) {
  # x ordered outward from the estimate, d the matching deviances
  ok <- !is.na(d)
  x <- x[ok]
  d <- d[ok]
  hit <- which(d >= cutoff)[1]
  if (is.na(hit)) return(NA_real_)
  if (hit == 1L) return(x[1])
  x[hit - 1] + (cutoff - d[hit - 1]) * (x[hit] - x[hit - 1]) / (d[hit] - d[hit - 1])
}

.summarise_profile <- function(row, lrow, curve, cutoff, level) {
  v0 <- row$value
  bnd <- .natural_bounds(row)
  up <- curve[curve$value >= v0, ]
  dn <- curve[curve$value <= v0, ]
  dn <- dn[order(-dn$value), ]
  upper <- .crossing(up$value, up$deviance, cutoff)
  lower <- .crossing(dn$value, dn$deviance, cutoff)
  side_type <- function(cross, edge, bound) {
    if (!is.na(cross)) return("closed")
    if (is.finite(bound) && abs(edge - bound) < 1e-8) return("bound")
    "open"
  }
  reached <- range(curve$value[!is.na(curve$loglik)])
  lower_type <- side_type(lower, reached[1], bnd[1])
  upper_type <- side_type(upper, reached[2], bnd[2])
  if (is.na(lower)) lower <- reached[1]
  if (is.na(upper)) upper <- reached[2]

  z <- stats::qnorm(1 - (1 - level) / 2)
  se <- if (nrow(lrow)) lrow$se else NA_real_
  if (.is_logit_par(row$name)) {
    wald <- stats::plogis(row$internal + c(-1, 1) * z * se)
  } else {
    wald <- row$value + c(-1, 1) * z * se
  }
  types <- c(lower_type, upper_type)
  status <- if (all(types == "closed")) "identified"
            else if (!any(types == "closed")) "flat"
            else if (any(types == "bound")) "bounded"
            else "partial"
  data.frame(param = row$param, item = row$item, name = row$name,
             est = v0, lower = lower, upper = upper,
             lower_type = lower_type, upper_type = upper_type,
             wald_lower = wald[1], wald_upper = wald[2],
             max_dev = max(curve$deviance, na.rm = TRUE),
             improved = FALSE,
             n_failed = sum(curve$point == "failed"),
             status = status, stringsAsFactors = FALSE)
}

.empty_profile_summary <- function() {
  data.frame(param = character(), item = character(), name = character(),
             est = numeric(), lower = numeric(), upper = numeric(),
             lower_type = character(), upper_type = character(),
             wald_lower = numeric(), wald_upper = numeric(),
             max_dev = numeric(), improved = logical(),
             n_failed = integer(), status = character(),
             stringsAsFactors = FALSE)
}
