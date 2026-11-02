#' Multi-start estimation with clustering of the solutions
#'
#' Refits the model from `n_starts` random starting values and groups the
#' solutions that agree (up to `tol_par` in every parameter). Each cluster is
#' classified as
#'
#' * **equivalent fit**: its best log-likelihood is within `tol_ll` of the
#'   overall best. Several such clusters mean that different parameter values
#'   fit the data equally well.
#' * **inferior optimum**: its best log-likelihood is lower.
#'
#' The reflection \eqn{\theta \to -\theta} (all slopes change sign) gives an
#' identical likelihood in a unidimensional model and is removed before
#' clustering.
#'
#' @param mod A unidimensional single-group model fitted with [mirt::mirt()].
#' @param n_starts Number of random starts.
#' @param tol_par Two solutions are put in the same cluster when no
#'   parameter (reporting scale) differs by more than this.
#' @param tol_ll Clusters whose best log-likelihood is within `tol_ll` of the
#'   overall best count as equivalent fit.
#' @param unstable_cut Named vector of spreads (reporting scale, by parameter
#'   type; graded intercepts use `d`) above which a parameter is flagged as
#'   unstable across equivalent-fit clusters.
#' @param parallel If `TRUE`, starts are run in parallel with `future.apply`.
#' @param ... Further arguments passed to [mirt::mirt()] for every refit.
#'
#' @return An object of class `irtident_multistart`, a list with
#'   \describe{
#'     \item{solutions}{one row per start: log-likelihood, convergence and
#'       cluster. Solutions that reached the iteration limit are included and
#'       marked in `converged`;}
#'     \item{estimates}{matrix of estimates, starts in rows (reporting
#'       scale);}
#'     \item{clusters}{one row per cluster with its size, best
#'       log-likelihood, gap to the overall best and type;}
#'     \item{params}{one row per parameter with its spread across the
#'       equivalent-fit clusters (`spread_top`) and across all solutions
#'       (`spread_all`), and an `unstable` flag;}
#'     \item{original_is_best}{whether the original fit is among the best
#'       solutions (in an equivalent-fit cluster).}
#'   }
#' @examples
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' fit <- mirt(dat, 1, itemtype = "3PL", verbose = FALSE)
#' set.seed(1)
#' ms <- id_multistart(fit, n_starts = 4)
#' ms
#' @export
id_multistart <- function(mod, n_starts = 20, tol_par = 0.05, tol_ll = 0.1,
                          unstable_cut = c(a1 = 0.5, d = 0.5, g = 0.1, u = 0.1,
                                           COV_11 = 0.25),
                          parallel = FALSE, ...) {
  .check_mod(mod)
  pt <- .par_table(mod)
  base <- mirt::mod2values(mod)
  est_rows <- match(pt$parnum, base$parnum)

  # draw all starting values up front so results do not depend on `parallel`
  starts <- lapply(seq_len(n_starts), function(i) .random_start(base))
  extra <- list(...)
  fits <- .lapply(starts, function(pv) {
    fit <- do.call(.refit, c(list(mod = mod, pars = pv), extra))
    if (is.null(fit)) return(NULL)
    vals <- mirt::mod2values(fit)$value[est_rows]
    list(loglik = mirt::extract.mirt(fit, "logLik"),
         converged = isTRUE(mirt::extract.mirt(fit, "converged")),
         values = .align_reflection(vals, pt$name))
  }, parallel = parallel)

  orig <- list(loglik = mirt::extract.mirt(mod, "logLik"),
               converged = isTRUE(mirt::extract.mirt(mod, "converged")),
               values = .align_reflection(pt$value, pt$name))
  all_fits <- c(list(orig), fits)
  ok <- !vapply(all_fits, is.null, logical(1))
  ll <- vapply(all_fits, function(f) if (is.null(f)) NA_real_ else f$loglik, 0)
  conv <- vapply(all_fits, function(f) if (is.null(f)) FALSE else f$converged, TRUE)
  est <- t(vapply(all_fits, function(f) {
    if (is.null(f)) rep(NA_real_, nrow(pt)) else f$values
  }, numeric(nrow(pt))))
  colnames(est) <- pt$param
  rownames(est) <- c("original", paste0("start", seq_len(n_starts)))

  # Keep returned estimates for comparison, with convergence reported separately.
  use <- which(ok)
  cluster <- rep(NA_integer_, length(all_fits))
  if (length(use) > 1L) {
    d <- stats::dist(est[use, , drop = FALSE], method = "maximum")
    cluster[use] <- stats::cutree(stats::hclust(d, method = "complete"), h = tol_par)
  } else if (length(use) == 1L) {
    cluster[use] <- 1L
  }
  # relabel clusters by decreasing best log-likelihood
  best_ll <- tapply(ll[use], cluster[use], max)
  relabel <- stats::setNames(seq_along(best_ll), names(sort(best_ll, decreasing = TRUE)))
  cluster[use] <- relabel[as.character(cluster[use])]

  ll_best <- max(ll[use])
  clusters <- do.call(rbind, lapply(sort(unique(cluster[use])), function(k) {
    m <- which(cluster == k)
    gap <- ll_best - max(ll[m])
    data.frame(cluster = k, size = length(m), n_converged = sum(conv[m]),
               loglik = max(ll[m]), gap = gap,
               type = if (gap <= tol_ll) "equivalent fit" else "inferior optimum",
               has_original = 1L %in% m, stringsAsFactors = FALSE)
  }))

  top <- clusters$cluster[clusters$type == "equivalent fit"]
  centre <- function(k) {
    m <- which(cluster == k)
    est[m[which.max(ll[m])], ]
  }
  top_centres <- do.call(rbind, lapply(top, centre))
  spread <- function(M) {
    if (is.null(M) || nrow(M) < 2L) return(rep(0, ncol(est)))
    apply(M, 2, function(x) diff(range(x, na.rm = TRUE)))
  }
  params <- data.frame(param = pt$param, est = pt$value,
                       spread_top = spread(top_centres),
                       spread_all = spread(est[use, , drop = FALSE]),
                       stringsAsFactors = FALSE)
  params$unstable <- params$spread_top >
    .target_for(pt$name, unstable_cut)
  rownames(params) <- NULL

  solutions <- data.frame(start = rownames(est), loglik = ll, converged = conv,
                          cluster = cluster, stringsAsFactors = FALSE)
  rownames(solutions) <- NULL
  structure(list(solutions = solutions, estimates = est, clusters = clusters,
                 params = params, original_is_best = isTRUE(clusters$type[clusters$has_original] ==
                                           "equivalent fit"),
                 improvement = ll_best - ll[1], n_starts = n_starts,
                 n_failed = sum(!ok[-1]), n_nonconv = sum(ok[-1] & !conv[-1]),
                 tol = c(par = tol_par, ll = tol_ll)),
            class = "irtident_multistart")
}

# Random starting values for the freely estimated parameters of a values table
.random_start <- function(pv) {
  est <- which(pv$est)
  for (i in est) {
    nm <- pv$name[i]
    pv$value[i] <- switch(
      nm,
      a1 = stats::runif(1, 0.3, 3),
      d = stats::runif(1, -2.5, 2.5),
      g = stats::runif(1, 0.01, 0.35),
      u = stats::runif(1, 0.75, 0.99),
      COV_11 = stats::runif(1, 0.3, 3),
      MEAN_1 = stats::runif(1, -1, 1),
      pv$value[i]
    )
  }
  # graded intercepts: random but strictly descending
  for (it in unique(pv$item[pv$class == "graded"])) {
    r <- which(pv$item == it & grepl("^d[0-9]+$", pv$name) & pv$est)
    if (!length(r)) next
    k <- as.integer(sub("^d", "", pv$name[r]))
    r <- r[order(k)]
    pv$value[r] <- sort(stats::runif(length(r), -3, 3), decreasing = TRUE)
  }
  pv
}
