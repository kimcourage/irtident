#' Item-level local identification diagnostics
#'
#' Uses the estimated covariance matrix of the parameters (the inverse of the
#' information matrix) to describe parameter precision and variance inflation.
#'
#' For each parameter \eqn{j} the sampling variance is split as
#' \deqn{SE_j^2 = SE_{cond,j}^2 \times VIF^{within}_j \times VIF^{between}_j}
#' where
#' * \eqn{SE_{cond,j} = 1/\sqrt{I_{jj}}} is the standard error the parameter
#'   would have if every other parameter were known (intrinsic information);
#' * \eqn{VIF^{within}_j} is the variance inflation caused by not knowing the
#'   other parameters *of the same item* (e.g. the association between
#'   the guessing parameter `g` and the intercept `d` in the 3PL);
#' * \eqn{VIF^{between}_j} is the extra inflation caused by not knowing the
#'   parameters of *other items* and of the latent distribution.
#'
#' A large total VIF indicates confounding in the local quadratic
#' approximation; a large `se_cond` indicates little conditional information.
#' Parameters sitting at a bound (e.g. `g` near
#' 0) are flagged separately because the usual asymptotics do not apply.
#'
#' Near-null directions of the correlation matrix of the estimates are
#' reported in `directions`: each lists the parameters that move together
#' along a weakly identified combination.
#'
#' Standard errors and VIFs are on mirt's internal scale, where `g` and `u`
#' are logits. Point estimates are reported on the usual scale.
#'
#' @param mod A unidimensional single-group model fitted with [mirt::mirt()].
#'   It does not need to have been fitted with `SE = TRUE`; if it was not,
#'   the information matrix is computed at the estimates.
#' @param vif_cut Total variance inflation factor above which a parameter is
#'   flagged as confounded.
#' @param se_cut Standard error (internal scale) above which a parameter is
#'   flagged as imprecise, and conditional standard error above which it is
#'   flagged as having little information.
#' @param boundary Distance from a bound (e.g. `g` < `boundary`) at which a
#'   parameter is flagged as being at the boundary.
#' @param eig_cut Eigenvalues of the correlation matrix of the estimates below
#'   this value define weakly identified directions. Parameters at a bound
#'   are left out of this eigen-analysis.
#' @param SE.type `NULL` (default) uses the covariance matrix stored in `mod`,
#'   or mirt's default information matrix if there is none. Any value
#'   accepted by the `SE.type` argument of [mirt::mirt()] recomputes it, e.g.
#'   `"Fisher"` for the expected information. The expected information
#'   depends only on the parameter values and N, so it can be used at the
#'   design stage: set the parameters to planned values (`pars` with
#'   `TOL = NaN`) and check identification before collecting data. It
#'   enumerates all response patterns, so it is only feasible for short
#'   tests.
#'
#' @return An object of class `irtident_local`, a list with
#'   \describe{
#'     \item{params}{one row per estimated parameter, with `se`, `se_cond`,
#'       `vif`, `vif_within`, `vif_between`, flags and a `status` of
#'       `"ok"`, `"weak"` or `"boundary"`;}
#'     \item{items}{one row per item, summarising the worst parameter;}
#'     \item{directions}{weakly identified parameter combinations;}
#'     \item{cor}{correlation matrix of the estimates;}
#'     \item{singular}{`TRUE` if the information matrix could not be
#'       inverted. This may reflect weak identification or numerical problems.}
#'   }
#' @seealso [id_profile()], [id_multistart()], [id_diagnose()]
#' @examples
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' fit <- mirt(dat, 1, itemtype = "3PL", verbose = FALSE)
#' loc <- id_local(fit)
#' loc
#' loc$params
#' @export
id_local <- function(mod, vif_cut = 10, se_cut = 1, boundary = 0.01,
                     eig_cut = 1 / vif_cut, SE.type = NULL) {
  .check_mod(mod)
  pt <- .par_table(mod)
  info <- .get_info(mod, SE.type = SE.type)
  V <- info$V
  mirt_condnum <- info$condnum

  at_bound <- (is.finite(pt$lbound) & pt$value - pt$lbound < boundary) |
              (is.finite(pt$ubound) & pt$ubound - pt$value < boundary)

  if (is.null(V)) {
    pt$se <- pt$se_cond <- pt$vif <- pt$vif_within <- pt$vif_between <- NA_real_
    pt$at_bound <- at_bound
    pt$flags <- ifelse(at_bound, "boundary", "")
    pt$status <- ifelse(at_bound, "boundary", "unknown")
    out <- list(params = pt, items = .local_items(pt, NULL),
                directions = list(), singular = TRUE,
                condnum = NA_real_, mirt_condnum = mirt_condnum,
                cuts = list(vif = vif_cut, se = se_cut, boundary = boundary,
                            eig = eig_cut))
    return(structure(out, class = "irtident_local"))
  }

  # align the covariance matrix with the parameter table via parnum
  pn <- as.integer(sub(".*\\.", "", colnames(V)))
  ord <- match(pt$parnum, pn)
  V <- V[ord, ord, drop = FALSE]
  dimnames(V) <- list(pt$param, pt$param)

  I <- .sym_inv(V)
  se <- sqrt(diag(V))
  se_cond <- 1 / sqrt(pmax(diag(I), .Machine$double.eps))

  # SE when the other parameters of the same item are known
  se_within <- numeric(nrow(pt))
  for (it in unique(pt$item)) {
    b <- which(pt$item == it)
    Ib <- .sym_inv(V[b, b, drop = FALSE])
    se_within[b] <- 1 / sqrt(pmax(diag(Ib), .Machine$double.eps))
  }
  vif <- se^2 / se_cond^2
  vif_within <- se^2 / se_within^2
  vif_between <- se_within^2 / se_cond^2

  pt$se <- se
  pt$se_cond <- se_cond
  pt$vif <- vif
  pt$vif_within <- vif_within
  pt$vif_between <- vif_between
  pt$at_bound <- at_bound

  # At a bound, report the boundary flag instead of quadratic diagnostics.
  confounded <- vif > vif_cut & !at_bound
  low_info <- se_cond > se_cut & !at_bound
  imprecise <- se > se_cut & !at_bound
  source <- ifelse(vif_within >= vif_between, "within", "between")
  flags <- character(nrow(pt))
  add <- function(f, cond, lab) ifelse(cond, ifelse(f == "", lab, paste(f, lab, sep = ", ")), f)
  flags <- add(flags, at_bound, "boundary")
  flags <- add(flags, confounded, paste0("confounded (", source, ")"))
  flags <- add(flags, low_info, "little information")
  flags <- add(flags, imprecise & !confounded & !low_info, "imprecise")
  pt$flags <- flags
  pt$status <- ifelse(at_bound, "boundary",
                      ifelse(confounded | low_info | imprecise, "weak", "ok"))

  # weakly identified directions
  R <- stats::cov2cor(V)
  inner <- which(!at_bound)
  e <- if (length(inner)) {
    eigen(R[inner, inner, drop = FALSE], symmetric = TRUE)
  } else {
    list(values = NA_real_, vectors = matrix(NA_real_, 0, 0))
  }
  small <- which(e$values < eig_cut)
  directions <- lapply(small, function(k) {
    v <- e$vectors[, k]
    keep <- order(-abs(v))
    keep <- keep[abs(v[keep]) >= 0.25]
    data.frame(param = pt$param[inner][keep], loading = round(v[keep], 3),
               stringsAsFactors = FALSE)
  })
  names(directions) <- if (length(small)) .fmt(e$values[small], 5) else NULL

  out <- list(params = pt, items = .local_items(pt, R), directions = directions,
              cor = R, singular = FALSE,
              condnum = max(e$values) / max(min(e$values), 1e-300),
              mirt_condnum = mirt_condnum,
              cuts = list(vif = vif_cut, se = se_cut, boundary = boundary,
                          eig = eig_cut))
  structure(out, class = "irtident_local")
}

.status_rank <- c(ok = 0, unknown = 1, weak = 2, boundary = 3,
                  `not identified` = 4)

.worst <- function(s) {
  s <- s[!is.na(s)]
  if (!length(s)) return(NA_character_)
  s[which.max(.status_rank[s])]
}

.local_items <- function(pt, R) {
  items <- unique(pt$item)
  rows <- lapply(items, function(it) {
    b <- which(pt$item == it)
    kappa <- NA_real_
    if (!is.null(R) && length(b) > 1L) {
      ev <- eigen(R[b, b, drop = FALSE], symmetric = TRUE, only.values = TRUE)$values
      kappa <- max(ev) / max(min(ev), 1e-300)
    } else if (!is.null(R)) {
      kappa <- 1
    }
    w <- if (all(is.na(pt$vif[b]))) b[1] else b[which.max(pt$vif[b])]
    fl <- unique(unlist(strsplit(pt$flags[b][pt$flags[b] != ""], ", ")))
    data.frame(item = it, itemtype = pt$itemtype[b[1]], n_par = length(b),
               max_vif = if (all(is.na(pt$vif[b]))) NA_real_ else max(pt$vif[b]),
               worst_param = pt$name[w], kappa = kappa,
               status = .worst(pt$status[b]),
               flags = paste(fl, collapse = "; "),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
