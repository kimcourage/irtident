# Internal helpers -----------------------------------------------------------

# Look up a threshold by parameter type; graded intercepts use "d".
.target_for <- function(name, values) {
  key <- ifelse(grepl("^d[0-9]+$", name), "d", name)
  out <- unname(values[key])
  out[is.na(out)] <- Inf
  out
}

.supported_itemtypes <- c("Rasch", "2PL", "3PL", "3PLu", "4PL", "graded")

#' Validate a fitted mirt model for use with irtident
#' @noRd
.check_mod <- function(mod) {
  if (!inherits(mod, "SingleGroupClass")) {
    stop("`mod` must be a single-group model fitted with mirt::mirt().",
         call. = FALSE)
  }
  if (mirt::extract.mirt(mod, "nfact") != 1L) {
    stop("Only unidimensional models are supported in this version.",
         call. = FALSE)
  }
  it <- mirt::extract.mirt(mod, "itemtype")
  bad <- setdiff(unique(it), .supported_itemtypes)
  if (length(bad)) {
    stop("Unsupported itemtype(s): ", paste(bad, collapse = ", "),
         ". Supported: ", paste(.supported_itemtypes, collapse = ", "), ".",
         call. = FALSE)
  }
  pv <- mirt::mod2values(mod)
  if (any(pv$const != "none")) {
    stop("Models with equality constraints are not supported in this version.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Is a parameter name estimated on the logit scale internally by mirt?
#' @noRd
.is_logit_par <- function(name) name %in% c("g", "u")

#' Map reporting-scale values to mirt's internal (optimisation) scale
#' @noRd
.to_internal <- function(value, name) {
  out <- value
  lg <- .is_logit_par(name)
  out[lg] <- stats::qlogis(pmin(pmax(value[lg], 1e-12), 1 - 1e-12))
  out
}

#' Table of freely estimated parameters
#'
#' Returns one row per estimated parameter with the reporting-scale value
#' (as in `mirt::mod2values()`), the internal value and a label such as
#' `"Item.2:g"`.
#' @noRd
.par_table <- function(mod) {
  pv <- mirt::mod2values(mod)
  it <- mirt::extract.mirt(mod, "itemtype")
  items <- colnames(mirt::extract.mirt(mod, "data"))
  pv <- pv[pv$est, , drop = FALSE]
  data.frame(
    param    = paste0(pv$item, ":", pv$name),
    item     = pv$item,
    name     = pv$name,
    itemtype = ifelse(pv$item %in% items, it[match(pv$item, items)], "group"),
    parnum   = pv$parnum,
    value    = pv$value,
    internal = .to_internal(pv$value, pv$name),
    lbound   = pv$lbound,
    ubound   = pv$ubound,
    stringsAsFactors = FALSE
  )
}

#' Refit a model with a modified parameter table
#'
#' Returns the fitted model, or `NULL` if mirt fails. Warnings and messages
#' are suppressed because refits are expected to visit awkward regions.
#' @noRd
.refit <- function(mod, pars, data = NULL, SE = FALSE, ...) {
  if (is.null(data)) data <- mirt::extract.mirt(mod, "data")
  args <- list(data = data, model = 1,
               itemtype = mirt::extract.mirt(mod, "itemtype"),
               pars = pars, verbose = FALSE, SE = SE)
  args <- utils::modifyList(args, list(...))
  tryCatch(
    suppressWarnings(suppressMessages(do.call(mirt::mirt, args))),
    error = function(e) NULL
  )
}

#' Covariance matrix of the estimated parameters (internal scale)
#'
#' Uses the stored matrix when the model was fitted with `SE = TRUE` and no
#' `SE.type` is requested; otherwise re-evaluates the model at its current
#' parameter values with `TOL = NaN`.
#' @noRd
.get_info <- function(mod, SE.type = NULL) {
  src <- mod
  V <- if (is.null(SE.type)) mirt::extract.mirt(mod, "vcov") else NULL
  if (is.null(V) || length(V) <= 1L || anyNA(V)) {
    args <- list(mod = mod, pars = mirt::mod2values(mod), SE = TRUE, TOL = NaN)
    if (!is.null(SE.type)) args$SE.type <- SE.type
    src <- do.call(.refit, args)
    V <- if (is.null(src)) NULL else mirt::extract.mirt(src, "vcov")
  }
  ext <- function(what) {
    if (is.null(src)) return(NA)
    tryCatch(mirt::extract.mirt(src, what), error = function(e) NA)
  }
  ok <- !(is.null(V) || length(V) <= 1L || anyNA(V))
  list(V = if (ok) V else NULL,
       condnum = ext("condnum"),
       secondorder = ext("secondordertest"))
}

#' Inverse of a symmetric matrix, falling back to a pseudo-inverse
#' @noRd
.sym_inv <- function(M) {
  out <- tryCatch(solve(M), error = function(e) NULL)
  if (!is.null(out)) return(out)
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > max(e$values) * 1e-12
  e$vectors[, keep, drop = FALSE] %*%
    (t(e$vectors[, keep, drop = FALSE]) / e$values[keep])
}

#' lapply that optionally runs in parallel through future.apply
#' @noRd
.lapply <- function(X, FUN, parallel = FALSE) {
  if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      warning("`parallel = TRUE` needs the 'future.apply' package; ",
              "running sequentially.", call. = FALSE)
    } else {
      return(future.apply::future_lapply(X, FUN, future.seed = TRUE))
    }
  }
  lapply(X, FUN)
}

#' Resolve a `params` argument to a vector of parameter labels
#' @noRd
.resolve_params <- function(params, pt, default) {
  if (is.null(params)) return(default)
  if (identical(params, "all")) return(pt$param[pt$item != "GROUP"])
  if (is.data.frame(params)) params <- paste0(params$item, ":", params$name)
  params <- as.character(params)
  unknown <- setdiff(params, pt$param)
  if (length(unknown)) {
    stop("Unknown or fixed parameter(s): ", paste(unknown, collapse = ", "),
         ". Use labels of the form 'item:name', e.g. '",
         pt$param[1], "'.", call. = FALSE)
  }
  params
}

#' Flip slopes when a solution is the reflection theta -> -theta
#' @noRd
.align_reflection <- function(values, names) {
  slope <- names == "a1"
  if (any(slope) && sum(values[slope]) < 0) values[slope] <- -values[slope]
  values
}

.fmt <- function(x, digits = 3) formatC(x, digits = digits, format = "f")
