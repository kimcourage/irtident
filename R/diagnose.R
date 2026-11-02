#' Run all identification diagnostics and summarise them by item
#'
#' Runs [id_local()], then [id_profile()] on the parameters it flags, and
#' [id_multistart()], and combines the evidence into one verdict per
#' parameter and per item:
#'
#' * `"ok"`: no diagnostic raised a concern;
#' * `"weak"`: the diagnostics indicate imprecision (high variance
#'   inflation, little information, or a profile interval open on one side);
#' * `"boundary"`: the estimate sits at a bound (e.g. `g` = 0), so the item
#'   behaves like a simpler model;
#' * `"not identified"`: neither side of the profile crosses the cutoff,
#'   or equivalent-fit solutions differ beyond the specified tolerance.
#'
#' These are diagnostic labels, not proofs of structural identification.
#' Each flagged parameter also gets a suggested follow-up check.
#' [id_samplesize()] is not run here because it is much slower; call it on
#' the parameters of interest afterwards.
#'
#' @param mod A unidimensional single-group model fitted with [mirt::mirt()].
#' @param profile,multistart Whether to run the profile and multi-start
#'   diagnostics.
#' @param profile_params Passed to `params` of [id_profile()].
#' @param n_grid Grid size for [id_profile()].
#' @param n_starts Number of random starts for [id_multistart()].
#' @param parallel Passed to the individual diagnostics.
#' @param ... Arguments passed to [id_local()] (`vif_cut`, `se_cut`,
#'   `boundary`).
#'
#' @return An object of class `irtident` with elements `items`, `params`,
#'   `local`, `profile`, `multistart` and `model_info`.
#' @examples
#' \donttest{
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' fit <- mirt(dat, 1, itemtype = "3PL", verbose = FALSE)
#' set.seed(1)
#' dx <- id_diagnose(fit, n_starts = 10)
#' dx
#' plot(dx)
#' }
#' @export
id_diagnose <- function(mod, profile = TRUE, multistart = TRUE,
                        profile_params = NULL, n_grid = 11, n_starts = 20,
                        parallel = FALSE, ...) {
  .check_mod(mod)
  local <- id_local(mod, ...)
  prof <- NULL
  flagged <- any(local$params$status %in% c("weak", "boundary"))
  if (profile && (flagged || !is.null(profile_params))) {
    prof <- id_profile(mod, params = profile_params, n_grid = n_grid,
                       local = local, parallel = parallel)
  }
  ms <- if (multistart) id_multistart(mod, n_starts = n_starts, parallel = parallel)
  params <- .combine_evidence(local, prof, ms)

  items <- do.call(rbind, lapply(unique(params$item), function(it) {
    p <- params[params$item == it, ]
    bad <- p[p$verdict != "ok", ]
    data.frame(item = it, itemtype = p$itemtype[1], verdict = .worst(p$verdict),
               issues = if (nrow(bad)) paste0(bad$name, ": ", bad$issue, collapse = "; ")
                        else "",
               stringsAsFactors = FALSE)
  }))
  info <- list(
    itemtype = mirt::extract.mirt(mod, "itemtype"),
    N = nrow(mirt::extract.mirt(mod, "data")),
    loglik = mirt::extract.mirt(mod, "logLik"),
    converged = isTRUE(mirt::extract.mirt(mod, "converged")),
    condnum = local$mirt_condnum
  )
  structure(list(items = items, params = params, local = local, profile = prof,
                 multistart = ms, model_info = info),
            class = "irtident")
}

.combine_evidence <- function(local, prof, ms) {
  lp <- local$params
  out <- data.frame(param = lp$param, item = lp$item, name = lp$name,
                    itemtype = lp$itemtype, est = lp$value, se = lp$se,
                    vif = lp$vif, local = lp$status, local_flags = lp$flags,
                    stringsAsFactors = FALSE)
  out$profile <- NA_character_
  out$lower <- out$upper <- NA_real_
  if (!is.null(prof) && nrow(prof$summary)) {
    m <- match(out$param, prof$summary$param)
    out$profile <- prof$summary$status[m]
    out$lower <- prof$summary$lower[m]
    out$upper <- prof$summary$upper[m]
    improved <- prof$summary$improved[m]
  } else {
    improved <- rep(FALSE, nrow(out))
  }
  out$unstable <- NA
  if (!is.null(ms)) out$unstable <- ms$params$unstable[match(out$param, ms$params$param)]

  verdict <- out$local
  verdict[verdict == "unknown"] <- "weak"
  verdict[out$profile %in% "partial" & verdict == "ok"] <- "weak"
  verdict[out$profile %in% "flat" & verdict != "boundary"] <- "not identified"
  verdict[out$unstable %in% TRUE & verdict != "boundary"] <- "not identified"
  out$verdict <- verdict

  partner <- .confounding_partner(local)
  out$issue <- ""
  out$advice <- ""
  for (i in which(out$verdict != "ok")) {
    r <- out[i, ]
    iss <- character()
    adv <- character()
    if (r$verdict == "boundary") {
      iss <- c(iss, paste0("at bound (", .fmt(r$est), ")"))
      adv <- c(adv, switch(r$name,
        g = "consider a simpler model (fix g = 0), or assess sensitivity to a prior on g",
        u = "consider fixing u = 1, or assess sensitivity to a prior on u",
        "check the model; consider fixing this parameter"))
    }
    if (!is.na(r$profile) && r$profile == "flat") {
      iss <- c(iss, "flat profile likelihood")
      if (r$verdict != "boundary") {
        adv <- c(adv, "check a wider profile grid and convergence; assess sensitivity to parameter constraints or priors")
      }
    } else if (!is.na(r$profile) && r$profile == "bounded" && r$verdict != "boundary") {
      iss <- c(iss, paste0("profile interval reaches the bound (", r$name, " = ",
                           .fmt(r$lower, 2), ")"))
      adv <- c(adv, paste0("the data do not rule out ", r$name, " = ", .fmt(r$lower, 2),
                           "; report the profile interval rather than the Wald interval"))
    } else if (!is.na(r$profile) && r$profile == "partial") {
      iss <- c(iss, "profile interval open on one side")
      if (r$name == "a1") adv <- c(adv, "slope may be unbounded (near-Guttman item); consider a prior on a1")
    }
    if (grepl("confounded", r$local_flags)) {
      iss <- c(iss, paste0("confounded with ", partner[r$param], " (VIF ", round(r$vif, 1), ")"))
      if (!length(adv)) adv <- c(adv, paste0("estimate ", r$name, " and ", partner[r$param],
                                             "; inspect their joint uncertainty and sensitivity to priors"))
    }
    if (grepl("little information", r$local_flags)) iss <- c(iss, "little information")
    if (grepl("imprecise", r$local_flags)) iss <- c(iss, "large SE")
    if (isTRUE(r$unstable)) {
      iss <- c(iss, "differs between equally good solutions")
      adv <- c(adv, "results depend on starting values")
    }
    if (isTRUE(improved[i])) {
      iss <- c(iss, "a better optimum exists")
      adv <- c(adv, "refit with more iterations or other starting values")
    }
    out$issue[i] <- paste(unique(iss), collapse = ", ")
    out$advice[i] <- paste(unique(adv), collapse = "; ")
  }
  out
}

# most strongly correlated other parameter, for each parameter
.confounding_partner <- function(local) {
  R <- local$cor
  if (is.null(R)) return(stats::setNames(rep(NA_character_, nrow(local$params)),
                                         local$params$param))
  diag(R) <- 0
  stats::setNames(colnames(R)[apply(abs(R), 1, which.max)], rownames(R))
}
