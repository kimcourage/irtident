# print and plot methods ------------------------------------------------------

.round_df <- function(df, digits = 3) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(x) round(x, digits))
  df
}

#' @export
print.irtident_local <- function(x, ...) {
  cat("<irtident_local>", nrow(x$params), "estimated parameters,",
      sum(x$items$item != "GROUP"), "items\n")
  if (x$singular) {
    cat("The information matrix could not be inverted; check identification",
        "and numerical stability.\n")
  } else {
    cat("Condition number of the correlation matrix of estimates:",
        format(signif(x$condnum, 3)), "\n")
  }
  cat("\nItems:\n")
  print(.round_df(x$items, 2), row.names = FALSE)
  if (length(x$directions)) {
    cat("\nWeakly identified directions (eigenvalue: parameters, loadings):\n")
    for (k in seq_along(x$directions)) {
      d <- x$directions[[k]]
      cat("  ", names(x$directions)[k], ": ",
          paste0(d$param, " (", .fmt(d$loading, 2), ")", collapse = ", "),
          "\n", sep = "")
    }
  }
  cat("\nSee $params for parameter-level results.\n")
  invisible(x)
}

#' @export
print.irtident_profile <- function(x, ...) {
  cat("<irtident_profile> ", nrow(x$summary), " parameter(s), ",
      round(100 * x$level), "% profile likelihood intervals\n\n", sep = "")
  if (!nrow(x$summary)) return(invisible(x))
  s <- x$summary
  tab <- data.frame(param = s$param, est = .fmt(s$est),
                    profile = paste0("[", .fmt(s$lower), ", ", .fmt(s$upper), "]"),
                    sides = paste(s$lower_type, s$upper_type, sep = "/"),
                    wald = paste0("[", .fmt(s$wald_lower), ", ", .fmt(s$wald_upper), "]"),
                    max_dev = .fmt(s$max_dev, 2), status = s$status,
                    stringsAsFactors = FALSE)
  print(tab, row.names = FALSE)
  if (any(s$improved)) {
    cat("\nNote: the profile of", paste(s$param[s$improved], collapse = ", "),
        "reached a log-likelihood more than", x$tol_ll,
        "above the original fit, which was therefore not the maximum.\n")
  }
  invisible(x)
}

#' @export
print.irtident_multistart <- function(x, ...) {
  cat("<irtident_multistart>", x$n_starts, "random starts")
  if (x$n_failed) cat(",", x$n_failed, "failed")
  if (x$n_nonconv) cat(",", x$n_nonconv, "hit the iteration limit")
  cat("\n\n")
  print(.round_df(x$clusters, 4), row.names = FALSE)
  eq <- sum(x$clusters$type == "equivalent fit")
  cat("\n")
  if (eq > 1) {
    cat(eq, "clusters reach the same log-likelihood (within",
        x$tol[["ll"]], ") with different parameter values.\n")
  }
  if (!x$original_is_best) {
    cat("The original fit is not among the best solutions (best log-likelihood is",
        format(signif(x$improvement, 3)), "higher).\n")
  }
  un <- x$params[x$params$unstable, ]
  if (nrow(un)) {
    cat("Parameters that differ between equally good solutions:\n")
    print(.round_df(un[, c("param", "est", "spread_top")], 3), row.names = FALSE)
  } else {
    cat("No parameter differs between equally good solutions.\n")
  }
  invisible(x)
}

#' @export
print.irtident_samplesize <- function(x, ...) {
  cat("<irtident_samplesize>", x$reps, "replications per sample size;",
      "current N =", x$current_N, "\n\n")
  s <- x$summary[, c("param", "N", "current", "mean", "sd", "mirt_se",
                     "n_used", "success", "conv", "se_failure")]
  print(.round_df(s, 3), row.names = FALSE)
  cat("\nmean, sd and mirt_se use the same n_used converged fits with usable SEs.\n",
      "mirt_se = sqrt(mean(SE^2)); se_failure is conditional on convergence.\n",
      "Estimates and SEs are on the reporting scale.\n", sep = "")
  invisible(x)
}

#' @export
print.irtident <- function(x, ...) {
  mi <- x$model_info
  types <- paste(unique(mi$itemtype), collapse = "/")
  cat("irtident identification diagnostics\n")
  cat("Model: ", types, ", ", length(mi$itemtype), " items, N = ", mi$N,
      ", logLik = ", format(round(mi$loglik, 2), nsmall = 2),
      ", converged: ", if (mi$converged) "yes" else "no", "\n", sep = "")
  if (!is.na(mi$condnum)) {
    cat("mirt information-matrix condition number: ",
        formatC(mi$condnum, format = "g", digits = 3),
        "\n",
        sep = "")
  }
  cat("\nItems:\n")
  w_item <- max(nchar(x$items$item))
  w_type <- max(nchar(x$items$itemtype))
  w_par <- max(c(1L, nchar(x$params$name[x$params$verdict != "ok"])))
  for (i in seq_len(nrow(x$items))) {
    it <- x$items[i, ]
    cat("  ", formatC(it$item, width = -w_item), "  ",
        formatC(it$itemtype, width = -w_type), "  ", it$verdict, "\n", sep = "")
    p <- x$params[x$params$item == it$item & x$params$verdict != "ok", ]
    if (nrow(p)) {
      for (j in seq_len(nrow(p))) {
        cat("      ", formatC(p$name[j], width = -w_par), "  ", p$issue[j], "\n",
            sep = "")
      }
    }
  }
  adv <- x$params[x$params$advice != "", c("param", "advice")]
  if (nrow(adv)) {
    cat("\nAdvice:\n")
    for (i in seq_len(nrow(adv))) cat("  ", adv$param[i], ": ", adv$advice[i], "\n", sep = "")
  }
  if (!is.null(x$multistart)) {
    ms <- x$multistart
    cat("\nMulti-start: ", ms$n_starts, " starts, ",
        sum(ms$clusters$type == "equivalent fit"), " equivalent-fit cluster(s), ",
        sum(ms$clusters$type == "inferior optimum"), " inferior optimum/optima; ",
        if (ms$original_is_best) "the original fit is among the best solutions."
        else "the original fit is NOT among the best solutions.", "\n", sep = "")
  }
  cat("\nDetails: $params, $local, $profile, $multistart; plot(x) shows the profiles.\n")
  invisible(x)
}

#' Plot profile likelihoods
#'
#' Draws the profile deviance of each profiled parameter with the
#' \eqn{\chi^2_1} cutoff (dashed) and the estimate (dotted).
#'
#' @param x An `irtident_profile` or `irtident` object.
#' @param params Optional subset of parameter labels to plot.
#' @param ... Further graphical parameters passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.irtident_profile <- function(x, params = NULL, ...) {
  s <- x$summary
  if (!is.null(params)) s <- s[s$param %in% params, ]
  if (!nrow(s)) {
    message("Nothing to plot.")
    return(invisible(x))
  }
  n <- nrow(s)
  nc <- ceiling(sqrt(n))
  op <- graphics::par(mfrow = c(ceiling(n / nc), nc), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(op))
  for (i in seq_len(n)) {
    cv <- x$curves[[s$param[i]]]
    cv <- cv[!is.na(cv$deviance), ]
    ylim <- c(min(0, cv$deviance), max(x$cutoff * 1.5, cv$deviance))
    graphics::plot(cv$value, cv$deviance, type = "b", pch = 19, cex = 0.6,
                   xlab = s$name[i], ylab = "profile deviance", ylim = ylim,
                   main = paste0(s$param[i], " (", s$status[i], ")"),
                   cex.main = 0.95, ...)
    graphics::abline(h = x$cutoff, lty = 2, col = "grey40")
    graphics::abline(v = s$est[i], lty = 3, col = "grey40")
  }
  invisible(x)
}

#' @rdname plot.irtident_profile
#' @export
plot.irtident <- function(x, params = NULL, ...) {
  if (is.null(x$profile)) {
    message("No profile diagnostics to plot.")
    return(invisible(x))
  }
  plot(x$profile, params = params, ...)
  invisible(x)
}

#' Plot multi-start solutions
#'
#' Shows the log-likelihood of every start relative to the best one, sorted,
#' with converged solutions filled and the original fit marked.
#'
#' @param x An `irtident_multistart` object.
#' @param ... Further graphical parameters passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.irtident_multistart <- function(x, ...) {
  s <- x$solutions[!is.na(x$solutions$loglik), ]
  s$gap <- max(s$loglik) - s$loglik
  s <- s[order(s$gap), ]
  graphics::plot(seq_len(nrow(s)), s$gap, pch = ifelse(s$converged, 19, 1),
                 xlab = "solution (sorted)",
                 ylab = "log-likelihood below best", ...)
  graphics::abline(h = x$tol[["ll"]], lty = 2, col = "grey40")
  o <- which(s$start == "original")
  graphics::points(o, s$gap[o], pch = 4, cex = 2, col = "red")
  graphics::legend("topleft", bty = "n", cex = 0.8,
                   pch = c(19, 1, 4), col = c("black", "black", "red"),
                   legend = c("converged", "iteration limit", "original fit"))
  invisible(x)
}

#' Plot sample-size simulation results
#'
#' Repeated-estimate SD and mirt SE against sample size for each parameter.
#' Both curves use the same replications and reporting scale. Inspect
#' `n_used`, `conv` and `se_failure` in the printed summary alongside the plot.
#'
#' @param x An `irtident_samplesize` object.
#' @param ... Further graphical parameters passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.irtident_samplesize <- function(x, ...) {
  s <- x$summary
  ps <- unique(s$param)
  n <- length(ps)
  nc <- ceiling(sqrt(n))
  op <- graphics::par(mfrow = c(ceiling(n / nc), nc), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(op))
  for (p in ps) {
    d <- s[s$param == p, ]
    d <- d[order(d$N), ]
    yl <- range(c(0, d$sd, d$mirt_se), na.rm = TRUE)
    if (diff(yl) == 0) yl <- c(0, 1)
    graphics::plot(d$N, d$sd, type = "b", pch = 19, log = "x", ylim = yl,
                   xlab = "N", ylab = "SD / SE", main = p, cex.main = 0.95, ...)
    graphics::lines(d$N, d$mirt_se, type = "b", pch = 17, lty = 2,
                    col = "#0072B2")
    graphics::legend("topright", legend = c("Repeated-estimate SD", "mirt SE (RMS)"),
                     col = c("black", "#0072B2"), pch = c(19, 17),
                     lty = c(1, 2), bty = "n", cex = 0.75)
  }
  invisible(x)
}
