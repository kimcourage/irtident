test_that("sample size simulation returns the expected structure", {
  set.seed(3)
  ss <- id_samplesize(get_fit("2PL"), params = c("Item.1:a1", "Item.1:d"),
                      N = c(300, 1000), reps = 3)
  expect_s3_class(ss, "irtident_samplesize")
  expect_equal(nrow(ss$summary), 4L)
  expect_named(ss, c("summary", "estimates", "standard_errors", "reps",
                     "current_N", "generating"))
  expect_false(any(c("target", "meets", "bias", "rmse", "true") %in%
                     names(ss$summary)))
  expect_equal(nrow(ss$estimates), 6L)
  expect_equal(nrow(ss$standard_errors), 6L)
  expect_equal(ss$current_N, nrow(lsat7))
  expect_true(all(ss$summary$success >= 0 & ss$summary$success <= 1))
  expect_true(all(ss$summary$n_used >= 2L))
  pt <- .par_table(get_fit("2PL"))
  for (i in seq_len(nrow(ss$summary))) {
    s <- ss$summary[i, ]
    use <- ss$estimates$N == s$N & ss$estimates$ok &
      ss$estimates$converged & is.finite(ss$estimates[[s$param]]) &
      is.finite(ss$standard_errors[[s$param]]) &
      ss$standard_errors[[s$param]] > 0
    z <- ss$estimates[use, s$param]
    se <- ss$standard_errors[use, s$param]
    expect_equal(s$current, pt$value[pt$param == s$param])
    expect_equal(s$mean, mean(z))
    expect_equal(s$sd, sd(z))
    expect_equal(s$mirt_se, sqrt(mean(se^2)))
    expect_equal(s$n_used, sum(use))
  }
  output <- paste(capture.output(print(ss)), collapse = "\n")
  expect_match(output, "mirt_se")
  expect_match(output, "se_failure")
  expect_false(grepl("Required sample size|target|meets|bias|rmse|true", output))
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(ss))
})

test_that("precision summaries compare the same usable replications", {
  z <- .samplesize_precision(c(1, 2, 100, -100), c(0.1, 0.2, NA, 0.4),
                             rep(TRUE, 4), c(TRUE, TRUE, TRUE, FALSE))
  expect_identical(z$used, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(z$mean, 1.5)
  expect_equal(z$sd, sqrt(0.5))
  expect_equal(z$mirt_se, sqrt(0.025))
  expect_equal(z$se_failure, 1 / 3)

  none <- .samplesize_precision(1, 0.1, TRUE, FALSE)
  expect_true(all(is.na(unlist(none[-1]))))
  one <- .samplesize_precision(1, 0.1, TRUE, TRUE)
  expect_true(is.na(one$sd))
  expect_equal(one$mirt_se, 0.1)
})

test_that("mirt SEs are matched by parameter number and converted to reporting scale", {
  pt <- data.frame(name = c("a1", "g", "u", "d"),
                   parnum = 1:4, value = c(1.5, 0.2, 0.8, 0))
  V <- diag(c(9, 4, 16, 1))
  colnames(V) <- c("g.2", "a1.1", "d.4", "u.3")
  expect_equal(.samplesize_reporting_se(pt, V), c(2, 0.48, 0.16, 4))
  pt$value[2:3] <- c(0, 1)
  expect_true(all(is.na(.samplesize_reporting_se(pt, V)[2:3])))
  diag(V) <- c(NA, -1, Inf, 0)
  expect_true(all(is.na(.samplesize_reporting_se(pt, V))))
  expect_true(all(is.na(.samplesize_reporting_se(pt, NULL))))
})

test_that("SE evaluation failure preserves estimates and reports exclusions", {
  original <- .samplesize_fit
  local_mocked_bindings(.samplesize_fit = function(args) {
    if (isTRUE(args$SE)) return(NULL)
    original(args)
  })
  set.seed(3)
  ss <- id_samplesize(get_fit("2PL"), params = "Item.1:a1", N = 1000, reps = 2)
  expect_true(all(is.finite(ss$estimates[["Item.1:a1"]])))
  expect_equal(ss$summary$success, 1)
  expect_equal(ss$summary$conv, 1)
  expect_equal(ss$summary$se_failure, 1)
  expect_equal(ss$summary$n_used, 0)
  expect_true(all(is.na(ss$summary[c("mean", "sd", "mirt_se")])))
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(ss))
})

test_that("an error is raised when there is nothing to study", {
  expect_error(id_samplesize(get_fit("2PL"), N = 300, reps = 1),
               "No parameters to study")
})

test_that("truth must match the model", {
  expect_error(id_samplesize(get_fit("2PL"), params = "Item.1:a1",
                             truth = get_fit("3PL"), N = 300, reps = 1),
               "same items and itemtypes")
})

test_that("removed cutoff arguments cannot be passed silently to mirt", {
  expect_error(id_samplesize(get_fit("2PL"), rmse_target = c(a1 = 0.3)),
               "Removed argument.*rmse_target")
  expect_error(id_samplesize(get_fit("2PL"), min_success = 0.9),
               "Removed argument.*min_success")
  expect_error(id_samplesize(get_fit("2PL"), SE = FALSE), "SE.*managed")
})
