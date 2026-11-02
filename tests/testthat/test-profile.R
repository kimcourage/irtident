test_that("profile of a well identified 2PL slope is closed and near Wald", {
  pr <- id_profile(get_fit("2PL"), params = "Item.1:a1", n_grid = 9)
  s <- pr$summary
  expect_identical(s$status, "identified")
  expect_true(s$lower < s$est && s$est < s$upper)
  expect_equal(s$lower, s$wald_lower, tolerance = 0.15)
  expect_equal(s$upper, s$wald_upper, tolerance = 0.15)
  cv <- pr$curves[["Item.1:a1"]]
  expect_equal(min(cv$deviance, na.rm = TRUE), 0, tolerance = 1e-3)
  expect_false(s$improved)
})

test_that("3PL item 2 guessing: profile interval reaches g = 0", {
  pr <- id_profile(get_fit("3PL"), params = "Item.2:g", n_grid = 9)
  s <- pr$summary
  expect_identical(s$lower_type, "bound")
  expect_identical(s$upper_type, "closed")
  expect_identical(s$status, "bounded")
  expect_equal(s$lower, 0)
  # Wald interval wrongly excludes 0
  expect_gt(s$wald_lower, 0.05)
})

test_that("the bound of g is judged with the chi-squared mixture cutoff", {
  expect_equal(irtident:::.boundary_cutoff(0.95), stats::qchisq(0.90, 1))
  row <- data.frame(param = "I:g", item = "I", name = "g", value = 0.3,
                    internal = stats::qlogis(0.3), lbound = 0, ubound = 1,
                    stringsAsFactors = FALSE)
  lrow <- data.frame(param = "I:g", se = 0.5)
  curve_with <- function(d0) {
    data.frame(value = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
               loglik = 0, deviance = c(d0, 1, 0.5, 0, 2, 5),
               point = "fitted", stringsAsFactors = FALSE)
  }
  summ <- function(d0) {
    irtident:::.summarise_profile(row, lrow, curve_with(d0),
                                  cutoff = stats::qchisq(0.95, 1), level = 0.95)
  }
  inside <- summ(2.0)   # below 2.71: g = 0 not rejected
  between <- summ(3.0)  # between 2.71 and 3.84: g = 0 rejected at the bound
  beyond <- summ(5.0)   # above 3.84: interior crossing
  expect_identical(inside$lower_type, "bound")
  expect_identical(inside$status, "bounded")
  expect_identical(between$lower_type, "closed")
  expect_equal(between$lower, 0)
  expect_identical(between$status, "identified")
  expect_identical(beyond$lower_type, "closed")
  expect_gt(beyond$lower, 0)
  expect_lt(beyond$lower, 0.1)
})

test_that("3PL boundary guessing parameter has a flat profile", {
  pr <- id_profile(get_fit("3PL"), params = "Item.5:g", n_grid = 7)
  expect_identical(pr$summary$status, "flat")
  expect_lt(pr$summary$max_dev, pr$cutoff)
})

test_that("graded intercepts can be profiled without ordering failures", {
  pr <- id_profile(get_fit("graded"), params = c("Future:d1", "Future:d3"),
                   n_grid = 7)
  expect_true(all(pr$summary$n_failed == 0))
  expect_true(all(pr$summary$status == "identified"))
})

test_that("user grids are respected", {
  g <- c(0, 0.1, 0.2, 0.4)
  pr <- id_profile(get_fit("3PL"), params = "Item.2:g",
                   grid = list(`Item.2:g` = g))
  vals <- pr$curves[["Item.2:g"]]$value
  expect_true(all(g %in% vals))
  expect_length(vals, length(g) + 1L) # the estimate is always added
})

test_that("nothing flagged means nothing profiled", {
  expect_message(pr <- id_profile(get_fit("2PL")), "nothing to profile")
  expect_equal(nrow(pr$summary), 0L)
})

test_that("print and plot methods run", {
  pr <- id_profile(get_fit("3PL"), params = "Item.2:g", n_grid = 5)
  expect_output(print(pr), "profile likelihood intervals")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(pr))
})
