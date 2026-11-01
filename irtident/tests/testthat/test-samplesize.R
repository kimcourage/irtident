test_that("sample size simulation returns the expected structure", {
  set.seed(3)
  ss <- id_samplesize(get_fit("2PL"), params = c("Item.1:a1", "Item.1:d"),
                      N = c(300, 1000), reps = 3)
  expect_s3_class(ss, "irtident_samplesize")
  expect_equal(nrow(ss$summary), 4L)
  expect_equal(nrow(ss$required), 2L)
  expect_equal(nrow(ss$estimates), 6L)
  expect_true(all(ss$summary$success >= 0 & ss$summary$success <= 1))
  expect_true(all(ss$summary$rmse >= 0, na.rm = TRUE))
  expect_output(print(ss), "Required sample size")
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

test_that("rmse targets are matched to graded intercepts", {
  out <- irtident:::.target_for(c("d2", "a1", "zz"), c(a1 = 0.3, d = 0.2))
  expect_equal(out, c(0.2, 0.3, Inf))
})
