test_that("sample size simulation returns the expected structure", {
  set.seed(3)
  ss <- id_samplesize(get_fit("2PL"), params = c("Item.1:a1", "Item.1:d"),
                      N = c(300, 1000), reps = 3)
  expect_s3_class(ss, "irtident_samplesize")
  expect_equal(nrow(ss$summary), 4L)
  expect_named(ss, c("summary", "estimates", "reps"))
  expect_false(any(c("target", "meets") %in% names(ss$summary)))
  expect_equal(nrow(ss$estimates), 6L)
  expect_true(all(ss$summary$success >= 0 & ss$summary$success <= 1))
  expect_true(all(ss$summary$rmse >= 0, na.rm = TRUE))
  for (i in seq_len(nrow(ss$summary))) {
    s <- ss$summary[i, ]
    z <- ss$estimates[ss$estimates$N == s$N & ss$estimates$ok, s$param]
    expect_equal(s$rmse, sqrt(mean((z - s$true)^2, na.rm = TRUE)))
    expect_equal(s$bias, mean(z, na.rm = TRUE) - s$true)
  }
  output <- paste(capture.output(print(ss)), collapse = "\n")
  expect_match(output, "rmse")
  expect_false(grepl("Required sample size|target|meets", output))
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
})
