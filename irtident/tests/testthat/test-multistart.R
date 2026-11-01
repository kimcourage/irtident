test_that("2PL: all starts agree and the original fit is best", {
  set.seed(11)
  ms <- id_multistart(get_fit("2PL"), n_starts = 4)
  expect_equal(nrow(ms$solutions), 5L)
  expect_equal(sum(ms$clusters$type == "equivalent fit"), 1L)
  expect_true(ms$original_is_best)
  expect_false(any(ms$params$unstable))
})

test_that("results do not depend on the RNG being used inside fits", {
  set.seed(5)
  a <- id_multistart(get_fit("2PL"), n_starts = 2)
  set.seed(5)
  b <- id_multistart(get_fit("2PL"), n_starts = 2)
  expect_equal(a$solutions$loglik, b$solutions$loglik)
})

test_that("reflected solutions are aligned", {
  v <- c(-1.2, 0.5, -0.8, 0.3)
  nm <- c("a1", "d", "a1", "d")
  expect_equal(irtident:::.align_reflection(v, nm), c(1.2, 0.5, 0.8, 0.3))
  expect_equal(irtident:::.align_reflection(-v, nm), -v)
})

test_that("random starts keep graded intercepts descending", {
  set.seed(1)
  pv <- irtident:::.random_start(mirt::mod2values(get_fit("graded")))
  for (it in unique(pv$item[pv$class == "graded"])) {
    d <- pv$value[pv$item == it & grepl("^d[0-9]+$", pv$name)]
    expect_true(all(diff(d) < 0))
  }
})

test_that("print and plot methods run", {
  set.seed(2)
  ms <- id_multistart(get_fit("2PL"), n_starts = 2)
  expect_output(print(ms), "random starts")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(ms))
})
