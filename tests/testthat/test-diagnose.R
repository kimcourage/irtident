test_that("3PL on LSAT7: item verdicts match the known structure", {
  skip_on_cran()
  set.seed(1)
  dx <- id_diagnose(get_fit("3PL"), n_grid = 7, n_starts = 4)
  expect_s3_class(dx, "irtident")
  v <- stats::setNames(dx$items$verdict, dx$items$item)
  expect_identical(unname(v["Item.1"]), "boundary")
  expect_true(v["Item.2"] %in% c("weak", "not identified"))
  adv <- dx$params$advice[dx$params$param == "Item.1:g"]
  expect_match(adv, "fix g = 0")
  expect_output(print(dx), "Advice")
})

test_that("components can be switched off", {
  dx <- id_diagnose(get_fit("2PL"), profile = FALSE, multistart = FALSE)
  expect_null(dx$profile)
  expect_null(dx$multistart)
  expect_true(all(dx$items$verdict == "ok"))
  expect_message(plot(dx), "No profile")
})
