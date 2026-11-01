test_that("non-mirt input is rejected", {
  expect_error(id_local(lm(dist ~ speed, cars)), "single-group model")
})

test_that("multidimensional models are rejected", {
  fit2 <- fit_quiet(lsat7, 2, itemtype = "2PL")
  expect_error(id_local(fit2), "unidimensional")
})

test_that("unsupported itemtypes are rejected", {
  fit <- fit_quiet(mirt::Science, 1, itemtype = "gpcm")
  expect_error(id_local(fit), "Unsupported itemtype")
})

test_that("unknown parameter labels give an informative error", {
  expect_error(id_profile(get_fit("2PL"), params = "Item.9:a1"),
               "Unknown or fixed parameter")
  expect_error(id_profile(get_fit("2PL"), params = "Item.1:g"),
               "Unknown or fixed parameter")
})

test_that("params can be given as a data frame", {
  pt <- irtident:::.par_table(get_fit("2PL"))
  out <- irtident:::.resolve_params(data.frame(item = "Item.1", name = "d"),
                                    pt, default = character())
  expect_identical(out, "Item.1:d")
})
