test_that("variance decomposition is internally consistent", {
  loc <- id_local(get_fit("3PL"))
  p <- loc$params
  inner <- !p$at_bound
  expect_equal(nrow(p), 15L)
  expect_true(all(p$vif[inner] >= 1 - 1e-6))
  expect_true(all(p$vif_within[inner] >= 1 - 1e-6))
  expect_true(all(p$vif_between[inner] >= 1 - 1e-6))
  expect_equal(p$vif, p$vif_within * p$vif_between, tolerance = 1e-8)
  expect_equal(p$se^2, p$se_cond^2 * p$vif, tolerance = 1e-8)
})

test_that("3PL on LSAT7: g at bound for item 1, g-d trade-off for item 2", {
  loc <- id_local(get_fit("3PL"))
  p <- loc$params
  expect_identical(p$status[p$param == "Item.1:g"], "boundary")
  expect_identical(p$status[p$param == "Item.2:g"], "weak")
  expect_match(p$flags[p$param == "Item.2:g"], "confounded \\(within\\)")
  expect_gt(p$vif_within[p$param == "Item.2:g"], p$vif_between[p$param == "Item.2:g"])
  # the weakest direction raises item 2's g while lowering its d
  expect_true(length(loc$directions) >= 1L)
  w <- loc$directions[[1]]
  expect_true(all(c("Item.2:g", "Item.2:d") %in% w$param))
  expect_true(sign(w$loading[w$param == "Item.2:g"]) !=
                sign(w$loading[w$param == "Item.2:d"]))
  expect_lt(as.numeric(names(loc$directions)[1]), loc$cuts$eig)
})

test_that("weak direction agrees with a constrained refit", {
  fit <- get_fit("3PL")
  w <- id_local(fit)$directions[[1]]
  pv <- mirt::mod2values(fit)
  r <- pv$item == "Item.2" & pv$name == "g"
  d0 <- pv$value[pv$item == "Item.2" & pv$name == "d"]
  pv$est[r] <- FALSE
  pv$value[r] <- pv$value[r] + 0.05
  refit <- fit_quiet(lsat7, 1, itemtype = "3PL", pars = pv)
  d1 <- mirt::mod2values(refit)$value[pv$item == "Item.2" & pv$name == "d"]
  expect_identical(sign(d1 - d0),
                   sign(w$loading[w$param == "Item.2:d"] /
                          w$loading[w$param == "Item.2:g"]))
})

test_that("2PL on LSAT7 is fine everywhere", {
  loc <- id_local(get_fit("2PL"))
  expect_true(all(loc$params$status == "ok"))
  expect_true(all(loc$items$status == "ok"))
})

test_that("models fitted without SE give the same answer as with SE", {
  with_se <- fit_quiet(lsat7, 1, itemtype = "2PL", SE = TRUE)
  a <- id_local(get_fit("2PL"))$params
  b <- id_local(with_se)$params
  expect_equal(a$se, b$se, tolerance = 1e-3)
})

test_that("Rasch: latent variance is included, between-item inflation only", {
  p <- id_local(get_fit("Rasch"))$params
  expect_true("GROUP:COV_11" %in% p$param)
  expect_equal(p$vif_within[p$item != "GROUP"], rep(1, 5), tolerance = 1e-8)
})

test_that("graded models work", {
  loc <- id_local(get_fit("graded"))
  expect_equal(nrow(loc$items), 4L)
  expect_true(all(loc$params$status == "ok"))
})

test_that("print method runs", {
  expect_output(print(id_local(get_fit("3PL"))), "Weakly identified directions")
})
