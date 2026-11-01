# Small models shared by the tests (fitted once per test run)
lsat7 <- mirt::expand.table(mirt::LSAT7)

fit_quiet <- function(...) suppressWarnings(suppressMessages(mirt::mirt(..., verbose = FALSE)))

.fits <- new.env()
get_fit <- function(type) {
  if (is.null(.fits[[type]])) {
    .fits[[type]] <- switch(type,
      `3PL` = fit_quiet(lsat7, 1, itemtype = "3PL"),
      `2PL` = fit_quiet(lsat7, 1, itemtype = "2PL"),
      Rasch = fit_quiet(lsat7, 1, itemtype = "Rasch"),
      graded = fit_quiet(mirt::Science, 1, itemtype = "graded")
    )
  }
  .fits[[type]]
}
