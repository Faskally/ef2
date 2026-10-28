fake_sdrep <- function(pdHess, v, log_sigma = c(-1, -10)) {
  par <- c(beta = 0.5, beta = 1, log_sigma = log_sigma[1], log_sigma = log_sigma[2])
  list(pdHess = pdHess, par.fixed = par, cov.fixed = diag(v, length(par)))
}
labels <- c("(Intercept)", "x", "s(year)", "s(site)")

test_that("no warning when the Hessian is positive definite", {
  expect_no_warning(check_sdreport(fake_sdrep(TRUE, 1), labels))
  expect_true(check_sdreport(fake_sdrep(TRUE, 1), labels))
})

test_that("warning names a standard deviation estimated close to zero", {
  sdrep <- fake_sdrep(FALSE, c(0.1, 0.1, 0.2, -5))
  expect_warning(
    check_sdreport(sdrep, labels),
    "sigma s\\(site\\) = 4.54e-05, close to zero"
  )
  expect_warning(check_sdreport(sdrep, labels), "removing or simplifying")
  # parameters with valid variances are not listed
  expect_false(grepl("s\\(year\\)", tryCatch(
    check_sdreport(sdrep, labels),
    warning = conditionMessage
  )))
})

test_that("warning names fixed effects without a valid variance", {
  sdrep <- fake_sdrep(FALSE, c(0.1, NaN, 0.2, 0.3))
  msg <- tryCatch(check_sdreport(sdrep, labels), warning = conditionMessage)
  expect_match(msg, "could not be computed for:\n  x$")
  expect_false(grepl("simplifying", msg))
})

test_that("generic warning when no variance is invalid", {
  sdrep <- fake_sdrep(FALSE, 0.1)
  msg <- tryCatch(check_sdreport(sdrep, labels), warning = conditionMessage)
  expect_equal(msg, "the Hessian is not positive definite, so standard errors are unreliable")
})

test_that("muffle_nan only removes the NaNs produced warning", {
  expect_no_warning(muffle_nan(sqrt(-1)))
  expect_warning(muffle_nan(warning("other")), "other")
})
