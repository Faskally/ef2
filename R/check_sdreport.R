# Warn when standard errors cannot be computed, naming the parameters
# affected. This is usually because a standard deviation is estimated at
# (close to) zero, for example when the data have no information on the
# penalised part of a smooth.
#   sdrep: the output of sdreport
#   labels: names for the fixed parameters, in the order of sdrep$par.fixed
check_sdreport <- function(sdrep, labels) {
  if (isTRUE(sdrep$pdHess)) {
    return(invisible(TRUE))
  }

  v <- diag(sdrep$cov.fixed)
  bad <- which(!is.finite(v) | v <= 0)
  is_sigma <- names(sdrep$par.fixed) == "log_sigma"

  msg <- "the Hessian is not positive definite, so standard errors are unreliable"
  if (length(bad)) {
    msg <- paste0(msg, ", and could not be computed for:")
    details <- vapply(bad, function(i) {
      if (is_sigma[i]) {
        sigma <- exp(sdrep$par.fixed[[i]])
        paste0(
          "\n  sigma ", labels[i], " = ", format(signif(sigma, 3)),
          if (sigma < 1e-3) {
            ", close to zero, so the part of the term it controls is shrunk to zero"
          }
        )
      } else {
        paste0("\n  ", labels[i])
      }
    }, character(1))
    msg <- paste0(msg, paste(details, collapse = ""))
    if (any(is_sigma[bad])) {
      msg <- paste0(msg, "\nconsider removing or simplifying the terms affected")
    }
  }

  warning(msg, call. = FALSE)
  invisible(FALSE)
}

# evaluate expr without the 'NaNs produced' warning from sqrt() of negative
# variances in sdreport, which check_sdreport explains
muffle_nan <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (identical(conditionMessage(w), "NaNs produced")) {
      invokeRestart("muffleWarning")
    }
  })
}
