#' Summary of an efn2 fit
#'
#' @param object,x an efn2 fitted model, or its summary
#' @param digits the number of significant digits to print
#' @param ... further arguments, ignored
#' @return \code{summary} returns an object of class \code{summary.efn2},
#'   with the coefficient tables for capture probability
#'   (\code{coefficients}) and log abundance (\code{coefficients_n}), and the
#'   standard deviations of the random effects (\code{sigma} and
#'   \code{sigma_n}).
#' @name efn2-methods
#' @rdname efn2-methods
#' @method summary efn2
#' @exportS3Method summary efn2
summary.efn2 <- function(object, ...) {
  out <- object[c("call", "formula", "abundance", "aic", "convergence")]
  out$logLik <- logLik(object)

  # NaN standard errors are explained by the warning from efn2
  fx <- muffle_nan(summary(object$sdrep, "fixed", p.value = TRUE))
  # for the cut model, capture probability is from stage 1, and the
  # standard errors for abundance allow for its uncertainty
  fx_capture <- fx
  re_capture <- NULL
  if (isTRUE(object$cut)) {
    fx_capture <- muffle_nan(summary(object$sdrep_capture, "fixed", p.value = TRUE))
    re_capture <- muffle_nan(summary(object$sdrep_capture, "report", p.value = FALSE))
  }
  out$coefficients <- fx_capture[rownames(fx_capture) == "beta", , drop = FALSE]
  rownames(out$coefficients) <- names(object$coefficients)
  out$coefficients_n <- fx[rownames(fx) == "beta_n", , drop = FALSE]
  rownames(out$coefficients_n) <- names(object$coefficients_n)
  if (isTRUE(object$cut) && is.matrix(object$Vb_n)) {
    # standard errors allowing for the uncertainty in capture probability
    se <- sqrt(diag(object$Vb_n))
    z <- out$coefficients_n[, 1] / se
    out$coefficients_n[, 2] <- se
    out$coefficients_n[, 3] <- z
    out$coefficients_n[, 4] <- 2 * stats::pnorm(-abs(z))
  }

  re <- muffle_nan(summary(object$sdrep, "report", p.value = FALSE))
  if (length(object$sigma)) {
    out$sigma <- (if (is.null(re_capture)) re else re_capture)
    out$sigma <- out$sigma[rownames(out$sigma) == "sigma", , drop = FALSE]
    rownames(out$sigma) <- paste("sigma", object$re_labels)
  }
  if (length(object$sigma_n)) {
    out$sigma_n <- re[rownames(re) == "sigma_n", , drop = FALSE]
    rownames(out$sigma_n) <- paste("sigma", object$re_labels_n)
  }

  if (identical(object$distribution, "negbin")) {
    out$theta <- re[rownames(re) == "theta", , drop = FALSE]
    rownames(out$theta) <- "theta"
  }
  out$distribution <- object$distribution
  out$cut <- isTRUE(object$cut)
  out$nsample <- length(object$catch)
  class(out) <- "summary.efn2"
  out
}

#' @rdname efn2-methods
#' @method print summary.efn2
#' @exportS3Method print summary.efn2
print.summary.efn2 <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n", sep = "")

  section <- function(title, formula, sigma, coefficients) {
    cat(title, "\n", sep = "")
    print(formula)
    cat("\n")
    if (!is.null(sigma)) {
      cat("Random effects:\n")
      print(sigma, digits = digits, ...)
      cat("\n")
    }
    if (nrow(coefficients)) {
      cat("Coefficients:\n")
      print(coefficients, digits = digits, ...)
      cat("\n")
    }
  }
  # [[ rather than $, which would partially match sigma to sigma_n
  section("Capture probability (logit scale):", x[["formula"]], x[["sigma"]], x[["coefficients"]])
  section("Abundance (log scale):", x[["abundance"]], x[["sigma_n"]], x[["coefficients_n"]])
  if (!is.null(x[["theta"]])) {
    cat("Negative binomial size:\n")
    print(x[["theta"]], digits = digits, ...)
    cat("\n")
  }

  cat("Samples:", x$nsample, "\n")
  if (isTRUE(x$cut)) {
    cat("Cut model: capture probability from the conditional likelihood, not\n")
    cat("influenced by the model for abundance; no AIC\n")
  } else {
    cat("Log likelihood:", format(signif(as.numeric(x$logLik), digits)),
      " AIC:", format(signif(x$aic, digits)), "\n")
    cat("(not comparable with efp2, which uses the conditional likelihood)\n")
  }
  invisible(x)
}

#' @rdname efn2-methods
#' @method print efn2
#' @exportS3Method print efn2
print.efn2 <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print(summary(x), digits = digits, ...)
  invisible(x)
}
