#' Methods for efp2 fits
#'
#' Methods to summarise and extract information from a model fitted
#' using \code{\link{efp2}}.
#'
#' @param object,x an efp2 fitted model
#' @param type the type of fitted value: \code{"p"} for the capture
#'   probability in each pass, \code{"pi"} for the probability that a fish
#'   is first caught in each pass, \code{"both"} for a list containing both
#'   and \code{"lpmatrix"} for the design matrix.
#' @param parm the parameters to compute confidence intervals for, either a
#'   vector of numbers or names. The default is all fixed effects.
#' @param level the confidence level required
#' @param digits the number of significant digits to print
#' @param ... further arguments, ignored unless stated
#'
#' @details
#' The variance matrix returned by \code{vcov} is for the fixed effects
#' only, the joint variance matrix of the fixed and random effects is
#' available as \code{object$Vp}. Confidence intervals are Wald intervals
#' based on this variance matrix.
#'
#' The log likelihood is the marginal log likelihood when the model
#' contains random effects or penalised smooths, and the degrees of freedom
#' are the number of fixed effects plus the number of standard deviations.
#'
#' Residuals, analysis of deviance tables and diagnostic plots are not
#' available for efp2 fits, and the corresponding methods return an error.
#'
#' @return \code{fitted} returns a vector, list or matrix depending on
#'   \code{type}; \code{model.matrix} a matrix; \code{vcov} a matrix;
#'   \code{confint} a matrix with columns giving the lower and upper
#'   confidence limits; \code{logLik} an object of class \code{logLik};
#'   \code{formula} the model formula; \code{summary} an object of class
#'   \code{summary.efp2}, a list containing the coefficient table
#'   (\code{coefficients}) and the table of random effect standard
#'   deviations (\code{sigma}). The \code{print} methods return their
#'   argument invisibly.
#'
#' @examples
#' data(ef2_data)
#' ef2_data$sampleID <- factor(paste(ef2_data$date, ef2_data$siteID,
#'                                  ef2_data$lifestage, ef2_data$species))
#' m1 <- efp2(count ~ lifestage, data = ef2_data, pass = pass, id = sampleID)
#'
#' summary(m1)
#' vcov(m1)
#' confint(m1)
#' logLik(m1)
#' head(fitted(m1))
#'
#' @name efp2-methods
NULL

#' @rdname efp2-methods
#' @export
vcov.efp2 <- function(object, ...) {
  if (!is.matrix(object$Vb)) {
    stop("the variance matrix is not available: refit the model with hessian = TRUE")
  }
  object$Vb
}

#' @rdname efp2-methods
#' @importFrom stats confint.default logLik
#' @export
confint.efp2 <- function(object, parm, level = 0.95, ...) {
  confint.default(object, parm, level = level, ...)
}

#' @rdname efp2-methods
#' @export
logLik.efp2 <- function(object, ...) {
  structure(
    object$llik,
    df = object$rank, nobs = nrow(object$G), class = "logLik"
  )
}

#' @rdname efp2-methods
#' @export
formula.efp2 <- function(x, ...) {
  x$formula
}

#' @rdname efp2-methods
#' @method summary efp2
#' @exportS3Method summary efp2
summary.efp2 <- function(object, ...) {
  out <- object[c("call", "formula", "aic", "df.null", "df.residual", "convergence")]
  out$logLik <- logLik(object)

  # NaN standard errors are explained by the warning from efp2
  fx <- muffle_nan(summary(object$sdrep, "fixed", p.value = TRUE))
  fx <- fx[rownames(fx) == "beta", , drop = FALSE]
  rownames(fx) <- names(object$coefficients)
  out$coefficients <- fx

  if (length(object$sigma)) {
    re <- muffle_nan(summary(object$sdrep, "report", p.value = FALSE))
    re <- re[rownames(re) == "sigma", , drop = FALSE]
    rownames(re) <- paste("sigma", object$re_labels)
    out$sigma <- re
  }

  class(out) <- "summary.efp2"
  out
}

#' @rdname efp2-methods
#' @method print summary.efp2
#' @exportS3Method print summary.efp2
print.summary.efp2 <- function(x, digits = max(3L, getOption("digits") - 3L),  ...) {

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  cat("Formula:\n")
  print(x$formula)
  cat("\n")

  if (!is.null(x$sigma)) {
    cat("Random effects:\n")
    print(x$sigma, digits = digits, ...)
    cat("\n")
  }
  else {
    cat("No random effects\n\n")
  }

  if (nrow(x$coefficients)) {
    cat("Coefficients:\n")
    print(x$coefficients, digits = digits, ...)
  }
  else {
    cat("No coefficients\n\n")
  }

  cat(
    "\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
    x$df.residual, "Residual\n"
  )

  cat("AIC:", format(signif(x$aic, digits)))

  cat("\n")
  invisible(x)
}

#' @rdname efp2-methods
#' @export
residuals.efp2 <- function(object, ...) {
  stop("residuals are not available for efp2 fits")
}

#' @rdname efp2-methods
#' @export
anova.efp2 <- function(object, ...) {
  stop("anova is not available for efp2 fits, models can be compared using AIC")
}

#' @rdname efp2-methods
#' @export
plot.efp2 <- function(x, ...) {
  stop("plot is not available for efp2 fits, see as.gam2 for plotting model terms")
}
