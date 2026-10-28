#' Cumulative capture probability of each sample
#'
#' The probability that a fish is caught on one of the passes of a sample,
#' \eqn{P = 1 - \prod_i (1 - p_i)}, where \eqn{p_i} is the capture
#' probability on pass \eqn{i}.
#'
#' Standard errors are computed on the logit scale by the delta method,
#' using the covariance of the fixed and random effects. The derivative of
#' \eqn{\textrm{logit}(P)} with respect to the linear predictor of pass
#' \eqn{i} is \eqn{p_i / P}. Confidence intervals are computed on the logit
#' scale and back transformed.
#'
#' @param object an efp2 fitted model
#' @param newdata an optional data.frame with one row per pass, containing
#'   the covariates in the model and the variables used for \code{pass} and
#'   \code{id} in the call to \code{\link{efp2}}. If omitted the data used
#'   to fit the model are used. An offset is only included for the data used
#'   to fit the model.
#' @param se.fit if TRUE standard errors and confidence intervals are
#'   returned
#' @param level the confidence level
#' @return a data.frame with one row per sample, in the order the samples
#'   first appear in the data, with columns \code{id}, \code{npass} the
#'   number of passes and \code{P}. If \code{se.fit = TRUE} there are also
#'   columns \code{se}, the standard error of \code{logit(P)}, and
#'   \code{lower} and \code{upper}, the confidence limits for \code{P}.
#' @references
#' Glover, R. S., Fryer, R. J., Soulsby, C., Bacon, P. J. and Malcolm, I. A.
#' (2018) Incorporating estimates of capture probability and river network
#' covariance in novel habitat - abundance models. Ecological Indicators,
#' 93, 302-315. \doi{10.1016/j.ecolind.2018.05.013}
#' @seealso \code{\link{abundance}}
#' @examples
#' ef2_data$pass12 <- factor(pmin(ef2_data$pass, 2))
#' ef2_data$sampleID <- as.numeric(factor(paste(ef2_data$date,
#'   ef2_data$siteID, ef2_data$lifestage)))
#' ef2_data <- ef2_data[order(ef2_data$sampleID, ef2_data$pass), ]
#' m1 <- efp2(count ~ lifestage + pass12, data = ef2_data,
#'   pass = pass, id = sampleID)
#' head(capture_prob(m1, se.fit = TRUE))
#' @importFrom stats qlogis qnorm
#' @export
capture_prob <- function(object, newdata = NULL, se.fit = FALSE, level = 0.95) {
  data <- if (is.null(newdata)) object$data else newdata
  id <- call_variable(object, "id", data)
  sample <- factor(id, levels = unique(id))

  eta <- predict(object, newdata = newdata)
  # log(1 - p), stable for large eta
  log_q <- plogis(eta, lower.tail = FALSE, log.p = TRUE)
  P <- -expm1(drop(rowsum(log_q, sample, reorder = FALSE)))

  out <- data.frame(
    id = unique(id),
    npass = as.vector(table(sample)),
    P = P
  )

  if (se.fit) {
    if (!is.matrix(object$Vp)) {
      stop("standard errors are not available: refit the model with hessian = TRUE")
    }
    X <- predict(object, newdata = newdata, type = "lpmatrix")
    X <- X[, object$gam_cols, drop = FALSE]
    # gradient of logit(P) for each sample: sum over passes of x_i p_i / P
    p <- plogis(eta)
    J <- rowsum(as.matrix(X) * (p / P[as.integer(sample)]), sample, reorder = FALSE)
    out$se <- sqrt(pmax(0, rowSums((J %*% object$Vp) * J)))
    z <- qnorm(1 - (1 - level) / 2)
    out$lower <- plogis(qlogis(P) - z * out$se)
    out$upper <- plogis(qlogis(P) + z * out$se)
  }

  rownames(out) <- NULL
  out
}

#' Estimate abundance and density from removal samples
#'
#' Abundance is estimated as \eqn{\hat N = T / P}, where \eqn{T} is the total
#' catch over the passes of a sample and \eqn{P} is the cumulative capture
#' probability, see \code{\link{capture_prob}}.
#'
#' The variance of \eqn{\hat N} combines the binomial sampling variation in
#' the catch, given \eqn{P}, and the uncertainty in the estimate of
#' \eqn{P} (Huggins, 1989):
#' \deqn{\textrm{var}(\hat N) \approx \hat N (1 - P) / P +
#'   \hat N^2 (1 - P)^2 \textrm{var}(\textrm{logit}(P)).}
#' Confidence intervals are computed on the log scale, and the lower limit
#' is not less than the total catch. When the total catch is zero the
#' estimated abundance is zero and the confidence limits are \code{NA}.
#'
#' @inheritParams capture_prob
#' @param newdata an optional data.frame as for \code{\link{capture_prob}},
#'   which must also contain the counts used as the response in the model
#' @param area an optional name of a column in the data giving the area
#'   fished, or a vector of areas with one value per row of the data. If
#'   given, density (abundance per unit area) is also returned.
#' @return a data.frame with one row per sample, in the order the samples
#'   first appear in the data, with columns \code{id}, \code{catch} the
#'   total catch, \code{P}, \code{N} the estimated abundance, its standard
#'   error \code{se}, and confidence limits \code{lower} and \code{upper}.
#'   If \code{area} is given there are also columns \code{area},
#'   \code{density}, \code{density_lower} and \code{density_upper}.
#' @references
#' Huggins, R. M. (1989) On the statistical analysis of capture experiments.
#' Biometrika, 76, 133-140. \doi{10.1093/biomet/76.1.133}
#' @seealso \code{\link{capture_prob}}
#' @examples
#' ef2_data$pass12 <- factor(pmin(ef2_data$pass, 2))
#' ef2_data$sampleID <- as.numeric(factor(paste(ef2_data$date,
#'   ef2_data$siteID, ef2_data$lifestage)))
#' ef2_data <- ef2_data[order(ef2_data$sampleID, ef2_data$pass), ]
#' m1 <- efp2(count ~ lifestage + pass12, data = ef2_data,
#'   pass = pass, id = sampleID)
#' head(abundance(m1, area = "area"))
#' @export
abundance <- function(object, newdata = NULL, area = NULL, level = 0.95) {
  data <- if (is.null(newdata)) object$data else newdata
  id <- call_variable(object, "id", data)
  sample <- factor(id, levels = unique(id))

  y <- tryCatch(
    eval(object$formula[[2]], data, environment(object$formula)),
    error = function(e) NULL
  )
  if (is.null(y) || length(y) != nrow(data)) {
    stop("the data must contain the counts used as the response in the model")
  }

  out <- capture_prob(object, newdata = newdata, se.fit = TRUE, level = level)
  out$npass <- NULL
  out$catch <- as.vector(rowsum(y, sample, reorder = FALSE))
  out <- out[c("id", "catch", "P", "se", "lower", "upper")]
  se_logitP <- out$se
  out$lower <- out$upper <- NULL

  P <- out$P
  N <- out$catch / P
  var_N <- N * (1 - P) / P + N^2 * (1 - P)^2 * se_logitP^2
  out$N <- N
  out$se <- sqrt(var_N)

  z <- qnorm(1 - (1 - level) / 2)
  cv <- out$se / N
  out$lower <- pmax(out$catch, N * exp(-z * cv))
  out$upper <- N * exp(z * cv)
  zero <- out$catch == 0
  out$lower[zero] <- out$upper[zero] <- NA

  if (!is.null(area)) {
    if (is.character(area) && length(area) == 1) {
      area <- data[[area]]
    }
    if (length(area) != nrow(data)) {
      stop("area must be a column name or have one value per row of the data")
    }
    a <- tapply(area, sample, unique, simplify = FALSE)
    if (any(lengths(a) != 1)) {
      stop("area must be the same for all passes of a sample")
    }
    out$area <- unlist(a)
    out$density <- out$N / out$area
    out$density_lower <- out$lower / out$area
    out$density_upper <- out$upper / out$area
  }

  rownames(out) <- NULL
  out
}

# evaluate the pass or id argument of the call to efp2 in data
call_variable <- function(object, name, data) {
  expr <- object$call[[name]]
  if (is.null(expr)) expr <- as.name(name)
  x <- tryCatch(eval(expr, data, environment(object$formula)),
    error = function(e) NULL
  )
  if (is.null(x) || length(x) != nrow(data)) {
    stop("the data must contain the variable '", deparse(expr), "' used for ",
      name, " in the call to efp2")
  }
  x
}
