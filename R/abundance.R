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
  cp <- capture_prob_parts(object, newdata, gradient = se.fit)

  out <- data.frame(
    id = unique(cp$id),
    npass = as.vector(table(cp$sample)),
    P = cp$P
  )

  if (se.fit) {
    out$se <- sqrt(pmax(0, rowSums((cp$J %*% object$Vp) * cp$J)))
    z <- qnorm(1 - (1 - level) / 2)
    out$lower <- plogis(qlogis(cp$P) - z * out$se)
    out$upper <- plogis(qlogis(cp$P) + z * out$se)
  }

  rownames(out) <- NULL
  out
}

# the cumulative capture probability of each sample, and optionally the
# gradient of logit(P) with respect to the coefficients in Vp
capture_prob_parts <- function(object, newdata = NULL, gradient = FALSE) {
  data <- if (is.null(newdata)) object$data else newdata
  id <- call_variable(object, "id", data)
  sample <- factor(id, levels = unique(id))

  eta <- predict(object, newdata = newdata)
  # log(1 - p), stable for large eta
  log_q <- plogis(eta, lower.tail = FALSE, log.p = TRUE)
  P <- -expm1(drop(rowsum(log_q, sample, reorder = FALSE)))

  J <- NULL
  if (gradient) {
    if (!is.matrix(object$Vp)) {
      stop("standard errors are not available: refit the model with hessian = TRUE")
    }
    X <- predict(object, newdata = newdata, type = "lpmatrix")
    X <- X[, object$gam_cols, drop = FALSE]
    # gradient of logit(P) for each sample: sum over passes of x_i p_i / P
    p <- plogis(eta)
    J <- rowsum(as.matrix(X) * (p / P[as.integer(sample)]), sample, reorder = FALSE)
  }

  list(data = data, id = id, sample = sample, P = P, J = J)
}

#' Estimate abundance and density from removal samples
#'
#' Abundance is estimated as \eqn{\hat N = T / P}, where \eqn{T} is the total
#' catch over the passes of a sample and \eqn{P} is the cumulative capture
#' probability, see \code{\link{capture_prob}}.
#'
#' Two targets are available, which have the same estimate but different
#' standard errors:
#' \describe{
#'   \item{\code{"realised"}}{the number of fish \eqn{N} present when the
#'     sample was taken. Given \eqn{N}, the catch is binomial, \eqn{T \sim
#'     \textrm{Bin}(N, P)}, and (Huggins, 1989)
#'     \deqn{\textrm{var}(\hat N) \approx \hat N (1 - P) / P +
#'       \hat N^2 (1 - P)^2 \textrm{var}(\textrm{logit}(P)).}
#'     The lower confidence limit is not less than the catch.}
#'   \item{\code{"expected"}}{the expected number of fish \eqn{\lambda},
#'     for example density times area, when the number present is
#'     \eqn{N \sim \textrm{Poisson}(\lambda)}. The catch is then
#'     \eqn{T \sim \textrm{Poisson}(\lambda P)}, and the first term of the
#'     variance is \eqn{\hat\lambda / P}, which also includes the natural
#'     variation of \eqn{N} about \eqn{\lambda}. The standard error is
#'     larger by a factor of about \eqn{1 / \sqrt{1 - P}} when \eqn{P} is
#'     well estimated.}
#' }
#' Use \code{"realised"} to say how many fish were present, and
#' \code{"expected"} to describe the underlying abundance, for example
#' when comparing sites or years. A Poisson model for the counts with
#' offset \eqn{\log(A P)}, where \eqn{A} is the area, is a model for the
#' expected abundance.
#'
#' Confidence intervals are computed on the log scale. When the total catch
#' is zero the estimated abundance is zero and the confidence limits are
#' \code{NA}.
#'
#' With \code{by}, abundance is summed over the samples in each group, for
#' example over the lifestages caught at a site in a year. The variation in
#' the catches is independent between samples, but the estimates of \eqn{P}
#' are correlated, as they depend on the same model coefficients, and this
#' is allowed for:
#' \deqn{\textrm{var}(\sum_s \hat N_s) \approx \sum_s v_s + w^\top \Sigma w,}
#' where \eqn{v_s} is \eqn{\hat N_s (1 - P_s) / P_s} for the realised
#' abundance and \eqn{\hat N_s / P_s} for the expected abundance,
#' \eqn{w_s = \hat N_s (1 - P_s)} and \eqn{\Sigma} is the covariance
#' matrix of \eqn{\textrm{logit}(P_s)}.
#'
#' @inheritParams capture_prob
#' @param newdata an optional data.frame as for \code{\link{capture_prob}},
#'   which must also contain the counts used as the response in the model
#' @param area an optional name of a column in the data giving the area
#'   fished, or a vector of areas with one value per row of the data. If
#'   given, density (abundance per unit area) is also returned. With
#'   \code{by}, the area must be the same for all samples in a group.
#' @param target the abundance to estimate, \code{"realised"} (the default)
#'   or \code{"expected"}, see Details
#' @param by an optional character vector of the names of variables in the
#'   data, which must be the same for all passes of a sample. Abundance is
#'   summed over the samples with each combination of their values.
#' @return a data.frame with one row per sample, in the order the samples
#'   first appear in the data, with columns \code{id}, \code{catch} the
#'   total catch, \code{P}, \code{N} the estimated abundance, its standard
#'   error \code{se}, and confidence limits \code{lower} and \code{upper}.
#'   If \code{area} is given there are also columns \code{area},
#'   \code{density}, \code{density_lower} and \code{density_upper}. With
#'   \code{by} there is one row per group, sorted by the \code{by}
#'   variables, with the \code{by} variables and \code{nsample}, the number
#'   of samples, in place of \code{id} and \code{P}.
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
#'
#' # summed over lifestages for each site and year
#' head(abundance(m1, area = "area", by = c("siteID", "year")))
#'
#' # expected abundance, for comparing years
#' abundance(m1, by = "year", target = "expected")
#' @export
abundance <- function(object, newdata = NULL, area = NULL, by = NULL,
  target = c("realised", "expected"), level = 0.95) {
  target <- match.arg(target)
  cp <- capture_prob_parts(object, newdata, gradient = TRUE)
  data <- cp$data
  sample <- cp$sample
  P <- cp$P

  y <- tryCatch(
    eval(object$formula[[2]], data, environment(object$formula)),
    error = function(e) NULL
  )
  if (is.null(y) || length(y) != nrow(data)) {
    stop("the data must contain the counts used as the response in the model")
  }
  catch <- as.vector(rowsum(y, sample, reorder = FALSE))
  N <- catch / P

  # one value per sample of a variable with one value per row of the data
  per_sample <- function(x, what) {
    x <- tapply(x, sample, unique, simplify = FALSE)
    if (any(lengths(x) != 1)) {
      stop(what, " must be the same for all passes of a sample")
    }
    unlist(x, use.names = FALSE)
  }

  # groups of samples
  if (is.null(by)) {
    out <- data.frame(id = unique(cp$id))
    group <- factor(seq_along(P))
  } else {
    missing <- setdiff(by, names(data))
    if (length(missing)) {
      stop("the data must contain the 'by' variables: ", paste(missing, collapse = ", "))
    }
    out <- as.data.frame(lapply(data[by], per_sample, what = "the 'by' variables"))
    group <- interaction(out, drop = TRUE, lex.order = TRUE)
    group <- factor(group, levels = unique(group))
    out <- out[!duplicated(group), , drop = FALSE]
    out$nsample <- as.vector(table(group))
  }

  # variance: variation in the catch of each sample given P, which is
  # independent between samples, and the uncertainty in P, which is
  # correlated between samples. The catch is binomial given the realised
  # abundance, and Poisson given the expected abundance. The derivative of
  # N with respect to logit(P) is -N (1 - P)
  var_count <- if (target == "realised") N * (1 - P) / P else N / P
  var_count <- rowsum(var_count, group, reorder = FALSE)
  G <- rowsum(cp$J * (N * (1 - P)), group, reorder = FALSE)
  var_P <- rowSums((G %*% object$Vp) * G)

  out$catch <- as.vector(rowsum(catch, group, reorder = FALSE))
  if (is.null(by)) out$P <- P
  out$N <- as.vector(rowsum(N, group, reorder = FALSE))
  out$se <- sqrt(as.vector(var_count) + pmax(0, var_P))

  z <- qnorm(1 - (1 - level) / 2)
  cv <- out$se / out$N
  out$lower <- out$N * exp(-z * cv)
  if (target == "realised") {
    # there are at least as many fish as were caught
    out$lower <- pmax(out$catch, out$lower)
  }
  out$upper <- out$N * exp(z * cv)
  zero <- out$catch == 0
  out$lower[zero] <- out$upper[zero] <- NA

  if (!is.null(area)) {
    if (is.character(area) && length(area) == 1) {
      area <- data[[area]]
    }
    if (length(area) != nrow(data)) {
      stop("area must be a column name or have one value per row of the data")
    }
    area <- tapply(per_sample(area, "area"), group, unique, simplify = FALSE)
    if (any(lengths(area) != 1)) {
      stop("area must be the same for all samples in a group: density is only ",
        "given when each group of samples is from a single area fished")
    }
    out$area <- unlist(area, use.names = FALSE)
    out$density <- out$N / out$area
    out$density_lower <- out$lower / out$area
    out$density_upper <- out$upper / out$area
  }

  if (!is.null(by)) {
    out <- out[do.call(order, unname(out[by])), , drop = FALSE]
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
