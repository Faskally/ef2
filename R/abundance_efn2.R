# Abundance for efn2 fits, called by abundance()
#
# The realised abundance is the catch plus the fish that were not caught,
# N = T + M, and given the parameters M is Poisson with mean
# lambda (1 - P). The estimate is T + M_hat, with M_hat = lambda (1 - P),
# and the variance of M combines this Poisson variance with the
# uncertainty in lambda and P, from the joint covariance of the
# coefficients of the two models. Confidence intervals are computed for M
# on the log scale, so the lower limit is at least the catch.
#
# The expected abundance is lambda, with confidence intervals on the log
# scale.
abundance_efn2 <- function(object, newdata = NULL, area = NULL, by = NULL,
  target = c("realised", "expected"), level = 0.95) {
  target <- match.arg(target)
  if (!is.matrix(object$Vp_joint)) {
    stop("standard errors are not available: refit the model with hessian = TRUE")
  }

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

  # expected abundance of each sample, in the order of the samples
  sdata <- data[!duplicated(sample), , drop = FALSE]
  Xn <- efp_lpmatrix(list(Gsetup = object$Gsetup_n, G = object$G_n), sdata)
  log_lambda <- drop(Xn %*% object$gam_coefficients_n) +
    abundance_offset(object$abundance, sdata)
  lambda <- exp(log_lambda)
  Xn <- Xn[, object$gam_cols_n, drop = FALSE]

  # derivatives with respect to the coefficients of both models, the
  # capture probability coefficients first, as in Vp_joint
  if (target == "realised") {
    est <- catch + lambda * (1 - P)
    missed <- lambda * (1 - P)
    var_process <- missed
    # d M / d logit(P) = -lambda P (1 - P), d M / d log(lambda) = M
    D <- cbind(cp$J * (-lambda * P * (1 - P)), Xn * missed)
  } else {
    est <- lambda
    var_process <- rep(0, length(lambda))
    D <- cbind(cp$J * 0, Xn * lambda)
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
    out <- as.data.frame(lapply(data[by], function(x) {
      x <- tapply(x, sample, unique, simplify = FALSE)
      if (any(lengths(x) != 1)) {
        stop("the 'by' variables must be the same for all passes of a sample")
      }
      unlist(x, use.names = FALSE)
    }))
    group <- interaction(out, drop = TRUE, lex.order = TRUE)
    group <- factor(group, levels = unique(group))
    out <- out[!duplicated(group), , drop = FALSE]
    out$nsample <- as.vector(table(group))
  }

  Dg <- rowsum(as.matrix(D), group, reorder = FALSE)
  var_param <- pmax(0, rowSums((Dg %*% object$Vp_joint) * Dg))
  var_total <- as.vector(rowsum(var_process, group, reorder = FALSE)) + var_param

  out$catch <- as.vector(rowsum(catch, group, reorder = FALSE))
  if (is.null(by)) {
    out$P <- P
    out$lambda <- lambda
  }
  out$N <- as.vector(rowsum(est, group, reorder = FALSE))
  out$se <- sqrt(var_total)

  z <- qnorm(1 - (1 - level) / 2)
  if (target == "realised") {
    M <- out$N - out$catch
    cv <- out$se / M
    out$lower <- out$catch + M * exp(-z * cv)
    out$upper <- out$catch + M * exp(z * cv)
  } else {
    cv <- out$se / out$N
    out$lower <- out$N * exp(-z * cv)
    out$upper <- out$N * exp(z * cv)
  }

  if (!is.null(area)) {
    if (is.character(area) && length(area) == 1) {
      area <- data[[area]]
    }
    if (length(area) != nrow(data)) {
      stop("area must be a column name or have one value per row of the data")
    }
    area <- tapply(area, sample, unique, simplify = FALSE)
    if (any(lengths(area) != 1)) {
      stop("area must be the same for all passes of a sample")
    }
    area <- tapply(unlist(area, use.names = FALSE), group, unique, simplify = FALSE)
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

# the offset in a model formula, evaluated in data
abundance_offset <- function(formula, data) {
  tt <- stats::terms(formula)
  off <- attr(tt, "offset")
  out <- rep(0, nrow(data))
  if (length(off)) {
    vars <- attr(tt, "variables")
    for (i in off) {
      out <- out + eval(vars[[i + 1]], data, environment(formula))
    }
  }
  out
}
