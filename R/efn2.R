#' Fit a joint model for capture probability and abundance
#'
#' Fits a model in which the number of fish in each sample is Poisson, with
#' its log mean (the expected abundance) given by covariates, smooths and
#' random effects, and capture probability is modelled as in
#' \code{\link{efp2}}. This is a multinomial N-mixture, or hierarchical
#' removal, model.
#'
#' The number of fish in sample \eqn{s} is \eqn{N_s \sim
#' \textrm{Poisson}(\lambda_s)}, and the catches on the passes are then
#' independent Poisson with means \eqn{\lambda_s \pi_{si}}, where
#' \eqn{\pi_{si}} is the probability that a fish is first caught on pass
#' \eqn{i}. The likelihood is the product of the conditional likelihood of
#' the catches given the total catch, used by \code{\link{efp2}}, and a
#' Poisson likelihood for the total catch \eqn{T_s} with mean \eqn{\lambda_s
#' P_s}, where \eqn{P_s} is the cumulative capture probability.
#'
#' In \code{efp2} the total catch carries no information about capture
#' probability: \code{efp2} is equivalent to this model with a separate
#' \eqn{\lambda_s} for every sample. When the model for \eqn{\lambda_s} has
#' structure, the total catches also contribute to the estimates of capture
#' probability.
#'
#' The realised abundance of each sample is estimated as the catch plus the
#' expected number of fish that were not caught, given the catch. For a
#' Poisson number of fish this is \eqn{T_s + \hat\lambda_s (1 - \hat P_s)}.
#'
#' With \code{distribution = "negbin"} the number of fish is negative
#' binomial with mean \eqn{\lambda_s} and size \eqn{\theta}, which allows for
#' more variation between samples than can be explained by the model for
#' abundance. The total catch is then negative binomial with mean
#' \eqn{\lambda_s P_s} and the same size, and the realised abundance is
#' estimated as \eqn{T_s + \hat\lambda_s (1 - \hat P_s) (\hat\theta + T_s) /
#' (\hat\theta + \hat\lambda_s \hat P_s)}. As \eqn{\theta} tends to zero this
#' tends to \eqn{T_s / \hat P_s}, the estimate from the conditional
#' likelihood, and as \eqn{\theta} tends to infinity to the Poisson estimate.
#'
#' The log likelihood includes the likelihood of the total catches, so it
#' cannot be compared with that of a model fitted by \code{efp2}.
#'
#' @inheritParams efp2
#' @param formula a formula for logit capture probability, with the counts
#'   on each pass as the response, as for \code{\link{efp2}}
#' @param abundance a one sided formula for log expected abundance, for
#'   example \code{~ lifestage + s(siteID, bs = "re") + offset(log(area))}.
#'   The variables must be the same for all passes of a sample. Smooth terms
#'   and random effects can be used as in \code{formula}.
#' @param distribution the distribution of the number of fish in a sample,
#'   \code{"poisson"} or \code{"negbin"}, the negative binomial with mean
#'   \eqn{\lambda_s} and size \eqn{\theta}, which is estimated
#' @param offset an optional offset for the linear predictor of capture
#'   probability. An offset for log abundance is given in \code{abundance}.
#' @return If \code{fit = TRUE}, an object of class \code{efn2}, which
#'   inherits from \code{efp2}, so the methods for \code{efp2} fits apply
#'   to the model for capture probability. In addition to the components of
#'   an \code{efp2} fit it contains
#'   \item{coefficients_n}{the fixed effects for log abundance}
#'   \item{sigma_n}{the standard deviations of the random effects for log
#'     abundance}
#'   \item{lambda}{the expected abundance of each sample}
#'   \item{catch}{the total catch of each sample}
#'   \item{N}{the estimated realised abundance of each sample}
#'   \item{theta}{for the negative binomial, the estimated size parameter}
#'   \item{Vp_joint}{the covariance matrix of the coefficients of the models
#'     for capture probability and log abundance}
#'   If \code{fit = FALSE}, the data that would be passed to the optimiser.
#' @references
#' Dorazio, R. M., Jelks, H. L. and Jordan, F. (2005) Improving removal-based
#' estimates of abundance by sampling a population of spatially distinct
#' subpopulations. Biometrics, 61, 1093-1101.
#' \doi{10.1111/j.1541-0420.2005.00360.x}
#'
#' Royle, J. A. (2004) Generalized estimators of avian abundance from count
#' survey data. Animal Biodiversity and Conservation, 27, 375-386.
#' \doi{10.32800/abc.2004.27.0375}
#' @seealso \code{\link{efp2}}
#' @examples
#' ef2_data$pass12 <- factor(pmin(ef2_data$pass, 2))
#' ef2_data$sampleID <- as.numeric(factor(paste(ef2_data$date,
#'   ef2_data$siteID, ef2_data$lifestage)))
#' ef2_data <- ef2_data[order(ef2_data$sampleID, ef2_data$pass), ]
#' m1 <- efn2(count ~ lifestage + pass12,
#'   abundance = ~ lifestage + siteID + offset(log(area)),
#'   data = ef2_data, pass = pass, id = sampleID)
#' m1
#' @export
efn2 <- function(formula, abundance, data = NULL, pass = pass, id = id,
  distribution = c("poisson", "negbin"),
  offset = NULL, verbose = FALSE, hessian = TRUE,
  fit = TRUE, control = list(), sparse = TRUE) {

  distribution <- match.arg(distribution)

  if (is.null(data)) stop("must supply data")
  if (!inherits(abundance, "formula") || length(abundance) != 2) {
    stop("abundance must be a one sided formula, for example ~ lifestage")
  }

  pass <- eval(substitute(pass), data, environment(formula))
  id <- eval(substitute(id), data, environment(formula))
  if (any(tapply(pass, id, is.unsorted))) {
    stop("data must be sorted by pass")
  }
  if (is.null(offset)) {
    offset <- rep(0, nrow(data))
  }

  # model for capture probability, one row per pass
  ms <- model_setup(formula, data, sparse = sparse)
  sample_id <- factor(id)
  efdat <- c(list(y = ms$Gsetup$y, sample_id = sample_id), model_data(ms, offset))

  # model for log abundance, one row per sample in the order of the levels
  # of sample_id
  vars <- intersect(all.vars(abundance), names(data))
  varies <- vapply(vars, function(v) {
    any(tapply(data[[v]], sample_id, function(x) length(unique(x))) > 1)
  }, TRUE)
  if (any(varies)) {
    stop("the variables in abundance must be the same for all passes of a ",
      "sample: ", paste(vars[varies], collapse = ", "))
  }
  first <- match(levels(sample_id), as.character(sample_id))
  sdata <- data[first, , drop = FALSE]
  catch <- as.vector(tapply(efdat$y, sample_id, sum))
  sdata$.ef_catch <- catch
  nformula <- stats::update(abundance, .ef_catch ~ .)
  msn <- model_setup(nformula, sdata, sparse = sparse)
  offset_n <- msn$Gsetup$offset
  if (length(offset_n) != nrow(sdata)) offset_n <- rep(0, nrow(sdata))
  ndat <- c(list(y = catch, distribution = distribution), model_data(msn, offset_n))

  if (!fit) {
    return(list(efdat = efdat, ndat = ndat))
  }

  # start log abundance near the mean catch
  beta_n <- rep(0, ncol(ndat$X))
  intercept <- colnames(ndat$X) == "(Intercept)"
  start <- log(mean(catch) + 0.5) - mean(offset_n)
  if (any(intercept)) beta_n[intercept] <- start else beta_n[] <- start

  params <- list(
    alpha = rep(0, ncol(efdat$Z)),
    beta = rep(0, ncol(efdat$X)),
    log_sigma = rep(0, length(ms$re_labels)),
    alpha_n = rep(0, ncol(ndat$Z)),
    beta_n = beta_n,
    log_sigma_n = rep(0, length(msn$re_labels))
  )
  if (distribution == "negbin") {
    params$log_theta <- log(10)
  }
  rand <- c("alpha", "alpha_n")[c(length(ms$blocks), length(msn$blocks)) > 0]
  if (!length(rand)) rand <- NULL

  obj <- MakeADFun(
    function(params) ef_nll_joint(params, efdat, ndat),
    params,
    random = rand,
    silent = !verbose
  )

  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
  if (opt$convergence != 0) {
    warning("model has not converged: ", opt$message)
  }

  out <- list()
  out$llik <- -1 * as.numeric(obj$fn(opt$par))
  par_full <- obj$env$last.par
  out$rep <- obj$report(par_full)

  # capture probability, in the order of the data
  ord_output <- order(id, pass)
  out$p <- unname(unlist(out$rep$ps))[order(ord_output)]
  out$pi <- as.vector(
    ave(out$p, sample_id, FUN = function(x) x * c(1, cumprod(1 - x)[-length(x)]))
  )
  out$P <- setNames_(out$rep$P, levels(sample_id))
  out$lambda <- setNames_(out$rep$lambda, levels(sample_id))
  out$catch <- setNames_(catch, levels(sample_id))
  out$distribution <- distribution
  if (distribution == "negbin") {
    out$theta <- exp(opt$par[["log_theta"]])
  }
  out$N <- out$catch + expected_missed(out$catch, out$lambda, out$P, out$theta)

  # if the model for abundance cannot describe the variation in the total
  # catches, the fit can explain it by capture probability instead, with
  # capture probability close to zero and abundance very large
  if (any(out$P < 1e-3)) {
    warning(
      "cumulative capture probability is estimated close to zero for ",
      sum(out$P < 1e-3), " sample(s), so the estimates of abundance are not ",
      "reliable. The model for abundance may be too simple to describe the ",
      "variation in the total catches: consider adding a random effect for ",
      "samples, s(f, bs = \"re\"), to the abundance formula.",
      call. = FALSE
    )
  }

  labels <- c(
    colnames(ms$X)[ms$fx_keep], ms$re_labels,
    colnames(msn$X)[msn$fx_keep], msn$re_labels,
    if (distribution == "negbin") "log(theta)"
  )
  out$sdrep <- muffle_nan(
    sdreport(obj, getJointPrecision = hessian && !is.null(rand))
  )
  check_sdreport(out$sdrep, labels)

  out$data <- data
  out$offset <- offset
  out$formula <- formula
  out$abundance <- abundance
  out$terms <- ms$Gsetup$terms
  out$call <- match.call()
  out$rank <- length(opt$par)
  out$aic <- -2 * out$llik + 2 * out$rank

  # capture probability, with the same components as an efp2 fit
  out$G <- ms$G
  out$Gfit <- ms$Gfit
  out$fx_cols <- ms$fx_cols
  out$re_cols <- ms$re_cols
  out$gam_cols <- ms$gam_cols
  out$re_labels <- ms$re_labels
  out$coefficients <- opt$par[names(opt$par) == "beta"]
  names(out$coefficients) <- colnames(ms$X)[ms$fx_keep]
  out$random <- par_full[names(par_full) == "alpha"]
  names(out$random) <- colnames(ms$Z)
  out$sigma <- exp(opt$par[names(opt$par) == "log_sigma"])
  names(out$sigma) <- ms$re_labels
  out$gam_coefficients <- drop(as.matrix(ms$to_gam %*% c(out$coefficients, out$random)))
  names(out$gam_coefficients) <- colnames(ms$G)
  out$Gsetup <- ms$Gsetup

  # log abundance
  out$data_n <- sdata
  out$G_n <- msn$G
  out$gam_cols_n <- msn$gam_cols
  out$re_labels_n <- msn$re_labels
  out$offset_n <- offset_n
  out$coefficients_n <- opt$par[names(opt$par) == "beta_n"]
  names(out$coefficients_n) <- colnames(msn$X)[msn$fx_keep]
  out$random_n <- par_full[names(par_full) == "alpha_n"]
  names(out$random_n) <- colnames(msn$Z)
  out$sigma_n <- exp(opt$par[names(opt$par) == "log_sigma_n"])
  names(out$sigma_n) <- msn$re_labels
  out$gam_coefficients_n <- drop(as.matrix(msn$to_gam %*% c(out$coefficients_n, out$random_n)))
  names(out$gam_coefficients_n) <- colnames(msn$G)
  out$Gsetup_n <- msn$Gsetup

  out$df.null <- nrow(ms$G)
  out$df.residual <- nrow(ms$G) - out$rank
  out$fitted <- out$p
  out$residuals <- rep(0, nrow(data))
  out$null.deviance <- NA
  out$deviance <- NA
  out$family <- binomial()

  out$hessian <- NULL
  out$Vb <- out$Vb_n <- out$Vp <- out$Vp_n <- out$Vp_joint <- NULL
  if (hessian) {
    # covariance of the fixed and random effects of both models
    if (is.null(rand)) {
      V <- try(solve(obj$he(opt$par)), silent = TRUE)
      nms <- names(opt$par)
    } else {
      V <- try(as.matrix(solve(out$sdrep$jointPrecision)), silent = TRUE)
      nms <- rownames(out$sdrep$jointPrecision)
    }
    if (inherits(V, "try-error")) {
      warning("the covariance matrix could not be computed")
    } else {
      b <- which(nms == "beta")
      bn <- which(nms == "beta_n")
      out$Vb <- V[b, b, drop = FALSE]
      dimnames(out$Vb) <- list(names(out$coefficients), names(out$coefficients))
      out$Vb_n <- V[bn, bn, drop = FALSE]
      dimnames(out$Vb_n) <- list(names(out$coefficients_n), names(out$coefficients_n))

      # covariance of the coefficients of the columns of Gfit and G_n
      # for the negative binomial, log(theta) is the last row and column
      lt <- which(nms == "log_theta")
      keep <- c(b, which(nms == "alpha"), bn, which(nms == "alpha_n"), lt)
      to_gam <- Matrix::bdiag(c(
        list(
          ms$to_gam[ms$gam_cols, , drop = FALSE],
          msn$to_gam[msn$gam_cols, , drop = FALSE]
        ),
        if (length(lt)) list(Matrix::Diagonal(1))
      ))
      Vj <- as.matrix(to_gam %*% V[keep, keep, drop = FALSE] %*% Matrix::t(to_gam))
      np <- length(ms$gam_cols)
      nn <- length(msn$gam_cols)
      jn <- c(colnames(ms$Gfit), paste0("n:", colnames(msn$G)[msn$gam_cols]),
        if (length(lt)) "log_theta")
      dimnames(Vj) <- list(jn, jn)
      out$Vp_joint <- Vj
      out$Vp <- Vj[seq_len(np), seq_len(np), drop = FALSE]
      dimnames(out$Vp) <- list(colnames(ms$Gfit), colnames(ms$Gfit))
      out$Vp_n <- Vj[np + seq_len(nn), np + seq_len(nn), drop = FALSE]
      dimnames(out$Vp_n) <- list(colnames(msn$G)[msn$gam_cols], colnames(msn$G)[msn$gam_cols])
    }
  }

  out$convergence <- opt$convergence
  class(out) <- c("efn2", "efp2", "glm", "lm")
  out
}

# the expected number of fish not caught, given the total catch. For a
# Poisson number of fish this is lambda (1 - P), whatever the catch. For a
# negative binomial number of fish, with size theta, the number of fish is
# Poisson with mean lambda G, G ~ Gamma(theta, theta), and given the catch G
# is Gamma(theta + T, theta + lambda P), so the expected number not caught
# is lambda (1 - P) (theta + T) / (theta + lambda P). This tends to the
# Poisson value as theta tends to infinity, and to T (1 - P) / P, so that
# the abundance is T / P as estimated from the conditional likelihood, as
# theta tends to zero.
expected_missed <- function(catch, lambda, P, theta = NULL) {
  if (is.null(theta)) {
    lambda * (1 - P)
  } else {
    lambda * (1 - P) * (theta + catch) / (theta + lambda * P)
  }
}

# names without stats::setNames
setNames_ <- function(x, nm) {
  names(x) <- nm
  x
}
