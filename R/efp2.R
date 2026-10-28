#' Estimate capture probabilities from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters
#'
#'
#' @param formula a formula object. Random effects are specified using
#'   \code{s(f, bs = "re")} terms, where \code{f} is a factor, and smooths
#'   using other \code{s()} terms, see Details.
#' @param data a data.frame containing all relevant info
#' @param pass a vector of integers giving the pass number of the observation
#' @param id a vector of integers identifying an observation (a set of electrofishing passes)
#' @param offset an optional offset for the linear predictor of capture probability
#' @param verbose if TRUE optimiser messages are printed to the screen
#' @param hessian if TRUE the hessian is computed and the covariance matrix of the parameters is returned via Vb
#' @param fit if TRUE model is fitted if FALSE the data that would be passed to the optimiser is returned
#' @param control a list of control commands to be passed to \code{\link[stats]{nlminb}}
#' @param sparse if TRUE (the default) penalised smooths are set up so that
#'   sparse bases and penalties stay sparse, see Details. If FALSE the
#'   constrained basis and penalties from \code{mgcv} are used, which are
#'   dense and so slower to fit.
#' @return If \code{fit = TRUE}, an object of class \code{efp2}, a list which
#'   inherits from \code{glm}, containing among other things the fixed effect
#'   estimates (\code{coefficients}), the random effect estimates
#'   (\code{random}), the standard deviations of the random effects and
#'   smooths (\code{sigma}),
#'   the fitted capture probabilities (\code{fitted}), the log likelihood
#'   (\code{llik}) and the variance matrix of the fixed effects (\code{Vb}).
#'   See \code{\link{efp2-methods}} for the methods available. If
#'   \code{fit = FALSE}, a list containing the data that would be passed to
#'   the optimiser.
#'
#' @references
#' Millar, C. P., Fryer, R. J., Millidine, K. J. and Malcolm, I. A. (2016)
#' Modelling capture probability of Atlantic salmon (Salmo salar) from a
#' diverse national electrofishing dataset: Implications for the estimation
#' of abundance. Fisheries Research, 177, 1-12.
#' \doi{10.1016/j.fishres.2016.01.001}
#'
#' @details
#' Random effects on the logit capture probability are included by adding
#' \code{s(f, bs = "re")} terms to the formula, where \code{f} is a factor
#' in \code{data}. For example, \code{s(visitID, bs = "re")} adds between
#' visit variation and \code{s(sampleID, bs = "re")} adds between sample
#' variation. Each term has its own standard deviation, and the
#' random effects are integrated out of the likelihood using the Laplace
#' approximation, so that \code{logLik} and \code{AIC} are based on the marginal
#' likelihood. Fitted values are conditional on the estimated random effects.
#'
#' Smooth terms can also be included in the formula. These are penalised,
#' and the amount of smoothing is estimated: a smooth is fitted as a random
#' effect, in which the part of the smooth that is not penalised (for example
#' a linear trend) is a fixed effect, and the standard deviation of the
#' random part controls how smooth the term is. A standard deviation close to
#' zero means that only the unpenalised part of the smooth is needed. As for
#' other random effects, \code{logLik} and \code{AIC} are based on the
#' marginal likelihood.
#'
#' The recommended smooths are first and second order random walks, using
#' the Markov random field basis of \code{mgcv} with a penalty from
#' \code{\link{rw_penalty}}:
#' \code{s(f, bs = "mrf", xt = list(penalty = rw_penalty(levels(f))))},
#' where \code{f} is a factor whose levels are the steps of the walk. The
#' basis and the penalty are sparse, and are kept sparse when fitting, so
#' these models are quick to fit even when there are many levels. The
#' standard deviation is that of the steps of the random walk, and the
#' starting values of the walk (the first level for a first order walk, the
#' first two for a second order walk) are fixed effects: this means that the
#' intercept is the logit capture probability at the first level of the
#' walk, and that a second order walk has an additional fixed effect for its
#' linear trend, which is on an arbitrary scale.
#'
#' Any other penalised smooth from \code{mgcv} can be used, such as
#' \code{s(x)} or \code{te(x, z)}, but their penalties are dense so they
#' will take longer to fit, particularly if the basis dimension \code{k} is
#' large. To fit an unpenalised regression spline, whose flexibility is set
#' entirely by \code{k}, use \code{fx = TRUE}, for example
#' \code{s(x, k = 3, fx = TRUE)}.
#'
#' The fixed effects are estimated by maximum likelihood, so the likelihood
#' depends slightly on how the unpenalised part of a smooth is separated
#' from the penalised part. With \code{sparse = TRUE} one coefficient of the
#' smooth is fixed for each function that is not penalised, as described
#' above for random walks. With \code{sparse = FALSE}, and for smooths with
#' \code{by} variables or models without an intercept, the identifiability
#' constraints of \code{mgcv} are applied to the smooth first, which makes
#' the basis and the penalty dense. The two approaches give very similar,
#' but not identical, fits.
#'
#' @examples
#' # create two electrofishing site visits with 3 and 4 passes and 2 lifestages
#' ef2_data <- data.frame(n      = c(100, 53, 24, 50, 26, 12,
#'                                  100, 53, 24, 50, 26, 12),
#'                       pass   = c(  1,  2,  3,  1,  2,  3,
#'                                    1,  2,  3,  1,  2,  3),
#'                       stage  = c(  1,  1,  1,  2,  2,  2,
#'                                    1,  1,  1,  2,  2,  2),
#'                       sample = c(  1,  1,  1,  2,  2,  2,
#'                                    3,  3,  3,  4,  4,  4))
#'
#' ef_data2 <- data.frame(n      = c(100, 53, 24, 50, 26, 12,
#'                                   100, 53, 24, 12, 50, 26, 12, 6),
#'                        pass   = c(  1,  2,  3,  1,  2,  3,
#'                                     1,  2,  3,  4,  1,  2,  3, 4),
#'                        stage  = c(  1,  1,  1,  2,  2,  2,
#'                                     1,  1,  1,  1,  2,  2,  2, 2),
#'                        sample = c(  1,  1,  1,  2,  2,  2,
#'                                     3,  3,  3,  3,  4,  4,  4, 4))
#'
#' ef_data3 <- data.frame(n      = c(100, 53, 24, 50, 26, 12, 40,
#'                                   100, 53, 24, 12, 50, 26, 12, 6, 40),
#'                        pass   = c(  1,  2,  3,  1,  2,  3, 1,
#'                                     1,  2,  3,  4,  1,  2,  3, 4, 1),
#'                        stage  = c(  1,  1,  1,  2,  2,  2, 1,
#'                                     1,  1,  1,  1,  2,  2,  2, 2, 2),
#'                        sample = c(  1,  1,  1,  2,  2,  2, 5,
#'                                     3,  3,  3,  3,  4,  4,  4, 4, 6))
#'
#' # Fit a simple model
#' m2 <- efp2(n ~ 1 + factor(stage), data = ef2_data, pass = pass, id = sample)
#' cbind(ef2_data, fit = fitted(m2))
#' m3 <- efp2(n ~ 1 + factor(stage), data = ef_data2, pass = pass, id = sample)
#' cbind(ef_data2, fit = fitted(m3))
#' m4 <- efp2(n ~ 1 + factor(stage), data = ef_data3, pass = pass, id = sample)
#' cbind(ef_data3, fit = fitted(m4))
#'
#' # create two electrofishing site visits with 3 and 4 passes and 2 lifestages
#' ef2_data <- data.frame(n      = c(200, 53, 24, 100, 26, 12,
#'                                  200, 53, 24, 100, 26, 12),
#'                       pass   = c(  1,  2,  3,  1,  2,  3,
#'                                    1,  2,  3,  1,  2,  3),
#'                       stage  = c(  1,  1,  1,  2,  2,  2,
#'                                    1,  1,  1,  2,  2,  2),
#'                       sample = c(  1,  1,  1,  2,  2,  2,
#'                                    3,  3,  3,  4,  4,  4))
#' # Fit a simple model
#' m2 <- efp2(n ~ 1 + factor(stage) + factor(replace(pass, pass> 2, 2)),
#'           data = ef2_data, pass = pass, id = sample)
#' out <- cbind(ef2_data, p = fitted(m2, type = "p"))
#' out
#'
#' # a smooth effect of a continuous covariate, day of the year, using a
#' # second order random walk
#' data(ef2_data)
#' ef2_data$sampleID <- factor(paste(ef2_data$date, ef2_data$siteID,
#'                                  ef2_data$lifestage, ef2_data$species))
#' ef2_data$pass12 <- factor(replace(ef2_data$pass, ef2_data$pass > 2, 2))
#' days <- min(ef2_data$doy):max(ef2_data$doy)
#' ef2_data$fdoy <- factor(ef2_data$doy, levels = days)
#' rw2 <- rw_penalty(days, order = 2)
#' m5 <- efp2(count ~ siteID + pass12 + lifestage +
#'             s(fdoy, bs = "mrf", xt = list(penalty = rw2)),
#'           data = ef2_data, pass = pass, id = sampleID)
#' m5
#'
#' # predicted capture probability over the smooth, with standard errors
#' newdata <- data.frame(siteID = "Site1", pass12 = "1", lifestage = "Fry",
#'                       fdoy = factor(days, levels = days))
#' pred <- predict(m5, newdata = newdata, type = "response", se.fit = TRUE)
#' cbind(newdata, fit = pred$fit, se = pred$se.fit)
#'
#' # an unpenalised regression spline
#' m6 <- efp2(count ~ siteID + pass12 + lifestage + s(doy, k = 3, fx = TRUE),
#'           data = ef2_data, pass = pass, id = sampleID)
#' m6
#'
#' @importFrom RTMB MakeADFun sdreport ADREPORT REPORT dnorm dgmrf
#' @importFrom RTMB logspace_add logspace_sub plogis
#' @importFrom mgcv gam
#' @importFrom methods as
#' @import Matrix
#' @importFrom stats binomial
#' @importFrom stats nlminb
#' @importFrom stats ave
#'
#' @export
efp2 <- function(formula, data = NULL, pass = pass, id = id,
  offset = NULL, verbose = FALSE, hessian = TRUE,
  fit = TRUE, control = list(), sparse = TRUE) {

  # some checks
  if (is.null(data)) stop("must supply data")
  
  # get pass information
  pass <- substitute(pass)
  pass <- eval(pass, data, environment(formula))
  # get id
  id <- substitute(id)
  id <- eval(id, data, environment(formula))
  
  # this is necessary to ensure that passes are ordered within samples so that
  # the order of input data and model outputs match. Likelihood function 
  # requires data ordered by pass.
  if (any(tapply(pass, id, is.unsorted))) {
    stop("data must be sorted by pass")
  }

  # set up offset
  if (is.null(offset)) {
    offset <- rep(0, nrow(data))
  }

  # set up model
  ms <- model_setup(formula, data, sparse = sparse)
  Gsetup <- ms$Gsetup
  G <- ms$G
  X <- ms$X
  Z <- ms$Z
  blocks <- ms$blocks
  re_labels <- ms$re_labels
  re_cols <- ms$re_cols
  fx_cols <- ms$fx_cols
  fx_keep <- ms$fx_keep
  to_gam <- ms$to_gam
  gam_cols <- ms$gam_cols
  Gfit <- ms$Gfit

  # get data for likelihood
  ord_output <- order(id, pass)

  efdat <- c(
    list(y = Gsetup$y, sample_id = factor(id)),
    model_data(ms, offset)
  )
  # the standard deviation of each random effect
  efdat$re_id <- as.integer(unlist(lapply(blocks, function(x) rep(x$sig[1], length(x$idx)))))

  if (!fit) {
    return(efdat)
  }

  params <-
    list(
      alpha = rep(0, ncol(efdat$Z)),
      beta = rep(0, ncol(efdat$X)),
      log_sigma = rep(0, length(re_labels))
    )

  if (length(blocks)) {
    rand <- "alpha"
  } else {
    rand <- NULL
  }

  obj <-
    MakeADFun(
      function(params) ef_nll(params, efdat),
      params,
      random = rand,
      silent = !verbose
    )

  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
  
  # check convergence
  
  if (opt$convergence != 0) {
    warning("model has not converged: ", opt$message)
  }
  
  # set up output object
  out <- list()

  # evaluate at the optimum, this also sets the random effects to their
  # estimates at the optimum
  out$llik <- -1 * as.numeric(obj$fn(opt$par))
  par_full <- obj$env$last.par
  out$rep <- obj$report(par_full)

  # extract transformed parameters
  out$p <- unname(unlist(out$rep$ps))
  # Re-order to match input data
  out$p <- out$p[order(ord_output)]
  # probability of capture in each pass, allowing for the fish removed in
  # previous passes
  out$pi <- as.vector(
    ave(out$p, factor(id), FUN = function(x) x * c(1, cumprod(1 - x)[-length(x)]))
  )
  # cumulative capture probability of each sample, named by id
  out$P <- out$rep$P
  names(out$P) <- levels(efdat$sample_id)
  
  fx_pars <- grep("beta", names(opt$par))
  out$sdrep <- muffle_nan(
    sdreport(obj, getJointPrecision = hessian && length(blocks) > 0)
  )
  check_sdreport(out$sdrep, c(colnames(X)[fx_keep], re_labels))
  # keep data for re-fitting
  # efdat needs to be reordered to be compatible with other structures 
  out$data <- data
  out$offset <- offset
  out$formula <- formula # for printing and summary
  out$terms <- Gsetup$terms
  out$call <- match.call()
  out$rank <- length(opt$par)
  out$aic <- -2 * out$llik + 2 * out$rank
  out$G <- G
  out$Gfit <- Gfit
  out$fx_cols <- fx_cols
  out$re_cols <- re_cols
  out$gam_cols <- gam_cols
  out$re_labels <- re_labels
  out$coefficients <- opt$par[fx_pars]
  names(out$coefficients) <- colnames(X)[fx_keep]
  out$random <- par_full[names(par_full) == "alpha"]
  names(out$random) <- colnames(Z)
  out$sigma <- exp(opt$par[names(opt$par) == "log_sigma"])
  names(out$sigma) <- re_labels
  # coefficients for the columns of G, used for prediction
  out$gam_coefficients <- drop(as.matrix(to_gam %*% c(out$coefficients, out$random)))
  names(out$gam_coefficients) <- colnames(G)
  out$df.null <- nrow(G)
  out$df.residual <- nrow(G) - out$rank
  out$fitted <- out$p

  out$residuals <- rep(0, nrow(data))
  out$null.deviance <- NA
  out$deviance <- NA
  out$family <- binomial()
  if (!hessian) {
    out$hessian <- NULL
    out$Vb <- NULL
    out$Vp <- NULL
  } else if (length(blocks)) {
    # covariance of fixed and random effects
    out$hessian <- NULL
    V <- try(as.matrix(solve(out$sdrep$jointPrecision)))
    if (inherits(V, "try-error")) {
      out$Vb <- out$Vp <- V
    } else {
      nms <- rownames(out$sdrep$jointPrecision)
      which_b <- which(nms == "beta")
      which_a <- which(nms == "alpha")
      out$Vb <- V[which_b, which_b, drop = FALSE]
      rownames(out$Vb) <- colnames(out$Vb) <- names(out$coefficients)
      # covariance of the coefficients of the columns of Gfit
      keep <- c(which_b, which_a)
      to_gam_fit <- to_gam[gam_cols, , drop = FALSE]
      out$Vp <- as.matrix(to_gam_fit %*% V[keep, keep, drop = FALSE] %*% t(to_gam_fit))
      rownames(out$Vp) <- colnames(out$Vp) <- colnames(Gfit)
    }
  } else {
    out$hessian <- obj$he(opt$par)[fx_pars, fx_pars, drop = FALSE]
    rownames(out$hessian) <- colnames(out$hessian) <- colnames(Gfit)
    out$Vb <- try(solve(out$hessian))
    out$Vp <- out$Vb
  }

  out$Gsetup <- Gsetup
  out$convergence <- opt$convergence

  class(out) <- c("efp2", "glm", "lm")
  out
}
