# log likelihood of a single multipass electrofishing sample, and the log
# of the cumulative capture probability over all passes, computed on
# the log scale from the linear predictor so that it stays finite when
# capture probabilities are close to 0 or 1
ll_sample <- function(eta, y) {
  npasses <- length(y)
  log_p <- -logspace_add(0 * eta, -eta)
  log_q <- -logspace_add(0 * eta, eta) # log(1 - p)
  # log probability of first capture on each pass
  log_pi <- log_p
  if (npasses > 1) {
    log_pi[-1] <- log_p[-1] + cumsum(log_q)[-npasses]
  }
  # log probability of capture over all passes, log(1 - prod(1 - p))
  log_P <- logspace_sub(0, sum(log_q))
  list(ll = sum(y * log_pi) - sum(y) * log_P, log_P = log_P)
}

# linear predictor of a model set up by model_setup
linear_predictor <- function(dat, beta, alpha) {
  eta <- drop(dat$X %*% beta) + dat$offset
  if (length(alpha)) {
    # add random effects
    eta <- eta + drop(dat$Z %*% alpha)
  }
  eta
}

# negative log density of the random effects and penalised smooths
re_nll <- function(alpha, sigma, blocks) {
  nll <- 0
  for (block in blocks) {
    u <- alpha[block$idx]
    if (is.null(block$S)) {
      # independent random effects
      nll <- nll - sum(dnorm(u, 0, sigma[block$sig], log = TRUE))
    } else {
      # gaussian with sparse precision matrix sum_j S_j / sigma_j^2, with
      # one term for each penalty
      Q <- block$S[[1]] / sigma[block$sig[1]]^2
      for (j in seq_along(block$S)[-1]) {
        Q <- Q + block$S[[j]] / sigma[block$sig[j]]^2
      }
      nll <- nll - dgmrf(u, 0, Q, log = TRUE)
    }
  }
  nll
}

# negative log likelihood of the catches conditional on the total catch of
# each sample, including the random effects for capture probability
#   returns the negative log likelihood, the log cumulative capture
#   probability and the capture probabilities of each sample, and the
#   standard deviations of the random effects
capture_nll <- function(alpha, beta, log_sigma, efdat) {
  eta <- linear_predictor(efdat, beta, alpha)

  p <- plogis(eta)
  idx <- split(seq_along(efdat$y), efdat$sample_id)
  ps <- unname(lapply(idx, function(i) p[i]))
  etas <- unname(lapply(idx, function(i) eta[i]))
  ys <- lapply(idx, function(i) efdat$y[i])

  sigma <- exp(log_sigma)
  nll <- re_nll(alpha, sigma, efdat$blocks)

  log_P <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    ll <- ll_sample(etas[[i]], ys[[i]])
    nll <- nll - ll$ll
    log_P[[i]] <- ll$log_P
  }

  list(nll = nll, log_P = do.call(c, log_P), ps = ps, sigma = sigma)
}

# negative log likelihood, the RTMB objective function
ef_nll <- function(params, efdat) {
  cap <- capture_nll(params$alpha, params$beta, params$log_sigma, efdat)

  # cumulative capture probability of each sample
  P <- exp(cap$log_P)
  REPORT(P)

  # delta method for sigma
  sigma <- cap$sigma
  if (length(sigma)) {
    ADREPORT(sigma)
  }
  REPORT(sigma)

  ps <- cap$ps
  REPORT(ps)

  cap$nll
}

# joint negative log likelihood of the catches and the number of fish, the
# RTMB objective function for efn2. The number of fish in sample s is
# Poisson (or negative binomial) with mean lambda_s, so the total catch is
# Poisson (or negative binomial with the same size) with mean lambda_s P_s,
# and the catches given the total catch are multinomial as in ef_nll. The
# likelihood is the product of the two.
#   ndat: the model for log(lambda), with one row per sample, in the order
#     of the levels of efdat$sample_id, y the total catch of each sample,
#     and distribution "poisson" or "negbin"
ef_nll_joint <- function(params, efdat, ndat) {
  cap <- capture_nll(params$alpha, params$beta, params$log_sigma, efdat)

  log_lambda <- linear_predictor(ndat, params$beta_n, params$alpha_n)
  sigma_n <- exp(params$log_sigma_n)
  nll <- cap$nll + re_nll(params$alpha_n, sigma_n, ndat$blocks)

  # the total catch has mean lambda P, and is Poisson if the number of fish
  # is Poisson, or negative binomial with size theta if the number of fish
  # is negative binomial
  y <- ndat$y
  log_mu <- log_lambda + cap$log_P
  if (identical(ndat$distribution, "negbin")) {
    theta <- exp(params$log_theta)
    log_theta_mu <- logspace_add(params$log_theta + 0 * log_mu, log_mu)
    nll <- nll - sum(
      lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
        theta * (params$log_theta - log_theta_mu) + y * (log_mu - log_theta_mu)
    )
    REPORT(theta)
    ADREPORT(theta)
  } else {
    nll <- nll - sum(y * log_mu - exp(log_mu) - lgamma(y + 1))
  }
  # the multinomial coefficients, which are left out of the conditional
  # likelihood, so that this is the full log likelihood of the catches
  nll <- nll - (sum(lgamma(y + 1)) - sum(lgamma(efdat$y + 1)))

  P <- exp(cap$log_P)
  lambda <- exp(log_lambda)
  REPORT(P)
  REPORT(lambda)

  sigma <- cap$sigma
  if (length(sigma)) {
    ADREPORT(sigma)
  }
  if (length(sigma_n)) {
    ADREPORT(sigma_n)
  }
  REPORT(sigma)
  REPORT(sigma_n)

  ps <- cap$ps
  REPORT(ps)

  nll
}
