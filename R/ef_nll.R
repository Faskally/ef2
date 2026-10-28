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

# negative log likelihood, the RTMB objective function
ef_nll <- function(params, efdat) {
  alpha <- params$alpha
  beta <- params$beta
  log_sigma <- params$log_sigma

  eta <- drop(efdat$X %*% beta) + efdat$offset
  if (length(alpha)) {
    # add random effects
    eta <- eta + drop(efdat$Z %*% alpha)
  }

  p <- plogis(eta)
  idx <- split(seq_along(efdat$y), efdat$sample_id)
  ps <- unname(lapply(idx, function(i) p[i]))
  etas <- unname(lapply(idx, function(i) eta[i]))
  ys <- lapply(idx, function(i) efdat$y[i])

  nll <- 0 # negative log likelihood

  sigma <- exp(log_sigma)

  # add random effects and penalised smooths
  for (block in efdat$blocks) {
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

  log_P <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    ll <- ll_sample(etas[[i]], ys[[i]])
    nll <- nll - ll$ll
    log_P[[i]] <- ll$log_P
  }

  # cumulative capture probability of each sample
  P <- exp(do.call(c, log_P))
  REPORT(P)

  # delta method for sigma
  if (length(sigma)) {
    ADREPORT(sigma)
  }
  REPORT(sigma)

  REPORT(ps)

  nll
}
