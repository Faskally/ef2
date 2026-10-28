# Bias in abundance estimates from removal models with random effects
#
# Replicates the simulation in vignettes/articles/simulate.R with different
# seeds, fits the model with visit and sample random effects, and compares
# estimators of abundance with the truth.
#
# Run from the package root: Rscript bias/bias_study.R [nrep] [ncores]

library(ef2)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
nrep <- if (length(args) > 0) as.integer(args[1]) else 100
ncores <- if (length(args) > 1) as.integer(args[2]) else detectCores()

sim_code <- readLines("vignettes/articles/simulate.R")
seed_line <- grep("^set.seed\\(", sim_code)
stopifnot(length(seed_line) == 1)

simulate <- function(seed) {
  code <- sim_code
  code[seed_line] <- paste0("set.seed(", seed, ")")
  env <- new.env()
  eval(parse(text = code), envir = env)
  env
}

one_rep <- function(seed) {
  sim <- simulate(seed)
  ef_sim <- sim$ef_sim
  ef_sim$fsampleID <- factor(ef_sim$sampleID)

  m <- efp2(
    count ~ lifestage + pass12 + s(year, k = 6) + s(doy, k = 6) +
      s(visitID, bs = "re") + s(fsampleID, bs = "re"),
    data = ef_sim, pass = pass, id = sampleID
  )

  a <- abundance(m)
  cp <- capture_prob(m, se.fit = TRUE)
  N <- sim$samples$N[match(a$id, sim$samples$sampleID)]
  true_P <- tapply(ef_sim$p, ef_sim$sampleID, function(p) 1 - prod(1 - p))
  true_P <- true_P[as.character(a$id)]

  # estimators of N
  l <- qlogis(cp$P)
  est <- list(
    plugin = a$N,
    # posterior mean of 1 / P, with logit(P) normal
    lognormal = a$catch * (1 + exp(-l + cp$se^2 / 2)),
    # with the true capture probabilities, for reference
    true_P = a$catch / true_P
  )

  q <- cut(true_P, quantile(true_P, 0:4 / 4), include.lowest = TRUE,
    labels = paste0("Q", 1:4))

  ratio <- sapply(est, function(x) c(
    total = sum(x) / sum(N),
    tapply(x, q, sum) / tapply(N, q, sum)
  ))

  # coverage of the realised intervals, centred on each estimator with the
  # plugin standard error
  z <- qnorm(0.975)
  coverage <- sapply(est[c("plugin", "lognormal")], function(x) {
    cv <- a$se / x
    lower <- pmax(a$catch, x * exp(-z * cv))
    upper <- x * exp(z * cv)
    mean(lower <= N & N <= upper, na.rm = TRUE)
  })

  list(
    seed = seed,
    ratio = ratio,
    coverage = coverage,
    sigma = m$sigma[c("s(visitID)", "s(fsampleID)")],
    slope = unname(coef(lm(qlogis(cp$P) ~ qlogis(true_P)))[2]),
    converged = m$convergence == 0
  )
}

res <- mclapply(seq_len(nrep), function(i) {
  tryCatch(one_rep(1000 + i), error = function(e) NULL)
}, mc.cores = ncores)
res <- Filter(Negate(is.null), res)

dir.create("bias/results", showWarnings = FALSE)
saveRDS(res, "bias/results/bias_study.rds")

cat("replicates:", length(res), " converged:", sum(sapply(res, `[[`, "converged")), "\n\n")

ratios <- simplify2array(lapply(res, `[[`, "ratio"))
cat("estimated / true abundance, mean (sd) over replicates\n")
m <- apply(ratios, 1:2, mean)
s <- apply(ratios, 1:2, sd)
print(matrix(sprintf("%.4f (%.4f)", m, s), nrow(m), dimnames = dimnames(m)), quote = FALSE)

cat("\ncoverage of 95% realised intervals\n")
print(round(rowMeans(sapply(res, `[[`, "coverage")), 3))

cat("\nrandom effect standard deviations (true 0.4 visit, 0.3 sample)\n")
print(round(rowMeans(sapply(res, `[[`, "sigma")), 3))

cat("\nslope of estimated on true logit(P):", round(mean(sapply(res, `[[`, "slope")), 3), "\n")
