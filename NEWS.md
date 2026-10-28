# ef2 2.0.0

* package renamed from ef to ef2 so it can be installed alongside ef 1.x;
  fitted models now have class `efp2`
* models are fitted with `nlminb()` rather than `optim()` (BFGS), which
  falsely converged on some platforms; `control` is passed to `nlminb()`
* the likelihood is computed on the log scale, so it stays finite when
  capture probabilities are close to 0 or 1
* `efp()` is renamed `efp2()`, so ef and ef2 can be attached together
* remove `logit()`, `invlogit()` and `transpar()`, use `qlogis()` and
  `plogis()` from stats instead
* `as.gam()` is renamed `as.gam2()` and `ef_data` is renamed `ef2_data`,
  so no names clash with ef
* fix the choice of coefficients removed for the null space of a smooth,
  which could give a singular penalty and NaN likelihood on macOS
* when standard errors cannot be computed, `efp2()` warns which parameters
  are affected, instead of "NaNs produced"
* new `capture_prob()` and `abundance()` give the cumulative capture
  probability, abundance and density of each sample with confidence
  intervals, optionally summed over groups of samples with `by`; the
  cumulative capture probability is also reported in `$P`
* `abundance(target = "expected")` gives confidence intervals for the
  expected abundance, when the number of fish present is Poisson

# ef 2.0.0 (2026-10-09)

* move likelihood computation to RTMB, the package no longer has compiled code
* random effects can be included using `s(f, bs = "re")` terms in the model
  formula, for example for between visit and between sample variation
* remove `overdispersion()` and `BICadj()`, these are replaced by visit and
  sample random effects
* remove the `sample_re` argument of `efp()`
* smooth terms are now penalised and the amount of smoothing is estimated,
  by fitting smooths as random effects. Previously smooths were unpenalised
  regression splines, which are still available using `fx = TRUE`, for
  example `s(x, k = 3, fx = TRUE)`
* the recommended smooths are first and second order random walks, using
  `s(f, bs = "mrf", xt = list(penalty = rw_penalty(levels(f))))`. The new
  function `rw_penalty()` creates the penalty, and these are fitted using
  sparse matrices so they are quick to fit even with many levels. Other
  smooths from `mgcv` can be used but are slower to fit
* add a `predict()` method for `efp` fits, with predictions for new data,
  standard errors and the option to exclude random effects
* `fitted(type = "pi")` now returns the probability of capture in each pass
* add day of the year (`doy`) to the `ef_data` example dataset, a continuous
  covariate which can be used to fit smooth terms
* add examples of smooth terms to the help and the README (#1)
* add `summary()`, `vcov()`, `confint()`, `logLik()` and `formula()` methods
  for `efp` fits. `residuals()`, `anova()` and `plot()` now give an
  informative error, rather than failing in the methods inherited from `glm`
* remove the unused `init` argument of `efp()`

# ef 1.2.0 (2020-10-18)

* move likelihood computation to TMB
* introduce random effects at the sample level

# ef 1.1.0 (2016-08-01)

* initial package developed using STAN for likelihood computations
