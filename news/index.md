# Changelog

## ef2 2.0.0

- package renamed from ef to ef2 so it can be installed alongside ef
  1.x; fitted models now have class `efp2`
- models are fitted with
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) rather than
  [`optim()`](https://rdrr.io/r/stats/optim.html) (BFGS), which falsely
  converged on some platforms; `control` is passed to
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html)
- the likelihood is computed on the log scale, so it stays finite when
  capture probabilities are close to 0 or 1
- `efp()` is renamed
  [`efp2()`](https://faskally.github.io/ef2/reference/efp2.md), so ef
  and ef2 can be attached together
- remove `logit()`, `invlogit()` and `transpar()`, use
  [`qlogis()`](https://rdrr.io/r/stats/Logistic.html) and
  [`plogis()`](https://rdrr.io/r/stats/Logistic.html) from stats instead
- `as.gam()` is renamed
  [`as.gam2()`](https://faskally.github.io/ef2/reference/as.gam2.md) and
  `ef_data` is renamed `ef2_data`, so no names clash with ef
- fix the choice of coefficients removed for the null space of a smooth,
  which could give a singular penalty and NaN likelihood on macOS
- when standard errors cannot be computed,
  [`efp2()`](https://faskally.github.io/ef2/reference/efp2.md) warns
  which parameters are affected, instead of “NaNs produced”
- new
  [`capture_prob()`](https://faskally.github.io/ef2/reference/capture_prob.md)
  and
  [`abundance()`](https://faskally.github.io/ef2/reference/abundance.md)
  give the cumulative capture probability, abundance and density of each
  sample with confidence intervals; the cumulative capture probability
  is also reported in `$P`
