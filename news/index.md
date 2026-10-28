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
