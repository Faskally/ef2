# Cumulative capture probability of each sample

The probability that a fish is caught on one of the passes of a sample,
\\P = 1 - \prod_i (1 - p_i)\\, where \\p_i\\ is the capture probability
on pass \\i\\.

## Usage

``` r
capture_prob(object, newdata = NULL, se.fit = FALSE, level = 0.95)
```

## Arguments

- object:

  an efp2 fitted model

- newdata:

  an optional data.frame with one row per pass, containing the
  covariates in the model and the variables used for `pass` and `id` in
  the call to
  [`efp2`](https://faskally.github.io/ef2/reference/efp2.md). If omitted
  the data used to fit the model are used. An offset is only included
  for the data used to fit the model.

- se.fit:

  if TRUE standard errors and confidence intervals are returned

- level:

  the confidence level

## Value

a data.frame with one row per sample, in the order the samples first
appear in the data, with columns `id`, `npass` the number of passes and
`P`. If `se.fit = TRUE` there are also columns `se`, the standard error
of `logit(P)`, and `lower` and `upper`, the confidence limits for `P`.

## Details

Standard errors are computed on the logit scale by the delta method,
using the covariance of the fixed and random effects. The derivative of
\\\textrm{logit}(P)\\ with respect to the linear predictor of pass \\i\\
is \\p_i / P\\. Confidence intervals are computed on the logit scale and
back transformed.

## References

Glover, R. S., Fryer, R. J., Soulsby, C., Bacon, P. J. and Malcolm, I.
A. (2018) Incorporating estimates of capture probability and river
network covariance in novel habitat - abundance models. Ecological
Indicators, 93, 302-315.
[doi:10.1016/j.ecolind.2018.05.013](https://doi.org/10.1016/j.ecolind.2018.05.013)

## See also

[`abundance`](https://faskally.github.io/ef2/reference/abundance.md)

## Examples

``` r
ef2_data$pass12 <- factor(pmin(ef2_data$pass, 2))
ef2_data$sampleID <- as.numeric(factor(paste(ef2_data$date,
  ef2_data$siteID, ef2_data$lifestage)))
ef2_data <- ef2_data[order(ef2_data$sampleID, ef2_data$pass), ]
m1 <- efp2(count ~ lifestage + pass12, data = ef2_data,
  pass = pass, id = sampleID)
head(capture_prob(m1, se.fit = TRUE))
#>   id npass         P        se     lower     upper
#> 1  1     3 0.9460482 0.2445284 0.9156713 0.9658903
#> 2  2     3 0.9858539 0.2609472 0.9766293 0.9914693
#> 3  3     3 0.9460482 0.2445284 0.9156713 0.9658903
#> 4  4     3 0.9858539 0.2609472 0.9766293 0.9914693
#> 5  5     3 0.9460482 0.2445284 0.9156713 0.9658903
#> 6  6     3 0.9858539 0.2609472 0.9766293 0.9914693
```
