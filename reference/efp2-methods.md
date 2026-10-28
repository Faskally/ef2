# Methods for efp2 fits

Methods to summarise and extract information from a model fitted using
[`efp2`](https://faskally.github.io/ef2/reference/efp2.md).

## Usage

``` r
# S3 method for class 'efp2'
fitted(object, type = "p", ...)

# S3 method for class 'efp2'
vcov(object, ...)

# S3 method for class 'efp2'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'efp2'
logLik(object, ...)

# S3 method for class 'efp2'
formula(x, ...)

# S3 method for class 'efp2'
summary(object, ...)

# S3 method for class 'summary.efp2'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'efp2'
residuals(object, ...)

# S3 method for class 'efp2'
anova(object, ...)

# S3 method for class 'efp2'
plot(x, ...)

# S3 method for class 'efp2'
model.matrix(object, ...)

# S3 method for class 'efp2'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object, x:

  an efp2 fitted model

- type:

  the type of fitted value: `"p"` for the capture probability in each
  pass, `"pi"` for the probability that a fish is first caught in each
  pass, `"both"` for a list containing both and `"lpmatrix"` for the
  design matrix.

- ...:

  further arguments, ignored unless stated

- parm:

  the parameters to compute confidence intervals for, either a vector of
  numbers or names. The default is all fixed effects.

- level:

  the confidence level required

- digits:

  the number of significant digits to print

## Value

`fitted` returns a vector, list or matrix depending on `type`;
`model.matrix` a matrix; `vcov` a matrix; `confint` a matrix with
columns giving the lower and upper confidence limits; `logLik` an object
of class `logLik`; `formula` the model formula; `summary` an object of
class `summary.efp2`, a list containing the coefficient table
(`coefficients`) and the table of random effect standard deviations
(`sigma`). The `print` methods return their argument invisibly.

## Details

The variance matrix returned by `vcov` is for the fixed effects only,
the joint variance matrix of the fixed and random effects is available
as `object$Vp`. Confidence intervals are Wald intervals based on this
variance matrix.

The log likelihood is the marginal log likelihood when the model
contains random effects or penalised smooths, and the degrees of freedom
are the number of fixed effects plus the number of standard deviations.

Residuals, analysis of deviance tables and diagnostic plots are not
available for efp2 fits, and the corresponding methods return an error.

## Examples

``` r
data(ef2_data)
ef2_data$sampleID <- factor(paste(ef2_data$date, ef2_data$siteID,
                                 ef2_data$lifestage, ef2_data$species))
m1 <- efp2(count ~ lifestage, data = ef2_data, pass = pass, id = sampleID)

summary(m1)
#> 
#> Call:
#> efp2(formula = count ~ lifestage, data = ef2_data, pass = pass, 
#>     id = sampleID)
#> 
#> Formula:
#> count ~ lifestage
#> <environment: 0x55cc212670b0>
#> 
#> No random effects
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z^2|)
#> (Intercept)     0.6820    0.07119   9.580  9.704e-22
#> lifestageParr   0.6212    0.10291   6.036  1.579e-09
#> 
#> Degrees of Freedom: 90 Total (i.e. Null);  88 Residual
#> AIC: 3192
vcov(m1)
#>                (Intercept) lifestageParr
#> (Intercept)    0.005067977  -0.005067977
#> lifestageParr -0.005067977   0.010590042
confint(m1)
#>                   2.5 %    97.5 %
#> (Intercept)   0.5424690 0.8215276
#> lifestageParr 0.4194697 0.8228614
logLik(m1)
#> 'log Lik.' -1594.226 (df=2)
head(fitted(m1))
#> [1] 0.7863670 0.6641846 0.7863670 0.6641846 0.7863670 0.6641846
```
