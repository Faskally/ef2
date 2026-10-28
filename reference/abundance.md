# Estimate abundance and density from removal samples

Abundance is estimated as \\\hat N = T / P\\, where \\T\\ is the total
catch over the passes of a sample and \\P\\ is the cumulative capture
probability, see
[`capture_prob`](https://faskally.github.io/ef2/reference/capture_prob.md).

## Usage

``` r
abundance(
  object,
  newdata = NULL,
  area = NULL,
  by = NULL,
  target = c("realised", "expected"),
  level = 0.95
)
```

## Arguments

- object:

  an efp2 fitted model

- newdata:

  an optional data.frame as for
  [`capture_prob`](https://faskally.github.io/ef2/reference/capture_prob.md),
  which must also contain the counts used as the response in the model

- area:

  an optional name of a column in the data giving the area fished, or a
  vector of areas with one value per row of the data. If given, density
  (abundance per unit area) is also returned. With `by`, the area must
  be the same for all samples in a group.

- by:

  an optional character vector of the names of variables in the data,
  which must be the same for all passes of a sample. Abundance is summed
  over the samples with each combination of their values.

- target:

  the abundance to estimate, `"realised"` (the default) or `"expected"`,
  see Details

- level:

  the confidence level

## Value

a data.frame with one row per sample, in the order the samples first
appear in the data, with columns `id`, `catch` the total catch, `P`, `N`
the estimated abundance, its standard error `se`, and confidence limits
`lower` and `upper`. If `area` is given there are also columns `area`,
`density`, `density_lower` and `density_upper`. With `by` there is one
row per group, sorted by the `by` variables, with the `by` variables and
`nsample`, the number of samples, in place of `id` and `P`.

## Details

Two targets are available, which have the same estimate but different
standard errors:

- `"realised"`:

  the number of fish \\N\\ present when the sample was taken. Given
  \\N\\, the catch is binomial, \\T \sim \textrm{Bin}(N, P)\\, and
  (Huggins, 1989) \$\$\textrm{var}(\hat N) \approx \hat N (1 - P) / P +
  \hat N^2 (1 - P)^2 \textrm{var}(\textrm{logit}(P)).\$\$ The lower
  confidence limit is not less than the catch.

- `"expected"`:

  the expected number of fish \\\lambda\\, for example density times
  area, when the number present is \\N \sim \textrm{Poisson}(\lambda)\\.
  The catch is then \\T \sim \textrm{Poisson}(\lambda P)\\, and the
  first term of the variance is \\\hat\lambda / P\\, which also includes
  the natural variation of \\N\\ about \\\lambda\\. The standard error
  is larger by a factor of about \\1 / \sqrt{1 - P}\\ when \\P\\ is well
  estimated.

Use `"realised"` to say how many fish were present, and `"expected"` to
describe the underlying abundance, for example when comparing sites or
years. A Poisson model for the counts with offset \\\log(A P)\\, where
\\A\\ is the area, is a model for the expected abundance.

Confidence intervals are computed on the log scale. When the total catch
is zero the estimated abundance is zero and the confidence limits are
`NA`.

With `by`, abundance is summed over the samples in each group, for
example over the lifestages caught at a site in a year. The variation in
the catches is independent between samples, but the estimates of \\P\\
are correlated, as they depend on the same model coefficients, and this
is allowed for: \$\$\textrm{var}(\sum_s \hat N_s) \approx \sum_s v_s +
w^\top \Sigma w,\$\$ where \\v_s\\ is \\\hat N_s (1 - P_s) / P_s\\ for
the realised abundance and \\\hat N_s / P_s\\ for the expected
abundance, \\w_s = \hat N_s (1 - P_s)\\ and \\\Sigma\\ is the covariance
matrix of \\\textrm{logit}(P_s)\\.

## References

Huggins, R. M. (1989) On the statistical analysis of capture
experiments. Biometrika, 76, 133-140.
[doi:10.1093/biomet/76.1.133](https://doi.org/10.1093/biomet/76.1.133)

## See also

[`capture_prob`](https://faskally.github.io/ef2/reference/capture_prob.md)

## Examples

``` r
ef2_data$pass12 <- factor(pmin(ef2_data$pass, 2))
ef2_data$sampleID <- as.numeric(factor(paste(ef2_data$date,
  ef2_data$siteID, ef2_data$lifestage)))
ef2_data <- ef2_data[order(ef2_data$sampleID, ef2_data$pass), ]
m1 <- efp2(count ~ lifestage + pass12, data = ef2_data,
  pass = pass, id = sampleID)
head(abundance(m1, area = "area"))
#>   id catch         P          N        se lower      upper area    density
#> 1  1    49 0.9460482  51.794403 1.8495074    49  55.549234  122 0.42454429
#> 2  2    42 0.9858539  42.602661 0.7975209    42  44.194802  122 0.34920214
#> 3  3    10 0.9460482  10.570286 0.7888318    10  12.235159  227 0.04656514
#> 4  4   132 0.9858539 133.894076 1.4715777   132 136.809604  227 0.58984174
#> 5  5     5 0.9460482   5.285143 0.5534132     5   6.489138  228 0.02318045
#> 6  6    64 0.9858539  64.918340 0.9944564    64  66.896993  228 0.28472956
#>   density_lower density_upper
#> 1    0.40163934    0.45532159
#> 2    0.34426230    0.36225248
#> 3    0.04405286    0.05389938
#> 4    0.58149780    0.60268548
#> 5    0.02192982    0.02846113
#> 6    0.28070175    0.29340787

# summed over lifestages for each site and year
head(abundance(m1, area = "area", by = c("siteID", "year")))
#>   siteID year nsample catch         N       se    lower     upper area
#> 1  Site1 2011       2    91  94.39706 2.044479  91.0000  98.49043  122
#> 2  Site1 2012       2   181 189.57232 3.635209 182.5797 196.83278  122
#> 3  Site1 2013       2   136 140.08545 2.261193 136.0000 144.58816  122
#> 4  Site1 2014       2   119 123.86583 2.531473 119.0023 128.92813  122
#> 5  Site1 2015       2    98 101.75358 2.167267  98.0000 106.09126  124
#> 6  Site2 2011       2   157 162.79521 2.804362 157.3905 168.38550  264
#>     density density_lower density_upper
#> 1 0.7737464     0.7459016     0.8072986
#> 2 1.5538715     1.4965547     1.6133835
#> 3 1.1482414     1.1147541     1.1851488
#> 4 1.0152937     0.9754286     1.0567880
#> 5 0.8205934     0.7903226     0.8555747
#> 6 0.6166485     0.5961762     0.6378239

# expected abundance, for comparing years
abundance(m1, by = "year", target = "expected")
#>   year nsample catch        N       se    lower    upper
#> 1 2011       6   355 368.2884 19.79013 331.4732 409.1925
#> 2 2012       6   688 719.7241 28.45972 666.0508 777.7226
#> 3 2013       6   473 483.5429 22.36209 441.6417 529.4196
#> 4 2014       6   493 511.5122 23.43563 467.5812 559.5707
#> 5 2015       6   273 280.2463 17.04151 248.7592 315.7189
```
