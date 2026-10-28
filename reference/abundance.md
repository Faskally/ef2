# Estimate abundance and density from removal samples

Abundance is estimated as \\\hat N = T / P\\, where \\T\\ is the total
catch over the passes of a sample and \\P\\ is the cumulative capture
probability, see
[`capture_prob`](https://faskally.github.io/ef2/reference/capture_prob.md).

## Usage

``` r
abundance(object, newdata = NULL, area = NULL, level = 0.95)
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
  (abundance per unit area) is also returned.

- level:

  the confidence level

## Value

a data.frame with one row per sample, in the order the samples first
appear in the data, with columns `id`, `catch` the total catch, `P`, `N`
the estimated abundance, its standard error `se`, and confidence limits
`lower` and `upper`. If `area` is given there are also columns `area`,
`density`, `density_lower` and `density_upper`.

## Details

The variance of \\\hat N\\ combines the binomial sampling variation in
the catch, given \\P\\, and the uncertainty in the estimate of \\P\\
(Huggins, 1989): \$\$\textrm{var}(\hat N) \approx \hat N (1 - P) / P +
\hat N^2 (1 - P)^2 \textrm{var}(\textrm{logit}(P)).\$\$ Confidence
intervals are computed on the log scale, and the lower limit is not less
than the total catch. When the total catch is zero the estimated
abundance is zero and the confidence limits are `NA`.

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
#>   id catch         P        se          N lower      upper area    density
#> 1  1    49 0.9460482 1.8495074  51.794403    49  55.549234  122 0.42454429
#> 2  2    42 0.9858539 0.7975209  42.602661    42  44.194802  122 0.34920214
#> 3  3    10 0.9460482 0.7888318  10.570286    10  12.235159  227 0.04656514
#> 4  4   132 0.9858539 1.4715777 133.894076   132 136.809604  227 0.58984174
#> 5  5     5 0.9460482 0.5534132   5.285143     5   6.489138  228 0.02318045
#> 6  6    64 0.9858539 0.9944564  64.918340    64  66.896993  228 0.28472956
#>   density_lower density_upper
#> 1    0.40163934    0.45532159
#> 2    0.34426230    0.36225248
#> 3    0.04405286    0.05389938
#> 4    0.58149780    0.60268548
#> 5    0.02192982    0.02846113
#> 6    0.28070175    0.29340787
```
