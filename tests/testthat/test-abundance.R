prep_data <- function() {
  data(ef2_data)
  ef2_data$pass12 <- factor(pmin(ef2_data$pass, 2))
  ef2_data$sampleID <- as.numeric(factor(paste(
    ef2_data$date, ef2_data$siteID, ef2_data$lifestage
  )))
  ef2_data[order(ef2_data$sampleID, ef2_data$pass), ]
}

# simulate removal samples with known abundance
sim_removal <- function(ns = 300, N = 40, p = 0.35, seed = 1) {
  set.seed(seed)
  data <- expand.grid(pass = 1:3, sample = 1:ns)
  Ntrue <- rpois(ns, N)
  left <- Ntrue
  data$n <- NA
  for (i in 1:3) {
    r <- data$pass == i
    data$n[r] <- rbinom(ns, left, p)
    left <- left - data$n[r]
  }
  list(data = data, N = Ntrue)
}

test_that("cumulative capture probability matches the fitted pass probabilities", {
  data <- prep_data()
  m1 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)
  cp <- capture_prob(m1)

  expect_equal(cp$id, unique(data$sampleID))
  expect_equal(cp$npass, as.vector(table(data$sampleID)[as.character(cp$id)]))
  by_hand <- tapply(fitted(m1), data$sampleID, function(p) 1 - prod(1 - p))
  expect_equal(cp$P, as.vector(by_hand[as.character(cp$id)]))
  # and the value reported by the likelihood
  expect_equal(cp$P, unname(m1$P[as.character(cp$id)]))
})

test_that("delta method standard errors match a finite difference gradient", {
  data <- prep_data()
  m1 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)
  cp <- capture_prob(m1, se.fit = TRUE)

  logit_P <- function(object) qlogis(capture_prob(object)$P)
  h <- 1e-6
  J <- sapply(m1$gam_cols, function(j) {
    up <- down <- m1
    up$gam_coefficients[j] <- up$gam_coefficients[j] + h
    down$gam_coefficients[j] <- down$gam_coefficients[j] - h
    (logit_P(up) - logit_P(down)) / (2 * h)
  })
  se_fd <- sqrt(rowSums((J %*% m1$Vp) * J))
  expect_equal(cp$se, se_fd, tolerance = 1e-6)

  expect_true(all(cp$lower < cp$P & cp$P < cp$upper))
  cp99 <- capture_prob(m1, se.fit = TRUE, level = 0.99)
  expect_true(all(cp99$lower < cp$lower & cp$upper < cp99$upper))
})

test_that("capture probability for new data", {
  data <- prep_data()
  m1 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)

  # single pass samples
  first <- data[data$pass == 1, ]
  cp1 <- capture_prob(m1, newdata = first)
  expect_equal(cp1$npass, rep(1, nrow(first)))
  expect_equal(cp1$P, unname(predict(m1, newdata = first, type = "response")))

  # the order of the rows does not matter
  set.seed(1)
  shuffled <- data[sample(nrow(data)), ]
  cp <- capture_prob(m1, se.fit = TRUE)
  cps <- capture_prob(m1, newdata = shuffled, se.fit = TRUE)
  cps <- cps[match(cp$id, cps$id), ]
  rownames(cps) <- NULL
  expect_equal(cps, cp)
})

test_that("abundance and density", {
  data <- prep_data()
  m1 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)
  a <- abundance(m1, area = "area")
  cp <- capture_prob(m1, se.fit = TRUE)

  expect_equal(a$id, cp$id)
  expect_equal(a$catch, as.vector(tapply(data$count, data$sampleID, sum)[as.character(a$id)]))
  expect_equal(a$N, a$catch / cp$P)
  expect_equal(a$se^2, a$N * (1 - cp$P) / cp$P + a$N^2 * (1 - cp$P)^2 * cp$se^2)
  expect_true(all(a$lower >= a$catch & a$lower <= a$N & a$N < a$upper))
  expect_equal(a$density, a$N / a$area)
  expect_equal(a$density_upper, a$upper / a$area)

  # area as a vector gives the same result
  expect_equal(abundance(m1, area = data$area), a)
})

test_that("abundance with a zero catch", {
  data <- prep_data()
  m1 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)
  zero <- data[data$sampleID == data$sampleID[1], ]
  zero$count <- 0
  expect_no_warning(a <- abundance(m1, newdata = zero))
  expect_equal(a$N, 0)
  expect_true(is.na(a$lower) && is.na(a$upper))
})

test_that("abundance confidence intervals have the nominal coverage", {
  sim <- sim_removal()
  m <- efp2(n ~ 1, data = sim$data, pass = pass, id = sample)
  a <- abundance(m)
  coverage <- mean(a$lower <= sim$N & sim$N <= a$upper)
  expect_gt(coverage, 0.9)
  expect_lt(coverage, 0.99)
})

test_that("informative errors", {
  data <- prep_data()
  m1 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)

  expect_error(capture_prob(m1, newdata = data[names(data) != "sampleID"]), "sampleID")
  expect_error(abundance(m1, newdata = data[names(data) != "count"]), "counts")
  expect_error(abundance(m1, area = 1:2), "area")
  bad_area <- data$area
  bad_area[2] <- bad_area[2] + 1
  expect_error(abundance(m1, area = bad_area), "same for all passes")

  m0 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID,
    hessian = FALSE)
  expect_error(capture_prob(m0, se.fit = TRUE), "hessian = TRUE")
})

test_that("abundance summed by group", {
  data <- prep_data()
  m1 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)
  a <- abundance(m1)

  # each sample its own group gives the same estimates
  by_id <- abundance(m1, by = "sampleID")
  by_id <- by_id[match(a$id, by_id$sampleID), ]
  expect_equal(by_id$N, a$N)
  expect_equal(by_id$se, a$se)
  expect_equal(by_id$lower, a$lower)

  # sums over the lifestages at each site in each year
  g <- abundance(m1, by = c("siteID", "year"), area = "area")
  samples <- unique(data[c("sampleID", "siteID", "year")])
  samples <- samples[match(a$id, samples$sampleID), ]
  key <- paste(samples$siteID, samples$year)
  gkey <- paste(g$siteID, g$year)
  expect_equal(g$N, as.vector(tapply(a$N, key, sum)[gkey]))
  expect_equal(g$catch, as.vector(tapply(a$catch, key, sum)[gkey]))
  expect_equal(g$nsample, rep(2, nrow(g)))
  expect_false(is.unsorted(order(g$siteID, g$year)))
  expect_equal(g$density, g$N / g$area)

  # variance from the full covariance of logit(P), with a finite
  # difference gradient
  cp <- capture_prob(m1)
  h <- 1e-6
  J <- sapply(m1$gam_cols, function(j) {
    up <- down <- m1
    up$gam_coefficients[j] <- up$gam_coefficients[j] + h
    down$gam_coefficients[j] <- down$gam_coefficients[j] - h
    (qlogis(capture_prob(up)$P) - qlogis(capture_prob(down)$P)) / (2 * h)
  })
  Sigma <- J %*% m1$Vp %*% t(J)
  w <- a$N * (1 - cp$P)
  var_g <- sapply(gkey, function(k) {
    s <- key == k
    sum(a$N[s] * (1 - cp$P[s]) / cp$P[s]) + drop(t(w[s]) %*% Sigma[s, s] %*% w[s])
  })
  expect_equal(g$se^2, unname(var_g), tolerance = 1e-6)

  # the estimates are positively correlated, so the variance of a sum is
  # more than the sum of the variances
  y <- abundance(m1, by = "year")
  year_of <- as.character(samples$year)
  expect_true(all(y$se^2 > tapply(a$se^2, year_of, sum)[as.character(y$year)]))
})

test_that("summed abundance confidence intervals have the nominal coverage", {
  sim <- sim_removal(ns = 600, N = 30, p = 0.3)
  sim$data$group <- (sim$data$sample - 1) %/% 3
  m <- efp2(n ~ 1, data = sim$data, pass = pass, id = sample)
  a <- abundance(m, by = "group")
  truth <- tapply(sim$N, (seq_along(sim$N) - 1) %/% 3, sum)
  coverage <- mean(a$lower <= truth & truth <= a$upper)
  expect_gt(coverage, 0.9)
  expect_lt(coverage, 0.99)
})

test_that("informative errors for summed abundance", {
  data <- prep_data()
  m1 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)

  expect_error(abundance(m1, by = "nothing"), "nothing")
  # pass varies within a sample
  expect_error(abundance(m1, by = "pass"), "same for all passes")
  # sites have different areas in different years
  expect_error(abundance(m1, by = "siteID", area = "area"), "single area")
})

test_that("expected abundance", {
  data <- prep_data()
  m1 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)
  r <- abundance(m1)
  e <- abundance(m1, target = "expected")
  cp <- capture_prob(m1, se.fit = TRUE)

  # same estimate, Poisson rather than binomial variance for the counts
  expect_equal(e$N, r$N)
  expect_equal(e$se^2, e$N / cp$P + e$N^2 * (1 - cp$P)^2 * cp$se^2)
  expect_true(all(e$se > r$se))
  # the lower limit is not truncated at the catch
  expect_true(any(e$lower < e$catch))
  expect_equal(e$lower, e$N * exp(-qnorm(0.975) * e$se / e$N))

  # and summed over groups
  eg <- abundance(m1, by = "year", target = "expected")
  rg <- abundance(m1, by = "year")
  expect_equal(eg$N, rg$N)
  samples <- unique(data[c("sampleID", "year")])
  year <- samples$year[match(r$id, samples$sampleID)]
  binom <- tapply(r$N * (1 - cp$P) / cp$P, year, sum)
  pois <- tapply(e$N / cp$P, year, sum)
  yr <- as.character(eg$year)
  expect_equal(eg$se^2 - rg$se^2, as.vector(pois[yr] - binom[yr]))

  expect_error(abundance(m1, target = "other"), "should be one of")
})

test_that("expected abundance confidence intervals have the nominal coverage", {
  lambda <- 40
  sim <- sim_removal(ns = 600, N = lambda, p = 0.35)
  m <- efp2(n ~ 1, data = sim$data, pass = pass, id = sample)
  e <- abundance(m, target = "expected")
  coverage <- mean(e$lower <= lambda & lambda <= e$upper)
  expect_gt(coverage, 0.9)
  expect_lt(coverage, 0.99)

  # the realised intervals are too narrow for the expected abundance
  r <- abundance(m)
  expect_lt(mean(r$lower <= lambda & lambda <= r$upper), 0.85)
})
