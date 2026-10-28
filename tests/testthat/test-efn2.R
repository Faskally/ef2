prep_efn2_data <- function() {
  data(ef2_data)
  ef2_data$pass12 <- factor(pmin(ef2_data$pass, 2))
  ef2_data$sampleID <- as.numeric(factor(paste(
    ef2_data$date, ef2_data$siteID, ef2_data$lifestage
  )))
  ef2_data$visitID <- factor(paste(ef2_data$date, ef2_data$siteID))
  ef2_data[order(ef2_data$sampleID, ef2_data$pass), ]
}

test_that("the joint likelihood is the product of Poisson catches", {
  data <- prep_efn2_data()
  dat <- efn2(count ~ lifestage + pass12,
    abundance = ~ lifestage + offset(log(area)),
    data = data, pass = pass, id = sampleID, fit = FALSE)

  set.seed(1)
  params <- list(
    alpha = numeric(0), beta = c(0.5, 0.3, -0.2), log_sigma = numeric(0),
    alpha_n = numeric(0), beta_n = c(-1, 0.2), log_sigma_n = numeric(0)
  )
  nll <- ef2:::ef_nll_joint(params, dat$efdat, dat$ndat)

  # the same likelihood computed directly
  X <- model.matrix(~ lifestage + pass12, data)
  p <- plogis(drop(X %*% params$beta))
  pi <- ave(p, data$sampleID, FUN = function(x) x * c(1, cumprod(1 - x)[-length(x)]))
  samples <- data[!duplicated(data$sampleID), ]
  lambda <- exp(params$beta_n[1] + params$beta_n[2] * (samples$lifestage == "Parr")) *
    samples$area
  names(lambda) <- samples$sampleID
  direct <- -sum(dpois(data$count, lambda[as.character(data$sampleID)] * pi, log = TRUE))

  expect_equal(nll, direct)
})

test_that("with a separate abundance for each sample the fit is that of efp2", {
  data <- prep_efn2_data()
  data$fsample <- factor(data$sampleID)
  m0 <- efp2(count ~ lifestage + pass12, data = data, pass = pass, id = sampleID)
  m <- efn2(count ~ lifestage + pass12, abundance = ~ 0 + fsample,
    data = data, pass = pass, id = sampleID)

  expect_equal(m$convergence, 0)
  expect_equal(coef(m), coef(m0), tolerance = 1e-5)
  expect_equal(m$P, m0$P, tolerance = 1e-5)
  # the estimated abundance is then T / P
  a <- abundance(m0)
  expect_equal(unname(m$N[as.character(a$id)]), a$N, tolerance = 1e-4)
})

test_that("estimated abundance and expected abundance", {
  data <- prep_efn2_data()
  m <- efn2(count ~ lifestage + pass12,
    abundance = ~ lifestage + siteID + offset(log(area)),
    data = data, pass = pass, id = sampleID)

  expect_equal(m$convergence, 0)
  expect_equal(m$N, m$catch + m$lambda * (1 - m$P))
  expect_equal(unname(m$catch), as.vector(tapply(data$count, data$sampleID, sum)))
  # lambda includes the offset
  samples <- data[match(names(m$lambda), data$sampleID), ]
  Xn <- model.matrix(~ lifestage + siteID, samples)
  expect_equal(unname(m$lambda), unname(exp(drop(Xn %*% m$coefficients_n)) * samples$area))

  expect_equal(dim(m$Vp_joint), rep(length(m$gam_cols) + length(m$gam_cols_n), 2))
  expect_true(all(diag(m$Vb_n) > 0))
})

test_that("random effects in both models", {
  data <- prep_efn2_data()
  data$fsample <- factor(data$sampleID)
  m <- efn2(count ~ lifestage + pass12 + s(visitID, bs = "re"),
    abundance = ~ lifestage + s(siteID, bs = "re") + s(fsample, bs = "re") +
      offset(log(area)),
    data = data, pass = pass, id = sampleID)

  expect_equal(m$convergence, 0)
  expect_equal(names(m$sigma), "s(visitID)")
  expect_equal(names(m$sigma_n), c("s(siteID)", "s(fsample)"))
  expect_true(is.matrix(m$Vp_joint))
  # capture probability is close to that from efp2
  m0 <- efp2(count ~ lifestage + pass12 + s(visitID, bs = "re"),
    data = data, pass = pass, id = sampleID)
  expect_equal(m$P, m0$P, tolerance = 0.01)

  out <- capture.output(print(m))
  capture <- out[seq_len(grep("^Abundance", out) - 1)]
  expect_true(any(grepl("sigma s\\(visitID\\)", capture)))
  expect_false(any(grepl("sigma s\\(siteID\\)", capture)))
  expect_true(any(grepl("not comparable with efp2", out)))
})

test_that("summary does not show abundance random effects for capture probability", {
  data <- prep_efn2_data()
  m <- efn2(count ~ lifestage + pass12,
    abundance = ~ lifestage + s(siteID, bs = "re"),
    data = data, pass = pass, id = sampleID)
  out <- capture.output(print(m))
  capture <- out[seq_len(grep("^Abundance", out) - 1)]
  expect_false(any(grepl("Random effects", capture)))
})

test_that("informative errors", {
  data <- prep_efn2_data()
  expect_error(
    efn2(count ~ pass12, abundance = ~ pass, data = data, pass = pass, id = sampleID),
    "same for all passes"
  )
  expect_error(
    efn2(count ~ pass12, abundance = count ~ lifestage, data = data, pass = pass, id = sampleID),
    "one sided"
  )
})

test_that("a warning when capture probability is estimated close to zero", {
  # without overdispersion in abundance, the variation in the total catches
  # is explained by the random effects for capture probability
  data <- prep_efn2_data()
  expect_warning(
    efn2(count ~ lifestage + pass12 + s(visitID, bs = "re"),
      abundance = ~ lifestage + s(siteID, bs = "re") + offset(log(area)),
      data = data, pass = pass, id = sampleID),
    "too simple"
  )
})

test_that("abundance from an efn2 fit", {
  data <- prep_efn2_data()
  data$fsample <- factor(data$sampleID)
  m <- efn2(count ~ lifestage + pass12 + s(visitID, bs = "re"),
    abundance = ~ lifestage + s(siteID, bs = "re") + s(fsample, bs = "re") +
      offset(log(area)),
    data = data, pass = pass, id = sampleID)

  a <- abundance(m, area = "area")
  ids <- as.character(a$id)
  expect_equal(a$N, unname(m$N[ids]))
  expect_equal(a$lambda, unname(m$lambda[ids]))
  expect_true(all(a$lower >= a$catch & a$lower < a$N & a$N < a$upper))
  expect_equal(a$density, a$N / a$area)

  e <- abundance(m, target = "expected")
  expect_equal(e$N, unname(m$lambda[ids]))
  expect_true(all(e$lower < e$N & e$N < e$upper))

  # the variance of the number not caught: Poisson plus the delta method,
  # with a finite difference gradient
  missed <- function(object) {
    cp <- capture_prob(object)
    sdata <- object$data[!duplicated(object$data$sampleID), ]
    Xn <- ef2:::efp_lpmatrix(list(Gsetup = object$Gsetup_n, G = object$G_n), sdata)
    lambda <- exp(drop(Xn %*% object$gam_coefficients_n) + log(sdata$area))
    lambda * (1 - cp$P)
  }
  h <- 1e-6
  grad <- function(component, cols) {
    sapply(cols, function(j) {
      up <- down <- m
      up[[component]][j] <- up[[component]][j] + h
      down[[component]][j] <- down[[component]][j] - h
      (missed(up) - missed(down)) / (2 * h)
    })
  }
  D <- cbind(grad("gam_coefficients", m$gam_cols), grad("gam_coefficients_n", m$gam_cols_n))
  M <- a$N - a$catch
  expect_equal(a$se^2, M + rowSums((D %*% m$Vp_joint) * D), tolerance = 1e-5)

  # sums over groups
  g <- abundance(m, by = c("siteID", "year"))
  samples <- data[match(ids, data$sampleID), ]
  expect_equal(g$N, as.vector(tapply(a$N, paste(samples$siteID, samples$year), sum)[paste(g$siteID, g$year)]))
  expect_true(all(g$se > 0))

  # a sample with no fish caught still has an estimate
  zero <- data[data$sampleID == data$sampleID[1], ]
  zero$count <- 0
  z <- abundance(m, newdata = zero)
  expect_equal(z$catch, 0)
  expect_gt(z$N, 0)
  expect_true(z$lower >= 0 && z$upper > z$N)
})

test_that("the negative binomial likelihood tends to the Poisson", {
  data <- prep_efn2_data()
  dat <- efn2(count ~ lifestage + pass12, abundance = ~ lifestage + offset(log(area)),
    data = data, pass = pass, id = sampleID, distribution = "negbin", fit = FALSE)
  params <- list(
    alpha = numeric(0), beta = c(0.5, 0.3, -0.2), log_sigma = numeric(0),
    alpha_n = numeric(0), beta_n = c(-1, 0.2), log_sigma_n = numeric(0),
    log_theta = 20
  )
  nb <- ef2:::ef_nll_joint(params, dat$efdat, dat$ndat)
  dat$ndat$distribution <- "poisson"
  pois <- ef2:::ef_nll_joint(params, dat$efdat, dat$ndat)
  expect_equal(nb, pois, tolerance = 1e-6)
})

test_that("the expected number not caught interpolates between efp2 and Poisson", {
  catch <- c(10, 40)
  lambda <- c(30, 50)
  P <- c(0.6, 0.9)
  expect_equal(ef2:::expected_missed(catch, lambda, P), lambda * (1 - P))
  expect_equal(ef2:::expected_missed(catch, lambda, P, theta = 1e8), lambda * (1 - P),
    tolerance = 1e-6)
  expect_equal(ef2:::expected_missed(catch, lambda, P, theta = 1e-8), catch * (1 - P) / P,
    tolerance = 1e-6)
})

test_that("negative binomial abundance", {
  data <- prep_efn2_data()
  m0 <- efp2(count ~ lifestage + pass12 + s(visitID, bs = "re"),
    data = data, pass = pass, id = sampleID)
  # the model that is degenerate with Poisson abundance
  expect_no_warning(
    m <- efn2(count ~ lifestage + pass12 + s(visitID, bs = "re"),
      abundance = ~ lifestage + s(siteID, bs = "re") + offset(log(area)),
      distribution = "negbin", data = data, pass = pass, id = sampleID)
  )
  expect_equal(m$convergence, 0)
  expect_gt(min(m$P), 0.5)
  expect_equal(m$P, m0$P, tolerance = 0.02)
  expect_equal(m$N, m$catch + ef2:::expected_missed(m$catch, m$lambda, m$P, m$theta))
  expect_true(any(grepl("Negative binomial size", capture.output(print(m)))))

  # the variance, with a finite difference gradient including log(theta)
  a <- abundance(m)
  missed <- function(object) {
    cp <- capture_prob(object)
    sdata <- object$data[!duplicated(object$data$sampleID), ]
    Xn <- ef2:::efp_lpmatrix(list(Gsetup = object$Gsetup_n, G = object$G_n), sdata)
    lambda <- exp(drop(Xn %*% object$gam_coefficients_n) + log(sdata$area))
    ef2:::expected_missed(a$catch, lambda, cp$P, object$theta)
  }
  h <- 1e-6
  grad <- function(component, cols) {
    sapply(cols, function(j) {
      up <- down <- m
      up[[component]][j] <- up[[component]][j] + h
      down[[component]][j] <- down[[component]][j] - h
      (missed(up) - missed(down)) / (2 * h)
    })
  }
  up <- down <- m
  up$theta <- m$theta * exp(h)
  down$theta <- m$theta * exp(-h)
  D <- cbind(
    grad("gam_coefficients", m$gam_cols),
    grad("gam_coefficients_n", m$gam_cols_n),
    (missed(up) - missed(down)) / (2 * h)
  )
  M <- a$N - a$catch
  b <- m$theta + a$catch
  cc <- m$theta + a$lambda * a$P
  var_process <- M + (a$lambda * (1 - a$P))^2 * b / cc^2
  expect_equal(a$se^2, var_process + rowSums((D %*% m$Vp_joint) * D), tolerance = 1e-5)
})
