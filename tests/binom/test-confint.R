local_edition(3)

# Source: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
# (Section 'Clopper-Pearson interval')
confint_exact <- function(n, x, clvl) {
  a <- 1 - clvl
  m <- x / n
  lowI = qbeta(a/2, x, n-x+1)
  uppI = qbeta(1-a/2, x+1, n-x)
  data.frame(m, lowI, uppI)
}

test_that("Test 'exact' method", {
  res <- binom::binom.confint(
      x = x_test,
      n = n_test,
      conf.level = clvl_test,
      methods = c('exact'))
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  results_exp <- confint_exact(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp) 
  expect_equal(lowI_res, lowI_exp) 
  expect_equal(uppI_res, uppI_exp) 
})


# Source: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
# (Section 'Agresti–Coull interval')
confint_ac <- function(n, x, clvl) {
  a <- 1 - clvl
  z_a <- qnorm(1 - a*0.5)
  n_tilde <- n + z_a*z_a
  x_tilde <- x + z_a*z_a*0.5
  p_tilde <- x_tilde / n_tilde
  lowI <- p_tilde - z_a * sqrt(p_tilde/n_tilde * (1-p_tilde))
  uppI <- p_tilde + z_a * sqrt(p_tilde/n_tilde * (1-p_tilde))
  m <- x / n
  data.frame(m, lowI, uppI)
}

test_that("Test Agresti-Coull method", {
  res <- binom::binom.confint(
      x = x_test,
      n = n_test,
      conf.level = clvl_test,
      methods = c('agresti-coull'))
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  
  results_exp <- confint_ac(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp) 
  expect_equal(lowI_res, lowI_exp) 
  expect_equal(uppI_res, uppI_exp) 
})


# Source: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
# (Section 'Problems with using a normal approximation or "Wald interval"')
confint_wald <- function(n, x, clvl) {
  a <- 1 - clvl
  z_a <- qnorm(1 - a/2)
  p_hat <- x / n
  lowI <- p_hat - z_a/sqrt(n) * sqrt(p_hat*(1-p_hat)) 
  uppI <- p_hat + z_a/sqrt(n) * sqrt(p_hat*(1-p_hat)) 
  m <- x / n
  data.frame(m, lowI, uppI)
}

test_that("Test Wald method", {
  res <- binom::binom.confint(
      x = x_test,
      n = n_test,
      conf.level = clvl_test,
      methods = c('asymptotic'))
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  results_exp <- confint_wald(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp) 
  expect_equal(lowI_res, lowI_exp) 
  expect_equal(uppI_res, uppI_exp) 
})


# Source: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
# (Section 'Wilson score interval')
confint_wilson <- function(n, x, clvl) {
  a <- 1 - clvl
  z_a <- qnorm(1 - a/2)
  p_hat <- x / n
  lowI <- 1/(1+z_a*z_a / n) *
    (p_hat + z_a*z_a/(2*n) - z_a/(2*n) * sqrt(4*n*p_hat*(1-p_hat)+z_a*z_a))
  uppI <- 1/(1+z_a*z_a / n) *
    (p_hat + z_a*z_a/(2*n) + z_a/(2*n) * sqrt(4*n*p_hat*(1-p_hat)+z_a*z_a))
  m <- x / n
  data.frame(m, lowI, uppI)
}

test_that("Test Wilson method", {
  res <- binom::binom.confint(
      x = x_test,
      n = n_test,
      conf.level = clvl_test,
      methods = c('wilson'))
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  results_exp <- confint_wilson(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp) 
  expect_equal(lowI_res, lowI_exp) 
  expect_equal(uppI_res, uppI_exp) 
})

#' allenr: This still fails.   
# Implementation according to package's Help page
confint_prop.test <- function(n, x, clvl) {
  extract_confint <- function(nn, xx, ll) {
    res <- prop.test(xx, nn, conf.level = ll)$conf.int
    data.frame(m=xx/nn, lowI=res[1], uppI=res[2])
  }
  list(n, x, clvl) |>
    pmap(extract_confint) |>
    list_rbind()
}

test_that("Test 'prop.test' method", {
  single_binom.confint <- function(x, n, clvl) 
    binom::binom.confint(x, n, conf.level = clvl, method = c('prop.test'))
  
  res <- list(x_test, n_test, clvl_test) |>
    pmap(single_binom.confint)|>
    list_rbind()
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  results_exp <- confint_prop.test(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp) 
  expect_equal(lowI_res, lowI_exp) 
  expect_equal(uppI_res, uppI_exp) 
})


# Source: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
# (Section 'Jeffreys interval')
confint_bayes <- function(n, x, clvl) {
  confint_bayes_single <- function(nn, xx, ll, prio_s1=0.5, prio_s2=0.5) {
    # calculate posterior shape parameters
    post_s1 <- xx + prio_s1 
    post_s2 <- nn - xx + prio_s2 
    # calculate posterior's mean
    m <- post_s1/(post_s1 + post_s2)
    
    # dirty solution to find the highest probability density (hpd) interval
    # within the specified clvl: use qbeta to find the value p such that
    # qbeta(p+clvl) - qbeta(p) is minimized. Then [qbeta(p), qbeta(p+clvl)]
    # is the hpd confidence interval.
    nsteps <- 1e5
    psteps <- c(1:nsteps)/nsteps
    qs <- qbeta(psteps, post_s1, post_s2)
    lvlgap <- round(ll * nsteps)
    amin <- (lead(qs, lvlgap) - qs) |> which.min()
    data.frame(m, lowI=qs[amin], uppI=qs[amin+lvlgap])
  }
  list(n, x, clvl) |>
    pmap(confint_bayes_single) |>
    list_rbind()
}

test_that("Test 'bayes' method", {suppressWarnings({
  single_confint <- function(x, n, clvl) 
    binom::binom.confint(x, n, conf.level = clvl, methods = c('bayes'))
  
  res <- list(x_test, n_test, clvl_test) |>
    pmap(single_confint) |>
    list_rbind()
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  results_exp <- confint_bayes(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp, tolerance = 1e-5) 
  expect_equal(lowI_res, lowI_exp, tolerance = 1e-5) 
  expect_equal(uppI_res, uppI_exp, tolerance = 1e-5) 
})})


# implementation according to
# https://github.com/cran/binom/blob/master/inst/doc/binom.pdf
confint_logit <- function(n, x, clvl) {
  lowI <- numeric(length(n))
  uppI <- numeric(length(n))
  
  # Filter non-degenerative cases (n=x or x=0, otherwise division by 0)
  i <- !(n==x)&!(x==0)
  xx <- x[i]
  nn <- n[i]
  ll <- clvl[i]
  
  a <- 1-ll
  b <- log(xx/(nn-xx))
  p_hat = xx/nn
  c <- qnorm(1-a/2)/sqrt(nn*p_hat*(1-p_hat))
  
  lowI[i] <- 1-1/(1+exp(b-c))
  uppI[i] <- 1-1/(1+exp(b+c))
  
  # In the degenerative case, we use the exact method
  lowI[!i] <- confint_exact(n[!i], x[!i], clvl[!i]) |> pull('lowI')
  uppI[!i] <- confint_exact(n[!i], x[!i], clvl[!i]) |> pull('uppI')
  
  data.frame(m=x/n, lowI, uppI)
}

test_that("Test 'logit' method", {
  res <- binom::binom.confint(
      x = x_test,
      n = n_test,
      conf.level = clvl_test,
      methods = c('logit'))
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  results_exp <- confint_logit(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp) 
  expect_equal(lowI_res, lowI_exp) 
  expect_equal(uppI_res, uppI_exp) 
})


# implementation according to
# https://github.com/cran/binom/blob/master/inst/doc/binom.pdf
confint_cloglog <- function(n, x, clvl) {
  lowI <- numeric(length(n))
  uppI <- numeric(length(n))
  
  # Filter non-degenerative cases (n=x or x=0, ow division by 0)
  i <- !(n==x)&!(x==0)
  xx <- x[i]
  nn <- n[i]
  ll <- clvl[i]
  
  cloglog <- function(p) log(-log(p))
  cloglog_inv <- function(z) exp(-exp(z))
  
  p_hat <- xx/nn
  a <- 1-ll
  b <- cloglog(p_hat)
  sd <- sqrt((1-p_hat)/(n*p_hat*log(p_hat)^2))
  
  m_lowI <- b + qnorm(1-a/2)*sd
  m_uppI <- b - qnorm(1-a/2)*sd
  lowI[i] <- cloglog_inv(m_lowI)
  uppI[i] <- cloglog_inv(m_uppI)
  
  # In the degenerative case, we use the exact method
  lowI[!i] <- confint_exact(n[!i], x[!i], clvl[!i]) |> pull('lowI')
  uppI[!i] <- confint_exact(n[!i], x[!i], clvl[!i]) |> pull('uppI')
  
  data.frame(m=x/n, lowI, uppI)
}

test_that("Test 'cloglog' method", {
  res <- binom::binom.confint(
      x = x_test,
      n = n_test,
      conf.level = clvl_test,
      methods = c('cloglog'))
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  results_exp <- confint_cloglog(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp) 
  expect_equal(lowI_res, lowI_exp) 
  expect_equal(uppI_res, uppI_exp) 
})


# implementation according to
# https://github.com/cran/binom/blob/master/inst/doc/binom.pdf
confint_probit <- function(n, x, clvl) {
  lowI <- numeric(length(n))
  uppI <- numeric(length(n))
  
  # Filter non-degenerative cases (n=x or x=0, ow division by 0)
  i <- !(n==x)&!(x==0)
  xx <- x[i]
  nn <- n[i]
  ll <- clvl[i]
  
  probit <- qnorm
  probit_inv <- pnorm
  a <- 1-ll
  p_hat <- xx/nn
  
  b <- probit(p_hat)
  sd <- sqrt(p_hat*(1-p_hat)/(n*dnorm(b)^2))
  
  p_lowI <- b - qnorm(1-a/2)*sd
  p_uppI <- b + qnorm(1-a/2)*sd
  lowI[i] <- probit_inv(p_lowI)
  uppI[i] <- probit_inv(p_uppI)
  
  # In the degenerative case, we use the exact method
  lowI[!i] <- confint_exact(n[!i], x[!i], clvl[!i]) |> pull('lowI')
  uppI[!i] <- confint_exact(n[!i], x[!i], clvl[!i]) |> pull('uppI')
  
  data.frame(m=x/n, lowI, uppI)
}

test_that("Test 'probit' method", {
  res <- binom::binom.confint(
      x = x_test,
      n = n_test,
      conf.level = clvl_test,
      methods = c('probit'))
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  results_exp <- confint_probit(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp) 
  expect_equal(lowI_res, lowI_exp) 
  expect_equal(uppI_res, uppI_exp) 
})

confint_profile <- function(n, x, clvl) {
  # Calculation for a single (nn, xx, ll) triple
  confint_profile_internal <- function(nn, xx, ll) {
    if (xx==0 || xx==nn) return(confint_exact(nn, xx, ll))
    
    p_hat <- xx / nn
    # log-likelihood function (up to a constant)
    loglik <- function(p) {
      if (p <= 0 || p >= 1) return(-inf)
      xx * log(p) + (nn - xx) * log(1 - p)
    }
    
    # likelihood ratio statistic (profile deviance)
    Lambda <- function(p) {
      2 * (loglik(p_hat) - loglik(p))
    }
    
    # Grid of p values
    p_grid <- seq(0.00001, 0.99999, length.out = 100000)
    lambda_values <- p_grid |> sapply(Lambda)
    
    # Find confidence interval
    threshold <- qchisq(ll, df = 1)
    inside_CI <- which(lambda_values <= threshold)
    lowI <- min(p_grid[inside_CI])
    uppI <- max(p_grid[inside_CI])
    
    data.frame(m=p_hat, lowI, uppI)
  }
   
  # Iterate over all triples and combine in single df
  list(n, x, clvl) |>
    pmap(confint_profile_internal) |>
    list_rbind()
}

test_that("Test 'profile' method", {
  single_confint <- function(x, n, clvl) 
    binom::binom.confint(x, n,
                         conf.level = clvl,
                         methods = c('profile'),
                         bayes = FALSE,
                         maxsteps = 100)
  
  res <- list(x_test, n_test, clvl_test) |>
    pmap(single_confint) |>
    list_rbind()
  
  mean_res <- res |> pull('mean')
  lowI_res <- res |> pull('lower')
  uppI_res <- res |> pull('upper')
  
  results_exp <- confint_profile(n_test, x_test, clvl_test)
  
  mean_exp <- results_exp |> pull('m')
  lowI_exp <- results_exp |> pull('lowI')
  uppI_exp <- results_exp |> pull('uppI')
  
  expect_equal(mean_res, mean_exp, tolerance = 1e-4) 
  expect_equal(lowI_res, lowI_exp, tolerance = 1e-4) 
  expect_equal(uppI_res, uppI_exp, tolerance = 1e-4) 
})
