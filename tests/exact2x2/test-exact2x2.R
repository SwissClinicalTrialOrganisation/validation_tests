# Write relevent tests for the function in here
# Consider the type of function:
#   - is it deterministic or statistic?
#   - is it worth checking for errors/warnings under particular conditions?
local_edition(3)


# reference:
# fay, m. p. (2009).
# supplementary material to confidence intervals that match
# fisher's exact or blaker's exact tests.
# biostatistics, table 2.

################################################################################
# Check p-values, odds ratios and confidence intervals against fay (2009)
################################################################################

test_that("exact2x2 agrees with published results from fay", {
  
  methods <- c("minlike", "central", "blaker")
  
  # Tremors
  
  tab <- matrix(c(1, 14, 4, 615), nrow = 2, byrow = TRUE)
  
  expected_p <- c(0.113, 0.226, 0.113)
  
  expected_ci <- list(
    c(0.42, 89.89),
    c(0.21, 119.89),
    c(0.42, 89.89)
  )
  
  for (i in seq_along(methods)) {
    
    result <- exact2x2::exact2x2(tab, tsmethod = methods[i])
    
    expect_true(
      abs(result$p.value - expected_p[i]) <= 0.0005 + 1e-10
    )
    
    expect_equal(unname(result$estimate), 10.85,
                 tolerance = 0.01)
    
    expect_equal(as.numeric(result$conf.int), expected_ci[[i]],
                 tolerance = 0.01)
  }
  
  # Vomiting / diarrhea
  
  tab <- matrix(c(5, 10, 78, 541), nrow = 2, byrow = TRUE)
  
  expected_p <- c(0.035, 0.071, 0.035)
  
  expected_ci <- list(
    c(1.11, 11.14),
    c(0.90, 11.46),
    c(1.11, 11.27)
  )
  
  for (i in seq_along(methods)) {
    
    result <- exact2x2::exact2x2(tab, tsmethod = methods[i])
    
    expect_true(
      abs(result$p.value - expected_p[i]) <= 0.0005 + 1e-10
    )
    
    expect_equal(unname(result$estimate), 3.46,
                 tolerance = 0.01)
    
    expect_equal(as.numeric(result$conf.int), expected_ci[[i]],
                 tolerance = 0.01)
  }
  
  # Abdominal pain
  
  tab <- matrix(c(4, 11, 50, 569), nrow = 2, byrow = TRUE)
  
  expected_p <- c(0.032, 0.063, 0.032)
  
  expected_ci <- list(
    c(1.17, 14.17),
    c(0.92, 14.58),
    c(1.17, 14.22)
  )
  
  for (i in seq_along(methods)) {
    
    result <- exact2x2::exact2x2(tab, tsmethod = methods[i])
    
    expect_true(
      abs(result$p.value - expected_p[i]) <= 0.0005 + 1e-10
    )
    
    expect_equal(unname(result$estimate), 4.12,
                 tolerance = 0.01)
    
    expect_equal(as.numeric(result$conf.int), expected_ci[[i]],
                 tolerance = 0.01)
  }
})


################################################################################
# Check vector input against contingency table input
################################################################################

test_that("exact2x2 gives the same results for vector and table input", {
  
  Treat <- c(rep("A", 15), rep("B", 20))
  
  Resp <- c(
    rep("Yes", 5), rep("No", 10),
    rep("Yes", 8), rep("No", 12)
  )
  
  tab <- table(Treat, Resp)
  
  methods <- c("minlike", "central", "blaker")
  
  for (method in methods) {
    
    result_vector <- exact2x2::exact2x2(
      x = Treat,
      y = Resp,
      tsmethod = method
    )
    
    result_table <- exact2x2::exact2x2(
      x = tab,
      tsmethod = method
    )
    
    expect_equal(result_vector$p.value, result_table$p.value)
    
    expect_equal(
      unname(result_vector$estimate),
      unname(result_table$estimate)
    )
    
    expect_equal(
      as.numeric(result_vector$conf.int),
      as.numeric(result_table$conf.int)
    )
  }
})
################################################################################
# Check Fisher exact p-value and odds ratio against stats::fisher.test
################################################################################

test_that("exact2x2 minlike agrees with stats::fisher.test", {
  
  #---------------------------
  # Regular case
  #---------------------------
  
  tab <- matrix(c(8, 2, 3, 7), nrow = 2, byrow = TRUE)
  
  result <- exact2x2::exact2x2(tab, tsmethod = "minlike")
  reference <- stats::fisher.test(tab)
  
  expect_equal(result$p.value, reference$p.value)
  expect_equal(unname(result$estimate), unname(reference$estimate))
  
  
  #---------------------------
  # Small sample
  #---------------------------
  
  tab <- matrix(c(1, 2, 3, 1), nrow = 2, byrow = TRUE)
  
  result <- exact2x2::exact2x2(tab, tsmethod = "minlike")
  reference <- stats::fisher.test(tab)
  
  expect_equal(result$p.value, reference$p.value)
  expect_equal(unname(result$estimate), unname(reference$estimate))
  
  
  #---------------------------
  # Zero cell
  #---------------------------
  
  tab <- matrix(c(0, 10, 5, 5), nrow = 2, byrow = TRUE)
  
  result <- exact2x2::exact2x2(tab, tsmethod = "minlike")
  reference <- stats::fisher.test(tab)
  
  expect_equal(result$p.value, reference$p.value)
  expect_equal(unname(result$estimate), unname(reference$estimate))
})



################################################################################
# Check equal proportions
################################################################################

test_that("exact2x2 minlike returns the expected results for equal proportions", {
  
  tab <- matrix(c(5, 5, 5, 5), nrow = 2, byrow = TRUE)
  
  result <- exact2x2::exact2x2(tab, tsmethod = "minlike")
  
  expect_equal(unname(result$estimate), 1, tolerance = 1e-6)
  expect_equal(result$p.value, 1)
  
  # The 95% confidence interval should contain the null odds ratio.
  
  expect_true(result$conf.int[1] <= 1)
  expect_true(result$conf.int[2] >= 1)
})


################################################################################
# Check exact p-value for a small sample
################################################################################

test_that("exact2x2 minlike returns the expected hypergeometric p-value", {
  
  tab <- matrix(c(1, 2, 3, 1), nrow = 2, byrow = TRUE)
  
  # The two-sided Fisher p-value is the sum of all
  # hypergeometric probabilities not greater than the
  # probability of the observed table.
  
  m <- sum(tab[1, ])
  n <- sum(tab[2, ])
  k <- sum(tab[, 1])
  a <- tab[1, 1]
  
  support <- max(0, k - n):min(m, k)
  probs <- stats::dhyper(support, m, n, k)
  p_obs <- stats::dhyper(a, m, n, k)
  
  expected_p <- sum(probs[probs <= p_obs * (1 + 1e-7)])
  
  result <- exact2x2::exact2x2(tab, tsmethod = "minlike")
  
  expect_equal(result$p.value, expected_p, tolerance = 1e-10)
})



################################################################################
# Check consistency between p-value and confidence interval
################################################################################

test_that("exact2x2 p-values and confidence intervals are consistent", {
  
  methods <- c("minlike", "central", "blaker")
  
  tables <- list(
    matrix(c(5, 5, 5, 5), nrow = 2, byrow = TRUE),
    matrix(c(5, 10, 78, 541), nrow = 2, byrow = TRUE),
    matrix(c(1, 14, 4, 615), nrow = 2, byrow = TRUE)
  )
  
  for (method in methods) {
    
    for (tab in tables) {
      
      result <- exact2x2::exact2x2(
        tab,
        tsmethod = method,
        conf.level = 0.95
      )
      
      p <- result$p.value
      ci <- as.numeric(result$conf.int)
      
      expect_true(is.finite(p))
      expect_true(p >= 0 && p <= 1)
      
      expect_length(ci, 2)
      expect_true(ci[1] <= ci[2])
      
      if (abs(p - 0.05) > 1e-8) {
        
        significant <- p < 0.05
        excludes_null <- ci[1] > 1 || ci[2] < 1
        
        expect_identical(significant, excludes_null)
      }
    }
  }
})

