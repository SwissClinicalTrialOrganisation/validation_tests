# Write relevant tests for the function in here
# Consider the type of function:
#   - is it deterministic or statistic?
#   - is it worth checking for errors/warnings under particular conditions?
library(testthat)
local_edition(3)
################################################################################
# Check Risk difference (RD)
################################################################################

test_that("twoby2 returns the expected probability difference", {

  # twoby2() calculates the probability difference as rd = r_a - r_b.
  
  #---------------------------
  # Regular case
  #---------------------------
  tab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  
  r_a <- 9 / 10
  r_b <- 3 / 10
  expected_rd <- r_a - r_b
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rd_row <- grep("Probability difference", rownames(result$measures))
  rd <- result$measures[rd_row, 1]
  
  expect_equal(unname(rd), expected_rd, tolerance = 1e-4)
  
  #---------------------------
  # Zero events in one group
  #--------------------------
  tab <- matrix(c(5, 51, 0, 29), nrow = 2, byrow = TRUE)
  
  r_a <- 5 / 56
  r_b <- 0 / 29
  expected_rd <- r_a - r_b
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rd_row <- grep("Probability difference", rownames(result$measures))
  rd <- result$measures[rd_row, 1]
  
  expect_equal(unname(rd), expected_rd, tolerance = 1e-4)
  
  #--------------------------
  # Zero events in both groups
  #--------------------------
  tab <- matrix(c(0, 10, 0, 10), nrow = 2, byrow = TRUE)
  
  expected_rd <- 0
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rd_row <- grep("Probability difference", rownames(result$measures))
  rd <- result$measures[rd_row, 1]
  
  expect_equal(unname(rd), expected_rd, tolerance = 1e-4)
  
})




test_that("twoby2 reverses the direction of the probability difference when groups are reversed", {
  
  # Reversing groups A and B changes
  #   rd = r_a - r_b to  rd_reversed = r_b - r_a = -rd.
  #
  # This test checks that twoby2() correctly preserves the direction
  # of the group comparison.
  
  tab_ab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  tab_ba <- matrix(c(3, 7, 9, 1), nrow = 2, byrow = TRUE)
  
  result_ab <- Epi::twoby2(tab_ab, print = FALSE, conf.level = 0.95)
  result_ba <- Epi::twoby2(tab_ba, print = FALSE, conf.level = 0.95)
  
  rd_row <- grep("Probability difference", rownames(result_ab$measures))
  rd_ab <- result_ab$measures[rd_row, 1]
  rd_ba <- result_ba$measures[rd_row, 1]
  
  expect_equal(unname(rd_ba), -unname(rd_ab), tolerance = 1e-4)
})


test_that("twoby2 returns B minus A probability difference after reversing rows and columns", {
  
  # This reproduces the table orientation used in the intended analysis.
  # The original table has groups A and B in rows and non-event/event
  # in columns. Reversing both dimensions places group B in the first
  # row and the event in the first column.
  #
  # twoby2() is therefore expected to return
  #   rd = r_b - r_a.
  
  tab <- matrix(c(7, 3, 1, 9), nrow = 2, byrow = TRUE)
  tab <- tab[2:1, 2:1]
  
  r_a <- 3 / 10
  r_b <- 9 / 10
  expected_rd <- r_b - r_a
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rd_row <- grep("Probability difference", rownames(result$measures))
  rd <- result$measures[rd_row, 1]
  
  expect_equal(unname(rd), expected_rd, tolerance = 1e-4)
})

#===============================================================================
# Check Risk difference confidence interval
#===============================================================================

test_that("twoby2 returns the expected probability difference confidence interval", {
  
  # Source:
  # Newcombe RG. Interval estimation for the difference between
  # independent proportions: comparison of eleven methods.
  # Statistics in Medicine. 1998;17:873-890, Table II.
  #
  # Epi::twoby2() obtains the confidence interval from Epi::ci.pd().
  # In Epi 2.66, ci.pd() uses method = "Nc" by default, identified in
  # its source code as method 10 from Newcombe which
  # combines Wilson score intervals without continuity correction.
  
  
  # Regular case: contrast (b), 9/10 versus 3/10
  # Published 95% CI: (0.1705, 0.8090)
  
  tab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rd_row <- grep("Probability difference", rownames(result$measures))
  rd_ci <- result$measures[rd_row, 2:3]
  
  expect_equal(unname(rd_ci), c(0.1705, 0.8090), tolerance = 0.001)
  
  
  # Zero events in one group: contrast (d), 5/56 versus 0/29
  # Published 95% CI: (-0.0381, 0.1926)
  
  tab <- matrix(c(5, 51, 0, 29), nrow = 2, byrow = TRUE)
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rd_row <- grep("Probability difference", rownames(result$measures))
  rd_ci <- result$measures[rd_row, 2:3]
  
  expect_equal(unname(rd_ci), c(-0.0381, 0.1926),  tolerance =  0.001)
  
  
  # Zero events in both groups: contrast (f), 0/10 versus 0/10
  # Published 95% CI: (-0.2775, 0.2775)
  
  tab <- matrix(c(0, 10, 0, 10), nrow = 2, byrow = TRUE)
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rd_row <- grep("Probability difference", rownames(result$measures))
  rd_ci <- result$measures[rd_row, 2:3]
  
  expect_equal(unname(rd_ci), c(-0.2775, 0.2775), tolerance = 0.001)
})  



test_that("twoby2 reverses the probability difference confidence limits when groups are reversed", {
  
  # Reversing groups transforms
  #   (lower, upper) into (-upper, -lower).
  
  tab_ab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  tab_ba <- matrix(c(3, 7, 9, 1), nrow = 2, byrow = TRUE)
  
  result_ab <- Epi::twoby2(tab_ab, print = FALSE, conf.level = 0.95)
  result_ba <- Epi::twoby2(tab_ba, print = FALSE, conf.level = 0.95)
  
  rd_row <- grep("Probability difference", rownames(result_ab$measures))
  rd_ci_ab <- result_ab$measures[rd_row, 2:3]
  rd_ci_ba <- result_ba$measures[rd_row, 2:3]
  
  expected_rd_ci_ba <- c(-rd_ci_ab[2], -rd_ci_ab[1])
  
  expect_equal(unname(rd_ci_ba), unname(expected_rd_ci_ba), tolerance = 1e-4)
})

################################################################################
# Check Relative risk (RR)
################################################################################

test_that("twoby2 returns the expected relative risk", {
  
  # The expected relative risk is calculated directly as
  #   rr = r_a / r_b.
  
  #--------------------------
  # Regular case
  #--------------------------
  tab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  
  r_a <- 9 / 10
  r_b <- 3 / 10
  expected_rr <- r_a / r_b
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rr_row <- grep("Relative Risk", rownames(result$measures))
  rr <- result$measures[rr_row, 1]
  
  expect_equal(unname(rr), expected_rr, tolerance = 1e-4)
  
  #--------------------------
  # Zero events in reference group
  #--------------------------
  tab <- matrix(c(5, 51, 0, 29), nrow = 2, byrow = TRUE)
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rr_row <- grep("Relative Risk", rownames(result$measures))
  rr <- result$measures[rr_row, 1]
  
  # The risk in the reference group is zero, therefore the
  # relative risk is infinite.
  expect_true(is.infinite(rr))
})


test_that("twoby2 returns the reciprocal relative risk when groups are reversed", {
  
  # Reversing groups changes
  #   rr = r_a / r_b
  # to
  #   rr_reversed = r_b / r_a = 1 / rr.
  
  tab_ab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  tab_ba <- matrix(c(3, 7, 9, 1), nrow = 2, byrow = TRUE)
  
  result_ab <- Epi::twoby2(tab_ab, print = FALSE, conf.level = 0.95)
  result_ba <- Epi::twoby2(tab_ba, print = FALSE, conf.level = 0.95)
  
  rr_row <- grep("Relative Risk", rownames(result_ab$measures))
  rr_ab <- result_ab$measures[rr_row, 1]
  rr_ba <- result_ba$measures[rr_row, 1]
  
  expect_equal(unname(rr_ba), 1 / unname(rr_ab), tolerance = 1e-4)
})


test_that("twoby2 returns B over A relative risk for the intended analysis orientation", {
  
  # This reproduces the table orientation used in the intended analysis.
  # Reversing both dimensions places group B in the first row and the
  # event in the first column.
  #
  # twoby2() is therefore expected to return  rr = r_b / r_a.
  
  tab <- matrix(c(7, 3, 1, 9), nrow = 2, byrow = TRUE)
  tab <- tab[2:1, 2:1]
  
  r_a <- 3 / 10
  r_b <- 9 / 10
  expected_rr <- r_b / r_a
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rr_row <- grep("Relative Risk", rownames(result$measures))
  rr <- result$measures[rr_row, 1]
  
  expect_equal(unname(rr), expected_rr, tolerance = 1e-4)
})


#===============================================================================
# Check Relative risk confidence interval
#===============================================================================

test_that("twoby2 returns the expected relative risk confidence interval", {
  
  # The 95% Wald confidence interval is calculated on the log scale:
  #
  #   se(log(rr)) = sqrt(1/a - 1/n_a + 1/c - 1/n_b)
  #   ci = exp(log(rr) +/- z_(0.975) * se(log(rr))).
  
  #--------------------------
  # Regular case
  #--------------------------
  tab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  
  a <- 9
  n_a <- 10
  c <- 3
  n_b <- 10
  
  r_a <- a / n_a
  r_b <- c / n_b
  expected_rr <- r_a / r_b
  
  se_log_rr <- sqrt(1 / a - 1 / n_a + 1 / c - 1 / n_b)
  z <- stats::qnorm(0.975)
  
  expected_rr_ci <- exp(
    log(expected_rr) + c(-1, 1) * z * se_log_rr
  )
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rr_row <- grep("Relative Risk", rownames(result$measures))
  rr_ci <- result$measures[rr_row, 2:3]
  
  expect_equal(unname(rr_ci), expected_rr_ci, tolerance = 1e-4)
  
  #--------------------------
  # Zero events in reference group
  #--------------------------
  tab <- matrix(c(5, 51, 0, 29), nrow = 2, byrow = TRUE)
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  rr_row <- grep("Relative Risk", rownames(result$measures))
  rr_ci <- result$measures[rr_row, 2:3]
  
  # The standard log-Wald confidence interval is not defined when
  # the reference group has zero events. Epi::twoby2() returns
  # NaN for the lower limit and Inf for the upper limit.
  expect_true(is.nan(rr_ci[1]))
  expect_true(is.infinite(rr_ci[2]))
})


test_that("twoby2 returns reciprocal relative risk confidence limits when groups are reversed", {
  
  # Reversing groups transforms the Wald confidence interval
  # (lower, upper) into (1 / upper, 1 / lower).
  
  tab_ab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  tab_ba <- matrix(c(3, 7, 9, 1), nrow = 2, byrow = TRUE)
  
  result_ab <- Epi::twoby2(tab_ab, print = FALSE, conf.level = 0.95)
  result_ba <- Epi::twoby2(tab_ba, print = FALSE, conf.level = 0.95)
  
  rr_row <- grep("Relative Risk", rownames(result_ab$measures))
  rr_ci_ab <- result_ab$measures[rr_row, 2:3]
  rr_ci_ba <- result_ba$measures[rr_row, 2:3]
  
  expected_rr_ci_ba <- c(1 / rr_ci_ab[2], 1 / rr_ci_ab[1])
  
  expect_equal(unname(rr_ci_ba), unname(expected_rr_ci_ba),  tolerance = 1e-4)
})


################################################################################
# Check Fisher exact p-value
################################################################################

test_that("twoby2 returns the expected Fisher exact p-value", {
  
  # Reference test: comparison with stats::fisher.test().
  #
  # This test verifies that twoby2() returns the same Fisher exact
  # p-value as the reference implementation in stats::fisher.test().
  
  #---------------------------
  # Regular case
  #---------------------------
  
  tab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  expected_p <- stats::fisher.test(tab)$p.value
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  exact_p <- result$p.value[2]
  
  expect_equal(unname(exact_p), unname(expected_p), tolerance = 1e-4)
  
  #---------------------------
  # Zero events in one group
  #---------------------------
  
  tab <- matrix(c(5, 51, 0, 29), nrow = 2, byrow = TRUE)
  expected_p <- stats::fisher.test(tab)$p.value
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  exact_p <- result$p.value[2]
  
  expect_equal(unname(exact_p),  unname(expected_p),tolerance = 1e-4)
})


################################################################################
# Check asymptotic p-value
################################################################################

test_that("twoby2 returns the expected asymptotic p-value", {
  
  # The asymptotic p-value is based on the large-sample Wald test for
  # the log odds ratio:
  #
  #   se(log(or)) = sqrt(1/a + 1/b + 1/c + 1/d)
  #   z = log(or) / se(log(or))
  #   p = 2 * P(Z >= abs(z))
  #---------------------------
  # Regular case
  #---------------------------
  tab <- matrix(c(9, 1, 3, 7), nrow = 2, byrow = TRUE)
  
  a <- 9
  b <- 1
  c <- 3
  d <- 7
  
  expected_or <- (a / b) / (c / d)
  se_log_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  z <- log(expected_or) / se_log_or
  expected_p <- 2 * stats::pnorm(-abs(z))
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  asymptotic_p <- result$p.value[1]
  
  expect_equal(unname(asymptotic_p), expected_p, tolerance = 1e-4)
  
  #---------------------------
  # Zero events in one group
  #---------------------------
  
  # With a zero cell, the uncorrected Wald test for the log odds ratio
  # is not defined. twoby2() returns NaN for the asymptotic p-value.
  
  tab <- matrix(c(5, 51, 0, 29), nrow = 2, byrow = TRUE)
  
  result <- Epi::twoby2(tab, print = FALSE, conf.level = 0.95)
  asymptotic_p <- result$p.value[1]
  
  expect_true(is.nan(asymptotic_p))
})

