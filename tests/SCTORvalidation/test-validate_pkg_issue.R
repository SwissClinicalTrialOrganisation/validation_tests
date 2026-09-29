# Write relevent tests for the function in here
# Consider the type of function:
#   - is it deterministic or statistic?
#   - is it worth checking for errors/warnings under particular conditions?
local_edition(3)

test_that("valid assessment passes", {
  res <- validate_pkg_issue(pkg_score())
  expect_true(res$score_ok)
})

test_that("invalid answer is reported", {
  res <- validate_pkg_issue(pkg_score(author = "Famous"))
  expect_false(res$score_ok)
  expect_match(res$message, "package author seems to be invalid")
})
