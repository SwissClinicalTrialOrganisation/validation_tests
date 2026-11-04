if(!require('exact2x2')) install.packages('exact2x2')
library(exact2x2)
library(testthat)
withr::defer({
  detach(package:exact2x2)
}, teardown_env())
