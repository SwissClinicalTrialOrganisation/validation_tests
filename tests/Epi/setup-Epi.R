if(!require('Epi')) install.packages('Epi')
library(Epi)
library(testthat)
withr::defer({
  detach(package:Epi)
}, teardown_env())
