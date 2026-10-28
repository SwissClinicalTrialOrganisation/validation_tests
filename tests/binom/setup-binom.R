if(!require('binom')) install.packages('binom')
if(!require('dplyr')) install.packages('dplyr')
if(!require('purrr')) install.packages('purrr')

library('binom')
library('testthat')

library('dplyr', include.only = c('pull', 'lead'))
library('purrr', include.only = c('pmap', 'list_rbind'))

# Create test data
set.seed(1234)
n_test <- c(10, 1e2, 1e4, 1e8) |> rep(5) |> sort()
x_test <- (runif(length(n_test)) * n_test) |> round()
clvl_test <- ((runif(length(n_test))+1) * 0.5) |> round(digits = 3)

withr::defer({
  detach(package:binom)
  detach(package:dplyr)
  detach(package:purrr)
}, teardown_env())
