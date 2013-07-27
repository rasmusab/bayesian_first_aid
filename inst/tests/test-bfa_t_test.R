context("BFA t-test")

test_that("bfa_t_test accepts the same format as t.test", {
  library(datasets)
  
  ## The examples from the t.test documentation. Uses Student's sleep data from datasets.
  fba_t_test(1:10, y = c(7:20))      # P = .00001855
  fba_t_test(1:10, y = c(7:20, 200)) # P = .1245    -- NOT significant anymore
  with(sleep, fba_t_test(extra[group == 1], extra[group == 2]))
  fba_t_test(extra ~ group, data = sleep)
})