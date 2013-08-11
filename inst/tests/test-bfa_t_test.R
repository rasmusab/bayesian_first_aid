context("BFA t-test")

test_that("bfa_t_test.default accepts the same format as t.test.default", {
  library(datasets)
  ## The examples from the t.test documentation. Uses Student's sleep data from datasets.
  bfa_t_test(1:10, y = c(7:20))      # P = .00001855
  bfa_t_test(1:10, y = c(7:20, 200)) # P = .1245    -- NOT significant anymore
  with(sleep, bfa_t_test(extra[group == 1], extra[group == 2]))
})

# test_that("bfa_t_test.formula accepts the same format as t.test.default", {
#   library(datasets)
#   bfa_t_test(extra ~ group, data = sleep)
# })

test_that("bfa_t_test returns bfa_t_test* object of the right class", {
  x <- rnorm(10)
  y <- rnorm(10)
  
  test_that(bfa_t_test(x, y), is_a("bfa_two_sample_t_test"))
  test_that(bfa_t_test(x), is_a("bfa_one_sample_t_test"))
  test_that(bfa_t_test(x, y, paired=TRUE), is_a("bfa_paired_t_test"))
})

test_that("the three jags_*_t_test functions returns 'mcmc.list' objects", {
  x <- rnorm(10)
  y <- rnorm(10, mean=1, sd=2)
  
  expect_that(jags_one_sample_t_test(x), is_a("mcmc.list"))
  expect_that(jags_two_sample_t_test(x, y), is_a("mcmc.list"))
  expect_that(jags_paired_t_test(x, y), is_a("mcmc.list"))
})