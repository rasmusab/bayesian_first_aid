context("BFA poisson-test")


test_that("the  jags_*_poisson_test functions returns an 'mcmc.list' objects", {
  
  expect_that(jags_one_sample_poisson_test(5, 10), is_a("mcmc.list"))
  expect_that(jags_two_sample_poisson_test(5, 10, 20, 25), is_a("mcmc.list"))
})

test_that("the bfa.poisson.test functions returns the correct 'bfa_*_poisson_test' objects", {
  
  expect_that(bayes.poisson.test(5, 10), is_a("bfa_one_sample_poisson_test"))
  expect_that(bayes.poisson.test(c(5, 20), c(10, 25)), is_a("bfa_two_sample_poisson_test"))
})