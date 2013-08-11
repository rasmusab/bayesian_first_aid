context("BFA poisson-test")


test_that("the  jags_*_poisson_test functions returns an 'mcmc.list' objects", {
  
  expect_that(jags_one_sample_poisson_test(5, 10), is_a("mcmc.list"))
  expect_that(jags_two_sample_poisson_test(5, 10, 20, 25), is_a("mcmc.list"))
})