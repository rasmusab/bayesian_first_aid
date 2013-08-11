context("BFA cor-test")


test_that("the  jags_cor_test function returns an 'mcmc.list' objects", {
  x1 <- rnorm(10)
  x2 <- rnorm(10, mean=1, sd=2)
  
  expect_that(jags_cor_test(x1, x2), is_a("mcmc.list"))
})