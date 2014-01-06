context("BFA cor-test")


test_that("the  jags_cor_test function returns an 'mcmc.list' objects", {
  x1 <- rnorm(10)
  x2 <- rnorm(10, mean=1, sd=2)
  
  expect_that(jags_cor_test(x1, x2, n.update=9, n.iter=9), is_a("mcmc.list"))
})

test_that("the  bfa.cor.test function returns a bfa_cor_test object", {
  x1 <- rnorm(10)
  x2 <- rnorm(10, mean=1, sd=2)
  d <- data.frame(x1 = x1, x2 = x2)
  
  
  expect_that(bayes.cor.test(x1, x2, n.iter=30), is_a("bfa_cor_test"))
  expect_that(bayes.cor.test(~ x1 + x2, data=d, n.iter=30), is_a("bfa_cor_test"))
})