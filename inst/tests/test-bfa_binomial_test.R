context("BFA binomial test")

test_that("bfa.binom.test accepts the same format as binom.test and returns a bfa_binom_test object", {
  
  ## These examples are from the binom.test documentation.
  expect_that(bfa.binom.test(c(682, 243), p = 3/4), is_a("bfa_binom_test"))
  expect_that(bfa.binom.test(682, 682 + 243, p = 3/4), is_a("bfa_binom_test"))
})

test_that("jags_binom_test returns an mcmc.list object", {
  expect_that(jags_binom_test(5, 10), is_a("mcmc.list"))
})

test_that("bfa.binom.test with generic function does not give an error when run", {
  fit <- bfa.binom.test(6, 11, p=0.84, alternative="two.sided")
  expect_output(print(fit), ".")
  expect_output(summary(fit), ".")
  expect_output(diagnostics(fit), ".")
  expect_true({plot(fit); TRUE})
})


