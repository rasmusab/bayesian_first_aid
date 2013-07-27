context("BFA binomial test")

test_that("bfa_binom_test accepts the same format as binom.test", {
  
  ## These examples are from the binom.test documentation.
  bfa_binom_test(c(682, 243), p = 3/4)
  bfa_binom_test(682, 682 + 243, p = 3/4)
})