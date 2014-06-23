context("bayes.binom.test")

test_that("bayes.binom.test accepts the same format as binom.test and returns a bayes_binom_test object", {
  
  ## These examples are from the binom.test documentation.
  expect_that(bayes.binom.test(c(682, 243), p = 3/4, n.iter=30), is_a("bayes_binom_test"))
  expect_that(bayes.binom.test(682, 682 + 243, p = 3/4, n.iter=30), is_a("bayes_binom_test"))
})

test_that("jags_binom_test returns an mcmc.list object", {
  expect_that(jags_binom_test(5, 10, n.iter=30), is_a("mcmc.list"))
})

test_that("bayes.binom.test with generic functions does not give an error when run", {
  fit <- bayes.binom.test(6, 11, p=0.84, alternative="two.sided", n.iter=30)
  expect_output(print(fit), ".")
  expect_output(summary(fit), ".")
  expect_output(diagnostics(fit), ".")
  expect_true({plot(fit); TRUE})
  expect_output( eval(parse(text=capture.output(model.code(fit)))) , ".")
  expect_is(as.data.frame(fit), "data.frame")
  expect_is(as.matrix(fit), "matrix")
})


