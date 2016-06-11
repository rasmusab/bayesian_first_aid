context("bayes.poisson.test")


test_that("the  jags_*_poisson_test functions returns an 'mcmc.list' objects", {
  expect_that(jags_one_sample_poisson_test(5, 10, n.chains = 3, n.iter=100, progress.bar = "none"), is_a("mcmc.list"))
  expect_that(jags_two_sample_poisson_test(5, 10, 20, 25, n.chains = 3, n.iter=100, progress.bar = "none"), is_a("mcmc.list"))
})

test_that("the bayes.poisson.test functions returns the correct 'bayes_*_poisson_test' objects", {
  expect_that(bayes.poisson.test(5, 10, n.iter=30), is_a("bayes_one_sample_poisson_test"))
  expect_that(bayes.poisson.test(c(5, 20), c(10, 25), n.iter=30), is_a("bayes_two_sample_poisson_test"))
})

test_that("bayes.poisson.test with generic functions does not give an error when run", {  
  # one sample
  fit <- bayes.poisson.test(5, 10, n.iter=100)
  expect_output(print(fit), ".")
  expect_output(summary(fit), ".")
  expect_is(summary(fit), "matrix")
  expect_output(diagnostics(fit), ".")
  expect_true({plot(fit); TRUE})
  # Not implemented yet
  expect_output( eval(parse(text=capture.output(model.code(fit)))) , ".")
  expect_is(as.data.frame(fit), "data.frame")
  expect_is(as.matrix(fit), "matrix")
  
  # two sample
  fit <- bayes.poisson.test(c(5, 1), c(10, 100), n.iter=100)
  expect_output(print(fit), ".")
  expect_output(summary(fit), ".")
  expect_output(diagnostics(fit), ".")
  expect_true({plot(fit); TRUE})
  # Not implemented yet
  expect_output( eval(parse(text=capture.output(model.code(fit)))) , ".")
  expect_is(as.data.frame(fit), "data.frame")
  expect_is(as.matrix(fit), "matrix")
  
})