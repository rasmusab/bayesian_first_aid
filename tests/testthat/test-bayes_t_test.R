context("bayes.t.test")

test_that("bayes.t.test.default accepts the same format as t.test.default", {
  library(datasets)
  ## The examples from the t.test documentation. Uses Student's sleep data from datasets.
  expect_that(bayes.t.test(1:10, y = c(7:20), n.iter=30), is_a("bayes_two_sample_t_test"))      # P = .00001855
  expect_that(bayes.t.test(1:10, y = c(7:20, 200), n.iter=30), is_a("bayes_two_sample_t_test")) # P = .1245    -- NOT significant anymore
  expect_that(with(sleep, bayes.t.test(extra[group == 1], extra[group == 2], n.iter=30)), is_a("bayes_two_sample_t_test"))
  expect_that(bayes.t.test(extra ~ group, data = sleep, n.iter=30), is_a("bayes_two_sample_t_test"))
})

test_that("bayes.t.test returns bayes_*_t_test object of the right class", {
  x <- rnorm(10)
  y <- rnorm(10)
  
  expect_that(bayes.t.test(x, y, n.iter=30), is_a("bayes_two_sample_t_test"))
  expect_that(bayes.t.test(x, n.iter=30), is_a("bayes_one_sample_t_test"))
  expect_that(bayes.t.test(x, y, paired=TRUE, n.iter=30), is_a("bayes_paired_t_test"))
})

test_that("the three jags_*_t_test functions returns 'mcmc.list' objects", {
  x <- rnorm(10)
  y <- rnorm(10, mean=1, sd=2)
  
  expect_that(jags_one_sample_t_test(x, n.update=30, n.iter=30), is_a("mcmc.list"))
  expect_that(jags_two_sample_t_test(x, y, n.update=30, n.iter=30), is_a("mcmc.list"))
  expect_that(jags_paired_t_test(x, y, n.update=30, n.iter=30), is_a("mcmc.list"))
})

test_that("bayes.t.test with generic functions does not give an error when run", {
  x <- rnorm(10)
  y <- rnorm(10)
  
  # Two sample
  fit <- bayes.t.test(x, y, n.iter=100)
  expect_output(print(fit), ".")
  expect_output(summary(fit), ".")
  expect_is(summary(fit), "matrix")
  expect_output(diagnostics(fit), ".")
  expect_true({plot(fit); TRUE})
  # Not implemented yet
  #expect_output( eval(parse(text=capture.output(model.code(fit)))) , ".")
  expect_is(as.data.frame(fit), "data.frame")
  expect_is(as.matrix(fit), "matrix")
  
  # One sample
  fit <- bayes.t.test(x, n.iter=100)
  expect_output(print(fit), ".")
  expect_output(summary(fit), ".")
  expect_output(diagnostics(fit), ".")
  expect_true({plot(fit); TRUE})
  # Not implemented yet
  expect_output( eval(parse(text=capture.output(model.code(fit)))) , ".")
  expect_is(as.data.frame(fit), "data.frame")
  expect_is(as.matrix(fit), "matrix")
  
  # Paired samples
  fit <- bayes.t.test(x, y, paired=TRUE, n.iter=100)
  expect_output(print(fit), ".")
  expect_output(summary(fit), ".")
  expect_output(diagnostics(fit), ".")
  expect_true({plot(fit); TRUE})
  # Not implemented yet
  expect_output( eval(parse(text=capture.output(model.code(fit)))) , ".")
  expect_is(as.data.frame(fit), "data.frame")
  expect_is(as.matrix(fit), "matrix")
  
})