context("bayes.prop.test")

test_that("bayes.prop.test accepts the same format as prop.test and returns a bayes_prop_test object", {
  
  ## These examples are from the prop.test documentation.
  # This should be a binom test since only one group is given.
  heads <- rbinom(1, size = 100, prob = .5)
  expect_that(bayes.prop.test(heads, 100, n.iter = 100) , is_a("bayes_binom_test"))
  
  smokers  <- c( 83, 90, 129, 70 )
  patients <- c( 86, 93, 136, 82 )  
  expect_that(bayes.prop.test(smokers, patients, n.iter = 100), is_a("bayes_prop_test"))
  expect_that(bayes.prop.test(smokers, patients, p = c(0.1, 0.2, 0.3, 0.4), conf.level = 0.9, correct = TRUE, n.iter = 100, 
                              alternative = "greater"), is_a("bayes_prop_test"))
})

test_that("jags_prop_test returns an mcmc.list object", {
  expect_that(jags_prop_test(c(5, 10), c(10, 10), n.iter=30), is_a("mcmc.list"))
})

test_that("bayes.prop.test with generic functions does not give an error when run", {
  smokers  <- c( 83, 90, 129, 70 )
  patients <- c( 86, 93, 136, 82 )  
  fit <- bayes.prop.test(smokers, patients, p = c(0.1, 0.2, 0.3, 0.4), conf.level = 0.9, n.iter = 100)
  expect_output(print(fit), ".")
  expect_output(summary(fit), ".")
  expect_is(summary(fit), "matrix")
  expect_output(diagnostics(fit), ".")
  expect_true({plot(fit); TRUE})
  expect_output( eval(parse(text=capture.output(model.code(fit)))) , ".")
  expect_is(as.data.frame(fit), "data.frame")
  expect_is(as.matrix(fit), "matrix")
})


