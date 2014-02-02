#' Bayesian First Aid Alternative to the Poisson Test
#' 
#' Descritions description description
#' 
#' Details details details
#' 
#' @param x
#' @param T
#' @param r
#' @param alternative
#' @param conf.level
#' @param n.iter
#' @param progress.bar The type of progress bar. Possible values are "text",
#'   "gui", and "none".
#' 
#' 
#' 
#' @return A list of class \code{bayes_one_sample_poisson_test} or
#'   \code{bayes_two_sample_poisson_test} that contains information about the
#'   analysis. It can be further inspected using the functions \code{summary},
#'   \code{plot}, \code{\link{diagnostics}} and \code{\link{model.code}}.
#' @export
bayes.poisson.test <- function (x, T = 1, r = 1, alternative = c("two.sided", "less", "greater"),
                                conf.level = 0.95, n.iter = 15000, progress.bar="none") 
{
  
  ### BEGIN code from poisson.test ###
  DNAME <- deparse(substitute(x))
  DNAME <- paste(DNAME, "time base:", deparse(substitute(T)))
  if ((l <- length(x)) != length(T)) 
    if (length(T) == 1L) 
      T <- rep(T, l)
  else stop("'x' and 'T' have incompatible length")
  xr <- round(x)
  if (any(!is.finite(x) | (x < 0)) || max(abs(x - xr)) > 1e-07) 
    stop("'x' must be finite, nonnegative, and integer")
  x <- xr
  if (any(is.na(T) | (T < 0))) 
    stop("'T' must be nonnegative")
  if ((k <- length(x)) < 1L) 
    stop("not enough data")
  if (k > 2L) 
    stop("the case k > 2 is unimplemented")
  if (!missing(r) && (length(r) > 1 || is.na(r) || r < 0)) 
    stop("'r' must be a single positive number")
  alternative <- match.arg(alternative)
  ### END code from poisson.test ###
  
  if (k == 2) {
    # two samle poison test
    mcmc_samples <- jags_two_sample_poisson_test(x[1], T[1], x[2], T[2], 
                                                 n.chains=3, n.iter= ceiling(n.iter / 3), progress.bar=progress.bar)
    bfa_object <- list(mcmc_samples = mcmc_samples, x = x, T = T)
    class(bfa_object) <- c("bayes_two_sample_poisson_test", "bayesian_first_aid")
  }
  else { # k == 1
    # one samle poison test
    mcmc_samples <- jags_one_sample_poisson_test(x, T, n.chains=3, n.iter= ceiling(n.iter / 3), progress.bar=progress.bar)
    bfa_object <- list(mcmc_samples = mcmc_samples, x = x, T = T)
    class(bfa_object) <- c("bayes_one_sample_poisson_test", "bayesian_first_aid")
  }
  bfa_object
}

one_sample_poisson_model_string <- "model {
  x ~ dpois(rate * t)
  rate ~ dgamma(0.5, 0.00001)
}"

jags_one_sample_poisson_test <- function(x, t, n.adapt= 500, n.chains=3, n.iter=5000, progress.bar="none") {  
  mcmc_samples <- run_jags(one_sample_poisson_model_string, data = list(x = x, t = t), inits = list(rate = x / t), 
                           params = c("rate"), n.chains = n.chains, n.adapt = n.adapt,
                           n.update = 0, n.iter = n.iter, thin = 1, progress.bar=progress.bar)
  mcmc_samples
}


two_sample_poisson_model_string <- "model {
  x1 ~ dpois(rate1 * t1)
  rate1 ~ dgamma(0.5, 0.00001)
  x2 ~ dpois(rate2 * t2)
  rate2 ~ dgamma(0.5, 0.00001)
  rate_diff <- rate2 - rate1
  rate_ratio <- rate1 / rate2
}"

jags_two_sample_poisson_test <- function(x1, t1, x2, t2, n.adapt= 500, n.chains=3, n.iter=5000, progress.bar="none") {
  data_list = list(x1 = x1, t1 = t1, x2 = x2, t2 = t2)
  init_list = list(rate1 = x1 / t1, rate2 = x2 / t2)
  mcmc_samples <- run_jags(two_sample_poisson_model_string, data = data_list, inits = init_list, 
                           params = c("rate1", "rate2", "rate_diff", "rate_ratio"), n.chains = n.chains, n.adapt = n.adapt,
                           n.update = 0, n.iter = n.iter, thin = 1, progress.bar=progress.bar)
  mcmc_samples
}

# ...


### One sample poisson test S3 methods ###

#' @export
print.bayes_one_sample_poisson_test <- function(x, ...) {
  cat("\n --- Bayesian first aid one sample poisson test ---\n\n")
  print(summary(x$mcmc_samples))
}

#' @export
summary.bayes_one_sample_poisson_test <- function(object, ...) {
  cat("\nSummary\n")
  print(object)
}

#' @export
plot.bayes_one_sample_poisson_test <- function(x, ...) {
  plot(x$mcmc_samples)
}

#' @export
diagnostics.bayes_one_sample_poisson_test <- function(fit) {
  cat("Not implemented\n")
  plot(fit$mcmc_samples)
}

#' @export
model.code.bayes_one_sample_poisson_test <- function(fit) {
  print(jags_one_sample_poisson_test)
}

### Two sample poisson test S3 methods ###

#' @export
print.bayes_two_sample_poisson_test <- function(x, ...) {
  cat("\n --- Bayesian first aid two sample poisson test ---\n\n")
  print(summary(x$mcmc_samples))
}

#' @export
summary.bayes_two_sample_poisson_test <- function(object, ...) {
  cat("\nSummary\n")
  print(object)
}

#' @export
plot.bayes_two_sample_poisson_test <- function(x, ...) {
  plot(x$mcmc_samples)
}

#' @export
diagnostics.bayes_two_sample_poisson_test <- function(fit) {
  cat("Not implemented\n")
  plot(fit$mcmc_samples)
}

#' @export
model.code.bayes_two_sample_poisson_test <- function(fit) {
  print(jags_two_sample_poisson_test)
}