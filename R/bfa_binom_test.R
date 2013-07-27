
#' Generate MCMC samples from the Bayesian First Aid binomial model using JAGS.
#' 
#' Descritions description description
#' 
#' Details details details
#' 
#' @param x The number of successes
#' @param n The number of trials
#' 
#' @return
#' An object of type \code{mcmc.list} (defined by the \code{coda} package) that contains the MCMC samples from the model.
jags_binom_test <- function(x, n) {
  model_string <- "
    model {
      x ~ dbinom(theta, n)
      theta ~ dbeta(0.5, 0.5)
    }
  "
  
  jags_model <- jags.model(textConnection(model_string), data=list(x = x, n = n), n.adapt= 1000)
  update(jags_model, 1000, progress.bar="none")
  mcmc_samples <- coda.samples(jags_model, c("theta"), n.iter=10000)
  mcmc_samples
}

bfa_binom_test <- function (x, n, p = 0.5, alternative = c("two.sided", "less", "greater"), conf.level = 0.95) {
  ### Begin code from binom.test 
  DNAME <- deparse(substitute(x))
  xr <- round(x)
  if (any(is.na(x) | (x < 0)) || max(abs(x - xr)) > 1e-07) 
    stop("'x' must be nonnegative and integer")
  x <- xr
  if (length(x) == 2L) {
    n <- sum(x)
    x <- x[1L]
  }
  else if (length(x) == 1L) {
    nr <- round(n)
    if ((length(n) > 1L) || is.na(n) || (n < 1) || abs(n - 
                                                         nr) > 1e-07 || (x > nr)) 
      stop("'n' must be a positive integer >= 'x'")
    DNAME <- paste(DNAME, "and", deparse(substitute(n)))
    n <- nr
  }
  else stop("incorrect length of 'x'")
  if (!missing(p) && (length(p) > 1L || is.na(p) || p < 0 || 
                        p > 1)) 
    stop("'p' must be a single number between 0 and 1")
  alternative <- match.arg(alternative)
  if (!((length(conf.level) == 1L) && is.finite(conf.level) && 
          (conf.level > 0) && (conf.level < 1))) 
    stop("'conf.level' must be a single number between 0 and 1")
  ### END code from binom.test
  
  mcmc_samples <- jags_binom_test(x, n)
  structure(list(x = x, n = n, p = p, data.name = DNAME, mcmc_samples = mcmc_samples), 
            class = "bfa_binom_test")
}

print.bfa_binom_test <- function(x) {
  cat("\n --- Bayesian first aid binomial test ---\n\n")
  print(summary(x$mcmc_samples))
}

summary.bfa_binom_test <- function(object) {
  cat("\nSummary\n")
  print(object)
}

plot.bfa_binom_test <- function(x) {
  plot(x$mcmc_samples)
}

model_diagram.bfa_binom_test <- function(x) {
  print(jags_binom_test)
}

model_code.bfa_binom_test <- function(x) {
  print(jags_binom_test)
}

