#'Bayesian First Aid Alternative to the Poisson Test
#'
#'\code{poisson.test} estimates the rate parameter for one or two groups of 
#'counts using Bayesian estimation assuming a Poisson distribution. This 
#'procedure is intended as a replacement for \code{\link{prop.test}}.
#'
#'Give data on the number of counts {\eqn{x}} during {\eqn{T}} periods (e.g., 
#'days, hours, years, etc.) the underlying rate \eqn{\lambda} is estimated 
#'assuming the following model:
#'
#'\deqn{x \sim \mathrm{Poisson}(\lambda \cdot T}{x ~ Poisson(\lambda·T)} 
#'\deqn{\lambda \sim \mathrm{Gamma}(0.5, 0.00001)}{\lambda ~ Gamma(0.5, 
#'0.00001)}
#'
#'Here the Gamma prior on \eqn{\lambda} is an approximation of Jeffrey's' prior 
#'for this model. In the case of two groups, both rate parameters are separately
#'estimated using the model above.
#'
#'@param x number of events. A vector of length one or two.
#'@param T time base for event count. A vector of length one or two.
#'@param r hypothesized rate or rate ratio
#'@param alternative ignored and is only retained in order to mantain 
#'  compatibility with \code{\link{poisson.test}}.
#'@param cred.mass the amount of probability mass that will be contained in 
#'  reported credible intervals. This argument fills a similar role as 
#'  \code{conf.level} in \code{\link{poisson.test}}.
#'@param n.iter
#'@param progress.bar The type of progress bar. Possible values are "text", 
#'  "gui", and "none".
#'@param conf.level same as \code{cred.mass} and is only retained in order to 
#'  mantain compatibility with \code{\link{poisson.test}}.
#'  
#'  
#'  
#'@return A list of class \code{bayes_one_sample_poisson_test} or 
#'  \code{bayes_two_sample_poisson_test} that contains information about the 
#'  analysis. It can be further inspected using the functions \code{summary}, 
#'  \code{plot}, \code{\link{diagnostics}} and \code{\link{model.code}}.
#'  
#' @examples
#' 
#' # Data from Boice, J. D., & Monson, R. R. (1977). 
#' # Breast cancer in women after repeated fluoroscopic examinations of the chest.
#' # Journal of the National Cancer Institute, 59(3), 823-832.
#' 
#' # 41 cases of breast cancer during 28,010 person-years in the treatment group
#' # of women receiving X-ray fluoroscopy and 15 cases of breast cancer during 
#' # 19 017 person-years in the control group.
#' 
#' no_cases <- c(41, 15)
#' no_years <- c(28010, 19017)
#' 
#' bayes.poisson.test(no_cases, no_years)
#' 
#' # Save the return value in order to inspect the model result further.
#' fit <- bayes.poisson.test(no_cases, no_years)
#' summary(fit)
#' plot(fit)
#' 
#' # MCMC diagnostics (should not be necessary for such a simple model)
#' diagnostics(fit)
#' 
#' # Print out the R code to run the model. This can be copy'n'pasted into
#' # an R-script and further modified.
#' model.code(fit)
#' 
#'  
#'@export
bayes.poisson.test <- function (x, T = 1, r = 1, alternative = c("two.sided", "less", "greater"),
                                cred.mass = 0.95, n.iter = 15000, progress.bar="none", conf.level) 
{
  
  if(! missing(alternative)) {
    warning("The argument 'alternative' is ignored by bayes.poisson.test")
  }
  
  if(! missing(conf.level)) {
    cred.mass <- conf.level
  }
  
  x_name <- deparse(substitute(x))
  t_name <- deparse(substitute(T))
  
  ### BEGIN Slightly modified code from  poisson.test ###
  DNAME <- deparse(substitute(x))
  DNAME <- paste(DNAME, ", time base: ", deparse(substitute(T)), sep = "")
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
    stop(paste("The case for more than two groups is unimplemented.",
               "Run the model with just two groups and use the model.code function",
               "to print the code that implements the model. This code can then be",
               "easily modified to work with many groups!", sep="\n"))
  if (!missing(r) && (length(r) > 1 || is.na(r) || r < 0)) 
    stop("'r' must be a single positive number")
  alternative <- match.arg(alternative)
  ### END code from poisson.test ###

  if (k == 2) {
    # two samle poison test
    mcmc_samples <- jags_two_sample_poisson_test(x[1], T[1], x[2], T[2], 
                                                 n.chains=3, n.iter= ceiling(n.iter / 3), progress.bar=progress.bar)
    stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = r)
    bfa_object <- list(x = x, t = T, r = r, cred_mass = cred.mass, x_name = x_name, t_name = t_name,
                       data_name = DNAME, mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- c("bayes_two_sample_poisson_test", "bayesian_first_aid")
  }
  else { # k == 1
    # one samle poison test
    mcmc_samples <- jags_one_sample_poisson_test(x, T, n.chains=3, n.iter= ceiling(n.iter / 3), progress.bar=progress.bar)
    stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = r)
    bfa_object <- list(x = x, t = T, r = r, cred_mass = cred.mass, x_name = x_name, t_name = t_name,
                       data_name = DNAME, mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- c("bayes_one_sample_poisson_test", "bayesian_first_aid")
  }
  bfa_object
}

one_sample_poisson_model_string <- "model {
  x ~ dpois(rate * t)
  rate ~ dgamma(0.5, 0.00001)
  x_rep ~ dpois(rate * t)
}"

jags_one_sample_poisson_test <- function(x, t, n.chains, n.iter, progress.bar) {  
  mcmc_samples <- run_jags(one_sample_poisson_model_string, data = list(x = x, t = t), inits = list(rate = (x + 0.5) / t), 
                           params = c("rate", "x_rep"), n.chains = n.chains, n.adapt = 0,
                           n.update = 100, n.iter = n.iter, thin = 1, progress.bar=progress.bar)
  mcmc_samples
}


two_sample_poisson_model_string <- "model {
  for(group_i in 1:2) {
    x[group_i] ~ dpois(rate[group_i] * t[group_i])
    rate[group_i] ~ dgamma(0.5, 0.00001)
    x_rep[group_i] ~ dpois(rate[group_i] * t[group_i])
  }
  rate_diff <- rate[1] - rate[2]
  rate_ratio <- rate[1] / rate[2]

}"

jags_two_sample_poisson_test <- function(x1, t1, x2, t2, n.chains, n.iter, progress.bar) {
  data_list = list(x = c(x1, x2), t = c(t1, t2))
  init_list = list(rate = ( c(x1, x2) + 0.5 ) / c(t1, t2) )
  mcmc_samples <- run_jags(two_sample_poisson_model_string, data = data_list, inits = init_list, 
                           params = c("rate", "x_rep","rate_diff", "rate_ratio"), n.chains = n.chains, n.adapt = 0,
                           n.update = 100, n.iter = n.iter, thin = 1, progress.bar=progress.bar)
  mcmc_samples
}

# ...


### One sample poisson test S3 methods ###

#' @export
print.bayes_one_sample_poisson_test <- function(x, ...) {
  s <- format_stats(x$stats)
  
  cat("\n")
  cat("\tBayesian Fist Aid poisson test - one sample\n")
  cat("\n")
  cat("number of events: ", x$x, ", time periods: ", x$t, sep="")
  cat("\n")
  cat("Estimated event rate:\n")
  cat(" ", s["rate", "median"], "\n")
  cat(s["rate","HDI%"],"% credible interval:\n", sep="")
  cat(" ", s["rate",  c("HDIlo", "HDIup")], "\n")
  cat("The event rate is more than", s["rate", "comp"] , "by a probability of", s["rate", "%>comp"], "\n")
  cat("and less than", s["rate", "comp"] , "by a probability of", s["rate", "%<comp"], " .\n")
  cat("\n")
  invisible(NULL)
}

#' @export
summary.bayes_one_sample_poisson_test <- function(object, ...) {
  cat("\nSummary\n")
  print(object)
  invisible(NULL)
}

#' @method plot bayes_one_sample_poisson_test
#' @export
plot.bayes_one_sample_poisson_test <- function(x, ...) {
  plot(x$mcmc_samples)
  invisible(NULL)
}

#' @export
diagnostics.bayes_one_sample_poisson_test <- function(fit) {
  cat("Not implemented\n")
  plot(fit$mcmc_samples)
  invisible(NULL)
}

#' @export
model.code.bayes_one_sample_poisson_test <- function(fit) {
  print(jags_one_sample_poisson_test)
  invisible(NULL)
}

### Two sample poisson test S3 methods ###


#' @export
print.bayes_two_sample_poisson_test <- function(x, ...) {
  s <- format_stats(x$stats)
  
  cat("\n")
  cat("\tBayesian Fist Aid poisson test - two sample\n")
  cat("\n")
  cat("number of events: ", paste(x$x[1], " and ", x$x[2], sep=""), ", time periods: ", paste(x$t[1], " and ", x$t[2], sep=""), "\n", sep="")
  cat("\n")
  cat("  Estimated event rate [", s[1, "HDI%"] ,"% credible interval]\n", sep="")
  cat("Group 1: ", s["rate[1]", "median"], " [", s["rate[1]", "HDIlo"], ", ", s["rate[1]", "HDIup"], "]","\n", sep="")
  cat("Group 2: ", s["rate[2]", "median"], " [", s["rate[2]", "HDIlo"], ", ", s["rate[2]", "HDIup"], "]","\n", sep="")
  cat("\n")
  cat("The event rate of group 1 is more than",  s["rate_ratio", "comp"], "times that of group 2 by a probability", "\n")
  cat("of", s["rate_ratio", "%>comp"], "and less than",  s["rate_ratio", "comp"], "times that of group 2 by a probability of", s["rate_ratio", "%<comp"], " .\n")
  cat("\n")
  invisible(NULL)
}

#' @export
summary.bayes_two_sample_poisson_test <- function(object, ...) {
  cat("\nSummary\n")
  print(object)
  invisible(NULL)
}

#' @method plot bayes_two_sample_poisson_test
#' @export
plot.bayes_two_sample_poisson_test <- function(x, ...) {
  plot(x$mcmc_samples)
  invisible(NULL)
}

#' @export
diagnostics.bayes_two_sample_poisson_test <- function(fit) {
  cat("Not implemented\n")
  plot(fit$mcmc_samples)
  invisible(NULL)
}

#' @export
model.code.bayes_two_sample_poisson_test <- function(fit) {
  print(jags_two_sample_poisson_test)
  invisible(NULL)
}