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
#'estimated using the model above. For two groups, the ratio of the rates is 
#'calculated where a ratio of, say, 2.5 would mean that the rate of group 1 is 
#'2.5 times that of group 2. Note that the mean and the highest desity interval
#'for the rate ratio are calculated on the log transformed samples and then
#'transformed back to the original scale.
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
    
    # Calculating the rate_ratio mean and HDI on the log scale.
    ratio_stats <- mcmc_stats(log(as.matrix(mcmc_samples)[,"rate_ratio"]), cred_mass = cred.mass, comp_val = log(r))
    stats["rate_ratio", c("mean", "HDIlo", "HDIup")] <- exp(ratio_stats[,c("mean", "HDIlo", "HDIup")])
    stats["rate_ratio", "sd"] <- NA # Setting sd to NA, as it is clear what should be presented
                                    # when the mean was calculated on the log transformed rate_ratio
    
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
  x ~ dpois(lambda * t)
  lambda ~ dgamma(0.5, 0.00001)
  x_pred ~ dpois(lambda * t)
}"

jags_one_sample_poisson_test <- function(x, t, n.chains, n.iter, progress.bar) {  
  mcmc_samples <- run_jags(one_sample_poisson_model_string, data = list(x = x, t = t), inits = list(lambda = (x + 0.5) / t), 
                           params = c("lambda", "x_pred"), n.chains = n.chains, n.adapt = 0,
                           n.update = 100, n.iter = n.iter, thin = 1, progress.bar=progress.bar)
  mcmc_samples
}


two_sample_poisson_model_string <- "model {
  for(group_i in 1:2) {
    x[group_i] ~ dpois(lambda[group_i] * t[group_i])
    lambda[group_i] ~ dgamma(0.5, 0.00001)
    x_pred[group_i] ~ dpois(lambda[group_i] * t[group_i])
  }
  rate_diff <- lambda[1] - lambda[2]
  rate_ratio <- lambda[1] / lambda[2]
}"

jags_two_sample_poisson_test <- function(x1, t1, x2, t2, n.chains, n.iter, progress.bar) {
  data_list = list(x = c(x1, x2), t = c(t1, t2))
  init_list = list(lambda = ( c(x1, x2) + 0.5 ) / c(t1, t2) )
  mcmc_samples <- run_jags(two_sample_poisson_model_string, data = data_list, inits = init_list, 
                           params = c("lambda", "x_pred","rate_diff", "rate_ratio"), n.chains = n.chains, n.adapt = 0,
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
  cat(" ", s["lambda", "median"], "\n")
  cat(s["lambda","HDI%"],"% credible interval:\n", sep="")
  cat(" ", s["lambda",  c("HDIlo", "HDIup")], "\n")
  cat("The event rate is more than", s["lambda", "comp"] , "by a probability of", s["lambda", "%>comp"], "\n")
  cat("and less than", s["lambda", "comp"] , "by a probability of", s["lambda", "%<comp"], " .\n")
  cat("\n")
  invisible(NULL)
}

#' @export
summary.bayes_one_sample_poisson_test <- function(object, ...) {
  s <- round_or_signif(object$stats, 3)
  
  cat("  Model parameters and generated quantities\n")
  cat("lambda: the rate of the process.\n")
  cat("x_pred: predicted event count during", object$t ,"periods.\n")
  cat("\n")
  cat("  Measures\n" )
  print(s[, c("mean", "sd", "HDIlo", "HDIup", "%<comp", "%>comp")])
  cat("\n")
  cat("'HDIlo' and 'HDIup' are the limits of a ", s[1, "HDI%"] ,"% HDI credible interval.\n", sep="")
  cat("'%<comp' and '%>comp' are the probabilities of the respective parameter being\n")
  cat("smaller or larger than ", s[1, "comp"] ,".\n", sep="")
  
  cat("\n")
  cat("  Quantiles\n" )
  print(s[, c("q2.5%", "q25%", "median","q75%", "q97.5%")] )
  invisible(NULL)
}

#' @method plot bayes_one_sample_poisson_test
#' @export
plot.bayes_one_sample_poisson_test <- function(x, ...) {
  old_par <- par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0), mfcol=c(2,1))
  sample_mat <- as.matrix(x$mcmc_samples)
  lambda <- sample_mat[, "lambda"]
  # Calculating the xlim to include the comparison rate and 0 *unless* they are too 
  # far away from the samples.
  xlim = range(lambda)
  if(0 > xlim[1] - diff(range(lambda)) / 2) {
    xlim[1] <- 0
  }
  if(x$r > xlim[1] - diff(range(lambda)) / 2) {
    xlim[1] <- min(x$r, xlim[1])
  }
  if(x$r < xlim[2] + diff(range(lambda)) / 2) {
    xlim[2] <- max(x$r, xlim[2])
  }
    
  plotPost(lambda, cred_mass= x$cred_mass, comp_val=x$r, cex=1, cex.lab=1.5,
           xlim=xlim, main = "Rate of occurence", xlab=expression(lambda), show_median= TRUE)
  hist_data <- discrete_hist(sample_mat[, "x_pred"], c(0, max(sample_mat[, "x_pred"])), ylab="Probability", x_marked= x$x,
                             xlab = paste("Event count during", x$t, "periods"), main="Data w. Post. Pred.")
  #legend("topright", legend="Data", col="red",  lty=1, lwd=3)
  par(old_par)
  invisible(NULL)
}

#' @export
diagnostics.bayes_one_sample_poisson_test <- function(fit) {
  print_mcmc_info(fit$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(fit$stats, 3))
  cat("\n")
  cat("  Model parameters and generated quantities\n")
  cat("lambda: the rate of the process.\n")
  cat("x_pred: predicted event count during", fit$t ,"periods.\n")
  cat("\n")
  
  old_par <- par( mar=c(3.5,2.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  plot(fit$mcmc_samples)
  par(old_par)
  invisible(NULL)
}

#' @export
model.code.bayes_one_sample_poisson_test <- function(fit) {
  cat("### Model code for the Bayesian First Aid one sample Poisson test ###\n")
  cat("require(rjags)\n\n")
  
  cat("# Setting up the data\n")
  cat("x <-", fit$x, "\n")
  cat("t <-", fit$t, "\n")
  cat("\n")
  pretty_print_function_body(one_sample_poisson_model_code)
  invisible(NULL)
}

# Not to be run, just to be printed
one_sample_poisson_model_code <- function(x, t) {
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string
  
  # Running the model
  model <- jags.model(textConnection(model_string), data = list(x = x, t = t), n.chains = 3)
  samples <- coda.samples(model, c("lambda", "x_pred"), n.iter=5000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)  
}
one_sample_poisson_model_code <- inject_model_string(one_sample_poisson_model_code, one_sample_poisson_model_string)

### Two sample poisson test S3 methods ###


#' @export
print.bayes_two_sample_poisson_test <- function(x, ...) {
  s <- format_stats(x$stats)
  
  cat("\n")
  cat("\tBayesian Fist Aid poisson test - two sample\n")
  cat("\n")
  cat("number of events: ", paste(x$x[1], " and ", x$x[2], sep=""), ", time periods: ", paste(x$t[1], " and ", x$t[2], sep=""), "\n", sep="")
  cat("\n")
  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep="")
  cat("Group 1 rate: ", s["lambda[1]", "median"], " [", s["lambda[1]", "HDIlo"], ", ", s["lambda[1]", "HDIup"], "]","\n", sep="")
  cat("Group 2 rate: ", s["lambda[2]", "median"], " [", s["lambda[2]", "HDIlo"], ", ", s["lambda[2]", "HDIup"], "]","\n", sep="")
  cat("Rate ratio (Group 1 rate / Group 2 rate):\n              ", s["rate_ratio", "median"], " [", s["rate_ratio", "HDIlo"], ", ", s["rate_ratio", "HDIup"], "]","\n", sep="")
  
  cat("\n")
  cat("The event rate of group 1 is more than",  s["rate_ratio", "comp"], "times that of group 2 by a probability", "\n")
  cat("of", s["rate_ratio", "%>comp"], "and less than",  s["rate_ratio", "comp"], "times that of group 2 by a probability of", s["rate_ratio", "%<comp"], ".\n")
  cat("\n")
  invisible(NULL)
}

print_bayes_two_sample_poisson_test_params <- function(x) {
  cat("  Model parameters and generated quantities\n")
  cat("lambda[1]: the rate of the process of group 1.\n")
  cat("lambda[2]: the rate of the process of group 2.\n")
  cat("x_pred[1]: predicted event count of group 1 during", x$t[1] ,"periods.\n")
  cat("x_pred[2]: predicted event count of group 2 during", x$t[2] ,"periods.\n")
  cat("rate_diff: The difference lambda[1] - lambda[2].\n")
  cat("rate_ratio: The ratio lambda[1] / lambda[2].\n")
  invisible(NULL)
}

#' @export
summary.bayes_two_sample_poisson_test <- function(object, ...) {
  s <- round_or_signif(object$stats, 3)
  
  print_bayes_two_sample_poisson_test_params(object)
  cat("\n")
  cat("  Measures\n" )
  print(s[, c("mean", "sd", "HDIlo", "HDIup", "%<comp", "%>comp")])
  cat("\n")
  cat("'HDIlo' and 'HDIup' are the limits of a ", s[1, "HDI%"] ,"% HDI credible interval.\n", sep="")
  cat("'%<comp' and '%>comp' are the probabilities of the respective parameter being\n")
  cat("smaller or larger than ", s[1, "comp"] ,".\n", sep="")
  cat("For rate_ratio the mean, 'HDIlo' and 'HDIup' are calculated on the log transformed\n", sep="")
  cat("samples and then transformed back to the original scale.\n", sep="")
  
  
  cat("\n")
  cat("  Quantiles\n" )
  print(s[, c("q2.5%", "q25%", "median","q75%", "q97.5%")] )
  invisible(NULL)
}

#' @method plot bayes_two_sample_poisson_test
#' @export
plot.bayes_two_sample_poisson_test <- function(x, ...) {
  old_par <- par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0), mfcol=c(3,1))
  sample_mat <- as.matrix(x$mcmc_samples)
  lambda1 <- sample_mat[, "lambda[1]"]
  lambda2 <- sample_mat[, "lambda[2]"]
  
  xlim <- range(lambda1, lambda2)
  if(0 > xlim[1] - diff(range(lambda1, lambda2)) / 2) {
    xlim[1] <- 0
  }
  plotPost(lambda1, cred_mass= x$cred_mass, cex=1.5, cex.lab=1.8, xlim=xlim,
           main = "Rate of occurence for group 1", xlab=expression(lambda[1]), show_median= TRUE)
  plotPost(lambda2, cred_mass= x$cred_mass, cex=1.5, cex.lab=1.8, xlim=xlim,
           main = "Rate of occurence for group 2", xlab=expression(lambda[2]), show_median= TRUE)
  
  # Centers the ratio plot on 1.0.
  xlim <- 2^(max(abs(range(log2(sample_mat[, "rate_ratio"])))) * c(-1, 1))
  plotPost(sample_mat[, "rate_ratio"], cred_mass= x$cred_mass, comp_val = x$r, cex=1.5, cex.lab=1.8, xlim=xlim,
           main = "Rate ratio between group 1 and group 2", xlab=expression(lambda[1] / lambda[2]), show_median= TRUE, log_base = 2)
  
  
  par(old_par)
  invisible(NULL)
}

#' @export
diagnostics.bayes_two_sample_poisson_test <- function(fit) {
  print_mcmc_info(fit$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(fit$stats, 3))
  cat("\n")
  print_bayes_two_sample_poisson_test_params(fit)
  cat("\n")
  
  old_par <- par( mar=c(3.5,2.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  plot(fit$mcmc_samples)
  par(old_par)
  invisible(NULL)
}

#' @export
model.code.bayes_two_sample_poisson_test <- function(fit) {
  cat("### Model code for the Bayesian First Aid two sample Poisson test ###\n")
  cat("require(rjags)\n\n")
  
  cat("# Setting up the data\n")
  cat("x <-", deparse(fit$x), "\n")
  cat("t <-", deparse(fit$t), "\n")
  cat("\n")
  pretty_print_function_body(two_sample_poisson_model_code)
  invisible(NULL)
}

# Not to be run, just to be printed
two_sample_poisson_model_code <- function(x, t) {
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string
  
  # Running the model
  model <- jags.model(textConnection(model_string), data = list(x = x, t = t), n.chains = 3)
  samples <- coda.samples(model, c("lambda", "x_pred", "rate_diff", "rate_ratio"), n.iter=5000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)  
}
two_sample_poisson_model_code <- inject_model_string(two_sample_poisson_model_code, two_sample_poisson_model_string)
