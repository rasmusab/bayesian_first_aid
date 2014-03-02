#'Bayesian First Aid Alternative to the Binomial Test
#'
#'\code{bayes.binom.test} estimates the relative frequency of success using 
#'Bayesian estimation and is intended as a replacement for 
#'\code{\link{binom.test}}.
#'
#'Given data on the number of successes and failures \code{bayes.binom.test}
#'estimates \eqn{\theta}, the relative frequency of success, assuming the
#'following model:
#'
#' \deqn{x \sim \mathrm{Binom}(\theta, n)}{x ~ Binomial(\theta, n)} 
#' \deqn{\theta \sim \mathrm{Beta}(1, 1)}{\theta ~ Beta(1, 1)}
#'
#'\if{html}{\figure{binom_diagram.svg}{options: height=250}} 
#'\if{latex}{\figure{binom_diagram.pdf}}
#' 
#' Here the prior on \eqn{\theta} is a non-informative \eqn{\mathrm{Beta}(1, 
#' 1)}{Beta(1, 1)} distribution which is identical to a \eqn{\mathrm{Uniform}(0,
#' 1)}{Uniform(0, 1)} distribution. By \code{plot}ing and looking at a
#' \code{summary} of the object returned by \code{bayes.binom.test} you can get
#' more information about the shape of the posterior and the posterior
#' predictive distribution. \code{\link{model.code}} prints out the
#' corresponding R code underlying \code{bayes.binom.test} which can be
#' copy-n-pasted into an R script and modified, for example, changing the prior
#' on \eqn{\theta}.
#' 
#' 
#' @param x number of successes, or a vector of length 2 giving the numbers of 
#'   successes and failures, respectively.
#' @param n number of trials; ignored if x has length 2.
#' @param comp.theta a fixed relative frequency of success to compare with the 
#'   estimated relative frequency of success. This argument fills a similar role
#'   as \code{p} in \code{\link{binom.test}}.
#' @param alternative ignored and is only retained in order to mantain 
#'   compatibility with \code{\link{binom.test}}.
#' @param cred.mass the amount of probability mass that will be contained in 
#'   reported credible intervals. This argument fills a similar role as 
#'   \code{conf.level} in \code{\link{binom.test}}.
#' @param n.iter The number of iterations to run the MCMC sampling.
#' @param progress.bar The type of progress bar. Possible values are "text",
#'   "gui", and "none".
#' @param p same as \code{comp.theta} and is only retained in order to mantain 
#'   compatibility with \code{\link{binom.test}}.
#' @param conf.level same as \code{cred.mass} and is only retained in order to 
#'   mantain compatibility with \code{\link{binom.test}}.
#'   
#'   
#' @return A list of class \code{bayes_binom_test} that contains information
#'   about the analysis. It can be further inspected using the functions 
#'   \code{summary}, \code{plot}, \code{\link{diagnostics}} and 
#'   \code{\link{model.code}}.
#' 
#' @examples
#' # A college dormitory recently sponsored a taste comparison between
#' # two major soft drinks. Of the 64 students who participated, 39 selected
#' # brand A, and only 25 selected brand B. Example from 
#' # http://www.elderlab.yorku.ca/~aaron/Stats2022/BinomialTest.htm
#' 
#' bayes.binom.test(x = 39, n = 64)
#' 
#' # Save the return value in order to inspect the model result further.
#' fit <- bayes.binom.test(x = 39, n = 64, cred.mass=0.8)
#' summary(fit)
#' plot(fit)
#' 
#' # MCMC diagnostics (should not be necessary for such a simple model)
#' diagnostics(fit)
#' 
#' # Print out the R code to run the model. This can be copy n' pasted into
#' # an R-script and further modified.
#' model.code(fit)
#' 
#' @export
bayes.binom.test <- function (x, n, comp.theta = 0.5, alternative = NULL, cred.mass = 0.95, n.iter=15000, progress.bar="none", p, conf.level) {
  
  if(! missing(alternative)) {
    warning("The argument 'alternative' is ignored by bayes.binom.test")
  }
  
  if(! missing(p)) {
    comp.theta <- p
  }
  
  if(! missing(conf.level)) {
    cred.mass <- conf.level
  }
  
  ### Begin code from binom.test 
  x_name <- deparse(substitute(x))
  DNAME <- x_name
  n_name <- deparse(substitute(n))
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
    DNAME <- paste(x_name, "and", n_name)
    n <- nr
  }
  else stop("incorrect length of 'x'")
  if (!missing(comp.theta) && (length(comp.theta) > 1L || is.na(comp.theta) || comp.theta < 0 || 
                                comp.theta > 1)) 
    stop("'comp.theta' or 'p' must be a single number between 0 and 1")
  if (!((length(cred.mass) == 1L) && is.finite(cred.mass) && 
          (cred.mass > 0) && (cred.mass < 1))) 
    stop("'cred.mass' or 'conf.level' must be a single number between 0 and 1")
  ### END code from binom.test
  
  mcmc_samples <- jags_binom_test(x, n, n.chains = 3, n.iter = ceiling(n.iter / 3) , progress.bar=progress.bar)
  stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = comp.theta)
  bfa_object <- list(x = x, n = n, comp_theta = comp.theta, cred_mass = cred.mass,
                     x_name = x_name, n_name = n_name, data_name = DNAME,
                     mcmc_samples = mcmc_samples, stats = stats) 
  class(bfa_object) <- c("bayes_binom_test", "bayesian_first_aid")
  bfa_object
}


binom_model_string <- "model {
  x ~ dbinom(theta, n)
  theta ~ dbeta(1, 1)
  x_pred ~ dbinom(theta, n)
}"

jags_binom_test <- function(x, n, n.chains=3, n.iter=5000, progress.bar="none") {
  mcmc_samples <- run_jags(binom_model_string, data = list(x = x, n = n), inits = list(theta = (x + 1) / (n + 2)), 
                           params = c("theta", "x_pred"), n.chains = n.chains, n.adapt = 0,
                           n.update = 0, n.iter = n.iter, thin = 1, progress.bar=progress.bar)
  mcmc_samples
}

### binom test S3 methods ###

#' @export
print.bayes_binom_test <- function(x, ...) {
  
  s <- format_stats(x$stats)["theta",]
  
  cat("\n")
  cat("\tBayesian First Aid binomial test\n")
  cat("\n")
  cat("data: ", x$data_name, "\n", sep="")
  cat("number of successes = ", x$x,", number of trials = ", x$n, "\n", sep="")
  cat("Estimated relative frequency of success:\n")
  cat(" ", s["mean"], "\n")
  cat(s["HDI%"],"% credible interval:\n", sep="")
  cat(" ", s[ c("HDIlo", "HDIup")], "\n")
  cat("The relative frequency of success is more than", s["comp"] , "by a probability of", s["%>comp"], "\n")
  cat("and less than", s["comp"] , "by a probability of", s["%<comp"], "\n")
  cat("\n")
}


#' @export
summary.bayes_binom_test <- function(object, ...) {
  s <- round(object$stats, 3)
  
  cat("  Data\n")
  cat("number of successes = ", object$x,", number of trials = ", object$n, "\n", sep="")
  cat("\n")
  
  cat("  Model parameters and generated quantities\n")
  cat("theta: the relative frequency of success\n")
  cat("x_pred: predicted number of successes in a replication\n")
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
  # Lägga till HDI
  # 'HDIlo' and 'HDIup' are the limits of a 95% HDI credible interval.
}


#' @export
plot.bayes_binom_test <- function(x, ...) {
  layout(matrix(c(1,2), nrow=2, ncol=1 , byrow=FALSE) )
  old_par <- par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  sample_mat <- as.matrix(x$mcmc_samples)
  plotPost(sample_mat[, "theta"], cred_mass= x$cred_mass, comp_val=x$comp_theta, xlim=c(0, 1), cex=1, cex.lab=1.5,
           main = "Relative Frequency of Success", xlab=expression(theta))
  hist_data <- discrete_hist(sample_mat[, "x_pred"], c(0, x$n), ylab="Probability", x_marked= x$x,
                             xlab = "Number of sucesses",main="Data w. Post. Pred.")
  #legend("topright", legend="Data", col="red",  lty=1, lwd=3)
  par(old_par)
}


#' @export
diagnostics.bayes_binom_test <- function(fit) {
  print_mcmc_info(fit$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(fit$stats, 3))
  cat("\n")
  
  cat("  Model parameters and generated quantities\n")
  cat("theta: The relative frequency of success\n")
  cat("x_pred: Predicted number of successes in a replication\n")
  old_par <- par( mar=c(3.5,2.5,2.5,0.6) , mgp=c(2.25,0.7,0) )
  plot(fit$mcmc_samples)
  par(old_par)
}


#' @export
model.code.bayes_binom_test <- function(fit) {
  cat("### Model code for the Bayesian First Aid alternative to the binomial test ###\n")
  cat("require(rjags)\n\n")
  
  cat("# Setting up the data\n")
  cat("x <-", fit$x, "\n")
  cat("n <-", fit$n, "\n")
  cat("\n")
  pretty_print_function_body(binom_model_code)
}

# Not to be run, just to be printed
binom_model_code <- function() {
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string
  
  # Running the model
  model <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                      n.chains = 3, n.adapt=1000)
  samples <- coda.samples(model, c("theta", "x_pred"), n.iter=5000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)  
}
binom_model_code <- inject_model_string(binom_model_code, binom_model_string)
