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
jags_binom_test <- function(x, n, n.chains=3, n.iter=1000) {
  model_string <- "
    model {
      x ~ dbinom(theta, n)
      theta ~ dbeta(1, 1)
      x_pred ~ dbinom(theta, n)
    }
  "
  mcmc_samples <- run_jags(model_string, data = list(x = x, n = n), inits = list(theta = (x + 1) / (n + 2)), 
                           params = c("theta", "x_pred"), n.chains = n.chains, n.adapt = 0,
                           n.update = 0, n.iter = n.iter, thin = 1)
  mcmc_samples
}

bfa.binom.test <- function (x, n, p = 0.5, alternative = c("two.sided", "less", "greater"), conf.level = 0.95) {
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
  stats <- mcmc_stats(mcmc_samples, cred_mass = conf.level, comp_val = p)
  bfa_object <- list(x = x, n = n, p = p, conf.level = conf.level,
                 data_name = DNAME, mcmc_samples = mcmc_samples, stats = stats) 
  class(bfa_object) <- "bfa_binom_test"
  bfa_object
}

### binom test S3 methods ###

print.bfa_binom_test <- function(x) {
  
  s <- round(x$stats["theta",], 3)
  
  cat("\n")
  cat("\tBayesian first aid binomial test\n")
  cat("\n")
  cat("data: ", x$data_name, "\n", sep="")
  cat("number of successes = ", x$x,", number of trials = ", x$n, "\n", sep="")
  cat("Estimated relative frequency of success:\n")
  cat(" ", s["mean"], "\n")
  cat(s["HDI%"],"percent credible interval:\n")
  cat(" ", s[ c("HDIlo", "HDIup")], "\n")
  cat("The relative frequency of success is more than", s["comp"] , "by a probability of", s["%>comp"], "\n")
  cat("and less than", s["comp"] , "by a probability of", s["%<comp"], "\n")
  cat("\n")
  
  # The output of binom.test
  
  #   Exact binomial test
  # 
  # data:  5 and 10
  # number of successes = 5, number of trials = 10, p-value = 1
  # alternative hypothesis: true probability of success is not equal to 0.5
  # 95 percent confidence interval:
  #  0.187 0.813
  # sample estimates:
  # probability of success 
  #                    0.5 
}

summary.bfa_binom_test <- function(x) {
  s <- round(x$stats, 3)
  
  cat("  Data\n")
  cat("number of successes = ", x$x,", number of trials = ", x$n, "\n", sep="")
  cat("\n")
  
  cat("  Model parameters and generated quantities\n")
  cat("theta: The relative frequency of success\n")
  cat("x_pred: Predicted number of successes in a replication\n")
  cat("\n")
  cat("  Measures\n" )
  print(s[, c("mean", "sd", "HDIlo", "HDIup", "%>comp", "%<comp")])
  cat("\n")
  cat("'HDIlo' and 'HDIup' are the limits of a ", s[1, "HDI%"] ,"% HDI credible interval.\n", sep="")
  cat("'%>comp' and '%<comp' are the probability of the respective parameter being larger\n")
  cat("than ", s[1, "comp"] ,". This comparison might not make sense for all parameters.\n", sep="")
  
  cat("\n")
  cat("  Quantiles\n" )
  print(s[, c("q2.5%", "q25%", "median","q75%", "q97.5%")] )
  # Lägga till HDI
  # 'HDIlo' and 'HDIup' are the limits of a 95% HDI credible interval.
}

plot.bfa_binom_test <- function(x) {
  old_par <- par(mfcol=c(2,1), mar=c(4.2,4,3,2))
  sample_mat <- as.matrix(x$mcmc_samples)
  plotPost(sample_mat[, "theta"], cred_mass= x$conf.level, comp_val=x$p, xlim=c(0, 1), cex=1, cex.lab=1,
           main = "Posterior distribution", xlab="The relative frequency of success")
  hist_data <- discrete_hist(sample_mat[, "x_pred"], c(0, x$n), ylab="Probability", x_marked= x$x,
                             xlab = "Number of sucesses",main="Posterior predictive distribution")
  legend("topright", legend="Data", col="red",  lty=1, lwd=3)
  par(old_par)
}

diagnostics.bfa_binom_test <- function(x) {
  s <- round(x$stats, 3)
  mcmc_samples <- x$mcmc_samples
  
  cat("\n", "Iterations = ", start(mcmc_samples), ":", end(mcmc_samples), "\n", sep = "")
  cat("Thinning interval =", thin(mcmc_samples), "\n")
  cat("Number of chains =", nchain(mcmc_samples), "\n")
  cat("Sample size per chain =", (end(mcmc_samples) - start(mcmc_samples))/thin(mcmc_samples) + 1, "\n")
  
  cat("\n")
  cat("  Diagnostic measures\n")
  print(s[, c("mean", "sd", "mcmc.se", "Rhat", "n.eff")])
  cat("\n")
  cat("MCMC SE: the estimated standard error of the MCMC approximation of the mean.\n")
  cat("n.eff: a crude measure of effective MCMC sample size.\n")
  cat("Rhat: the potential scale reduction factor (at convergence, Rhat=1).\n")
  
  cat("\n")
  cat("  Model parameters and generated quantities\n")
  cat("theta: The relative frequency of success\n")
  cat("x_pred: Predicted number of successes in a replication\n")
  
  plot(x$mcmc_samples)
}

model_diagram.bfa_binom_test <- function(x) {
  print(jags_binom_test)
}

model_code.bfa_binom_test <- function(x) {
  print(jags_binom_test)
}

