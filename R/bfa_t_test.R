# ...


jags_one_sample_t_test <- function(x, comp_mu = 0,n.adapt= 100, n.chains=3, n.update = 100, n.iter=500, thin=1) {
  model_string <- "
    model {
      for(i in 1:length(x)) {
        x[i] ~ dt( mu , tau , nu )
      }
      x_pred ~ dt( mu , tau , nu )
      eff_size <- (mu - comp_mu) / sigma
  
      mu ~ dnorm( mean_mu , precision_mu )
      tau <- 1/pow( sigma , 2 )
      sigma ~ dunif( sigmaLow , sigmaHigh )
      nu <- nuMinusOne+1
      nuMinusOne ~ dexp(1/29)
    }"

  data_list <- list(
    x = x,
    mean_mu = mean(x) ,
    precision_mu = 1 / (sd(x)^2 * 1000000),
    sigmaLow = sd(x) / 1000 ,
    sigmaHigh = sd(x) * 1000 ,
    comp_mu = comp_mu
  )
  
  
  inits_list <- list(mu= mean(x, trim=0.2), sigma = mad(x), nuMinusOne = 4)
  params <- c("mu", "sigma", "nu", "eff_size", "x_pred")
  mcmc_samples <- run_jags(model_string, data = data_list, inits = inits_list, 
                           params = params, n.chains = n.chains, n.adapt = n.adapt,
                           n.update = n.update, n.iter = n.iter, thin = thin)   
  mcmc_samples
}

#' Adapted from John Kruschke's original BEST code.
jags_two_sample_t_test <- function(x, y, n.adapt= 100, n.chains=3, n.update = 100, n.iter=500, thin=1) {
  model_string <- "
    model {
      for(i in 1:length(x)) {
        x[i] ~ dt( mu_x , tau_x , nu )
      }
      x_pred ~ dt( mu_x , tau_x , nu )
      for(i in 1:length(y)) {
        y[i] ~ dt( mu_y , tau_y , nu )
      }
      y_pred ~ dt( mu_y , tau_y , nu )
      eff_size <- (mu_x - mu_y) / sqrt((pow(sigma_x, 2) + pow(sigma_y, 2)) / 2)
      mu_diff <- mu_x - mu_y
      sigma_diff <-sigma_x - sigma_y 
      
      mu_x ~ dnorm( mean_mu , precision_mu )
      tau_x <- 1/pow( sigma_x , 2 )
      sigma_x ~ dunif( sigmaLow , sigmaHigh )

      mu_y ~ dnorm( mean_mu , precision_mu )
      tau_y <- 1/pow( sigma_y , 2 )
      sigma_y ~ dunif( sigmaLow , sigmaHigh )

      nu <- nuMinusOne+1
      nuMinusOne ~ dexp(1/29)
    }"

  # Specify the data in a list, for later shipment to JAGS:
  data_list <- list(
    x = x ,
    y = y ,
    mean_mu = mean(c(x, y)) ,
    precision_mu = 1 / (sd(c(x, y))^2 * 1000000),
    sigmaLow = sd(c(x, y)) / 1000 ,
    sigmaHigh = sd(c(x, y)) * 1000 
  )
  
  # Initial values of MCMC chains based on data:
  inits_list <- list(
    mu_x = mean(x, trim=0.2),
    mu_y = mean(y, trim=0.2),
    sigma_x = mad(x),
    sigma_y = mad(y),
    nuMinusOne = 4 
  )
  
  params <- c("mu_x", "sigma_x", "mu_y", "sigma_y", "mu_diff", "sigma_diff","nu", "eff_size", "x_pred", "y_pred")  
  mcmc_samples <- run_jags(model_string, data = data_list, inits = inits_list, 
                           params = params, n.chains = n.chains, n.adapt = n.adapt,
                           n.update = n.update, n.iter = n.iter, thin = thin)
  mcmc_samples
}

# Right now, this is basically just calling jags_one_sample_t_test but I'm 
# keeping it in case I would want to change it in the future.
jags_paired_t_test <- function(x, y, comp_mu = 0,n.adapt= 100, n.chains=3, n.update = 100, n.iter=500, thin=1) {
  if(is.null(y)) { # assume x is the aldread calculated difference between the two groups
    x_diff <- x
  } else {
    x_diff <- y - x
  }
  jags_one_sample_t_test(x_diff, comp_mu=comp_mu,n.adapt=n.adapt, n.chains=n.chains, n.update=n.update, 
                         n.iter = n.iter, thin=thin)  
}


bfa.t.test.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
                               mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95, ...) {
  
  if(var.equal) {
    var.equal <- FALSE
    warning("To assume equal variance of 'x' and 'y' is not supported. Continuing by estimating the variance of 'x' and 'y' separately.")
  }
  
  ### Original (but slighly modified) code from t.test.default ###
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                                 conf.level < 0 || conf.level > 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  
  # removing incomplete cases and preparing the data vectors (x & y)
  if (!is.null(y)) {
    xname <- deparse(substitute(x))
    yname <- deparse(substitute(y))
    dname <- paste(xname, "and", yname)
    if (paired) 
      xok <- yok <- complete.cases(x, y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else {
    xname <- deparse(substitute(x))
    dname <- xname
    if (paired) 
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  
  # Checking that there is enough data. Even though BEST handles the case with
  # one data point it is still usefull to do these checks.
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  if (is.null(y)) {
    if (nx < 2) 
      stop("not enough 'x' observations")
    df <- nx - 1
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx)) 
      stop("data are essentially constant")
  }
  else {
    ny <- length(y)
    if (nx < 2) 
      stop("not enough 'x' observations")
    if (ny < 2) 
      stop("not enough 'y' observations")
    my <- mean(y)
    vy <- var(y)
    stderrx <- sqrt(vx/nx)
    stderry <- sqrt(vy/ny)
    stderr <- sqrt(stderrx^2 + stderry^2)
    df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), 
                                                abs(my))) 
      stop("data are essentially constant")
  }
  
  ### Own code starts here ###
  
  if(paired) {
    mcmc_samples <- jags_paired_t_test(x, y)
    stats <- mcmc_stats(mcmc_samples, cred_mass = conf.level, comp_val = mu)
    bfa_object <- list(x = x, y = y, x_diff = y - x, comp = mu,
                        x_name = xname, y_name = yname, data_name = dname,
                        mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- "bfa_paired_t_test"
  
  } else if(is.null(y)) {
    mcmc_samples <- jags_one_sample_t_test(x, comp_mu = mu)
    stats <- mcmc_stats(mcmc_samples, cred_mass = conf.level, comp_val = mu)
    bfa_object <- list(x = x, comp = mu, x_name = xname, data_name = dname,
                       mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- "bfa_one_sample_t_test"
    
  } else { # is two sample t.test
    mcmc_samples <- jags_two_sample_t_test(x, y)
    stats <- mcmc_stats(mcmc_samples, cred_mass = conf.level, comp_val = mu)
    bfa_object <- list(x = x, y = y, comp = mu,
                       x_name = xname, y_name = yname, data_name = dname,
                       mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- "bfa_two_sample_t_test"
  }
  bfa_object
}

bfa.t.test.formula <- function(formula, data, subset, na.action, ...) {
  
  ### Original code from t.test.formula ###
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), "term.labels")) != 1L)) 
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L) 
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  
  ### Own code starts here ###
  bfa_object <- do.call("bfa.t.test", c(DATA, list(...)))
  bfa_object$data_name <- DNAME
  if (length(levels(g)) == 2L) {
    bfa_object$x_name <- paste("group", levels(g)[1])
    bfa_object$y_name <- paste("group", levels(g)[2])
  }
  bfa_object
  
}

### One sample t-test S3 methods ###

print.bfa_one_sample_t_test <- function(x) {
  
  s <- round(x$stats, 3)
  
  cat("\n")
  cat("\tBayesian estimation supersedes the t test (BEST)\n")
  cat("\n")
  cat("data: ", x$data_name, "\n", sep="")
  cat("\n")
  cat("  Estimates [", s[1, "HDI%"] ," percent credible interval]\n", sep="")
  cat("mean of ",  x$x_name, ": ", s["mu", "mean"], " [", s["mu", "HDIlo"], ", ", s["mu", "HDIup"] , "]\n",sep="")
  cat("sd of ",  x$x_name, ": ", s["sigma", "mean"], " [", s["sigma", "HDIlo"], ", ", s["sigma", "HDIup"] , "]\n",sep="")
  
  cat("\n")
  cat("The mean is more than", s["mu","comp"] , "by a probability of", s["mu","%>comp"], "\n")
  cat("and less than", s["mu", "comp"] , "by a probability of", s["mu", "%<comp"], "\n")
  cat("\n")
  
  # Output from t.test , one sample
  #
  #  One Sample t-test
  #
  #data:  my_x
  #t = 0.223, df = 9, p-value = 0.8284
  #alternative hypothesis: true mean is not equal to 0
  #95 percent confidence interval:
  # -0.365  0.445
  #sample estimates:
  #mean of x 
  #   0.0399 
}

summary.bfa_one_sample_t_test <- function(x) {
  s <- round(x$stats, 3)
  print(s)
#   cat("  Data\n")
#   cat(x$data_name, ", n = ", length(x$x), sep="")
#   cat("\n")
#   
#   cat("  Model parameters and generated quantities\n")
#   cat("mu: The relative frequency of success\n")
#   cat("sigma: The relative frequency of success\n")
#   cat("nu: The relative frequency of success\n")
#   cat("x_pred: Predicted number of successes in a replication\n")
#   cat("\n")
#   cat("  Measures\n" )
#   print(stats[, c("mean", "sd", "HDIlo", "HDIup", "%>comp", "%<comp")])
#   cat("\n")
#   cat("'HDIlo' and 'HDIup' are the limits of a ", stats[1, "HDI%"] ,"% HDI credible interval.\n", sep="")
#   cat("'%>comp' and '%<comp' are the probability of the respective parameter being larger\n")
#   cat("than ", stats[1, "comp"] ,". This comparison might not make sense for all parameters.\n", sep="")
#   
#   cat("\n")
#   cat("  Quantiles\n" )
#   print(stats[, c("q2.5%", "q25%", "median","q75%", "q97.5%")] )

}

plot.bfa_one_sample_t_test <- function(bfa_result) {
  plot(bfa_result$mcmc_samples)
}

diagnostics.bfa_one_sample_t_test <- function(bfa_result) {
  plot(bfa_result$mcmc_samples)
}

model_diagram.bfa_one_sample_t_test <- function(bfa_result) {
  print(jags_binom_test)
}

model_code.bfa_one_sample_t_test <- function(bfa_result) {
  print(jags_binom_test)
}

### Two sample t-test S3 methods ###

print.bfa_two_sample_t_test <- function(bfa_result) {
  cat("\n --- Bayesian first aid two sample t-test ---\n\n")
  
  # Define matrix for storing summary info:
  summaryInfo = matrix( 0 , nrow=9 , ncol=6 , dimnames=list(
    PARAMETER=c( "mu1" , "mu2" , "muDiff" , "sigma1" , "sigma2" , "sigmaDiff" ,
                 "nu" , "nuLog10" , "effSz" ),
    SUMMARY.INFO=c( "mean" , "median" , "mode" , "HDIlow" , "HDIhigh" ,
                    "pcgtZero" ) 
  ))
  
  mcmcChain <- as.matrix(bfa_result$mcmc_samples)
  summaryInfo[ "mu1" , ] = mcmcSummary( mcmcChain[,"mu_x"] )
  summaryInfo[ "mu2" , ] = mcmcSummary( mcmcChain[,"mu_y"] )
  summaryInfo[ "muDiff" , ] = mcmcSummary( mcmcChain[,"mu_x"]
                                           - mcmcChain[,"mu_y"] , 
                                           compVal=0 )
  summaryInfo[ "sigma1" , ] = mcmcSummary( mcmcChain[,"sigma_x"] )
  summaryInfo[ "sigma2" , ] = mcmcSummary( mcmcChain[,"sigma_y"] )
  summaryInfo[ "sigmaDiff" , ] = mcmcSummary( mcmcChain[,"sigma_x"]
                                              - mcmcChain[,"sigma_y"] , 
                                              compVal=0 )
  summaryInfo[ "nu" , ] = mcmcSummary( mcmcChain[,"nu"] )
  summaryInfo[ "nuLog10" , ] = mcmcSummary( log10(mcmcChain[,"nu"]) )
  
  N1 = length(bfa_result$x)
  N2 = length(bfa_result$y)
  effSzChain = ( ( mcmcChain[,"mu_x"] - mcmcChain[,"mu_y"] ) 
                 / sqrt( ( mcmcChain[,"sigma_x"]^2 + mcmcChain[,"sigma_y"]^2 ) / 2 ) ) 
  summaryInfo[ "effSz" , ] = mcmcSummary( effSzChain , compVal=0 )
  # Or, use sample-size weighted version:
  # effSz = ( mu1 - mu2 ) / sqrt( ( sigma1^2 *(N1-1) + sigma2^2 *(N2-1) ) 
  #                               / (N1+N2-2) )
  # Be sure also to change plot label in BESTplot function, below.
  print(summaryInfo)
}

summary.bfa_two_sample_t_test <- function(bfa_result) {
  s <- round(x$stats, 3)
  print(s)
}

plot.bfa_two_sample_t_test <- function(bfa_result) {
  plot(bfa_result$mcmc_samples)
}

diagnostics.bfa_two_sample_t_test <- function(bfa_result) {
  plot(bfa_result$mcmc_samples)
}

model_diagram.bfa_two_sample_t_test <- function(bfa_result) {
  print(jags_two_sample_t_test)
}

model_code.bfa_two_sample_t_test <- function(bfa_result) {
  print(jags_two_sample_t_test)
}

### Paired samples t-test S3 methods ###

print.bfa_paired_t_test <- function(bfa_result) {
  cat("\n --- Bayesian first aid binomial test ---\n\n")
  print(summary(bfa_result$mcmc_samples))
}

summary.bfa_paired_t_test <- function(bfa_result) {
  cat("\nSummary\n")
  print(bfa_result)
}

plot.bfa_paired_t_test <- function(bfa_result) {
  plot(bfa_result$mcmc_samples)
}

model_diagram.bfa_paired_t_test <- function(bfa_result) {
  print(jags_binom_test)
}

diagnostics.bfa_paired_t_test <- function(bfa_result) {
  plot(bfa_result$mcmc_samples)
}

model_code.bfa_paired_t_test <- function(bfa_result) {
  print(jags_binom_test)
}