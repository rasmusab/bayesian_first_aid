# ...

jags_one_sample_t_test <- function(x, n.adapt= 1000, n.chains=3, n.update = 1000, n.iter=1000, thin=1) {
  model_string <- "
    model {
      for(i in 1:length(x)) {
        x[i] ~ dt( mu , tau , nu )
      }
  
    mu ~ dnorm( mean_mu , precision_mu )
    tau <- 1/pow( sigma , 2 )
    sigma ~ dunif( sigmaLow , sigmaHigh )
    nu <- nuMinusOne+1
    nuMinusOne ~ dexp(1/29)
  }"

  # Specify the data in a list, for later shipment to JAGS:
  data_list <- list(
    x = x,
    mean_mu = mean(x) ,
    precision_mu = 1 / (sd(x)^2 * 1000000),
    sigmaLow = sd(x) / 1000 ,
    sigmaHigh = sd(x) * 1000 
  )
  
  # Initial values of MCMC chains based on data:
  inits_list <- list(mu= mean(x, trim=0.2), sigma = mad(x), nuMinusOne = 4)
  
  jagsModel = jags.model(textConnection(model_string) , data=data_list , inits=inits_list , 
                         n.chains=n.chains , n.adapt=n.adapt)
  # Burn-in
  update( jagsModel , n.update, progress.bar="none")
  # running the model.
  # Later increase the number of n.iter steps to 33333
  codaSamples = coda.samples( jagsModel , variable.names=c("mu", "sigma", "nu"),
                              n.iter=n.iter, thin=thin)
  codaSamples
}

#' Adapted from John Kruschke's original BEST code.
jags_two_sample_t_test <- function(x, y, n.adapt= 1000, n.chains=3, n.update = 1000, n.iter=1000, thin=1) {
  model_string <- "
    model {
      for(i in 1:length(x)) {
        x[i] ~ dt( mu[1] , tau[1] , nu )
      }
      for(i in 1:length(y)) {
        y[i] ~ dt( mu[2] , tau[2] , nu )
      }

      mu[1] ~ dnorm( mean_mu , precision_mu )
      tau[1] <- 1/pow( sigma[1] , 2 )
      sigma[1] ~ dunif( sigmaLow , sigmaHigh )
      mu[2] ~ dnorm( mean_mu , precision_mu )
      tau[2] <- 1/pow( sigma[2] , 2 )
      sigma[2] ~ dunif( sigmaLow , sigmaHigh )
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
    mu = c(mean(x), mean(y)),
    sigma = c(sd(x), sd(y)),
    nuMinusOne = 4 
  )
  
  jagsModel = jags.model(textConnection(model_string) , data=data_list , inits=inits_list , 
                          n.chains=n.chains , n.adapt=n.adapt)
  # Burn-in
  update( jagsModel , n.update, progress.bar="none")
  # running the model.
  # Later increase the number of n.iter steps to 33333
  codaSamples = coda.samples( jagsModel , variable.names=c("mu", "sigma", "nu"),
                              n.iter=n.iter , thin=thin)
  codaSamples
}

jags_paired_t_test <- function(x, y, n.adapt= 1000, n.chains=3, n.update = 1000, n.iter=1000, thin=1) {
  x_diff <- y - x
  jags_one_sample_t_test(y - x, n.adapt, n.chains, n.update, n.iter, thin)  
}


# Stand in function untill I include the code from the original t.test function
bfa.t.test.default <- function(x, y = NULL, paired = FALSE) {
  if(is.null(y)) {
    mcmc_samples <- jags_one_sample_t_test(x)
    bfa_object <- list(mcmc_samples = mcmc_samples, x = x)
    class(bfa_object) <- "bfa_one_sample_t_test"
    return()
  } else if(paired) {
    mcmc_samples <- jags_paired_t_test(x, y)
    bfa_object <- list(mcmc_samples = mcmc_samples, x = x, y = y, 
                       x_diff = y - x)
    class(bfa_object) <- "bfa_paired_t_test"
  } else { # is two sample t.test
    mcmc_samples <- jags_two_sample_t_test(x, y)
    bfa_object <- list(mcmc_samples = mcmc_samples, x = x, y = y)
    class(bfa_object) <- "bfa_two_sample_t_test"
  }
  
  return(bfa_object)
}

bfa.t.test.formula <- function() {
  
}

### One sample t-test S3 methods ###

print.bfa_one_sample_t_test <- function(bfa_result) {
  cat("\n --- Bayesian first aid binomial test ---\n\n")
  print(summary(bfa_result$mcmc_samples))
}

summary.bfa_one_sample_t_test <- function(bfa_result) {
  cat("\nSummary\n")
  print(bfa_result)
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
  summaryInfo[ "mu1" , ] = mcmcSummary( mcmcChain[,"mu[1]"] )
  summaryInfo[ "mu2" , ] = mcmcSummary( mcmcChain[,"mu[2]"] )
  summaryInfo[ "muDiff" , ] = mcmcSummary( mcmcChain[,"mu[1]"]
                                           - mcmcChain[,"mu[2]"] , 
                                           compVal=0 )
  summaryInfo[ "sigma1" , ] = mcmcSummary( mcmcChain[,"sigma[1]"] )
  summaryInfo[ "sigma2" , ] = mcmcSummary( mcmcChain[,"sigma[2]"] )
  summaryInfo[ "sigmaDiff" , ] = mcmcSummary( mcmcChain[,"sigma[1]"]
                                              - mcmcChain[,"sigma[2]"] , 
                                              compVal=0 )
  summaryInfo[ "nu" , ] = mcmcSummary( mcmcChain[,"nu"] )
  summaryInfo[ "nuLog10" , ] = mcmcSummary( log10(mcmcChain[,"nu"]) )
  
  N1 = length(x$x)
  N2 = length(x$y)
  effSzChain = ( ( mcmcChain[,"mu[1]"] - mcmcChain[,"mu[2]"] ) 
                 / sqrt( ( mcmcChain[,"sigma[1]"]^2 + mcmcChain[,"sigma[2]"]^2 ) / 2 ) ) 
  summaryInfo[ "effSz" , ] = mcmcSummary( effSzChain , compVal=0 )
  # Or, use sample-size weighted version:
  # effSz = ( mu1 - mu2 ) / sqrt( ( sigma1^2 *(N1-1) + sigma2^2 *(N2-1) ) 
  #                               / (N1+N2-2) )
  # Be sure also to change plot label in BESTplot function, below.
  return( summaryInfo )
}

summary.bfa_two_sample_t_test <- function(bfa_result) {
  print(bfa_result)
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