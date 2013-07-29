# ...

jags_one_sample_t_test <- function(x) {
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
  inits_list <- list(mu= mean(x), sigma = sd(x), nuMinusOne = 4)
  
  jagsModel = jags.model(textConnection(model_string) , data=data_list , inits=inits_list , 
                         n.chains=3 , n.adapt=500)
  # Burn-in
  update( jagsModel , 1000, progress.bar="none")
  # running the model.
  # Later increase the number of n.iter steps to 33333
  codaSamples = coda.samples( jagsModel , variable.names=c("mu", "sigma", "nu"),
                              n.iter=333 , thin=3)
  codaSamples
}

#' Adapted from John Kruschke's original BEST code.
jags_two_sample_t_test <- function(x1, x2) {
  model_string <- "
    model {
      for(i in 1:length(x1)) {
        x1[i] ~ dt( mu[1] , tau[1] , nu )
      }
      for(i in 1:length(x2)) {
        x2[i] ~ dt( mu[2] , tau[2] , nu )
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
    x1 = x1 ,
    x2 = x2 ,
    mean_mu = mean(c(x1, x2)) ,
    precision_mu = 1 / (sd(c(x1, x2))^2 * 1000000),
    sigmaLow = sd(c(x1, x2)) / 1000 ,
    sigmaHigh = sd(c(x1, x2)) * 1000 
  )
  
  # Initial values of MCMC chains based on data:
  inits_list <- list(
    mu = c(mean(x1), mean(x2)),
    sigma = c(sd(x1), sd(x2)),
    nuMinusOne = 4 
  )
  
  jagsModel = jags.model(textConnection(model_string) , data=data_list , inits=inits_list , 
                          n.chains=3 , n.adapt=500)
  # Burn-in
  update( jagsModel , 1000, progress.bar="none")
  # running the model.
  # Later increase the number of n.iter steps to 33333
  codaSamples = coda.samples( jagsModel , variable.names=c("mu", "sigma", "nu"),
                              n.iter=333 , thin=3)
  codaSamples
}

jags_paired_t_test <- function(x, y) {
  x_diff <- y - x
  jags_one_sample_t_test(y - x)  
}


# Stand in function untill I include the code from the original t.test function
bfa_t_test.default <- function(x, y = NULL, paired = FALSE) {
  if(is.null(y)) {
    mcmc_samples <- jags_one_sample_t_test(x)
    bfa_object <- list(mcmc_samples = mcmc_samples, x = x)
    class(bfa_object) <- "bfa_one_sample_t_test"
    return()
  } else if(paired) {
    mcmc_samples <- jags_paired_t_test(x, y)
    bfa_object <- list(mcmc_samples = mcmc_samples, x1 = x, x2 = y, 
                       x_diff = y - x)
    class(bfa_object) <- "bfa_paired_t_test"
  } else { # is two sample t.test
    mcmc_samples <- jags_two_sample_t_test(x, y)
    bfa_object <- list(mcmc_samples = mcmc_samples, x1 = x, x2 = y)
    class(bfa_object) <- "bfa_two_sample_t_test"
  }
  
  return(bfa_object)
}

bfa_t_test.formula <- function() {
  
}

### One sample t-test S3 methods ###

print.bfa_one_sample_t_test <- function(x) {
  cat("\n --- Bayesian first aid binomial test ---\n\n")
  print(summary(x$mcmc_samples))
}

summary.bfa_one_sample_t_test <- function(object) {
  cat("\nSummary\n")
  print(object)
}

plot.bfa_one_sample_t_test <- function(x) {
  plot(x$mcmc_samples)
}

model_diagram.bfa_one_sample_t_test <- function(x) {
  print(jags_binom_test)
}

model_code.bfa_one_sample_t_test <- function(x) {
  print(jags_binom_test)
}

### Two sample t-test S3 methods ###

print.bfa_two_sample_t_test <- function(x) {
  cat("\n --- Bayesian first aid two sample t-test ---\n\n")
  
  # Define matrix for storing summary info:
  summaryInfo = matrix( 0 , nrow=9 , ncol=6 , dimnames=list(
    PARAMETER=c( "mu1" , "mu2" , "muDiff" , "sigma1" , "sigma2" , "sigmaDiff" ,
                 "nu" , "nuLog10" , "effSz" ),
    SUMMARY.INFO=c( "mean" , "median" , "mode" , "HDIlow" , "HDIhigh" ,
                    "pcgtZero" ) 
  ))
  
  mcmcChain <- as.matrix(x$mcmc_samples)
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
  
  N1 = length(x$x1)
  N2 = length(x$x2)
  effSzChain = ( ( mcmcChain[,"mu[1]"] - mcmcChain[,"mu[2]"] ) 
                 / sqrt( ( mcmcChain[,"sigma[1]"]^2 + mcmcChain[,"sigma[2]"]^2 ) / 2 ) ) 
  summaryInfo[ "effSz" , ] = mcmcSummary( effSzChain , compVal=0 )
  # Or, use sample-size weighted version:
  # effSz = ( mu1 - mu2 ) / sqrt( ( sigma1^2 *(N1-1) + sigma2^2 *(N2-1) ) 
  #                               / (N1+N2-2) )
  # Be sure also to change plot label in BESTplot function, below.
  return( summaryInfo )
}

summary.bfa_two_sample_t_test <- function(object) {
  print(object)
}

plot.bfa_two_sample_t_test <- function(x) {
  plot(x$mcmc_samples)
}

model_diagram.bfa_two_sample_t_test <- function(x) {
  print(jags_two_sample_t_test)
}

model_code.bfa_two_sample_t_test <- function(x) {
  print(jags_two_sample_t_test)
}

### Paired samples t-test S3 methods ###

print.bfa_paired_t_test <- function(x) {
  cat("\n --- Bayesian first aid binomial test ---\n\n")
  print(summary(x$mcmc_samples))
}

summary.bfa_paired_t_test <- function(object) {
  cat("\nSummary\n")
  print(object)
}

plot.bfa_paired_t_test <- function(x) {
  plot(x$mcmc_samples)
}

model_diagram.bfa_paired_t_test <- function(x) {
  print(jags_binom_test)
}

model_code.bfa_paired_t_test <- function(x) {
  print(jags_binom_test)
}