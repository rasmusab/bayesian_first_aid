jags_cor_test <- function(x1, x2) {
  model_string <- "
    model {
      for(i in 1:n) {
        #x[i,1:2] ~ dmnorm(mu[], prec[ , ])
        x[i,1:2] ~ dmt(mu[], prec[ , ], nu) 
      }

      ## priors for elements of precision matrix
      prec[1:2,1:2] <- inverse(cov[,])
      
      cov[1,1] <- sigma[1] * sigma[1]
      cov[1,2] <- sigma[1] * sigma[2] * rho
      cov[2,1] <- sigma[1] * sigma[2] * rho
      cov[2,2] <- sigma[2] * sigma[2]
      
      sigma[1] ~ dunif(0, 1000) 
      sigma[2] ~ dunif(0, 1000)
      tau[1] <- 1 / pow(sigma[1], 2)
      tau[2] <- 1 / pow(sigma[2], 2)
      rho ~ dunif(-1, 1)

      mu[1] ~ dnorm(0, 0.0001)
      mu[2] ~ dnorm(0, 0.0001)

      nu <- nuMinusOne+1
      nuMinusOne ~ dexp(1/29)
    }
  "
  
  jags_model <- jags.model(textConnection(model_string), data=list(x = cbind(x1, x2), n = length(x1)), n.adapt= 1000, 
                           inits = list(mu=c(median(x1), median(x2)), rho=cor(x1, x2, method="spearman"), sigma = c(mad(x1), mad(x2))))
  update(jags_model, 1000)
  mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "nu"), n.iter=1000)
  mcmc_samples
}