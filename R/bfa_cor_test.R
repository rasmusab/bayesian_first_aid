jags_cor_test <- function(x1, x2) {
  model_string <- "
    model {
      for(i in 1:n) {
        x[i,1:2] ~ dmnorm(mu[], prec[ , ])
        #x[i,1:2] ~ dmt(mu[], prec[ , ], nu)  # <- not working
      }

      ## priors for elements of precision matrix
      prec[1:2,1:2] ~ dwish(R[,],k)
      R[1,1] <- .01
      R[1,2] <- 0
      R[2,1] <- 0
      R[2,2] <- .01
      k <- 2
      mu[1] ~ dnorm(0, 0.0001)
      mu[2] ~ dnorm(0, 0.0001)
      sigma[1:2,1:2] <- inverse(prec[ , ])
      rho <- sigma[1,2]/sqrt(sigma[1,1]*sigma[2,2])
      nu ~ dunif(2, 250)
    }
  "
  
  jags_model <- jags.model(textConnection(model_string), data=list(x = cbind(x1, x2), n = length(x1)), n.adapt= 1000)
  update(jags_model, 1000)
  mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "nu"), n.iter=1000)
  mcmc_samples
}