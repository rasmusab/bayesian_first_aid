#' Generate MCMC samples from the Bayesian First Aid one sample poisson model using JAGS.
#' 
#' Descritions description description
#' 
#' Details details details
#'
#'  @param x the number of events
#'  @param t The number of time steps during the \code{x} events were observed.
#'  
#'  #' @return
#' An object of type \code{mcmc.list} (defined by the \code{coda} package) that contains the MCMC samples from the model.
jags_one_sample_poisson_test <- function(x, t) {
  model_string <- "
    model {
      x ~ dpois(rate * t)
      rate ~ dgamma(0.5, 0.00001)
    }
  "
  
  jags_model <- jags.model(textConnection(model_string), data=list(x = x, t = t), n.adapt= 200, 
                          inits = list(rate = x / t))
  mcmc_samples <- coda.samples(jags_model, c("rate"), n.iter=10000)
  mcmc_samples
}


jags_two_sample_poisson_test <- function(x1, t1, x2, t2) {
  model_string <- "
    model {
      x1 ~ dpois(rate1 * t1)
      rate1 ~ dgamma(0.5, 0.00001)
      x2 ~ dpois(rate2 * t2)
      rate2 ~ dgamma(0.5, 0.00001)
      rate_diff <- rate2 - rate1
      rate_ratio <- rate1 / rate2
    }
  "
  data_list = list(x1 = x1, t1 = t1, x2 = x2, t2 = t2)
  init_list = list(rate1 = x1 / t1, rate2 = x2 / t2)
  jags_model <- jags.model(textConnection(model_string), data = data_list, n.adapt= 200, 
                           inits = init_list)
  mcmc_samples <- coda.samples(jags_model, c("rate1", "rate2", "rate_diff", "rate_ratio"), n.iter=10000)
  mcmc_samples
}