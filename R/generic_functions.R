

#' Plots and prints diagnostics regarding the convergence of the model.
#' 
#' Descritions description description
#' 
#' Details details details
#' 
#' @param fit
#' 
#' @return
#' An object of type something...
#' @export
diagnostics <- function(fit) {
  UseMethod("diagnostics")
}

#' Prints code that replicates the model you just ran. 
#' 
#' This is good if you better want to understand how the model is 
#' implemented or if you want to run a modified version of the code.
#' 
#' Details details details
#' 
#' @param fit
#' 
#' @return
#' An object of type something...
#' @export
model.code <- function(fit) {
  UseMethod("model.code")
}

#' Returns the MCMC samples from a Bayesian First Aid fit as a data frame.
#' 
#' Returns a dataframe with the MCMC samples from a Bayesian First Aid Object. The samples does not include the burn in phase and all chains are collapsed in the resulting data frame.
#' 
#' Details details details
#' 
#' @param fit
#' 
#' @return
#' A data frame
#' @export
as.data.frame.bayesian_first_aid <- function(x, ...) {
  as.data.frame(as.matrix(fit$mcmc_samples))
}

#' Returns the MCMC samples from a Bayesian First Aid fit as a matrix.
#' 
#' Returns a matrix with the MCMC samples from a Bayesian First Aid Object. The samples does not include the burn in phase and all chains are collapsed in the resulting matrix.
#' 
#' Details details details
#' 
#' @param fit
#' 
#' @return
#' A matrix
#' @export
as.matrix.bayesian_first_aid <- function(x, ...) {
  as.matrix(fit$mcmc_samples)
}
