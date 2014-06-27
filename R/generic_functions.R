#' Plots and prints diagnostics regarding the convergence of the model.
#' 
#' @param fit The output from a Bayesian First Aid model.
#' 
#' @export
diagnostics <- function(fit) {
  UseMethod("diagnostics")
}

#' Prints code that replicates the model you just ran. 
#' 
#' This is good if you better want to understand how the model is 
#' implemented or if you want to run a modified version of the code.
#' 
#' 
#' @param fit The output from a Bayesian First Aid model.
#' 
#' @export
model.code <- function(fit) {
  UseMethod("model.code")
}

#' Returns the MCMC samples from a Bayesian First Aid fit as a data frame.
#' 
#' Returns a dataframe with the MCMC samples from a Bayesian First Aid Object. 
#' The samples does not include the burn in phase and all chains are collapsed 
#' in the resulting data frame.
#' 
#' In the case where there are many parameters with the same name but different 
#' index (for example, mu[1], mu[2], mu[3], etc.) these names will be rewritten 
#' by removing the brackets (that is, mu[1] -> mu1, mu[2] -> mu2, sigma[1,2] ->
#' sigma12, etc.).
#' 
#' @param x The output from a Bayesian First Aid model.
#' @param ... Not used.
#'   
#' @return A data frame with one column per parameter.
#' 
#' @examples
#' 
#' fit <- bayes.t.test(rnorm(10), rnorm(10))
#' d <- as.data.frame(fit)
#' mean(d$mu_x)
#' @export
as.data.frame.bayesian_first_aid <- function(x, ...) {
  d <- as.data.frame(as.matrix(x$mcmc_samples))
  colnames(d) <- str_replace_all(colnames(d), "\\[|\\]", "")
  d
}

#' Returns the MCMC samples from a Bayesian First Aid fit as a matrix.
#' 
#' Returns a matrix with the MCMC samples from a Bayesian First Aid Object. The
#' samples does not include the burn in phase and all chains are collapsed in
#' the resulting matrix.
#' 
#' This is just a wrapper for the \code{as.matrix} from the coda package.
#' Therefor, paramameters with the same name but different index will be column
#' names that includes brackets (for example, "mu[1]", "mu[2]", "sigma[1,2],
#' etc.)
#' 
#' @param x The output from a Bayesian First Aid model.
#' @param ... Not used
#'   
#' @return A matrix with one column per parameter.
#' 
#' @examples
#' 
#' fit <- bayes.t.test(rnorm(10), rnorm(10))
#' d <- as.matrix(fit)
#' mean(d[,"mu_x"])
#' @export
as.matrix.bayesian_first_aid <- function(x, ...) {
  as.matrix(x$mcmc_samples)
}
