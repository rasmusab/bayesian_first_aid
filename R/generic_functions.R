
#' Plots and prints diagnostics regarding the convergence of the model.
diagnostics <- function(fit) {
  UseMethod("diagnostics")
}

#' Plots the model diagram in the style of Kruschke(2010)
#' 
#' @references
#' Kruschke, J. (2010) Doing Bayesian Data Analysis: A Tutorial Introduction with R. Academic Press.
model_diagram <- function(fit) {
  UseMethod("model_diagram")
}

#' Prints code that replicates the model you just ran. 
#' 
#' This is good if you better want to understand how the model is 
#' implemented or if you want to run a modified version of the code.
model_code <- function(fit) {
  UseMethod("model_code")
}

bfa_t_test <- function(x, ...) {
  UseMethod("bfa_t_test")
}