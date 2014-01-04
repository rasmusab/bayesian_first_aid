

#' Plots and prints diagnostics regarding the convergence of the model.
#' 
#' Descritions description description
#' 
#' Details details details
#' 
#' @param x What is this param
#' 
#' @return
#' An object of type something...
#' @export
diagnostics <- function(fit) {
  UseMethod("diagnostics")
}



#' 

#' Prints code that replicates the model you just ran. 
#' 
#' This is good if you better want to understand how the model is 
#' implemented or if you want to run a modified version of the code.
#' 
#' Details details details
#' 
#' @param x What is this param
#' 
#' @return
#' An object of type something...
#' @export
model.code <- function(fit) {
  UseMethod("model.code")
}