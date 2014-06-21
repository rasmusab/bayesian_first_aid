# Todo
# comp.theta kan inte vara olika för de olika grupperna, men det kan p i prop.test


#'Bayesian First Aid Alternative to a Test of Proportions
#'
#'\code{bayes.prop.test} estimates the relative frequency of success for two or 
#'more groups using Bayesian estimation and is intended as a replacement for 
#'\code{\link{prop.test}}.
#'
#'Given data on the number of successes and failures \code{bayes.prop.test} 
#'estimates \eqn{\theta}₁…ₙ, the relative frequencies of success for each of the
#'\eqn{n} groups. The following model is assumed for each for each group:
#'
#'\deqn{x \sim \mathrm{Binom}(\theta, n)}{x ~ Binomial(\theta, n)} \deqn{\theta
#'\sim \mathrm{Beta}(1, 1)}{\theta ~ Beta(1, 1)}
#'
#'
#'Here the prior on the \eqn{\theta}s is a non-informative \eqn{\mathrm{Beta}(1,
#'1)}{Beta(1, 1)} distribution which is identical to a \eqn{\mathrm{Uniform}(0, 
#'1)}{Uniform(0, 1)} distribution. By \code{plot}ing and looking at a 
#'\code{summary} of the object returned by \code{bayes.prop.test} you can get 
#'more information about the shape of the posterior and the posterior predictive
#'distribution. \code{\link{model.code}} prints out the corresponding R code
#'underlying \code{bayes.prop.test} which can be copy-n-pasted into an R script
#'and modified, for example, changing the prior on \eqn{\theta}.
#'
#'
#'@param x a vector of counts of successes, a one-dimensional table with two
#'  entries, or a two-dimensional table (or matrix) with 2 columns, giving the
#'  counts of successes and failures, respectively.
#'@param n a vector of counts of trials; ignored if x is a matrix or a table.
#'@param comp.theta a vector of fixed relative frequencies of success to compare
#'  with the estimated relative frequency of success. The length of \code{comp.theta} must be
#'  the same as the number of groups specified by \code{x}, and its elements must be
#'  greater than 0 and less than 1. This argument fills a similar role as
#'  \code{p} in \code{\link{prop.test}}.
#'@param alternative ignored and is only retained in order to mantain 
#'  compatibility with \code{\link{prop.test}}.
#'@param cred.mass the amount of probability mass that will be contained in 
#'  reported credible intervals. This argument fills a similar role as 
#'  \code{conf.level} in \code{\link{prop.test}}.
#'@param correct ignored and is only retained in order to mantain 
#'  compatibility with \code{\link{prop.test}}
#'@param n.iter The number of iterations to run the MCMC sampling.
#'@param progress.bar The type of progress bar. Possible values are "text", 
#'  "gui", and "none".
#'@param p same as \code{comp.theta} and is only retained in order to mantain 
#'  compatibility with \code{\link{prop.test}}.
#'@param conf.level same as \code{cred.mass} and is only retained in order to 
#'  mantain compatibility with \code{\link{prop.test}}.
#'  
#'
#'  
#'@return A list of class \code{bayes_prop_test} that contains information 
#'  about the analysis. It can be further inspected using the functions 
#'  \code{summary}, \code{plot}, \code{\link{diagnostics}} and 
#'  \code{\link{model.code}}.
#'  
#' @examples
#' 
#' TODO, fixa!
#' 
#' # A college dormitory recently sponsored a taste comparison between
#' # two major soft drinks. Of the 64 students who participated, 39 selected
#' # brand A, and only 25 selected brand B. Example from 
#' # http://www.elderlab.yorku.ca/~aaron/Stats2022/BinomialTest.htm
#' 
#' bayes.binom.test(x = 39, n = 64)
#' 
#' # Save the return value in order to inspect the model result further.
#' fit <- bayes.binom.test(x = 39, n = 64, cred.mass=0.8)
#' summary(fit)
#' plot(fit)
#' 
#' # MCMC diagnostics (should not be necessary for such a simple model)
#' diagnostics(fit)
#' 
#' # Print out the R code to run the model. This can be copy n' pasted into
#' # an R-script and further modified.
#' model.code(fit)
#' 
#' @seealso \code{\link{bayes.binom.test}} for when you want to estimate the relative frequency for only one group.
#' 
#'@export


bayes.prop.test <- function (x, n, comp.theta = NULL, alternative = NULL, cred.mass = 0.95, correct = NULL, n.iter=15000, progress.bar="none", p, conf.level) {
  
  if(! missing(alternative)) {
    warning("The argument 'alternative' is ignored by bayes.prop.test")
  }
  
  if(! missing(correct)) {
    warning("The argument 'correct' is ignored by bayes.prop.test")
  }
  
  if(! missing(p)) {
    comp.theta <- p
  }
  
  if(! missing(conf.level)) {
    cred.mass <- conf.level
  }
  
  x_name <- deparse(substitute(x))
  n_name <- deparse(substitute(n))
  
  ### Begin slightly modified code from prop.test 
  DNAME <- deparse(substitute(x))
  if (is.table(x) && length(dim(x)) == 1L) {
    if (dim(x) != 2L) 
      stop("table 'x' should have 2 entries")
    l <- 1
    n <- sum(x)
    x <- x[1L]
  }
  else if (is.matrix(x)) {
    if (ncol(x) != 2L) 
      stop("'x' must have 2 columns")
    l <- nrow(x)
    n <- rowSums(x)
    x <- x[, 1L]
  }
  else {
    DNAME <- paste(DNAME, "out of", deparse(substitute(n)))
    if ((l <- length(x)) != length(n)) 
      stop("'x' and 'n' must have the same length")
  }
  OK <- complete.cases(x, n)
  x <- x[OK]
  n <- n[OK]
  if ((k <- length(x)) < 1L) 
    stop("not enough data")
  if (any(n <= 0)) 
    stop("elements of 'n' must be positive")
  if (any(x < 0)) 
    stop("elements of 'x' must be nonnegative")
  if (any(x > n)) 
    stop("elements of 'x' must not be greater than those of 'n'")
  if(length(comp.theta) == 1) {
    comp.theta <- rep(comp.theta, length(x))
  }
  if (is.null(comp.theta) && (k == 1)) 
    comp.theta <- 0.5
  if (!is.null(comp.theta)) {
    if (length(comp.theta) != l) 
      stop("'comp.theta' must have the same length as 'x' and 'n' or be a single number")
    comp.theta <- comp.theta[OK]
    if (any((comp.theta <= 0) | (comp.theta >= 1))) 
      stop("elements of 'comp.theta' must be in (0,1)")
  }
  if ((length(cred.mass) != 1L) || is.na(cred.mass) ||
        (cred.mass <= 0) || (cred.mass >= 1)) 
    stop("'cred.mass' must be a single number between 0 and 1")

  ### END code from prop.test 
  
  if(length(x) == 1) {
    return(bayes.binom.test(x, n, comp.theta, cred.mass = ifelse(is.null(cred.mass), 0.5, cred.mass),
                            n.iter = n.iter, progress.bar = progress.bar))
  }
  
  mcmc_samples <- jags_prop_test(x, n, n.chains = 3, n.iter = ceiling(n.iter / 3) , progress.bar=progress.bar)
  if(is.null(comp.theta)) {
    stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = 0.5)
  } else {
    stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = comp.theta)
  }
 
  bfa_object <- list(x = x, n = n, comp_theta = comp.theta, cred_mass = cred.mass,
                     x_name = x_name, n_name = n_name, data_name = DNAME,
                     mcmc_samples = mcmc_samples, stats = stats) 
  class(bfa_object) <- c("bayes_prop_test", "bayesian_first_aid")
  bfa_object
}

prop_model_string <- "model {
  for(i in 1:length(x)) {
    x[i] ~ dbinom(theta[i], n[i])
    theta[i] ~ dbeta(1, 1)
    x_pred[i] ~ dbinom(theta[i], n[i])
  }
}"

jags_prop_test <- function(x, n, n.chains=3, n.iter=5000, progress.bar="none") {
  mcmc_samples <- run_jags(prop_model_string, data = list(x = x, n = n), inits = list(theta = (x + 1) / (n + 2)), 
                           params = c("theta", "x_pred"), n.chains = n.chains, n.adapt = 0,
                           n.update = 0, n.iter = n.iter, thin = 1, progress.bar=progress.bar)
  mcmc_samples
}


format_group_diffs <- function(bfa_object) {
  samples <- as.matrix(bfa_object$mcmc_samples)
  param_names <- colnames(samples)[ str_detect(colnames(samples), "theta")]
  diff_mat <- matrix("", nrow = (length(param_names) - 1) * 2, ncol = length(param_names) - 1)
  for(i in 1:(length(param_names) - 1)) {
    for(j in 2:length(param_names)) {
      if(! i >= j) { 
        param1 <- paste("theta[", i, "]",sep="")
        param2 <- paste("theta[", j, "]", sep="")
        theta_diff <- samples[,param1] - samples[,param2]
        med_theta <- round(median(theta_diff), 2)
        hdi <- signif(HDIofMCMC(theta_diff, bfa_object$cred_mass), 2)
        diff_mat[1 + (i - 1) * 2, j - 1] <- med_theta
        diff_mat[2 + (i - 1) * 2, j - 1] <- paste("[", hdi[1], ", ", hdi[2], "]", sep="")
      }
    }
  }
  diff_mat <- format(diff_mat, width = max(nchar(diff_mat)), justify = "centre")
  rownames(diff_mat) <- rep(1:ceiling(nrow(diff_mat) / 2), each = 2)
  rownames(diff_mat)[1:nrow(diff_mat) %% 2 == 0] <- ""
  rownames(diff_mat) <- paste0("  ", rownames(diff_mat))
  colnames(diff_mat) <- format(as.character(1:ncol(diff_mat) + 1), width = max(nchar(diff_mat)), justify = "centre")
  diff_mat
}


#' @export
print.bayes_prop_test <- function(x, ...) {
  
  s <- format_stats(x$stats)
  
  cat("\n")
  cat("\tBayesian First Aid propotion test\n")
  cat("\n")
  cat("data: ", x$data_name, "\n", sep="")
  pad_width <- max(nchar(as.character(c(x$x, x$n)))) + 1
  cat("number of successes: ", paste(str_pad(x$x, pad_width), collapse = ","), "\n", sep="")
  cat("number of trials:    ", paste(str_pad(x$n, pad_width), collapse = ","), "\n", sep="")
  cat("Estimated relative frequency of success [", s[1, "HDI%"] ,"% credible interval]:\n", sep="")
  for(param_i in which(str_detect(rownames(s), "theta"))) {
    param <- paste("theta[", param_i, "]", sep="")
    cat("  Group ", param_i,": " ,s[param, "median"], " [", paste(s[param, c("HDIlo", "HDIup")], collapse = ", "),"]\n", sep = "")
  }
  
  group_diffs <- format_group_diffs(x)
  if(ncol(group_diffs) > 1) {
    cat("Estimated pairwise group differences (row - column) with", s[1, "HDI%"] ,"% cred. intervals:\n")
    cat(format("Group", width = 2 + nchar(rownames(group_diffs)[1]) * 2 + sum(nchar(colnames(group_diffs))),
               justify = "centre"), "\n", sep="") 
    print(format_group_diffs(x), quote=FALSE)
  } else {
    cat("Estimated group difference (Group 1 - Group 2):\n")  
    cat("  ", str_trim(group_diffs[1,1]), " ",group_diffs[2,1], "\n", sep="")
  }

  if(! is.null(fit$comp_theta)) {
    cat("The prob. that the relative frequency of success is less/more than comp. val:\n")
    comp_table <- s[str_detect(rownames(s), "theta"), c("comp", "%<comp", "%>comp")]
    rownames(comp_table) <- paste("  Group ", 1:nrow(comp_table), ":", sep="")
    colnames(comp_table) <- c("comp. val.",  " <", " >") 
    print(format(comp_table, justify="centre"), quote=FALSE)
  }
  cat("\n")
}

