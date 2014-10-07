#'Bayesian First Aid Alternative to a Test of Proportions
#'
#'\code{bayes.prop.test} estimates the relative frequency of success for two or 
#'more groups using Bayesian estimation and is intended as a replacement for 
#'\code{\link{prop.test}}.
#'
#'Given data on the number of successes and failures \code{bayes.prop.test} 
#'estimates \ifelse{latex}{\eqn{\theta_{1...m}}}{\eqn{\theta}₁…ₘ}, the relative
#'frequencies of success for each of the \eqn{m} groups. The following model is
#'assumed for each group:
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
#'The \code{print} and \code{plot} function will only work well with a small 
#'number of groups (2 to 6). If you have more groups you might want to run the 
#'model with a small number of groups and then print the model code using 
#'\code{\link{model.code}} and fit the model using that code. The current model does
#'not assume any dependency between the groups, if this is an unreasonable assumption
#'you might want to modify the model code (from \code{\link{model.code}}) to 
#'include a dependency between the groups (see 
#'\href{http://lingpipe-blog.com/2009/09/23/bayesian-estimators-for-the-beta-binomial-model-of-batting-ability/}{here}
#'for an example).
#'
#'@param x a vector of counts of successes, a one-dimensional table with two 
#'  entries, or a two-dimensional table (or matrix) with 2 columns, giving the 
#'  counts of successes and failures, respectively.
#'@param n a vector of counts of trials; ignored if x is a matrix or a table.
#'@param comp.theta a vector of fixed relative frequencies of success to compare
#'  with the estimated relative frequency of success. The length of 
#'  \code{comp.theta} must be the same as the number of groups specified by 
#'  \code{x}, and its elements must be greater than 0 and less than 1. This 
#'  argument fills a similar role as \code{p} in \code{\link{prop.test}}.
#'@param alternative ignored and is only retained in order to mantain 
#'  compatibility with \code{\link{prop.test}}.
#'@param cred.mass the amount of probability mass that will be contained in 
#'  reported credible intervals. This argument fills a similar role as 
#'  \code{conf.level} in \code{\link{prop.test}}.
#'@param correct ignored and is only retained in order to mantain compatibility 
#'  with \code{\link{prop.test}}
#'@param n.iter The number of iterations to run the MCMC sampling.
#'@param progress.bar The type of progress bar. Possible values are "text", 
#'  "gui", and "none".
#'@param p same as \code{comp.theta} and is only retained in order to mantain 
#'  compatibility with \code{\link{prop.test}}.
#'@param conf.level same as \code{cred.mass} and is only retained in order to 
#'  mantain compatibility with \code{\link{prop.test}}.
#'  
#'  
#'@return A list of class \code{bayes_prop_test} that contains information about
#'  the analysis. It can be further inspected using the functions 
#'  \code{summary}, \code{plot}, \code{\link{diagnostics}} and 
#'  \code{\link{model.code}}.
#'  
#' @examples
#' 
#' 
#' # Data from Muller, F. H., Tobakmissbrauch und Lungencarcinom,
#' # Zeit. f. Krebsforsch. 49, 57-85, 1939. One of the early papers
#' # investigating the relation between smoking and lung cancer.
#'
#' # Number of heavy smokers in one group of 86 lung cancer patients
#' # and one group of 86 healthy individuals.
#' no_heavy_smokers <- c(56, 31)
#' no_cases <- c(86, 86)
#'
#' bayes.prop.test(no_heavy_smokers, no_cases)
#' 
#' # Save the return value in order to inspect the model result further.
#' fit <- bayes.prop.test(no_heavy_smokers, no_cases)
#' summary(fit)
#' plot(fit)
#' 
#' # MCMC diagnostics (should not be necessary for such a simple model)
#' diagnostics(fit)
#' 
#' # Print out the R code to run the model. This can be copy'n'pasted into
#' # an R-script and further modified.
#' model.code(fit)
#' 
#' 
#'@seealso \code{\link{bayes.binom.test}} for when you want to estimate the 
#'  relative frequency for only one group.
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
    temp_comp_val <- 0.5
  } else {
    temp_comp_val <- comp.theta
  }
  stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = temp_comp_val)
  diff_stats <- mcmc_stats(create_theta_diff_matrix(as.matrix(mcmc_samples)))
  stats <- rbind(stats, diff_stats)
  bfa_object <- list(x = x, n = n, comp_theta = comp.theta, cred_mass = cred.mass,
                     x_name = x_name, n_name = n_name, data_name = DNAME,
                     mcmc_samples = mcmc_samples, stats = stats) 
  class(bfa_object) <- c("bayes_prop_test", "bayesian_first_aid")
  bfa_object
}

create_theta_diff_matrix <- function(samples_mat) {
  n_groups <- sum(str_count(colnames(samples_mat), "theta\\["))
  combs <- combn(n_groups, 2)
  theta_diffs <- sapply(1:ncol(combs), function(comb_i) {
    i <- combs[1, comb_i]
    j <- combs[2, comb_i]
    theta_diff <- samples_mat[,paste0("theta[", i,"]")] - samples_mat[,paste0("theta[", j,"]")] 
    theta_diff <- matrix(theta_diff, nrow = 1, dimnames = NULL)
    theta_diff
  })
  
  colnames(theta_diffs) <- apply(combs, 2, function(comb) {paste0("theta_diff[", comb[1], ",", comb[2], "]")})
  theta_diffs
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
  s <- bfa_object$stats
  n_groups <- length(bfa_object$x)
  med_diff_mat <- matrix("", nrow = n_groups, ncol = n_groups)
  hdi_diff_mat <- matrix("", nrow = n_groups, ncol = n_groups)
  diff_names <- rownames(s)[ str_detect(rownames(s), "theta_diff\\[")]
  for(diff_name in diff_names) {
    indices_match <- str_match(diff_name, "\\[(\\d+),(\\d+)\\]$")
    i <- as.numeric(indices_match[1,2])
    j <- as.numeric(indices_match[1,3])
    med_diff_mat[i, j] <- as.character(round(s[diff_name, "median"], 2) )
    hdi_diff_mat[i, j] <- paste("[", signif(s[diff_name, "HDIlo"], 2), ", ", signif(s[diff_name, "HDIup"], 2), "]", sep="")
  }
  diff_mat <- matrix("", nrow = n_groups * 2, ncol = n_groups)
  for(i in seq_len(n_groups)) {
    diff_mat[1 + (i - 1) * 2,] <- med_diff_mat[i,] 
    diff_mat[2 + (i - 1) * 2,] <- hdi_diff_mat[i,] 
  }
  rownames(diff_mat) <- rep(seq_len(n_groups), each = 2)
  rownames(diff_mat)[1:nrow(diff_mat) %% 2 == 0] <- ""
  rownames(diff_mat) <- paste0("  ", rownames(diff_mat))
  diff_mat <- format(diff_mat, width = max(nchar(diff_mat)), justify = "centre")
  colnames(diff_mat) <- format(as.character(1:ncol(diff_mat)), width = max(nchar(diff_mat)), justify = "centre")
  diff_mat <- diff_mat[-c(nrow(diff_mat) - 1, nrow(diff_mat)), -1, drop=FALSE]
  diff_mat
}


# TODO
#' @method plot bayes_prop_test
#' @export
plot.bayes_prop_test <- function(x, ...) {
  samples <- as.matrix(x$mcmc_samples)
  # Throw away everything except the what we want to plot, the theta samples.
  samples <- samples[,str_detect(colnames(samples), "^theta\\[")]
  n_groups <- length(x$x)
  diff_samples <- create_theta_diff_matrix(as.matrix(x$mcmc_samples)) 
  layout_mat <- matrix( 0 , nrow=n_groups, ncol=n_groups)
  #layout_mat[,1] <- seq_len(n_groups)
  diag(layout_mat) <- seq_len(n_groups)
  
  old_par <- par(no.readonly = TRUE)
  layout_mat <- t(layout_mat)
  layout_mat[lower.tri(layout_mat)] <- seq(n_groups + 1, by = 2,length.out = (ncol(diff_samples)))
  layout_mat <- t(layout_mat)
  layout_mat[lower.tri(layout_mat)] <- seq(n_groups + 2, by = 2,length.out = (ncol(diff_samples)))
  layout(layout_mat)
  par( mar=c(3.5,2,2,2) , mgp=c(2.25,0.7,0) )
  post_xlim <- range(apply(samples, 2, quantile, probs = c(0.001, 0.999)))
  # Some rules for making the post_xlim nice, with a preference for showing endpoints of the scale
  xlim_length <- abs(diff(post_xlim))
  if( post_xlim[1] - xlim_length < 0) {
    post_xlim[1] <- 0
  }
  if(post_xlim[2] + xlim_length > 1) {
    post_xlim[2] <- 1
  }
  plotPost(samples[,"theta[1]"], cex.lab = 1.5, xlab=bquote(theta[1]), main=paste("Rel. Freq. Group 1"),  
           cred_mass= x$cred_mass, col="#5DE293" , show_median=TRUE, comp_val=x$comp_theta[1], xlim=post_xlim)
  for(i in 2:n_groups) {
    plotPost(samples[,paste0("theta[",i, "]")], cex.lab = 1.5, xlab=bquote(theta[.(i)]), main=paste("Group", i),  
             cred_mass= x$cred_mass, col="#5DE293" , show_median=TRUE, comp_val=x$comp_theta[i], xlim=post_xlim, show_labels = FALSE)
  }
  diff_xlim <- range(apply(diff_samples, 2, quantile, probs = c(0.001, 0.999)))
  if(all(diff_xlim < 0)) {
    diff_xlim[2] <- 0
  } else if(all(diff_xlim > 0)) {
    diff_xlim[1] <- 0
  }
  for(i in 1:ncol(diff_samples)) {
    diff_name <- colnames(diff_samples)[i]
    indices_match <- str_match(diff_name, "\\[(\\d+),(\\d+)\\]$")
    group_i <- as.numeric(indices_match[1,2])
    group_j <- as.numeric(indices_match[1,3])
    plotPost(diff_samples[,i], cex.lab = 1.5, xlab=bquote(theta[.(group_i)] - theta[.(group_j)]),
             main="", cred_mass= x$cred_mass, col="skyblue" , show_median=TRUE, 
             comp_val=0, xlim=diff_xlim, show_labels = FALSE)
    plotPost(-diff_samples[,i], cex.lab = 1.5, xlab=bquote(theta[.(group_j)] - theta[.(group_i)]),
             main="", cred_mass= x$cred_mass, col="skyblue" , show_median=TRUE, 
             comp_val=0, xlim=sort(-diff_xlim), show_labels = FALSE)
  }

  par(old_par)
  invisible(NULL)
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
  for(param_i in which(str_detect(rownames(s), "theta\\["))) {
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

  if(! is.null(x$comp_theta)) {
    cat("The prob. that the relative frequency of success is less/more than comp. val:\n")
    comp_table <- s[str_detect(rownames(s), "theta\\["), c("comp", "%<comp", "%>comp")]
    rownames(comp_table) <- paste("  Group ", 1:nrow(comp_table), ":", sep="")
    colnames(comp_table) <- c("comp. val.",  " <", " >") 
    print(format(comp_table, justify="centre"), quote=FALSE)
  }
  if(ncol(group_diffs) == 1) {
    cat("The relative frequency of success is larger for Group 1 by a probability\n")
    cat("of", s["theta_diff[1,2]", "%>comp"], "and larger for Group 2 by a probability of", s["theta_diff[1,2]", "%<comp"], ".\n")
  }
  cat("\n")
  invisible(NULL)
}

#' @export
summary.bayes_prop_test <- function(object, ...) {
  
  s <- round(object$stats, 3)
  
  cat("  Data\n")
  pad_width <- max(nchar(as.character(c(object$x, object$n)))) + 1
  cat("number of successes: ", paste(str_pad(object$x, pad_width), collapse = ","), "\n", sep="")
  cat("number of trials:    ", paste(str_pad(object$n, pad_width), collapse = ","), "\n", sep="")
  cat("\n")
  
  cat("  Model parameters and generated quantities\n")
  cat("theta[i]: the relative frequency of success for Group i\n")
  cat("x_pred[i]: predicted number of successes in a replication for Group i\n")
  cat("theta_diff[i,j]: the difference between two groups (theta[i] - theta[j])\n")
  cat("\n")
  cat("  Measures\n" )
  print(s[, c("mean", "sd", "HDIlo", "HDIup", "%<comp", "%>comp")])
  cat("\n")
  cat("'HDIlo' and 'HDIup' are the limits of a ", s[1, "HDI%"] ,"% HDI credible interval.\n", sep="")
  cat("'%<comp' and '%>comp' are the probabilities of the respective parameter being\n")
  cat("smaller or larger than ", s[1, "comp"] ,".\n", sep="")
  
  cat("\n")
  cat("  Quantiles\n" )
  print(s[, c("q2.5%", "q25%", "median","q75%", "q97.5%")] )
  invisible(NULL)
}

#' @export
diagnostics.bayes_prop_test <- function(fit) {
  print_mcmc_info(fit$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(fit$stats, 3))
  cat("\n")
  
  cat("  Model parameters and generated quantities\n")
  cat("theta: The relative frequency of success\n")
  cat("x_pred: Predicted number of successes in a replication\n")
  cat("theta_diff[i,j]: the difference between two groups (theta[i] - theta[j])\n")
  old_par <- par( mar=c(3.5,2.5,2.5,0.6) , mgp=c(2.25,0.7,0) )
  plot(fit$mcmc_samples)
  par(old_par)
  invisible(NULL)
}

# Model code for the Bayesian First Aid alternative to the test of proportions #

#' @export
model.code.bayes_prop_test <- function(fit) {
  cat("### Model code for the Bayesian First Aid  ###\n### alternative to the test of proportions ###\n")
  cat("require(rjags)\n\n")
  
  cat("# Setting up the data\n")
  cat("x <-", deparse(fit$x, ), "\n")
  cat("n <-", deparse(fit$n), "\n")
  cat("\n")
  pretty_print_function_body(prop_model_code)
  invisible(NULL)
}

# Not to be run, just to be printed
prop_model_code <- function(x, n) {
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string
  
  # Running the model
  model <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                      n.chains = 3, n.adapt=1000)
  samples <- coda.samples(model, c("theta", "x_pred"), n.iter=5000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)
  
  # You can extract the mcmc samples as a matrix and compare the thetas 
  # of the groups. For example, the following shows the median and 95%
  # credible interval for the difference between Group 1 and Group 2.
  samp_mat <- as.matrix(samples)
  quantile(samp_mat[, "theta[1]"] - samp_mat[, "theta[2]"], c(0.025, 0.5, 0.975))
}
prop_model_code <- inject_model_string(prop_model_code, prop_model_string)