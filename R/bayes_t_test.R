# ...

#'Bayesian First Aid alternative to the t-test
#'
#'\code{bayes.t.test} estimates the mean of one group, or the difference in 
#'means between two groups, using Bayesian estimation and is intended as a 
#'replacement for \code{\link{t.test}}. Is based on Bayesian Estimation
#'Supersedes the t-test (BEST) (Kruschke, 2012).
#'
#'As with the \code{\link{t.test}} function \code{bayes.t.test} estimates one of
#'three models depending on the arguments given. All three models are based on 
#'the \emph{Bayesian Estimation Supersedes the t test} (BEST) model developed by
#'Kruschke (2013).
#'
#'If one vecor is supplied a one sample BEST is run. BEST assumes the data (\eqn{x})
#'is distributed as a t distribution, a more robust alternative to the normal 
#'distribution due to its wider tails. Except for the mean (\eqn{\mu}) and the
#'scale (\eqn{\sigma}) the t has one additional parameter, the
#'degree-of-freedoms (\eqn{\nu}), where the lower \eqn{\nu} is the wider the
#'tails become. When \eqn{\nu} gets larger the t distribution approaches the
#'normal distribution. While it would be possible to fix \eqn{\nu} to a single
#'value BEST instead estimates \eqn{\nu} allowing the t-distribution to become
#'more or less normal depending on the data. Here is the full model for the one
#'sample BEST:
#'
#'\deqn{x_i \sim \mathrm{t}(\mu, \sigma, \nu)}{x[i] ~ t(\mu, \sigma, \nu)} 
#'\deqn{\mu \sim \mathrm{Normal}(M_\mu, S_\mu)}{\mu ~ Normal(M[\mu], S[\mu])} 
#'\deqn{\sigma \sim \mathrm{Uniform}(L_\sigma, H__\sigma)}{\sigma ~
#'Uniform(L[\sigma], H[\sigma])} \deqn{\nu \sim \mathrm{Shifted-Exp}(1/29,
#'\mathrm{shift}=1)}{\nu ~ Shifted-Exp(1/29, shift=1)}
#'
#'\figure{one_sample_best_diagram.png}{A graphical diagram of the one sample the
#'BEST model}
#'
#'The constants \eqn{M_\mu, S_\mu, L_\sigma, H__\sigma}{M[\mu], S[\mu],
#'L[\sigma]} and \eqn{H__\sigma}{H[\sigma]} are set so that the priors on
#'\eqn{\mu} and \eqn{\sigma} are essentially flat.
#'
#'If two vectors are supplied a two sample BEST is run. This is essentially the 
#'same as estimaiting two separate one sample BEST except for that both groups 
#'are assumed to have the same \eqn{\nu}. Here is a Kruschke style diagram 
#'showing the two sample BEST model:
#'
#'\figure{two_sample_best_diagram.png}{A graphical diagram of the two sample the
#'BEST model}
#'
#'If two vectors are supplied and \code{paired=TRUE} then the paired difference
#'between \code{x - y} is modeled using the one sample BEST.
#'
#'
#'@param x a (non-empty) numeric vector of data values.
#'@param y an optional (non-empty) numeric vector of data values.
#'@param alternative ignored and is only retained in order to mantain 
#'  compatibility with \code{\link{t.test}}.
#'@param mu a fixed relative mean value to compare with the estimated mean (or 
#'  the difference in means when performing a two sample BEST).
#'@param paired a logical indicating whether you want to estimate a paired samples BEST.
#'@param var.equal ignored and is only retained in order to mantain 
#'  compatibility with \code{\link{t.test}}.
#'@param cred.mass the amount of probability mass that will be contained in 
#'  reported credible intervals. This argument fills a similar role as 
#'  \code{conf.level} in \code{\link{t.test}}.
#'@param n.iter The number of iterations to run the MCMC sampling.
#'@param progress.bar The type of progress bar. Possible values are "text", 
#'  "gui", and "none".
#'@param conf.level same as \code{cred.mass} and is only retained in order to 
#'  mantain compatibility with \code{\link{binom.test}}.
#'@param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} is a 
#'  numeric variable giving the data values and \code{rhs} a factor with two 
#'  levels giving the corresponding groups.
#'@param data an optional matrix or data frame (or similar: see 
#'  \code{\link{model.frame}}) containing the variables in the formula formula. 
#'  By default the variables are taken from \code{environment(formula)}.
#'@param subset an optional vector specifying a subset of observations to be 
#'  used.
#'@param na.action a function which indicates what should happen when the data 
#'  contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#'@param ... further arguments to be passed to or from methods.
#'  
#'  
#'@return A list of class \code{bayes_paired_t_test}, 
#'  \code{bayes_one_sample_t_test} or \code{bayes_two_sample_t_test} that 
#'  contains information about the analysis. It can be further inspected using 
#'  the functions \code{summary}, \code{plot}, \code{\link{diagnostics}} and 
#'  \code{\link{model.code}}.
#'  
#' @examples
#' # Using Student's sleep data as in the t.test example
#' # bayes.t.test can be called in the same way as t.test 
#' # so you can both supply two vectors...
#' 
#' bayes.t.test(sleep$extra[sleep$group == 1], sleep$extra[sleep$group == 2])
#' 
#' # ... or use the formula interface.
#' 
#' bayes.t.test(extra ~ group, data = sleep)
#' 
#' # Save the return value in order to inspect the model result further.
#' fit <- bayes.t.test(extra ~ group, data = sleep)
#' summary(fit)
#' plot(fit)
#' 
#' # MCMC diagnostics
#' diagnostics(fit)
#' 
#' # Print out the R code to run the model. This can be copy n' pasted into
#' # an R-script and further modified.
#' model.code(fit)
#'  
#'@references Kruschke, J. K. (2013). Bayesian estimation supersedes the t test.
#'  \emph{Journal of Experimental Psychology: General}, 142(2), 573.
#'  
#'@export
#'@rdname bayes.t.test
bayes.t.test <- function(x, ...) {
  UseMethod("bayes.t.test")
}


#' @export
#' @rdname bayes.t.test
bayes.t.test.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
                                 mu = 0, paired = FALSE, var.equal = FALSE, cred.mass = 0.95, n.iter = 30000, progress.bar="text", conf.level,...) {
  
  if(! missing(conf.level)) {
    cred.mass <- conf.level
  }
  
  if(var.equal) {
    var.equal <- FALSE
    warning("To assume equal variance of 'x' and 'y' is not supported. Continuing by estimating the variance of 'x' and 'y' separately.")
  }
  
  if(! missing(alternative)) {
    warning("The argument 'alternative' is ignored by bayes.binom.test")
  } 
  
  ### Original (but slighly modified) code from t.test.default ###
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  if (!missing(cred.mass) && (length(cred.mass) != 1 || !is.finite(cred.mass) || 
                                cred.mass < 0 || cred.mass > 1)) 
    stop("'cred.mass' or 'conf.level' must be a single number between 0 and 1")
  
  # removing incomplete cases and preparing the data vectors (x & y)
  if (!is.null(y)) {
    x_name <- deparse(substitute(x))
    y_name <- deparse(substitute(y))
    data_name <- paste(x_name, "and", y_name)
    if (paired) 
      xok <- yok <- complete.cases(x, y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else {
    x_name <- deparse(substitute(x))
    data_name <- x_name
    if (paired) 
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  
  # Checking that there is enough data. Even though BEST handles the case with
  # one data point it is still usefull to do these checks.
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  if (is.null(y)) {
    if (nx < 2) 
      stop("not enough 'x' observations")
    df <- nx - 1
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx)) 
      stop("data are essentially constant")
  }
  else {
    ny <- length(y)
    if (nx < 2) 
      stop("not enough 'x' observations")
    if (ny < 2) 
      stop("not enough 'y' observations")
    my <- mean(y)
    vy <- var(y)
    stderrx <- sqrt(vx/nx)
    stderry <- sqrt(vy/ny)
    stderr <- sqrt(stderrx^2 + stderry^2)
    df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), 
                                                abs(my))) 
      stop("data are essentially constant")
  }
  
  ### Own code starts here ###
  
  if(paired) {
    mcmc_samples <- jags_paired_t_test(x, y, n.chains= 3, n.iter = ceiling(n.iter / 3), progress.bar=progress.bar)
    stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = mu)
    bfa_object <- list(x = x, y = y, pair_diff = x - y, comp = mu, cred_mass = cred.mass,
                       x_name = x_name, y_name = y_name, data_name = data_name,
                       x_data_expr = x_name, y_data_expr = y_name,
                       mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- c("bayes_paired_t_test", "bayesian_first_aid")
    
  } else if(is.null(y)) {
    mcmc_samples <- jags_one_sample_t_test(x, comp_mu = mu, n.chains= 3, n.iter = ceiling(n.iter / 3), progress.bar=progress.bar)
    stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = mu)
    bfa_object <- list(x = x, comp = mu, cred_mass = cred.mass, x_name = x_name, x_data_expr = x_name,
                       data_name = data_name, mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- c("bayes_one_sample_t_test", "bayesian_first_aid")
    
  } else { # is two sample t.test
    mcmc_samples <- jags_two_sample_t_test(x, y, n.chains= 3, n.iter = ceiling(n.iter / 3), progress.bar=progress.bar)
    stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = mu)
    bfa_object <- list(x = x, y = y, comp = mu, cred_mass = cred.mass,
                       x_name = x_name, y_name = y_name, data_name = data_name,
                       x_data_expr = x_name, y_data_expr = y_name,
                       mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- c("bayes_two_sample_t_test", "bayesian_first_aid")
  }
  bfa_object
}


#' @export
#' @rdname bayes.t.test
bayes.t.test.formula <- function(formula, data, subset, na.action, ...) {
  
  ### Original code from t.test.formula ###
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), "term.labels")) != 1L)) 
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  data_name <- paste(names(mf), collapse = " by ")
  response_name <- names(mf)[1]
  group_name <- names(mf)[2]
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L) 
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  
  ### Own code starts here ###
  bfa_object <- do.call("bayes.t.test", c(DATA, list(...)))
  bfa_object$data_name <- data_name
  bfa_object$x_name <- paste("group", levels(g)[1])
  bfa_object$y_name <- paste("group", levels(g)[2])
  if(!missing(data)) {
    data_expr <- deparse(substitute(data))
    bfa_object$x_data_expr <- 
      paste("subset(", data_expr, ", as.factor(", group_name, ") == ",
            deparse(levels(g)[1]), ", ", response_name, ", drop = TRUE)", sep="")
    bfa_object$y_data_expr <- 
      paste("subset(", data_expr, ", as.factor(", group_name, ") == ",
            deparse(levels(g)[2]), ", ", response_name, ", drop = TRUE)", sep="")
  } else {
    bfa_object$x_data_expr <- 
      paste(response_name, "[", "as.factor(", group_name, ") == ",
            deparse(levels(g)[1]),"]",sep="")
    bfa_object$y_data_expr <- 
      paste(response_name, "[", "as.factor(", group_name, ") == ",
            deparse(levels(g)[2]),"]",sep="")
  }
  bfa_object  
}


one_sample_t_model_string <- "model {
  for(i in 1:length(x)) {
    x[i] ~ dt( mu , tau , nu )
  }
  x_pred ~ dt( mu , tau , nu )
  eff_size <- (mu - comp_mu) / sigma

  mu ~ dnorm( mean_mu , precision_mu )
  tau <- 1/pow( sigma , 2 )
  sigma ~ dunif( sigma_low , sigma_high )
  # A trick to get an exponentially distributed prior on nu that starts at 1.
  nu <- nuMinusOne + 1 
  nuMinusOne ~ dexp(1/29)
}"

jags_one_sample_t_test <- function(x, comp_mu = 0,n.adapt= 500, n.chains=3, n.update = 100, n.iter=5000, thin=1, progress.bar="text") {
  data_list <- list(
    x = x,
    mean_mu = mean(x, trim=0.2) ,
    precision_mu = 1 / (mad0(x)^2 * 1000000),
    sigma_low = mad0(x) / 1000 ,
    sigma_high = mad0(x) * 1000 ,
    comp_mu = comp_mu
  )
  
  inits_list <- list(mu = mean(x, trim=0.2), sigma = mad0(x), nuMinusOne = 4)
  params <- c("mu", "sigma", "nu", "eff_size", "x_pred")
  mcmc_samples <- run_jags(one_sample_t_model_string, data = data_list, inits = inits_list, 
                           params = params, n.chains = n.chains, n.adapt = n.adapt,
                           n.update = n.update, n.iter = n.iter, thin = thin, progress.bar=progress.bar)   
  mcmc_samples
}

two_sample_t_model_string <- "model {
  for(i in 1:length(x)) {
    x[i] ~ dt( mu_x , tau_x , nu )
  }
  x_pred ~ dt( mu_x , tau_x , nu )
  for(i in 1:length(y)) {
    y[i] ~ dt( mu_y , tau_y , nu )
  }
  y_pred ~ dt( mu_y , tau_y , nu )
  eff_size <- (mu_x - mu_y) / sqrt((pow(sigma_x, 2) + pow(sigma_y, 2)) / 2)
  mu_diff <- mu_x - mu_y
  sigma_diff <-sigma_x - sigma_y 
  
  # The priors
  mu_x ~ dnorm( mean_mu , precision_mu )
  tau_x <- 1/pow( sigma_x , 2 )
  sigma_x ~ dunif( sigma_low , sigma_high )

  mu_y ~ dnorm( mean_mu , precision_mu )
  tau_y <- 1/pow( sigma_y , 2 )
  sigma_y ~ dunif( sigma_low , sigma_high )

  # A trick to get an exponentially distributed prior on nu that starts at 1.
  nu <- nuMinusOne+1
  nuMinusOne ~ dexp(1/29)
}"

# Adapted from John Kruschke's original BEST code.
jags_two_sample_t_test <- function(x, y, n.adapt= 500, n.chains=3, n.update = 100, n.iter=5000, thin=1, progress.bar="text") {
  data_list <- list(
    x = x ,
    y = y ,
    mean_mu = mean(c(x, y), trim=0.2) ,
    precision_mu = 1 / (mad0(c(x, y))^2 * 1000000),
    sigma_low = mad0(c(x, y)) / 1000 ,
    sigma_high = mad0(c(x, y)) * 1000 
  )
  
  inits_list <- list(
    mu_x = mean(x, trim=0.2),
    mu_y = mean(y, trim=0.2),
    sigma_x = mad0(x),
    sigma_y = mad0(y),
    nuMinusOne = 4
  )
  
  params <- c("mu_x", "sigma_x", "mu_y", "sigma_y", "mu_diff", "sigma_diff","nu", "eff_size", "x_pred", "y_pred")  
  mcmc_samples <- run_jags(two_sample_t_model_string, data = data_list, inits = inits_list, 
                           params = params, n.chains = n.chains, n.adapt = n.adapt,
                           n.update = n.update, n.iter = n.iter, thin = thin, progress.bar=progress.bar)
  mcmc_samples
}

# Right now, this is basically just calling jags_one_sample_t_test but I'm 
# keeping it in case I would want to change it in the future.

paired_samples_t_model_string <- "model {
  for(i in 1:length(pair_diff)) {
    pair_diff[i] ~ dt( mu_diff , tau_diff , nu )
  }
  diff_pred ~ dt( mu_diff , tau_diff , nu )
  eff_size <- (mu_diff - comp_mu) / sigma_diff
  
  mu_diff ~ dnorm( mean_mu , precision_mu )
  tau_diff <- 1/pow( sigma_diff , 2 )
  sigma_diff ~ dunif( sigma_low , sigma_high )
  # A trick to get an exponentially distributed prior on nu that starts at 1.
  nu <- nuMinusOne + 1 
  nuMinusOne ~ dexp(1/29)
}"


jags_paired_t_test <- function(x, y, comp_mu = 0, n.adapt= 500, n.chains=3, n.update = 100, n.iter=5000, thin=1, progress.bar="text") {
  pair_diff <- x - y
  data_list <- list(
    pair_diff = pair_diff,
    mean_mu = mean(pair_diff, trim=0.2) ,
    precision_mu = 1 / (mad0(pair_diff)^2 * 1000000),
    sigma_low = mad0(pair_diff) / 1000 ,
    sigma_high = mad0(pair_diff) * 1000 ,
    comp_mu = comp_mu
  )
  
  inits_list <- list(mu_diff = mean(pair_diff, trim=0.2), 
                     sigma_diff = mad0(pair_diff), 
                     nuMinusOne = 4)
  params <- c("mu_diff", "sigma_diff", "nu", "eff_size", "diff_pred")
  mcmc_samples <- run_jags(paired_samples_t_model_string, data = data_list, inits = inits_list, 
                           params = params, n.chains = n.chains, n.adapt = n.adapt,
                           n.update = n.update, n.iter = n.iter, thin = thin, progress.bar=progress.bar) 
  mcmc_samples
}

####################################
### One sample t-test S3 methods ###
####################################

#' @export
print.bayes_one_sample_t_test <- function(x, ...) {
  
  s <- format_stats(x$stats)
  
  cat("\n")
  cat("\tBayesian estimation supersedes the t test (BEST) - one sample\n")
  cat("\n")
  cat("data: ", x$data_name,  ", n = ", length(x$x),"\n", sep="")
  cat("\n")
  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep="")
  cat("mean of ",  x$x_name, ": ", s["mu", "median"], " [", s["mu", "HDIlo"], ", ", s["mu", "HDIup"] , "]\n",sep="")
  cat("sd of ",  x$x_name, ": ", s["sigma", "median"], " [", s["sigma", "HDIlo"], ", ", s["sigma", "HDIup"] , "]\n",sep="")
  
  cat("\n")
  cat("The mean is more than", s["mu","comp"] , "by a probability of", s["mu","%>comp"], "\n")
  cat("and less than", s["mu", "comp"] , "by a probability of", s["mu", "%<comp"], "\n")
  cat("\n")
  invisible(NULL)
}

print_bayes_one_sample_t_test_params <- function(x) {
  cat("  Model parameters and generated quantities\n")
  cat("mu: the mean of", x$data_name, "\n")
  cat("sigma: the scale of", x$data_name,", a consistent\n  estimate of SD when nu is large.\n")
  cat("nu: the degrees-of-freedom for the t distribution fitted to",x$data_name , "\n")
  cat("eff_size: the effect size calculated as (mu - ", x$comp ,") / sigma\n", sep="")
  cat("x_pred: predicted distribution for a new datapoint generated as",x$data_name , "\n")
}

#' @export
summary.bayes_one_sample_t_test <- function(object, ...) {
  s <- round(object$stats, 3)
  
  cat("  Data\n")
  cat(object$data_name, ", n = ", length(object$x), "\n", sep="")
  cat("\n")
  
  print_bayes_one_sample_t_test_params(object)
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

#' @method plot bayes_one_sample_t_test
#' @export
plot.bayes_one_sample_t_test <- function(x, ...) {
  stats <- x$stats
  mcmc_samples <- x$mcmc_samples
  samples_mat <- as.matrix(mcmc_samples)
  mu = samples_mat[,"mu"]
  sigma = samples_mat[,"sigma"]
  nu = samples_mat[,"nu"]
  
  #layout( matrix( c(3,3,4,4,5,5, 1,1,1,1,2,2) , nrow=6, ncol=2 , byrow=FALSE ) )
  layout( matrix( c(2,2,3,3,4,4, 1,1) , nrow=4, ncol=2 , byrow=FALSE ) )
  old_par <- par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  
  n_curves <- 30
  rand_i <- sample(nrow(samples_mat), n_curves)
  hist_with_t_curves(x$x, stats["mu", "mean"], stats["sigma", "median"], mu[rand_i], sigma[rand_i],
                      nu[rand_i], x$data_name, main= "Data w. Post. Pred.", x_range= range(x$x))
  
  # Plot posterior distribution of parameter nu:
  #paramInfo = plotPost( log10(nu) , col="skyblue" ,
  #                      xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , show_mode=TRUE ,
  #                      main="Normality" ) #  (<0.7 suggests kurtosis)

  # distribution of mu:
  xlim = range( c( mu , x$comp ) )
  plotPost(mu ,  xlim=xlim , cex.lab = 1.75 , comp_val = x$comp, cred_mass= x$cred_mass,
           xlab=bquote(mu) , main=paste("Mean") , col="skyblue", show_median=TRUE )
  
  # distribution of sigma:
  plotPost(sigma, cex.lab = 1.75, xlab=bquote(sigma), main=paste("Std. Dev."),  
           cred_mass= x$cred_mass, col="skyblue" , show_median=TRUE )

  # effect size:
  plotPost(samples_mat[, "eff_size"] , comp_val=0 , xlab=bquote( (mu-.(x$comp)) / sigma ),
           cred_mass= x$cred_mass, show_median=TRUE , cex.lab=1.75 , main="Effect Size" , col="skyblue" )
  par(old_par)
  invisible(NULL)
}

#' @export
diagnostics.bayes_one_sample_t_test <- function(fit) {
    
  print_mcmc_info(fit$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(fit$stats, 3))
  cat("\n")
  print_bayes_one_sample_t_test_params(fit)
  cat("\n")
  
  old_par <- par( mar=c(3.5,2.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  plot(fit$mcmc_samples)
  par(old_par)
  invisible(NULL)
}

#' @export
model.code.bayes_one_sample_t_test <- function(fit) {  
  cat("### Model code for Bayesian estimation supersedes the t test - one sample ###\n")
  
  cat("require(rjags)\n\n")
  cat("# Setting up the data\n")
  cat("x <-", fit$x_data_expr, "\n")
  cat("comp_mu <- ", fit$comp, "\n")
  cat("\n")
  pretty_print_function_body(one_sample_t_test_model_code)
  invisible(NULL)
}

# Not to be run, just to be printed
one_sample_t_test_model_code <- function(x, comp_mu) {
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string
  
  # Setting parameters for the priors that in practice will result
  # in flat priors on mu and sigma.
  mean_mu = mean(x, trim=0.2)
  precision_mu = 1 / (mad(x)^2 * 1000000)
  sigma_low = mad(x) / 1000 
  sigma_high = mad(x) * 1000
  
  # Initializing parameters to sensible starting values helps the convergence
  # of the MCMC sampling. Here using robust estimates of the mean (trimmed)
  # and standard deviation (MAD).
  inits_list <- list(mu = mean(x, trim=0.2), sigma = mad(x), nuMinusOne = 4)
  
  data_list <- list(
    x = x,
    comp_mu = comp_mu,
    mean_mu = mean_mu,
    precision_mu = precision_mu,
    sigma_low = sigma_low,
    sigma_high = sigma_high)
  
  # The parameters to monitor.
  params <- c("mu", "sigma", "nu", "eff_size", "x_pred")
  
  # Running the model
  model <- jags.model(textConnection(model_string), data = data_list,
                      inits = inits_list, n.chains = 3, n.adapt=1000)
  update(model, 500) # Burning some samples to the MCMC gods....
  samples <- coda.samples(model, params, n.iter=10000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)  
}
one_sample_t_test_model_code <- inject_model_string(one_sample_t_test_model_code, one_sample_t_model_string)



####################################
### Two sample t-test S3 methods ###
####################################
#' @export
print.bayes_two_sample_t_test <- function(x, ...) {
  s <- format_stats(x$stats)
  
  cat("\n")
  cat("\tBayesian estimation supersedes the t test (BEST) - two sample\n")
  cat("\n")
  cat("data: ", x$x_name, " (n = ", length(x$x) ,") and ", x$y_name," (n = ", length(x$y) ,")\n", sep="")
  cat("\n")
  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep="")
  cat("mean of ",  x$x_name, ": ", s["mu_x", "median"], " [", s["mu_x", "HDIlo"], ", ", s["mu_x", "HDIup"] , "]\n",sep="")
  cat("mean of ",  x$y_name, ": ", s["mu_y", "median"], " [", s["mu_y", "HDIlo"], ", ", s["mu_y", "HDIup"] , "]\n",sep="")
  cat("difference of the means: ", s["mu_diff", "median"], " [", s["mu_diff", "HDIlo"], ", ", s["mu_diff", "HDIup"] , "]\n",sep="")
  cat("sd of ",  x$x_name, ": ", s["sigma_x", "median"], " [", s["sigma_x", "HDIlo"], ", ", s["sigma_x", "HDIup"] , "]\n",sep="")
  cat("sd of ",  x$y_name, ": ", s["sigma_y", "median"], " [", s["sigma_y", "HDIlo"], ", ", s["sigma_y", "HDIup"] , "]\n",sep="")
  
  cat("\n")
  cat("The difference of the means is greater than", s["mu_diff","comp"] , "by a probability of", s["mu_diff","%>comp"], "\n")
  cat("and less than", s["mu_diff", "comp"] , "by a probability of", s["mu_diff", "%<comp"], "\n")
  cat("\n")
  invisible(NULL)
}
 

print_bayes_two_sample_t_test_params <- function(x) {
  cat("  Model parameters and generated quantities\n")
  cat("mu_x: the mean of", x$x_name, "\n")
  cat("sigma_x: the scale of", x$x_name,", a consistent\n  estimate of SD when nu is large.\n")
  cat("mu_y: the mean of", x$y_name, "\n")
  cat("sigma_y: the scale of", x$y_name,"\n")
  cat("mu_diff: the difference in means (mu_x - mu_y)\n")
  cat("sigma_diff: the difference in scale (sigma_x - sigma_y)\n")
  cat("nu: the degrees-of-freedom for the t distribution\n")
  cat("  fitted to",x$data_name , "\n")
  cat("eff_size: the effect size calculated as \n", sep="")
  cat("  (mu_x - mu_y) / sqrt((sigma_x^2 + sigma_y^2) / 2)\n", sep="")
  cat("x_pred: predicted distribution for a new datapoint\n")
  cat("  generated as",x$x_name , "\n")
  cat("y_pred: predicted distribution for a new datapoint\n")
  cat("  generated as",x$y_name , "\n")
  invisible(NULL)
}

#' @export
summary.bayes_two_sample_t_test <- function(object, ...) {
  s <- round(object$stats, 3)
  
  cat("  Data\n")
  cat(object$x_name, ", n = ", length(object$x), "\n", sep="")
  cat(object$y_name, ", n = ", length(object$y), "\n", sep="")
  cat("\n")
  
  print_bayes_two_sample_t_test_params(object)
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

#' @method plot bayes_two_sample_t_test
#' @export
plot.bayes_two_sample_t_test <- function(x, ...) {
  stats <- x$stats
  mcmc_samples <- x$mcmc_samples
  samples_mat <- as.matrix(mcmc_samples)
  mu_x = samples_mat[,"mu_x"]
  sigma_x = samples_mat[,"sigma_x"]
  mu_y = samples_mat[,"mu_y"]
  sigma_y = samples_mat[,"sigma_y"]
  nu = samples_mat[,"nu"]
  
  #layout( matrix( c(4,5,7,8,3,1,2,6,9,10) , nrow=5, byrow=FALSE ) )
  layout( matrix( c(3,4,5,1,2,6) , nrow=3, byrow=FALSE ) )
  old_par <- par( mar=c(3.5,3.5,2.5,0.51) , mgp=c(2.25,0.7,0) )
  
  
  # Plot data with post predictive distribution
  n_curves <- 30
  data_range <- range(c(x$x, x$y))
  rand_i <- sample(nrow(samples_mat), n_curves)
  hist_with_t_curves(x$x, stats["mu_x", "mean"], stats["sigma_x", "median"], mu_x[rand_i], sigma_x[rand_i],
                     nu[rand_i], x$x_name, main= paste("Data",  x$x_name, "w. Post. Pred."), x_range= data_range)
  hist_with_t_curves(x$y, stats["mu_y", "mean"], stats["sigma_y", "median"], mu_y[rand_i], sigma_y[rand_i],
                     nu[rand_i], x$y_name, main= paste("Data",  x$y_name, "w. Post. Pred."), x_range= data_range)
  
  # Plot posterior distribution of parameter nu:
#   plotPost( log10(nu) , col="skyblue" , cred_mass= x$cred_mass,
#             xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , show_mode=TRUE ,
#             main="Normality" ) #  (<0.7 suggests kurtosis)
  
  # Plot posterior distribution of parameters mu_x, mu_y, and their difference:
  xlim = range( c( mu_x , mu_y ) )
  plotPost( mu_x ,  xlim=xlim , cex.lab = 1.75 , cred_mass= x$cred_mass, show_median=TRUE,
            xlab=bquote(mu[x]) , main=paste(x$x_name,"Mean") , col="skyblue" )
  plotPost( mu_y ,  xlim=xlim , cex.lab = 1.75 ,  cred_mass= x$cred_mass, show_median=TRUE,
            xlab=bquote(mu[y]) , main=paste(x$y_name,"Mean") , col="skyblue" )
  plotPost( samples_mat[,"mu_diff"] , comp_val= x$comp , cred_mass= x$cred_mass,
            xlab=bquote(mu[x] - mu[y]) , cex.lab = 1.75 , show_median=TRUE,
            main="Difference of Means" , col="skyblue" )
  
  # Save this to var.test
  #
  # Plot posterior distribution of param's sigma_x, sigma_y, and their difference:
#   xlim=range( c( sigma_x , sigma_y ) )
#   plotPost( sigma_x ,  xlim=xlim , cex.lab = 1.75 ,cred_mass= x$cred_mass,
#             xlab=bquote(sigma[x]) , main=paste(x$x_name, "Std. Dev.") , 
#             col="skyblue" , show_mode=TRUE )
#   plotPost( sigma_y ,  xlim=xlim , cex.lab = 1.75 ,cred_mass= x$cred_mass,
#             xlab=bquote(sigma[y]) , main=paste(x$y_name, "Std. Dev.") , 
#             col="skyblue" , show_mode=TRUE )
#   plotPost( samples_mat[, "sigma_diff"] , comp_val= x$comp , cred_mass= x$cred_mass, 
#             xlab=bquote(sigma[x] - sigma[y]) , cex.lab = 1.75 ,
#             main="Difference of Std. Dev.s" , col="skyblue" , show_mode=TRUE )
  
  # Plot of estimated effect size. Effect size is d-sub-a from 
  # Macmillan & Creelman, 1991; Simpson & Fitter, 1973; Swets, 1986a, 1986b.
  plotPost( samples_mat[, "eff_size"] , comp_val=0 , cred_mass= x$cred_mass,
            xlab=bquote( (mu[x]-mu[y]) / sqrt((sigma[x]^2 +sigma[y]^2 )/2 ) ),
            show_median=TRUE , cex.lab=1.0 , main="Effect Size" , col="skyblue" )
  
  par(old_par)
}

#' @export
diagnostics.bayes_two_sample_t_test <- function(fit) {
  print_mcmc_info(fit$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(fit$stats, 3))
  cat("\n")
  print_bayes_two_sample_t_test_params(fit)
  cat("\n")
  
  old_par <- par( mar=c(3.5,2.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  plot(fit$mcmc_samples)
  par(old_par)
  invisible(NULL)
}

#' @export
model.code.bayes_two_sample_t_test <- function(fit) {
  cat("## Model code for Bayesian estimation supersedes the t test - two sample ##\n")
  
  cat("require(rjags)\n\n")
  cat("# Setting up the data\n")
  cat("x <-", fit$x_data_expr, "\n")
  cat("y <-", fit$y_data_expr, "\n")
  cat("\n")
  pretty_print_function_body(two_sample_t_test_model_code)
  invisible(NULL)
}

# Not to be run, just to be printed
two_sample_t_test_model_code <- function(x, y) {
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string
  
  # Setting parameters for the priors that in practice will result
  # in flat priors on the mu's and sigma's.
  mean_mu = mean( c(x, y), trim=0.2)
  precision_mu = 1 / (mad( c(x, y) )^2 * 1000000)
  sigma_low = mad( c(x, y) ) / 1000 
  sigma_high = mad( c(x, y) ) * 1000
  
  # Initializing parameters to sensible starting values helps the convergence
  # of the MCMC sampling. Here using robust estimates of the mean (trimmed)
  # and standard deviation (MAD).
  inits_list <- list(
    mu_x = mean(x, trim=0.2), mu_y = mean(y, trim=0.2),
    sigma_x = mad(x), sigma_y = mad(y),
    nuMinusOne = 4)
  
  data_list <- list(
    x = x, y = y,    
    mean_mu = mean_mu,
    precision_mu = precision_mu,
    sigma_low = sigma_low,
    sigma_high = sigma_high)
  
  # The parameters to monitor.
  params <- c("mu_x", "mu_y", "mu_diff", "sigma_x", "sigma_y", "sigma_diff",
              "nu", "eff_size", "x_pred", "y_pred")
  
  # Running the model
  model <- jags.model(textConnection(model_string), data = data_list,
                      inits = inits_list, n.chains = 3, n.adapt=1000)
  update(model, 500) # Burning some samples to the MCMC gods....
  samples <- coda.samples(model, params, n.iter=10000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)  
}
two_sample_t_test_model_code <- inject_model_string(two_sample_t_test_model_code, two_sample_t_model_string)

########################################
### Paired samples t-test S3 methods ###
########################################


#' @export
print.bayes_paired_t_test <- function(x, ...) {
  s <- format_stats(x$stats)
  cat("\n")
  cat("\tBayesian estimation supersedes the t test (BEST) - paired samples\n")
  cat("\n")
  cat("data: ", x$data_name,  ", n = ", length(x$pair_diff),"\n", sep="")
  cat("\n")
  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep="")
  cat("mean paired difference: ", s["mu_diff", "median"], " [", s["mu_diff", "HDIlo"], ", ", s["mu_diff", "HDIup"] , "]\n",sep="")
  cat("sd of the paired differences: ", s["sigma_diff", "median"], " [", s["sigma_diff", "HDIlo"], ", ", s["sigma_diff", "HDIup"] , "]\n",sep="")
  
  cat("\n")
  cat("The mean difference is more than", s["mu_diff","comp"] , "by a probability of", s["mu_diff","%>comp"], "\n")
  cat("and less than", s["mu_diff", "comp"] , "by a probability of", s["mu_diff", "%<comp"], "\n")
  cat("\n")
  invisible(NULL)
}


print_bayes_paired_t_test_params <- function(x) {

  cat("  Model parameters and generated quantities\n")
  cat("mu_diff: the mean pairwise difference between", x$x_name, "and", x$y_name, "\n")
  cat("sigma_diff: the scale of the pairwise difference, a consistent\n  estimate of SD when nu is large.\n")
  cat("nu: the degrees-of-freedom for the t distribution fitted to the pairwise difference\n")
  cat("eff_size: the effect size calculated as (mu_diff - ", x$comp ,") / sigma_diff\n", sep="")
  cat("diff_pred: predicted distribution for a new datapoint generated\n  as the pairwise difference between", x$x_name, "and", x$y_name,"\n")
}

#' @export
summary.bayes_paired_t_test <- function(object, ...) {
  s <- round(object$stats, 3)
  
  cat("  Data\n")
  cat(object$x_name, ", n = ", length(object$x), "\n", sep="")
  cat(object$y_name, ", n = ", length(object$y), "\n", sep="")
  cat("\n")
  
  print_bayes_paired_t_test_params(object)
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

#' @method plot bayes_paired_t_test
#' @export
plot.bayes_paired_t_test <- function(x, y, ...) {
  stats <- x$stats
  mcmc_samples <- x$mcmc_samples
  samples_mat <- as.matrix(mcmc_samples)
  mu_diff = samples_mat[,"mu_diff"]
  sigma_diff = samples_mat[,"sigma_diff"]
  nu = samples_mat[,"nu"]
  
  #layout( matrix( c(3,3,4,4,5,5, 1,1,1,1,2,2) , nrow=6, ncol=2 , byrow=FALSE ) )
  layout( matrix( c(2,2,3,3,4,4, 1,1) , nrow=4, ncol=2 , byrow=FALSE ) )
  old_par <- par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  
  n_curves <- 30
  rand_i <- sample(nrow(samples_mat), n_curves)
  hist_with_t_curves(x$pair_diff, stats["mu_diff", "mean"], stats["sigma_diff", "median"], mu_diff[rand_i], sigma_diff[rand_i],
                     nu[rand_i], x$data_name, main= "Data w. Post. Pred.", x_range= range(x$pair_diff))
  
  # Plot posterior distribution of parameter nu:
  #paramInfo = plotPost( log10(nu) , col="skyblue" ,
  #                      xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , show_mode=TRUE ,
  #                      main="Normality" ) #  (<0.7 suggests kurtosis)
  
  # distribution of mu_diff:
  xlim = range( c( mu_diff , x$comp ) )
  plotPost(mu_diff ,  xlim=xlim , cex.lab = 1.75 , comp_val = x$comp, cred_mass= x$cred_mass,
           xlab=bquote(mu[diff]) , main=paste("Mean difference") , col="skyblue", show_median=TRUE )
  
  # distribution of sigma_diff:
  plotPost(sigma_diff, cex.lab = 1.75, xlab=bquote(sigma[diff]), main=paste("Std. Dev. of difference"),  
           cred_mass= x$cred_mass, col="skyblue" , show_median=TRUE )
  
  # effect size:
  plotPost(samples_mat[, "eff_size"] , comp_val=0 , xlab=bquote( (mu[diff] -.(x$comp)) / sigma[diff] ),
           cred_mass= x$cred_mass, show_median=TRUE , cex.lab=1.75 , main="Effect Size" , col="skyblue" )
  par(old_par)
  invisible(NULL)
}

#' @export
diagnostics.bayes_paired_t_test <- function(fit) {
  print_mcmc_info(fit$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(fit$stats, 3))
  cat("\n")
  print_bayes_paired_t_test_params(fit)
  cat("\n")
  
  old_par <- par( mar=c(3.5,2.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  plot(fit$mcmc_samples)
  par(old_par)
  invisible(NULL)
}

#' @export
model.code.bayes_paired_t_test <- function(fit) {
  cat("## Model code for Bayesian estimation supersedes the t test - paired samples ##\n")
  
  cat("require(rjags)\n\n")
  cat("# Setting up the data\n")
  cat("x <-", fit$x_data_expr, "\n")
  cat("y <-", fit$y_data_expr, "\n")
  cat("pair_diff <- x - y\n")
  cat("comp_mu <- ", fit$comp, "\n")
  cat("\n")
  pretty_print_function_body(paired_samples_t_test_model_code)
  invisible(NULL)
}

# Not to be run, just to be printed
paired_samples_t_test_model_code <- function(pair_diff, comp_mu) {
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string
  
  # Setting parameters for the priors that in practice will result
  # in flat priors on mu and sigma.
  mean_mu = mean(pair_diff, trim=0.2)
  precision_mu = 1 / (mad(pair_diff)^2 * 1000000)
  sigma_low = mad(pair_diff) / 1000 
  sigma_high = mad(pair_diff) * 1000
  
  # Initializing parameters to sensible starting values helps the convergence
  # of the MCMC sampling. Here using robust estimates of the mean (trimmed)
  # and standard deviation (MAD).
  inits_list <- list(
    mu_diff = mean(pair_diff, trim=0.2),
    sigma_diff = mad(pair_diff),
    nuMinusOne = 4)
  
  data_list <- list(
    pair_diff = pair_diff,
    comp_mu = comp_mu,
    mean_mu = mean_mu,
    precision_mu = precision_mu,
    sigma_low = sigma_low,
    sigma_high = sigma_high)
  
  # The parameters to monitor.
  params <- c("mu_diff", "sigma_diff", "nu", "eff_size", "diff_pred")
  
  # Running the model
  model <- jags.model(textConnection(model_string), data = data_list,
                      inits = inits_list, n.chains = 3, n.adapt=1000)
  update(model, 500) # Burning some samples to the MCMC gods....
  samples <- coda.samples(model, params, n.iter=10000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)  
}
paired_samples_t_test_model_code <- inject_model_string(paired_samples_t_test_model_code, paired_samples_t_model_string)
