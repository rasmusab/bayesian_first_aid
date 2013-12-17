# ...


jags_one_sample_t_test <- function(x, comp_mu = 0,n.adapt= 1000, n.chains=3, n.update = 1000, n.iter=5000, thin=1) {
  model_string <- "
    model {
      for(i in 1:length(x)) {
        x[i] ~ dt( mu , tau , nu )
      }
      x_pred ~ dt( mu , tau , nu )
      eff_size <- (mu - comp_mu) / sigma
  
      mu ~ dnorm( mean_mu , precision_mu )
      tau <- 1/pow( sigma , 2 )
      sigma ~ dunif( sigmaLow , sigmaHigh )
      nu <- nuMinusOne+1
      nuMinusOne ~ dexp(1/29)
    }"

  data_list <- list(
    x = x,
    mean_mu = mean(x) ,
    precision_mu = 1 / (sd(x)^2 * 1000000),
    sigmaLow = sd(x) / 1000 ,
    sigmaHigh = sd(x) * 1000 ,
    comp_mu = comp_mu
  )
  
  
  inits_list <- list(mu= mean(x, trim=0.2), sigma = mad(x), nuMinusOne = 4)
  params <- c("mu", "sigma", "nu", "eff_size", "x_pred")
  mcmc_samples <- run_jags(model_string, data = data_list, inits = inits_list, 
                           params = params, n.chains = n.chains, n.adapt = n.adapt,
                           n.update = n.update, n.iter = n.iter, thin = thin)   
  mcmc_samples
}

#' Adapted from John Kruschke's original BEST code.
jags_two_sample_t_test <- function(x, y, n.adapt= 1000, n.chains=3, n.update = 1000, n.iter=5000, thin=1) {
  model_string <- "
    model {
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
      
      mu_x ~ dnorm( mean_mu , precision_mu )
      tau_x <- 1/pow( sigma_x , 2 )
      sigma_x ~ dunif( sigmaLow , sigmaHigh )

      mu_y ~ dnorm( mean_mu , precision_mu )
      tau_y <- 1/pow( sigma_y , 2 )
      sigma_y ~ dunif( sigmaLow , sigmaHigh )

      nu <- nuMinusOne+1
      nuMinusOne ~ dexp(1/29)
    }"

  data_list <- list(
    x = x ,
    y = y ,
    mean_mu = mean(c(x, y)) ,
    precision_mu = 1 / (sd(c(x, y))^2 * 1000000),
    sigmaLow = sd(c(x, y)) / 1000 ,
    sigmaHigh = sd(c(x, y)) * 1000 
  )
  
  inits_list <- list(
    mu_x = mean(x, trim=0.2),
    mu_y = mean(y, trim=0.2),
    sigma_x = mad(x),
    sigma_y = mad(y),
    nuMinusOne = 4 
  )
  
  params <- c("mu_x", "sigma_x", "mu_y", "sigma_y", "mu_diff", "sigma_diff","nu", "eff_size", "x_pred", "y_pred")  
  mcmc_samples <- run_jags(model_string, data = data_list, inits = inits_list, 
                           params = params, n.chains = n.chains, n.adapt = n.adapt,
                           n.update = n.update, n.iter = n.iter, thin = thin)
  mcmc_samples
}

# Right now, this is basically just calling jags_one_sample_t_test but I'm 
# keeping it in case I would want to change it in the future.
jags_paired_t_test <- function(x, y, comp_mu = 0,n.adapt= 1000, n.chains=3, n.update = 1000, n.iter=5000, thin=1) {
  if(is.null(y)) { # assume x is the aldread calculated difference between the two groups
    pair_diff <- x
  } else {
    pair_diff <- y - x
  }
  mcmc_samples <- 
    jags_one_sample_t_test(pair_diff, comp_mu=comp_mu,n.adapt=n.adapt, n.chains=n.chains, 
                           n.update=n.update, n.iter = n.iter, thin=thin)
  # Renaming the parameters to match a paired test
  for(i in seq_along(mcmc_samples)) {
    cnames <- colnames(mcmc_samples[[i]])
    cnames[cnames %in% c("mu", "sigma", "x_pred")] <- c("mu_diff", "sigma_diff", "diff_pred")
    colnames(mcmc_samples[[i]]) <- cnames
  }
  mcmc_samples
}


bfa.t.test.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
                               mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95, ...) {
  
  if(var.equal) {
    var.equal <- FALSE
    warning("To assume equal variance of 'x' and 'y' is not supported. Continuing by estimating the variance of 'x' and 'y' separately.")
  }
  
  ### Original (but slighly modified) code from t.test.default ###
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                                 conf.level < 0 || conf.level > 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  
  # removing incomplete cases and preparing the data vectors (x & y)
  if (!is.null(y)) {
    xname <- deparse(substitute(x))
    yname <- deparse(substitute(y))
    dname <- paste(xname, "and", yname)
    if (paired) 
      xok <- yok <- complete.cases(x, y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else {
    xname <- deparse(substitute(x))
    dname <- xname
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
    mcmc_samples <- jags_paired_t_test(x, y)
    stats <- mcmc_stats(mcmc_samples, cred_mass = conf.level, comp_val = mu)
    bfa_object <- list(x = x, y = y, pair_diff = y - x, comp = mu, cred_mass = conf.level,
                        x_name = xname, y_name = yname, data_name = dname,
                        mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- "bfa_paired_t_test"
  
  } else if(is.null(y)) {
    mcmc_samples <- jags_one_sample_t_test(x, comp_mu = mu)
    stats <- mcmc_stats(mcmc_samples, cred_mass = conf.level, comp_val = mu)
    bfa_object <- list(x = x, comp = mu, cred_mass = conf.level, x_name = xname, 
                       data_name = dname, mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- "bfa_one_sample_t_test"
    
  } else { # is two sample t.test
    mcmc_samples <- jags_two_sample_t_test(x, y)
    stats <- mcmc_stats(mcmc_samples, cred_mass = conf.level, comp_val = mu)
    bfa_object <- list(x = x, y = y, comp = mu, cred_mass = conf.level,
                       x_name = xname, y_name = yname, data_name = dname,
                       mcmc_samples = mcmc_samples, stats = stats)
    class(bfa_object) <- "bfa_two_sample_t_test"
  }
  bfa_object
}

bfa.t.test.formula <- function(formula, data, subset, na.action, ...) {
  
  ### Original code from t.test.formula ###
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), "term.labels")) != 1L)) 
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L) 
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  
  ### Own code starts here ###
  bfa_object <- do.call("bfa.t.test", c(DATA, list(...)))
  bfa_object$data_name <- DNAME
  if (length(levels(g)) == 2L) {
    bfa_object$x_name <- paste("group", levels(g)[1])
    bfa_object$y_name <- paste("group", levels(g)[2])
  }
  bfa_object
  
}

### One sample t-test S3 methods ###

print.bfa_one_sample_t_test <- function(x) {
  
  s <- round(x$stats, 3)
  
  cat("\n")
  cat("\tBayesian estimation supersedes the t test (BEST)\n")
  cat("\n")
  cat("data: ", x$data_name, "\n", sep="")
  cat("\n")
  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep="")
  cat("mean of ",  x$x_name, ": ", s["mu", "mean"], " [", s["mu", "HDIlo"], ", ", s["mu", "HDIup"] , "]\n",sep="")
  cat("sd of ",  x$x_name, ": ", s["sigma", "mean"], " [", s["sigma", "HDIlo"], ", ", s["sigma", "HDIup"] , "]\n",sep="")
  
  cat("\n")
  cat("The mean is more than", s["mu","comp"] , "by a probability of", s["mu","%>comp"], "\n")
  cat("and less than", s["mu", "comp"] , "by a probability of", s["mu", "%<comp"], "\n")
  cat("\n")
}

print_bfa_one_sample_t_test_params <- function(x) {
  cat("  Model parameters and generated quantities\n")
  cat("mu: The mean of", x$data_name, "\n")
  cat("sigma: The standard deviation of", x$data_name,"\n")
  cat("nu: The degrees-of-freedom for the t distribution fitted to",x$data_name , "\n")
  cat("eff_size: The effect size calculated as (mu - ", x$comp ,") / sigma\n", sep="")
  cat("x_pred: Predicted distribution for a new datapoint generated as",x$data_name , "\n")
}

summary.bfa_one_sample_t_test <- function(x) {
  s <- round(x$stats, 3)
  
  cat("  Data\n")
  cat(x$data_name, ", n = ", length(x$x), "\n", sep="")
  cat("\n")
  
  print_bfa_one_sample_t_test_params(x)
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
}

plot.bfa_one_sample_t_test <- function(x) {
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
           xlab=bquote(mu) , main=paste("Mean") , col="skyblue" )
  
  # distribution of sigma:
  plotPost(sigma, cex.lab = 1.75, xlab=bquote(sigma), main=paste("Std. Dev."),  
           cred_mass= x$cred_mass, col="skyblue" , show_mode=TRUE )

  # effect size:
  plotPost(samples_mat[, "eff_size"] , comp_val=0 , xlab=bquote( (mu-.(x$comp)) / sigma ),
           cred_mass= x$cred_mass, show_mode=TRUE , cex.lab=1.75 , main="Effect Size" , col="skyblue" )
  par(old_par)
  invisible(NULL)
}

diagnostics.bfa_one_sample_t_test <- function(x) {
    
  print_mcmc_info(x$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(x$stats, 3))
  cat("\n")
  print_bfa_one_sample_t_test_params(x)
  cat("\n")
  
  old_par <- par( mar=c(3.5,2.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  plot(x$mcmc_samples)
  par(old_par)
  
}

model_diagram.bfa_one_sample_t_test <- function(x) {
  print(jags_binom_test)
}

model_code.bfa_one_sample_t_test <- function(x) {
  print(jags_binom_test)
}

### Two sample t-test S3 methods ###

print.bfa_two_sample_t_test <- function(x) {
  s <- round(x$stats, 3)
  
  cat("\n")
  cat("\tBayesian estimation supersedes the t test (BEST)\n")
  cat("\n")
  cat("data: ", x$data_name, "\n", sep="")
  cat("\n")
  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep="")
  cat("mean of ",  x$x_name, ": ", s["mu_x", "mean"], " [", s["mu_x", "HDIlo"], ", ", s["mu_x", "HDIup"] , "]\n",sep="")
  cat("mean of ",  x$y_name, ": ", s["mu_y", "mean"], " [", s["mu_y", "HDIlo"], ", ", s["mu_y", "HDIup"] , "]\n",sep="")
  cat("difference of the means: ", s["mu_diff", "mean"], " [", s["mu_diff", "HDIlo"], ", ", s["mu_diff", "HDIup"] , "]\n",sep="")
  cat("sd of ",  x$x_name, ": ", s["sigma_x", "mean"], " [", s["sigma_x", "HDIlo"], ", ", s["sigma_x", "HDIup"] , "]\n",sep="")
  cat("sd of ",  x$y_name, ": ", s["sigma_y", "mean"], " [", s["sigma_y", "HDIlo"], ", ", s["sigma_y", "HDIup"] , "]\n",sep="")
  
  cat("\n")
  cat("The difference of the means is more than", s["mu_diff","comp"] , "by a probability of", s["mu_diff","%>comp"], "\n")
  cat("and less than", s["mu_diff", "comp"] , "by a probability of", s["mu_diff", "%<comp"], "\n")
  cat("\n")
}
 
print_bfa_two_sample_t_test_params <- function(x) {
  cat("  Model parameters and generated quantities\n")
  cat("mu_x: The mean of", x$x_name, "\n")
  cat("sigma_x: The standard deviation of", x$x_name,"\n")
  cat("mu_y: The mean of", x$y_name, "\n")
  cat("sigma_y: The standard deviation of", x$y_name,"\n")
  cat("mu_diff: the difference in means (mu_x - mu_y)\n")
  cat("sigma_diff: the difference in SD (sigma_x - sigma_y)\n")
  cat("nu: The degrees-of-freedom for the t distribution\n")
  cat("  fitted to",x$data_name , "\n")
  cat("eff_size: The effect size calculated as \n", sep="")
  cat("  (mu_x - mu_y) / sqrt((sigma_x^2 + sigma_y^2) / 2)\n", sep="")
  cat("x_pred: Predicted distribution for a new datapoint\n")
  cat("  generated as",x$x_name , "\n")
  cat("y_pred: Predicted distribution for a new datapoint\n")
  cat("  generated as",x$y_name , "\n")
}

summary.bfa_two_sample_t_test <- function(x) {
  s <- round(x$stats, 3)
  
  cat("  Data\n")
  cat(x$x_name, ", n = ", length(x$x), "\n", sep="")
  cat(x$y_name, ", n = ", length(x$y), "\n", sep="")
  cat("\n")
  
  print_bfa_two_sample_t_test_params(x)
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
}

plot.bfa_two_sample_t_test <- function(x) {
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
  plotPost( mu_x ,  xlim=xlim , cex.lab = 1.75 , cred_mass= x$cred_mass,
            xlab=bquote(mu[x]) , main=paste(x$x_name,"Mean") , col="skyblue" )
  plotPost( mu_y ,  xlim=xlim , cex.lab = 1.75 ,  cred_mass= x$cred_mass,
            xlab=bquote(mu[y]) , main=paste(x$y_name,"Mean") , col="skyblue" )
  plotPost( samples_mat[,"mu_diff"] , comp_val= x$comp , cred_mass= x$cred_mass,
            xlab=bquote(mu[x] - mu[y]) , cex.lab = 1.75 , 
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
            show_mode=TRUE , cex.lab=1.0 , main="Effect Size" , col="skyblue" )
  
  par(old_par)
}

diagnostics.bfa_two_sample_t_test <- function(x) {
  print_mcmc_info(x$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(x$stats, 3))
  cat("\n")
  print_bfa_two_sample_t_test_params(x)
  cat("\n")
  
  old_par <- par( mar=c(3.5,2.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  plot(x$mcmc_samples)
  par(old_par)
  
}

model_diagram.bfa_two_sample_t_test <- function(bfa_result) {
  print(jags_two_sample_t_test)
}

model_code.bfa_two_sample_t_test <- function(bfa_result) {
  print(jags_two_sample_t_test)
}

### Paired samples t-test S3 methods ###

print.bfa_paired_t_test <- function(x) {
  s <- round(x$stats, 3)
  # Todo: ändra för att passa paired test.
  cat("\n")
  cat("\tBayesian estimation supersedes the t test (BEST)\n")
  cat("\n")
  cat("data: ", x$data_name, "\n", sep="")
  cat("\n")
  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep="")
  cat("mean paired difference: ", s["mu_diff", "mean"], " [", s["mu_diff", "HDIlo"], ", ", s["mu_diff", "HDIup"] , "]\n",sep="")
  cat("sd of the paired differences: ", s["sigma_diff", "mean"], " [", s["sigma_diff", "HDIlo"], ", ", s["sigma_diff", "HDIup"] , "]\n",sep="")
  
  cat("\n")
  cat("The mean difference is more than", s["mu_diff","comp"] , "by a probability of", s["mu_diff","%>comp"], "\n")
  cat("and less than", s["mu_diff", "comp"] , "by a probability of", s["mu_diff", "%<comp"], "\n")
  cat("\n")
  
  #   Paired t-test
  #   
  #   data:  runif(10) and rnorm(10)
  #   t = 1.9, df = 9, p-value = 0.0896
  #   alternative hypothesis: true difference in means is not equal to 0
  #   95 percent confidence interval:
  #     -0.106  1.230
  #   sample estimates:
  #     mean of the differences 
  #   0.562 
  
}

print_bfa_paired_t_test_params <- function(x) {

  cat("  Model parameters and generated quantities\n")
  cat("mu_diff: The mean pairwise difference between", x$x_name, "and", x$y_name, "\n")
  cat("sigma_diff: The standard deviation of the pairwise\n  difference between", x$x_name, "and", x$y_name,"\n")
  cat("nu: The degrees-of-freedom for the t distribution fitted to", x$data_name , "\n")
  cat("eff_size: The effect size calculated as (mu_diff - ", x$comp ,") / sigma_diff\n", sep="")
  cat("diff_pred: Predicted distribution for a new datapoint generated\n  as the pairwise difference between", x$x_name, "and", x$y_name,"\n")
}

summary.bfa_paired_t_test <- function(x) {
  s <- round(x$stats, 3)
  
  cat("  Data\n")
  cat(x$x_name, ", n = ", length(x$x), "\n", sep="")
  cat(x$y_name, ", n = ", length(x$y), "\n", sep="")
  cat("\n")
  
  print_bfa_paired_t_test_params(x)
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
}

plot.bfa_paired_t_test <- function(x) {
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
           xlab=bquote(mu[diff]) , main=paste("Mean difference") , col="skyblue" )
  
  # distribution of sigma_diff:
  plotPost(sigma_diff, cex.lab = 1.75, xlab=bquote(sigma[diff]), main=paste("Std. Dev. of difference"),  
           cred_mass= x$cred_mass, col="skyblue" , show_mode=TRUE )
  
  # effect size:
  plotPost(samples_mat[, "eff_size"] , comp_val=0 , xlab=bquote( (mu[diff] -.(x$comp)) / sigma[diff] ),
           cred_mass= x$cred_mass, show_mode=TRUE , cex.lab=1.75 , main="Effect Size" , col="skyblue" )
  par(old_par)
  invisible(NULL)
}


diagnostics.bfa_paired_t_test <- function(x) {
  print_mcmc_info(x$mcmc_samples)  
  cat("\n")
  print_diagnostics_measures(round(x$stats, 3))
  cat("\n")
  print_bfa_paired_t_test_params(x)
  cat("\n")
  
  old_par <- par( mar=c(3.5,2.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  plot(x$mcmc_samples)
  par(old_par)
}

model_diagram.bfa_paired_t_test <- function(x) {
  print(jags_binom_test)
}

model_code.bfa_paired_t_test <- function(x) {
  print(jags_binom_test)
}