#' Title title title
#' 
#' Descritions description description
#' 
#' Details details details
#' 
#' @param x 
#' @param ... 
#' @param y 
#' @param alternative 
#' @param method 
#' @param exact 
#' @param cred.mass
#' @param continuity 
#' @param n.iter
#' @param progress.bar The type of progress bar. Possible values are "text",
#'   "gui", and "none".
#' @param formula 
#' @param data 
#' @param subset 
#' @param na.action
#' @param conf.level 
#' 
#' @return A list of class \code{bayes_cor_test} that contains information about
#'   the analysis. It can be further inspected using the functions
#'   \code{summary}, \code{plot}, \code{\link{diagnostics}} and 
#'   \code{\link{model.code}}.
#' @export
#' @rdname bayes.cor.test
bayes.cor.test <- function(x, ...) {
  UseMethod("bayes.cor.test")
}

cor_model_string <- "model {
  for(i in 1:n) {
    xy[i,1:2] ~ dmt(mu[], prec[ , ], nu) 
  }

  ## priors for elements of the precision matrix
  prec[1:2,1:2] <- inverse(cov[,])
  
  cov[1,1] <- sigma[1] * sigma[1]
  cov[1,2] <- sigma[1] * sigma[2] * rho
  cov[2,1] <- sigma[1] * sigma[2] * rho
  cov[2,2] <- sigma[2] * sigma[2]
  
  sigma[1] ~ dunif(0, 1000) 
  sigma[2] ~ dunif(0, 1000)
  tau[1] <- 1 / pow(sigma[1], 2)
  tau[2] <- 1 / pow(sigma[2], 2)
  rho ~ dunif(-1, 1)

  mu[1] ~ dnorm(0, 0.0001)
  mu[2] ~ dnorm(0, 0.0001)

  nu <- nuMinusOne+1
  nuMinusOne ~ dexp(1/29)
}"

jags_cor_test <- function(x, y, n.adapt= 500, n.chains=3, n.update = 100, n.iter=5000, thin=1, progress.bar="text") {
  data_list = list(xy = cbind(x, y), n = length(x))
  # Use robust estimates of the parameters as initial values
  inits_list = list(mu=c(mean(x, trim=0.2), mean(y, trim=0.2)), rho=cor(x, y, method="spearman"), 
                    sigma = c(mad(x), mad(y)), nuMinusOne = 5)
  mcmc_samples <- run_jags(cor_model_string, data = data_list, inits = inits_list, 
                           params = c("rho", "mu", "sigma", "nu"), n.chains = n.chains, n.adapt = n.adapt,
                           n.update = n.update, n.iter = n.iter, thin = thin, progress.bar=progress.bar)
  mcmc_samples
}

#' @method bayes.cor.test default
#' @export
#' @rdname bayes.cor.test
bayes.cor.test.default <- function (x, y, alternative = c("two.sided", "less", "greater"), 
                                  method = c("pearson", "kendall", "spearman"), exact = NULL, 
                                  cred.mass = 0.95, continuity = FALSE, n.iter = 15000, progress.bar="text",
                                  conf.level, ...) 
{
  
  if(! missing(conf.level)) {
    cred.mass <- conf.level
  }
  
  ### BEGIN code from cor.test.default ###
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  x_name <- deparse(substitute(x))
  y_name <- deparse(substitute(y))
  data_name <- paste(x_name, "and", y_name)
  if (length(x) != length(y)) 
    stop("'x' and 'y' must have the same length")
  if (!is.numeric(x)) 
    stop("'x' must be a numeric vector")
  if (!is.numeric(y)) 
    stop("'y' must be a numeric vector")
  # removes uncomplete pairs, this shouldn't be neccessary if JAGS could handle missing data in dmvt
  OK <- complete.cases(x, y)
  x <- x[OK]
  y <- y[OK]
  n <- length(x)
  if (n < 3L) 
    stop("not enough observations. Need at least three complete observation.")
  ### END code from cor.test.default
  if (method == "kendall" || method == "spearman") {
    stop("no non-parametric correlation comparable to Kendall's tau or Spearman's rho has been implemented yet.")
  }
  mcmc_samples <- jags_cor_test(x, y, n.chains=3, n.iter=ceiling(n.iter / 3), progress.bar=progress.bar)
  stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = 0)
  bfa_result <- list(x = x, y = y, cred_mass = cred.mass, x_name = x_name, y_name = y_name, 
                     data_name = data_name, mcmc_samples = mcmc_samples, stats = stats)
  class(bfa_result) <- c("bayes_cor_test", "bayesian_first_aid")
  bfa_result
  
}

#' @method bayes.cor.test formula
#' @export
#' @rdname bayes.cor.test
bayes.cor.test.formula <- function (formula, data, subset, na.action, ...) 
{
  ### BEGIN code from cor.test.formula ###
  if (missing(formula) || !inherits(formula, "formula") || 
        length(formula) != 2L) 
    stop("'formula' missing or invalid")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m, environment(formula))
  if (length(mf) != 2L) 
    stop("invalid formula")
  DNAME <- paste(names(mf), collapse = " and ")
  names(mf) <- c("x", "y")
  ### END code from cor.test.formula ###
  
  bfa_result <- do.call("bayes.cor.test.default", c(mf, list(...)))
  bfa_result$data_name <- DNAME
  bfa_result
}

### Cor test S3 methods ###

#' @export
print.bayes_cor_test <- function(x, ...) {
  cat("\n --- Bayesian first aid cor test ---\n\n")
  print(summary(x$mcmc_samples))
}

#' @export
summary.bayes_cor_test <- function(object, ...) {
  cat("\nSummary\n")
  print(object)
}

#' @export
plot.bayes_cor_test <- function(x, ...) {
  stats <- x$stats
  mcmc_samples <- x$mcmc_samples
  samples_mat <- as.matrix(mcmc_samples)
  
  rho <- samples_mat[, "rho"]
    
  old_par <- par(no.readonly = TRUE)
  
  zones <- matrix(c(1,1,1,
                    0,6,0,
                    3,2,5,
                    0,4,0), ncol = 3, byrow = TRUE)
  layout(zones, widths=c(0.5,4,1), heights = c(6,3,10,1.5))
  
  ### fig 1, the posterior of rho ###
  par(mar = c(2,2,2,2))
  plotPost(rho, comp_val=0, cred_mass=0.95, xlim=c(-1, 1), xlab="", main=expression(Correlation ~ (rho)))
  
  
  ### fig 2, the scatterplot ###
  
  # Sampling from the posterior predictive distribution 
  xy_rep <- do.call(rbind, lapply(sample(1:nrow(samples_mat), 1000), function(i) {
    sigma1 <- samples_mat[i, "sigma[1]"]
    sigma2 <- samples_mat[i, "sigma[2]"]
    rho <- samples_mat[i, "rho"]
    cov_mat <- cbind(c(sigma1^2, sigma1 * sigma2 * rho),
                     c(sigma1 * sigma2 * rho, sigma2^2))
    rmt(100, samples_mat[i, c("mu[1]", "mu[2]")], cov_mat, samples_mat[i, "nu"])
  }))
  
  x_rep <- xy_rep[,1]
  y_rep <- xy_rep[,2]
  
  # Calculating the 2d density of the posterior predictive distribution x_rep and y_rep
  dens_limits <- c(median(x_rep) - IQR(x_rep) * 5, median(x_rep) + IQR(x_rep) * 5,
                   median(y_rep) - IQR(y_rep) * 5, median(y_rep) + IQR(y_rep) * 5)
  bandwidth <- c(IQR(x_rep), IQR(y_rep)) / 1.349 * 0.5
  dens_2d <- kde2d(x_rep, y_rep, n = 40, lims=dens_limits, h=bandwidth)
  sorted_z <- sort(dens_2d$z, decreasing=TRUE)
  post_95_limit <- sorted_z[which.min(abs(cumsum(sorted_z / sum(sorted_z)) - 0.95))]
  post_50_limit <- sorted_z[which.min(abs(cumsum(sorted_z / sum(sorted_z)) - 0.5))]
  
  # These messy lines calculates limits of the plot that makes sure both all the data
  # and the density estimate is visible. Also centers the plot on the median of the data.
  plot_xlim <- c(median(x$x) - max( abs(x$x - median(x$x))), median(x$x) + max( abs(x$x - median(x$x))))
  plot_xlim[1] <- min(plot_xlim[1], dens_2d$x[apply(dens_2d$z > post_95_limit, 2, any)] - diff(dens_2d$x[1:2]) / 2)
  plot_xlim[2] <- max(plot_xlim[2], dens_2d$x[apply(dens_2d$z > post_95_limit, 2, any)] + diff(dens_2d$x[1:2])/ 2)
  plot_ylim <- c(median(x$y) - max( abs(x$y - median(x$y))), median(x$y) + max( abs(x$y - median(x$y))))
  plot_ylim[1] <- min(plot_ylim[1], dens_2d$y[apply(dens_2d$z > post_95_limit, 1, any)] - diff(dens_2d$y[1:2])/ 2)
  plot_ylim[2] <- max(plot_ylim[2], dens_2d$y[apply(dens_2d$z > post_95_limit, 1, any)] + diff(dens_2d$y[1:2])/ 2)
  
  par(mar = c(2,2,0,0), xaxt="s", yaxt="s", bty="o")
  plot(x$x, x$y, col="white", xlim=plot_xlim, ylim=plot_ylim)
  .filled.contour(dens_2d$x, dens_2d$y, dens_2d$z, levels=c(0, post_95_limit, post_50_limit, max(sorted_z)), col=c(rgb(1,1,1, 0), "#bddeeb", "#87ceeb"))
  points(x$x, x$y, pch=1, col=c(rgb(0,0,0, 1)))
  
  # Save the limits for use with the marginal plots
  scatter_lim <- par("usr")
  
  ### The titles of the scatter plot ###
  
  par(xaxt="n", yaxt="n",bty="n",  mar = c(.3,2,.3,0) +.05)
  # fig 3
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))
  text(0,0, x$y_name, cex=1.5, srt=90)
  # fig 4
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))
  text(0,0, x$x_name, cex=1.5)
  
  ### fig 5, the y histogram ###
  par(mar = c(2,0,0,1))
  
  i <- sample(1:nrow(samples_mat), 20)
  y_mu <-  samples_mat[i, "mu[2]"]
  y_sigma <-  samples_mat[i, "sigma[2]"]
  y_nu <- samples_mat[i, "nu"]
  
  hist_with_t_curves(x$y, stats["mu[2]", "mean"], stats["sigma[2]", "mean"], y_mu, y_sigma, y_nu, 
                     axes = FALSE, horiz=TRUE, plot_n=FALSE, x_lim=scatter_lim[3:4], axs="i")

  ### fig 6, The x histogram ###
  par(mar = c(0,2,1,0))
  
  i <- sample(1:nrow(samples_mat), 20)
  x_mu <-  samples_mat[i, "mu[1]"]
  x_sigma <-  samples_mat[i, "sigma[1]"]
  x_nu <- samples_mat[i, "nu"]
  
  hist_with_t_curves(x$x,stats["mu[1]", "mean"], stats["sigma[1]", "mean"], x_mu, x_sigma, y_nu,
                     horiz=FALSE, plot_n=TRUE, axes = FALSE, x_lim=scatter_lim[1:2], axs="i")
  
  par(old_par)
  invisible(NULL)
}

#' @export
diagnostics.bayes_cor_test <- function(fit) {
  cat("Not implemented\n")
  plot(fit$mcmc_samples)
}

#' @export
model.code.bayes_cor_test <- function(fit) {
  print(jags_cor_test)
}
