\name{as.data.frame.bayesian_first_aid}
\alias{as.data.frame.bayesian_first_aid}
\title{Returns the MCMC samples from a Bayesian First Aid fit as a data frame.}
\usage{
\method{as.data.frame}{bayesian_first_aid}(x, ...)
}
\arguments{
  \item{x}{The output from a Bayesian First Aid model.}
}
\value{
A data frame with one column per parameter.
}
\description{
Returns a dataframe with the MCMC samples from a Bayesian
First Aid Object. The samples does not include the burn in
phase and all chains are collapsed in the resulting data
frame.
}
\details{
In the case where there are many parameters with the same
name but different index (for example, mu[1], mu[2], mu[3],
etc.) these names will be rewritten by removing the
brackets (that is, mu[1] -> mu1, mu[2] -> mu2, sigma[1,2]
-> sigma12, etc.).
}
\examples{
fit <- bayes.t.test(rnorm(10), rnorm(10))
d <- as.data.frame(fit)
mean(d$mu_x)
}

