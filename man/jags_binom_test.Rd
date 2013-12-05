\name{jags_binom_test}
\alias{jags_binom_test}
\title{Generate MCMC samples from the Bayesian First Aid binomial model using JAGS.}
\usage{
  jags_binom_test(x, n, n.adapt = 1000, n.chains = 3,
    n.iter = 1000)
}
\arguments{
  \item{x}{The number of successes}

  \item{n}{The number of trials}
}
\value{
  An object of type \code{mcmc.list} (defined by the
  \code{coda} package) that contains the MCMC samples from
  the model.
}
\description{
  Descritions description description
}
\details{
  Details details details
}

