\name{bayes.binom.test}
\alias{bayes.binom.test}
\title{Bayesian First Aid Alternative to the Binomial Test}
\usage{
bayes.binom.test(x, n, comp.theta = 0.5, alternative = NULL,
  cred.mass = 0.95, n.iter = 15000, p, conf.level)
}
\arguments{
  \item{x}{number of successes, or a vector of length 2
  giving the numbers of successes and failures,
  respectively.}

  \item{n}{number of trials; ignored if x has length 2.}

  \item{comp.theta}{a fixed relative frequency of success
  to compare with the estimated relative frequency of
  success. This argument fills a similar role as \code{p}
  in \code{\link{binom.test}}.}

  \item{alternative}{ignored and is only retained in order
  to mantain compatibility with \code{\link{binom.test}}.}

  \item{cred.mass}{the amount of probability mass that will
  be contained in reported credible intervals. This
  argument fills a similar role as \code{conf.level} in
  \code{\link{binom.test}}.}

  \item{n.iter}{The number of iterations to run the MCMC
  sampling.}

  \item{p}{same as \code{comp.theta} and is only retained
  in order to mantain compatibility with
  \code{\link{binom.test}}.}

  \item{conf.level}{same as \code{cred.mass} and is only
  retained in order to mantain compatibility with
  \code{\link{binom.test}}.}
}
\value{
A list of class \code{bfa_binom_test} that contains
information about the analysis. It can be further inspected
using the functions \code{summary}, \code{plot},
\code{\link{diagnostics}} and \code{\link{model.code}}.
}
\description{
Descritions description description
}
\details{
\deqn{x \sim \mathrm{Binom}(\theta, n)}{x ~ Binom(\theta,
n)} \deqn{\theta \sim \mathrm{Beta}(1, 1)}{\theta ~ Beta(1,
1)}

\if{html}{\figure{binom_diagram.svg}{options: height=250}}
\if{latex}{binom_diagram.svg}
}
\examples{
bayes.binom.test(5, 10)

# Save the return value in order to inspect the model result further.
fit <- bayes.binom.test(x=67, n=100, cred.mass=0.8)
summary(fit)
plot(fit)
# MCMC diagnostics (should not be necessary for such a simple model)
diagnostics(fit)
# Print out the R code to run the model. This can be copy n' pasted into
# an R-script and further modified.
model.code(fit)
}

