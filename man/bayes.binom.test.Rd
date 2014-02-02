\name{bayes.binom.test}
\alias{bayes.binom.test}
\title{Bayesian First Aid Alternative to the Binomial Test}
\usage{
bayes.binom.test(x, n, comp.theta = 0.5, alternative = NULL,
  cred.mass = 0.95, n.iter = 15000, progress.bar = "none", p, conf.level)
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

  \item{progress.bar}{The type of progress bar. Possible
  values are "text", "gui", and "none".}

  \item{p}{same as \code{comp.theta} and is only retained
  in order to mantain compatibility with
  \code{\link{binom.test}}.}

  \item{conf.level}{same as \code{cred.mass} and is only
  retained in order to mantain compatibility with
  \code{\link{binom.test}}.}
}
\value{
A list of class \code{bayes_binom_test} that contains
information about the analysis. It can be further inspected
using the functions \code{summary}, \code{plot},
\code{\link{diagnostics}} and \code{\link{model.code}}.
}
\description{
\code{bayes.binom.test} estimates the relative frequency of
success using Bayesian estimation and is intended as a
replacement for \code{\link{binom.test}}.
}
\details{
Given data on the number of successes and failures
\code{bayes.binom.test} estimates \eqn{\theta}, the
relative frequency of success, assuming the following
model:

\deqn{x \sim \mathrm{Binom}(\theta, n)}{x ~
Binomial(\theta, n)} \deqn{\theta \sim \mathrm{Beta}(1,
1)}{\theta ~ Beta(1, 1)}

\if{html}{\figure{binom_diagram.svg}{options: height=250}}
\if{latex}{\figure{binom_diagram.pdf}}

Here the prior on \eqn{\theta} is a non-informative
\eqn{\mathrm{Beta}(1, 1)}{Beta(1, 1)} distribution which is
identical to a \eqn{\mathrm{Uniform}(0, 1)}{Uniform(0, 1)}
distribution. By \code{plot}ing and looking at a
\code{summary} of the object returned by
\code{bayes.binom.test} you can get more information about
the shape of the posterior and the posterior predictive
distribution. \code{\link{model.code}} prints out the
corresponding R code underlying \code{bayes.binom.test}
which can be copy-n-pasted into an R script and modified,
for example, changing the prior on \eqn{\theta}.
}
\examples{
# A college dormitory recently sponsored a taste comparison between
# two major soft drinks. Of the 64 students who participated, 39 selected
# brand A, and only 25 selected brand B. Example from
# http://www.elderlab.yorku.ca/~aaron/Stats2022/BinomialTest.htm

bayes.binom.test(x = 39, n = 64)

# Save the return value in order to inspect the model result further.
fit <- bayes.binom.test(x = 39, n = 64, cred.mass=0.8)
summary(fit)
plot(fit)

# MCMC diagnostics (should not be necessary for such a simple model)
diagnostics(fit)

# Print out the R code to run the model. This can be copy n' pasted into
# an R-script and further modified.
model.code(fit)
}

