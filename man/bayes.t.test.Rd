\name{bayes.t.test}
\alias{bayes.t.test}
\alias{bayes.t.test.default}
\alias{bayes.t.test.formula}
\title{Title title title}
\usage{
bayes.t.test(x, ...)

\method{bayes.t.test}{default}(x, y = NULL, alternative = c("two.sided",
  "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE,
  cred.mass = 0.95, n.iter = 15000, conf.level, ...)

\method{bayes.t.test}{formula}(formula, data, subset, na.action, ...)
}
\arguments{
  \item{x}{a (non-empty) numeric vector of data values.}

  \item{y}{an optional (non-empty) numeric vector of data
  values.}

  \item{alternative}{ignored and is only retained in order
  to mantain compatibility with \code{\link{t.test}}.}

  \item{mu}{a fixed relative mean value to compare with the
  estimated mean (or the difference in means when
  performing a two sample BEST).}

  \item{paired}{a logical indicating whether you want to
  assume the data is paired and you want to run a one
  sample BEST of the pair difference.}

  \item{var.equal}{ignored and is only retained in order to
  mantain compatibility with \code{\link{t.test}}.}

  \item{cred.mass}{the amount of probability mass that will
  be contained in reported credible intervals. This
  argument fills a similar role as \code{conf.level} in
  \code{\link{t.test}}.}

  \item{n.iter}{The number of iterations to run the MCMC
  sampling.}

  \item{conf.level}{same as \code{cred.mass} and is only
  retained in order to mantain compatibility with
  \code{\link{binom.test}}.}

  \item{formula}{a formula of the form \code{lhs ~ rhs}
  where \code{lhs} is a numeric variable giving the data
  values and \code{rhs} a factor with two levels giving the
  corresponding groups.}

  \item{data}{an optional matrix or data frame (or similar:
  see \code{\link{model.frame}}) containing the variables
  in the formula formula.  By default the variables are
  taken from \code{environment(formula)}.}

  \item{subset}{an optional vector specifying a subset of
  observations to be used.}

  \item{na.action}{a function which indicates what should
  happen when the data contain \code{NA}s. Defaults to
  \code{getOption("na.action")}.}

  \item{...}{further arguments to be passed to or from
  methods.}
}
\value{
An object of type something...
}
\description{
Descritions description description
}
\details{
Details details details

\figure{bestmodel.png}{A graphical diagram of the BEST
model}

\deqn{y \sim \mathrm{Norm}(\mu, \sigma)}{y ~ Norm(mu,
sigma)}
}

