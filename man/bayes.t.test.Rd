\name{bayes.t.test}
\alias{bayes.t.test}
\alias{bayes.t.test.default}
\alias{bayes.t.test.formula}
\title{Title title title}
\usage{
bayes.t.test(x, ...)

\method{bayes.t.test}{default}(x, y = NULL, alternative = c("two.sided",
  "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE,
  conf.level = 0.95, n.iter = 15000, ...)

\method{bayes.t.test}{formula}(formula, data, subset, na.action, ...)
}
\arguments{
  \item{x}{}

  \item{...}{}

  \item{y}{}

  \item{alternative}{}

  \item{mu}{}

  \item{paired}{}

  \item{var.equal}{}

  \item{conf.level}{}

  \item{n.iter}{}

  \item{formula}{}

  \item{data}{}

  \item{subset}{}

  \item{na.action}{}
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

