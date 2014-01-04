
## ----echo=FALSE----------------------------------------------------------
set.seed(123)


## ------------------------------------------------------------------------
library(BayesianFirstAid)
binom.test(x=16, n=20)


## ------------------------------------------------------------------------
bayes.binom.test(x=16, n=20)


## ------------------------------------------------------------------------
binom.test(x=16, n=20, p=0.75, conf.level=0.8)


## ------------------------------------------------------------------------
bayes.binom.test(x=16, n=20, p=0.75, conf.level=0.8)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bayes.binom.test(x=16, n=20))
plot(fit)
summary(fit)
diagnostics(fit)


## ------------------------------------------------------------------------
model.code(fit)


## ------------------------------------------------------------------------
x <- rnorm(n=30, mean=2, sd=5)
t.test(x)


## ------------------------------------------------------------------------
bayes.t.test(x)


## ------------------------------------------------------------------------
y <- rnorm(n=30, mean=10, sd=5)
bayes.t.test(x, y, mu=1, conf.level=0.5)


## ------------------------------------------------------------------------
# Creating "paired" data
d <- data.frame(x = rnorm(n=30, mean=10, sd=5), group = rep(c("A", "B"), 15))
                
bayes.t.test(x ~ group, paired=T, data=d)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bayes.t.test(x, mu=1) )
plot(fit)
summary(fit)
diagnostics(fit)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bayes.t.test(x, y, conf.level=0.8) )
plot(fit)
summary(fit)
diagnostics(fit)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bayes.t.test(x ~ group, paired=T, data=d))
plot(fit)
summary(fit)
diagnostics(fit)


