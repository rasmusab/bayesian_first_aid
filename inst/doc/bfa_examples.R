
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
bayes.t.test(x, progress.bar="none")


## ------------------------------------------------------------------------
y <- rnorm(n=30, mean=10, sd=5)
bayes.t.test(x, y, mu=1, conf.level=0.5, progress.bar="none")


## ------------------------------------------------------------------------
# Creating "paired" data
d <- data.frame(x = rnorm(n=30, mean=10, sd=5), group = rep(c("A", "B"), 15))
                
bayes.t.test(x ~ group, paired=T, data=d, progress.bar="none")


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bayes.t.test(x, mu=1, progress.bar="none") )
plot(fit)
summary(fit)
diagnostics(fit)
model.code(fit)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bayes.t.test(x, y, conf.level=0.8, progress.bar="none") )
plot(fit)
summary(fit)
diagnostics(fit)
model.code(fit)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bayes.t.test(x ~ group, paired=T, data=d, progress.bar="none"))
plot(fit)
summary(fit)
diagnostics(fit)
model.code(fit)


## ------------------------------------------------------------------------
poisson.test(15, 3, r=3)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bayes.poisson.test(15, 3, r=3))
plot(fit)
summary(fit)
diagnostics(fit)
model.code(fit)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bayes.poisson.test(c(15, 20), c(3, 3), r=1))
plot(fit)
summary(fit)
diagnostics(fit)
model.code(fit)


## ------------------------------------------------------------------------
a = c(-1.5, -1.6, -0.5, -1.5, 0.7, 2.1, -1.3, 0.8)
b = c(-1.9, -1.8, -1.5, -0.6, 0.4, 0.7, -1.7, -1.8)
cor.test(a, b, conf.level=0.8)


## ------------------------------------------------------------------------
(fit <- bayes.cor.test(a, b, conf.level=0.8, progress.bar="none") )
plot(fit)
summary(fit)
diagnostics(fit)
model.code(fit)


