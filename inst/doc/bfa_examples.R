
## ----echo=FALSE----------------------------------------------------------
set.seed(123)


## ------------------------------------------------------------------------
binom.test(x=16, n=20)


## ------------------------------------------------------------------------
bfa.binom.test(x=16, n=20)


## ------------------------------------------------------------------------
binom.test(x=16, n=20, p=0.75, conf.level=0.8)


## ------------------------------------------------------------------------
bfa.binom.test(x=16, n=20, p=0.75, conf.level=0.8)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bfa.binom.test(x=16, n=20))
plot(fit)
summary(fit)
diagnostics(fit)



## ------------------------------------------------------------------------
x <- rnorm(n=30, mean=2, sd=5)
t.test(x)


## ------------------------------------------------------------------------
bfa.t.test(x)


## ------------------------------------------------------------------------
y <- rnorm(n=30, mean=10, sd=5)
bfa.t.test(x, y, mu=1, conf.level=0.5)


## ------------------------------------------------------------------------
# Creating "paired" data
d <- data.frame(x = rnorm(n=30, mean=10, sd=5), group = rep(c("A", "B"), 15))
                
bfa.t.test(x ~ group, paired=T, data=d)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bfa.t.test(x, mu=1) )
plot(fit)
summary(fit)
diagnostics(fit)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bfa.t.test(x, y, conf.level=0.8) )
plot(fit)
summary(fit)
diagnostics(fit)


## ----fig.width=6, fig.height=5, dpi=96-----------------------------------
(fit <- bfa.t.test(x ~ group, paired=T, data=d))
plot(fit)
summary(fit)
diagnostics(fit)


