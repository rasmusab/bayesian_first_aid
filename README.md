Bayesian First Aid
======================

Bayesian First Aid is an R package that implements Bayesian alternatives to the most commonly used statistical tests. The package aims at being as easy as possible to pick up and use, especially if you are already used to the classical `.test` functions. The main gimmick is that the Bayesian alternatives will have the same calling semantics as the corresponding classical test functions save for the addition of `bayes.` to the beginning of the function name. That is, going from a classical binomial test:

``` S
binom.test(x=7, n=10, p=0.33)
```

```
    Exact binomial test

data:  7 and 10
number of successes = 7, number of trials = 10, p-value = 0.01855
alternative hypothesis: true probability of success is not equal to 0.33
95 percent confidence interval:
 0.3475 0.9333
sample estimates:
probability of success 
                   0.7 
```

... to the Bayesian estimation version is as easy as prepending `bayes.` to the function call:

``` S
bayes.binom.test(x=7, n=10, p=0.33)
```
```
    Bayesian first aid binomial test

data: 7 and 10
number of successes = 7, number of trials = 10
Estimated relative frequency of success:
  0.669 
95% credible interval:
  0.42 0.91 
The relative frequency of success is more than 0.33 by a probability of 0.993 
and less than 0.33 by a probability of 0.007 
```

Hence the name of the package: Bayesian First Aid.

<center>
<img src="http://i.imgur.com/gQAiFkD.png"/>
</center>

Read more about the rational behind Bayesian First aid in [the announcement on my blog](http://sumsar.net/blog/2014/01/bayesian-first-aid/). If you use Bayesian First Aid and want to cite it, please use the following citation:

Bååth, R., (2014) Bayesian First Aid: A Package that Implements Bayesian Alternatives to the Classical \*.test Functions in R. In the proceedings of *UseR! 2014 - the International R User Conference*. [pdf](http://sumsar.net/files/academia/baath_user14_abstract.pdf)

```
@inproceedings{baath2014bayesian,
  title={Bayesian First Aid: A Package that Implements Bayesian Alternatives to the Classical *.test Functions in R},
  author={B{\aa}{\aa}th, Rasmus},
  booktitle={UseR! 2014 - the International R User Conference},
  year={2014}
}
```
*Be aware that Bayesian First Aid is under heavy development and is still in a beta stage.*

Installation
--------------------

This package requires a working instalation of [JAGS](http://mcmc-jags.sourceforge.net/). To install the package from github you also need the [devtools package](http://cran.r-project.org/web/packages/devtools/index.html):

```
install.packages("devtools")
```

To install Bayesian First Aid then run:

```
library(devtools)
install_github("rasmusab/bayesian_first_aid")
```

The installation will take a couple of minutes due to the vignettes taking quite a long time to build, so please be patient :)

The Implemented Models
--------------------------------

Currently Bayesian First Aid implements five alternative models to classical tests. Each of these alternative models are described in detail on my blog (see links below):

* Binomial Test: `bayes.binom.test(x,n)` . [Code example and details](http://sumsar.net/blog/2014/01/bayesian-first-aid-binomial-test/)
* One Sample and Paired Samples t-test: `bayes.t.test(x)` . [Code example and details](http://sumsar.net/blog/2014/02/bayesian-first-aid-one-sample-t-test/)
* Two Sample t-test: `bayes.t.test(x, y)` . [Code example and details](http://sumsar.net/blog/2014/02/bayesian-first-aid-two-sample-t-test/)
* Pearson Correlation Test: `bayes.cor.test(x,y)` . [Code example and details](http://sumsar.net/blog/2014/03/bayesian-first-aid-pearson-correlation-test/)
* Test of Proportions: `prop.test` . [Code example and details](http://sumsar.net/blog/2014/06/bayesian-first-aid-prop-test/)

It should be easy to start tinkering with the models underlying Bayesian First Aid. The generic function `model.code` takes a Bayesian First Aid object and prints out the underlying model which is ready to be copy-n-pasted into an R script and tinkered with from there. In addition to `model.code` all Bayesian First Aid objects will be `plot`able and `summary`iseable. A call to `diagnostics` will show some MCMC diagnostics (even if this shouldn't be necessary to look at for the simple models). For an example of how this would work, see the example further down.

Except for the already implemented Bayesian alternatives "replacing" the following tests are currently on the todo list:

* `var.test`
* `poisson.test`
* `chisq.test`
* `oneway.test`
* `wilcox.test`

An Example of Bayesian First Aid in Action
-------------------------------------------

Just a quick look at how Bayesian First Aid works. The following simple problem is from a statistical methods course (*[source](http://www.elderlab.yorku.ca/~aaron/Stats2022/BinomialTest.htm)*):

> A college dormitory recently sponsored a taste comparison between two major soft drinks. Of the 64 students who participated, 39 selected brand A, and only 25 selected brand B.  Do these data indicate a significant preference? 

This problem is written in order to be tested by a binomial test, so let's do that:


``` S
binom.test(c(39, 25))
```

```

  Exact binomial test

data:  c(39, 25)
number of successes = 39, number of trials = 64, p-value = 0.1034
alternative hypothesis: true probability of success is not equal to 0.5
95 percent confidence interval:
 0.4793 0.7290
sample estimates:
probability of success 
                0.6094 
```


Bummer, seems like there is no statistically significant difference between the two brands. Now we're going to run the Bayesian First Aid alternative simply by prepending `bayes.` to our function call:


``` S
library(BayesianFirstAid)
bayes.binom.test(c(39, 25))
```

```

	Bayesian first aid binomial test

data: c(39, 25)
number of successes = 39, number of trials = 64
Estimated relative frequency of success:
  0.606 
95% credible interval:
  0.492 0.727 
The relative frequency of success is more than 0.5 by a probability of 0.959 
and less than 0.5 by a probability of 0.041 
```


Great, we just estimated the relative frequency θ of choosing brand A assuming the following model:

<img src="http://i.imgur.com/4YZxHDt.png"/>

In this simple example the estimates from `binom.test` and `bayes.binom.test` are close to identical. Still we get to know that `The relative frequency of success is more than 0.5 by a probability of 0.956` which indicates that brand A might be more popular. At the end of the day we are not that interested in *whether* brand A is more popular that brand B but rather *how much more* popular brand A is. This is easier to see if we look at the posterior distribution of θ: 


``` S
plot(bayes.binom.test(c(39, 25)))
```

![plot of chunk unnamed-chunk-3](http://i.imgur.com/GAoNH6s.png) 


So it seems that brand A might be more popular than brand B but not by too much as the posterior has most of its mass at a relative frequency between 0.5 and 0.7. The college cafeteria should probably keep both brands in stock if possible. If they need to choose one brand, pick A.

Perhaps you want to tinker with the model above, maybe you have some prior knowledge that you want to incorporate. By using the `model.code` function we get a nice printout of the model that we can copy-n-paste into an R script.


``` S
model.code(bayes.binom.test(c(39, 25)))
```

``` S
### Model code for the Bayesian First Aid alternative to the binomial test ###
require(rjags)

# Setting up the data
x <- 39 
n <- 64 

# The model string written in the JAGS language
model_string <-"model {
  x ~ dbinom(theta, n)
  theta ~ dbeta(1, 1)
  x_pred ~ dbinom(theta, n)
}"

# Running the model
model <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                    n.chains = 3, n.adapt=1000)
samples <- coda.samples(model, c("theta", "x_pred"), n.iter=5000)

#Inspecting the posterior
plot(samples)
summary(samples)  
```
