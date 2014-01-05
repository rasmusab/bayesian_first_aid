*This is still in development and currently not fit for any purpose.*

The idea with this R package is to make "replacements" for the most commonly used tests in R such as `t.test`, `binom.test` and `cor.test`. These replacements will be based on Bayesian estimation and will in that sense neither be "null hypothesis" nor "tests". They will be replacements in that they will have similar assumptions as the original tests and will answer the the same type of question as one probably have when using the corresponding null hypothesis test.

The gimmick of the package is that the Bayesian versions have functions calls that are compatible with the classical tests' functions. That is going from a classical binomial test:

``` S
binom.test(x=7, n=10)
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

... to the Bayesian estimation version is as easy as prepending bayes. to the function call:

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
