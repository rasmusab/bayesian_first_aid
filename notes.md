Bayesian First Aid
========================================================

Tests to include:

  * t.test 
    - Using BEST
  * cor.test
    - Maybe fit multivariate t-dist. Or something non-parametric?
  * binom.test
    - Just fit regular Bernouli distribution
  * chisq.test & fisher.test
    - Pointers to how to do this are here: http://andrewgelman.com/2009/10/13/what_is_the_bay/ and http://lingpipe-blog.com/2009/10/13/bayesian-counterpart-to-fisher-exact-test-on-contingency-tables/. These are just for the 2x2 case and I'm not sure how to generalize it. chisq.test can be used for many  things. Maybe start by implementing it just for the common 2x2 case.
  * aov
  * lm
  * var.test
    - Same as t.test but change the output to reflect the comparison between the variances.
  * poisson.test
    - Robustify the poisson test by using a negative-binomial distribution but construct it by pasting together a Poisson with a gamma distributed mean as on p. 253 in the Bugs book.
  

