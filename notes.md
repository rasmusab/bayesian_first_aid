Bayesian First Aid
========================================================

Tests to include:

  * t.test 
    - Using BEST
  * cor.test
    - Maybe fit multivariate t-dist. Or something non-parametric?
  * binom.test
    - Just fit regular Bernouli distribution
  * prop.test 
  * chisq.test & fisher.test
    - Pointers to how to do this are here: http://andrewgelman.com/2009/10/13/what_is_the_bay/ and http://lingpipe-blog.com/2009/10/13/bayesian-counterpart-to-fisher-exact-test-on-contingency-tables/. These are just for the 2x2 case and I'm not sure how to generalize it. chisq.test can be used for many  things. Maybe start by implementing it just for the common 2x2 case.
  * aov
  * oneway.test
  * lm
  * var.test
    - Same as t.test but change the output to reflect the comparison between the variances.
  * poisson.test
    - Just fit a standard poisson distribution. The poisson.test is quite limited...
  
I should write somewhere that the package uses modified code from the coda and the BEST packages.
