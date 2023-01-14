
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aae.pop.inference

<!-- badges: start -->

![R-CMD-check](https://github.com/aae-stats/aae.pop.inference/actions/workflows/check-standard.yaml/badge.svg)
[![Codecov test
coverage](https://codecov.io/github/aae-stats/aae.pop.inference/main/graph/badge.svg)](https://app.codecov.io/github/aae-stats/aae.pop.inference/)
<!-- badges: end -->

`aae.pop.inference` supports inference on population dynamics models
implemented in `aae.pop`. These models are relatively complex, making
standard inference methods (e.g. maximum likelihood, MCMC) challenging,
especially as the models become realistically complex (many population
classes, temporal variability). `aae.pop.inference` uses a sequential
approximate Bayesian computation approach to support fast inference on
these complex models.

## Installation

You can install the released version of aae.pop.inference from
[GitHub](https://github.com/aae-stats/aae.pop.inference) with:

``` r
remotes::install_github("aae-stats/aae.pop.inference")
```

## Example: estimating carrying capacity and density dependence

This example uses data simulated from a model as a target to estimate
the parameters of this model. This relatively simple example (estimating
parameters of a known model) illustrates how to perform inference on
these models but also highlights a case in which the target parameters
are highly correlated.

``` r
# some packages needed
library(aae.pop)
library(aae.pop.inference)
#> Loading required package: EasyABC
#> Loading required package: abc
#> Loading required package: abc.data
#> Loading required package: nnet
#> Loading required package: quantreg
#> Loading required package: SparseM
#> 
#> Attaching package: 'SparseM'
#> The following object is masked from 'package:base':
#> 
#>     backsolve
#> Loading required package: MASS
#> Loading required package: locfit
#> locfit 1.5-9.7    2023-01-02

# define a basic population model with aae.pop
nclass <- 5
popmat <- matrix(0, ncol = nclass, nrow = nclass)
popmat[transition(popmat)] <- c(0.3, 0.5, 0.7, 0.8)
popmat[reproduction(popmat)] <- c(0, 1, 4, 10)
dd <- density_dependence(
  funs = ricker(100, exclude = 1),
  masks = reproduction(popmat)
)
pop_fn <- dynamics(matrix = popmat, dd)

# simulate some data from this model to use as a target data set
obs <- simulate(
  pop_fn,
  nsim = 1000,
  args = list(density_dependence = list(theta = 0.4))
)

# define wrapper so we can simulate from the model based on
#   the inference parameters alone
popsim_fn <- function(par) {

  tmp <- update(
    pop_fn,
    density_dependence(
      funs = ricker(par[1], exclude = 1),
      masks = reproduction(popmat)
    )
  )

  sim <- simulate(tmp, nsim = 100, args = list(density_dependence = list(theta = par[2])))

  stat_abundance_trend(sim)

}

# calculate the target from the simulated data set
target <- stat_abundance_trend(obs)

# estimate parameters by comparing the simulation model to
#   the simulated data set with uniform priors on 
#   carrying capacity (U[50, 500]) and the theta parameter
#   of the Ricker density dependence (U[0.01, 1])
abc_est <- inference(
  model = popsim_fn,
  prior = list(
    c("unif", 50, 500),
    c("unif", 0.01, 1)
  ),
  target = target,
  nb_simul = 100,
  progress_bar = FALSE
)
```

Plotting the estimated posteriors for the two parameters suggests
relatively broad distributions with high uncertainty. However, a plot of
these two parameters against one another highlights their close
relationship, with a near perfect linear correlation between the two.

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-2-2.png" width="100%" />

## Issues

Please leave feedback, bug reports or feature requests at the GitHub
[issues page](https://github.com/aae-stats/aae.pop.inference/issues)
