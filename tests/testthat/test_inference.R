context("inference")

# set a seed so the RNG aspects of this test don't cause random
#   failures
set.seed(251521)

# setup: define a pop model to simulate and use in ABC
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
  nsim = 100,
  args = list(density_dependence = list(theta = 0.4))
)

# repeat but with many more replicates
obs_long <- simulate(
  pop_fn,
  nsim = 1000,
  args = list(density_dependence = list(theta = 0.4))
)

# define wrapper for a quick simulation from parameters alone
popsim_fn <- function(par) {

  tmp <- update(
    pop_fn,
    aae.pop::density_dependence(
      funs = aae.pop::ricker(par[1], exclude = 1),
      masks = aae.pop::reproduction(popmat)
    )
  )

  sim <- simulate(tmp, nsim = 10, args = list(density_dependence = list(theta = par[2])))

  stat_abundance_trend(sim)

}

# and repeat this but for more replicates
popsim_fn_long <- function(par) {

  tmp <- update(
    pop_fn,
    aae.pop::density_dependence(
      funs = aae.pop::ricker(par[1], exclude = 1),
      masks = aae.pop::reproduction(popmat)
    )
  )

  sim <- simulate(tmp, nsim = 100, args = list(density_dependence = list(theta = par[2])))

  stat_abundance_trend(sim)

}

# calculate the target from the simulated data set
target <- stat_abundance_trend(obs)
target_long <- stat_abundance_trend(obs_long)

test_that("inference creates an ABC output", {

  # quick (inaccurate) inference
  inference_test <- inference(
    model = popsim_fn,
    prior = list(
      c("unif", 50, 500),
      c("unif", 0.01, 1)
    ),
    target = target,
    nb_simul = 10
  )

  # quick inference using full function call
  easyabc_test <- ABC_sequential(
    method = "Lenormand",
    model = popsim_fn,
    prior = list(
      c("unif", 50, 500),
      c("unif", 0.01, 1)
    ),
    nb_simul = 10,
    summary_stat_target = target,
    progress_bar = TRUE
  )

  # check inference has the same output structure as an EasyABC call
  expect_equal(names(inference_test), names(easyabc_test))

})


test_that("inference generates an approximately correct output", {

  # longer (more accurate) inference
  inference_test <- inference(
    model = popsim_fn_long,
    prior = list(
      c("unif", 50, 500),
      c("unif", 0.01, 1)
    ),
    target = target_long,
    nb_simul = 100
  )

  # linear regression of scaled parameter estimates should have
  #   slope near 1
  test_val <- coef(
    lm(
      scale(inference_test$param[, 1]) ~ -1 + scale(inference_test$param[, 2])
    )
  )

  # check inference has the same output structure as an EasyABC call
  expect_gt(test_val, 0.95)
  expect_lt(test_val, 1.05)

})

test_that("inference errors informatively", {

  # with incorrect model input
  expect_error(
    inference(
      model = target,
      prior = list(
        c("unif", 50, 500),
        c("unif", 0.01, 1)
      ),
      target = target,
      nb_simul = 10
    ),
    "must be a function"
  )

  # with incorrect prior class
  expect_error(
    inference(
      model = popsim_fn,
      prior = c(
        c("unif", 50, 500),
        c("unif", 0.01, 1)
      ),
      target = target,
      nb_simul = 10
    ),
    "must be a list"
  )

  # with incorrect prior values for a given distn
  expect_error(
    inference(
      model = popsim_fn,
      prior = list(
        c("unif", 50),
        c("unif", 0.01, 1)
      ),
      target = target,
      nb_simul = 10
    ),
    "at least one elemenet of prior has the wrong length"
  )

  # with un-implemented prior distribution
  expect_error(
    inference(
      model = popsim_fn,
      prior = list(
        c("labradoodle", 50, 500),
        c("unif", 0.01, 1)
      ),
      target = target,
      nb_simul = 10
    ),
    "prior contains unknown distribution"
  )

  # with mismatched model outputs and target
  expect_error(
    inference(
      model = popsim_fn,
      prior = list(
        c("unif", 50, 500),
        c("unif", 0.01, 1)
      ),
      target = c(target, 1),
      nb_simul = 10
    ),
    "must return an output with the same length as target"
  )

})
