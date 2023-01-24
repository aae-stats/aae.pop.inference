context("define")

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
obs_long <- simulate(
  pop_fn,
  nsim = 1000,
  args = list(density_dependence = list(theta = 0.4))
)

test_that("define works with fixed priors", {

  # define the inference model
  pop_inf <- define(x = pop_fn, nsim = 5)

  # and run it
  pars <- aae.pop.inference::inference(
    model = pop_inf$model,
    prior = pop_inf$prior,
    target = pop_inf$stat(obs),
    nb_simul = 100,
    progress_bar = FALSE,
    p_acc_min = 0.3
  )

  # just expect an output (not accurate) with correct number of pars and sims
  expect_equal(13L, ncol(pars$param))

})

test_that("define works with custom masks", {

  # define the inference model
  pop_inf <- define(
    x = pop_fn,
    masks = list(
      survival = aae.pop::transition(pop_fn$matrix),
      reproduction = aae.pop::reproduction(pop_fn$matrix)
    ),
    nsim = 5
  )

  # and run it
  pars <- aae.pop.inference::inference(
    model = pop_inf$model,
    prior = pop_inf$prior,
    target = pop_inf$stat(obs),
    nb_simul = 100,
    progress_bar = FALSE,
    p_acc_min = 0.3
  )

  # just expect an output (not accurate) with correct number of pars and sims
  expect_equal(8L, ncol(pars$param))

})

test_that("define works with custom priors", {

  # define the inference model
  pop_inf <- define(
    x = pop_fn,
    masks = list(
      survival = aae.pop::transition(pop_fn$matrix),
      reproduction = aae.pop::reproduction(pop_fn$matrix)
    ),
    priors = list(
      survival = list(
        c("unif", 0.2, 0.4),
        c("unif", 0.3, 0.6),
        c("unif", 0.4, 0.9),
        c("unif", 0.5, 1)
      ),
      reproduction = list(
        c("lognormal", 0, 1),
        c("lognormal", 0, 1),
        c("lognormal", 0, 3),
        c("lognormal", 0, 3)
      )
    ),
    nsim = 5
  )

  # and run it
  pars <- aae.pop.inference::inference(
    model = pop_inf$model,
    prior = pop_inf$prior,
    target = pop_inf$stat(obs),
    nb_simul = 100,
    progress_bar = FALSE,
    p_acc_min = 0.3
  )

  # just expect an output (not accurate) with correct number of pars and sims
  expect_equal(8L, ncol(pars$param))

})

test_that("define works accurately with custom priors", {

  # define the inference model
  pop_inf <- define(
    x = pop_fn,
    masks = list(
      survival = aae.pop::transition(pop_fn$matrix),
      reproduction = aae.pop::reproduction(pop_fn$matrix)
    ),
    priors = list(
      survival = list(
        c("unif", 0.2, 0.4),
        c("unif", 0.3, 0.6),
        c("unif", 0.4, 0.9),
        c("unif", 0.5, 1)
      ),
      reproduction = list(
        c("lognormal", 0, 1),
        c("lognormal", 0, 1),
        c("lognormal", 0, 3),
        c("lognormal", 0, 3)
      )
    ),
    nsim = 100,
    args = list(density_dependence = list(theta = 0.4))
  )

  # and run it
  pars <- aae.pop.inference::inference(
    model = pop_inf$model,
    prior = pop_inf$prior,
    target = pop_inf$stat(obs_long),
    nb_simul = 200,
    progress_bar = FALSE,
    p_acc_min = 0.05
  )

  # get true params
  target <- c(
    pop_fn$matrix[transition(pop_fn$matrix)],
    pop_fn$matrix[reproduction(pop_fn$matrix)]
  )

  # get modelled quantiles to see if true values are correlated
  value <- apply(pars$param, 2, quantile, p = 0.5)

  # just expect an output (not accurate) with correct number of pars and sims
  expect_gt(coef(lm(value ~ -1 + target)), 0.95)
  expect_lt(coef(lm(value ~ -1 + target)), 1.05)

})

test_that("define works with aae.pop args and options", {

  # define the inference model
  pop_inf <- define(
    x = pop_fn,
    masks = list(
      survival = aae.pop::transition(pop_fn$matrix),
      reproduction = aae.pop::reproduction(pop_fn$matrix)
    ),
    nsim = 5,
    args = list(density_dependence = list(theta = 0.4))
  )

  # and run it
  pars <- aae.pop.inference::inference(
    model = pop_inf$model,
    prior = pop_inf$prior,
    target = pop_inf$stat(obs),
    nb_simul = 100,
    progress_bar = FALSE,
    p_acc_min = 0.3
  )

  # just expect an output (not accurate) with correct number of pars and sims
  expect_equal(8L, ncol(pars$param))

  # define the inference model
  pop_inf <- define(
    x = pop_fn,
    masks = list(
      survival = aae.pop::transition(pop_fn$matrix),
      reproduction = aae.pop::reproduction(pop_fn$matrix)
    ),
    nsim = 5,
    args = list(density_dependence = list(theta = 0.4)),
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )

  # and run it
  pars <- aae.pop.inference::inference(
    model = pop_inf$model,
    prior = pop_inf$prior,
    target = pop_inf$stat(obs),
    nb_simul = 100,
    progress_bar = FALSE,
    p_acc_min = 0.3
  )

  # just expect an output (not accurate) with correct number of pars and sims
  expect_equal(8L, ncol(pars$param))

})

test_that("define works correctly with different summary functions", {

  # define the inference model
  pop_inf <- define(
    x = pop_fn,
    masks = list(
      survival = aae.pop::transition(pop_fn$matrix),
      reproduction = aae.pop::reproduction(pop_fn$matrix)
    ),
    stat = aae.pop.inference::stat_abundance_moment,
    nsim = 5,
    args = list(density_dependence = list(theta = 0.4))
  )

  # and run it
  pars <- aae.pop.inference::inference(
    model = pop_inf$model,
    prior = pop_inf$prior,
    target = pop_inf$stat(obs),
    nb_simul = 100,
    progress_bar = FALSE,
    p_acc_min = 0.3
  )

  # just expect an output (not accurate) with correct number of pars and sims
  expect_equal(8L, ncol(pars$param))

})

test_that("define errors informatively with incorrect priors", {

  # check survival as a list of incorrect dims
  expect_error(
    define(
      x = pop_fn,
      priors = list(
        survival = lapply(seq_len(25), function(i) c("unif", 1, 10))
      )
    ),
    "survival prior is a list with 25"
  )

  # check reproduction as a list of incorrect dims
  expect_error(
    define(
      x = pop_fn,
      priors = list(
        reproduction = lapply(seq_len(25), function(i) c("unif", 1, 10))
      )
    ),
    "reproduction prior is a list with 25"
  )

  # check survival as a vector of wrong class
  expect_error(
    define(
      x = pop_fn,
      priors = list(
        survival = c(1, 2, 3)
      )
    ),
    "priors must be character vectors or lists of character vectors"
  )

  # check reproduction as a vector of wrong class
  expect_error(
    define(
      x = pop_fn,
      priors = list(
        reproduction = c(1, 2, 3)
      )
    ),
    "priors must be character vectors or lists of character vectors"
  )

})

test_that("define errors informatively without dynamics input", {

  # check survival as a list of incorrect dims
  expect_error(
    define(x = rnorm(10)),
    "must be a dynamics object"
  )

})

test_that("print and is methods work correctly", {

  # define the inference model
  pop_inf <- define(
    x = pop_fn,
    masks = list(
      survival = aae.pop::transition(pop_fn$matrix),
      reproduction = aae.pop::reproduction(pop_fn$matrix)
    ),
    stat = aae.pop.inference::stat_abundance_moment,
    nsim = 5,
    args = list(density_dependence = list(theta = 0.4))
  )

  # check print returns an output
  expect_output(print(pop_inf), "Specified model and priors for use")

  # check the `is` method works
  expect_true(is.definition(pop_inf))

})
