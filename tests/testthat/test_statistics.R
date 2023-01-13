context("statistics")

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

test_that("stat_abundance_trend returns correct values under defaults", {

  # calculate values from data
  value <- stat_abundance_trend(obs)

  # and calculate targets for each set of classes
  target1 <- apply(obs, c(1, 3), sum)
  target1 <- c(apply(target1, 2, median), apply(target1, 2, sd))
  target2 <- apply(obs[, 2:5, ], c(1, 3), sum)
  target2 <- c(apply(target2, 2, median), apply(target2, 2, sd))
  target3 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target3 <- c(apply(target3, 2, median), apply(target3, 2, sd))

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))

})

test_that("stat_abundance_trend works with different classes", {

  # calculate values from data
  value <- stat_abundance_trend(obs, classes = list(c(1), c(4:5), c(2)))

  # and calculate targets for each set of classes
  target1 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target1 <- c(apply(target1, 2, median), apply(target1, 2, sd))
  target2 <- apply(obs[, 4:5, ], c(1, 3), sum)
  target2 <- c(apply(target2, 2, median), apply(target2, 2, sd))
  target3 <- apply(obs[, 2, , drop = FALSE], c(1, 3), sum)
  target3 <- c(apply(target3, 2, median), apply(target3, 2, sd))

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))

})

test_that("stat_abundance_trend works correctly with different fns", {

  # calculate values from data
  value <- stat_abundance_trend(obs, fn = mean)

  # and calculate targets for each set of classes
  target1 <- apply(obs, c(1, 3), sum)
  target1 <- apply(target1, 2, mean)
  target2 <- apply(obs[, 2:5, ], c(1, 3), sum)
  target2 <- apply(target2, 2, mean)
  target3 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target3 <- apply(target3, 2, mean)

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))

  # calculate values from data
  value <- stat_abundance_trend(obs, fn = var)

  # and calculate targets for each set of classes
  target1 <- apply(obs, c(1, 3), sum)
  target1 <- apply(target1, 2, var)
  target2 <- apply(obs[, 2:5, ], c(1, 3), sum)
  target2 <- apply(target2, 2, var)
  target3 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target3 <- apply(target3, 2, var)

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))

  # calculate values from data
  value <- stat_abundance_trend(obs, fn = list(mean, sd, var, median))

  # and calculate targets for each set of classes
  target1 <- apply(obs, c(1, 3), sum)
  target1 <- c(
    apply(target1, 2, mean),
    apply(target1, 2, sd),
    apply(target1, 2, var),
    apply(target1, 2, median)
  )
  target2 <- apply(obs[, 2:5, ], c(1, 3), sum)
  target2 <- c(
    apply(target2, 2, mean),
    apply(target2, 2, sd),
    apply(target2, 2, var),
    apply(target2, 2, median)
  )
  target3 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target3 <- c(
    apply(target3, 2, mean),
    apply(target3, 2, sd),
    apply(target3, 2, var),
    apply(target3, 2, median)
  )

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))

})

test_that("stat_abundance_moment returns correct values under defaults", {

  # calculate values from data
  value <- stat_abundance_moment(obs)

  # and calculate targets for each set of classes
  target1 <- apply(obs, c(1, 3), sum)
  target1 <- apply(
    target1,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target1 <- c(apply(target1, 1, median), apply(target1, 1, sd))
  target2 <- apply(obs[, 2:5, ], c(1, 3), sum)
  target2 <- apply(
    target2,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target2 <- c(apply(target2, 1, median), apply(target2, 1, sd))
  target3 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target3 <- apply(
    target3,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target3 <- c(apply(target3, 1, median), apply(target3, 1, sd))

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))

})

test_that("stat_abundance_moment works with different classes", {

  # calculate values from data
  value <- stat_abundance_moment(obs, classes = list(c(1), c(4:5), c(2)))

  # and calculate targets for each set of classes
  target1 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target1 <- apply(
    target1,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target1 <- c(apply(target1, 1, median), apply(target1, 1, sd))
  target2 <- apply(obs[, 4:5, ], c(1, 3), sum)
  target2 <- apply(
    target2,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target2 <- c(apply(target2, 1, median), apply(target2, 1, sd))
  target3 <- apply(obs[, 2, , drop = FALSE], c(1, 3), sum)
  target3 <- apply(
    target3,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target3 <- c(apply(target3, 1, median), apply(target3, 1, sd))

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))

})

test_that("stat_abundance_moment works with different fns", {

  # calculate values from data
  value <- stat_abundance_moment(obs, fn = var)

  # and calculate targets for each set of classes
  target1 <- apply(obs, c(1, 3), sum)
  target1 <- apply(
    target1,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target1 <- apply(target1, 1, var)
  target2 <- apply(obs[, 2:5, ], c(1, 3), sum)
  target2 <- apply(
    target2,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target2 <- apply(target2, 1, var)
  target3 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target3 <- apply(
    target3,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target3 <- apply(target3, 1, var)

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))


  # calculate values from data
  value <- stat_abundance_moment(obs, fn = mean)

  # and calculate targets for each set of classes
  target1 <- apply(obs, c(1, 3), sum)
  target1 <- apply(
    target1,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target1 <- apply(target1, 1, mean)
  target2 <- apply(obs[, 2:5, ], c(1, 3), sum)
  target2 <- apply(
    target2,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target2 <- apply(target2, 1, mean)
  target3 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target3 <- apply(
    target3,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target3 <- apply(target3, 1, mean)

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))

  # calculate values from data
  value <- stat_abundance_moment(obs, fn = list(sum, mean, var))

  # and calculate targets for each set of classes
  target1 <- apply(obs, c(1, 3), sum)
  target1 <- apply(
    target1,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target1 <- c(apply(target1, 1, sum), apply(target1, 1, mean), apply(target1, 1, var))
  target2 <- apply(obs[, 2:5, ], c(1, 3), sum)
  target2 <- apply(
    target2,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target2 <- c(apply(target2, 1, sum), apply(target2, 1, mean), apply(target2, 1, var))
  target3 <- apply(obs[, 1, , drop = FALSE], c(1, 3), sum)
  target3 <- apply(
    target3,
    1,
    function(x) c(mean(x), var(x), moments::skewness(x), moments::kurtosis(x))
  )
  target3 <- c(apply(target3, 1, sum), apply(target3, 1, mean), apply(target3, 1, var))

  # check statistic calcs match manual calcs
  expect_equal(value, c(target1, target2, target3))

})

test_that("stat_abundance_trend and _moment error informatively", {

  # when input is not a simulation object
  expect_error(
    stat_abundance_moment(rnorm(10)),
    "must be a simulation object"
  )

  # when fn is not a function
  expect_error(
    stat_abundance_moment(obs, fn = rnorm(10)),
    "fn must be a function"
  )

  # or a list of functions
  expect_error(
    stat_abundance_moment(obs, fn = list(mean, rnorm(10))),
    "all elements must be functions"
  )

})
