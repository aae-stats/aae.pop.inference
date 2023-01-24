#' @name inference
#' @title Estimate parameters of a simulation model with sequential
#'   approximate Bayesian computation
#' @description Generic function to estimate parameters using sequential
#'   ABC, targeted specifically towards models implemented with the
#'   \code{\link[aae.pop:simulate]{aae.pop::simulate()}} function
#'   in the \code{aae.pop} package.
NULL

#' @rdname inference
#'
#' @export
#'
#' @importFrom EasyABC ABC_sequential
#' @importFrom stats runif rlnorm rnorm rexp quantile
#' @importFrom graphics hist par
#'
#' @param model a function to simulate a population dynamics object
#'   from a set of parameters. This function must return a
#'   vector of summary statistics with the same length as \code{target}
#' @param prior a list of prior distributions for each parameter, specified
#'   as described in
#'   \code{\link[EasyABC:ABC_sequential]{EasyABC::ABC_sequential()}}.
#'   Currently supported options are the uniform, normal, lognormal,
#'   and exponential distributions, and each list element should be a
#'   vector specifying the distribution ("unif", "normal", "lognormal",
#'   or "exponential"), followed by one (for the exponential distribution)
#'   or two parameters (all other distributions)
#' @param target values to be compared to model outputs simulated
#'   with \code{model}
#' @param method sequential ABC method to use for sampling, see
#'   \code{\link[EasyABC:ABC_sequential]{EasyABC::ABC_sequential()}} for
#'   details
#' @param x output from \code{inference}
#' @param y ignored, included for consistency with generic plot method
#' @param which the parameter to plot, defaults to all estimated parameters
#' @param object output from \code{inference}
#' @param \dots additional arguments passed to the sequential ABC sampler
#'   (see \code{\link[EasyABC:ABC_sequential]{EasyABC::ABC_sequential()}})
#'
#' @details Approximate Bayesian computation is a fast method for parameter
#'   inference in cases where model likelihoods are intractable (or at
#'   least difficult to compute). The \code{inference} function is a wrapper
#'   for the \code{\link[EasyABC:ABC_sequential]{EasyABC::ABC_sequential()}}
#'   function intended to simplify parameter inference for population
#'   dynamics models created with
#'   \code{\link[aae.pop:dynamics]{aae.pop::dynamics()}}.
#'
#'   Although standard MCMC algorithms are feasible for many of these models,
#'   the use of ABC allows fast and flexible specification of models based on
#'   existing simulation tools such as those in the \code{aae.pop} package.
#'
#' @examples
#' # use the tools in aae.pop to set up an example
#' library(aae.pop)
#'
#' # setup: define a pop model to simulate and use in ABC
#' nclass <- 5
#' popmat <- matrix(0, ncol = nclass, nrow = nclass)
#' popmat[transition(popmat)] <- c(0.3, 0.5, 0.7, 0.8)
#' popmat[reproduction(popmat)] <- c(0, 1, 4, 10)
#' dd <- density_dependence(
#'   funs = ricker(100, exclude = 1),
#'   masks = reproduction(popmat)
#' )
#' pop_fn <- dynamics(matrix = popmat, dd)
#'
#' # simulate some data from this model to use as a target data set
#' obs <- simulate(
#'   pop_fn,
#'   nsim = 100,
#'   args = list(density_dependence = list(theta = 0.4))
#' )
#'
#' # define wrapper for a simulation from K and DD parameters
#' popsim_fn <- function(par) {
#'
#'   tmp <- update(
#'     pop_fn,
#'     aae.pop::density_dependence(
#'       funs = aae.pop::ricker(par[1], exclude = 1),
#'       masks = aae.pop::reproduction(popmat)
#'     )
#'   )
#'
#'   sim <- simulate(
#'     tmp,
#'     nsim = 10,
#'     args = list(density_dependence = list(theta = par[2]))
#'   )
#'
#'   stat_abundance_trend(sim)
#'
#' }
#'
#' # calculate the target from the simulated data set
#' target <- stat_abundance_trend(obs)
#'
#' # and run it with uniform priors for both parameters (quick run,
#' #   nb_simul too low and p_acc_min too high for accurate results)
#' pars <- inference(
#'   model = popsim_fn,
#'   prior = list(
#'     c("unif", 50, 500),
#'     c("unif", 0.01, 1)
#'   ),
#'   target = target,
#'   nb_simul = 10,
#'   progress_bar = TRUE
#' )
#'
#' # plot estimated parameters (not expected to be accurate)
#' plot(pars)
inference <- function(
    model,
    prior,
    target,
    method = NULL,
    ...
) {

  # the simulation method must be a function and prior
  #   must be a list
  if (!is.function(model))
    stop("model must be a function", call. = FALSE)
  if (!is.list(prior))
    stop("prior must be a list", call. = FALSE)

  # and then check the prior is correctly specified for ABC_sequential
  check_prior(prior)

  # and the output from model should match the length of target
  prior_draw <- sapply(prior, sample_once)
  sim_test <- model(prior_draw)
  if (length(sim_test) != length(target)) {
    stop(
      "model must return an output with the same length as target",
      call. = FALSE
    )
  }

  # set baseline arguments and overwrite any that are set by user
  args <- list(
    model = model,
    prior = prior,
    summary_stat_target = target,
    method = "Lenormand",
    nb_simul = 1000,
    progress_bar = TRUE
  )
  args_dots <- list(...)
  args[names(args_dots)] <- args_dots

  # define and return inference output
  as_inference(do.call(EasyABC::ABC_sequential, args))

}

# internal function to check that the specified prior is defined
#   as required by EasyABC::ABC_sequential
check_prior <- function(x) {

  # check prior distn is implemented in ABC_sequential
  prior_distn <- sapply(x, function(x) x[1])
  if (!all(prior_distn %in% c("unif", "lognormal", "normal", "exponential"))) {
    stop(
      "prior contains unknown distribution; must be one of unif, normal, ",
      "lognormal, or exponential",
      call. = FALSE
    )
  }

  # check prior has the correct length for each type
  prior_ok <- TRUE
  for (i in seq_along(x)) {
    if (x[[i]][1] == "exponential") {
      args_ok <- length(x[[i]]) == 2 &
        !is.na(as.numeric(x[[i]][2]))
    } else {
      args_ok <- length(x[[i]]) == 3 &
        !any(is.na(as.numeric(x[[i]][2:3])))
    }
    prior_ok <- prior_ok & args_ok
  }

  # error if not OK
  if (!prior_ok) {
    stop(
      "at least one elemenet of prior has the wrong length or values; each ",
      " element should contain a distribution followed by one (for ",
      "exponential) or two numeric values (all other distributions)",
      call. = FALSE
    )
  }

}

# internal function to draw a single set of values from a prior defined
#   as described in EasyABC::ABC_sequential
sample_once <- function(x) {

  # specify a RNG for the appropriate distribution, defaulting
  #   to runif if not otherwise specified
  distn <- switch(
    x[1],
    "unif" = runif,
    "lognormal" = rlnorm,
    "normal" = rnorm,
    "exponential" = rexp,
    runif
  )

  # draw a single value from the relevant distribution
  do.call(distn, c(1, lapply(seq_along(x)[-1], function(i) as.numeric(x[i]))))

}

# S3 is method
#' @rdname inference
#' @export
# nolint start
is.inference <- function(x) {
  # nolint end
  inherits(x, "inference")
}

# S3 print method
#' @rdname inference
#' @export
# nolint start
print.inference <- function(x, ...) {
  # nolint end
  cat(
    paste0(
    "Outputs from sequential ABC inference ",
    "on population dynamics model\n"
    )
  )
}

# S3 summary method
#' @rdname inference
#' @export
# nolint start
summary.inference <- function(object, ...) {
  # nolint end

  # calculate quantiles of the estimated posterior distribution
  param_dist <- apply(
    object$param,
    2,
    quantile,
    pr = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
  )
  colnames(param_dist) <- paste0("Parameter ", seq_len(ncol(object$param)))
  rownames(param_dist) <- paste0(
    c("5th", "10th", "25th", "50th", "75th", "90th", "95th"),
    " percentile"
  )

  # print these quantiles with info on compute time
  cat(paste0(
    "Inference on ", ncol(object$param),
    " model parameters took ", round(object$computime, 3),
    " seconds. Posterior distributions for each parameter ",
    " had the following quantiles:\n"
  ))
  print(round(param_dist, 4))

}

# S3 plot method
#' @rdname inference
#' @export
# nolint start
plot.inference <- function(x, y, ..., which = NULL) {
  # nolint end

  # grab plotting pars and reset at end
  old_mfrow <- par()$mfrow
  on.exit(par(mfrow = old_mfrow))

  # how many parameters are there?
  nparam <- ncol(x$param)

  # plot all if a subset isn't specified
  if (is.null(which))
    which <- seq_len(nparam)

  # set some basic arguments to keep the plot clean
  args <- list(
    main = "",
    xlab = "Parameter",
    las = 1
  )

  # overwrite with `...` if required
  arg_list <- list(...)
  args[names(arg_list)] <- arg_list

  # work out plotting dimensions
  sqrt_npar <- sqrt(length(which))
  ncol <- floor(sqrt_npar)
  nrow <- ceiling(sqrt_npar)
  par(mfrow = c(nrow, ncol))
  for (i in which)
    do.call(hist, c(list(x$param[, i]), args))

}

# internal function: set inference class
as_inference <- function(x) {
  as_class(x, name = "inference", type = "list")
}
