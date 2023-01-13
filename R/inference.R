#' @name inference
#' @title Estimate parameters of a simulation model with sequential
#'   approximate Bayesian computation
#' @description Generic function to estimate parameters using sequential
#'   ABC, targeted specifically towards models implemented with the
#'   \code{\link{simulate}} function in the \code{aae.pop} package.
NULL

#' @rdname inference
#'
#' @export
#'
#' @importFrom EasyABC ABC_sequential
#' @importFrom stats runif rlnorm rnorm rexp quantile
#'
#' @param model a function to simulate a population dynamics object
#'   from a set of parameters. This function must return a
#'   vector of summary statistics with the same length as \code{target}
#' @param prior a list of prior distributions for each parameter, specified
#'   as described in
#'   \code{\link[EasyABC:ABC_sequential]{EasyABC::ABC_sequential()}}.
#'   Currently supported options are the uniform, normal, lognormal,
#'   and exponential distributions
#' @param target values to be compared to model outputs simulated
#'   with \code{model}
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
#' # to add
inference <- function(
    model,
    prior,
    target,
    method = NULL,
    ...
) {

  # the simulation method must be a function
  if(!is.function(model))
    stop("model must be a function", call. = FALSE)

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
  do.call(distn, lapply(seq_along(x)[-1], function(i) x[i]))

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
#' @export
# nolint start
summary.inference <- function(x, ...) {
  # nolint end

  # calculate quantiles of the estimated posterior distribution
  param_dist <- apply(
    x$param,
    2,
    quantile,
    pr = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
  )
  colnames(param_dist) <- paste0("Parameter ", seq_len(ncol(x$param)))
  rownames(param_dist) <- paste0(
    c("5th", "10th", "25th", "50th", "75th", "90th", "95th"),
    " percentile"
  )

  # print these quantiles with info on compute time
  cat(paste0(
    "Inference on ", ncol(x$param),
    " model parameters took ", x$computime,
    " seconds. Posterior distributions for each parameter ",
    " had the following quantiles:\n"
  ))
  print(param_dist)

  # and return silently
  outputs <- list(
    quantiles = param_dist
  )

}

# S3 plot method
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

  # work out plotting dimensions
  sqrt_npar <- sqrt(length(which))
  ncol <- floor(sqrt_npar)
  nrow <- ceiling(sqrt_npar)
  par(mfrow = c(nrow, ncol))
  for (i in which)
    hist(x$param[[i]], ...)

}

# internal function: set inference class
as_inference <- function(x) {
  type <- "array"
  if (is.matrix(x))
    type <- "matrix"
  as_class(x, name = "inference", type = type)
}
