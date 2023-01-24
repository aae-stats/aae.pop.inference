#' @name statistics
#' @title Calculate summary statistics from population trajectories
#' @description Functions to define summary statistics based on
#'   observed or simulated population dynamics, which can be used to
#'   support inference on population dynamics models using
#'   approximate Bayesian computation.
NULL

#' @rdname statistics
#'
#' @export
#'
#' @importFrom stats sd median quantile var
#'
#' @param x an object returned from
#'   \code{\link[aae.pop:simulate]{aae.pop::simulate()}} containing
#'   simulated population trajectories
#' @param classes a list of integer vectors specifying which how classes
#'   should be grouped to calculate abundances. Groupings do not need to
#'   be disjoint, that is, classes can be included in multiple groups.
#'   Defaults to \code{list(seq_len(nclass), c(2:nclass), c(1))}, which
#'   returns summed abundances of all individuals, all individuals except
#'   the first class, and the first class only
#' @param zscale logical determining whether summary statistics are z-scaled.
#'   Defaults to \code{TRUE}
#' @param fn function or vector of functions defining how abundances are
#'   summarised over all replicates. Defaults to \code{c(median, sd)},
#'   which calculates the median and standard deviation of abundances
#'   over replicate trajectories
#' @param \dots additional arguments passed to fn
#'
#' @details \code{stat_abundance_trend} returns the abundances of different
#'   groupings of population classes in each year, averaged over replicate
#'   trajectories according to \code{fn}. The default is to calculate the
#'   median and standard deviation of abundances over all replicates, with
#'   groupings representing all individuals, all individuals except the first
#'   class (often young-of-year recruits), and the first class only.
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
#' # calculate summary statistics by year
#' stat_abundance_trend(obs)
#'
#' # calculate summary statistics as moments over all years
#' stat_abundance_moment(obs)
stat_abundance_trend <- function(
  x,
  classes = NULL,
  zscale = TRUE,
  fn = c(median, sd),
  ...
) {

  # need to summarise matrices and arrays differently, so check
  #   dims of x
  if (!aae.pop::is.simulation(x)) {
    stop(
      "x must be a simulation object returned by aae.pop::simulate()",
      call. = FALSE
    )
  }

  # check fn and wrap it in a list if it's an isolated function to make
  #   the loop below work cleanly
  if (!is.list(fn)) {
    if (!is.function(fn))
      stop("fn must be a function or list of functions", call. = FALSE)
  }
  if (is.list(fn)) {
    if (!all(sapply(fn, class) == "function"))
      stop("if fn is a list, all elements must be functions", call. = FALSE)
  }
  if (is.function(fn))
    fn <- list(fn)

  # set default classes if not specified
  if (is.null(classes))
    classes <- list(seq_len(ncol(x)), seq_len(ncol(x))[-1], c(1))

  # calculate abundances of each grouping of classes
  abund <- lapply(
    classes,
    function(idx, .x) apply(.x[, idx, , drop = FALSE], c(1, 3), sum), .x = x
  )

  # and then summarise over replicates based on fn
  stat <- vector("list", length = length(fn))
  for (i in seq_along(fn)) {
    stat[[i]] <- unlist(
      lapply(abund, function(.x, .y) apply(.x, 2, .y, ...), .y = fn[[i]])
    )
  }

  # need this as a single vector for ABC_sequential
  stat <- unlist(stat)

  # optionally zscale to keep scale of stats similar
  if (zscale)
    stat <- (stat - mean(stat)) / sd(stat)

  # return
  stat

}

#' @rdname statistics
#'
#' @export
#'
#' @importFrom moments skewness kurtosis
#'
#' @param x an object returned from
#'   \code{\link[aae.pop:simulate]{aae.pop::simulate()}} containing
#'   simulated population trajectories
#' @param classes a list of integer vectors specifying which how classes
#'   should be grouped to calculate abundances. Groupings do not need to
#'   be disjoint, that is, classes can be included in multiple groups.
#'   Defaults to \code{list(seq_len(nclass), c(2:nclass), c(1))}, which
#'   returns summed abundances of all individuals, all individuals except
#'   the first class, and the first class only
#' @param zscale logical determining whether summary statistics are z-scaled.
#'   Defaults to \code{TRUE}
#' @param fn function or vector of functions defining how abundances are
#'   summarised over all replicates. Defaults to \code{c(median, sd)},
#'   which calculates the median and standard deviation of abundances
#'   over replicate trajectories
#' @param \dots additional arguments passed to fn
#'
#' @details \code{stat_abundance_moment} returns the moments over all years
#'   of population trajectories. Populations are grouped according to
#'   \code{classes}, as for \code{stat_abundance_trend}. Moments are
#'   averaged over all replicates according to \code{fn}, which defaults
#'   to the median and standard deviation. The first four moments are
#'   included (mean, variance, skewness, kurtosis).
#'
stat_abundance_moment <- function(
  x,
  classes = NULL,
  zscale = TRUE,
  fn = c(median, sd),
  ...
) {

  # need to summarise matrices and arrays differently, so check
  #   dims of x
  if (!aae.pop::is.simulation(x)) {
    stop(
      "x must be a simulation object returned by aae.pop::simulate()",
      call. = FALSE
    )
  }

  # check fn and wrap it in a list if it's an isolated function to make
  #   the loop below work cleanly
  if (!is.list(fn)) {
    if (!is.function(fn))
      stop("fn must be a function or list of functions", call. = FALSE)
  }
  if (is.list(fn)) {
    if (!all(sapply(fn, class) == "function"))
      stop("if fn is a list, all elements must be functions", call. = FALSE)
  }
  if (is.function(fn))
    fn <- list(fn)

  # set default classes if not specified
  if (is.null(classes))
    classes <- list(seq_len(ncol(x)), seq_len(ncol(x))[-1], c(1))

  # calculate abundances of each grouping of classes
  abund <- lapply(
    classes,
    function(idx, .x) apply(.x[, idx, , drop = FALSE], c(1, 3), sum), .x = x
  )

  # and then summarise to give final summary stat
  stat <- vector("list", length = length(fn))
  for (i in seq_along(fn)) {

    # first summarise over time steps (keep replicates)
    stat[[i]] <- lapply(abund, function(.x) apply(.x, 1, moment_fn))

    # then flatten over replicates
    stat[[i]] <- unlist(
      lapply(stat[[i]], function(.x, .y) apply(.x, 1, .y, ...), .y = fn[[i]])
    )
  }

  # need this as a single vector for ABC_sequential
  stat <- unlist(stat)

  # optionally zscale to keep scale of stats similar
  if (zscale)
    stat <- (stat - mean(stat)) / sd(stat)

  # return
  stat

}

# internal function: wrapper to calculate moments in a single function
moment_fn <- function(x, ...) {
  c(
    mean(x, ...),
    var(x, ...),
    moments::skewness(x, ...),
    moments::kurtosis(x, ...)
  )
}
