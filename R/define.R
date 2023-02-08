#' @name definition
#' @title Specify a population model and priors for inference on vital rates
#' @description Return a simulation function and suitable priors to perform
#'   inference on the vital rates component of a stochastic matrix population
#'   model with the \code{aae.pop.inference} package.
NULL

#' @rdname definition
#'
#' @export
#'
#' @importFrom stats simulate
#' @import aae.pop
#'
#' @param x an object returned from
#'   \code{\link[aae.pop:dynamics]{aae.pop::dynamics()}} containing
#'   a stochastic matrix population model
#' @param masks a logical matrix or vector (or list of these) defining cells
#'   to be estimated. See \code{\link[aae.pop:masks]{aae.pop::masks()}} for
#'   details
#' @param priors a named list of prior distributions, with one element for
#'   survival and another for reproduction. Each element can be a single
#'   character vector (see Details) or a list of character vectors. The
#'   former assumes all survival/reproduction priors are identical, while
#'   the latter allows different priors for each parameter
#' @param stat a function to summarise the output of a population dynamics
#'   model run with \code{\link[aae.pop:simulate]{aae.pop::simulate()}}.
#'   Defaults to
# nolint start
#'   \code{\link[aae.pop.inference:stat_abundance_trend]{aae.pop.inference::stat_abundance_trend()}}
# nolint end
#' @param sim_args additional arguments passed to
#'   \code{\link[aae.pop:simulate]{aae.pop::simulate()}}
#' @param \dots additional arguments passed to \code{stat}
#'
#' @details \code{define} returns a pre-defined template that can be used with
#'   \code{\link[aae.pop.inference:inference]{aae.pop.inference::inference()}}
#'   to estimate the demographic parameters of a stochastic population model.
#'   More complex models will require a custom specification but \code{define}
#'   is a useful placeholder in simpler cases where model arguments are fixed
#'   and only the demographic vital rates are variable.
#'
#'   Priors should be specified as described in
#'   \code{\link[EasyABC:ABC_sequential]{EasyABC::ABC_sequential()}}.
#'   Currently supported options are the uniform, normal, lognormal,
#'   and exponential distributions, and each list element should be a
#'   vector specifying the distribution ("unif", "normal", "lognormal",
#'   or "exponential"), followed by one (for the exponential distribution)
#'   or two parameters (all other distributions)
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
#' # define the inference model for all non-zero parameters
#' #   (nsim too low for accurate results)
#' pop_inf <- define(
#'   x = pop_fn,
#'   masks = list(
#'     survival = aae.pop::transition(pop_fn$matrix),
#'     reproduction = aae.pop::reproduction(pop_fn$matrix)
#'   ),
#'   stat = stat_abundance_trend,
#'   sim_args = list(nsim = 5),
#'   classes = list(c(1:3), c(3:5))  # argument passed to stat
#' )
#'
#' # and run it (quick run, nb_simul too low and p_acc_min too high for
#' #   accurate results)
#' pars <- aae.pop.inference::inference(
#'   model = pop_inf$model,
#'   prior = pop_inf$prior,
#'   target = do.call(pop_inf$stat, c(list(obs), pop_inf$stat_args)),
#'   nb_simul = 100,
#'   progress_bar = TRUE,
#'   p_acc_min = 0.3
#' )
#'
#' # plot estimated parameters (not expected to be accurate)
#' plot(pars)
define <- function(
    x,
    masks = list(),
    priors = list(),
    stat = aae.pop.inference::stat_abundance_trend,
    sim_args = list(),
    ...
) {

  # check that x is a dynamics object
  if (!aae.pop::is.dynamics(x)) {
    stop(
      "x must be a dynamics object created with aae.pop::dynamics",
      call. = FALSE
    )
  }

  # pull out matrix
  dyn <- x$matrix

  # set default args for simulations
  args <- list(nsim = 100)
  args[names(sim_args)] <- sim_args

  # grab any args for stat
  stat_args <- list(...)

  # set defaults for masks
  mask_set <- list(
    survival = aae.pop::combine(
      aae.pop::survival(dyn),
      aae.pop::transition(dyn)
    ),
    reproduction = aae.pop::reproduction(dyn)
  )
  mask_set[names(masks)] <- masks

  # set defaults for priors
  prior_set <- list(
    survival = c("unif", 0, 1),
    reproduction = c("lognormal", 0, 2)
  )
  prior_set[names(priors)] <- priors

  # set up a parameter list (survival first, then reproduction terms)
  nsurv <- sum(mask_set$survival)
  nrepr <- sum(mask_set$reproduction)
  idx_surv <- seq_len(nsurv)
  idx_repr <- nsurv + seq_len(nrepr)

  # check, expand, and combine priors
  prior <- c(
    expand_prior(
      x = prior_set$survival,
      nx = nsurv,
      type = "survival"
    ),
    expand_prior(
      x = prior_set$reproduction,
      nx = nrepr,
      type = "reproduction"
    )
  )

  # and define the model
  model <- function(pars) {

    # use the model as specified during initial call
    xset <- x

    # update the model matrix
    xset$matrix[mask_set$survival] <- pars[idx_surv]
    xset$matrix[mask_set$reproduction] <- pars[idx_repr]

    # simulate with provided args
    sims <- do.call(simulate, c(list(xset), args))

    # and return summary stat
    do.call(stat, c(list(sims), stat_args))

  }

  # return
  as_definition(
    list(
      model = model,
      prior = prior,
      stat = stat,
      sim_args = args,
      stat_args = stat_args
    )
  )

}

# S3 is method
#' @rdname definition
#' @export
# nolint start
is.definition <- function(x) {
  # nolint end
  inherits(x, "definition")
}

# S3 print method
#' @rdname definition
#' @export
# nolint start
print.definition <- function(x, ...) {
  # nolint end
  cat(
    paste0(
      "Specified model and priors for use with sequential ABC inference ",
      "with aae.pop.inference::inference()\n"
    )
  )
}

# internal function: check and expand priors
expand_prior <- function(x, nx, type) {

  if (is.list(x)) {

    # all good if it's the correct length
    if (length(x) != nx) {
      stop(
        paste0(type, " prior is a list with ", length(x),
        " elements but the ", type, " mask has ", nx, " non-zero elements"),
        call. = FALSE
      )
    }

  } else {

    # otherwise it needs to be a single character vector
    if (!is.character(x)) {
      stop(
        "priors must be character vectors or lists of character vectors",
        call. = FALSE
      )
    }

    # if so, can be expanded to a list with nx elements
    x <- lapply(
      seq_len(nx),
      function(i) x
    )

  }

  # and return
  x

}

# internal function: set definition class
as_definition <- function(x) {
  as_class(x, name = "definition", type = "list")
}
