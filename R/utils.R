# set an object class
as_class <- function(
    object,
    name,
    type = c("function", "list", "matrix", "array", "dynamics")
) {

  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))

  object

}
