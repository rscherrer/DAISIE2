# Wrapper around a function to avoid argument collision between environments
call_fun <- function(fun, x, extra = list(), untrans = NULL) {

  # fun: the function to call
  # x: the first argument
  # extra: list with all the other arguments
  # untrans: function to untransform the parameters if needed

  # Check
  testit::assert(is.function(fun))
  testit::assert(is_numeric_vector(x))
  testit::assert(is.list(extra))
  testit::assert(is.null(untrans) | is.function(untrans))

  # Untransform the (re-scaled) parameters if needed
  if (!is.null(untrans)) x <- untrans(x)

  # Call the function
  do.call(fun, c(list(x), extra))

}
