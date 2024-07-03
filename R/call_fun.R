# Wrapper around a function to avoid argument collision between environments
call_fun <- function(fun, x, extra = list()) {

  # fun: the function to call
  # x: the first argument
  # extra: list with all the other arguments

  # Call the function
  do.call(fun, c(list(x), extra))

}
