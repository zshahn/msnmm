#' Build a q-function from a formula (returns model matrix)
#' @export
make_q <- function(formula) {
  force(formula)
  function(data) model.matrix(formula, data)
}
