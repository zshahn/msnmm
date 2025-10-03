#' @keywords internal
#' @noRd
`%nin%` <- function(x, table) !match(x, table, nomatch = 0L)
