#' Coefficients for msnmm_pairs_fit
#' @param object An msnmm_pairs_fit object.
#' @export
coef.msnmm_pairs_fit <- function(object, ...) {
  object$coef
}

#' Variance-covariance for msnmm_pairs_fit
#' @export
vcov.msnmm_pairs_fit <- function(object, ...) {
  object$vcov
}

#' Summarize msnmm_pairs_fit
#'
#' Produces a matrix with Estimate, SE, z, and p-values.
#' @export
summary.msnmm_pairs_fit <- function(object, ...) {
  beta <- object$coef
  V    <- object$vcov
  se   <- suppressWarnings(sqrt(diag(V)))
  z    <- beta / se
  p    <- 2 * stats::pnorm(-abs(z))
  tab  <- cbind(Estimate = beta, SE = se, z = z, `Pr(>|z|)` = p)
  class(tab) <- c("summary.msnmm_pairs_fit", class(tab))
  attr(tab, "call") <- object$call
  tab
}

#' @export
print.summary.msnmm_pairs_fit <- function(x, digits = 4, ...) {
  cat("Call:\n")
  print(attr(x, "call"))
  cat("\nCoefficients:\n")
  printCoefmat(x, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}
