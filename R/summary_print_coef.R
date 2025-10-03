# R/summary-print-coef.R

# minimal summary for msnmm_pairs_fit
summary.msnmm_pairs_fit <- function(object, ...) {
  est <- object$coef
  V   <- object$vcov
  se  <- if (is.null(V)) rep(NA_real_, length(est)) else sqrt(diag(V))
  z   <- est / se
  p   <- 2 * stats::pnorm(-abs(z))
  out <- cbind(Estimate = est, SE = se, z = z, `Pr(>|z|)` = p)
  class(out) <- c("summary.msnmm_pairs_fit", "matrix")
  out
}

print.summary.msnmm_pairs_fit <- function(x, ...) {
  printCoefmat(x, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}

