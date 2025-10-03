#' @export
summary.msnmm_pairs_fit <- function(object, ...) {
  est <- object$coef
  if (is.null(est)) stop("No coefficients found in object$coef.")

  V <- object$vcov
  if (is.null(V) || !is.matrix(V) || anyNA(diag(V))) {
    out <- cbind(Estimate = est)
  } else {
    se <- sqrt(pmax(diag(V), 0))
    z  <- ifelse(se > 0, est / se, NA_real_)
    p  <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    out <- cbind(Estimate = est, SE = se, z = z, `Pr(>|z|)` = p)
  }

  class(out) <- c("summary.msnmm_pairs_fit", "matrix")
  out
}

#' @export
print.summary.msnmm_pairs_fit <- function(x, digits = 4, ...) {
  stats::printCoefmat(x, digits = digits,
                      has.Pvalue = "Pr(>|z|)" %in% colnames(x))
  invisible(x)
}

# Optional: if you also use class "msnmm_fit", point its summary to the same code
#' @export
summary.msnmm_fit <- function(object, ...) {
  class(object) <- "msnmm_pairs_fit"
  summary(object, ...)
}
