coef.msnmm_fit <- function(object, ...) object$coef

vcov.msnmm_fit <- function(object, ...) object$vcov

summary.msnmm_fit <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  z  <- object$coef / se
  p  <- 2 * stats::pnorm(-abs(z))
  out <- cbind(Estimate = object$coef, SE = se, z = z, `Pr(>|z|)` = p)
  class(out) <- c("summary.msnmm_fit", "matrix")
  out
}
