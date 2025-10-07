# R/zzz-stubs.R

#' @export
fit_glm <- function(...) stop("fit_glm() not yet implemented in this build.")

#' @export
fit_superlearner <- function(...) stop("fit_superlearner() not yet implemented in this build.")

#' @export
fit_xgboost <- function(...) stop("fit_xgboost() not yet implemented in this build.")

#' @export
create_wide_format <- function(...) stop("create_wide_format() not yet implemented here (use msnmm_tv_pairs_min workflow).")

#' @export
make_diffs <- function(...) stop("make_diffs() not yet implemented in this build.")

#' @export
make_pastA <- function(...) stop("make_pastA() not yet implemented in this build.")

#' @export
make_q <- function(...) stop("make_q() not yet implemented in this build.")

#' @export
msnmm_point <- function(...) stop("msnmm_point() not yet included in this build.")

#' @export
msnmm_tv2 <- function(...) stop("msnmm_tv2() not yet included in this build.")

# If msnmm_fit isn’t ready, either remove its export or stub it:
#' @export
msnmm_fit <- function(...) stop("msnmm_fit() not yet included in this build (use msnmm_tv_pairs_min).")
