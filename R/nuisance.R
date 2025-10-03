#' Specify a GLM learner (nuisance)
#' @param formula optional; if NULL, we build from vars
#' @param family glm family
#' @param vars character vector of column names
#' @export
fit_glm <- function(formula = NULL, family = gaussian(), vars = NULL) {
  structure(list(kind = "glm", formula = formula, family = family, vars = vars),
            class = "ms_nuisance_spec")
}

#' Specify an xgboost learner (nuisance)
#' @param params list passed to xgboost::xgb.train
#' @param nrounds number of boosting rounds
#' @param vars character vector of column names
#' @export
fit_xgboost <- function(params = list(), nrounds = 200, vars = NULL) {
  structure(list(kind = "xgb", params = params, nrounds = nrounds, vars = vars),
            class = "ms_nuisance_spec")
}

#' Specify a SuperLearner learner (nuisance)
#' @export
fit_superlearner <- function(SL.library, family = gaussian(), vars = NULL) {
  structure(list(kind = "sl", SL.library = SL.library, family = family, vars = vars),
            class = "ms_nuisance_spec")
}

#' Low-level fitter used by msnmm
#' @keywords internal
.ms_fit <- function(spec, X, y, family = NULL) {
  stopifnot(inherits(spec, "ms_nuisance_spec"))
  fam <- spec$family %||% family
  
  if (spec$kind == "glm") {
    if (!is.null(spec$formula)) {
      df <- as.data.frame(X); df$.y <- y
      fit <- stats::glm(spec$formula, data = df, family = fam)
      pred <- function(newX) as.numeric(stats::predict(fit, newdata = as.data.frame(newX), type = if (identical(fam$family, "binomial")) "response" else "link"))
    } else {
      fit <- stats::glm.fit(x = cbind(1, as.matrix(X)), y = y, family = fam)
      pred <- function(newX) {
        eta <- cbind(1, as.matrix(newX)) %*% stats::coef(fit)
        if (identical(fam$family, "binomial")) 1 / (1 + exp(-eta)) else as.numeric(eta)
      }
    }
    return(list(object = fit, predict = pred))
  }
  
  if (spec$kind == "xgb") {
    if (!requireNamespace("xgboost", quietly = TRUE))
      stop("xgboost not installed; install.packages('xgboost')")
    obj <- if (!is.null(fam) && identical(fam$family, "binomial")) "binary:logistic" else "reg:squarederror"
    d <- xgboost::xgb.DMatrix(as.matrix(X), label = y)
    bst <- xgboost::xgb.train(params = utils::modifyList(list(objective = obj), spec$params),
                              data = d, nrounds = spec$nrounds, verbose = 0)
    pred <- function(newX) as.numeric(stats::predict(bst, as.matrix(newX)))
    return(list(object = bst, predict = pred))
  }
  
  if (spec$kind == "sl") {
    if (!requireNamespace("SuperLearner", quietly = TRUE))
      stop("SuperLearner not installed; install.packages('SuperLearner')")
    fit <- SuperLearner::SuperLearner(Y = y, X = as.data.frame(X),
                                      SL.library = spec$SL.library, family = fam)
    pred <- function(newX) as.numeric(stats::predict(fit, newdata = as.data.frame(newX))$pred)
    return(list(object = fit, predict = pred))
  }
  
  stop("Unknown nuisance kind.")
}
