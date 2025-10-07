# --- tiny xgboost wrappers with a uniform predict() interface ---------------

ms_xgb_bin <- function(vars, nrounds = 200,
                       params = list(objective = "binary:logistic",
                                     eval_metric = "logloss",
                                     max_depth = 4, eta = 0.3, nthread = 1)) {
  force(vars); force(nrounds); force(params)
  structure(list(
    type = "bin",
    vars = vars,
    fit = function(data, y, rows = rep(TRUE, nrow(data))) {
      X <- as.matrix(data[rows, vars, drop = FALSE])
      yy <- as.numeric(y[rows])
      mdl <- xgboost::xgboost(params = params,
                              data = xgboost::xgb.DMatrix(data = X, label = yy),
                              nrounds = nrounds, verbose = 0)
      list(
        predict = function(newdata) {
          as.numeric(predict(mdl, as.matrix(newdata[, vars, drop = FALSE])))
        }
      )
    }
  ), class = "ms_learner")
}

ms_xgb_reg <- function(vars, nrounds = 300,
                       params = list(objective = "reg:squarederror",
                                     max_depth = 4, eta = 0.3, nthread = 1)) {
  force(vars); force(nrounds); force(params)
  structure(list(
    type = "reg",
    vars = vars,
    fit = function(data, y, rows = rep(TRUE, nrow(data))) {
      X <- as.matrix(data[rows, vars, drop = FALSE])
      yy <- as.numeric(y[rows])
      mdl <- xgboost::xgboost(params = params,
                              data = xgboost::xgb.DMatrix(data = X, label = yy),
                              nrounds = nrounds, verbose = 0)
      list(
        predict = function(newdata) {
          as.numeric(predict(mdl, as.matrix(newdata[, vars, drop = FALSE])))
        }
      )
    }
  ), class = "ms_learner")
}

#' Fit a GLM (wrapper)
#' @param formula model formula
#' @param data    data.frame
#' @param family  a GLM family, e.g. stats::binomial()
#' @return a fitted "glm" object
#' @export
fit_glm <- function(formula, data, family) {
  stats::glm(formula = formula, data = data, family = family)
}


#' Fit with SuperLearner (generic)
#' @param SL.library character vector (e.g. c("SL.ranger","SL.xgboost"))
#' @export
fit_superlearner <- function(formula, data, family, SL.library, ...) {
  stopifnot(requireNamespace("SuperLearner", quietly = TRUE))
  mf <- stats::model.frame(formula, data = data)
  y  <- stats::model.response(mf)
  x  <- mf[, setdiff(colnames(mf), all.vars(stats::terms(stats::update(formula, . ~ . -1)))) , drop = FALSE] # robust-ish
  x  <- mf[, setdiff(colnames(mf), as.character(formula[[2]])), drop = FALSE]
  SuperLearner::SuperLearner(Y = y, X = x, SL.library = SL.library,
                             family = if (identical(family$family, "binomial")) stats::binomial() else gaussian(),
                             ...)
}

#' Fit xgboost (simple recipe)
#' @export
fit_xgboost <- function(formula, data, family, nrounds = 100, params = list(), ...) {
  stopifnot(requireNamespace("xgboost", quietly = TRUE))
  mf <- stats::model.frame(formula, data = data)
  y  <- stats::model.response(mf)
  X  <- model.matrix(stats::terms(formula), data = data)
  if (identical(family$family, "binomial")) {
    params <- utils::modifyList(list(objective = "binary:logistic", eval_metric = "logloss"), params)
  } else {
    params <- utils::modifyList(list(objective = "reg:squarederror"), params)
  }
  m <- xgboost::xgboost(data = xgboost::xgb.DMatrix(X, label = y),
                        params = params, nrounds = nrounds, verbose = 0)
  structure(list(model = m, x_cols = colnames(X), family = family),
            class = "xgb_fit")
}

# Predictors for wrappers
#' @export
predict.xgb_fit <- function(object, newdata, type = "response", ...) {
  X <- model.matrix(delete.response(terms( ~ . )), newdata[, object$x_cols, drop = FALSE])
  p <- stats::predict(object$model, newdata = X)
  if (identical(object$family$family, "binomial") && type != "link") return(p)
  p
}

#' @export
predict.SuperLearner <- function(object, newdata, type = "response", ...) {
  pr <- stats::predict(object, newdata = newdata)$pred
  if (!is.null(dim(pr))) pr <- as.numeric(pr)
  pr
}

