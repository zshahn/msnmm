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
#' Lightweight GLM wrapper that returns a predictor function
#' @export
fit_glm <- function(formula, data, family = gaussian(), ...) {
  fit <- stats::glm(formula, data = data, family = family, ...)
  list(
    fit = fit,
    predict = function(newdata, type = NULL) {
      # default type: response for both gaussian and binomial
      t <- if (!is.null(type)) type else "response"
      as.numeric(stats::predict(fit, newdata = newdata, type = t))
    }
  )
}




#' Fit with SuperLearner (generic)
#' @param SL.library character vector (e.g. c("SL.ranger","SL.xgboost"))
#' @export
# SuperLearner wrapper: formula -> design matrix, train SL, return predict() closure
fit_superlearner <- function(formula, data, family = gaussian(),
                             learners = c("SL.glm","SL.ranger"), ...) {
  # Build model frame & design matrix
  mf  <- stats::model.frame(formula = formula, data = data, drop.unused.levels = TRUE)
  trm <- stats::terms(mf)
  X   <- stats::model.matrix(trm, mf)
  Y   <- stats::model.response(mf)

  # 1) DROP INTERCEPT (ranger hates "(Intercept)" in data)
  if (ncol(X) && "(Intercept)" %in% colnames(X)) {
    X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  }

  # 2) SANITIZE COLUMN NAMES (make syntactically valid)
  safe_train <- make.names(colnames(X), unique = TRUE, allow_ = TRUE)
  colnames(X) <- safe_train

  # 3) Fit SL on safe, intercept-free X
  sl_fit <- SuperLearner::SuperLearner(
    Y          = Y,
    X          = as.data.frame(X),
    family     = family,
    SL.library = learners
  )

  # 4) Return a predictor that repeats the exact same steps
  list(
    model = sl_fit,
    terms = trm,
    safe_names = safe_train,
    predict = function(newdata, type = NULL) {
      mf_new <- stats::model.frame(trm, data = newdata, drop.unused.levels = TRUE)
      Xnew   <- stats::model.matrix(trm, mf_new)

      # drop intercept at predict time too
      if (ncol(Xnew) && "(Intercept)" %in% colnames(Xnew)) {
        Xnew <- Xnew[, setdiff(colnames(Xnew), "(Intercept)"), drop = FALSE]
      }

      # sanitize & align columns to training cols
      colnames(Xnew) <- make.names(colnames(Xnew), unique = TRUE, allow_ = TRUE)
      missing_cols <- setdiff(safe_train, colnames(Xnew))
      if (length(missing_cols)) {
        Xnew <- cbind(
          Xnew,
          matrix(0, nrow(Xnew), length(missing_cols), dimnames = list(NULL, missing_cols))
        )
      }
      Xnew <- Xnew[, safe_train, drop = FALSE]

      as.numeric(SuperLearner::predict.SuperLearner(
        sl_fit, newdata = as.data.frame(Xnew)
      )$pred)
    }
  )
}




#' Fit xgboost (simple recipe)
#' Internal: XGBoost learner (formula interface)
#' Accepts family either as character ("gaussian"/"binomial") or as a stats::family object.
#' Returns a list with $predict(newdata) -> numeric vector (response scale).
#' @keywords internal
fit_xgboost <- function(formula, data, family = "gaussian",
                        nrounds = 200,
                        params = list(),
                        ...) {

  # ---- normalize family to a character name ----
  fam_name <- if (is.character(family)) {
    tolower(family)
  } else if (is.list(family) && !is.null(family$family)) {
    tolower(family$family)
  } else {
    stop("fit_xgboost(): 'family' must be \"gaussian\"/\"binomial\" or a stats::family object.")
  }
  if (!fam_name %in% c("gaussian","binomial")) {
    stop("fit_xgboost(): unsupported family: ", fam_name)
  }

  # ---- model.frame / design matrices ----
  mf <- stats::model.frame(formula, data = data, drop.unused.levels = TRUE)
  y  <- model.response(mf)
  X  <- stats::model.matrix(stats::delete.response(stats::terms(mf)), data = mf)

  # ---- default params by family ----
  if (fam_name == "gaussian") {
    params_default <- list(objective = "reg:squarederror", eval_metric = "rmse", eta = 0.1, max_depth = 6)
  } else {
    params_default <- list(objective = "binary:logistic",  eval_metric = "logloss", eta = 0.1, max_depth = 6)
  }
  # user params override defaults
  params <- utils::modifyList(params_default, params, keep.null = TRUE)

  # ---- train ----
  dtr <- xgboost::xgb.DMatrix(data = X, label = y)
  bst <- xgboost::xgb.train(
    params  = params,
    data    = dtr,
    nrounds = nrounds,
    verbose = 0
  )

  # ---- return predict wrapper on response scale ----
  pred_fun <- function(newdata) {
    mf_new <- stats::model.frame(stats::terms(mf), data = newdata, drop.unused.levels = TRUE)
    X_new  <- stats::model.matrix(stats::delete.response(stats::terms(mf)), data = mf_new)
    pr     <- stats::predict(bst, newdata = X_new)
    # gaussian: already on response; binomial: already probability due to objective
    as.numeric(pr)
  }

  list(predict = pred_fun)
}

# Predictors for wrappers
#' @export
predict.xgb_fit <- function(object, newdata, type = "response", ...) {
  X <- model.matrix(delete.response(terms( ~ . )), newdata[, object$x_cols, drop = FALSE])
  p <- stats::predict(object$model, newdata = X)
  if (identical(object$family$family, "binomial") && type != "link") return(p)
  p
}



