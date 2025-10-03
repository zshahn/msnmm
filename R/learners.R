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
