#' Make V-fold indices
#' @export
make_folds <- function(n, V = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  split(sample.int(n), rep(1:V, length.out = n))
}

# Generic CF fit/predict for a single regression
.cf_fit_predict <- function(formula, data, family, V, fit_fun, pred_fun = predict, ...) {
  n <- nrow(data); out <- numeric(n)
  folds <- make_folds(n, V = V)
  for (v in seq_along(folds)) {
    idx_te <- folds[[v]]
    idx_tr <- setdiff(seq_len(n), idx_te)
    fit <- do.call(fit_fun, c(list(formula = formula, data = data[idx_tr, , drop=FALSE], family = family), list(...)))
    out[idx_te] <- do.call(pred_fun, list(object = fit, newdata = data[idx_te, , drop=FALSE]))
  }
  out
}
