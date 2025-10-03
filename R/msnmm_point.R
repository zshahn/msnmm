#' Point-exposure marginal SNMM with cross-fitting
#' @param data data.frame with columns y0,y1,a and specified covariates
#' @param y0,y1,a column names
#' @param L0 adjustment set at baseline
#' @param Ltilde0 effect-modifier set (for the blip)
#' @param blip formula for γ(ltilde; ψ), linear in ψ
#' @param q "score" (default) or a function(data)-> n x p matrix or a formula
#' @param nuisance list(pi=..., mu=...) each created with fit_glm/xgboost/SuperLearner
#' @param folds integer >=2
#' @param trim upper-trim on propensity (1 - π); 0 means none
#' @export
msnmm_point <- function(data, y0, y1, a, L0, Ltilde0,
                        blip = ~ 1,
                        q = "score",
                        nuisance = list(pi = fit_glm(vars = L0, family = binomial()),
                                        mu = fit_glm(vars = L0, family = gaussian())),
                        folds = 2, trim = 0, solver = c("nleqslv","oneshot")[1],
                        control = list(), seed = NULL) {

  dat <- data
  .ms_stopifnot_cols(dat, c(y0, y1, a))
  .ms_stopifnot_cols(dat, L0, "L0")
  .ms_stopifnot_cols(dat, Ltilde0, "Ltilde0")

  D <- dat[[y1]] - dat[[y0]]
  A <- dat[[a]]
  if (!.ms_is_binary(A)) stop("Treatment column '", a, "' must be binary 0/1.")

  Qfun <- if (identical(q, "score")) function(d) model.matrix(blip, d[, Ltilde0, drop = FALSE])
  else if (inherits(q, "function")) q
  else if (inherits(q, "formula")) make_q(q)
  else stop("q must be 'score', a function, or a formula.")

  n <- nrow(dat)
  fold_id <- .ms_make_folds(n, folds, seed = seed)
  requireNamespace("nleqslv")

  psi_list <- IF_list <- D_list <- vector("list", folds)

  for (s in seq_len(folds)) {
    idx_tr <- which(fold_id != s)
    idx_te <- which(fold_id == s)

    Xpi <- dat[idx_tr, L0, drop = FALSE]
    ypi <- A[idx_tr]
    pi_fit <- .ms_fit(nuisance$pi, X = Xpi, y = ypi, family = binomial())

    Xmu <- dat[idx_tr, L0, drop = FALSE]
    rows_mu <- idx_tr[which(A[idx_tr] == 0)]
    mu_fit <- .ms_fit(nuisance$mu, X = Xmu[A[idx_tr] == 0, , drop = FALSE],
                      y = D[rows_mu], family = gaussian())

    # predict on test
    pi_hat <- .ms_clip(pi_fit$predict(dat[idx_te, L0, drop = FALSE]))
    if (trim > 0) pi_hat <- pmin(pi_hat, 1 - trim)
    mu_hat <- mu_fit$predict(dat[idx_te, L0, drop = FALSE])

    Q <- Qfun(dat[idx_te, , drop = FALSE])
    p <- ncol(Q)

    estfun <- function(psi) {
      gamma <- as.numeric(Q %*% psi)
      eps <- ((A[idx_te] - pi_hat) / (1 - pi_hat)) * (D[idx_te] - mu_hat) - A[idx_te] * gamma
      as.numeric(crossprod(Q, eps))
    }

    if (solver == "nleqslv") {
      sol <- nleqslv::nleqslv(x = rep(0, p), fn = estfun, control = control)
      psi_hat <- sol$x
    } else psi_hat <- rep(0, p)

    gamma <- as.numeric(Q %*% psi_hat)
    eps <- ((A[idx_te] - pi_hat) / (1 - pi_hat)) * (D[idx_te] - mu_hat) - A[idx_te] * gamma

    Dhat <- - crossprod(Q, A[idx_te] * Q) / length(idx_te)
    IF <- sweep(Q, 1, eps, `*`)

    psi_list[[s]] <- psi_hat
    IF_list[[s]] <- IF
    D_list[[s]]  <- Dhat
  }

  psi_hat <- Reduce(`+`, psi_list) / folds
  Dbar <- Reduce(`+`, D_list) / folds
  IF_all <- do.call(rbind, IF_list)
  Sigma <- solve(Dbar) %*% stats::cov(IF_all) %*% t(solve(Dbar)) / n

  out <- list(
    call = match.call(),
    coef = setNames(psi_hat, colnames(Qfun(dat))),
    vcov = Sigma,
    folds = folds,
    IF = IF_all,
    D = Dbar,
    blip = blip,
    Ltilde0 = Ltilde0
  )
  class(out) <- c("msnmm_point", "msnmm_fit")
  out
}
