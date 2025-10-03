#' Two-time (K=2) marginal SNMM with cross-fitting
#' @export
msnmm_tv2 <- function(data,
                      y = c("Y0","Y1","Y2"),
                      a = c("A0","A1"),
                      L = list(L0 = NULL, L1 = NULL),
                      Ltilde = list(L0 = NULL, L1 = NULL),
                      blip = list(g01 = ~ 1, g12 = ~ 1, g02 = ~ 1),
                      q = "score",
                      nuisance = list(
                        pi0  = fit_glm(vars = L$L0, family = binomial()),
                        pi1  = fit_glm(vars = c(L$L1, a[1]), family = binomial()),
                        mu01 = fit_glm(vars = L$L0, family = gaussian()),
                        mu12 = fit_glm(vars = c(L$L1, a[1]), family = gaussian())
                      ),
                      folds = 2, trim = 0, solver = "nleqslv", control = list(), seed = NULL) {
  
  dat <- data
  Y0 <- dat[[y[1]]]; Y1 <- dat[[y[2]]]; Y2 <- dat[[y[3]]]
  A0 <- dat[[a[1]]];  A1 <- dat[[a[2]]]
  if (!.ms_is_binary(A0) || !.ms_is_binary(A1)) stop("A0/A1 must be 0/1.")
  D01 <- Y1 - Y0
  D12 <- Y2 - Y1
  
  MM01 <- model.matrix(blip$g01, dat[, Ltilde$L0, drop = FALSE])
  MM12 <- model.matrix(blip$g12, dat[, Ltilde$L1, drop = FALSE])
  MM02 <- model.matrix(blip$g02, dat[, Ltilde$L0, drop = FALSE])
  p01 <- ncol(MM01); p12 <- ncol(MM12); p02 <- ncol(MM02); p <- p01 + p12 + p02
  
  fold_id <- .ms_make_folds(nrow(dat), folds, seed = seed)
  requireNamespace("nleqslv")
  
  psi_list <- IF_list <- D_list <- vector("list", folds)
  
  for (s in seq_len(folds)) {
    idx_tr <- which(fold_id != s)
    idx_te <- which(fold_id == s)
    
    pi0_fit <- .ms_fit(nuisance$pi0,  X = dat[idx_tr, L$L0, drop = FALSE], y = A0[idx_tr], family = binomial())
    pi1_fit <- .ms_fit(nuisance$pi1,  X = dat[idx_tr, c(L$L1, a[1]), drop = FALSE], y = A1[idx_tr], family = binomial())
    
    mu01_fit <- .ms_fit(nuisance$mu01,
                        X = dat[idx_tr, L$L0, drop = FALSE][A0[idx_tr] == 0, , drop = FALSE],
                        y = D01[idx_tr][A0[idx_tr] == 0], family = gaussian())
    mu12_fit <- .ms_fit(nuisance$mu12,
                        X = dat[idx_tr, c(L$L1, a[1]), drop = FALSE][A1[idx_tr] == 0, , drop = FALSE],
                        y = D12[idx_tr][A1[idx_tr] == 0], family = gaussian())
    
    pi0_hat <- .ms_clip(pi0_fit$predict(dat[idx_te, L$L0, drop = FALSE]))
    pi1_hat <- .ms_clip(pi1_fit$predict(dat[idx_te, c(L$L1, a[1]), drop = FALSE]))
    if (trim > 0) { pi0_hat <- pmin(pi0_hat, 1 - trim); pi1_hat <- pmin(pi1_hat, 1 - trim) }
    
    mu01_hat <- mu01_fit$predict(dat[idx_te, L$L0, drop = FALSE])
    mu12_hat <- mu12_fit$predict(dat[idx_te, c(L$L1, a[1]), drop = FALSE])
    
    # mu02: regress predicted mu12 on L0 among A0=0 (train), then predict on test
    mu02_reg <- .ms_fit(fit_glm(vars = L$L0, family = gaussian()),
                        X = dat[idx_tr, L$L0, drop = FALSE][A0[idx_tr] == 0, , drop = FALSE],
                        y = mu12_fit$predict(dat[idx_tr, c(L$L1, a[1]), drop = FALSE][A0[idx_tr] == 0, , drop = FALSE]),
                        family = gaussian())
    mu02_hat <- mu02_reg$predict(dat[idx_te, L$L0, drop = FALSE])
    
    Q01 <- MM01[idx_te, , drop = FALSE]
    Q12 <- MM12[idx_te, , drop = FALSE]
    Q02 <- MM02[idx_te, , drop = FALSE]
    
    estfun <- function(psi) {
      psi01 <- psi[seq_len(p01)]
      psi12 <- psi[p01 + seq_len(p12)]
      psi02 <- psi[p01 + p12 + seq_len(p02)]
      
      g01 <- as.numeric(Q01 %*% psi01)
      g12 <- as.numeric(Q12 %*% psi12)
      g02 <- as.numeric(Q02 %*% psi02)
      
      eps01 <- ((A0[idx_te] - pi0_hat) / (1 - pi0_hat)) * (D01[idx_te] - mu01_hat) - A0[idx_te] * g01
      eps12 <- ((A1[idx_te] - pi1_hat) / (1 - pi1_hat)) * (D12[idx_te] - mu12_hat) - A1[idx_te] * g12
      
      term1 <- (1 - (1 - A1[idx_te]) * (1 - A0[idx_te]) / ((1 - pi1_hat) * (1 - pi0_hat))) * D12[idx_te]
      term2 <- - ((1 - A0[idx_te]) / (1 - pi0_hat)) * (1 - ((1 - A1[idx_te]) / (1 - pi1_hat))) * mu12_hat
      term3 <- - (1 - (1 - A0[idx_te]) / (1 - pi0_hat)) * mu02_hat
      eps02 <- term1 + term2 + term3 - A1[idx_te] * g12 - A0[idx_te] * g02 + A0[idx_te] * g01
      
      c(colSums(Q01 * eps01), colSums(Q12 * eps12), colSums(Q02 * eps02))
    }
    
    sol <- nleqslv::nleqslv(x = rep(0, p), fn = estfun, control = control)
    psi_hat <- sol$x
    
    g01 <- as.numeric(Q01 %*% psi_hat[seq_len(p01)])
    g12 <- as.numeric(Q12 %*% psi_hat[p01 + seq_len(p12)])
    g02 <- as.numeric(Q02 %*% psi_hat[p01 + p12 + seq_len(p02)])
    
    eps01 <- ((A0[idx_te] - pi0_hat) / (1 - pi0_hat)) * (D01[idx_te] - mu01_hat) - A0[idx_te] * g01
    eps12 <- ((A1[idx_te] - pi1_hat) / (1 - pi1_hat)) * (D12[idx_te] - mu12_hat) - A1[idx_te] * g12
    term1 <- (1 - (1 - A1[idx_te]) * (1 - A0[idx_te]) / ((1 - pi1_hat) * (1 - pi0_hat))) * D12[idx_te]
    term2 <- - ((1 - A0[idx_te]) / (1 - pi0_hat)) * (1 - ((1 - A1[idx_te]) / (1 - pi1_hat))) * mu12_hat
    term3 <- - (1 - (1 - A0[idx_te]) / (1 - pi0_hat)) * mu02_hat
    eps02 <- term1 + term2 + term3 - A1[idx_te] * g12 - A0[idx_te] * g02 + A0[idx_te] * g01
    
    Dhat <- matrix(0, p, p)
    Dhat[1:p01, 1:p01] <- - crossprod(Q01, A0[idx_te] * Q01) / length(idx_te)
    Dhat[p01 + (1:p12), p01 + (1:p12)] <- - crossprod(Q12, A1[idx_te] * Q12) / length(idx_te)
    Dhat[p01 + p12 + (1:p02), p01 + p12 + (1:p02)] <- - crossprod(Q02, A0[idx_te] * Q02) / length(idx_te)
    
    IF <- cbind(Q01 * eps01, Q12 * eps12, Q02 * eps02)
    
    psi_list[[s]] <- psi_hat
    IF_list[[s]]  <- IF
    D_list[[s]]   <- Dhat
  }
  
  psi_hat <- Reduce(`+`, psi_list) / folds
  Dbar <- Reduce(`+`, D_list) / folds
  IF_all <- do.call(rbind, IF_list)
  Sigma <- solve(Dbar) %*% stats::cov(IF_all) %*% t(solve(Dbar)) / nrow(dat)
  
  nm <- c(paste0("psi01_", colnames(MM01)),
          paste0("psi12_", colnames(MM12)),
          paste0("psi02_", colnames(MM02)))
  
  out <- list(coef = setNames(psi_hat, nm), vcov = Sigma, folds = folds)
  class(out) <- c("msnmm_tv", "msnmm_fit")
  out
}
