#' Minimal unified (m,k) SNMM with optional ML + cross-fitting
#'
#' @param data,id,time,ntimes,... same as before
#' @param model one of "glm","xgb","sl" (SuperLearner)
#' @param cf logical; if TRUE, K-fold cross-fitting for all nuisances
#' @param folds integer; number of folds (default 2)
#' @param seed optional integer for reproducible folds
#' @param learners_outcome,learners_treat used when model="sl"
#'
#' @export
msnmm_tv_pairs_min <- function(
    data,
    id,
    time = "yearms",
    ntimes,
    exposure = "A", outcome = "Y",
    time_vars_to_widen = c("A","Y"),
    wide_ready = FALSE,
    outcome_nuisance_formulas,
    treatment_nuisance_formulas,
    blips,
    initiation = FALSE,
    verbose = TRUE,
    model = c("glm","xgb","sl"),
    cf = FALSE,
    folds = 2,
    seed = NULL,
    learners_outcome = c("SL.glm","SL.ranger"),
    learners_treat   = c("SL.glm","SL.ranger")
){
  model <- match.arg(model)
  dat <- data

  # ---- helpers (use your exported helpers when available) ----------------
  if (!all(paste0("diffs", 1:ntimes) %in% names(dat))) {
    dat <- make_diffs(dat, outcome = outcome, K = ntimes)   # exported helper
  }
  if (!all(paste0("past", exposure, "_", 2:ntimes) %in% names(dat))) {
    dat <- make_pastA(dat, exposure = exposure, K = ntimes) # exported helper
  }
  if (!isTRUE(wide_ready)) {
    dat <- create_wide_format(
      dat, id = id,
      time_vars = c(time_vars_to_widen, setdiff(colnames(dat), c(id, time_vars_to_widen, time))),
      max_timepoints = ntimes
    )
    dat <- dat[dat[[time]] == 1, , drop = FALSE]
  }

  # if "first treatment only": set future A_i = 0 (data coding only)
  if (isTRUE(initiation)) {
    for (i in 2:ntimes) {
      hit <- isTRUE(dat[[paste0("past", exposure, "_", i)]])
      dat[hit, paste0(exposure, "_", i)] <- 0
    }
  }

  # convenience lists
  Y <- lapply(0:ntimes, function(t) dat[[paste0(outcome, "_", t)]])
  A <- lapply(1:ntimes, function(t) dat[[paste0(exposure, "_", t)]])

  # fold assignment
  if (isTRUE(cf)) {
    if (!is.null(seed)) set.seed(seed)
    fold_id <- sample(rep(1:folds, length.out = nrow(dat)))
  } else {
    fold_id <- rep(1, nrow(dat))
    folds <- 1
  }

  # ------------------ 2) Outcome nuisances: preds{i}_{k} ------------------
  # seed diagonal like your script
  for (i in 1:ntimes) {
    dat[[paste0("preds", i+1, "_", i)]] <- dat[[paste0("diffs", i)]]
  }

  fit_and_predict_outcome <- function(formula, i, k, te_rows) {
    # fit only on (A_i==0) AND (no past A for initiation settings when relevant)
    base_fit <- dat[[paste0(exposure, "_", i)]] == 0
    # (no extra filter here; that's already the structural conditioning)
    df_fit <- dat[base_fit, , drop = FALSE]
    if (nrow(df_fit) == 0L) {
      stop("No rows to fit outcome nuisance at (i=", i, ", k=", k, ").")
    }

    fam <- gaussian()

    # choose backend
    if (!isTRUE(cf)) {
      if (model == "glm") {
        mod <- fit_glm(formula, data = df_fit, family = fam)
        pred <- mod$predict(dat)
      } else if (model == "xgb") {
        mod <- fit_xgboost(formula, data = df_fit, family = "gaussian")
        pred <- mod$predict(dat)
      } else { # "sl"
        mod <- fit_superlearner(formula, data = df_fit, family = fam,
                                learners = learners_outcome)
        pred <- mod$predict(dat)
      }
      dat[[paste0("preds", i, "_", k)]] <<- as.numeric(pred)
      return(invisible(NULL))
    }

    # cross-fitting branch
    pred_all <- numeric(nrow(dat))
    for (s in 1:folds) {
      tr <- fold_id != s
      te <- fold_id == s

      # training rows must also satisfy A_i==0
      tr_rows <- which(tr & (dat[[paste0(exposure, "_", i)]] == 0))
      df_tr <- dat[tr_rows, , drop = FALSE]
      if (nrow(df_tr) == 0L) next

      if (model == "glm") {
        mod <- fit_glm(formula, data = df_tr, family = fam)
      } else if (model == "xgb") {
        mod <- fit_xgboost(formula, data = df_tr, family = "gaussian")
      } else {
        mod <- fit_superlearner(formula, data = df_tr, family = fam,
                                learners = learners_outcome)
      }
      pred_all[te] <- as.numeric(mod$predict(dat[te, , drop = FALSE]))
    }
    dat[[paste0("preds", i, "_", k)]] <<- pred_all
    invisible(NULL)
  }

  for (k in ntimes:1) {
    for (i in k:1) {
      form <- outcome_nuisance_formulas[[i]][[k - i + 1]]
      fit_and_predict_outcome(form, i, k, te_rows = NULL)
    }
  }

  # ------------------ 3) Treatment nuisances: A_i -------------------------
  fit_and_predict_treat <- function(formula, i) {
    fam <- binomial()
    # for initiation analyses, fit on "no past treatment" rows
    if (isTRUE(initiation) && i > 1) {
      fit_rows <- which(!dat[[paste0("past", exposure, "_", i)]])
    } else {
      fit_rows <- seq_len(nrow(dat))
    }

    if (!isTRUE(cf)) {
      if (model == "glm") {
        mod <- fit_glm(formula, data = dat[fit_rows, , drop = FALSE], family = fam)
        pred <- mod$predict(dat, type = "response")
      } else if (model == "xgb") {
        mod <- fit_xgboost(formula, data = dat[fit_rows, , drop = FALSE], family = "binomial")
        pred <- mod$predict(dat)
      } else {
        mod <- fit_superlearner(formula, data = dat[fit_rows, , drop = FALSE], family = fam,
                                learners = learners_treat)
        pred <- mod$predict(dat)
      }
      dat[[paste0(exposure, "_", i, "_hat")]] <<- as.numeric(pred)
      return(invisible(NULL))
    }

    # cross-fitting branch
    pred_all <- numeric(nrow(dat))
    for (s in 1:folds) {
      tr <- (fold_id != s)
      te <- (fold_id == s)

      tr_rows <- intersect(which(tr), fit_rows)
      df_tr <- dat[tr_rows, , drop = FALSE]
      if (nrow(df_tr) == 0L) next

      if (model == "glm") {
        mod <- fit_glm(formula, data = df_tr, family = fam)
        pred_all[te] <- as.numeric(mod$predict(dat[te, , drop = FALSE], type = "response"))
      } else if (model == "xgb") {
        mod <- fit_xgboost(formula, data = df_tr, family = "binomial")
        pred_all[te] <- as.numeric(mod$predict(dat[te, , drop = FALSE]))
      } else {
        mod <- fit_superlearner(formula, data = df_tr, family = fam,
                                learners = learners_treat)
        pred_all[te] <- as.numeric(mod$predict(dat[te, , drop = FALSE]))
      }
    }
    dat[[paste0(exposure, "_", i, "_hat")]] <<- pred_all
    invisible(NULL)
  }

  for (i in 1:ntimes) {
    fit_and_predict_treat(treatment_nuisance_formulas[[i]], i)
  }

  # if initiation: zero out predicted treatment after first treatment
  if (isTRUE(initiation)) {
    for (i in 2:ntimes) {
      hit <- isTRUE(dat[[paste0("past", exposure, "_", i)]])
      dat[hit, paste0(exposure, "_", i, "_hat")] <- 0
    }
  }

  # ------------------ 4) Weights and term2 -------------------------------
  for (i in 1:ntimes) {
    dat[[paste0("weightfactor_", i)]] <-
      (1 - dat[[paste0(exposure, "_", i)]]) / (1 - dat[[paste0(exposure, "_", i, "_hat")]])
  }
  for (m in 1:ntimes) {
    for (k in m:ntimes) {
      cols <- paste0("weightfactor_", m:k)
      dat[[paste0("weightprod_", m, k)]] <- apply(as.matrix(dat[, cols, drop = FALSE]), 1, prod)
    }
  }
  for (m in 1:ntimes) {
    for (k in m:ntimes) {
      for (j in m:k) {
        dat[[paste0("term2summands_", m, k, "_", j)]] <-
          dat[[paste0("weightprod_", m, j)]] *
          (dat[[paste0("preds", j+1, "_", k)]] - dat[[paste0("preds", j, "_", k)]])
      }
      cols <- paste0("term2summands_", m, k, "_", m:k)
      dat[[paste0("term2_", m, k)]] <- rowSums(as.matrix(dat[, cols, drop = FALSE]))
    }
  }

  # ------------------ 5) Drop inactive m, build Q ------------------------
  if (!("intercept" %in% names(dat))) dat$intercept <- 1
  active_m <- which(vapply(1:ntimes, function(m) any(A[[m]] == 1, na.rm = TRUE), FALSE))
  keep <- vapply(blips, function(b) b$m %in% active_m, FALSE)
  blips <- blips[keep]
  if (!length(blips)) stop("No blip blocks remain after dropping inactive m.")

  for (blk in blips) {
    miss <- setdiff(blk$preds, names(dat))
    if (length(miss)) stop("Missing columns in blip preds: ", paste(miss, collapse = ", "))
  }
  Q_list  <- vector("list", length(blips))
  p_block <- integer(length(blips))
  for (i in seq_along(blips)) {
    Q_list[[i]] <- as.matrix(dat[, blips[[i]]$preds, drop = FALSE])
    p_block[i]  <- ncol(Q_list[[i]])
  }
  p_total <- sum(p_block)

  # ------------------ 6) G-estimation system ----------------------------
  Yt <- function(t) Y[[t + 1L]]

  block_offsets <- c(0L, cumsum(p_block))[seq_along(p_block)]
  block_slice <- function(i) {
    off <- block_offsets[i]; (off + 1L):(off + p_block[i])
  }

  est_eq <- function(psi) {
    # blip contributions B_{j->k}
    B_contrib <- list()
    for (i in seq_along(blips)) {
      m_i <- blips[[i]]$m
      k_i <- blips[[i]]$k
      g   <- as.numeric(Q_list[[i]] %*% psi[block_slice(i)])
      key <- paste0(m_i, "_", k_i)
      B_contrib[[key]] <- A[[m_i]] * g
    }

    mom <- numeric(p_total)
    for (i in seq_along(blips)) {
      m_i <- blips[[i]]$m; k_i <- blips[[i]]$k

      sum_blips_mk <- 0
      for (j in m_i:k_i) {
        key <- paste0(j, "_", k_i)
        if (!is.null(B_contrib[[key]])) sum_blips_mk <- sum_blips_mk + B_contrib[[key]]
      }
      H_mk <- Yt(k_i) - sum_blips_mk

      if (k_i == m_i) {
        H_mk_1 <- Yt(m_i - 1)
      } else {
        sum_blips_mk1 <- 0
        for (j in m_i:(k_i - 1)) {
          key <- paste0(j, "_", k_i - 1)
          if (!is.null(B_contrib[[key]])) sum_blips_mk1 <- sum_blips_mk1 + B_contrib[[key]]
        }
        H_mk_1 <- Yt(k_i - 1) - sum_blips_mk1
      }

      term1 <- (H_mk - H_mk_1) - dat[[paste0("preds", m_i, "_", k_i)]]
      term2 <- dat[[paste0("term2_", m_i, k_i)]]
      eps   <- as.numeric(term1 - term2)

      mom[block_slice(i)] <- as.numeric(crossprod(Q_list[[i]], eps))
    }
    mom
  }

  sol  <- nleqslv::nleqslv(x = rep(0, p_total), fn = est_eq)
  psi_hat <- sol$x

  # crude IF covariance
  IF_blocks <- vector("list", length(blips))
  for (i in seq_along(blips)) {
    # recompute eps at psi_hat
    B_contrib <- list()
    for (ii in seq_along(blips)) {
      g   <- as.numeric(Q_list[[ii]] %*% psi_hat[block_slice(ii)])
      key <- paste0(blips[[ii]]$m, "_", blips[[ii]]$k)
      B_contrib[[key]] <- A[[blips[[ii]]$m]] * g
    }
    m_i <- blips[[i]]$m; k_i <- blips[[i]]$k
    sum_blips_mk <- 0
    for (j in m_i:k_i) {
      key <- paste0(j, "_", k_i)
      if (!is.null(B_contrib[[key]])) sum_blips_mk <- sum_blips_mk + B_contrib[[key]]
    }
    H_mk <- Yt(k_i) - sum_blips_mk
    if (k_i == m_i) {
      H_mk_1 <- Yt(m_i - 1)
    } else {
      sum_blips_mk1 <- 0
      for (j in m_i:(k_i - 1)) {
        key <- paste0(j, "_", k_i - 1)
        if (!is.null(B_contrib[[key]])) sum_blips_mk1 <- sum_blips_mk1 + B_contrib[[key]]
      }
      H_mk_1 <- Yt(k_i - 1) - sum_blips_mk1
    }
    term1 <- (H_mk - H_mk_1) - dat[[paste0("preds", m_i, "_", k_i)]]
    term2 <- dat[[paste0("term2_", m_i, k_i)]]
    eps   <- as.numeric(term1 - term2)
    IF_blocks[[i]] <- sweep(Q_list[[i]], 1, eps, `*`)
  }
  IFmat <- do.call(cbind, IF_blocks)
  V <- tryCatch(stats::cov(IFmat) / nrow(dat),
                error = function(e) diag(NA_real_, length(psi_hat)))

  nm_blocks <- unlist(mapply(function(b, Q) paste0("psi", b$m, "_", b$k, "_", colnames(Q)),
                             blips, Q_list, SIMPLIFY = FALSE))

  out <- list(
    call = match.call(),
    coef = stats::setNames(psi_hat, nm_blocks),
    vcov = V,
    solver = sol,
    ntimes = ntimes,
    blips = blips
  )
  class(out) <- "msnmm_pairs_fit"
  if (isTRUE(verbose)) {
    cat("Convergence code:", sol$termcd, "\n")
    if (!is.null(sol$message)) cat(sol$message, "\n")
  }
  out
}
