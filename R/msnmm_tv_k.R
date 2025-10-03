# Unified-(m,k) marginal SNMM g-estimation (no numerical guardrails)
# - Wide data: one row per unit after preprocessing
# - Nuisances: lm() as in your script (so you can see convergence / identifiability plainly)
# - Blips: one list of (m, k, preds) entries; "immediate" is just k = m+1
# - Initiation option: if TRUE, zero out future A once A turns on

suppressPackageStartupMessages({
  library(dplyr)
  library(nleqslv)
})

`%nin%` <- Negate(`%in%`)

# --------------------------
# helpers: wide transforms
# --------------------------

create_wide_format <- function(data, id, time_vars, max_timepoints) {
  time_fns <- lapply(0:max_timepoints, function(i) {
    if (i == 0) as.formula("~first(.)") else as.formula(paste0("~nth(., ", i + 1, ")"))
  })
  names(time_fns) <- as.character(0:max_timepoints)
  data %>%
    group_by(!!rlang::sym(id)) %>%
    mutate(across(all_of(time_vars), time_fns, .names = "{.col}_{.fn}")) %>%
    ungroup()
}

diffs <- function(data, outcome, max_timepoints){
  for(i in 1:max_timepoints) {
    curr <- rlang::sym(paste0(outcome, "_", i))
    prev <- rlang::sym(paste0(outcome, "_", i-1))
    data <- data %>% mutate(!!paste0("diffs", i) := !!curr - !!prev)
  }
  data
}

pastA <- function(data, exposure, max_timepoints){
  for(i in 2:max_timepoints) {
    data <- data %>%
      mutate(!!paste0("past", exposure, "_", i) :=
               .data[[paste0(exposure, "_", i-1)]] == 1)
  }
  data
}

# --------------------------
# unified blip estimator
# --------------------------

# blips: a list of entries, each a list(m = <int>, k = <int>, preds = <character vec>)
#   - m, k are 1-based time indices to match your script: Y_0, A_1,...,A_ntimes, Y_ntimes
#   - Example entry: list(m = 1, k = 2, preds = c("intercept", "lag_pov_1"))
#   - Immediate blip is just k = m + 1
#
# outcome_nuisance_formulas:
#   list of length ntimes; for each i, a list of length (ntimes - i + 1)
#   giving formula for preds{i+1}_{i+j-1} ~ f(L_i^{(j)}) exactly as in your script
#
# treatment_nuisance_formulas:
#   list of length ntimes; for each i, a single formula A_i ~ f(L_i)
#
# initiation: if TRUE, force A_i = 0 and A_i_hat = 0 after first initiation
#
# mtimes: optional int vector of m’s you want to include in the stack (defaults to all m that appear in `blips`)
#
msnmm_tv_pairs <- function(
    data,                     # long data (we will wide-ify), or pass wide_ready=TRUE if already wide per-unit
    id,                       # id column name
    time = "yearms",          # time column name (0..ntimes)
    ntimes,                   # number of Y increments and A’s (times 1..ntimes)
    exposure = "A", outcome = "Y",
    time_vars_to_widen = c("A","Y"),  # add others (lag_* etc.) before calling this if desired
    wide_ready = FALSE,       # set TRUE if `data` already has one row per unit (time==1) with A_i, Y_i, etc.
    outcome_nuisance_formulas,
    treatment_nuisance_formulas,
    blips,                    # unified (m,k,preds) list
    initiation = TRUE,
    verbose = TRUE
){
  dat <- data

  # 1) Wide format + basic derived columns -----------------------------------
  if (!wide_ready) {
    dat <- create_wide_format(dat, id = id,
                              time_vars = c(time_vars_to_widen, setdiff(colnames(dat), c(id, time_vars_to_widen, time))),
                              max_timepoints = ntimes)

    # make diffs and pastA on the wide version (but still has multiple rows per unit)
    dat <- diffs(dat, outcome, ntimes)
    dat <- pastA(dat, exposure, ntimes)

    # restrict to one row per unit (time==1, like your script)
    dat <- dat[dat[[time]] == 1, , drop = FALSE]
  } else {
    # assume caller already constructed:
    #  - A_1..A_ntimes, Y_0..Y_ntimes, diffs1..diffs_ntimes, pastA_2..pastA_ntimes
    if (!all(paste0("A_", 1:ntimes) %in% names(dat))) stop("wide_ready=TRUE but A_i not found")
    if (!all(paste0("Y_", 0:ntimes) %in% names(dat))) stop("wide_ready=TRUE but Y_i not found")
    if (!all(paste0("diffs", 1:ntimes) %in% names(dat))) dat <- diffs(dat, outcome, ntimes)
    if (!all(paste0("pastA_", 2:ntimes) %in% names(dat))) dat <- pastA(dat, exposure, ntimes)
  }

  # Optional: enforce "initiation" (first-treatment-only coding)
  if (isTRUE(initiation)) {
    for(i in 2:ntimes){
      hit <- dat[[paste0("pastA_", i)]] == TRUE
      dat[hit, paste0("A_", i)] <- 0
    }
  }

  # 2) Outcome regressions for m=k (preds_{i+1,i}) and then forward-fill preds_{i,k} ----
  # Initialize preds{i+1}_{i} as observed diffs_i (as in your script)
  for(i in 1:ntimes){
    dat[[paste0("preds", i+1, "_", i)]] <- dat[[paste0("diffs", i)]]
  }

  # Fit lm for each exposure_time=i and outcome_time=k >= i:
  # filter: A_i==0 and (if i>1) pastA_i==0
  for (outcome_time in ntimes:1) {
    for (exposure_time in outcome_time:1) {
      form <- outcome_nuisance_formulas[[exposure_time]][[outcome_time - exposure_time + 1]]
      filt <- paste0(exposure, "_", exposure_time, " == 0",
                     if (exposure_time > 1) paste0(" & past", exposure, "_", exposure_time, " == 0") else "")
      df_fit <- dat %>% filter(!!rlang::parse_expr(filt))
      if (nrow(df_fit) == 0L) stop("No rows to fit outcome nuisance at (i=", exposure_time, ", k=", outcome_time, ").")
      mod <- lm(form, data = df_fit)
      pred_name <- paste0("preds", exposure_time, "_", outcome_time)
      dat[[pred_name]] <- as.numeric(predict(mod, newdata = dat))
    }
  }

  # 3) Treatment models A_i ~ L_i  ------------------------------------------
  for (exposure_time in 1:ntimes) {
    form <- treatment_nuisance_formulas[[exposure_time]]
    df_fit <- dat  # your script filtered on time==i originally; here each row is per-unit
    mod <- lm(form, data = df_fit)
    dat[[paste0("A_", exposure_time, "_hat")]] <- as.numeric(predict(mod, newdata = dat))
  }

  if (isTRUE(initiation)) {
    for(i in 2:ntimes){
      hit <- dat[[paste0("pastA_", i)]] == TRUE
      dat[hit, paste0("A_", i, "_hat")] <- 0
    }
  }

  # 4) Weight factors and products ------------------------------------------
  for(i in 1:ntimes){
    dat[[paste0("weightfactor_", i)]] <-
      (1 - dat[[paste0("A_", i)]]) / (1 - dat[[paste0("A_", i, "_hat")]])
  }
  for(m in 1:ntimes){
    for(k in m:ntimes){
      cols <- paste0("weightfactor_", m:k)
      dat[[paste0("weightprod_", m, k)]] <- apply(as.matrix(dat[, cols, drop = FALSE]), 1, prod)
    }
  }

  # 5) Term2: sum_j=m..k weightprod_{m,j} * (preds_{j+1,k} - preds_{j,k})
  for(m in 1:ntimes){
    for(k in m:ntimes){
      for(j in m:k){
        dat[[paste0("term2summands_", m, k, "_", j)]] <-
          dat[[paste0("weightprod_", m, j)]] *
          (dat[[paste0("preds", j+1, "_", k)]] - dat[[paste0("preds", j, "_", k)]])
      }
      cols <- paste0("term2summands_", m, k, "_", m:k)
      dat[[paste0("term2_", m, k)]] <- rowSums(as.matrix(dat[, cols, drop = FALSE]))
    }
  }

  # 6) Build Q_{m,k} design blocks from unified blip spec -------------------
  #    Also define an "intercept" column if requested.
  if (!("intercept" %in% names(dat))) dat$intercept <- 1

  # ensure all referenced predictor columns exist
  for (blk in blips) {
    miss <- setdiff(blk$preds, names(dat))
    if (length(miss)) stop("Missing columns referenced in blip preds: ", paste(miss, collapse = ", "))
  }

  # stack layout: order is the given order in `blips`
  Q_list <- vector("list", length(blips))
  p_block <- integer(length(blips))
  for (i in seq_along(blips)) {
    Q_list[[i]] <- as.matrix(dat[, blips[[i]]$preds, drop = FALSE])
    p_block[i]  <- ncol(Q_list[[i]])
  }
  p_total <- sum(p_block)

  # Optional: limit the set of m’s included (defaults to all m in blips)
  if (exists("mtimes", inherits = FALSE)) {
    keep_idx <- which(vapply(blips, function(b) b$m %in% mtimes, FALSE))
  } else {
    keep_idx <- seq_along(blips)
  }

  # 7) Estimating function ---------------------------------------------------
  #    H_{m,k} = Y_k - sum_{j=m..k} [ A_j * (Q_{j,k} %*% psi_{j,k}) ]
  #    eps_{m,k} = (H_{m,k} - H_{m,k-1} - preds_{m,k}) - term2_{m,k}
  #    moment stack = sum_i Q_{m,k}^T eps_{m,k}

  # fast accessors
  Y <- lapply(0:ntimes, function(t) dat[[paste0("Y_", t)]] )
  A <- lapply(1:ntimes, function(t) dat[[paste0("A_", t)]] )

  # Map (m,k) blocks to their positions in psi
  block_offsets <- c(0L, cumsum(p_block))[seq_along(p_block)]
  block_slice <- function(i) {
    off <- block_offsets[i]
    idx <- (off + 1L):(off + p_block[i])
    idx
  }

  est_eq <- function(psi) {
    # Precompute all blip contributions B_{j->k} = A_j * (Q_{j,k} %*% psi_{j,k})
    # We need these for every (m,k); build a list keyed by (j,k) across all blocks.
    # We'll use names like "j_k".
    B_contrib <- list()
    for (i in seq_along(blips)) {
      m_i <- blips[[i]]$m
      k_i <- blips[[i]]$k
      q   <- Q_list[[i]]
      psi_i <- psi[block_slice(i)]
      g    <- as.numeric(q %*% psi_i)
      key  <- paste0(m_i, "_", k_i)
      B_contrib[[key]] <- A[[m_i]] * g
    }

    # Now compute moments for the kept blocks
    mom <- numeric(p_total)
    for (i in keep_idx) {
      m_i <- blips[[i]]$m
      k_i <- blips[[i]]$k
      # H_{m_i,k_i}
      sum_blips_mk <- 0
      for (j in m_i:k_i) {
        key <- paste0(j, "_", k_i)
        if (!is.null(B_contrib[[key]])) sum_blips_mk <- sum_blips_mk + B_contrib[[key]]
      }
      H_mk   <- Y[[k_i]] - sum_blips_mk
      if (k_i == m_i) {
        H_mk_1 <- Y[[m_i - 1]]  # since Y_0 exists; we assume m_i>=1
      } else {
        # H_{m,k-1} uses sum of blips up to (m..k-1)
        sum_blips_mk1 <- 0
        for (j in m_i:(k_i-1)) {
          key <- paste0(j, "_", k_i-1)
          if (!is.null(B_contrib[[key]])) sum_blips_mk1 <- sum_blips_mk1 + B_contrib[[key]]
        }
        H_mk_1 <- Y[[k_i - 1]] - sum_blips_mk1
      }

      term1 <- (H_mk - H_mk_1) - dat[[paste0("preds", m_i, "_", k_i)]]
      term2 <- dat[[paste0("term2_", m_i, k_i)]]
      eps   <- as.numeric(term1 - term2)

      q <- Q_list[[i]]
      sl <- block_slice(i)
      mom[sl] <- as.numeric(crossprod(q, eps))
    }
    mom
  }

  # 8) Solve
  psi0 <- rep(0, p_total)
  sol  <- nleqslv::nleqslv(x = psi0, fn = est_eq)
  psi_hat <- sol$x

  # 9) (Optional) crude sandwich — block-separable IF for reporting SEs ------
  # Not essential if you only need point estimates right now; keep simple & transparent.
  IF_blocks <- list()
  for (i in keep_idx) {
    m_i <- blips[[i]]$m
    k_i <- blips[[i]]$k

    # recompute with psi_hat
    # (reuse code; brevity over micro-optim here)
    B_contrib <- list()
    for (ii in seq_along(blips)) {
      q   <- Q_list[[ii]]
      psi_ii <- psi_hat[block_slice(ii)]
      g  <- as.numeric(q %*% psi_ii)
      key <- paste0(blips[[ii]]$m, "_", blips[[ii]]$k)
      B_contrib[[key]] <- A[[blips[[ii]]$m]] * g
    }
    sum_blips_mk <- 0
    for (j in m_i:k_i) {
      key <- paste0(j, "_", k_i)
      if (!is.null(B_contrib[[key]])) sum_blips_mk <- sum_blips_mk + B_contrib[[key]]
    }
    H_mk   <- Y[[k_i]] - sum_blips_mk
    if (k_i == m_i) {
      H_mk_1 <- Y[[m_i - 1]]
    } else {
      sum_blips_mk1 <- 0
      for (j in m_i:(k_i-1)) {
        key <- paste0(j, "_", k_i-1)
        if (!is.null(B_contrib[[key]])) sum_blips_mk1 <- sum_blips_mk1 + B_contrib[[key]]
      }
      H_mk_1 <- Y[[k_i - 1]] - sum_blips_mk1
    }
    term1 <- (H_mk - H_mk_1) - dat[[paste0("preds", m_i, "_", k_i)]]
    term2 <- dat[[paste0("term2_", m_i, k_i)]]
    eps   <- as.numeric(term1 - term2)

    q <- Q_list[[i]]
    IF_blocks[[i]] <- sweep(q, 1, eps, `*`)
  }
  IFmat <- do.call(cbind, IF_blocks)
  V <- tryCatch({
    # very rough block-diagonal D (per-block curvature) + empirical IF covariance
    # Here we just use identity curvature for transparency; SEs are for rough checks.
    stats::cov(IFmat) / nrow(dat)
  }, error = function(e) diag(NA_real_, p_total))

  # 10) Names & return -------------------------------------------------------
  nm_blocks <- unlist(mapply(function(b, Q) paste0("psi", b$m, "_", b$k, "_", colnames(Q)),
                             blips, Q_list, SIMPLIFY = FALSE))
  out <- list(
    call = match.call(),
    coef = setNames(psi_hat, nm_blocks),
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

# --------------------------
# Example usage skeleton
# --------------------------
# # 1) Prepare your long data: (id, yearms, A, Y, lag_* vars)
# # 2) Create wide:
# data_wide <- create_wide_format(
#   data = ses_sub_6tp_tvall,
#   id = "stcofips",
#   time_vars = c("A","Y","lag_pop","lag_black","lag_hisp","lag_pov"),
#   max_timepoints = 6
# )
# data_wide <- diffs(data_wide, "Y", ntimes = 6)
# data_wide <- pastA(data_wide, "A", ntimes = 6)
# data_wide <- subset(data_wide, yearms == 1)
#
# # 3) Build nuisance formulas (exactly like your script constructed)
# #    outcome_nuisance_formulas[[i]][[j]] is formula for preds{i}_{i+j-1}
#
# # 4) Define unified blips (examples):
# blips <- list(
#   list(m = 1, k = 2, preds = c("intercept","lag_pov_1")),  # immediate at t=1
#   list(m = 1, k = 3, preds = c("intercept","lag_pov_1")),  # skip 1->3
#   list(m = 2, k = 3, preds = c("intercept","lag_pov_2"))   # immediate at t=2
# )
#
# # 5) Fit
# fit <- msnmm_tv_pairs(
#   data = data_wide, id = "stcofips", time = "yearms", ntimes = 6,
#   exposure = "A", outcome = "Y",
#   wide_ready = TRUE,
#   outcome_nuisance_formulas = outcome_nuisance_formulas,
#   treatment_nuisance_formulas = treatment_nuisance_formulas,
#   blips = blips,
#   initiation = TRUE
# )
# fit$coef
