#' Unified (m,k) SNMM g-estimation (minimal)
#' @param data wide data
#' @param id ignored when wide_ready=TRUE
#' @param ntimes integer horizon
#' @param outcome_nuisance_formulas list[[i]][[k-i+1]] of formulas
#' @param treatment_nuisance_formulas list[[i]] of formulas
#' @param blips list of blocks: list(m=,k=,preds=)
#' @param initiation logical
#' @return object of class \code{msnmm_pairs_fit}
#' @export
msnmm_tv_pairs_min <- function(
    data, id, time = "yearms", ntimes,
    exposure = "A", outcome = "Y",
    time_vars_to_widen = c("A","Y"),
    wide_ready = FALSE,
    outcome_nuisance_formulas,
    treatment_nuisance_formulas,
    blips,
    initiation = TRUE,
    verbose = TRUE
){
  dat <- data
  if (isTRUE(wide_ready)) {
    # diffs
    need_d <- paste0("diffs", 1:ntimes)
    if (!all(need_d %in% names(dat))) {
      for (k in 1:ntimes)
        dat[[paste0("diffs", k)]] <- dat[[paste0(outcome, "_", k)]] - dat[[paste0(outcome, "_", k-1)]]
    }
    # pastA
    need_p <- paste0("pastA_", 2:ntimes)
    if (ntimes >= 2 && !all(need_p %in% names(dat))) {
      for (i in 2:ntimes) {
        prev <- paste0(exposure, "_", 1:(i-1))
        dat[[paste0("pastA_", i)]] <- as.integer(Reduce(`|`, lapply(prev, \(v) dat[[v]] == 1)))
      }
    }
  }
  # ----- helpers (same as before) -------------------------------------------
  create_wide_format <- function(data, id, time_vars, max_timepoints) {
    time_fns <- lapply(0:max_timepoints, function(i) {
      if (i == 0) as.formula("~first(.)") else as.formula(paste0("~nth(., ", i + 1, ")"))
    })
    names(time_fns) <- as.character(0:max_timepoints)
    dplyr::group_by(data, !!rlang::sym(id)) |>
      dplyr::mutate(dplyr::across(dplyr::all_of(time_vars), time_fns, .names = "{.col}_{.fn}")) |>
      dplyr::ungroup()
  }
  diffs <- function(data, outcome, max_timepoints){
    for(i in 1:max_timepoints) {
      curr <- rlang::sym(paste0(outcome, "_", i))
      prev <- rlang::sym(paste0(outcome, "_", i-1))
      data <- dplyr::mutate(data, !!paste0("diffs", i) := !!curr - !!prev)
    }
    data
  }
  pastA_flag <- function(df, exposure, i) {
    cols_prev <- .prevA_cols(i, exposure)
    if (!length(cols_prev)) return(rep(FALSE, nrow(df)))
    rowSums(df[, cols_prev, drop = FALSE] != 0) > 0
  }

  # ----- 1) Wide format ------------------------------------------------------
  if (!isTRUE(wide_ready)) {
    dat <- create_wide_format(
      dat, id = id,
      time_vars = c(time_vars_to_widen, setdiff(colnames(dat), c(id, time_vars_to_widen, time))),
      max_timepoints = ntimes
    )
    dat <- diffs(dat, outcome, ntimes)
    dat <- dat[dat[[time]] == 1, , drop = FALSE]
  } else {
    if (!all(paste0("A_", 1:ntimes) %in% names(dat))) stop("wide_ready=TRUE but A_i not found")
    if (!all(paste0("Y_", 0:ntimes) %in% names(dat))) stop("wide_ready=TRUE but Y_i not found")
    if (!all(paste0("diffs", 1:ntimes) %in% names(dat))) dat <- diffs(dat, outcome, ntimes)
  }

  # ----------------- 2) Outcome regressions: preds{i}_{k} --------------------
  for (i in 1:ntimes) {
    dat[[paste0("preds", i+1, "_", i)]] <- dat[[paste0("diffs", i)]]
  }
  for (k in ntimes:1) {
    for (i in k:1) {
      obj <- outcome_nuisance_formulas[[i]][[k - i + 1]]  # can be formula OR learner
      # rows: A_i == 0 & (if i>1) no past A_i
      filt <- paste0(exposure, "_", i, " == 0",
                     if (i > 1) paste0(" & past", exposure, "_", i, " == 0") else "")
      rows <- rlang::eval_tidy(rlang::parse_expr(filt), data = dat)
      if (!any(rows)) stop("No rows to fit outcome nuisance at (i=", i, ", k=", k, ").")

      y_vec <- dat[[paste0("diffs", k)]]

      if (inherits(obj, "formula")) {
        mod <- stats::lm(obj, data = dat[rows, , drop = FALSE])
        dat[[paste0("preds", i, "_", k)]] <- as.numeric(stats::predict(mod, newdata = dat))
      } else if (inherits(obj, "ms_learner")) {
        # Build the learner's X from the vars it captured
        mdl <- obj$fit(dat, y = y_vec, rows = rows)
        dat[[paste0("preds", i, "_", k)]] <- mdl$predict(dat)
      } else {
        stop("Unknown outcome nuisance spec at (i=", i, ", k=", k, ").")
      }
    }
  }

  # ----------------- 3) Treatment models A_i ~ L_i ---------------------------
  for (i in 1:ntimes) {
    obj <- treatment_nuisance_formulas[[i]]  # can be formula OR learner

    # (Optional) if you added past-A columns into the stems, they'll already be present in dat
    y_vec <- dat[[paste0("A_", i)]]

    if (inherits(obj, "formula")) {
      # logistic if you want: stats::glm(obj, data=dat, family=binomial())
      mod <- stats::glm(obj, data = dat, family = binomial())
      dat[[paste0("A_", i, "_hat")]] <- as.numeric(stats::predict(mod, newdata = dat, type = "response"))
    } else if (inherits(obj, "ms_learner")) {
      mdl <- obj$fit(dat, y = y_vec, rows = rep(TRUE, nrow(dat)))
      dat[[paste0("A_", i, "_hat")]] <- mdl$predict(dat)
    } else {
      stop("Unknown treatment nuisance spec at i=", i, ".")
    }
  }

  # # ----- 2) Seed diagonal preds and fit outcome nuisances with 'no past A' ---
  # for(i in 1:ntimes){
  #   dat[[paste0("preds", i+1, "_", i)]] <- dat[[paste0("diffs", i)]]
  # }
  # for (k in ntimes:1) {
  #   for (i in k:1) {
  #     form <- outcome_nuisance_formulas[[i]][[k - i + 1]]
  #     # fit on rows with A_i == 0 AND no past A before i (initiation context)
  #     no_past <- !pastA_flag(dat, exposure, i)
  #     df_fit  <- dat[no_past & dat[[paste0(exposure, "_", i)]] == 0, , drop = FALSE]
  #     if (nrow(df_fit) == 0L) stop(sprintf("No rows to fit outcome nuisance at (i=%d,k=%d).", i, k))
  #     mod <- stats::lm(form, data = df_fit)
  #     dat[[paste0("preds", i, "_", k)]] <- as.numeric(stats::predict(mod, newdata = dat))
  #   }
  # }
  #
  # # ----- 3) Treatment models A_i ~ (stems at i) + A_{i-1} if i>1 -------------
  # for (i in 1:ntimes) {
  #   form <- treatment_nuisance_formulas[[i]]
  #   mod  <- stats::glm(form, data = dat, family = stats::binomial())
  #   dat[[paste0("A_", i, "_hat")]] <- as.numeric(stats::predict(mod, newdata = dat, type = "response"))
  # }

  # ----- 4) Initiation coding (first-treatment only) -------------------------
  if (isTRUE(initiation)) {
    for(i in 2:ntimes){
      has_past <- pastA_flag(dat, exposure, i)
      dat[has_past, paste0("A_", i)]      <- 0
      dat[has_past, paste0("A_", i, "_hat")] <- 0
    }
  }

  # ----- 5) Weights and products (unchanged) ---------------------------------
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

  # ----- 6) Build Y/A lists BEFORE dropping inactive m -----------------------
  Y <- lapply(0:ntimes, function(t) dat[[paste0("Y_", t)]])
  A <- lapply(1:ntimes, function(t) dat[[paste0("A_", t)]])
  Yt <- function(t) Y[[t + 1L]]

  # drop inactive m
  active_m <- which(vapply(1:ntimes, function(m) any(A[[m]] == 1, na.rm = TRUE), FALSE))
  keep <- vapply(blips, function(b) b$m %in% active_m, FALSE)
  if (!all(keep)) blips <- blips[keep]
  if (length(blips) == 0L) stop("No blip blocks remain after dropping inactive m.")

  # ----- 7) Blip design ------------------------------------------------------
  if (!("intercept" %in% names(dat))) dat$intercept <- 1
  for (blk in blips) {
    miss <- setdiff(blk$preds, names(dat))
    if (length(miss)) stop("Missing columns referenced in blip preds: ", paste(miss, collapse = ", "))
  }
  Q_list  <- lapply(blips, function(b) as.matrix(dat[, b$preds, drop = FALSE]))
  p_block <- vapply(Q_list, ncol, 0L)
  p_total <- sum(p_block)
  block_offsets <- c(0L, cumsum(p_block))[seq_along(p_block)]
  block_slice <- function(i) { off <- block_offsets[i]; (off + 1L):(off + p_block[i]) }

  # ----- 8) Estimating equations --------------------------------------------
  est_eq <- function(psi) {
    # B_{j->k}
    B_contrib <- list()
    for (i in seq_along(blips)) {
      m_i <- blips[[i]]$m; k_i <- blips[[i]]$k
      g   <- as.numeric(Q_list[[i]] %*% psi[block_slice(i)])
      B_contrib[[paste0(m_i, "_", k_i)]] <- A[[m_i]] * g
    }

    mom <- numeric(p_total)
    for (i in seq_along(blips)) {
      m_i <- blips[[i]]$m; k_i <- blips[[i]]$k
      # H_{m,k}
      sum_blips_mk <- 0
      for (j in m_i:k_i) {
        key <- paste0(j, "_", k_i)
        if (!is.null(B_contrib[[key]])) sum_blips_mk <- sum_blips_mk + B_contrib[[key]]
      }
      H_mk <- Yt(k_i) - sum_blips_mk
      # H_{m,k-1}
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

  psi0 <- rep(0, p_total)
  sol  <- nleqslv::nleqslv(x = psi0, fn = est_eq)
  psi_hat <- sol$x

  nm_blocks <- unlist(mapply(function(b, Q) paste0("psi", b$m, "_", b$k, "_", colnames(Q)),
                             blips, Q_list, SIMPLIFY = FALSE))

  out <- list(
    call = match.call(),
    coef = stats::setNames(psi_hat, nm_blocks),
    vcov = NULL,
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
