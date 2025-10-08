#' Create wide-style lagged columns Y_t, A_t, and other time-varying vars
#' @param time_vars character vector of columns to expand
#' @param max_timepoints numeric K
#' @export
create_wide_format <- function(data, id, time_vars, max_timepoints) {
  time_fns <- lapply(0:max_timepoints, function(i) {
    if (i == 0) stats::as.formula("~first(.)") else stats::as.formula(paste0("~nth(., ", i + 1, ")"))
  })
  names(time_fns) <- as.character(0:max_timepoints)

  data |>
    dplyr::group_by(.data[[id]]) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(time_vars), time_fns, .names = "{.col}_{.fn}")) |>
    dplyr::ungroup()
}


# R/data_helpers.R

#' Create diffs1..diffsK from wide Y_0..Y_K
#' @export
make_diffs <- function(dat, outcome = "Y", K) {
  # sanity
  need <- paste0(outcome, "_", 0:K)
  if (!all(need %in% names(dat))) {
    stop("make_diffs(): missing columns: ", paste(setdiff(need, names(dat)), collapse = ", "))
  }
  for (i in 1:K) {
    dat[[paste0("diffs", i)]] <- dat[[paste0(outcome, "_", i)]] - dat[[paste0(outcome, "_", i-1)]]
  }
  dat
}

#' Create pastA_2..pastA_K flags from A_1..A_K (ever-treated before time i)
#'
#' pastA_i = 1 if any A_1..A_{i-1} == 1, else 0. This supports the
#' "no-past-A rows" filtering you use in the outcome nuisance regressions.
#'
#' @export
make_pastA <- function(dat, exposure = "A", K) {
  need <- paste0(exposure, "_", 1:K)
  if (!all(need %in% names(dat))) {
    stop("make_pastA(): missing columns: ", paste(setdiff(need, names(dat)), collapse = ", "))
  }
  for (i in 2:K) {
    prior <- dat[, paste0(exposure, "_", 1:(i-1)), drop = FALSE]
    # treat NA as 0 so 'any' stays robust
    prior[is.na(prior)] <- 0
    dat[[paste0("past", exposure, "_", i)]] <- as.integer(rowSums(prior != 0) > 0)
  }
  dat
}

