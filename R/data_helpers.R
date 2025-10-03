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

#' Make diffs ΔY_{t+1} = Y_{t+1} - Y_t for t=0..K-1 (wide names)
#' @export
make_diffs <- function(data, outcome, K) {
  for (i in 1:K) {
    cur <- paste0(outcome, "_", i)
    prev <- paste0(outcome, "_", i - 1)
    data[[paste0("diffs", i)]] <- data[[cur]] - data[[prev]]
  }
  data
}

#' Create past-treatment indicators (for "initiation" setups)
#' @export
make_pastA <- function(data, exposure, K) {
  for (i in 2:K) {
    data[[paste0("past", exposure, "_", i)]] <- as.integer(data[[paste0(exposure, "_", i - 1)]] == 1)
  }
  data
}
