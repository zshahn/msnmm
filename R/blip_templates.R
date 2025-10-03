# Expand tokens {t},{k},{d} in a character vector
.ms_glue_mk <- function(x, m, k) {
  x <- gsub("\\{t\\}", as.character(m), x, perl = TRUE)
  x <- gsub("\\{k\\}", as.character(k), x, perl = TRUE)
  x <- gsub("\\{d\\}", as.character(k - m), x, perl = TRUE)
  x
}

#' Build unified (m,k) blip blocks for ALL pairs (1 <= m <= k <= K)
#'
#' @param K integer number of time points
#' @param preds_template character vector of predictor column names with optional tokens:
#'   - "{t}" = time m; "{k}" = time k; "{d}" = (k-m)
#'   - "intercept" is allowed (handled by the estimator)
#' @param data_wide optional data.frame (wide) to warn early about missing columns
#' @return list(list(m=.., k=.., preds=..)) suitable for msnmm_tv_pairs()
#' @export
build_blips_all_from_template <- function(K,
                                          preds_template = c("intercept","Ltilde_{t}"),
                                          data_wide = NULL) {
  mk <- subset(expand.grid(m = 1:K, k = 1:K), k >= m)

  if (!is.null(data_wide)) {
    ex <- .ms_glue_mk(preds_template, 1, min(2, K))
    miss <- setdiff(setdiff(ex, "intercept"), names(data_wide))
    if (length(miss)) warning("Potential missing columns: ", paste(miss, collapse = ", "))
  }

  out <- vector("list", nrow(mk))
  for (i in seq_len(nrow(mk))) {
    m <- mk$m[i]; k <- mk$k[i]
    out[[i]] <- list(m = m, k = k, preds = .ms_glue_mk(preds_template, m, k))
  }
  out
}
