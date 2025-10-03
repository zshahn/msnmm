# ---- RHS helpers ---------------------------------------------------------

rhs_linear <- function(stems, t, ...) {
  paste(c(
    sprintf("%s_%d", stems, t),
    if (t > 1) sprintf("A_%d", 1:(t-1)) else character(0)
  ), collapse = " + ")
}

rhs_spline <- function(stems, t, df = 5, ...) {
  pieces <- c(sprintf("bs(%s_%d, df=%d)", stems, t, df))
  if (t > 1) pieces <- c(pieces, sprintf("A_%d", 1:(t-1)))
  paste(pieces, collapse = " + ")
}

# ---- Outcome builders ----------------------------------------------------

#' Build outcome nuisance formulas (minimal)
build_outcome_formulas_min <- function(K, stems, rhs_builder = rhs_linear, ...) {
  out <- vector("list", K)
  for (i in 1:K) {
    out[[i]] <- vector("list", K - i + 1)
    rhs <- rhs_builder(stems, i, ...)
    for (k in i:K) {
      out[[i]][[k - i + 1]] <- stats::as.formula(sprintf("preds%d_%d ~ %s", i + 1, k, rhs))
    }
  }
  out
}

# ---- Treatment builders --------------------------------------------------

#' Build treatment nuisance formulas (minimal)
build_treatment_formulas_min <- function(K, stems, rhs_builder = rhs_linear, include_pastA = TRUE, ...) {
  out <- vector("list", K)
  for (i in 1:K) {
    rhs <- rhs_builder(stems, i, ...)
    if (!include_pastA && i > 1) {
      # strip prior A terms if rhs_builder added them
      rhs <- gsub("(^|\\+| )A_\\d+", "", rhs)
      rhs <- gsub("\\++", "+", rhs)
      rhs <- gsub("^\\+|\\+$| ", "", rhs)
      if (rhs == "") rhs <- "1"
    }
    out[[i]] <- stats::as.formula(sprintf("A_%d ~ %s", i, rhs))
  }
  out
}

# ---- Blip template expander (for completeness) ---------------------------

build_blips_all_from_template <- function(K, preds_template = c("intercept"), data_wide = NULL) {
  .ms_glue_mk <- function(x, m, k) {
    x <- gsub("\\{t\\}", as.character(m), x, perl = TRUE)
    x <- gsub("\\{k\\}", as.character(k), x, perl = TRUE)
    x <- gsub("\\{d\\}", as.character(k - m), x, perl = TRUE)
    x
  }
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
