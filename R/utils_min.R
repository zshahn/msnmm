# R/utils_min.R

# Ensure diagonal seeds exist: preds{i+1}_{i} := diffs{i}
.seed_preds_diagonal <- function(dat, K, outcome = "Y") {
  for (i in 1:K) {
    nm <- paste0("preds", i + 1, "_", i)
    if (!nm %in% names(dat)) dat[[nm]] <- dat[[paste0("diffs", i)]]
  }
  dat
}

