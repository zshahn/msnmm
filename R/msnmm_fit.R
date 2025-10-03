#' Friendly front door to SNMM g-estimation
#' @export
msnmm_fit <- function(
    data,
    id = NULL,
    K = NULL,
    wide_ready = TRUE,
    # --- outcomes & treatments (wide names) ---
    y_cols = NULL,         # e.g., c("Y_0","Y_1","Y_2")
    a_cols = NULL,         # e.g., c("A_1","A_2")
    # --- easy mode: templates ---
    stems_outcome = NULL,  # e.g., c("S1","S2")
    stems_treat   = NULL,  # e.g., c("S1","S2")
    blip_template = c("intercept","S2_{t}"),   # tokens: {t},{k},{d}
    model = c("glm","spline","ml"),
    df_spline = 5,
    include_pastA = TRUE,
    initiation = FALSE,
    # --- power mode: direct lists ---
    blips = NULL,
    outcome_nuisance_formulas = NULL,
    treatment_nuisance_formulas = NULL,
    verbose = TRUE
) {
  model <- match.arg(model)

  # infer K if not provided
  if (is.null(K)) {
    if (!is.null(y_cols)) {
      K <- length(y_cols) - 1L
    } else {
      # heuristic: look for Y_0..Y_K in data
      y_detect <- sort(grep("^Y_\\d+$", names(data), value = TRUE))
      if (length(y_detect) < 2) stop("Cannot infer K; please supply K or y_cols.")
      idxs <- as.integer(sub("Y_", "", y_detect))
      K <- max(idxs)
    }
  }

  # if user passed direct lists, just forward to engine
  if (!is.null(blips) && !is.null(outcome_nuisance_formulas) && !is.null(treatment_nuisance_formulas)) {
    return(msnmm_tv_pairs_min(
      data = data, id = if (is.null(id)) "id" else id,
      ntimes = K, wide_ready = wide_ready,
      outcome_nuisance_formulas = outcome_nuisance_formulas,
      treatment_nuisance_formulas = treatment_nuisance_formulas,
      blips = blips, initiation = initiation, verbose = verbose
    ))
  }

  # template path: require stems
  if (is.null(stems_outcome) || is.null(stems_treat)) {
    stop("Template mode: please provide stems_outcome and stems_treat (e.g., c('S1','S2')).")
  }

  # choose RHS builder
  rhs_builder <-
    switch(model,
           glm    = rhs_linear,
           spline = function(stems, t, ...) rhs_spline(stems, t, df = df_spline),
           ml     = rhs_linear  # RHS string same; 'ml' will be honored by msnmm_tv_pairs_min later
    )

  # build nuisances
  outcome_fml <- build_outcome_formulas_min(K, stems = stems_outcome, rhs_builder = rhs_builder)
  treat_fml   <- build_treatment_formulas_min(K, stems = stems_treat,
                                              rhs_builder = rhs_builder,
                                              include_pastA = include_pastA)

  # build blips from template (all pairs m<=k)
  blips_built <- build_blips_all_from_template(
    K, preds_template = blip_template, data_wide = if (wide_ready) data else NULL
  )

  # run engine
  msnmm_tv_pairs_min(
    data = data, id = if (is.null(id)) "id" else id,
    ntimes = K, wide_ready = wide_ready,
    outcome_nuisance_formulas = outcome_fml,
    treatment_nuisance_formulas = treat_fml,
    blips = blips_built,
    initiation = initiation,
    verbose = verbose
  )
}
