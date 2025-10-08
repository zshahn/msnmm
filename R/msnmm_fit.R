#' High-level SNMM fit with templates, ML backends, and cross-fitting
#'
#' @param data data.frame in wide format if wide_ready=TRUE; otherwise long
#' @param K number of treatment times
#' @param wide_ready logical
#' @param stems_outcome character stems for outcome RHS builder ("S1","S2",...)
#' @param stems_treat   character stems for treatment RHS builder
#' @param blip_template character, e.g. c("intercept","S2_{t}")
#' @param rhs_outcome builder function or "linear"/"spline"
#' @param rhs_treat   builder function or "linear"/"spline"
#' @param df_outcome,df_treat df for rhs_spline() if used
#' @param model,cf,folds,seed,learners_outcome,learners_treat passed through
#' @param initiation logical; see core
#'
#' @export
msnmm_fit <- function(
    data,
    K,
    id = NULL,
    time = "yearms",
    wide_ready = TRUE,
    stems_outcome,
    stems_treat,
    blip_template = c("intercept","S2_{t}"),
    rhs_outcome = "linear",
    rhs_treat   = "linear",
    df_outcome = 5,
    df_treat   = 5,
    model = c("glm","xgb","sl"),
    cf = FALSE,
    folds = 2,
    seed = NULL,
    learners_outcome = c("SL.glm","SL.ranger"),
    learners_treat   = c("SL.glm","SL.ranger"),
    initiation = FALSE
){
  model <- match.arg(model)

  # pick RHS builders
  builder_from <- function(x, df) {
    if (is.character(x)) {
      x <- match.arg(x, c("linear","spline"))
      if (x == "linear") {
        function(stems, t, ...) rhs_linear(stems, t, ...)
      } else {
        function(stems, t, ...) rhs_spline(stems, t, df = df, ...)
      }
    } else if (is.function(x)) {
      x
    } else stop("rhs_* must be 'linear', 'spline', or a function.")
  }
  rhs_out_fn <- builder_from(rhs_outcome, df_outcome)
  rhs_trt_fn <- builder_from(rhs_treat,   df_treat)

  outcome_fml <- build_outcome_formulas_min(K, stems = stems_outcome, rhs_builder = rhs_out_fn)
  treat_fml   <- build_treatment_formulas_min(K, stems = stems_treat,   rhs_builder = rhs_trt_fn,
                                              include_pastA = TRUE)

  blips <- build_blips_all_from_template(K, preds_template = blip_template,
                                         data_wide = if (isTRUE(wide_ready)) data else NULL)

  msnmm_tv_pairs_min(
    data = data, id = if (is.null(id)) "id" else id, time = time,
    ntimes = K, wide_ready = wide_ready,
    outcome_nuisance_formulas   = outcome_fml,
    treatment_nuisance_formulas = treat_fml,
    blips = blips,
    initiation = initiation,
    verbose = TRUE,
    model = model,
    cf = cf,
    folds = folds,
    seed = seed,
    learners_outcome = learners_outcome,
    learners_treat   = learners_treat
  )
}
