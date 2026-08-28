# ============================================================
# multinma_class_regression.R
#
# PURPOSE
# -------
# Run a class-level network meta-regression (NMR) using the
# multinma package (Dias et al.).  Uses nmr_dataset_long.csv
# produced by prepare_nmr_dataset.R.
#
# The model is arm-based, continuous outcome (SMD), with
# class-level random effects and optional study-level
# covariates as effect modifiers on the treatment contrasts.
#
# DATASET VARIANTS (informed by data_type_consistency_check.R)
# ------------------------------------------------------------
# The BF (Baseline+Followup) block contains ~24 studies with
# extreme within-arm SMDs caused by yB << yF (i.e. the follow-up
# score is much larger than baseline, suggesting scale-switching
# or data-entry issues in those specific studies).  The BIN
# (Binary) block has a median SMD ~2× more negative than CFB,
# reflecting strong assumptions in the Furukawa conversion.
#
# We run NMR on three variants:
#
#   VARIANT A – ALL data (CFB + BF + BIN, no filtering)
#               Primary pre-specified analysis; most powered.
#
#   VARIANT B – BF-winsorised
#               Within each BF study-arm, approx_smd (= mean_change
#               / sd_change) is capped at WINSОR_LO / WINSОR_HI
#               quantiles of the CFB distribution.  mean_change is
#               rescaled accordingly; sd_change is unchanged.
#               Advantage: retains all 172 BF studies in the
#               network while damping outlier influence.
#               Disadvantage: arbitrary threshold; alters raw data.
#
#   VARIANT C – CFB + BF only (no BIN)
#               Excludes 34 BIN studies whose median SMD is
#               systematically 2× that of CFB.  Cleaner
#               methodologically if the Furukawa assumptions are
#               suspect.  Disadvantage: smaller, potentially
#               less-connected network.
#
# MODELS FIT (per variant)
# ------------------------
#  M0 : Base-case NMA (no covariates)
#  M1 : + dur_c        (trial duration, centred)
#  M2 : + age_c        (mean age, centred)
#  M3 : + sex_c        (% female, centred)
#  M4 : + dur_c + age_c + sex_c  (main NMR model)
#  M5 : + base_c       (baseline severity, centred)
#  M6 : + rob_flag     (overall ROB flag, binary)
#
# OUTPUTS (all in out_dir / variant subfolder)
# ---------------------------------------------
#  bf_outlier_diagnostic.csv
#  nmr_<label>_rel_eff.csv          (per model)
#  nmr_covariate_effects_<variant>.csv
#  nmr_model_comparison_dic_<variant>.csv
#  multinma_fitted_<variant>.rds
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(multinma)   # install via: devtools::install_github("dmphillippo/multinma")

# ------------------------------------------------------------------
# 0.  Paths and settings
# ------------------------------------------------------------------
base_dir <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli/netmeta_class_ms"
dat_file <- file.path(base_dir, "nmr_dataset_long.csv")
out_dir  <- file.path(base_dir, "multinma_regression")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Stan / MCMC settings
CHAINS  <- 4L
ITER    <- 2000L   # total iterations per chain (warmup = ITER/2)
SEED    <- 20240101L

# Reference class
REF_CLASS <- "Placebo"

# Winsorising thresholds for BF arms (applied to approx_smd = mean_change/sd_change)
# Defined relative to the CFB arm distribution (2.5th and 97.5th percentiles).
# Change WINSOR_QUANTILE to adjust strictness (e.g. 0.01 for 1st/99th).
WINSOR_QUANTILE <- 0.025

# ------------------------------------------------------------------
# 1.  Load data
# ------------------------------------------------------------------
message("Loading data ...")
dat <- read_csv(dat_file, show_col_types = FALSE) %>%
  filter(
    !is.na(classcode), !is.na(class),
    !is.na(n), !is.na(mean_change), !is.na(sd_change),
    n > 0, sd_change > 0
  )

message(sprintf("  Loaded: %d studies, %d arms", n_distinct(dat$studyid), nrow(dat)))
message(sprintf("  Data types: %s",
  paste(names(table(dat$data_type)),
        as.integer(table(dat$data_type)), sep = "=", collapse = ", ")))

# ------------------------------------------------------------------
# 2.  BF outlier diagnostic (informed by consistency check)
#     approx_smd = mean_change / sd_change (within-arm signal-to-noise)
# ------------------------------------------------------------------
cfb_smd_lo <- quantile(
  (dat %>% filter(data_type == "CFB"))$mean_change /
  (dat %>% filter(data_type == "CFB"))$sd_change,
  WINSOR_QUANTILE, na.rm = TRUE
)
cfb_smd_hi <- quantile(
  (dat %>% filter(data_type == "CFB"))$mean_change /
  (dat %>% filter(data_type == "CFB"))$sd_change,
  1 - WINSOR_QUANTILE, na.rm = TRUE
)
message(sprintf(
  "  CFB approx-SMD %.1f%%/%.1f%% quantiles: [%.3f, %.3f]  (winsorising bounds for BF)",
  WINSOR_QUANTILE * 100, (1 - WINSOR_QUANTILE) * 100,
  cfb_smd_lo, cfb_smd_hi
))

bf_outliers <- dat %>%
  filter(data_type == "BF") %>%
  mutate(approx_smd = mean_change / sd_change) %>%
  filter(approx_smd < cfb_smd_lo | approx_smd > cfb_smd_hi) %>%
  arrange(approx_smd) %>%
  select(studyid, arm, treatment, class, data_type,
         n, mean_change, sd_change, approx_smd)

message(sprintf("  BF arms outside CFB bounds: %d / %d",
                nrow(bf_outliers), sum(dat$data_type == "BF")))

write_csv(bf_outliers, file.path(out_dir, "bf_outlier_diagnostic.csv"))
message("  Saved: bf_outlier_diagnostic.csv")

# ------------------------------------------------------------------
# 3.  Build the three dataset variants
# ------------------------------------------------------------------

## VARIANT A – all data unchanged
dat_A <- dat

## VARIANT B – BF arms winsorised
dat_B <- dat %>%
  mutate(
    approx_smd     = mean_change / sd_change,
    approx_smd_win = pmin(pmax(approx_smd, cfb_smd_lo), cfb_smd_hi),
    # Rescale mean_change to match winsorised approx_smd (sd_change unchanged)
    mean_change    = ifelse(
      data_type == "BF",
      approx_smd_win * sd_change,
      mean_change
    )
  ) %>%
  select(-approx_smd, -approx_smd_win)

n_winsorised <- sum(
  dat %>% filter(data_type == "BF") %>%
    mutate(a = mean_change / sd_change) %>%
    pull(a) < cfb_smd_lo |
  dat %>% filter(data_type == "BF") %>%
    mutate(a = mean_change / sd_change) %>%
    pull(a) > cfb_smd_hi,
  na.rm = TRUE
)
message(sprintf("  Variant B: %d BF arms will have mean_change rescaled", n_winsorised))

## VARIANT C – exclude BIN studies
dat_C <- dat %>% filter(data_type != "BIN")
message(sprintf("  Variant C (no BIN): %d studies, %d arms",
                n_distinct(dat_C$studyid), nrow(dat_C)))

variants <- list(
  A_all         = dat_A,
  B_bf_winsor   = dat_B,
  C_no_bin      = dat_C
)

# ------------------------------------------------------------------
# 4.  Helper: prepare a dataset variant for multinma
#     (centre covariates → collapse to class level → build network)
# ------------------------------------------------------------------
prepare_net <- function(dat_v, variant_label) {

  # Centre covariates on the grand mean of THIS variant's studies
  study_covs <- dat_v %>%
    distinct(studyid, duration_weeks, mean_age, pct_female,
             baseline_mean, rob_high_any) %>%
    mutate(
      dur_c  = duration_weeks - mean(duration_weeks, na.rm = TRUE),
      age_c  = mean_age       - mean(mean_age,       na.rm = TRUE),
      sex_c  = pct_female     - mean(pct_female,     na.rm = TRUE),
      base_c = baseline_mean  - mean(baseline_mean,  na.rm = TRUE),
      dur_c  = ifelse(is.na(dur_c),  0, dur_c),
      age_c  = ifelse(is.na(age_c),  0, age_c),
      sex_c  = ifelse(is.na(sex_c),  0, sex_c),
      base_c = ifelse(is.na(base_c), 0, base_c),
      rob_flag = ifelse(is.na(rob_high_any), 0L, as.integer(rob_high_any))
    ) %>%
    select(studyid, dur_c, age_c, sex_c, base_c, rob_flag)

  dat_v <- dat_v %>% left_join(study_covs, by = "studyid")

  # Collapse multiple arms per study×class
  dat_class <- dat_v %>%
    group_by(studyid, classcode, class,
             dur_c, age_c, sex_c, base_c, rob_flag) %>%
    summarise(
      n           = sum(n),
      mean_change = sum(mean_change * n) / sum(n),
      sd_change   = sqrt(sum((n - 1) * sd_change^2) / sum(n - 1)),
      .groups     = "drop"
    ) %>%
    group_by(studyid) %>%
    filter(n_distinct(class) >= 2) %>%
    ungroup()

  if (!REF_CLASS %in% dat_class$class) {
    stop(sprintf("[%s] Reference class '%s' not found after filtering.",
                 variant_label, REF_CLASS))
  }

  message(sprintf("  [%s] %d studies, %d class-arms after collapse",
                  variant_label,
                  n_distinct(dat_class$studyid), nrow(dat_class)))

  set_agd_arm(
    data       = dat_class,
    study      = studyid,
    trt        = class,
    y          = mean_change,
    se         = sd_change / sqrt(n),
    trt_ref    = REF_CLASS,
    covariates = c("dur_c", "age_c", "sex_c", "base_c", "rob_flag")
  )
}

# ------------------------------------------------------------------
# 5.  Helper: fit one NMR model
# ------------------------------------------------------------------
fit_nmr <- function(net, regression_formula = NULL, label = "model") {
  message(sprintf("    Fitting %s ...", label))

  fit <- tryCatch(
    nma(
      net,
      trt_effects   = "random",
      regression    = regression_formula,
      prior_trt     = normal(scale = 10),
      prior_het     = half_normal(scale = 0.5),
      prior_reg     = normal(scale = 1),
      chains        = CHAINS,
      iter          = ITER,
      seed          = SEED,
      show_messages = FALSE
    ),
    error = function(e) {
      message(sprintf("    [ERROR %s] %s", label, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  rel_eff <- tryCatch(
    relative_effects(fit, trt_ref = REF_CLASS) %>%
      as.data.frame() %>%
      mutate(model = label),
    error = function(e) NULL
  )

  dic_val <- tryCatch({
    d <- dic(fit)
    # multinma returns an nma_dic object; extract scalar total DIC
    if (is.numeric(d)) d[[1]] else if (!is.null(d$dic)) d$dic else as.numeric(d)[[1]]
  }, error = function(e) NA_real_)

  list(fit = fit, rel_eff = rel_eff, dic = dic_val, label = label)
}

# ------------------------------------------------------------------
# 6.  Helper: extract posterior regression coefficients
# ------------------------------------------------------------------
extract_reg_coefs <- function(res) {
  if (is.null(res)) return(NULL)
  fit   <- res$fit
  label <- res$label

  # multinma stores regression params as "beta[covariate,treatment]"
  all_pars <- tryCatch(
    posterior::as_draws_df(fit),
    error = function(e) NULL
  )
  if (is.null(all_pars)) return(NULL)

  reg_pars <- grep("^beta\\[", names(all_pars), value = TRUE)
  if (length(reg_pars) == 0) return(NULL)

  all_pars %>%
    dplyr::select(dplyr::any_of(reg_pars)) %>%
    tidyr::pivot_longer(everything(), names_to = "parameter") %>%
    dplyr::group_by(parameter) %>%
    dplyr::summarise(
      mean  = mean(value),
      sd    = sd(value),
      q2.5  = quantile(value, 0.025),
      q50   = quantile(value, 0.50),
      q97.5 = quantile(value, 0.975),
      rhat  = posterior::rhat(as.vector(value)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(model = label)
}

# ------------------------------------------------------------------
# 7.  Run all models across all variants
# ------------------------------------------------------------------
model_specs <- list(
  M0 = list(formula = NULL,                         label = "M0_base"),
  M1 = list(formula = ~ dur_c,                      label = "M1_duration"),
  M2 = list(formula = ~ age_c,                      label = "M2_age"),
  M3 = list(formula = ~ sex_c,                      label = "M3_sex"),
  M4 = list(formula = ~ dur_c + age_c + sex_c,      label = "M4_dur_age_sex"),
  M5 = list(formula = ~ base_c,                     label = "M5_baseline"),
  M6 = list(formula = ~ rob_flag,                   label = "M6_rob")
)

all_results <- list()   # all_results[[variant]][[model]]

for (vname in names(variants)) {
  message(sprintf("\n=== Variant %s ===", vname))

  net_v <- tryCatch(
    prepare_net(variants[[vname]], vname),
    error = function(e) { message("  [ERROR] ", conditionMessage(e)); NULL }
  )
  if (is.null(net_v)) next

  v_results <- list()
  for (mname in names(model_specs)) {
    spec <- model_specs[[mname]]
    full_label <- paste0(vname, "_", spec$label)
    v_results[[mname]] <- fit_nmr(
      net_v,
      regression_formula = spec$formula,
      label              = full_label
    )
  }
  all_results[[vname]] <- v_results
}

# ------------------------------------------------------------------
# 8.  Save outputs per variant
# ------------------------------------------------------------------
message("\nSaving outputs ...")

for (vname in names(all_results)) {
  v_results <- all_results[[vname]]
  vdir <- file.path(out_dir, vname)
  if (!dir.exists(vdir)) dir.create(vdir)

  # Relative effects (all models combined)
  rel_all <- dplyr::bind_rows(lapply(v_results, function(r) r$rel_eff))
  if (nrow(rel_all) > 0)
    write_csv(rel_all, file.path(vdir, "nmr_relative_effects_all_models.csv"))

  # Per-model relative effects
  for (mname in names(v_results)) {
    r <- v_results[[mname]]
    if (!is.null(r) && !is.null(r$rel_eff))
      write_csv(r$rel_eff,
                file.path(vdir, paste0("nmr_", r$label, "_rel_eff.csv")))
  }

  # Covariate effects
  coef_tbl <- dplyr::bind_rows(lapply(v_results, extract_reg_coefs))
  if (!is.null(coef_tbl) && nrow(coef_tbl) > 0)
    write_csv(coef_tbl, file.path(vdir, "nmr_covariate_effects.csv"))

  # DIC comparison
  dic_tbl <- tibble::tibble(
    model     = vapply(v_results, function(r) if (!is.null(r)) r$label else NA_character_, character(1)),
    DIC       = vapply(v_results, function(r) {
      if (is.null(r)) return(NA_real_)
      d <- r$dic
      # d may be: scalar, named numeric, nma_dic object (list), or NA
      tryCatch(as.numeric(unlist(d))[[1]], error = function(e) NA_real_)
    }, numeric(1))
  ) %>%
    filter(!is.na(DIC)) %>%
    mutate(delta_DIC = DIC - min(DIC, na.rm = TRUE)) %>%
    arrange(DIC)
  write_csv(dic_tbl, file.path(vdir, "nmr_model_comparison_dic.csv"))

  # Fitted objects
  saveRDS(
    lapply(v_results, function(r) r$fit),
    file.path(vdir, paste0("multinma_fitted_", vname, ".rds"))
  )

  message(sprintf("  Saved outputs for variant %s to: %s", vname, vdir))
}

# ------------------------------------------------------------------
# 9.  Cross-variant comparison table
#     Compares M0 (base NMA) DIC and heterogeneity across variants
# ------------------------------------------------------------------
message("\n=== Cross-variant summary (M0 base NMA) ===")

for (vname in names(all_results)) {
  m0 <- all_results[[vname]][["M0"]]
  if (!is.null(m0)) {
    dic_scalar <- tryCatch(as.numeric(unlist(m0$dic))[[1]], error = function(e) NA_real_)
    message(sprintf("  %-20s  DIC = %s", vname,
                    if (is.na(dic_scalar)) "NA" else sprintf("%.1f", dic_scalar)))
  }
}

# ------------------------------------------------------------------
# 10.  Print covariate effects for M4 across variants
# ------------------------------------------------------------------
message("\n=== M4 covariate effects (dur + age + sex) across variants ===")
for (vname in names(all_results)) {
  m4 <- all_results[[vname]][["M4"]]
  if (!is.null(m4)) {
    coefs <- extract_reg_coefs(m4)
    if (!is.null(coefs) && nrow(coefs) > 0) {
      message(sprintf("  --- %s ---", vname))
      print(as.data.frame(coefs %>% select(parameter, mean, q2.5, q97.5)))
    }
  }
}

message(sprintf("\nAll outputs written to: %s", out_dir))

