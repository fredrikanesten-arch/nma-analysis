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
# MODELS FIT
# ----------
#  M0 : Base-case NMA (no covariates)
#  M1 : + duration_weeks  (trial length as effect modifier)
#  M2 : + mean_age        (patient age)
#  M3 : + pct_female      (sex distribution, % female)
#  M4 : + all three simultaneously
#  (Optionally) M5 : + baseline_mean  (baseline severity)
#  (Optionally) M6 : + rob_high_any   (overall ROB flag)
#
# NOTE: covariates in multinma NMR are study-level (centred
# across all studies) and modify treatment contrasts relative
# to the reference class (Placebo).
#
# OUTPUTS (all written to out_dir)
# ---------------------------------
#  nma_base_case_summary.csv
#  nmr_duration_summary.csv
#  nmr_age_summary.csv
#  nmr_sex_summary.csv
#  nmr_all3_summary.csv
#  model_comparison_dic.csv
#  covariate_effects_all_models.csv
#  multinma_class_objects.rds   (fitted model objects)
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(multinma)   # install via: devtools::install_github("dmphillippo/multinma")

# ------------------------------------------------------------------
# 0.  Paths
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

# ------------------------------------------------------------------
# 1.  Load data
# ------------------------------------------------------------------
message("Loading data ...")
dat <- read_csv(dat_file, show_col_types = FALSE)

# Keep only arms that have a mapped class and non-missing outcome
dat <- dat %>%
  filter(
    !is.na(classcode),
    !is.na(class),
    !is.na(n),
    !is.na(mean_change),
    !is.na(sd_change),
    n > 0,
    sd_change > 0
  )

message(sprintf("  Loaded: %d studies, %d arms", n_distinct(dat$studyid), nrow(dat)))

# ------------------------------------------------------------------
# 2.  Centre and scale covariates
#     Centering is done at the study level (one value per study)
#     using the grand mean across all included studies.
# ------------------------------------------------------------------
study_covs <- dat %>%
  distinct(studyid, duration_weeks, mean_age, pct_female,
           baseline_mean, rob_high_any) %>%
  mutate(
    # Centre around grand mean (missing → imputed with grand mean, i.e. 0 after centring)
    dur_c    = duration_weeks - mean(duration_weeks, na.rm = TRUE),
    age_c    = mean_age       - mean(mean_age,       na.rm = TRUE),
    sex_c    = pct_female     - mean(pct_female,     na.rm = TRUE),
    base_c   = baseline_mean  - mean(baseline_mean,  na.rm = TRUE),
    # Replace NA with 0 (grand-mean imputation; keeps study in the analysis)
    dur_c    = ifelse(is.na(dur_c),  0, dur_c),
    age_c    = ifelse(is.na(age_c),  0, age_c),
    sex_c    = ifelse(is.na(sex_c),  0, sex_c),
    base_c   = ifelse(is.na(base_c), 0, base_c),
    rob_flag = ifelse(is.na(rob_high_any), 0L, as.integer(rob_high_any))
  ) %>%
  select(studyid, dur_c, age_c, sex_c, base_c, rob_flag)

message(sprintf(
  "  Covariate grand-mean centres: duration=%.1f wk, age=%.1f yr, pct_female=%.1f%%",
  mean(dat$duration_weeks, na.rm = TRUE),
  mean(dat$mean_age, na.rm = TRUE),
  mean(dat$pct_female, na.rm = TRUE)
))

# Merge centred covariates back into arm-level data
dat <- dat %>%
  left_join(study_covs, by = "studyid")

# ------------------------------------------------------------------
# 3.  Collapse arms within study×class (same as netmeta_class_ms.R)
# ------------------------------------------------------------------
dat_class <- dat %>%
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

stopifnot(REF_CLASS %in% dat_class$class)
message(sprintf("  After class collapse: %d studies, %d class-arms",
                n_distinct(dat_class$studyid), nrow(dat_class)))

# ------------------------------------------------------------------
# 4.  Build multinma network object
#     multinma expects arm-based data with columns:
#       studyid, treatment, y (mean), sd, n
# ------------------------------------------------------------------
net <- set_agd_arm(
  data       = dat_class,
  study      = studyid,
  trt        = class,
  y          = mean_change,
  se         = sd_change / sqrt(n),
  # Pass covariates as extra columns for use in regression models
  trt_ref    = REF_CLASS
)

# ------------------------------------------------------------------
# 5.  Helper: fit a model, extract summary, compute DIC
# ------------------------------------------------------------------
fit_nmr <- function(net, regression_formula = NULL, label = "model") {
  message(sprintf("  Fitting %s ...", label))

  fit <- nma(
    net,
    trt_effects  = "random",
    regression   = regression_formula,
    prior_trt    = normal(scale = 10),
    prior_het    = half_normal(scale = 0.5),
    prior_reg    = normal(scale = 1),   # vague prior on regression coefficients
    chains       = CHAINS,
    iter         = ITER,
    seed         = SEED,
    show_messages = FALSE
  )

  # Summary of relative effects vs Placebo
  rel_eff <- relative_effects(fit, trt_ref = REF_CLASS) %>%
    as.data.frame() %>%
    mutate(model = label)

  # DIC
  dic_val <- tryCatch(dic(fit), error = function(e) NA_real_)

  list(fit = fit, rel_eff = rel_eff, dic = dic_val, label = label)
}

# ------------------------------------------------------------------
# 6.  Run models
# ------------------------------------------------------------------
message("\nFitting models ...")

results <- list()

# M0 – base case
results[["M0"]] <- fit_nmr(net, regression_formula = NULL, label = "M0_base")

# M1 – duration
results[["M1"]] <- fit_nmr(
  net,
  regression_formula = ~ dur_c,
  label = "M1_duration"
)

# M2 – age
results[["M2"]] <- fit_nmr(
  net,
  regression_formula = ~ age_c,
  label = "M2_age"
)

# M3 – sex
results[["M3"]] <- fit_nmr(
  net,
  regression_formula = ~ sex_c,
  label = "M3_sex"
)

# M4 – all three simultaneously
results[["M4"]] <- fit_nmr(
  net,
  regression_formula = ~ dur_c + age_c + sex_c,
  label = "M4_duration_age_sex"
)

# M5 – baseline severity (optional; comment out if convergence issues)
results[["M5"]] <- fit_nmr(
  net,
  regression_formula = ~ base_c,
  label = "M5_baseline_severity"
)

# M6 – ROB flag (optional)
results[["M6"]] <- fit_nmr(
  net,
  regression_formula = ~ rob_flag,
  label = "M6_rob"
)

# ------------------------------------------------------------------
# 7.  Extract and save covariate (regression) effects
# ------------------------------------------------------------------
message("\nExtracting covariate effects ...")

extract_reg_coefs <- function(res) {
  fit   <- res$fit
  label <- res$label
  coef_names <- names(fit$stanfit@sim$samples[[1]])
  reg_pars  <- grep("^beta\\[", coef_names, value = TRUE)
  if (length(reg_pars) == 0) return(NULL)

  posterior::as_draws_df(fit$stanfit) %>%
    select(dplyr::any_of(reg_pars)) %>%
    tidyr::pivot_longer(everything(), names_to = "parameter") %>%
    group_by(parameter) %>%
    summarise(
      mean   = mean(value),
      sd     = sd(value),
      q2.5   = quantile(value, 0.025),
      q50    = quantile(value, 0.50),
      q97.5  = quantile(value, 0.975),
      rhat   = posterior::rhat(as.vector(value)),
      .groups = "drop"
    ) %>%
    mutate(model = label)
}

coef_table <- dplyr::bind_rows(lapply(results, extract_reg_coefs))

# ------------------------------------------------------------------
# 8.  Model comparison (DIC)
# ------------------------------------------------------------------
dic_table <- tibble::tibble(
  model       = sapply(results, `[[`, "label"),
  DIC         = sapply(results, `[[`, "dic"),
  delta_DIC   = DIC - min(DIC, na.rm = TRUE)
) %>%
  arrange(DIC)

# ------------------------------------------------------------------
# 9.  Combine relative-effect tables
# ------------------------------------------------------------------
all_rel_eff <- dplyr::bind_rows(lapply(results, `[[`, "rel_eff"))

# ------------------------------------------------------------------
# 10.  Save outputs
# ------------------------------------------------------------------
message("Saving outputs ...")

write_csv(
  all_rel_eff,
  file.path(out_dir, "nmr_relative_effects_all_models.csv")
)
write_csv(
  coef_table,
  file.path(out_dir, "nmr_covariate_effects.csv")
)
write_csv(
  dic_table,
  file.path(out_dir, "nmr_model_comparison_dic.csv")
)

# Save individual model summaries
for (nm in names(results)) {
  tryCatch(
    write_csv(
      results[[nm]]$rel_eff,
      file.path(out_dir, paste0("nmr_", results[[nm]]$label, "_rel_eff.csv"))
    ),
    error = function(e) message(sprintf("  Could not save %s: %s", nm, e$message))
  )
}

# Save fitted objects (large; remove to save disk space)
saveRDS(
  lapply(results, `[[`, "fit"),
  file = file.path(out_dir, "multinma_fitted_objects.rds")
)

# ------------------------------------------------------------------
# 11.  Print summary to console
# ------------------------------------------------------------------
message("\n=== Model comparison (DIC) ===")
print(as.data.frame(dic_table))

message("\n=== Covariate effects summary ===")
if (nrow(coef_table) > 0) {
  print(as.data.frame(coef_table %>% select(model, parameter, mean, q2.5, q97.5)))
} else {
  message("  No regression coefficients extracted (check Stan parameter names).")
}

message(sprintf("\nAll outputs written to: %s", out_dir))
