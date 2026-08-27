# ============================================================
# data_type_consistency_check.R
#
# PURPOSE
# -------
# Analyse whether the three data formats in ms_depression_data
# (CFB = direct change-from-baseline, BF = baseline+followup
# converted to CFB, BIN = binary response converted to CFB)
# yield consistent SMD estimates.
#
# APPROACH
# --------
# 1. Separate NMAs by data_type (CFB-only, CFB+BF, ALL).
# 2. Compare class-vs-Placebo SMDs across the three networks
#    (scatter plots + correlation).
# 3. Compute per-arm "approximate within-arm SMD" (mean_change /
#    sd_change) and compare distributions by data_type.
# 4. Run a simple meta-regression with data_type as covariate
#    to test whether BF/BIN studies differ systematically from CFB.
# 5. Export tables and plots.
#
# OUTPUTS (written to out_dir)
# -----------------------------
#  consistency_arm_smd_by_type.csv      – per-arm approx SMD
#  consistency_network_estimates.csv    – NMA estimates per type
#  consistency_comparison_plot.pdf      – scatter + violin plots
#  consistency_datatype_metareg.csv     – meta-reg coefficient table
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(netmeta)

# ------------------------------------------------------------------
# 0.  Paths
# ------------------------------------------------------------------
base_dir <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli/netmeta_class_ms"
dat_file <- file.path(base_dir, "nmr_dataset_long.csv")
out_dir  <- file.path(base_dir, "consistency_check")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

REF_CLASS <- "Placebo"

# ------------------------------------------------------------------
# 1.  Load and prepare data
# ------------------------------------------------------------------
message("Loading data ...")
dat <- read_csv(dat_file, show_col_types = FALSE) %>%
  filter(
    !is.na(classcode), !is.na(class),
    !is.na(n), !is.na(mean_change), !is.na(sd_change),
    n > 0, sd_change > 0
  )

message(sprintf("  %d studies, %d arms loaded", n_distinct(dat$studyid), nrow(dat)))
message(sprintf("  Data type breakdown:"))
dat %>% count(data_type) %>% { message(sprintf("    %s", paste(.$data_type, .$n, sep = ": ", collapse = ", "))) }

# ------------------------------------------------------------------
# 2.  Approximate within-arm SMD
#     mean_change / sd_change  gives the within-arm effect size
#     (signed; negative = improvement on typical depression scales)
# ------------------------------------------------------------------
arm_smd <- dat %>%
  mutate(
    approx_smd = mean_change / sd_change,
    data_type  = factor(data_type, levels = c("CFB", "BF", "BIN"))
  )

# Save
write_csv(arm_smd %>%
            select(studyid, arm, treatment, class, data_type,
                   n, mean_change, sd_change, approx_smd),
          file.path(out_dir, "consistency_arm_smd_by_type.csv"))

# Summary statistics by data_type
type_summary <- arm_smd %>%
  group_by(data_type) %>%
  summarise(
    n_arms    = n(),
    n_studies = n_distinct(studyid),
    mean_smd  = mean(approx_smd, na.rm = TRUE),
    sd_smd    = sd(approx_smd, na.rm = TRUE),
    median_smd= median(approx_smd, na.rm = TRUE),
    q25       = quantile(approx_smd, 0.25, na.rm = TRUE),
    q75       = quantile(approx_smd, 0.75, na.rm = TRUE),
    .groups   = "drop"
  )
message("\n--- Within-arm approx SMD by data_type ---")
print(as.data.frame(type_summary))

# ------------------------------------------------------------------
# 3.  Helper: class-collapse → pairwise → NMA
# ------------------------------------------------------------------
run_nma <- function(dat_subset, label = "") {
  collapsed <- dat_subset %>%
    group_by(studyid, classcode, class) %>%
    summarise(
      n           = sum(n),
      mean_change = sum(mean_change * n) / sum(n),
      sd_change   = sqrt(sum((n - 1) * sd_change^2) / sum(n - 1)),
      .groups     = "drop"
    ) %>%
    group_by(studyid) %>%
    filter(n_distinct(class) >= 2) %>%
    ungroup()

  if (!REF_CLASS %in% collapsed$class) {
    message(sprintf("  [SKIP %s] Reference class '%s' not present.", label, REF_CLASS))
    return(NULL)
  }

  pw <- tryCatch(
    pairwise(
      treat   = class,
      mean    = mean_change,
      sd      = sd_change,
      n       = n,
      studlab = studyid,
      data    = collapsed,
      sm      = "SMD"
    ),
    error = function(e) { message(sprintf("  [ERROR %s] pairwise: %s", label, e$message)); NULL }
  )
  if (is.null(pw)) return(NULL)

  # Check connectivity; restrict to largest component containing REF_CLASS
  nc <- tryCatch(netconnection(pw), error = function(e) NULL)
  if (!is.null(nc) && nc$n.subnets > 1) {
    message(sprintf(
      "  [INFO %s] Network has %d sub-networks; restricting to component containing '%s'.",
      label, nc$n.subnets, REF_CLASS
    ))
    # Find studies in the component that contains REF_CLASS
    # netconnection stores study-to-subnet assignment in nc$subnet
    ref_subnet <- nc$subnet[nc$treat1 == REF_CLASS | nc$treat2 == REF_CLASS]
    ref_subnet <- unique(ref_subnet)
    if (length(ref_subnet) == 0) {
      message(sprintf("  [SKIP %s] Could not locate '%s' in any sub-network.", label, REF_CLASS))
      return(NULL)
    }
    ref_subnet <- ref_subnet[1]
    keep_studies <- unique(nc$studlab[nc$subnet == ref_subnet])
    pw <- pw[pw$studlab %in% keep_studies, ]
    if (nrow(pw) == 0) {
      message(sprintf("  [SKIP %s] No rows left after filtering to connected component.", label))
      return(NULL)
    }
  }

  nma_fit <- tryCatch(
    netmeta(
      TE      = pw$TE,
      seTE    = pw$seTE,
      treat1  = pw$treat1,
      treat2  = pw$treat2,
      studlab = pw$studlab,
      data    = pw,
      sm      = "SMD",
      random  = TRUE,
      common  = FALSE,
      reference.group = REF_CLASS
    ),
    error = function(e) { message(sprintf("  [ERROR %s] netmeta: %s", label, e$message)); NULL }
  )
  if (is.null(nma_fit)) return(NULL)

  # Extract estimates vs Placebo
  TE_mat    <- nma_fit$TE.random
  lower_mat <- nma_fit$lower.random
  upper_mat <- nma_fit$upper.random

  if (REF_CLASS %in% colnames(TE_mat)) {
    col <- REF_CLASS
    trt_names <- rownames(TE_mat)
    smd_vals  <- TE_mat[, col]
    lo_vals   <- lower_mat[, col]
    hi_vals   <- upper_mat[, col]
  } else {
    row <- REF_CLASS
    trt_names <- colnames(TE_mat)
    smd_vals  <- TE_mat[row, ]
    lo_vals   <- lower_mat[row, ]
    hi_vals   <- upper_mat[row, ]
  }

  tibble::tibble(
    class    = trt_names,
    smd      = as.numeric(smd_vals),
    lower    = as.numeric(lo_vals),
    upper    = as.numeric(hi_vals),
    data_subset = label,
    n_studies = n_distinct(collapsed$studyid),
    n_arms    = nrow(collapsed)
  )
}

# ------------------------------------------------------------------
# 4.  Run NMA on three subsets
# ------------------------------------------------------------------
message("\nRunning subset NMAs ...")
nma_cfb     <- run_nma(dat %>% filter(data_type == "CFB"),            "CFB_only")
nma_cfb_bf  <- run_nma(dat %>% filter(data_type %in% c("CFB","BF")), "CFB_plus_BF")
nma_all     <- run_nma(dat,                                           "ALL")

nma_compare <- dplyr::bind_rows(nma_cfb, nma_cfb_bf, nma_all)
write_csv(nma_compare, file.path(out_dir, "consistency_network_estimates.csv"))

# ------------------------------------------------------------------
# 5.  Scatter: ALL vs CFB_only estimates
# ------------------------------------------------------------------
message("\nBuilding comparison plots ...")

if (!is.null(nma_cfb) && !is.null(nma_all)) {
  scatter_dat <- nma_cfb %>%
    select(class, smd_cfb = smd) %>%
    inner_join(nma_all %>% select(class, smd_all = smd), by = "class")

  r_pearson <- cor(scatter_dat$smd_cfb, scatter_dat$smd_all, use = "complete.obs")

  p_scatter <- ggplot(scatter_dat, aes(x = smd_cfb, y = smd_all, label = class)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 2.5, colour = "#2c7bb6") +
    geom_text(size = 2.2, hjust = -0.1, vjust = 0.5, check_overlap = TRUE) +
    labs(
      title    = "NMA class estimates: CFB-only vs ALL data types",
      subtitle = sprintf("Pearson r = %.3f  (dashed = identity line)", r_pearson),
      x        = "SMD vs Placebo (CFB studies only)",
      y        = "SMD vs Placebo (all data types)"
    ) +
    theme_bw(base_size = 11)
} else {
  p_scatter <- ggplot() +
    labs(title = "Not enough data for CFB-only vs ALL scatter")
}

# ------------------------------------------------------------------
# 6.  Violin: distribution of within-arm approx SMD by data_type
# ------------------------------------------------------------------
p_violin <- ggplot(arm_smd, aes(x = data_type, y = approx_smd, fill = data_type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 1.2, fill = "white") +
  scale_fill_manual(values = c(CFB = "#1a9641", BF = "#fdae61", BIN = "#d7191c")) +
  labs(
    title = "Within-arm approx SMD (mean_change / sd_change) by data type",
    x     = "Data type",
    y     = "Approx SMD (mean_change / sd_change)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# Forest-style comparison: class estimates by subset
if (nrow(nma_compare) > 0) {
  top_classes <- nma_compare %>%
    filter(data_subset == "ALL", class != REF_CLASS) %>%
    arrange(smd) %>%
    slice_head(n = 20) %>%
    pull(class)

  p_forest <- nma_compare %>%
    filter(class %in% top_classes) %>%
    mutate(class = factor(class, levels = top_classes)) %>%
    ggplot(aes(x = smd, y = class, colour = data_subset,
               xmin = lower, xmax = upper)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.3, position = position_dodge(0.6), orientation = "y") +
    geom_point(size = 2, position = position_dodge(0.6)) +
    scale_colour_manual(
      values = c(CFB_only = "#1a9641", CFB_plus_BF = "#fdae61", ALL = "#d7191c"),
      name   = "Data subset"
    ) +
    labs(
      title = "Class SMD vs Placebo: sensitivity to data type inclusion",
      x     = "SMD vs Placebo (95% CI)",
      y     = NULL
    ) +
    theme_bw(base_size = 10)
} else {
  p_forest <- ggplot() + labs(title = "No NMA results to plot")
}

# Save all plots to one PDF
pdf(file.path(out_dir, "consistency_comparison_plots.pdf"), width = 12, height = 7)
print(p_violin)
if (!is.null(nma_cfb) && !is.null(nma_all)) print(p_scatter)
print(p_forest)
dev.off()

# ------------------------------------------------------------------
# 7.  Meta-regression: test data_type as covariate
#     Uses pairwise contrasts; data_type_BF and data_type_BIN are
#     indicators (reference = CFB).  A significant coefficient would
#     suggest the conversion introduces systematic bias.
# ------------------------------------------------------------------
message("\nRunning data_type meta-regression ...")

dat_class_full <- dat %>%
  group_by(studyid, classcode, class, data_type) %>%
  summarise(
    n           = sum(n),
    mean_change = sum(mean_change * n) / sum(n),
    sd_change   = sqrt(sum((n - 1) * sd_change^2) / sum(n - 1)),
    .groups     = "drop"
  ) %>%
  group_by(studyid) %>%
  filter(n_distinct(class) >= 2) %>%
  ungroup()

pw_full <- tryCatch(
  pairwise(
    treat   = class,
    mean    = mean_change,
    sd      = sd_change,
    n       = n,
    studlab = studyid,
    data    = dat_class_full,
    sm      = "SMD"
  ),
  error = function(e) { message("  pairwise failed: ", e$message); NULL }
)

if (!is.null(pw_full)) {
  # Attach data_type to pairwise contrasts (from study-level flag).
  # pairwise() drops extra columns, so we join back from the original data.
  # Each study belongs to exactly one block, so data_type is study-level.
  study_type <- dat %>%
    distinct(studyid, data_type) %>%
    group_by(studyid) %>% slice(1) %>% ungroup() %>%
    rename(study_data_type = data_type)   # avoid name collision with any pw_full column

  pw_full <- pw_full %>%
    left_join(study_type, by = c("studlab" = "studyid")) %>%
    mutate(
      type_BF  = as.integer(study_data_type == "BF"),
      type_BIN = as.integer(study_data_type == "BIN"),
      # Replace NA (studies not matched) with 0 = CFB reference
      type_BF  = ifelse(is.na(type_BF),  0L, type_BF),
      type_BIN = ifelse(is.na(type_BIN), 0L, type_BIN)
    )

  # Fit metareg (uses metafor internally via netmeta's metareg wrapper,
  # or directly via metafor if preferred)
  if (requireNamespace("metafor", quietly = TRUE)) {
    mr_fit <- metafor::rma(
      yi   = pw_full$TE,
      sei  = pw_full$seTE,
      mods = ~ type_BF + type_BIN,
      data = pw_full,
      method = "REML"
    )
    mr_coefs <- data.frame(
      parameter = rownames(coef(summary(mr_fit))),
      estimate  = coef(summary(mr_fit))[, "estimate"],
      se        = coef(summary(mr_fit))[, "se"],
      pval      = coef(summary(mr_fit))[, "pval"],
      row.names = NULL
    )
    message("\n--- Data-type meta-regression coefficients ---")
    print(mr_coefs)
    write_csv(mr_coefs, file.path(out_dir, "consistency_datatype_metareg.csv"))

    if (abs(mr_coefs$estimate[mr_coefs$parameter == "type_BF"])  > 0.2 |
        abs(mr_coefs$estimate[mr_coefs$parameter == "type_BIN"]) > 0.2) {
      message("  [NOTE] Coefficient(s) > 0.2 SMD: consider sensitivity analysis",
              " restricting to CFB-only studies.")
    } else {
      message("  [OK] Data-type coefficients < 0.2 SMD: conversion appears consistent.")
    }
  } else {
    message("  metafor not available – skipping meta-regression.",
            " Install with: install.packages('metafor')")
  }
}

message(sprintf("\nConsistency check complete. Outputs in: %s", out_dir))
