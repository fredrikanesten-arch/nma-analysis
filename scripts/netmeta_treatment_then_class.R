library(dplyr)
library(readr)
library(netmeta)

# Usage:
# Rscript scripts/netmeta_treatment_then_class.R <combined_long_csv> <trt_to_class_csv> <output_dir> [reference_treatment]

args <- commandArgs(trailingOnly = TRUE)
in_data <- if (length(args) >= 1) args[[1]] else "combined_long_mean_change_dataset_ms_smd_bias_adj.csv"
in_map <- if (length(args) >= 2) args[[2]] else "trt_to_class_ms.csv"
out_dir <- if (length(args) >= 3) args[[3]] else "netmeta_treatment_then_class_outputs"
reference_treatment_arg <- if (length(args) >= 4) args[[4]] else NA_character_

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

read_csv_robust <- function(path) {
  dat_csv <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(dat_csv) && ncol(dat_csv) > 1) return(dat_csv)

  dat_scsv <- tryCatch(readr::read_delim(path, delim = ";", show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(dat_scsv) && ncol(dat_scsv) > 1) return(dat_scsv)

  stop("Could not parse file: ", path)
}

pool_arm_sd <- function(n, mean, sd) {
  n_total <- sum(n)
  mean_total <- sum(mean * n) / n_total
  ss_within <- sum((n - 1) * sd^2)
  ss_between <- sum(n * (mean - mean_total)^2)
  sqrt((ss_within + ss_between) / (n_total - 1))
}

extract_vs_ref <- function(mat, ref) {
  if (is.null(mat)) return(NULL)
  if (ref %in% colnames(mat)) {
    return(data.frame(Treatment = rownames(mat), value = as.numeric(mat[, ref]), row.names = NULL))
  }
  if (ref %in% rownames(mat)) {
    return(data.frame(Treatment = colnames(mat), value = as.numeric(mat[ref, ]), row.names = NULL))
  }
  stop("Reference treatment not found in matrix dimnames.")
}

dat <- read_csv_robust(in_data)
trt_map <- read_csv_robust(in_map)

required_data_cols <- c("studyid", "treatment", "n", "mean_change", "sd_change")
required_map_cols <- c("trtcode", "classcode", "class")
missing_data_cols <- setdiff(required_data_cols, names(dat))
missing_map_cols <- setdiff(required_map_cols, names(trt_map))
if (length(missing_data_cols) > 0) stop("Missing required data columns: ", paste(missing_data_cols, collapse = ", "))
if (length(missing_map_cols) > 0) stop("Missing required map columns: ", paste(missing_map_cols, collapse = ", "))

name_col_candidates <- c("treatment_name", "trt_name", "treatment", "trt")
name_col <- name_col_candidates[name_col_candidates %in% names(trt_map)][1]

dat <- dat %>%
  mutate(
    treatment = as.numeric(treatment),
    n = as.numeric(n),
    mean_change = as.numeric(mean_change),
    sd_change = as.numeric(sd_change)
  )

trt_map <- trt_map %>%
  mutate(
    trtcode = as.numeric(trtcode),
    classcode = as.numeric(classcode),
    class = as.character(class)
  )

if (!is.na(name_col)) {
  trt_map <- trt_map %>% mutate(treatment_label = as.character(.data[[name_col]]))
} else {
  trt_map <- trt_map %>% mutate(treatment_label = sprintf("trt_%03d", trtcode))
}

trt_map <- trt_map %>%
  mutate(
    treatment_label = ifelse(is.na(treatment_label) | !nzchar(treatment_label),
                             sprintf("trt_%03d", trtcode),
                             treatment_label)
  )

trt_map_primary <- trt_map %>%
  filter(!is.na(trtcode), !is.na(classcode), !is.na(class), nzchar(class), !is.na(treatment_label), nzchar(treatment_label)) %>%
  arrange(trtcode, classcode, class, treatment_label) %>%
  group_by(trtcode) %>%
  slice(1) %>%
  ungroup()

dat_joined <- dat %>%
  left_join(trt_map_primary %>% select(treatment = trtcode, classcode, class, treatment_label), by = "treatment")

dat_audit <- dat_joined %>%
  mutate(
    drop_reason = case_when(
      is.na(studyid) ~ "missing_studyid",
      is.na(treatment) ~ "missing_treatment",
      is.na(classcode) | is.na(class) ~ "unmapped_class",
      is.na(n) ~ "missing_n",
      is.na(mean_change) ~ "missing_mean_change",
      is.na(sd_change) ~ "missing_sd_change",
      n <= 1 ~ "invalid_n",
      sd_change <= 0 ~ "non_positive_sd_change",
      TRUE ~ "kept_pre_collapse"
    )
  )

write_csv(count(dat_audit, drop_reason, sort = TRUE), file.path(out_dir, "drop_reason_counts_treatment_first.csv"))

dat_kept <- dat_audit %>%
  filter(drop_reason == "kept_pre_collapse")

# Collapse only exact duplicate treatment arms in a study (not class-level collapse).
collapsed_treatment <- dat_kept %>%
  group_by(studyid, treatment_label, classcode, class) %>%
  summarise(
    n = sum(n),
    mean_change = sum(mean_change * n) / sum(n),
    sd_change = pool_arm_sd(n, mean_change, sd_change),
    .groups = "drop"
  )

study_level_after_collapse <- collapsed_treatment %>%
  group_by(studyid) %>%
  summarise(
    n_treatments = n_distinct(treatment_label),
    .groups = "drop"
  )

write_csv(study_level_after_collapse, file.path(out_dir, "study_treatment_counts_after_collapse.csv"))

collapsed_treatment_nma <- collapsed_treatment %>%
  group_by(studyid) %>%
  filter(n_distinct(treatment_label) >= 2) %>%
  ungroup()

write_csv(collapsed_treatment_nma, file.path(out_dir, "treatment_level_arms_used_for_nma.csv"))

if (nrow(collapsed_treatment_nma) == 0) {
  stop("No analyzable treatment-level arms available after filtering.")
}

if (!is.na(reference_treatment_arg)) {
  reference_treatment <- reference_treatment_arg
} else {
  placebo_candidates <- trt_map_primary %>%
    filter(tolower(class) == "placebo") %>%
    pull(treatment_label) %>%
    unique()
  if (length(placebo_candidates) >= 1) {
    reference_treatment <- placebo_candidates[[1]]
  } else {
    reference_treatment <- unique(collapsed_treatment_nma$treatment_label)[[1]]
  }
}

if (!(reference_treatment %in% collapsed_treatment_nma$treatment_label)) {
  stop("Reference treatment not present in analyzed data: ", reference_treatment)
}

pw <- pairwise(
  treat = treatment_label,
  mean = mean_change,
  sd = sd_change,
  n = n,
  studlab = studyid,
  data = collapsed_treatment_nma,
  sm = "MD"
)

nma <- netmeta(
  TE = pw$TE,
  seTE = pw$seTE,
  treat1 = pw$treat1,
  treat2 = pw$treat2,
  studlab = pw$studlab,
  data = pw,
  sm = "MD",
  random = TRUE,
  common = FALSE,
  reference.group = reference_treatment
)

TE_mat <- if (!is.null(nma$TE.random)) nma$TE.random else nma$TE.common
lower_mat <- if (!is.null(nma$lower.random)) nma$lower.random else nma$lower.common
upper_mat <- if (!is.null(nma$upper.random)) nma$upper.random else nma$upper.common
se_mat <- if (!is.null(nma$seTE.random)) nma$seTE.random else nma$seTE.common

trt_eff <- extract_vs_ref(TE_mat, reference_treatment) %>% rename(SMD = value)
trt_low <- extract_vs_ref(lower_mat, reference_treatment) %>% rename(lower_CI = value)
trt_up <- extract_vs_ref(upper_mat, reference_treatment) %>% rename(upper_CI = value)
trt_se <- extract_vs_ref(se_mat, reference_treatment) %>% rename(seTE = value)

treatment_n_totals <- collapsed_treatment_nma %>%
  group_by(treatment_label) %>%
  summarise(N_total_treatment = sum(n), n_studies_treatment = n_distinct(studyid), .groups = "drop") %>%
  rename(Treatment = treatment_label)

treatment_class_lookup <- trt_map_primary %>%
  distinct(treatment_label, classcode, class) %>%
  rename(Treatment = treatment_label)

treatment_vs_ref <- trt_eff %>%
  left_join(trt_low, by = "Treatment") %>%
  left_join(trt_up, by = "Treatment") %>%
  left_join(trt_se, by = "Treatment") %>%
  left_join(treatment_n_totals, by = "Treatment") %>%
  left_join(treatment_class_lookup, by = "Treatment") %>%
  mutate(
    SMD = ifelse(Treatment == reference_treatment, 0, SMD),
    lower_CI = ifelse(Treatment == reference_treatment, 0, lower_CI),
    upper_CI = ifelse(Treatment == reference_treatment, 0, upper_CI),
    seTE = ifelse(Treatment == reference_treatment, NA_real_, seTE),
    z = ifelse(!is.na(seTE) & seTE > 0, SMD / seTE, NA_real_),
    p_value = ifelse(!is.na(z), 2 * stats::pnorm(-abs(z)), NA_real_)
  ) %>%
  arrange(classcode, Treatment)

class_vs_ref <- treatment_vs_ref %>%
  group_by(classcode, class) %>%
  summarise(
    n_treatments = n_distinct(Treatment),
    N_total_class = sum(N_total_treatment, na.rm = TRUE),
    n_studies_class = sum(n_studies_treatment, na.rm = TRUE),
    SMD_class_vs_ref = {
      idx <- which(!is.na(SMD) & !is.na(seTE) & seTE > 0)
      if (length(idx) == 0) {
        if (all(Treatment == reference_treatment)) 0 else NA_real_
      } else {
        weighted.mean(SMD[idx], w = 1 / (seTE[idx]^2))
      }
    },
    seTE_class_vs_ref = {
      idx <- which(!is.na(SMD) & !is.na(seTE) & seTE > 0)
      if (length(idx) == 0) NA_real_ else sqrt(1 / sum(1 / (seTE[idx]^2)))
    },
    .groups = "drop"
  ) %>%
  mutate(
    lower_CI = ifelse(!is.na(seTE_class_vs_ref), SMD_class_vs_ref - 1.96 * seTE_class_vs_ref, NA_real_),
    upper_CI = ifelse(!is.na(seTE_class_vs_ref), SMD_class_vs_ref + 1.96 * seTE_class_vs_ref, NA_real_),
    z = ifelse(!is.na(seTE_class_vs_ref) & seTE_class_vs_ref > 0, SMD_class_vs_ref / seTE_class_vs_ref, NA_real_),
    p_value = ifelse(!is.na(z), 2 * stats::pnorm(-abs(z)), NA_real_)
  ) %>%
  arrange(classcode)

write_csv(treatment_vs_ref, file.path(out_dir, "nma_treatment_level_vs_reference.csv"))
write_csv(class_vs_ref, file.path(out_dir, "nma_class_level_aggregated_from_treatments.csv"))
write_csv(as.data.frame(pw), file.path(out_dir, "pairwise_treatment_level.csv"))

saveRDS(
  list(
    input_data = dat,
    treatment_class_map = trt_map_primary,
    dat_audit = dat_audit,
    collapsed_treatment = collapsed_treatment,
    collapsed_treatment_nma = collapsed_treatment_nma,
    pairwise = pw,
    nma = nma,
    treatment_vs_ref = treatment_vs_ref,
    class_vs_ref = class_vs_ref,
    reference_treatment = reference_treatment
  ),
  file = file.path(out_dir, "treatment_first_nma_objects.rds")
)

message("Done. Reference treatment: ", reference_treatment)
message("Outputs written to: ", normalizePath(out_dir))
