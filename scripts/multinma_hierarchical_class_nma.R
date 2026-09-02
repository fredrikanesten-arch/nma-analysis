library(dplyr)
library(readr)

if (!requireNamespace("multinma", quietly = TRUE)) {
  stop("Package 'multinma' is required. Install with install.packages('multinma').")
}

# Fixed paths requested by user
base_dir <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli"
in_data  <- file.path(base_dir, "binfixed_class_ms", "combined_long_mean_change_dataset_ms_smd_bias_adj.csv")
in_map   <- file.path(base_dir, "clean_data", "trt_to_class_ms.csv")
out_dir  <- file.path(base_dir, "binfixed_class_ms")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

reference_class <- "Placebo"
chains <- 4
iter <- 8000
warmup <- 4000
seed <- 20260902
adapt_delta <- 0.99
max_treedepth <- 15

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

build_fallback_map <- function() {
  data.frame(trtcode = 1:99) %>%
    mutate(
      classcode = case_when(
        trtcode == 1 ~ 1, trtcode == 2 ~ 2, trtcode == 3 ~ 3, trtcode == 4 ~ 4, trtcode == 5 ~ 5,
        trtcode == 6 ~ 6, trtcode == 7 ~ 7, trtcode %in% c(8, 9) ~ 8,
        trtcode %in% c(10, 11, 12, 13) ~ 9, trtcode == 14 ~ 10, trtcode == 15 ~ 11, trtcode == 16 ~ 12,
        trtcode == 17 ~ 13, trtcode == 18 ~ 14, trtcode == 19 ~ 15,
        trtcode %in% c(20, 21, 22, 23) ~ 16, trtcode %in% c(24, 25, 26, 27) ~ 17,
        trtcode %in% c(28, 29) ~ 18, trtcode == 30 ~ 19, trtcode == 31 ~ 20, trtcode == 32 ~ 21,
        trtcode == 33 ~ 22, trtcode %in% c(34, 35) ~ 23, trtcode == 36 ~ 24, trtcode == 37 ~ 25,
        trtcode == 38 ~ 26, trtcode %in% c(39, 40, 41, 42, 43, 44) ~ 27,
        trtcode %in% c(45, 46, 47, 48, 49, 50) ~ 28, trtcode %in% c(51, 52) ~ 29,
        trtcode == 53 ~ 30, trtcode %in% c(54, 55, 56) ~ 31, trtcode %in% c(57, 58, 59) ~ 32,
        trtcode %in% c(60, 61, 62) ~ 33, trtcode %in% c(63, 64) ~ 34, trtcode == 65 ~ 35,
        trtcode == 66 ~ 36, trtcode %in% c(67, 68) ~ 37,
        trtcode %in% c(69, 70, 71, 72, 73, 74, 75) ~ 38, trtcode == 76 ~ 39,
        trtcode %in% c(77, 78) ~ 40, trtcode %in% c(79, 80, 81) ~ 41, trtcode %in% c(82, 83) ~ 42,
        trtcode == 84 ~ 43, trtcode == 85 ~ 44, trtcode == 86 ~ 45,
        trtcode %in% c(87, 88, 89) ~ 46, trtcode %in% c(90, 91) ~ 47, trtcode == 92 ~ 48,
        trtcode %in% c(93, 94, 95, 96, 97) ~ 49, trtcode %in% c(98, 99) ~ 50,
        TRUE ~ NA_real_
      ),
      class = case_when(
        trtcode == 1 ~ "Placebo",
        trtcode == 2 ~ "Attention placebo",
        trtcode == 3 ~ "No treatment",
        trtcode == 4 ~ "Waitlist",
        trtcode == 5 ~ "TAU",
        trtcode == 6 ~ "Mirtazapine",
        trtcode == 7 ~ "Trazodone",
        trtcode %in% c(8, 9) ~ "Behavioural therapies individual",
        trtcode %in% c(10, 11, 12, 13) ~ "Cognitive and cognitive behavioural therapies individual",
        trtcode == 14 ~ "Cognitive and cognitive behavioural therapies group",
        trtcode == 15 ~ "Problem solving individual",
        trtcode == 16 ~ "Problem solving group",
        trtcode == 17 ~ "Counselling individual",
        trtcode == 18 ~ "Interpersonal psychotherapy (IPT) individual",
        trtcode == 19 ~ "Psychoeducation group",
        trtcode %in% c(20, 21, 22, 23) ~ "Self-help",
        trtcode %in% c(24, 25, 26, 27) ~ "Self-help with support",
        trtcode %in% c(28, 29) ~ "Short-term psychodynamic psychotherapies individual",
        trtcode == 30 ~ "Music therapy group",
        trtcode == 31 ~ "Mindfulness or meditation group",
        trtcode == 32 ~ "Peer support group",
        trtcode == 33 ~ "Any psychotherapy",
        trtcode %in% c(34, 35) ~ "Cognitive and cognitive behavioural therapies individual + placebo",
        trtcode == 36 ~ "Interpersonal psychotherapy (IPT) individual + placebo",
        trtcode == 37 ~ "Counselling individual + placebo",
        trtcode == 38 ~ "Relaxation individual + placebo",
        trtcode %in% c(39, 40, 41, 42, 43, 44) ~ "SSRIs",
        trtcode %in% c(45, 46, 47, 48, 49, 50) ~ "TCAs",
        trtcode %in% c(51, 52) ~ "SNRIs",
        trtcode == 53 ~ "Any AD",
        trtcode %in% c(54, 55, 56) ~ "Sham acupuncture",
        trtcode %in% c(57, 58, 59) ~ "Acupuncture",
        trtcode %in% c(60, 61, 62) ~ "Exercise individual",
        trtcode %in% c(63, 64) ~ "Exercise group",
        trtcode == 65 ~ "Yoga group",
        trtcode == 66 ~ "Light therapy",
        trtcode %in% c(67, 68) ~ "Behavioural therapies individual + AD",
        trtcode %in% c(69, 70, 71, 72, 73, 74, 75) ~ "Cognitive and cognitive behavioural therapies individual + AD",
        trtcode == 76 ~ "Cognitive and cognitive behavioural therapies group + AD",
        trtcode %in% c(77, 78) ~ "Interpersonal psychotherapy (IPT) individual + AD",
        trtcode %in% c(79, 80, 81) ~ "Counselling individual + AD",
        trtcode %in% c(82, 83) ~ "Short-term psychodynamic psychotherapies individual + AD",
        trtcode == 84 ~ "Psychoeducation group + AD",
        trtcode == 85 ~ "Peer support group + AD",
        trtcode == 86 ~ "Relaxation individual + AD",
        trtcode %in% c(87, 88, 89) ~ "Exercise individual + AD",
        trtcode %in% c(90, 91) ~ "Exercise group + AD",
        trtcode == 92 ~ "Yoga group + AD",
        trtcode %in% c(93, 94, 95, 96, 97) ~ "Acupuncture + AD",
        trtcode %in% c(98, 99) ~ "Light therapy + AD",
        TRUE ~ NA_character_
      ),
      treatment_label = sprintf("trt_%03d", trtcode)
    ) %>%
    filter(!is.na(classcode), !is.na(class), nzchar(class))
}

dat <- read_csv_robust(in_data) %>%
  mutate(
    treatment = as.numeric(treatment),
    n = as.numeric(n),
    mean_change = as.numeric(mean_change),
    sd_change = as.numeric(sd_change)
  )

trt_map <- read_csv_robust(in_map) %>%
  mutate(
    trtcode = as.numeric(trtcode),
    classcode = as.numeric(classcode),
    class = as.character(class)
  )

required_data_cols <- c("studyid", "treatment", "n", "mean_change", "sd_change")
required_map_cols <- c("trtcode", "classcode", "class")
missing_data_cols <- setdiff(required_data_cols, names(dat))
missing_map_cols <- setdiff(required_map_cols, names(trt_map))
if (length(missing_data_cols) > 0) stop("Missing required data columns: ", paste(missing_data_cols, collapse = ", "))
if (length(missing_map_cols) > 0) stop("Missing required map columns: ", paste(missing_map_cols, collapse = ", "))

name_col_candidates <- c("treatment_name", "trt_name", "treatment", "trt")
name_col <- name_col_candidates[name_col_candidates %in% names(trt_map)][1]
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

trt_map_fallback <- build_fallback_map()
trt_map_dedup <- bind_rows(
  trt_map_primary,
  trt_map_fallback %>% anti_join(trt_map_primary %>% select(trtcode), by = "trtcode")
)

dat_joined <- dat %>%
  left_join(trt_map_dedup %>% select(treatment = trtcode, classcode, class, treatment_label), by = "treatment") %>%
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
      TRUE ~ "kept"
    )
  )

write_csv(count(dat_joined, drop_reason, sort = TRUE), file.path(out_dir, "multinma_drop_reason_counts.csv"))

unmapped <- dat %>%
  distinct(treatment) %>%
  left_join(trt_map_dedup %>% select(treatment = trtcode, classcode, class), by = "treatment") %>%
  filter(is.na(classcode) | is.na(class))
write_csv(unmapped, file.path(out_dir, "multinma_unmapped_treatments.csv"))

dat_kept <- dat_joined %>%
  filter(drop_reason == "kept")

# Collapse only exact duplicate treatment arms in each study.
arm_dat <- dat_kept %>%
  group_by(studyid, treatment_label, classcode, class) %>%
  summarise(
    n_total = sum(n),
    mean_change_total = sum(mean_change * n) / sum(n),
    sd_change_total = pool_arm_sd(n = n, mean = mean_change, sd = sd_change),
    .groups = "drop"
  ) %>%
  transmute(
    studyid = studyid,
    treatment_label = treatment_label,
    class = class,
    classcode = classcode,
    n = n_total,
    mean_change = mean_change_total,
    sd_change = sd_change_total,
    se_change = sd_change_total / sqrt(n_total)
  )

write_csv(arm_dat, file.path(out_dir, "multinma_treatment_level_arms.csv"))
study_trt_counts <- arm_dat %>%
  group_by(studyid) %>%
  summarise(n_treatments = n_distinct(treatment_label), .groups = "drop")

write_csv(study_trt_counts, file.path(out_dir, "multinma_study_treatment_counts.csv"))

single_arm_studies <- study_trt_counts %>%
  filter(n_treatments < 2)
write_csv(single_arm_studies, file.path(out_dir, "multinma_single_arm_studies_dropped.csv"))

arm_dat_nma <- arm_dat %>%
  filter(studyid %in% (study_trt_counts %>% filter(n_treatments >= 2) %>% pull(studyid)))

if (nrow(arm_dat_nma) == 0) {
  stop("No multi-arm studies available after removing single-arm studies.")
}

if (!reference_class %in% arm_dat_nma$class) {
  stop("Reference class not present in retained data: ", reference_class)
}
reference_treatment <- arm_dat_nma %>%
  filter(class == reference_class) %>%
  arrange(treatment_label) %>%
  pull(treatment_label) %>%
  unique() %>%
  .[[1]]

network <- multinma::set_agd_arm(
  arm_dat_nma,
  study = studyid,
  trt = treatment_label,
  y = mean_change,
  se = se_change,
  trt_class = class
)

fit <- multinma::nma(
  network,
  trt_effects = "random",
  class_effects = "exchangeable",
  likelihood = "normal",
  link = "identity",
  prior_intercept = multinma::normal(location = 0, scale = 10),
  prior_trt = multinma::normal(location = 0, scale = 5),
  prior_het = multinma::half_normal(scale = 1),
  prior_class_mean = multinma::normal(location = 0, scale = 5),
  prior_class_sd = multinma::half_normal(scale = 1),
  chains = chains,
  iter = iter,
  warmup = warmup,
  seed = seed,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth
)

rel <- multinma::relative_effects(fit, trt_ref = reference_treatment)
rel_sum <- tryCatch(
  as.data.frame(summary(rel)),
  error = function(e) {
    stop(
      "Could not summarize relative effects with summary(). ",
      "Please check multinma version and object structure. Original error: ", e$message
    )
  }
)
if (!("Treatment" %in% names(rel_sum))) {
  rel_sum$Treatment <- rownames(rel_sum)
}
if (!("mean" %in% names(rel_sum))) {
  mean_col <- intersect(c("mean", "Mean", "Estimate", "estimate", "median", "Median"), names(rel_sum))[1]
  if (!is.na(mean_col)) rel_sum <- rel_sum %>% rename(mean = all_of(mean_col))
}
if (!("sd" %in% names(rel_sum))) {
  sd_col <- intersect(c("sd", "SD", "Est.Error", "est.error", "se", "SE"), names(rel_sum))[1]
  if (!is.na(sd_col)) rel_sum <- rel_sum %>% rename(sd = all_of(sd_col))
}
if (!all(c("mean", "sd") %in% names(rel_sum))) {
  stop("Could not identify mean/sd columns in relative effects summary.")
}
rownames(rel_sum) <- NULL

treatment_n <- arm_dat_nma %>%
  group_by(treatment_label, class) %>%
  summarise(N_total_treatment = sum(n), n_studies_treatment = n_distinct(studyid), .groups = "drop") %>%
  rename(Treatment = treatment_label)

treatment_vs_ref <- rel_sum %>%
  left_join(treatment_n, by = "Treatment") %>%
  arrange(class, Treatment)

class_vs_ref <- treatment_vs_ref %>%
  group_by(class) %>%
  summarise(
    n_treatments = n_distinct(Treatment),
    N_total_class = sum(N_total_treatment, na.rm = TRUE),
    n_studies_class = sum(n_studies_treatment, na.rm = TRUE),
    mean_class_vs_ref = {
      idx <- which(!is.na(mean) & !is.na(sd) & sd > 0)
      if (length(idx) == 0) NA_real_ else weighted.mean(mean[idx], w = 1 / (sd[idx]^2))
    },
    .groups = "drop"
  ) %>%
  arrange(class)

write_csv(treatment_vs_ref, file.path(out_dir, "multinma_treatment_vs_reference_summary.csv"))
write_csv(class_vs_ref, file.path(out_dir, "multinma_class_summary_from_treatment_posteriors.csv"))

fit_sum <- capture.output(print(multinma::summary(fit)))
writeLines(fit_sum, con = file.path(out_dir, "multinma_model_summary.txt"))
diag_sum <- capture.output(print(summary(fit)))
writeLines(diag_sum, con = file.path(out_dir, "multinma_sampling_diagnostics.txt"))

saveRDS(
  list(
    arm_dat = arm_dat,
    arm_dat_nma = arm_dat_nma,
    study_treatment_counts = study_trt_counts,
    single_arm_studies = single_arm_studies,
    network = network,
    fit = fit,
    rel = rel,
    treatment_vs_ref = treatment_vs_ref,
    class_vs_ref = class_vs_ref,
    reference_treatment = reference_treatment
  ),
  file = file.path(out_dir, "multinma_hierarchical_class_model_objects.rds")
)

message("Done.")
message("Reference treatment: ", reference_treatment)
message("Outputs written to: ", normalizePath(out_dir))
