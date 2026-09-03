if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required. Install with: install.packages('dplyr')")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required. Install with: install.packages('readr')")
}
if (!requireNamespace("multinma", quietly = TRUE)) {
  stop("Package 'multinma' is required. Install with: install.packages('multinma')")
}

library(dplyr)
library(readr)
library(multinma)

# -----------------------------
# 0. Paths
# -----------------------------
base_dir <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli"
in_data <- file.path(base_dir, "binfixed_class_ms", "combined_long_mean_change_dataset_ms_smd_bias_adj.csv")
in_map <- file.path(base_dir, "clean_data", "trt_to_class_ms.csv")
out_dir <- file.path(base_dir, "binfixed_class_ms")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 1. Robust CSV reader
# -----------------------------
read_combined_long_robust <- function(path) {
  lines <- readr::read_lines(path, n_max = 5)

  dat_csv <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (!is.null(dat_csv) && ncol(dat_csv) > 1) {
    message("Read data as standard comma-separated CSV.")
    return(dat_csv)
  }

  dat_scsv <- tryCatch(
    readr::read_delim(path, delim = ";", show_col_types = FALSE),
    error = function(e) NULL
  )
  if (!is.null(dat_scsv) && ncol(dat_scsv) > 1) {
    message("Read data as semicolon-separated CSV.")
    return(dat_scsv)
  }

  repaired_lines <- readr::read_lines(path)
  repaired_lines <- sub('^"', "", repaired_lines)
  repaired_lines <- sub('"$', "", repaired_lines)
  repaired_lines <- gsub('""', '"', repaired_lines, fixed = TRUE)

  dat_repaired <- tryCatch(
    readr::read_csv(I(paste(repaired_lines, collapse = "\n")), show_col_types = FALSE),
    error = function(e) NULL
  )
  if (!is.null(dat_repaired) && ncol(dat_repaired) > 1) {
    message("Read data after repairing double-encoded CSV format.")
    return(dat_repaired)
  }

  stop(
    "Could not parse input data file. First lines were:\n",
    paste(lines, collapse = "\n")
  )
}

# -----------------------------
# 2. Read data and convert types
# -----------------------------
dat <- read_combined_long_robust(in_data)

trt_map <- read_delim(
  in_map,
  delim = ";",
  col_types = cols(.default = col_character())
)

dat <- dat %>%
  mutate(
    across(any_of(c("na", "arm", "treatment", "n", "mean_change", "sd_change")), as.numeric)
  )

trt_map <- trt_map %>%
  mutate(
    trtcode = as.numeric(trtcode),
    classcode = as.numeric(classcode)
  )

# -----------------------------
# 3. Deterministic treatment -> class mapping
# -----------------------------
trt_map <- trt_map %>%
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
    )
  )

# -----------------------------
# 4. Join + keep analyzable rows
# -----------------------------
dat_class <- dat %>%
  left_join(
    trt_map %>% select(treatment = trtcode, classcode, class),
    by = "treatment"
  ) %>%
  filter(
    !is.na(studyid),
    !is.na(classcode),
    !is.na(class),
    !is.na(n),
    !is.na(mean_change),
    !is.na(sd_change)
  )

# -----------------------------
# 5. Collapse duplicate same-class arms within study
# -----------------------------
collapsed_class <- dat_class %>%
  group_by(studyid, classcode, class) %>%
  summarise(
    n = sum(n),
    mean_change = sum(mean_change * n) / sum(n),
    sd_change = sqrt(sum((n - 1) * sd_change^2) / sum(n - 1)),
    .groups = "drop"
  ) %>%
  group_by(studyid) %>%
  filter(n_distinct(class) >= 2) %>%
  ungroup() %>%
  mutate(
    se_change = sd_change / sqrt(n)
  ) %>%
  filter(is.finite(se_change), se_change > 0)

stopifnot("Placebo" %in% collapsed_class$class)

# -----------------------------
# 6. Bayesian NMA (multinma)
# -----------------------------
network <- multinma::set_agd_arm(
  data = collapsed_class,
  study = studyid,
  trt = class,
  y = mean_change,
  se = se_change
)

fit <- multinma::nma(
  network,
  trt_effects = "random",
  consistency = "consistency",
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 20260903,
  refresh = 0
)

# -----------------------------
# 7. Summary table vs Placebo
# -----------------------------
extract_vs_placebo <- function(fit, collapsed_class, ref = "Placebo", digits = 3) {
  rel <- multinma::relative_effects(fit, trt_ref = ref)
  rel_sum <- as.data.frame(summary(rel))
  rel_sum <- tibble::as_tibble(rel_sum)

  n_tbl <- collapsed_class %>%
    group_by(class) %>%
    summarise(N = sum(n, na.rm = TRUE), .groups = "drop")

  # Try to support common multinma summary column names
  trt_col <- intersect(c("trt", "treatment", "class"), names(rel_sum))[1]
  mean_col <- intersect(c("mean", "Estimate", "estimate"), names(rel_sum))[1]
  lcl_col <- grep("^2\\.5%$|^lcl$|lower", names(rel_sum), value = TRUE, ignore.case = TRUE)[1]
  ucl_col <- grep("^97\\.5%$|^ucl$|upper", names(rel_sum), value = TRUE, ignore.case = TRUE)[1]

  if (any(is.na(c(trt_col, mean_col, lcl_col, ucl_col)))) {
    stop("Could not detect expected columns in relative_effects() summary output.")
  }

  out <- rel_sum %>%
    transmute(
      Treatment = .data[[trt_col]],
      SMD = as.numeric(.data[[mean_col]]),
      `lower CrI` = as.numeric(.data[[lcl_col]]),
      `upper CrI` = as.numeric(.data[[ucl_col]])
    ) %>%
    right_join(tibble::tibble(Treatment = sort(unique(collapsed_class$class))), by = "Treatment") %>%
    mutate(
      SMD = ifelse(Treatment == ref & is.na(SMD), 0, SMD),
      `lower CrI` = ifelse(Treatment == ref & is.na(`lower CrI`), 0, `lower CrI`),
      `upper CrI` = ifelse(Treatment == ref & is.na(`upper CrI`), 0, `upper CrI`)
    ) %>%
    left_join(n_tbl %>% rename(Treatment = class), by = "Treatment") %>%
    mutate(
      `SMD vs Placebo (mean, 95% CrI)` = ifelse(
        is.na(SMD),
        NA_character_,
        paste0(
          formatC(SMD, digits = digits, format = "f"),
          " (",
          formatC(`lower CrI`, digits = digits, format = "f"),
          " to ",
          formatC(`upper CrI`, digits = digits, format = "f"),
          ")"
        )
      )
    ) %>%
    select(Treatment, N, SMD, `lower CrI`, `upper CrI`, `SMD vs Placebo (mean, 95% CrI)`) %>%
    arrange(SMD)

  out
}

summary_table_full <- extract_vs_placebo(fit, collapsed_class, ref = "Placebo")

# Optional ranking summaries if available
rank_summary <- tryCatch({
  rank_probs <- multinma::posterior_ranks(fit, lower_better = TRUE)
  as.data.frame(summary(rank_probs)) %>%
    tibble::rownames_to_column("Treatment")
}, error = function(e) {
  tibble::tibble(note = paste("Ranking summary not produced:", conditionMessage(e)))
})

# -----------------------------
# 8. Save key outputs
# -----------------------------
write_csv(summary_table_full, file.path(out_dir, "bayes_nma_summary_table_full_multinma.csv"))
write_csv(rank_summary, file.path(out_dir, "bayes_nma_rank_summary_multinma.csv"))

saveRDS(
  list(
    dat = dat,
    trt_map = trt_map,
    collapsed_class = collapsed_class,
    network = network,
    fit = fit
  ),
  file = file.path(out_dir, "bayes_nma_class_level_analysis_objects_multinma.rds")
)

# -----------------------------
# 9. Class mapping consistency checks
# -----------------------------
dat_join_check <- dat %>%
  left_join(trt_map %>% select(treatment = trtcode, classcode, class), by = "treatment")
write_csv(dat_join_check, file.path(out_dir, "bayes_dat_join_check.csv"))

unmapped <- dat %>%
  distinct(treatment) %>%
  left_join(trt_map %>% select(treatment = trtcode, classcode, class), by = "treatment") %>%
  filter(is.na(classcode) | is.na(class))
write_csv(unmapped, file.path(out_dir, "bayes_unmapped_treatments.csv"))

map_consistency <- trt_map %>%
  group_by(trtcode) %>%
  summarise(n_class = n_distinct(classcode), .groups = "drop") %>%
  filter(n_class != 1)
write_csv(map_consistency, file.path(out_dir, "bayes_mapping_inconsistencies.csv"))

trt_class_counts <- dat_class %>%
  count(treatment, classcode, class, sort = TRUE)
write_csv(trt_class_counts, file.path(out_dir, "bayes_trt_class_counts_used.csv"))

dat_audit <- dat %>%
  left_join(trt_map %>% select(treatment = trtcode, classcode, class), by = "treatment") %>%
  mutate(
    drop_reason = case_when(
      is.na(studyid) ~ "missing_studyid",
      is.na(classcode) | is.na(class) ~ "unmapped_class",
      is.na(n) ~ "missing_n",
      is.na(mean_change) ~ "missing_mean_change",
      is.na(sd_change) ~ "missing_sd_change",
      TRUE ~ "kept"
    )
  )
write_csv(count(dat_audit, drop_reason), file.path(out_dir, "bayes_drop_reason_counts.csv"))

message("Bayesian multinma class-level analysis complete.")
