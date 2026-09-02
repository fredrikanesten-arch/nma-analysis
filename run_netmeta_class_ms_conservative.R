library(dplyr)
library(readr)
library(netmeta)

# -----------------------------
# 0. Paths
# -----------------------------
base_dir <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli"
in_data  <- file.path(base_dir, "binfixed_class_ms", "combined_long_mean_change_dataset_ms_smd_bias_adj.csv")
in_map   <- file.path(base_dir, "clean_data", "trt_to_class_ms.csv")
out_dir  <- file.path(base_dir, "binfixed_class_ms")

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

read_map_robust <- function(path) {
  dat_scsv <- tryCatch(
    readr::read_delim(path, delim = ";", show_col_types = FALSE),
    error = function(e) NULL
  )
  if (!is.null(dat_scsv) && ncol(dat_scsv) > 1) {
    message("Read mapping as semicolon-separated CSV.")
    return(dat_scsv)
  }

  dat_csv <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (!is.null(dat_csv) && ncol(dat_csv) > 1) {
    message("Read mapping as comma-separated CSV.")
    return(dat_csv)
  }

  stop("Could not parse treatment mapping file: ", path)
}

# -----------------------------
# 2. Read data and convert types
# -----------------------------
dat <- read_combined_long_robust(in_data)
trt_map <- read_map_robust(in_map)

required_map_cols <- c("trtcode", "classcode", "class")
missing_map_cols <- setdiff(required_map_cols, names(trt_map))
if (length(missing_map_cols) > 0) {
  stop("Mapping file missing required columns: ", paste(missing_map_cols, collapse = ", "))
}

dat <- dat %>%
  mutate(
    across(any_of(c("na", "arm", "treatment", "n", "mean_change", "sd_change")), as.numeric)
  )

trt_map <- trt_map %>%
  mutate(
    trtcode = as.numeric(trtcode),
    classcode = as.numeric(classcode),
    class = as.character(class)
  )

# Keep one deterministic mapping row per treatment (first non-missing class tuple).
trt_map_dedup <- trt_map %>%
  filter(!is.na(trtcode), !is.na(classcode), !is.na(class), nzchar(class)) %>%
  arrange(trtcode, classcode, class) %>%
  group_by(trtcode) %>%
  slice(1) %>%
  ungroup()

# -----------------------------
# 3. Join + keep analyzable rows
# -----------------------------
dat_class <- dat %>%
  left_join(
    trt_map_dedup %>% select(treatment = trtcode, classcode, class),
    by = "treatment"
  ) %>%
  filter(
    !is.na(studyid),
    !is.na(classcode),
    !is.na(class),
    !is.na(n),
    !is.na(mean_change),
    !is.na(sd_change),
    n > 1,
    sd_change > 0
  )

# -----------------------------
# 4. Collapse duplicate same-class arms within study (conservative)
# -----------------------------
# Uses full combined variance including between-arm mean dispersion:
# var = (sum((n-1)*sd^2) + sum(n*(mean-mean_total)^2)) / (sum(n)-1)
collapsed_class <- dat_class %>%
  group_by(studyid, classcode, class) %>%
  summarise(
    n = sum(n),
    mean_change = sum(mean_change * n) / sum(n),
    sd_change = {
      n_total <- sum(n)
      mean_total <- sum(mean_change * n) / n_total
      ss_within <- sum((n - 1) * sd_change^2)
      ss_between <- sum(n * (mean_change - mean_total)^2)
      sqrt((ss_within + ss_between) / (n_total - 1))
    },
    .groups = "drop"
  ) %>%
  group_by(studyid) %>%
  filter(n_distinct(class) >= 2) %>%
  ungroup()

stopifnot("Placebo" %in% collapsed_class$class)

# -----------------------------
# 5. Pairwise contrasts + NMA (conservative scale handling)
# -----------------------------
# Input arm means are already on a standardized scale from preprocessing,
# so we keep sm = "MD" to avoid re-standardizing in pairwise().
pw <- pairwise(
  treat = class,
  mean = mean_change,
  sd = sd_change,
  n = n,
  studlab = studyid,
  data = collapsed_class,
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
  reference.group = "Placebo",
  verbose = TRUE
)

# -----------------------------
# 6. Summary table vs Placebo
# -----------------------------
make_nma_excel_table_full <- function(nma, collapsed_class, ref = "Placebo",
                                      small.values = "good", digits = 2) {
  if (!ref %in% nma$trts) stop("Reference treatment '", ref, "' not found in nma$trts")

  n_tbl <- collapsed_class %>%
    group_by(class) %>%
    summarise(N = sum(n, na.rm = TRUE), .groups = "drop") %>%
    rename(Treatment = class) %>%
    mutate(`Treatment class` = Treatment)

  TE_mat <- if (!is.null(nma$TE.random)) nma$TE.random else nma$TE.common
  lower_mat <- if (!is.null(nma$lower.random)) nma$lower.random else nma$lower.common
  upper_mat <- if (!is.null(nma$upper.random)) nma$upper.random else nma$upper.common
  se_mat <- if (!is.null(nma$seTE.random)) nma$seTE.random else nma$seTE.common

  extract_against_ref <- function(mat, ref) {
    if (ref %in% colnames(mat)) {
      data.frame(Treatment = rownames(mat), value = as.numeric(mat[, ref]), row.names = NULL)
    } else if (ref %in% rownames(mat)) {
      data.frame(Treatment = colnames(mat), value = as.numeric(mat[ref, ]), row.names = NULL)
    } else stop("Reference treatment not found in matrix dimnames.")
  }

  eff <- extract_against_ref(TE_mat, ref)
  low <- extract_against_ref(lower_mat, ref)
  high <- extract_against_ref(upper_mat, ref)
  se <- extract_against_ref(se_mat, ref)

  effects_tbl <- eff %>%
    rename(SMD = value) %>%
    left_join(rename(low, lower_CI = value), by = "Treatment") %>%
    left_join(rename(high, upper_CI = value), by = "Treatment") %>%
    left_join(rename(se, seTE = value), by = "Treatment") %>%
    mutate(
      SMD = ifelse(Treatment == ref, 0, SMD),
      lower_CI = ifelse(Treatment == ref, 0, lower_CI),
      upper_CI = ifelse(Treatment == ref, 0, upper_CI),
      seTE = ifelse(Treatment == ref, NA, seTE),
      z = ifelse(!is.na(seTE) & seTE > 0, SMD / seTE, NA_real_),
      p_value = ifelse(!is.na(z), 2 * stats::pnorm(-abs(z)), NA_real_)
    )

  rnk <- netmeta::netrank(nma, small.values = small.values)
  ps <- if (!is.null(rnk$ranking.random)) rnk$ranking.random else rnk$ranking.common

  rank_tbl <- data.frame(
    Treatment = names(ps),
    P_score = as.numeric(ps),
    Rank = rank(-as.numeric(ps), ties.method = "average"),
    row.names = NULL
  )

  out <- data.frame(Treatment = nma$trts) %>%
    left_join(n_tbl, by = "Treatment") %>%
    left_join(effects_tbl, by = "Treatment") %>%
    left_join(rank_tbl, by = "Treatment") %>%
    mutate(`Treatment class` = ifelse(is.na(`Treatment class`), Treatment, `Treatment class`))

  fmt <- function(x) formatC(x, digits = digits, format = "f")
  out[[paste0("SMD vs ", ref, " (mean, 95% CI)")]] <- ifelse(
    is.na(out$SMD), NA_character_,
    paste0(fmt(out$SMD), " (", fmt(out$lower_CI), " to ", fmt(out$upper_CI), ")")
  )

  out %>%
    arrange(Rank) %>%
    select(
      Treatment, `Treatment class`, N, SMD, lower_CI, upper_CI,
      all_of(paste0("SMD vs ", ref, " (mean, 95% CI)")),
      z, p_value, P_score, Rank
    ) %>%
    rename(`lower CI` = lower_CI, `upper CI` = upper_CI, `p-value` = p_value, `P-score` = P_score)
}

summary_table_full <- make_nma_excel_table_full(
  nma = nma,
  collapsed_class = collapsed_class,
  ref = "Placebo"
)

# -----------------------------
# 7. Save key outputs
# -----------------------------
write_csv(summary_table_full, file.path(out_dir, "nma_summary_table_full_v2_conservative.csv"))

saveRDS(
  list(dat = dat, trt_map = trt_map_dedup, collapsed_class = collapsed_class, pw = pw, nma = nma),
  file = file.path(out_dir, "nma_class_level_analysis_objects_conservative.rds")
)

# -----------------------------
# 8. Class mapping consistency checks
# -----------------------------
dat_join_check <- dat %>%
  left_join(trt_map_dedup %>% select(treatment = trtcode, classcode, class), by = "treatment")
write_csv(dat_join_check, file.path(out_dir, "dat_join_check_conservative.csv"))

unmapped <- dat %>%
  distinct(treatment) %>%
  left_join(trt_map_dedup %>% select(treatment = trtcode, classcode, class), by = "treatment") %>%
  filter(is.na(classcode) | is.na(class))
write_csv(unmapped, file.path(out_dir, "unmapped_treatments_conservative.csv"))

map_consistency <- trt_map %>%
  filter(!is.na(trtcode), !is.na(classcode)) %>%
  group_by(trtcode) %>%
  summarise(n_class = n_distinct(classcode), .groups = "drop") %>%
  filter(n_class != 1)
write_csv(map_consistency, file.path(out_dir, "mapping_inconsistencies_conservative.csv"))

trt_class_counts <- dat_class %>%
  count(treatment, classcode, class, sort = TRUE)
write_csv(trt_class_counts, file.path(out_dir, "trt_class_counts_used_conservative.csv"))

dat_audit <- dat %>%
  left_join(trt_map_dedup %>% select(treatment = trtcode, classcode, class), by = "treatment") %>%
  mutate(
    drop_reason = case_when(
      is.na(studyid) ~ "missing_studyid",
      is.na(classcode) | is.na(class) ~ "unmapped_class",
      is.na(n) ~ "missing_n",
      is.na(mean_change) ~ "missing_mean_change",
      is.na(sd_change) ~ "missing_sd_change",
      n <= 1 ~ "invalid_n",
      sd_change <= 0 ~ "non_positive_sd_change",
      TRUE ~ "kept"
    )
  )
write_csv(count(dat_audit, drop_reason), file.path(out_dir, "drop_reason_counts_conservative.csv"))
