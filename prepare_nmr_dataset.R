# ============================================================
# prepare_nmr_dataset.R
#
# PURPOSE
# -------
# Merge ms_depression_data.xlsx (tab: "MS SMD base-case") with
# mmc3_included_studies.xlsx (tab: "MS depression-included
# studies") to produce a single long-format arm-level dataset
# suitable for class-level network meta-regression (NMR).
#
# The merged dataset contains the continuous outcome data
# (mean change, SD, n) together with study-level covariates
# from mmc3: country, trial duration, mean age, % female,
# baseline depression scale/score, and risk-of-bias items.
#
# THREE DATA FORMATS in ms_depression_data (MS SMD base-case)
# -----------------------------------------------------------
#  Block 1 – Change-from-baseline (CFB)
#    Columns: na | t[1:5] | yCFB[1:5] | sdCFB[1:5] | nCFB[1:5]
#    studyid column: col 23 (1-indexed R)
#    Excel rows 3–152  (150 studies)
#
#  Block 2 – Baseline + Follow-up (B+F)
#    Columns: na | t[1:5] | yB[1:5] | sdB[1:5] | yF[1:5]
#             | sdF[1:5] | n[1:5]
#    studyid column: col 33 (1-indexed R)
#    Excel rows 154–328  (175 studies)
#
#  Block 3 – Binary response + Baseline (BIN)
#    Columns: na | t[1:5] | r[1:5] | n[1:5] | yBR[1:5]
#             | sdBR[1:5] | q
#    studyid column: col 29 (1-indexed R)
#    Excel rows 331–364  (34 studies)
#
# CONVERSION ASSUMPTIONS (vague priors)
# -------------------------------------
#  Block 2 → CFB:
#    mean_change = yF - yB
#    sd_change   = sqrt(sdB^2 + sdF^2 - 2 * rho * sdB * sdF)
#    rho = 0.5   (standard weak-assumption; sensitivity in 0.3–0.7
#                 can be examined by changing RHO_BF below)
#
#  Block 3 → CFB (Furukawa normal-distribution method):
#    response rate p = r / n
#    response threshold c = yBR * (1 - q)
#      (e.g. q=0.5 → 50 % reduction from baseline = response)
#    Assuming post-treatment scores ~ N(mu_F, sdBR^2):
#      mu_F = c - sdBR * qnorm(p)        [solve P(yF < c) = p]
#    mean_change = mu_F - yBR = -yBR*q - sdBR * qnorm(p)
#    sd_change   = sdBR * sqrt(2 * (1 - rho))
#    rho = 0.5   (same vague prior; change via RHO_BIN below)
#
#  Both conversions tag the row with a data_type flag so you
#  can run sensitivity analyses excluding imputed studies.
#
# OUTPUT
# ------
#  <out_dir>/nmr_dataset_long.csv  – long arm-level data + covariates
#  <out_dir>/nmr_studylevel.csv    – one row per study (covariates only)
#  <out_dir>/nmr_missing_mmc3.csv  – studies without mmc3 match
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# ------------------------------------------------------------------
# 0.  User-configurable paths and assumptions
# ------------------------------------------------------------------
base_dir  <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli/article_supplements"
ms_file   <- file.path(base_dir, "ms_depression_data.xlsx")
mmc3_file <- file.path(base_dir, "mmc3_included_studies.xlsx")
out_dir   <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli/netmeta_class_ms"

# Imputed correlation for sd_change (vague prior; change to explore sensitivity)
RHO_BF  <- 0.5   # Block 2 (Baseline + Follow-up)
RHO_BIN <- 0.5   # Block 3 (Binary response)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------
# 1.  Class mapping (identical to netmeta_class_ms.R)
# ------------------------------------------------------------------
map_classcode <- function(trtcode) {
  dplyr::case_when(
    trtcode == 1  ~  1L, trtcode == 2  ~  2L, trtcode == 3  ~  3L,
    trtcode == 4  ~  4L, trtcode == 5  ~  5L, trtcode == 6  ~  6L,
    trtcode == 7  ~  7L, trtcode %in% c(8, 9) ~ 8L,
    trtcode %in% c(10, 11, 12, 13) ~  9L, trtcode == 14 ~ 10L,
    trtcode == 15 ~ 11L, trtcode == 16 ~ 12L, trtcode == 17 ~ 13L,
    trtcode == 18 ~ 14L, trtcode == 19 ~ 15L,
    trtcode %in% c(20, 21, 22, 23) ~ 16L,
    trtcode %in% c(24, 25, 26, 27) ~ 17L,
    trtcode %in% c(28, 29) ~ 18L, trtcode == 30 ~ 19L,
    trtcode == 31 ~ 20L, trtcode == 32 ~ 21L, trtcode == 33 ~ 22L,
    trtcode %in% c(34, 35) ~ 23L, trtcode == 36 ~ 24L,
    trtcode == 37 ~ 25L, trtcode == 38 ~ 26L,
    trtcode %in% c(39, 40, 41, 42, 43, 44) ~ 27L,
    trtcode %in% c(45, 46, 47, 48, 49, 50) ~ 28L,
    trtcode %in% c(51, 52) ~ 29L, trtcode == 53 ~ 30L,
    trtcode %in% c(54, 55, 56) ~ 31L, trtcode %in% c(57, 58, 59) ~ 32L,
    trtcode %in% c(60, 61, 62) ~ 33L, trtcode %in% c(63, 64) ~ 34L,
    trtcode == 65 ~ 35L, trtcode == 66 ~ 36L,
    trtcode %in% c(67, 68) ~ 37L,
    trtcode %in% c(69, 70, 71, 72, 73, 74, 75) ~ 38L,
    trtcode == 76 ~ 39L, trtcode %in% c(77, 78) ~ 40L,
    trtcode %in% c(79, 80, 81) ~ 41L, trtcode %in% c(82, 83) ~ 42L,
    trtcode == 84 ~ 43L, trtcode == 85 ~ 44L, trtcode == 86 ~ 45L,
    trtcode %in% c(87, 88, 89) ~ 46L, trtcode %in% c(90, 91) ~ 47L,
    trtcode == 92 ~ 48L,
    trtcode %in% c(93, 94, 95, 96, 97) ~ 49L,
    trtcode %in% c(98, 99) ~ 50L,
    TRUE ~ NA_integer_
  )
}

map_classname <- function(trtcode) {
  dplyr::case_when(
    trtcode == 1  ~ "Placebo",
    trtcode == 2  ~ "Attention placebo",
    trtcode == 3  ~ "No treatment",
    trtcode == 4  ~ "Waitlist",
    trtcode == 5  ~ "TAU",
    trtcode == 6  ~ "Mirtazapine",
    trtcode == 7  ~ "Trazodone",
    trtcode %in% c(8, 9) ~ "Behavioural therapies individual",
    trtcode %in% c(10, 11, 12, 13) ~
      "Cognitive and cognitive behavioural therapies individual",
    trtcode == 14 ~ "Cognitive and cognitive behavioural therapies group",
    trtcode == 15 ~ "Problem solving individual",
    trtcode == 16 ~ "Problem solving group",
    trtcode == 17 ~ "Counselling individual",
    trtcode == 18 ~ "Interpersonal psychotherapy (IPT) individual",
    trtcode == 19 ~ "Psychoeducation group",
    trtcode %in% c(20, 21, 22, 23) ~ "Self-help",
    trtcode %in% c(24, 25, 26, 27) ~ "Self-help with support",
    trtcode %in% c(28, 29) ~
      "Short-term psychodynamic psychotherapies individual",
    trtcode == 30 ~ "Music therapy group",
    trtcode == 31 ~ "Mindfulness or meditation group",
    trtcode == 32 ~ "Peer support group",
    trtcode == 33 ~ "Any psychotherapy",
    trtcode %in% c(34, 35) ~
      "Cognitive and cognitive behavioural therapies individual + placebo",
    trtcode == 36 ~
      "Interpersonal psychotherapy (IPT) individual + placebo",
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
    trtcode %in% c(69, 70, 71, 72, 73, 74, 75) ~
      "Cognitive and cognitive behavioural therapies individual + AD",
    trtcode == 76 ~
      "Cognitive and cognitive behavioural therapies group + AD",
    trtcode %in% c(77, 78) ~
      "Interpersonal psychotherapy (IPT) individual + AD",
    trtcode %in% c(79, 80, 81) ~ "Counselling individual + AD",
    trtcode %in% c(82, 83) ~
      "Short-term psychodynamic psychotherapies individual + AD",
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
}

# ------------------------------------------------------------------
# 2.  Helper: pivot a wide (per-study) row to arm-level long rows
# ------------------------------------------------------------------
# wide_df   : one row per study, columns named t1..t5, x1..x5, y1..y5, n1..n5
# id_col    : name of the study-id column
# n_arms_col: column giving number of arms (na[])
# Returns a long tibble with columns: studyid, arm, treatment, <x>, <y>, n
pivot_to_arms <- function(wide_df, id_col, n_arms_col,
                          t_cols,   # e.g. c("t1","t2","t3","t4","t5")
                          val_cols, # named list: outcome_name -> c("v1",..,"v5")
                          n_cols) { # e.g. c("n1","n2","n3","n4","n5")
  MAX_ARMS <- 5L

  wide_df %>%
    dplyr::filter(!is.na(.data[[id_col]]),
                  !is.na(.data[[n_arms_col]])) %>%
    dplyr::mutate(
      na_int = suppressWarnings(as.integer(.data[[n_arms_col]]))
    ) %>%
    dplyr::filter(!is.na(na_int)) %>%
    dplyr::rowwise() %>%
    dplyr::do({
      row   <- .
      na    <- row$na_int
      sid   <- row[[id_col]]
      arms  <- seq_len(min(na, MAX_ARMS))

      arm_list <- lapply(arms, function(k) {
        trt <- suppressWarnings(as.numeric(row[[t_cols[k]]]))
        n   <- suppressWarnings(as.numeric(row[[n_cols[k]]]))

        vals <- lapply(val_cols, function(vcols) {
          suppressWarnings(as.numeric(row[[vcols[k]]]))
        })
        names(vals) <- names(val_cols)

        c(list(studyid = sid, arm = k, treatment = trt, n = n), vals)
      })

      dplyr::bind_rows(lapply(arm_list, as.data.frame, stringsAsFactors = FALSE))
    }) %>%
    dplyr::ungroup()
}

# ------------------------------------------------------------------
# 3.  Read Block 1 – Change-from-baseline (CFB)
#     Excel rows 3–152: skip 2 rows (title + header), read 150 data rows.
#     Use col_names = FALSE to avoid NA/empty names from blank header cells.
# ------------------------------------------------------------------
message("Reading Block 1: CFB format ...")

# Build a full 50-column name vector; unused slots get a generic name
b1_colnames <- c(
  "na_arms",
  paste0("t",     1:5),   # cols 2–6
  paste0("yCFB",  1:5),   # cols 7–11
  paste0("sdCFB", 1:5),   # cols 12–16
  paste0("nCFB",  1:5),   # cols 17–21
  "hash_col",              # col 22  (#)
  "studyid",               # col 23
  paste0("spare_", 24:50) # cols 24–50  (lookup table + blank cells)
)

b1_raw <- read_excel(
  ms_file,
  sheet     = "MS SMD base-case",
  skip      = 2,           # skip title row AND header row
  n_max     = 150,
  col_names = b1_colnames,
  .name_repair = "minimal"
)

b1_raw <- b1_raw %>%
  dplyr::filter(
    !is.na(studyid),
    is.character(studyid) | is.character(as.character(studyid)),
    !grepl("^na\\[", as.character(studyid), ignore.case = TRUE)  # drop any stray header
  ) %>%
  dplyr::mutate(
    na_arms = suppressWarnings(as.integer(na_arms))
  ) %>%
  dplyr::filter(!is.na(na_arms))

b1_long <- pivot_to_arms(
  b1_raw,
  id_col     = "studyid",
  n_arms_col = "na_arms",
  t_cols     = paste0("t", 1:5),
  val_cols   = list(mean_change = paste0("yCFB",  1:5),
                    sd_change   = paste0("sdCFB", 1:5)),
  n_cols     = paste0("nCFB", 1:5)
) %>%
  dplyr::mutate(data_type = "CFB")

message(sprintf("  Block 1: %d studies, %d arms", n_distinct(b1_long$studyid), nrow(b1_long)))

# ------------------------------------------------------------------
# 4.  Read Block 2 – Baseline + Follow-up (B+F)
#     Excel rows 154–328: skip 153 rows, read 175 data rows.
# ------------------------------------------------------------------
message("Reading Block 2: Baseline+Follow-up format ...")

# Layout: na[] | t[1:5] | yB[1:5] | sdB[1:5] | yF[1:5] | sdF[1:5] | n[1:5] | # | studyid
b2_colnames <- c(
  "na_arms",
  paste0("t",   1:5),   # cols 2–6
  paste0("yB",  1:5),   # cols 7–11
  paste0("sdB", 1:5),   # cols 12–16
  paste0("yF",  1:5),   # cols 17–21
  paste0("sdF", 1:5),   # cols 22–26
  paste0("n",   1:5),   # cols 27–31
  "hash_col",           # col 32
  "studyid",            # col 33
  paste0("spare_", 34:50)   # cols 34–50  (50 total to match sheet width)
)

b2_raw <- read_excel(
  ms_file,
  sheet     = "MS SMD base-case",
  skip      = 153,        # title + B1 header + B1 data (150) + B2 header = 153
  n_max     = 175,
  col_names = b2_colnames,
  .name_repair = "minimal"
)

b2_raw <- b2_raw %>%
  dplyr::filter(
    !is.na(studyid),
    !grepl("^na\\[|^studyid", as.character(studyid), ignore.case = TRUE)
  ) %>%
  dplyr::mutate(na_arms = suppressWarnings(as.integer(na_arms))) %>%
  dplyr::filter(!is.na(na_arms))

# For each arm compute: mean_change = yF - yB
#                       sd_change   = sqrt(sdB^2 + sdF^2 - 2*rho*sdB*sdF)
b2_long <- pivot_to_arms(
  b2_raw,
  id_col     = "studyid",
  n_arms_col = "na_arms",
  t_cols     = paste0("t", 1:5),
  val_cols   = list(yB  = paste0("yB",  1:5),
                    sdB = paste0("sdB", 1:5),
                    yF  = paste0("yF",  1:5),
                    sdF = paste0("sdF", 1:5)),
  n_cols     = paste0("n", 1:5)
) %>%
  dplyr::mutate(
    mean_change = yF - yB,
    # sd_change with vague-prior correlation rho = RHO_BF
    sd_change   = sqrt(sdB^2 + sdF^2 - 2 * RHO_BF * sdB * sdF),
    data_type   = "BF"
  ) %>%
  dplyr::select(-yB, -sdB, -yF, -sdF)

message(sprintf("  Block 2: %d studies, %d arms", n_distinct(b2_long$studyid), nrow(b2_long)))

# ------------------------------------------------------------------
# 5.  Read Block 3 – Binary response + Baseline (BIN)
#     Excel rows 331–364: skip 330 rows, read 34 data rows.
# ------------------------------------------------------------------
message("Reading Block 3: Binary response format ...")

# Layout: na[] | t[1:5] | r[1:5] | n[1:5] | yBR[1:5] | sdBR[1:5] | q | # | studyid
b3_colnames <- c(
  "na_arms",
  paste0("t",    1:5),   # cols 2–6
  paste0("r",    1:5),   # cols 7–11  (responders)
  paste0("n",    1:5),   # cols 12–16 (total)
  paste0("yBR",  1:5),   # cols 17–21 (baseline mean)
  paste0("sdBR", 1:5),   # cols 22–26 (baseline SD)
  "q_thresh",            # col 27
  "hash_col",            # col 28
  "studyid"              # col 29  (last column — no spares needed)
)

b3_raw <- read_excel(
  ms_file,
  sheet     = "MS SMD base-case",
  skip      = 330,        # title + B1 h + B1 data(150) + B2 h + B2 data(175) + B3 h = 330
  n_max     = 34,
  col_names = b3_colnames,
  .name_repair = "minimal"
)

b3_raw <- b3_raw %>%
  dplyr::filter(
    !is.na(studyid),
    !grepl("^na\\[|^studyid", as.character(studyid), ignore.case = TRUE)
  ) %>%
  dplyr::mutate(na_arms = suppressWarnings(as.integer(na_arms))) %>%
  dplyr::filter(!is.na(na_arms))

# Furukawa normal-distribution conversion:
#   response rate p = r / n
#   threshold on follow-up scale c = yBR * (1 - q)  [e.g. 50% reduction]
#   Assuming yF ~ N(mu_F, sdBR^2):
#     mu_F = c - sdBR * qnorm(p)
#   mean_change = mu_F - yBR = -yBR*q - sdBR * qnorm(p)
#   sd_change   = sdBR * sqrt(2 * (1 - rho))   [rho = RHO_BIN]
b3_long <- pivot_to_arms(
  b3_raw,
  id_col     = "studyid",
  n_arms_col = "na_arms",
  t_cols     = paste0("t", 1:5),
  val_cols   = list(
    r    = paste0("r",    1:5),
    yBR  = paste0("yBR",  1:5),
    sdBR = paste0("sdBR", 1:5)
  ),
  n_cols = paste0("n", 1:5)
) %>%
  # q_thresh is study-level; join back from b3_raw
  dplyr::left_join(
    b3_raw %>%
      dplyr::select(studyid,
                    q_thresh = q_thresh) %>%
      dplyr::mutate(q_thresh = suppressWarnings(as.numeric(q_thresh))),
    by = "studyid"
  ) %>%
  dplyr::mutate(
    p           = r / n,
    # Guard against p=0 or p=1 (qnorm undefined at boundaries)
    p_clipped   = pmin(pmax(p, 0.001), 0.999),
    mean_change = -yBR * q_thresh - sdBR * stats::qnorm(p_clipped),
    sd_change   = sdBR * sqrt(2 * (1 - RHO_BIN)),
    data_type   = "BIN"
  ) %>%
  dplyr::select(-r, -yBR, -sdBR, -p, -p_clipped, -q_thresh)

message(sprintf("  Block 3: %d studies, %d arms", n_distinct(b3_long$studyid), nrow(b3_long)))

# ------------------------------------------------------------------
# 6.  Combine all blocks
# ------------------------------------------------------------------
message("Combining blocks ...")

nma_long <- dplyr::bind_rows(b1_long, b2_long, b3_long) %>%
  dplyr::mutate(
    treatment   = suppressWarnings(as.integer(treatment)),
    n           = suppressWarnings(as.numeric(n)),
    mean_change = suppressWarnings(as.numeric(mean_change)),
    sd_change   = suppressWarnings(as.numeric(sd_change))
  ) %>%
  # Remove arms with missing core data
  dplyr::filter(
    !is.na(treatment),
    !is.na(n),
    !is.na(mean_change),
    !is.na(sd_change),
    n > 0,
    sd_change > 0
  ) %>%
  # Apply class mapping
  dplyr::mutate(
    classcode = map_classcode(treatment),
    class     = map_classname(treatment)
  )

message(sprintf("  Combined: %d studies, %d arms", n_distinct(nma_long$studyid), nrow(nma_long)))

# ------------------------------------------------------------------
# 7.  Read mmc3 study-level covariates
# ------------------------------------------------------------------
message("Reading mmc3 covariates ...")

mmc3_raw <- read_excel(
  mmc3_file,
  sheet     = "MS depression-included studies",
  col_names = TRUE,
  .name_repair = "minimal"
)

# Expected columns (by name from header row):
#  Study ID | Arm1..Arm5 | Incl/Excl | N randomised | Country |
#  Duration (weeks) | Mean age (SD) | Sex (% female) |
#  Baseline depression scale | Mean baseline depression score |
#  SD of baseline depression score |
#  [ROB items x7] | Full reference

mmc3 <- mmc3_raw %>%
  dplyr::rename_with(~ dplyr::case_when(
    .x == "Study ID"                                       ~ "studyid",
    .x == "N randomised"                                   ~ "n_randomised",
    .x == "Country"                                        ~ "country",
    .x == "Duration (weeks)"                               ~ "duration_weeks",
    .x == "Mean age (SD)"                                  ~ "mean_age_raw",
    .x == "Sex (% female)"                                 ~ "pct_female",
    .x == "Baseline depression scale"                      ~ "baseline_scale",
    .x == "Mean baseline depression score"                 ~ "baseline_mean",
    .x == "SD of baseline depression score"                ~ "baseline_sd",
    .x == "Random sequence generation (selection bias)"    ~ "rob_sequence",
    .x == "Allocation concealment (selection bias)"        ~ "rob_allocation",
    .x == "Blinding of participants and personnel (performance bias)" ~
      "rob_blind_pp",
    .x == "Blinding of outcome assessment (detection bias)" ~
      "rob_blind_outcome",
    .x == "Incomplete outcome data (attrition bias)"       ~ "rob_attrition",
    .x == "Selective reporting (reporting bias)"           ~ "rob_reporting",
    .x == "Other bias"                                     ~ "rob_other",
    .x == "Full reference"                                 ~ "full_reference",
    TRUE ~ .x
  )) %>%
  dplyr::filter(!is.na(studyid)) %>%
  dplyr::mutate(
    # Extract numeric mean age from strings like "41.0 (11.2)" or "39 (SD NR)"
    mean_age = suppressWarnings(
      as.numeric(stringr::str_extract(as.character(mean_age_raw), "^[0-9.]+"))
    ),
    pct_female     = suppressWarnings(as.numeric(pct_female)),
    duration_weeks = suppressWarnings(as.numeric(duration_weeks)),
    n_randomised   = suppressWarnings(as.integer(n_randomised)),
    baseline_mean  = suppressWarnings(as.numeric(baseline_mean)),
    baseline_sd    = suppressWarnings(as.numeric(baseline_sd))
  ) %>%
  dplyr::select(
    studyid, country, duration_weeks, mean_age, pct_female,
    baseline_scale, baseline_mean, baseline_sd,
    rob_sequence, rob_allocation, rob_blind_pp,
    rob_blind_outcome, rob_attrition, rob_reporting, rob_other,
    full_reference
  )

message(sprintf("  mmc3: %d studies with covariates", nrow(mmc3)))

# ------------------------------------------------------------------
# 8.  Merge outcome data with covariates
# ------------------------------------------------------------------
message("Merging ...")

nmr_data <- nma_long %>%
  dplyr::left_join(mmc3, by = "studyid")

# Studies without mmc3 match
missing_mmc3 <- nmr_data %>%
  dplyr::filter(is.na(country)) %>%
  dplyr::distinct(studyid, data_type)

if (nrow(missing_mmc3) > 0) {
  message(sprintf(
    "  WARNING: %d studies could not be matched to mmc3 covariates:",
    n_distinct(missing_mmc3$studyid)
  ))
  message(paste("   ", missing_mmc3$studyid, collapse = "\n"))
}

# ------------------------------------------------------------------
# 9.  Derive convenience / NMR-ready variables
# ------------------------------------------------------------------
nmr_data <- nmr_data %>%
  dplyr::mutate(
    # Proportion male (complement of pct_female expressed as %)
    pct_male = 100 - pct_female,

    # Numeric ROB score: Low=0, Unclear=1, High=2 (for regression sensitivity)
    rob_score_sequence    = dplyr::case_when(
      grepl("Low",    rob_sequence,    ignore.case = TRUE) ~ 0L,
      grepl("Unclear",rob_sequence,    ignore.case = TRUE) ~ 1L,
      grepl("High",   rob_sequence,    ignore.case = TRUE) ~ 2L,
      TRUE ~ NA_integer_
    ),
    rob_score_allocation  = dplyr::case_when(
      grepl("Low",    rob_allocation,  ignore.case = TRUE) ~ 0L,
      grepl("Unclear",rob_allocation,  ignore.case = TRUE) ~ 1L,
      grepl("High",   rob_allocation,  ignore.case = TRUE) ~ 2L,
      TRUE ~ NA_integer_
    ),
    rob_score_blind_pp    = dplyr::case_when(
      grepl("Low",    rob_blind_pp,    ignore.case = TRUE) ~ 0L,
      grepl("Unclear",rob_blind_pp,    ignore.case = TRUE) ~ 1L,
      grepl("High",   rob_blind_pp,    ignore.case = TRUE) ~ 2L,
      TRUE ~ NA_integer_
    ),
    rob_score_blind_out   = dplyr::case_when(
      grepl("Low",    rob_blind_outcome,ignore.case = TRUE) ~ 0L,
      grepl("Unclear",rob_blind_outcome,ignore.case = TRUE) ~ 1L,
      grepl("High",   rob_blind_outcome,ignore.case = TRUE) ~ 2L,
      TRUE ~ NA_integer_
    ),
    rob_score_attrition   = dplyr::case_when(
      grepl("Low",    rob_attrition,   ignore.case = TRUE) ~ 0L,
      grepl("Unclear",rob_attrition,   ignore.case = TRUE) ~ 1L,
      grepl("High",   rob_attrition,   ignore.case = TRUE) ~ 2L,
      TRUE ~ NA_integer_
    ),

    # Overall ROB flag: "High" if any domain is high risk
    rob_high_any = dplyr::if_else(
      rob_score_sequence  == 2L |
        rob_score_allocation == 2L |
        rob_score_blind_pp   == 2L |
        rob_score_blind_out  == 2L |
        rob_score_attrition  == 2L,
      1L, 0L
    )
  ) %>%
  dplyr::arrange(studyid, arm)

# ------------------------------------------------------------------
# 10.  Study-level summary (one row per study)
# ------------------------------------------------------------------
nmr_studylevel <- nmr_data %>%
  dplyr::distinct(
    studyid, data_type,
    country, duration_weeks, mean_age, pct_female,
    baseline_scale, baseline_mean, baseline_sd,
    rob_sequence, rob_allocation, rob_blind_pp,
    rob_blind_outcome, rob_attrition, rob_reporting, rob_other,
    rob_score_sequence, rob_score_allocation,
    rob_score_blind_pp, rob_score_blind_out, rob_score_attrition,
    rob_high_any,
    full_reference
  )

# ------------------------------------------------------------------
# 11.  Save outputs
# ------------------------------------------------------------------
message("Saving outputs ...")

# Final column order for the long dataset
nmr_out <- nmr_data %>%
  dplyr::select(
    studyid, arm, treatment, classcode, class,
    n, mean_change, sd_change, data_type,
    country, duration_weeks, mean_age, pct_female,
    baseline_scale, baseline_mean, baseline_sd,
    rob_sequence, rob_allocation, rob_blind_pp,
    rob_blind_outcome, rob_attrition, rob_reporting, rob_other,
    rob_score_sequence, rob_score_allocation,
    rob_score_blind_pp, rob_score_blind_out, rob_score_attrition,
    rob_high_any
  )

readr::write_csv(nmr_out,        file.path(out_dir, "nmr_dataset_long.csv"))
readr::write_csv(nmr_studylevel, file.path(out_dir, "nmr_studylevel.csv"))
readr::write_csv(missing_mmc3,   file.path(out_dir, "nmr_missing_mmc3.csv"))

message(sprintf(
  "\nDone. Output: %d studies, %d arms written to %s",
  n_distinct(nmr_out$studyid), nrow(nmr_out), out_dir
))
message("  nmr_dataset_long.csv   – arm-level data + covariates (use for NMR)")
message("  nmr_studylevel.csv     – study-level covariates only")
message("  nmr_missing_mmc3.csv   – studies without covariate match")

# ------------------------------------------------------------------
# 12.  Quick diagnostics
# ------------------------------------------------------------------
message("\n--- Coverage summary ---")
message(sprintf("Studies with mean_age    : %d / %d",
                sum(!is.na(nmr_studylevel$mean_age)),     nrow(nmr_studylevel)))
message(sprintf("Studies with pct_female  : %d / %d",
                sum(!is.na(nmr_studylevel$pct_female)),   nrow(nmr_studylevel)))
message(sprintf("Studies with country     : %d / %d",
                sum(!is.na(nmr_studylevel$country)),      nrow(nmr_studylevel)))
message(sprintf("Studies with duration    : %d / %d",
                sum(!is.na(nmr_studylevel$duration_weeks)),nrow(nmr_studylevel)))
message(sprintf("Arms by data_type (CFB / BF / BIN): %s / %s / %s",
                sum(nmr_out$data_type == "CFB"),
                sum(nmr_out$data_type == "BF"),
                sum(nmr_out$data_type == "BIN")))
message(sprintf("Unmapped treatment codes: %d",
                sum(is.na(nmr_out$classcode))))

# ------------------------------------------------------------------
# 13.  Merge sanity checks
# ------------------------------------------------------------------
message("\n--- Merge sanity checks ---")

# (a) No duplicate study × arm combinations
dup_arms <- nmr_out %>%
  dplyr::count(studyid, arm) %>%
  dplyr::filter(n > 1)
if (nrow(dup_arms) == 0) {
  message("  [OK] No duplicate study × arm combinations")
} else {
  message(sprintf("  [WARN] %d duplicate study × arm combinations:", nrow(dup_arms)))
  print(dup_arms)
}

# (b) Plausible covariate ranges
age_range <- range(nmr_studylevel$mean_age, na.rm = TRUE)
sex_range <- range(nmr_studylevel$pct_female, na.rm = TRUE)
dur_range <- range(nmr_studylevel$duration_weeks, na.rm = TRUE)
message(sprintf("  mean_age range    : %.1f – %.1f  (expect ~20–80)", age_range[1], age_range[2]))
message(sprintf("  pct_female range  : %.1f – %.1f  (expect 0–100)",  sex_range[1], sex_range[2]))
message(sprintf("  duration_weeks    : %.1f – %.1f  (expect 1–104)",   dur_range[1], dur_range[2]))
if (age_range[1] < 10 | age_range[2] > 100)
  message("  [WARN] mean_age outside expected range – check mmc3 parsing")
if (sex_range[1] < 0 | sex_range[2] > 100)
  message("  [WARN] pct_female outside 0–100 – check mmc3 parsing")

# (c) Spot-check specific well-known studies (present in both files)
spot_ids <- c("Blumenthal 2007/Hoffman 2011", "Elkin 1989/Imber 1990", "Kasper 2005a")
message("\n  Spot-check studies (studyid | country | duration | age | pct_female):")
spot <- nmr_studylevel %>%
  dplyr::filter(studyid %in% spot_ids) %>%
  dplyr::select(studyid, country, duration_weeks, mean_age, pct_female)
if (nrow(spot) > 0) {
  for (i in seq_len(nrow(spot))) {
    message(sprintf("    %-40s | %-8s | %5.1f wk | age %4.1f | sex %4.1f%%",
                    spot$studyid[i],
                    ifelse(is.na(spot$country[i]), "NA", spot$country[i]),
                    ifelse(is.na(spot$duration_weeks[i]), NA_real_, spot$duration_weeks[i]),
                    ifelse(is.na(spot$mean_age[i]),     NA_real_, spot$mean_age[i]),
                    ifelse(is.na(spot$pct_female[i]),   NA_real_, spot$pct_female[i])))
  }
} else {
  message("    (none of the spot-check studies found in merged data)")
}

# (d) Class and treatment code consistency
class_check <- nmr_out %>%
  dplyr::distinct(treatment, classcode, class) %>%
  dplyr::group_by(treatment) %>%
  dplyr::filter(n_distinct(classcode) > 1 | n_distinct(class) > 1) %>%
  dplyr::ungroup()
if (nrow(class_check) == 0) {
  message("\n  [OK] All treatment codes map to a single class")
} else {
  message(sprintf("\n  [WARN] %d treatment codes with inconsistent class mapping:", nrow(class_check)))
  print(class_check)
}

# (e) SMD plausibility: flag extreme values
extreme_smd <- nmr_out %>%
  # Approximate SMD as mean_change / sd_change (within-arm; used only as sanity signal)
  dplyr::mutate(approx_smd = abs(mean_change / sd_change)) %>%
  dplyr::filter(approx_smd > 5)
if (nrow(extreme_smd) == 0) {
  message("  [OK] No arms with |mean_change / sd_change| > 5")
} else {
  message(sprintf("  [WARN] %d arms with |mean_change/sd_change| > 5 (possible scaling issue):",
                  nrow(extreme_smd)))
  print(extreme_smd %>% dplyr::select(studyid, arm, treatment, data_type,
                                       mean_change, sd_change, approx_smd))
}

message("\nSanity checks complete.")

