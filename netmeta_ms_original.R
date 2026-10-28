library(readxl)
library(readr)
library(dplyr)
library(tibble)
library(netmeta)

# -----------------------------
# 1) File paths and settings
# -----------------------------
base_dir <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli"
out_dir <- file.path(base_dir, "netmeta_ms_analyses")

input_xlsx <- file.path(out_dir, "mmc5_frequentist.xlsx")
input_mapping <- file.path(base_dir, "clean_data", "trt_to_class_ms.csv")

assumed_correlation <- 0.31
responder_sd <- 7.5
reference_code <- "1"  # Placebo

clip_response_proportions <- FALSE
response_p_clip <- 1e-6

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2) Read the treatment-class mapping
# -----------------------------
trt_map <- read_delim(
  input_mapping,
  delim = ";",
  locale = locale(encoding = "UTF-8"),
  show_col_types = FALSE
) %>%
  transmute(
    treatment = as.numeric(trtcode),
    classcode = as.character(as.numeric(classcode)),
    class = trimws(class)
  )

if (
  anyNA(trt_map) ||
  any(trt_map$class == "") ||
  anyDuplicated(trt_map$treatment) > 0
) {
  stop("The mapping contains missing values or duplicate treatment codes.")
}

class_labels <- trt_map %>%
  distinct(classcode, class)

if (anyDuplicated(class_labels$classcode) > 0) {
  stop("Each class code must have one class label.")
}

# -----------------------------
# 3) Read the workbook and locate the blocks
# -----------------------------
raw <- as.matrix(read_excel(
  input_xlsx,
  sheet = "MS SMD bias-adj",
  col_names = FALSE,
  col_types = "text",
  .name_repair = "minimal"
))

cfb_header <- which(raw[, 1] == "na[]" & raw[, 7] == "yCFB[,1]")[1]
bf_header <- which(raw[, 1] == "na[]" & raw[, 7] == "yB[,1]")[1]
resp_header <- which(raw[, 1] == "na[]" & raw[, 7] == "r[,1]")[1]

if (
  anyNA(c(cfb_header, bf_header, resp_header)) ||
  !(cfb_header < bf_header && bf_header < resp_header)
) {
  stop("Expected CFB, BF, and RESP headers in that order.")
}

block_starts <- c(cfb_header, bf_header, resp_header) + 1L
block_ends <- c(bf_header - 1L, resp_header - 1L, nrow(raw))

study_id_columns <- c(26, 37, 33)
sources <- c("change_from_baseline", "baseline_followup", "responders")

# -----------------------------
# 4) Extract arm-level data
# -----------------------------
output_rows <- list()

for (block in 1:3) {
  if (block_starts[block] > block_ends[block]) next
  
  for (row in seq.int(block_starts[block], block_ends[block])) {
    values <- suppressWarnings(as.numeric(raw[row, ]))
    study_id <- trimws(raw[row, study_id_columns[block]])
    
    if (
      is.na(values[1]) ||
      is.na(study_id) ||
      study_id %in% c("", "#")
    ) next
    
    for (arm in 1:5) {
      treatment <- values[1 + arm]
      if (is.na(treatment)) next
      
      # Zero placeholders allow the three sources to share one structure.
      y <- b <- s <- i <- r <- q <- 0
      n <- NA_real_
      
      if (block == 1) {
        # Reported change mean, change SD, and sample size.
        y <- values[6 + arm]
        s <- values[11 + arm]
        n <- values[16 + arm]
      }
      
      if (block == 2) {
        # Baseline and follow-up means and SDs, and sample size.
        b <- values[6 + arm]
        i <- values[11 + arm]
        y <- values[16 + arm]
        s <- values[21 + arm]
        n <- values[26 + arm]
      }
      
      if (block == 3) {
        # Responders, sample size, baseline mean and SD, and threshold.
        r <- values[6 + arm]
        n <- values[11 + arm]
        b <- values[16 + arm]
        s <- values[21 + arm]
        q <- values[27]
      }
      
      # Retain the original rule of skipping incomplete arms.
      if (anyNA(c(y, b, s, i, r, q, n))) next
      
      if (
        any(!is.finite(c(y, b, s, i, r, q, n))) ||
        n <= 1 ||
        s <= 0 ||
        (block == 2 && i <= 0) ||
        (block == 3 && (r < 0 || r > n || q <= 0 || q >= 1))
      ) {
        stop("Invalid data at spreadsheet row ", row, ", arm ", arm)
      }
      
      output_rows[[length(output_rows) + 1L]] <- data.frame(
        study_id = study_id,
        source = sources[block],
        arm = arm,
        treatment = treatment,
        n = n,
        y = y,
        b = b,
        s = s,
        i = i,
        r = r,
        q = q
      )
    }
  }
}

if (length(output_rows) == 0) stop("No usable study arms found.")

# Build a fresh dataset from the extracted rows on every run.
arm_data <- bind_rows(output_rows)

unmapped_treatments <- setdiff(arm_data$treatment, trt_map$treatment)

if (length(unmapped_treatments) > 0) {
  stop(
    "Treatment codes missing from the mapping: ",
    paste(unmapped_treatments, collapse = ", ")
  )
}

# Prevent the same study from entering through multiple outcome sources.
multiple_sources <- arm_data %>%
  distinct(study_id, source) %>%
  count(study_id) %>%
  filter(n > 1)

if (nrow(multiple_sources) > 0) {
  print(multiple_sources)
  stop("Resolve studies appearing in multiple source blocks before analysis.")
}

# -----------------------------
# 5) Pool within study and treatment
# -----------------------------
# Preserve the original pooling rules:
# unweighted means, within-arm pooled SDs, and averaged responder counts.
# Sample size is calculated last so earlier expressions use the arm sizes.
pooled_treatment <- arm_data %>%
  group_by(source, study_id, treatment) %>%
  summarise(
    y = mean(y),
    b = mean(b),
    s = sqrt(sum((n - 1) * s^2) / (sum(n) - n())),
    i = sqrt(sum((n - 1) * i^2) / (sum(n) - n())),
    r = mean(r),
    q = first(q),
    n = sum(n),
    .groups = "drop"
  ) %>%
  left_join(
    select(trt_map, treatment, classcode),
    by = "treatment"
  )

# -----------------------------
# 6) Pool within study and class
# -----------------------------
pooled_class <- pooled_treatment %>%
  group_by(source, study_id, classcode) %>%
  summarise(
    y = mean(y),
    b = mean(b),
    s = sqrt(sum((n - 1) * s^2) / (sum(n) - n())),
    i = sqrt(sum((n - 1) * i^2) / (sum(n) - n())),
    r = mean(r),
    q = first(q),
    n = sum(n),
    .groups = "drop"
  )

# -----------------------------
# 7) Calculate mean changes and analysis SDs
# -----------------------------
response_rows <- which(pooled_class$source == "responders")
pooled_class$response_mean <- NA_real_

if (length(response_rows) > 0) {
  response_proportion <- pooled_class$r[response_rows] /
    pooled_class$n[response_rows]
  
  if (clip_response_proportions) {
    response_proportion <- pmin(
      pmax(response_proportion, response_p_clip),
      1 - response_p_clip
    )
  }
  
  pooled_class$response_mean[response_rows] <-
    -pooled_class$q[response_rows] * pooled_class$b[response_rows] -
    responder_sd * qnorm(response_proportion)
}

analysis_data <- pooled_class %>%
  mutate(
    mean_change = case_when(
      source == "change_from_baseline" ~ y,
      source == "baseline_followup" ~ y - b,
      source == "responders" ~ response_mean
    ),
    sd_change = case_when(
      source == "baseline_followup" ~ sqrt(
        s^2 + i^2 - 2 * assumed_correlation * s * i
      ),
      TRUE ~ s
    )
  ) %>%
  group_by(study_id) %>%
  filter(n_distinct(classcode) >= 2) %>%
  ungroup()

if (
  any(!is.finite(analysis_data$mean_change)) ||
  any(!is.finite(analysis_data$sd_change)) ||
  any(analysis_data$sd_change <= 0)
) {
  stop("Invalid analysis means or SDs. Check the data and response proportions.")
}

if (!reference_code %in% analysis_data$classcode) {
  stop("Placebo is absent from the analysis data.")
}

# -----------------------------
# 8) Create comparisons and retain the Placebo-connected network
# -----------------------------
pairwise_data <- meta::pairwise(
  treat = analysis_data$classcode,
  mean = analysis_data$mean_change,
  sd = analysis_data$sd_change,
  n = analysis_data$n,
  studlab = analysis_data$study_id,
  sm = "SMD"
)

network_connection <- netconnection(pairwise_data)

if (network_connection$n.subnets > 1) {
  reference_subnet <- unique(
    network_connection$subnet[
      pairwise_data$treat1 == reference_code |
        pairwise_data$treat2 == reference_code
    ]
  )
  
  pairwise_data <- pairwise_data[
    network_connection$subnet == reference_subnet,
  ]
  
  message("Network disconnected: retaining only the Placebo-connected component.")
}

# Match participant and study counts to the retained comparisons.
retained_classes <- unique(c(
  as.character(pairwise_data$treat1),
  as.character(pairwise_data$treat2)
))

analysis_data <- analysis_data %>%
  filter(
    study_id %in% pairwise_data$studlab,
    classcode %in% retained_classes
  )

# -----------------------------
# 9) Fit the base random-effects NMA
# -----------------------------
nma <- netmeta(
  TE = TE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = studlab,
  data = pairwise_data,
  sm = "SMD",
  common = FALSE,
  random = TRUE,
  reference.group = reference_code,
  small.values = "good"
)

# -----------------------------
# 10) Extract effects, rankings, and sample sizes
# -----------------------------
effect_tbl <- tibble(
  classcode = rownames(nma$TE.random),
  SMD = as.numeric(nma$TE.random[, reference_code]),
  lower_ci = as.numeric(nma$lower.random[, reference_code]),
  upper_ci = as.numeric(nma$upper.random[, reference_code]),
  p_value_vs_reference = as.numeric(nma$pval.random[, reference_code])
) %>%
  mutate(
    SMD = if_else(classcode == reference_code, 0, SMD),
    lower_ci = if_else(classcode == reference_code, 0, lower_ci),
    upper_ci = if_else(classcode == reference_code, 0, upper_ci),
    p_value_vs_reference = if_else(
      classcode == reference_code,
      NA_real_,
      p_value_vs_reference
    )
  )

ranking_results <- netrank(
  nma,
  method = "P-score",
  small.values = "good"
)

p_scores <- ranking_results$ranking.random

rank_tbl <- tibble(
  classcode = names(p_scores),
  P_score = as.numeric(p_scores),
  Rank = rank(-as.numeric(p_scores), ties.method = "average")
)

sample_size_tbl <- analysis_data %>%
  group_by(classcode) %>%
  summarise(
    N_total = sum(n),
    k_total = n_distinct(study_id),
    .groups = "drop"
  )

# -----------------------------
# 11) Save the results table
# -----------------------------
out_tbl <- effect_tbl %>%
  left_join(class_labels, by = "classcode") %>%
  left_join(sample_size_tbl, by = "classcode") %>%
  left_join(rank_tbl, by = "classcode") %>%
  arrange(Rank) %>%
  select(
    classcode,
    Treatment = class,
    N_total,
    k_total,
    SMD,
    lower_ci,
    upper_ci,
    p_value_vs_reference,
    P_score,
    Rank
  )

output_file <- file.path(out_dir, "nma_ms_base_random_results.csv")

write_csv(out_tbl, output_file)
message("Saved: ", output_file)