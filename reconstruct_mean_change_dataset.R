library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(netmeta)

# =========================
# 0. Paths
# =========================

base_dir <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli"

input_cfb_path <- file.path(base_dir, "clean_data", "mmc5_cfb_ms.xlsx")
input_bf_path <- file.path(base_dir, "clean_data", "mmc5_bf_ms.xlsx")
input_resp_path <- file.path(base_dir, "clean_data", "mmc5_resp_ms.xlsx")

output_dir <- file.path(base_dir, "netmeta_class_ms")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_long_path <- file.path(output_dir, "combined_long_mean_change_dataset.csv")
output_wide_path <- file.path(output_dir, "combined_wide_mean_change_dataset.csv")

# =========================
# 1. Read Excel files
# =========================

dat_cfb_raw <- readxl::read_excel(input_cfb_path, col_names = TRUE)
dat_bf_raw <- readxl::read_excel(input_bf_path, col_names = TRUE)
dat_resp_raw <- readxl::read_excel(input_resp_path, col_names = TRUE)

names(dat_cfb_raw)
names(dat_bf_raw)
names(dat_resp_raw)

# =========================
# 2. Helper functions
# =========================

comma_num <- function(x) {
  if (is.numeric(x)) return(x)

  x <- as.character(x)
  x <- trimws(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  x[x %in% c("", "NA", "Na", "N/A", "n/a", "-", "–", ".", "..")] <- NA

  suppressWarnings(as.numeric(x))
}

clean_weird_names <- function(nms) {
  out <- nms
  out <- str_replace_all(out, fixed("[,"), "")
  out <- str_replace_all(out, fixed("]"), "")
  out <- str_replace_all(out, fixed("["), "")
  out <- str_replace_all(out, ",", "")
  out <- str_replace_all(out, "\\s+", "")
  out <- ifelse(out == "#", "hash", out)
  out
}

clean_block <- function(dat) {
  names(dat) <- clean_weird_names(names(dat))

  dat <- dat %>%
    mutate(across(everything(), ~ if (is.character(.x)) trimws(.x) else .x)) %>%
    filter(if_any(everything(), ~ !is.na(.x) & .x != ""))

  num_cols <- setdiff(names(dat), c("studyid", "hash"))
  dat %>%
    mutate(across(all_of(num_cols), comma_num))
}

fill_merged_cells <- function(dat) {
  out <- dat

  if ("studyid" %in% names(out)) {
    out <- out %>% tidyr::fill(studyid, .direction = "down")
  }

  if ("na" %in% names(out)) {
    out <- out %>% tidyr::fill(na, .direction = "down")
  }

  out
}

normalize_rho_columns <- function(dat) {
  out <- dat

  if ("q" %in% names(out) && !"rho" %in% names(out)) {
    out <- out %>% rename(rho = q)
  }

  rho_arm_cols <- intersect(names(out), paste0("q", seq_len(5)))
  for (col_name in rho_arm_cols) {
    new_name <- sub("^q", "rho", col_name)
    if (!new_name %in% names(out)) {
      names(out)[names(out) == col_name] <- new_name
    }
  }

  out
}

find_non_numeric_values <- function(dat, exclude = c("studyid", "hash")) {
  num_cols <- setdiff(names(dat), exclude)

  purrr::map_dfr(num_cols, function(col) {
    vals <- dat[[col]]
    vals <- as.character(vals)
    vals <- trimws(vals)
    vals2 <- gsub(",", ".", vals, fixed = TRUE)

    bad <- vals[
      !is.na(vals) &
        vals != "" &
        is.na(suppressWarnings(as.numeric(vals2))) &
        !vals %in% c("NA", "Na", "N/A", "n/a", "-", "–", ".", "..")
    ]

    tibble::tibble(
      column = col,
      bad_value = unique(bad)
    )
  })
}

wide_to_long_arms <- function(dat, prefixes, max_arms = 5) {
  long_list <- lapply(seq_len(max_arms), function(k) {
    cols_k <- paste0(prefixes, k)
    out <- dat %>%
      transmute(
        studyid = studyid,
        na = na,
        arm = k,
        across(any_of(cols_k))
      )
    names(out) <- c("studyid", "na", "arm", prefixes)
    out
  })

  bind_rows(long_list) %>%
    filter(!is.na(t))
}

extract_study_rho_lookup <- function(dat, fill_missing_with = NA_real_) {
  if (!"rho" %in% names(dat)) {
    stop("Expected a study-level correlation column named 'rho'.")
  }

  out <- dat %>%
    filter(!is.na(studyid)) %>%
    group_by(studyid) %>%
    tidyr::fill(rho, .direction = "downup") %>%
    summarise(
      rho_values = list(sort(unique(stats::na.omit(rho)))),
      .groups = "drop"
    ) %>%
    mutate(n_rho = lengths(rho_values))

  rho_conflicts <- out %>%
    filter(n_rho > 1)

  if (nrow(rho_conflicts) > 0) {
    stop(
      "Some studies have more than one non-missing rho value in responder data: ",
      paste(rho_conflicts$studyid, collapse = ", ")
    )
  }

  out %>%
    transmute(
      studyid = studyid,
      rho = dplyr::case_when(
        n_rho == 1 ~ purrr::map_dbl(rho_values, 1),
        is.finite(fill_missing_with) ~ fill_missing_with,
        TRUE ~ NA_real_
      )
    )
}

attach_arm_or_study_rho <- function(long_dat,
                                    dat_source,
                                    max_arms = 5,
                                    rho_assumed = NA_real_) {
  rho_arm_cols <- intersect(names(dat_source), paste0("rho", seq_len(max_arms)))

  if (length(rho_arm_cols) > 0) {
    rho_long <- wide_to_long_arms(
      dat_source,
      prefixes = "rho",
      max_arms = max_arms
    ) %>%
      select(studyid, arm, rho)

    return(
      long_dat %>%
        left_join(rho_long, by = c("studyid", "arm")) %>%
        mutate(rho = dplyr::coalesce(rho, rho_assumed))
    )
  }

  if ("rho" %in% names(dat_source)) {
    rho_lookup <- extract_study_rho_lookup(
      dat_source,
      fill_missing_with = rho_assumed
    )

    return(
      long_dat %>%
        left_join(rho_lookup, by = "studyid")
    )
  }

  long_dat %>%
    mutate(rho = rho_assumed)
}

sd_change_from_baseline_followup <- function(sd_baseline, sd_followup, rho) {
  sqrt(sd_baseline^2 + sd_followup^2 - 2 * rho * sd_baseline * sd_followup)
}

sd_change_equal_sd <- function(sd_baseline, rho) {
  sqrt(2 * sd_baseline^2 * (1 - rho))
}

# Mavranezouli et al. supplementary derivation:
# R_ik = Pr(Y_ik - X_ik <= -q X_ik)
# theta_ik = -q * mu_X,ik - Phi^{-1}(R_ik) * sigma_X,ik *
#            sqrt(q^2 + 2 * (1 - q) * (1 - rho_ik))
response_to_mean_change_mavranezouli <- function(p_response,
                                                 mu_baseline,
                                                 sd_baseline,
                                                 rho,
                                                 response_fraction) {
  p_response <- pmin(pmax(p_response, 1e-6), 1 - 1e-6)

  scale_term <- sd_baseline * sqrt(
    response_fraction^2 + 2 * (1 - response_fraction) * (1 - rho)
  )

  -response_fraction * mu_baseline - qnorm(p_response) * scale_term
}

assert_unique_study_source_arm <- function(dat, label) {
  dup_keys <- dat %>%
    count(studyid, source, arm, name = "n_rows") %>%
    filter(n_rows > 1)

  if (nrow(dup_keys) > 0) {
    stop(
      label,
      " created duplicate (studyid, source, arm) rows. ",
      "Inspect: count(", label, ", studyid, source, arm) |> filter(n_rows > 1)"
    )
  }

  dat
}

# =========================
# 3. User assumptions
# =========================

response_fraction <- 0.50
rho_assumed_bf <- 0.50
rho_assumed_resp <- 0.50

# =========================
# 4. Clean imported data
# =========================

dat_cfb <- dat_cfb_raw %>%
  clean_block() %>%
  fill_merged_cells()

dat_bf <- dat_bf_raw %>%
  clean_block() %>%
  fill_merged_cells() %>%
  normalize_rho_columns()

dat_resp <- dat_resp_raw %>%
  clean_block() %>%
  fill_merged_cells() %>%
  normalize_rho_columns()

# Optional diagnostics:
# find_non_numeric_values(dat_cfb_raw)
# find_non_numeric_values(dat_bf_raw)
# find_non_numeric_values(dat_resp_raw)

# =========================
# 5. Process change-from-baseline block
# =========================

process_cfb <- function(dat_cfb) {
  wide_to_long_arms(
    dat_cfb,
    prefixes = c("t", "yCFB", "sdCFB", "nCFB"),
    max_arms = 5
  ) %>%
    transmute(
      studyid = studyid,
      na = na,
      arm = arm,
      treatment = t,
      n = nCFB,
      mean_change = yCFB,
      sd_change = sdCFB,
      source = "change_from_baseline"
    ) %>%
    assert_unique_study_source_arm("process_cfb()")
}

# =========================
# 6. Process baseline/follow-up block
# =========================

process_bf <- function(dat_bf, rho_assumed = 0.5) {
  bf_prefixes <- c("t", "yB", "sdB", "yF", "sdF", "n")
  if (any(paste0("rho", seq_len(5)) %in% names(dat_bf))) {
    bf_prefixes <- c(bf_prefixes, "rho")
  }

  out <- wide_to_long_arms(
    dat_bf,
    prefixes = bf_prefixes,
    max_arms = 5
  )

  if (!"rho" %in% names(out)) {
    out <- attach_arm_or_study_rho(
      long_dat = out,
      dat_source = dat_bf,
      max_arms = 5,
      rho_assumed = rho_assumed
    )
  } else {
    out <- out %>%
      mutate(rho = dplyr::coalesce(rho, rho_assumed))
  }

  out %>%
    mutate(
      mean_change = yF - yB,
      sd_change = sd_change_from_baseline_followup(
        sd_baseline = sdB,
        sd_followup = sdF,
        rho = rho
      ),
      rho_used = rho
    ) %>%
    transmute(
      studyid = studyid,
      na = na,
      arm = arm,
      treatment = t,
      n = n,
      mean_change = mean_change,
      sd_change = sd_change,
      source = "baseline_followup",
      y_baseline = yB,
      sd_baseline = sdB,
      y_followup = yF,
      sd_followup = sdF,
      rho_used = rho_used
    ) %>%
    assert_unique_study_source_arm("process_bf()")
}

# =========================
# 7. Process responder block
# =========================

process_resp <- function(dat_resp,
                         response_fraction = 0.50,
                         rho_assumed = 0.50) {
  resp_prefixes <- c("t", "r", "n", "yBR", "sdBR")
  if (any(paste0("rho", seq_len(5)) %in% names(dat_resp))) {
    resp_prefixes <- c(resp_prefixes, "rho")
  }

  out <- wide_to_long_arms(
    dat_resp,
    prefixes = resp_prefixes,
    max_arms = 5
  )

  if (!"rho" %in% names(out)) {
    out <- attach_arm_or_study_rho(
      long_dat = out,
      dat_source = dat_resp,
      max_arms = 5,
      rho_assumed = rho_assumed
    )
  } else {
    out <- out %>%
      group_by(studyid) %>%
      tidyr::fill(rho, .direction = "downup") %>%
      ungroup() %>%
      mutate(rho = dplyr::coalesce(rho, rho_assumed))
  }

  out %>%
    mutate(
      p_response = r / n,
      sd_change = sd_change_equal_sd(sdBR, rho),
      mean_change = response_to_mean_change_mavranezouli(
        p_response = p_response,
        mu_baseline = yBR,
        sd_baseline = sdBR,
        rho = rho,
        response_fraction = response_fraction
      )
    ) %>%
    transmute(
      studyid = studyid,
      na = na,
      arm = arm,
      treatment = t,
      n = n,
      responders = r,
      p_response = p_response,
      mean_change = mean_change,
      sd_change = sd_change,
      source = "binary_response_converted",
      y_baseline = yBR,
      sd_baseline = sdBR,
      rho = rho,
      response_fraction = response_fraction
    ) %>%
    assert_unique_study_source_arm("process_resp()")
}

# =========================
# 8. Run all transformations
# =========================

long_cfb <- process_cfb(dat_cfb)

long_bf <- process_bf(
  dat_bf,
  rho_assumed = rho_assumed_bf
)

long_resp <- process_resp(
  dat_resp,
  response_fraction = response_fraction,
  rho_assumed = rho_assumed_resp
)

combined_long <- bind_rows(long_cfb, long_bf, long_resp)

# =========================
# 9. Optional: reshape back to one row per study/source
# =========================

make_arm_wide <- function(long_dat) {
  long_dat %>%
    assert_unique_study_source_arm("make_arm_wide() input") %>%
    select(studyid, source, arm, treatment, n, mean_change, sd_change) %>%
    pivot_wider(
      id_cols = c(studyid, source),
      names_from = arm,
      values_from = c(treatment, n, mean_change, sd_change),
      names_sep = ""
    ) %>%
    arrange(studyid, source)
}

combined_wide <- make_arm_wide(combined_long)

# =========================
# 10. Inspect / save
# =========================

print(dat_cfb)
print(dat_bf)
print(dat_resp)

print(long_cfb)
print(long_bf)
print(long_resp)

print(combined_long)
print(combined_wide)

readr::write_csv(combined_long, output_long_path)
readr::write_csv(combined_wide, output_wide_path)
