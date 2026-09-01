library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(netmeta)

# =========================
# 1. Read Excel files
# =========================

input_cfb_path <- "/Users/Fredrik/Desktop/nma_project/mavranezouli/clean_data/mmc5_cfb_ls.xlsx"
input_bf_path <- "/Users/Fredrik/Desktop/nma_project/mavranezouli/clean_data/mmc5_bf_ls.xlsx"
input_resp_path <- "/Users/Fredrik/Desktop/nma_project/mavranezouli/clean_data/mmc5_resp_ls.xlsx"

output_long_path <- "/Users/Fredrik/Desktop/nma_project/mavranezouli/ls_files/combined_long_mean_change_dataset.csv"
output_wide_path <- "/Users/Fredrik/Desktop/nma_project/mavranezouli/ls_files/combined_wide_mean_change_dataset.csv"

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
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "Na", "N/A")] <- NA
  as.numeric(x)
}

clean_weird_names <- function(nms) {
  out <- nms
  out <- str_replace_all(out, fixed("[,"), "")
  out <- str_replace_all(out, fixed("]"), "")
  out <- str_replace_all(out, fixed("["), "")
  out <- str_replace_all(out, ",", "")
  out <- str_replace_all(out, "\\s+", "")
  out <- ifelse(out == "na", "na", out)
  out <- ifelse(out == "rho", "rho", out)
  out <- ifelse(out == "q", "q", out)
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

normalize_rho_column <- function(dat) {
  if ("rho" %in% names(dat)) {
    return(dat)
  }

  if ("q" %in% names(dat)) {
    return(dat %>% rename(rho = q))
  }

  stop("Expected a study-level correlation column named 'rho' or 'q'.")
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

sd_change_from_baseline_followup <- function(sd_baseline, sd_followup, rho) {
  sqrt(sd_baseline^2 + sd_followup^2 - 2 * rho * sd_baseline * sd_followup)
}

sd_change_equal_sd <- function(sd_baseline, rho) {
  sqrt(2 * sd_baseline^2 * (1 - rho))
}

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

# =========================
# 3. User assumptions
# =========================

response_fraction <- 0.50
rho_assumed_bf <- 0.50

# =========================
# 4. Clean imported data
# =========================

dat_cfb <- clean_block(dat_cfb_raw)
dat_bf <- clean_block(dat_bf_raw)
dat_resp <- dat_resp_raw %>%
  clean_block() %>%
  normalize_rho_column()

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
    )
}

# =========================
# 6. Process baseline/follow-up block
# =========================

process_bf <- function(dat_bf, rho_assumed = 0.5) {
  wide_to_long_arms(
    dat_bf,
    prefixes = c("t", "yB", "sdB", "yF", "sdF", "n"),
    max_arms = 5
  ) %>%
    mutate(
      mean_change = yF - yB,
      sd_change = sd_change_from_baseline_followup(
        sd_baseline = sdB,
        sd_followup = sdF,
        rho = rho_assumed
      ),
      rho_used = rho_assumed
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
    )
}

# =========================
# 7. Process responder block
# =========================

process_resp <- function(dat_resp, response_fraction = 0.50) {
  wide_to_long_arms(
    dat_resp,
    prefixes = c("t", "r", "n", "yBR", "sdBR"),
    max_arms = 5
  ) %>%
    left_join(
      dat_resp %>%
        select(studyid, rho),
      by = "studyid"
    ) %>%
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
    )
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
  response_fraction = response_fraction
)

combined_long <- bind_rows(long_cfb, long_bf, long_resp)

# =========================
# 9. Optional: reshape back to one row per study
# =========================

make_arm_wide <- function(long_dat) {
  long_dat %>%
    select(studyid, arm, treatment, n, mean_change, sd_change, source) %>%
    pivot_wider(
      names_from = arm,
      values_from = c(treatment, n, mean_change, sd_change),
      names_sep = ""
    ) %>%
    arrange(studyid)
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

write.csv(combined_long, output_long_path, row.names = FALSE)
write.csv(combined_wide, output_wide_path, row.names = FALSE)
