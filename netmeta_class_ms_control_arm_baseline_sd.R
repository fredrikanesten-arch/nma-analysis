#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(netmeta)
})

read_csv_robust <- function(path) {
  dat_csv <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(dat_csv) && ncol(dat_csv) > 1) return(dat_csv)
  dat_scsv <- tryCatch(readr::read_delim(path, delim = ";", show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(dat_scsv) && ncol(dat_scsv) > 1) return(dat_scsv)
  stop("Could not parse CSV: ", path)
}

make_nma_summary <- function(nma, collapsed_class, ref = "Placebo", digits = 3) {
  if (!ref %in% nma$trts) stop("Reference treatment not found: ", ref)
  te <- if (!is.null(nma$TE.random)) nma$TE.random else nma$TE.common
  lo <- if (!is.null(nma$lower.random)) nma$lower.random else nma$lower.common
  hi <- if (!is.null(nma$upper.random)) nma$upper.random else nma$upper.common
  se <- if (!is.null(nma$seTE.random)) nma$seTE.random else nma$seTE.common
  
  pick_vs_ref <- function(m) {
    if (ref %in% colnames(m)) return(data.frame(class = rownames(m), value = as.numeric(m[, ref]), row.names = NULL))
    if (ref %in% rownames(m)) return(data.frame(class = colnames(m), value = as.numeric(m[ref, ]), row.names = NULL))
    stop("Reference not found in matrix")
  }
  
  n_tbl <- collapsed_class %>% group_by(class) %>% summarise(N = sum(n), .groups = "drop")
  out <- pick_vs_ref(te) %>%
    rename(estimate = value) %>%
    left_join(rename(pick_vs_ref(lo), lower = value), by = "class") %>%
    left_join(rename(pick_vs_ref(hi), upper = value), by = "class") %>%
    left_join(rename(pick_vs_ref(se), se = value), by = "class") %>%
    left_join(n_tbl, by = "class") %>%
    mutate(
      estimate = ifelse(class == ref, 0, estimate),
      lower = ifelse(class == ref, 0, lower),
      upper = ifelse(class == ref, 0, upper),
      se = ifelse(class == ref, NA_real_, se),
      z = ifelse(!is.na(se) & se > 0, estimate / se, NA_real_),
      p_value = ifelse(!is.na(z), 2 * pnorm(-abs(z)), NA_real_)
    ) %>%
    arrange(estimate)
  
  fmt <- function(x) formatC(x, digits = digits, format = "f")
  out$effect_95CI <- paste0(fmt(out$estimate), " (", fmt(out$lower), " to ", fmt(out$upper), ")")
  out
}

collapse_same_class <- function(dat) {
  dat %>%
    group_by(studyid, class) %>%
    summarise(
      n = sum(n),
      mean_change = sum(mean_change * n) / sum(n),
      sd_change = {
        n_total <- sum(n)
        mu <- sum(mean_change * n) / n_total
        ss_within <- sum((n - 1) * (sd_change^2))
        ss_between <- sum(n * (mean_change - mu)^2)
        sqrt((ss_within + ss_between) / (n_total - 1))
      },
      source = dplyr::first(source),
      control_class = dplyr::first(control_class),
      control_sd_ref = dplyr::first(control_sd_ref),
      control_sd_source = dplyr::first(control_sd_source),
      .groups = "drop"
    ) %>%
    group_by(studyid) %>%
    filter(n_distinct(class) >= 2) %>%
    ungroup()
}

main <- function() {
  base_dir <- "/home/runner/work/nma-analysis/nma-analysis"
  default_in_data <- file.path(base_dir, "combined_long_mean_change_dataset_ms_smd_bias_adj.csv")
  default_in_map <- file.path(base_dir, "trt_to_class_ms.csv")
  default_out_dir <- base_dir
  
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[[1]] %in% c("-h", "--help")) {
    cat("Usage: Rscript netmeta_class_ms_control_arm_baseline_sd.R [input_csv] [mapping_csv] [output_dir] [reference] [fallback]\n")
    cat("fallback: 'drop' (default) or 'sd_change' when control baseline SD is missing.\n")
    return(invisible(NULL))
  }
  
  in_data <- if (length(args) >= 1) args[[1]] else default_in_data
  in_map <- if (length(args) >= 2) args[[2]] else default_in_map
  out_dir <- if (length(args) >= 3) args[[3]] else default_out_dir
  reference_group <- if (length(args) >= 4) args[[4]] else "Placebo"
  fallback_mode <- if (length(args) >= 5) args[[5]] else "drop"
  
  if (!fallback_mode %in% c("drop", "sd_change")) {
    stop("Invalid fallback mode. Use 'drop' or 'sd_change'.")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!file.exists(in_data)) stop("Input data not found: ", in_data)
  
  dat <- read_csv_robust(in_data)
  
  if (!"class" %in% names(dat)) {
    if (!file.exists(in_map)) stop("Mapping file required when class column is absent: ", in_map)
    trt_map <- read_csv_robust(in_map)
    req_map <- c("trtcode", "class")
    if (!all(req_map %in% names(trt_map))) stop("Mapping file missing columns: trtcode, class")
    trt_map_use <- trt_map %>%
      mutate(trtcode = suppressWarnings(as.numeric(trtcode)), class = as.character(class)) %>%
      filter(!is.na(trtcode), !is.na(class), nzchar(class)) %>%
      group_by(trtcode) %>% slice(1) %>% ungroup()
    dat <- dat %>%
      mutate(treatment = suppressWarnings(as.numeric(treatment))) %>%
      left_join(trt_map_use %>% rename(treatment = trtcode), by = "treatment")
  }
  
  req <- c("studyid", "class", "n", "mean_change", "sd_change")
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) stop("Input data missing required columns: ", paste(miss, collapse = ", "))
  
  if (!"sd_baseline" %in% names(dat)) dat$sd_baseline <- NA_real_
  if (!"source" %in% names(dat)) dat$source <- "unknown"
  
  dat <- dat %>%
    mutate(
      studyid = as.character(studyid),
      class = as.character(class),
      n = suppressWarnings(as.numeric(n)),
      mean_change = suppressWarnings(as.numeric(mean_change)),
      sd_change = suppressWarnings(as.numeric(sd_change)),
      sd_baseline = suppressWarnings(as.numeric(sd_baseline)),
      source = as.character(source)
    ) %>%
    filter(!is.na(studyid), !is.na(class), nzchar(class))
  
  control_priority <- c("Placebo", "Attention placebo", "No treatment", "Waitlist", "TAU")
  
  control_pick <- dat %>%
    filter(class %in% control_priority, !is.na(n), n > 1, !is.na(sd_change), sd_change > 0) %>%
    mutate(
      control_order = match(class, control_priority),
      usable_baseline = !is.na(sd_baseline) & sd_baseline > 0
    ) %>%
    arrange(studyid, control_order, desc(usable_baseline), desc(n)) %>%
    group_by(studyid) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      control_sd_ref = dplyr::case_when(
        usable_baseline ~ sd_baseline,
        fallback_mode == "sd_change" ~ sd_change,
        TRUE ~ NA_real_
      ),
      control_sd_source = dplyr::case_when(
        usable_baseline ~ "sd_baseline",
        fallback_mode == "sd_change" ~ "sd_change_fallback",
        TRUE ~ "missing"
      )
    ) %>%
    select(studyid, control_class = class, control_sd_ref, control_sd_source)
  
  dat_std <- dat %>%
    left_join(control_pick, by = "studyid") %>%
    mutate(
      drop_reason = case_when(
        is.na(control_class) ~ "no_control_arm_found",
        is.na(control_sd_ref) ~ "control_sd_not_available",
        is.na(n) | n <= 1 ~ "invalid_n",
        is.na(mean_change) ~ "missing_mean_change",
        is.na(sd_change) | sd_change <= 0 ~ "invalid_sd_change",
        TRUE ~ "kept"
      )
    )
  
  write_csv(
    dat_std %>% count(drop_reason, sort = TRUE),
    file.path(out_dir, "drop_reason_counts_control_arm_baseline_sd.csv")
  )
  write_csv(
    dat_std %>% distinct(studyid, control_class, control_sd_ref, control_sd_source) %>% arrange(studyid),
    file.path(out_dir, "control_arm_selection_control_arm_baseline_sd.csv")
  )
  
  dat_std <- dat_std %>%
    filter(drop_reason == "kept") %>%
    mutate(
      mean_change = mean_change / control_sd_ref,
      sd_change = sd_change / control_sd_ref
    )
  
  collapsed_class <- collapse_same_class(dat_std)
  write_csv(collapsed_class, file.path(out_dir, "collapsed_class_control_arm_baseline_sd.csv"))
  
  if (!reference_group %in% collapsed_class$class) {
    stop("Reference group '", reference_group, "' is not present after filtering.")
  }
  
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
    reference.group = reference_group
  )
  
  summary_tbl <- make_nma_summary(nma, collapsed_class, ref = reference_group)
  write_csv(summary_tbl, file.path(out_dir, "nma_summary_table_control_arm_baseline_sd.csv"))
  
  saveRDS(
    list(input = dat, control_pick = control_pick, standardized = dat_std, collapsed_class = collapsed_class, pw = pw, nma = nma),
    file.path(out_dir, "nma_control_arm_baseline_sd_objects.rds")
  )
  
  cat("Done.\n")
  cat("Wrote outputs to: ", out_dir, "\n", sep = "")
}

main()
