suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

read_csv_robust <- function(path) {
  dat_csv <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(dat_csv) && ncol(dat_csv) > 1) return(dat_csv)
  dat_scsv <- tryCatch(readr::read_delim(path, delim = ";", show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(dat_scsv) && ncol(dat_scsv) > 1) return(dat_scsv)
  stop("Could not parse CSV: ", path)
}

derive_arm_flags <- function(source_chr) {
  src <- tolower(as.character(source_chr))
  data.frame(
    src_change_from_baseline = as.integer(src %in% c("change_from_baseline", "cfb", "smd")),
    src_baseline_followup = as.integer(src %in% c("baseline_followup", "basefu")),
    src_binary_response = as.integer(src %in% c("binary_response_converted", "responders", "response"))
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[[1]] %in% c("-h", "--help")) {
    cat("Usage: Rscript build_mmc5_ms_smd_bias_adj_dataset.R [input_csv] [mapping_csv] [output_csv]\n")
    cat("Adds class mapping + treatment-level + bridge eligibility flags to long dataset.\n")
    return(invisible(NULL))
  }

  in_data <- if (length(args) >= 1) args[[1]] else "combined_long_mean_change_dataset.csv"
  in_map <- if (length(args) >= 2) args[[2]] else "trt_to_class_ms.csv"
  out_csv <- if (length(args) >= 3) args[[3]] else "combined_long_mean_change_dataset_ms_smd_bias_adj.csv"

  if (!file.exists(in_data)) stop("Input data not found: ", in_data)
  if (!file.exists(in_map)) stop("Mapping file not found: ", in_map)

  dat <- read_csv_robust(in_data)
  trt_map <- read_csv_robust(in_map)

  req_data <- c("studyid", "treatment", "n", "mean_change", "sd_change", "source")
  miss_data <- setdiff(req_data, names(dat))
  if (length(miss_data) > 0) stop("Input data missing columns: ", paste(miss_data, collapse = ", "))

  req_map <- c("trtcode", "class")
  miss_map <- setdiff(req_map, names(trt_map))
  if (length(miss_map) > 0) stop("Mapping file missing columns: ", paste(miss_map, collapse = ", "))

  map_clean <- trt_map %>%
    transmute(
      treatment = suppressWarnings(as.numeric(trtcode)),
      class = as.character(class),
      treatment_label = dplyr::coalesce(
        dplyr::if_else("treatment_label" %in% names(trt_map), as.character(trt_map$treatment_label), NA_character_),
        dplyr::if_else("trtname" %in% names(trt_map), as.character(trt_map$trtname), NA_character_),
        dplyr::if_else("name" %in% names(trt_map), as.character(trt_map$name), NA_character_)
      )
    ) %>%
    filter(!is.na(treatment), !is.na(class), nzchar(class)) %>%
    group_by(treatment) %>%
    slice(1) %>%
    ungroup()

  out <- dat %>%
    mutate(
      treatment = suppressWarnings(as.numeric(treatment)),
      studyid = as.character(studyid),
      source = as.character(source)
    ) %>%
    left_join(map_clean, by = "treatment") %>%
    mutate(
      treatment_node = dplyr::coalesce(
        dplyr::na_if(treatment_label, ""),
        ifelse(!is.na(treatment), paste0("trt_", sprintf("%03d", as.integer(treatment))), NA_character_)
      ),
      class = dplyr::coalesce(class, "UNMAPPED"),
      bridge_any_ad = as.integer(class == "Any AD"),
      bridge_any_psychotherapy = as.integer(class == "Any psychotherapy"),
      bridge_either = as.integer(bridge_any_ad == 1 | bridge_any_psychotherapy == 1),
      bridge_set_any_ad_only = as.integer(bridge_any_ad == 1),
      bridge_set_any_ad_any_psychotherapy = as.integer(bridge_either == 1)
    )

  arm_flags <- derive_arm_flags(out$source)
  out <- bind_cols(out, arm_flags)

  write_csv(out, out_csv)

  summary_path <- paste0(tools::file_path_sans_ext(out_csv), "_bridge_flag_summary.csv")
  write_csv(
    out %>%
      summarise(
        n_rows = n(),
        n_studies = n_distinct(studyid),
        n_unmapped_rows = sum(class == "UNMAPPED"),
        n_bridge_any_ad_rows = sum(bridge_any_ad),
        n_bridge_any_psychotherapy_rows = sum(bridge_any_psychotherapy),
        n_bridge_either_rows = sum(bridge_either)
      ),
    summary_path
  )

  cat("Done.\n")
  cat("Wrote: ", out_csv, "\n", sep = "")
  cat("Wrote: ", summary_path, "\n", sep = "")
}

main()
