#!/usr/bin/env Rscript

DEFAULT_INPUT_DIR <- "C:\\Users\\fredr\\OneDrive\\Desktop\\nma_project\\mavranezouli\\clean_data"
DEFAULT_OUTPUT_DIR <- "C:\\Users\\fredr\\OneDrive\\Desktop\\nma_project\\mavranezouli\\placebo_reclassified"
DEFAULT_SHEET <- "MS SMD bias-adj"
PLACEBO_CODE <- 1
PARTIAL_PLACEBO_CODE <- 100
PARTIAL_PLACEBO_NAME <- "Partial placebo"
PLACEBO_CLASS_CODE <- 1
PLACEBO_CLASS_NAME <- "Placebo"
ARM_COLUMNS <- sprintf("Arm %d intervention", 1:5)
CONTROL_ARMS <- c("pill placebo", "attention placebo", "no treatment", "waitlist", "tau")
# Heuristic list of non-pharmacological intervention terms used to identify
# placebo comparisons that may need partial-placebo reclassification; extend
# this list when new non-drug intervention labels appear in the source workbooks.
NONPHARMA_KEYWORDS <- c(
  "acupuncture",
  "attentional bias",
  "behavioral",
  "behavioural",
  "bibliotherapy",
  "bright light",
  "cbt",
  "cognitive bias",
  "computerised",
  "computerized",
  "counselling",
  "counseling",
  "dbt",
  "dialectical",
  "exercise",
  "interpersonal counselling",
  "interpersonal psychotherapy",
  "light therapy",
  "mbct",
  "meditation",
  "mindfulness",
  "music therapy",
  "peer support",
  "positive psychological",
  "problem solving",
  "psychoeducational",
  "psychodynamic",
  "psychotherapy",
  "relaxation",
  "self-help",
  "supportive",
  "website",
  "yoga"
)

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0 || (length(left) == 1 && is.na(left))) right else left
}

normalize_string <- function(value) {
  if (length(value) == 0 || is.na(value)) {
    return("")
  }
  gsub("\\s+", " ", trimws(tolower(as.character(value))))
}

is_blank_cell <- function(value) {
  if (length(value) == 0 || is.na(value)) {
    return(TRUE)
  }
  identical(trimws(as.character(value)), "")
}

is_placeholder_missing <- function(value) {
  is_blank_cell(value) || identical(normalize_string(value), "na")
}

as_numeric_code <- function(value) {
  if (is_placeholder_missing(value)) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(value))
}

split_arm_components <- function(arm) {
  components <- trimws(unlist(strsplit(as.character(arm), "\\+")))
  components[nzchar(components)]
}

is_nonpharma_component <- function(component) {
  lowered <- normalize_string(component)
  if (lowered %in% CONTROL_ARMS) {
    return(FALSE)
  }
  any(vapply(NONPHARMA_KEYWORDS, function(keyword) grepl(keyword, lowered, fixed = TRUE), logical(1)))
}

has_nonpharma_arm <- function(arms) {
  active_arms <- arms[!vapply(arms, is_blank_cell, logical(1))]
  any(vapply(active_arms, function(arm) {
    any(vapply(split_arm_components(arm), is_nonpharma_component, logical(1)))
  }, logical(1)))
}

has_pill_placebo_arm <- function(arms) {
  any(vapply(arms, function(arm) normalize_string(arm) == "pill placebo", logical(1)))
}

has_blinding_issue <- function(performance_bias, detection_bias) {
  normalize_string(performance_bias) != "low risk" || normalize_string(detection_bias) != "low risk"
}

infer_mmc3_sheet <- function(mmc5_sheet) {
  prefix <- strsplit(as.character(mmc5_sheet), " ")[[1]][1]
  if (identical(prefix, "MS")) {
    return("MS depression-included studies")
  }
  if (identical(prefix, "LS")) {
    return("LS depression-included studies")
  }
  stop(sprintf("Cannot infer mmc3 sheet from '%s'. Supply --mmc3-sheet.", mmc5_sheet), call. = FALSE)
}

resolve_mmc3_sheet_name <- function(available_sheets, requested_sheet) {
  aliases <- c(
    "LS depression-included studies" = "LS depression -included studies",
    "LS depression -included studies" = "LS depression-included studies"
  )
  if (requested_sheet %in% available_sheets) {
    return(requested_sheet)
  }
  alias <- aliases[[requested_sheet]]
  if (!is.null(alias) && alias %in% available_sheets) {
    return(alias)
  }
  requested_sheet
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- list(
    input_dir = DEFAULT_INPUT_DIR,
    output_dir = DEFAULT_OUTPUT_DIR,
    mmc5 = NULL,
    mmc3 = NULL,
    lookup = NULL,
    sheet = DEFAULT_SHEET,
    mmc3_sheet = NULL
  )

  index <- 1
  while (index <= length(args)) {
    token <- args[[index]]
    if (!startsWith(token, "--")) {
      stop(sprintf("Unexpected argument '%s'. Use --name=value or --name value.", token), call. = FALSE)
    }

    if (grepl("=", token, fixed = TRUE)) {
      parts <- strsplit(substring(token, 3), "=", fixed = TRUE)[[1]]
      name <- parts[[1]]
      value <- paste(parts[-1], collapse = "=")
    } else {
      name <- substring(token, 3)
      index <- index + 1
      if (index > length(args) || startsWith(args[[index]], "--")) {
        stop(sprintf("Missing value for --%s.", name), call. = FALSE)
      }
      value <- args[[index]]
    }

    normalized_name <- gsub("-", "_", name)
    if (!normalized_name %in% names(options)) {
      stop(sprintf("Unknown option --%s.", name), call. = FALSE)
    }
    options[[normalized_name]] <- value
    index <- index + 1
  }

  if (is.null(options$mmc5)) {
    options$mmc5 <- file.path(options$input_dir, "mmc5_fixed.xlsx")
  }
  if (is.null(options$mmc3)) {
    options$mmc3 <- file.path(options$input_dir, "mmc3_included_studies.xlsx")
  }

  options
}

assert_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
}

load_study_sheet <- function(path, sheet_name) {
  assert_file_exists(path)
  available_sheets <- openxlsx::getSheetNames(path)
  sheet_name <- resolve_mmc3_sheet_name(available_sheets, sheet_name)
  if (!sheet_name %in% available_sheets) {
    stop(sprintf("Sheet '%s' not found in %s.", sheet_name, path), call. = FALSE)
  }

  study_df <- openxlsx::read.xlsx(path, sheet = sheet_name, colNames = TRUE)
  canonical_required_columns <- c(
    "Study ID",
    "Blinding of participants and personnel (performance bias)",
    "Blinding of outcome assessment (detection bias)"
  )
  normalized_names <- make.names(names(study_df), unique = TRUE)
  canonical_map <- stats::setNames(names(study_df), normalized_names)
  required_columns <- make.names(canonical_required_columns, unique = TRUE)
  missing_columns <- required_columns[!required_columns %in% normalized_names]
  if (length(missing_columns) > 0) {
    stop(
      sprintf(
        "Sheet '%s' in %s is missing required columns: %s",
        sheet_name,
        path,
        paste(canonical_required_columns[match(missing_columns, required_columns)], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  study_id_column <- canonical_map[[make.names("Study ID", unique = TRUE)]]
  performance_column <- canonical_map[[make.names("Blinding of participants and personnel (performance bias)", unique = TRUE)]]
  detection_column <- canonical_map[[make.names("Blinding of outcome assessment (detection bias)", unique = TRUE)]]
  arm_columns <- vapply(
    make.names(ARM_COLUMNS, unique = TRUE),
    function(name) if (name %in% names(canonical_map)) canonical_map[[name]] else NA_character_,
    character(1)
  )

  records <- list()
  for (row_index in seq_len(nrow(study_df))) {
    study_id <- trimws(as.character(study_df[[study_id_column]][[row_index]]))
    if (!nzchar(study_id)) {
      next
    }

    arms <- rep(NA_character_, length(ARM_COLUMNS))
    names(arms) <- ARM_COLUMNS
    for (arm_index in seq_along(ARM_COLUMNS)) {
      column_name <- arm_columns[[arm_index]]
      if (is.na(column_name)) {
        next
      }
      arm_value <- study_df[[column_name]][[row_index]]
      if (!is_blank_cell(arm_value)) {
        arms[[arm_index]] <- trimws(as.character(arm_value))
      }
    }

    records[[study_id]] <- list(
      study_id = study_id,
      arms = arms,
      performance_bias = trimws(as.character(study_df[[performance_column]][[row_index]] %||% "")),
      detection_bias = trimws(as.character(study_df[[detection_column]][[row_index]] %||% ""))
    )
  }

  list(records = records, sheet_name = sheet_name)
}

read_raw_sheet <- function(path, sheet_name) {
  assert_file_exists(path)
  available_sheets <- openxlsx::getSheetNames(path)
  if (!sheet_name %in% available_sheets) {
    stop(
      sprintf("Sheet '%s' not found in %s. Available sheets: %s", sheet_name, path, paste(available_sheets, collapse = ", ")),
      call. = FALSE
    )
  }

  openxlsx::readWorkbook(
    path,
    sheet = sheet_name,
    colNames = FALSE,
    rowNames = FALSE,
    skipEmptyRows = FALSE,
    skipEmptyCols = FALSE
  )
}

find_block_headers <- function(raw_sheet) {
  which(vapply(seq_len(nrow(raw_sheet)), function(row_index) {
    values <- unlist(raw_sheet[row_index, , drop = TRUE], use.names = FALSE)
    any(vapply(values, function(value) identical(as.character(value), "na[]"), logical(1)))
  }, logical(1)))
}

build_column_map <- function(raw_sheet, header_row) {
  headers <- unlist(raw_sheet[header_row, , drop = TRUE], use.names = FALSE)
  indices <- which(!vapply(headers, is_blank_cell, logical(1)))
  stats::setNames(indices, as.character(headers[indices]))
}

row_is_blank <- function(raw_sheet, row_index) {
  values <- unlist(raw_sheet[row_index, , drop = TRUE], use.names = FALSE)
  all(vapply(values, is_blank_cell, logical(1)))
}

list_target_blocks <- function(raw_sheet) {
  header_rows <- find_block_headers(raw_sheet)
  blocks <- list()
  block_index <- 0

  for (header_pos in seq_along(header_rows)) {
    header_row <- header_rows[[header_pos]]
    column_map <- build_column_map(raw_sheet, header_row)
    if (any(startsWith(names(column_map), "r["))) {
      next
    }

    next_header <- if (header_pos < length(header_rows)) header_rows[[header_pos + 1]] else nrow(raw_sheet) + 1
    start_row <- header_row + 1
    end_row <- next_header - 1
    while (end_row >= start_row && row_is_blank(raw_sheet, end_row)) {
      end_row <- end_row - 1
    }

    block_index <- block_index + 1
    blocks[[block_index]] <- list(
      block_index = block_index,
      header_row = header_row,
      start_row = start_row,
      end_row = end_row,
      column_map = column_map
    )
  }

  blocks
}

study_status <- function(study_record) {
  pill_placebo <- has_pill_placebo_arm(study_record$arms)
  nonpharma <- has_nonpharma_arm(study_record$arms)
  blinding_issue <- has_blinding_issue(study_record$performance_bias, study_record$detection_bias)

  if (!pill_placebo) {
    return(list(status = "not_reclassified", reason = "mmc3_has_no_pill_placebo_arm", pill_placebo = pill_placebo, nonpharma = nonpharma, blinding_issue = blinding_issue))
  }
  if (!nonpharma) {
    return(list(status = "not_reclassified", reason = "no_nonpharmacological_component_detected", pill_placebo = pill_placebo, nonpharma = nonpharma, blinding_issue = blinding_issue))
  }
  if (!blinding_issue) {
    return(list(status = "not_reclassified", reason = "no_blinding_issue_detected", pill_placebo = pill_placebo, nonpharma = nonpharma, blinding_issue = blinding_issue))
  }

  list(status = "reclassified", reason = "reclassified", pill_placebo = pill_placebo, nonpharma = nonpharma, blinding_issue = blinding_issue)
}

build_audit_row <- function(sheet_name, block_index, row_index, study_id, mmc3_sheet_name, placebo_columns, performance_bias, detection_bias, has_pill_placebo_arm, has_nonpharmacological_component, has_blinding_issue_value, status, reason, arms) {
  data.frame(
    sheet_name = sheet_name,
    block_index = block_index,
    worksheet_row = row_index,
    study_id = study_id,
    matched_sheet = mmc3_sheet_name,
    treat_columns_with_code_1 = paste(placebo_columns, collapse = " | "),
    performance_bias = performance_bias,
    detection_bias = detection_bias,
    has_pill_placebo_arm = has_pill_placebo_arm,
    has_nonpharmacological_component = has_nonpharmacological_component,
    has_blinding_issue = has_blinding_issue_value,
    status = status,
    reason = reason,
    arms = arms,
    stringsAsFactors = FALSE
  )
}

empty_flagged_frame <- function() {
  data.frame(
    sheet_name = character(),
    block_index = integer(),
    worksheet_row = integer(),
    study_id = character(),
    matched_sheet = character(),
    replaced_treat_columns = character(),
    performance_bias = character(),
    detection_bias = character(),
    arms = character(),
    original_treat_code = integer(),
    replacement_treat_code = integer(),
    stringsAsFactors = FALSE
  )
}

empty_audit_frame <- function() {
  data.frame(
    sheet_name = character(),
    block_index = integer(),
    worksheet_row = integer(),
    study_id = character(),
    matched_sheet = character(),
    treat_columns_with_code_1 = character(),
    performance_bias = character(),
    detection_bias = character(),
    has_pill_placebo_arm = logical(),
    has_nonpharmacological_component = logical(),
    has_blinding_issue = logical(),
    status = character(),
    reason = character(),
    arms = character(),
    stringsAsFactors = FALSE
  )
}

collect_reclassification_results <- function(mmc5_path, mmc3_records, sheet_name, mmc3_sheet_name) {
  workbook <- openxlsx::loadWorkbook(mmc5_path)
  raw_sheet <- read_raw_sheet(mmc5_path, sheet_name)
  target_blocks <- list_target_blocks(raw_sheet)
  flagged_rows <- list()
  audit_rows <- list()

  for (block in target_blocks) {
    column_map <- block$column_map
    treat_column_names <- sprintf("t[,%d]", 1:5)
    treat_column_names <- treat_column_names[treat_column_names %in% names(column_map)]
    treat_columns <- unname(column_map[treat_column_names])
    study_column <- unname(column_map[["studyid"]])
    if (length(treat_columns) == 0 || length(study_column) == 0 || is.na(study_column)) {
      next
    }

    for (row_index in seq.int(block$start_row, block$end_row)) {
      row_codes <- vapply(treat_columns, function(column_index) as_numeric_code(raw_sheet[[column_index]][[row_index]]), numeric(1))
      placebo_columns <- treat_column_names[!is.na(row_codes) & row_codes == PLACEBO_CODE]
      if (length(placebo_columns) == 0) {
        next
      }

      study_id <- trimws(as.character(raw_sheet[[study_column]][[row_index]] %||% ""))
      if (!nzchar(study_id) || is_placeholder_missing(study_id)) {
        audit_rows[[length(audit_rows) + 1]] <- build_audit_row(
          sheet_name = sheet_name,
          block_index = block$block_index,
          row_index = row_index,
          study_id = NA_character_,
          mmc3_sheet_name = mmc3_sheet_name,
          placebo_columns = placebo_columns,
          performance_bias = NA_character_,
          detection_bias = NA_character_,
          has_pill_placebo_arm = NA,
          has_nonpharmacological_component = NA,
          has_blinding_issue_value = NA,
          status = "manual_review",
          reason = "missing_study_id_in_mmc5",
          arms = NA_character_
        )
        next
      }

      study_record <- mmc3_records[[study_id]]
      if (is.null(study_record)) {
        audit_rows[[length(audit_rows) + 1]] <- build_audit_row(
          sheet_name = sheet_name,
          block_index = block$block_index,
          row_index = row_index,
          study_id = study_id,
          mmc3_sheet_name = mmc3_sheet_name,
          placebo_columns = placebo_columns,
          performance_bias = NA_character_,
          detection_bias = NA_character_,
          has_pill_placebo_arm = NA,
          has_nonpharmacological_component = NA,
          has_blinding_issue_value = NA,
          status = "manual_review",
          reason = "study_not_found_in_mmc3",
          arms = NA_character_
        )
        next
      }

      current_status <- study_status(study_record)
      arms_text <- paste(study_record$arms[!vapply(study_record$arms, is_blank_cell, logical(1))], collapse = " | ")

      audit_rows[[length(audit_rows) + 1]] <- build_audit_row(
        sheet_name = sheet_name,
        block_index = block$block_index,
        row_index = row_index,
        study_id = study_id,
        mmc3_sheet_name = mmc3_sheet_name,
        placebo_columns = placebo_columns,
        performance_bias = study_record$performance_bias,
        detection_bias = study_record$detection_bias,
        has_pill_placebo_arm = current_status$pill_placebo,
        has_nonpharmacological_component = current_status$nonpharma,
        has_blinding_issue_value = current_status$blinding_issue,
        status = current_status$status,
        reason = current_status$reason,
        arms = arms_text
      )

      if (current_status$status != "reclassified") {
        next
      }

      for (column_name in placebo_columns) {
        openxlsx::writeData(
          workbook,
          sheet = sheet_name,
          x = PARTIAL_PLACEBO_CODE,
          startCol = unname(column_map[[column_name]]),
          startRow = row_index,
          colNames = FALSE,
          rowNames = FALSE
        )
      }

      flagged_rows[[length(flagged_rows) + 1]] <- data.frame(
        sheet_name = sheet_name,
        block_index = block$block_index,
        worksheet_row = row_index,
        study_id = study_id,
        matched_sheet = mmc3_sheet_name,
        replaced_treat_columns = paste(placebo_columns, collapse = " | "),
        performance_bias = study_record$performance_bias,
        detection_bias = study_record$detection_bias,
        arms = arms_text,
        original_treat_code = PLACEBO_CODE,
        replacement_treat_code = PARTIAL_PLACEBO_CODE,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    workbook = workbook,
    flagged = bind_rows(flagged_rows, empty_flagged_frame()),
    audit = bind_rows(audit_rows, empty_audit_frame())
  )
}

bind_rows <- function(rows, empty_frame) {
  if (length(rows) == 0) {
    return(empty_frame)
  }
  do.call(rbind, rows)
}

default_lookup_path <- function(base_dir, sheet_name) {
  preferred_lookup <- if (startsWith(sheet_name, "LS")) file.path(base_dir, "trt_to_class_ls.csv") else file.path(base_dir, "trt_to_class_ms.csv")
  if (startsWith(sheet_name, "LS") && !file.exists(preferred_lookup)) {
    stop(sprintf("Expected LS lookup file was not found: %s. Supply --lookup explicitly.", preferred_lookup), call. = FALSE)
  }
  preferred_lookup
}

update_lookup <- function(lookup_path, output_path) {
  assert_file_exists(lookup_path)
  lookup <- utils::read.csv2(lookup_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  if (!"trtcode" %in% names(lookup)) {
    names(lookup)[[1]] <- sub("^\\ufeff", "", names(lookup)[[1]])
  }
  required_columns <- c("trtcode", "trt", "classcode", "class")
  names(lookup) <- sub("^\\ufeff", "", names(lookup))
  missing_columns <- required_columns[!required_columns %in% names(lookup)]
  if (length(missing_columns) > 0) {
    stop(
      sprintf("Lookup file %s is missing required columns: %s", lookup_path, paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  has_partial_placebo <- any(as.character(lookup$trtcode) == as.character(PARTIAL_PLACEBO_CODE))
  if (!has_partial_placebo) {
    partial_row <- stats::setNames(as.list(rep("", ncol(lookup))), names(lookup))
    partial_row <- as.data.frame(partial_row, stringsAsFactors = FALSE)
    partial_row$trtcode <- as.character(PARTIAL_PLACEBO_CODE)
    partial_row$trt <- PARTIAL_PLACEBO_NAME
    partial_row$classcode <- as.character(PLACEBO_CLASS_CODE)
    partial_row$class <- PLACEBO_CLASS_NAME
    lookup <- rbind(lookup, partial_row)
  }

  utils::write.table(lookup, output_path, sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE, fileEncoding = "UTF-8")
}

write_reports <- function(results, output_dir, mmc5_path, lookup_path, sheet_name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  safe_sheet_name <- gsub("[/\\\\:*?\"<>|]", "_", sheet_name)

  workbook_output <- file.path(output_dir, sprintf("%s_partial_placebo.xlsx", tools::file_path_sans_ext(basename(mmc5_path))))
  flagged_output <- file.path(output_dir, sprintf("flagged_partial_placebo_%s.csv", gsub(" ", "_", safe_sheet_name)))
  review_output <- file.path(output_dir, sprintf("manual_review_partial_placebo_%s.csv", gsub(" ", "_", safe_sheet_name)))
  lookup_output <- file.path(output_dir, sprintf("%s_partial_placebo.csv", tools::file_path_sans_ext(basename(lookup_path))))

  update_lookup(lookup_path, lookup_output)
  utils::write.csv(results$flagged, flagged_output, row.names = FALSE, quote = TRUE, na = "")
  if (nrow(results$audit) == 0) {
    review_rows <- results$audit
  } else {
    review_rows <- results$audit[results$audit$status != "reclassified", , drop = FALSE]
  }
  utils::write.csv(review_rows, review_output, row.names = FALSE, quote = TRUE, na = "")
  openxlsx::saveWorkbook(results$workbook, workbook_output, overwrite = TRUE)

  list(
    workbook_output = workbook_output,
    flagged_output = flagged_output,
    review_output = review_output,
    lookup_output = lookup_output
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_args(args)
  sheet_name <- options$sheet
  mmc3_sheet_name <- options$mmc3_sheet %||% infer_mmc3_sheet(sheet_name)
  study_sheet <- load_study_sheet(options$mmc3, mmc3_sheet_name)
  lookup_path <- options$lookup %||% default_lookup_path(options$input_dir, sheet_name)
  results <- collect_reclassification_results(options$mmc5, study_sheet$records, sheet_name, study_sheet$sheet_name)
  outputs <- write_reports(results, options$output_dir, options$mmc5, lookup_path, sheet_name)

  message(sprintf("Resolved mmc5 sheet: %s", sheet_name))
  message(sprintf("Matched mmc3 sheet: %s", study_sheet$sheet_name))
  message(sprintf("Flagged studies: %d", nrow(results$flagged)))
  message(sprintf("Manual review rows: %d", nrow(results$audit[results$audit$status != "reclassified", , drop = FALSE])))
  if (nrow(results$flagged) > 0) {
    for (row_index in seq_len(nrow(results$flagged))) {
      message(sprintf(
        "- %s (block %s, row %s)",
        results$flagged$study_id[[row_index]],
        results$flagged$block_index[[row_index]],
        results$flagged$worksheet_row[[row_index]]
      ))
    }
  }
  message(sprintf("Recoded workbook: %s", outputs$workbook_output))
  message(sprintf("Flag report: %s", outputs$flagged_output))
  message(sprintf("Manual review report: %s", outputs$review_output))
  message(sprintf("Updated lookup: %s", outputs$lookup_output))

  invisible(outputs)
}

if (sys.nframe() == 0) {
  main()
}
