script_args <- commandArgs(trailingOnly = FALSE)
file_args <- grep("^--file=", script_args, value = TRUE)
test_file <- if (length(file_args) > 0) sub("^--file=", "", file_args[[1]]) else file.path(getwd(), "tests", "test_reclassify_partial_placebo.R")
repo_root <- dirname(dirname(normalizePath(test_file)))
source(file.path(repo_root, "reclassify_partial_placebo.R"))

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

withr_tempdir <- function(code) {
  temp_dir <- tempfile("partial-placebo-test-")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  force(code(temp_dir))
}

make_test_workbook <- function(path) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "MS SMD bias-adj")
  headers <- c("na[]", "t[,1]", "t[,2]", "t[,3]", "t[,4]", "t[,5]", "studyid")
  openxlsx::writeData(wb, "MS SMD bias-adj", t(headers), startRow = 1, startCol = 1, colNames = FALSE)
  rows <- list(
    c(2, 1, 42, NA, NA, NA, "S1"),
    c(2, 1, 42, NA, NA, NA, "S2"),
    c(3, 1, 1, 42, NA, NA, "S3")
  )
  for (idx in seq_along(rows)) {
    openxlsx::writeData(wb, "MS SMD bias-adj", t(rows[[idx]]), startRow = idx + 1, startCol = 1, colNames = FALSE)
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

make_manual_review_workbook <- function(path) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "MS SMD bias-adj")
  headers <- c("na[]", "t[,1]", "t[,2]", "t[,3]", "t[,4]", "t[,5]", "studyid")
  openxlsx::writeData(wb, "MS SMD bias-adj", t(headers), startRow = 1, startCol = 1, colNames = FALSE)
  rows <- list(
    c(2, 1, 42, NA, NA, NA, NA),
    c(2, 1, 42, NA, NA, NA, "S4")
  )
  for (idx in seq_along(rows)) {
    openxlsx::writeData(wb, "MS SMD bias-adj", t(rows[[idx]]), startRow = idx + 1, startCol = 1, colNames = FALSE)
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

make_no_audit_workbook <- function(path) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "MS SMD bias-adj")
  headers <- c("na[]", "t[,1]", "t[,2]", "studyid")
  openxlsx::writeData(wb, "MS SMD bias-adj", t(headers), startRow = 1, startCol = 1, colNames = FALSE)
  openxlsx::writeData(wb, "MS SMD bias-adj", t(c(2, 42, 43, "S1")), startRow = 2, startCol = 1, colNames = FALSE)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

make_test_lookup <- function(path) {
  writeLines(
    c(
      "trtcode;trt;classcode;class",
      "1;Pill placebo;1;Placebo",
      "42;Fluoxetine;42;Fluoxetine"
    ),
    con = path,
    useBytes = TRUE
  )
}

make_lookup_with_partial_placebo <- function(path) {
  writeLines(
    c(
      "trtcode;trt;classcode;class",
      "1;Pill placebo;1;Placebo",
      "100;Partial placebo;1;Placebo"
    ),
    con = path,
    useBytes = TRUE
  )
}

test_should_flag <- function() {
  flagged <- list(
    study_id = "flagged",
    arms = c("Pill placebo", "Bright light therapy", NA, NA, NA),
    performance_bias = "Unclear risk",
    detection_bias = "Low risk"
  )
  not_flagged <- list(
    study_id = "not_flagged",
    arms = c("Pill placebo", "Fluoxetine", NA, NA, NA),
    performance_bias = "Low risk",
    detection_bias = "Low risk"
  )

  assert_true(study_status(flagged)$status == "reclassified", "Expected flagged study to be reclassified.")
  assert_true(study_status(not_flagged)$status != "reclassified", "Expected pharmacological-only study not to be reclassified.")
}

test_infer_and_resolve_ls_sheet <- function() {
  inferred <- infer_mmc3_sheet("LS SMD bias-adj")
  resolved <- resolve_mmc3_sheet_name(c("LS depression -included studies"), inferred)

  assert_true(inferred == "LS depression-included studies", "Expected LS inference to use the standard sheet title.")
  assert_true(resolved == "LS depression -included studies", "Expected LS sheet resolution to fall back to the workbook's actual sheet name.")
}

test_load_study_sheet_resolves_ls_alias <- function() {
  withr_tempdir(function(temp_dir) {
    workbook_path <- file.path(temp_dir, "mmc3_included_studies.xlsx")
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "LS depression -included studies")
    sheet_data <- data.frame(
      "Study ID" = "LS1",
      "Arm 1 intervention" = "Pill placebo",
      "Arm 2 intervention" = "Bright light therapy",
      "Blinding of participants and personnel (performance bias)" = "High risk",
      "Blinding of outcome assessment (detection bias)" = "Low risk",
      check.names = FALSE
    )
    openxlsx::writeData(wb, "LS depression -included studies", sheet_data)
    openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)

    study_sheet <- load_study_sheet(workbook_path, "LS depression-included studies")

    assert_true(identical(study_sheet$sheet_name, "LS depression -included studies"), "Expected LS alias resolution to surface the actual workbook sheet name.")
    assert_true(identical(study_sheet$records$LS1$study_id, "LS1"), "Expected LS alias resolution to load the study sheet.")
  })
}

test_collect_reclassification_results <- function() {
  withr_tempdir(function(temp_dir) {
    mmc5_path <- file.path(temp_dir, "mmc5_fixed.xlsx")
    make_test_workbook(mmc5_path)
    mmc3_records <- list(
      S1 = list(study_id = "S1", arms = c("Pill placebo", "Bright light therapy", NA, NA, NA), performance_bias = "High risk", detection_bias = "Low risk"),
      S2 = list(study_id = "S2", arms = c("Pill placebo", "Fluoxetine", NA, NA, NA), performance_bias = "Low risk", detection_bias = "Low risk"),
      S3 = list(study_id = "S3", arms = c("Pill placebo", "Pill placebo", "Bright light therapy", NA, NA), performance_bias = "High risk", detection_bias = "Low risk")
    )

    results <- collect_reclassification_results(mmc5_path, mmc3_records, "MS SMD bias-adj", "MS depression-included studies")
    output_path <- file.path(temp_dir, "out.xlsx")
    openxlsx::saveWorkbook(results$workbook, output_path, overwrite = TRUE)
    recoded <- openxlsx::readWorkbook(output_path, sheet = "MS SMD bias-adj", colNames = FALSE)

    assert_true(recoded[[2]][[2]] == PARTIAL_PLACEBO_CODE, "Expected S1 placebo code to be rewritten.")
    assert_true(recoded[[2]][[3]] == PLACEBO_CODE, "Expected S2 placebo code to remain unchanged.")
    assert_true(recoded[[2]][[4]] == PARTIAL_PLACEBO_CODE, "Expected first S3 placebo code to be rewritten.")
    assert_true(recoded[[3]][[4]] == PARTIAL_PLACEBO_CODE, "Expected second S3 placebo code to be rewritten.")
    assert_true(identical(results$flagged$study_id, c("S1", "S3")), "Expected S1 and S3 to be flagged.")
    assert_true(identical(results$audit$study_id, c("S1", "S2", "S3")), "Expected audit output for all placebo-coded rows.")
    assert_true(identical(results$audit$reason, c("reclassified", "no_nonpharmacological_component_detected", "reclassified")), "Expected audit reasons to capture the non-reclassified S2 branch.")
  })
}

test_write_reports <- function() {
  withr_tempdir(function(temp_dir) {
    mmc5_path <- file.path(temp_dir, "mmc5_fixed.xlsx")
    lookup_path <- file.path(temp_dir, "trt_to_class_ms.csv")
    make_test_workbook(mmc5_path)
    make_test_lookup(lookup_path)
    mmc3_records <- list(
      S1 = list(study_id = "S1", arms = c("Pill placebo", "Bright light therapy", NA, NA, NA), performance_bias = "High risk", detection_bias = "Low risk"),
      S2 = list(study_id = "S2", arms = c("Pill placebo", "Fluoxetine", NA, NA, NA), performance_bias = "Low risk", detection_bias = "Low risk"),
      S3 = list(study_id = "S3", arms = c("Pill placebo", "Pill placebo", "Bright light therapy", NA, NA), performance_bias = "High risk", detection_bias = "Low risk")
    )

    results <- collect_reclassification_results(mmc5_path, mmc3_records, "MS SMD bias-adj", "MS depression-included studies")
    outputs <- write_reports(results, file.path(temp_dir, "outputs"), mmc5_path, lookup_path, "MS SMD bias-adj")
    lookup <- read.csv2(outputs$lookup_output, stringsAsFactors = FALSE)
    flagged <- read.csv(outputs$flagged_output, stringsAsFactors = FALSE)
    review <- read.csv(outputs$review_output, stringsAsFactors = FALSE)

    assert_true(file.exists(outputs$workbook_output), "Expected recoded workbook to be written.")
    assert_true(any(as.character(lookup$trtcode) == as.character(PARTIAL_PLACEBO_CODE)), "Expected updated lookup to include Partial placebo.")
    assert_true(identical(flagged$study_id, c("S1", "S3")), "Expected flagged report to contain only reclassified studies.")
    assert_true(identical(review$study_id, "S2"), "Expected manual review report to contain only non-reclassified placebo-coded studies.")
  })
}

test_collect_reclassification_manual_review_branches <- function() {
  withr_tempdir(function(temp_dir) {
    mmc5_path <- file.path(temp_dir, "mmc5_fixed.xlsx")
    lookup_path <- file.path(temp_dir, "trt_to_class_ms.csv")
    make_manual_review_workbook(mmc5_path)
    make_test_lookup(lookup_path)
    mmc3_records <- list()

    results <- collect_reclassification_results(mmc5_path, mmc3_records, "MS SMD bias-adj", "MS depression-included studies")
    outputs <- write_reports(results, file.path(temp_dir, "outputs"), mmc5_path, lookup_path, "MS SMD bias-adj")
    review <- read.csv(outputs$review_output, stringsAsFactors = FALSE)
    assert_true(nrow(results$flagged) == 0, "Expected no flagged rows for manual-review-only workbook.")
    assert_true(identical(results$audit$status, c("manual_review", "manual_review")), "Expected both audit rows to require manual review.")
    assert_true(identical(results$audit$reason, c("missing_study_id_in_mmc5", "study_not_found_in_mmc3")), "Expected audit reasons for missing study ID and missing mmc3 match.")
    assert_true(identical(review$reason, c("missing_study_id_in_mmc5", "study_not_found_in_mmc3")), "Expected manual review CSV to preserve manual-review audit rows.")
  })
}

test_write_reports_with_empty_audit <- function() {
  withr_tempdir(function(temp_dir) {
    mmc5_path <- file.path(temp_dir, "mmc5_fixed.xlsx")
    lookup_path <- file.path(temp_dir, "trt_to_class_ms.csv")
    make_no_audit_workbook(mmc5_path)
    make_test_lookup(lookup_path)
    empty_results <- list(
      workbook = openxlsx::loadWorkbook(mmc5_path),
      flagged = empty_flagged_frame(),
      audit = empty_audit_frame()
    )

    outputs <- write_reports(empty_results, file.path(temp_dir, "outputs"), mmc5_path, lookup_path, "MS SMD bias-adj")
    review <- read.csv(outputs$review_output, stringsAsFactors = FALSE)

    assert_true(nrow(review) == 0, "Expected empty manual review report when no placebo-coded rows were audited.")
  })
}

test_update_lookup_does_not_duplicate_partial_placebo <- function() {
  withr_tempdir(function(temp_dir) {
    lookup_path <- file.path(temp_dir, "trt_to_class_ms.csv")
    output_path <- file.path(temp_dir, "trt_to_class_ms_partial_placebo.csv")
    make_lookup_with_partial_placebo(lookup_path)

    update_lookup(lookup_path, output_path)
    lookup <- read.csv2(output_path, stringsAsFactors = FALSE)

    assert_true(sum(as.character(lookup$trtcode) == as.character(PARTIAL_PLACEBO_CODE)) == 1, "Expected Partial placebo row not to be duplicated.")
  })
}

test_write_reports_sanitizes_sheet_name <- function() {
  withr_tempdir(function(temp_dir) {
    mmc5_path <- file.path(temp_dir, "mmc5_fixed.xlsx")
    lookup_path <- file.path(temp_dir, "trt_to_class_ms.csv")
    make_no_audit_workbook(mmc5_path)
    make_test_lookup(lookup_path)
    empty_results <- list(
      workbook = openxlsx::loadWorkbook(mmc5_path),
      flagged = empty_flagged_frame(),
      audit = empty_audit_frame()
    )

    outputs <- write_reports(empty_results, file.path(temp_dir, "outputs"), mmc5_path, lookup_path, "MS/\\SMD bias-adj")

    assert_true(file.exists(outputs$flagged_output), "Expected flagged output to be written for sanitized sheet names.")
    assert_true(!grepl("[/\\\\]", basename(outputs$flagged_output)), "Expected output filename to remove path separators from sheet names.")
  })
}

test_should_flag()
test_infer_and_resolve_ls_sheet()
test_load_study_sheet_resolves_ls_alias()
test_collect_reclassification_results()
test_write_reports()
test_collect_reclassification_manual_review_branches()
test_write_reports_with_empty_audit()
test_update_lookup_does_not_duplicate_partial_placebo()
test_write_reports_sanitizes_sheet_name()
cat("All R placebo reclassification tests passed.\n")
