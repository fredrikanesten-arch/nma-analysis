script_args <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
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

test_should_flag()
test_collect_reclassification_results()
test_write_reports()
cat("All R placebo reclassification tests passed.\n")
