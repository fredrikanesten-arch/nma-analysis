if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install with: install.packages('readxl')")
}

OUTPUT_COLUMNS <- c(
  "studyid",
  "na",
  "arm",
  "treatment",
  "n",
  "mean_change",
  "sd_change",
  "source",
  "y_baseline",
  "sd_baseline",
  "y_followup",
  "sd_followup",
  "r_used",
  "responders",
  "p_response",
  "q",
  "cutoff"
)

# -----------------------------
# Strict assumptions aligned to equations
# -----------------------------
ASSUMED_R_BF <- 0.20
P_CLIP <- 1e-6
INCLUDE_RESPONDERS <- TRUE

as_num <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "NA")) return(NA_real_)
  as.numeric(x)
}

as_int <- function(x) {
  v <- as_num(x)
  if (is.na(v)) return(NA_integer_)
  as.integer(v)
}

cell_value <- function(df, r, c) {
  if (r < 1 || c < 1 || r > nrow(df) || c > ncol(df)) return(NA)
  df[[c]][r]
}

find_header_row <- function(df, marker) {
  for (r in seq_len(nrow(df))) {
    if (identical(cell_value(df, r, 1), "na[]") && identical(cell_value(df, r, 7), marker)) {
      return(r)
    }
  }
  stop(sprintf("Could not find header row with marker %s", marker))
}

empty_output <- function() {
  data.frame(
    studyid = character(0),
    na = integer(0),
    arm = integer(0),
    treatment = integer(0),
    n = integer(0),
    mean_change = numeric(0),
    sd_change = numeric(0),
    source = character(0),
    y_baseline = numeric(0),
    sd_baseline = numeric(0),
    y_followup = numeric(0),
    sd_followup = numeric(0),
    r_used = numeric(0),
    responders = integer(0),
    p_response = numeric(0),
    q = numeric(0),
    cutoff = numeric(0),
    stringsAsFactors = FALSE
  )
}

append_row <- function(out, row) {
  out[nrow(out) + 1, OUTPUT_COLUMNS] <- row[OUTPUT_COLUMNS]
  out
}

pooled_sd_from_arms <- function(n, sd) {
  dfree <- sum(n) - length(n)
  if (dfree <= 0) stop("Invalid pooled SD degrees of freedom")
  pooled_var <- sum((n - 1) * (sd^2)) / dfree
  if (!is.finite(pooled_var) || pooled_var <= 0) stop("Invalid pooled SD variance")
  sqrt(pooled_var)
}

hedges_j <- function(dfree) {
  if (is.na(dfree) || dfree <= 1) return(1)
  1 - (3 / (4 * dfree - 1))
}

clip_p <- function(p) {
  pmin(pmax(p, P_CLIP), 1 - P_CLIP)
}

parse_cfb_block <- function(df, start_row, end_row) {
  out <- empty_output()

  for (r in start_row:end_row) {
    studyid <- cell_value(df, r, 26)
    na <- as_int(cell_value(df, r, 1))

    if (is.na(na) || is.na(studyid) || identical(studyid, "#")) next

    arms <- list()
    for (arm in 1:5) {
      trt <- as_int(cell_value(df, r, 1 + arm))
      n <- as_int(cell_value(df, r, 16 + arm))
      y <- as_num(cell_value(df, r, 6 + arm))
      sd <- as_num(cell_value(df, r, 11 + arm))

      if (is.na(trt)) next
      if (any(is.na(c(n, y, sd)))) stop(sprintf("Missing CFB values in row %d, arm %d", r, arm))

      arms[[length(arms) + 1]] <- list(arm = arm, trt = trt, n = n, y = y, sd = sd)
    }

    arm_n <- vapply(arms, function(a) a$n, numeric(1))
    arm_sd <- vapply(arms, function(a) a$sd, numeric(1))
    pooled_sd <- pooled_sd_from_arms(n = arm_n, sd = arm_sd)
    j <- hedges_j(sum(arm_n) - length(arm_n))

    for (a in arms) {
      row <- list(
        studyid = as.character(studyid),
        na = na,
        arm = a$arm,
        treatment = a$trt,
        n = a$n,
        mean_change = (a$y / pooled_sd) * j,
        sd_change = a$sd / pooled_sd,
        source = "change_from_baseline",
        y_baseline = NA_real_,
        sd_baseline = NA_real_,
        y_followup = NA_real_,
        sd_followup = NA_real_,
        r_used = NA_real_,
        responders = NA_integer_,
        p_response = NA_real_,
        q = NA_real_,
        cutoff = NA_real_
      )
      out <- append_row(out, row)
    }
  }

  out
}

parse_bf_block <- function(df, start_row, end_row, r_assumed = ASSUMED_R_BF) {
  out <- empty_output()

  for (r in start_row:end_row) {
    studyid <- cell_value(df, r, 37)
    na <- as_int(cell_value(df, r, 1))

    if (is.na(na) || is.na(studyid) || identical(studyid, "#")) next

    arms <- list()
    for (arm in 1:5) {
      trt <- as_int(cell_value(df, r, 1 + arm))
      n <- as_int(cell_value(df, r, 26 + arm))
      yb <- as_num(cell_value(df, r, 6 + arm))
      sdb <- as_num(cell_value(df, r, 11 + arm))
      yf <- as_num(cell_value(df, r, 16 + arm))
      sdf <- as_num(cell_value(df, r, 21 + arm))

      if (is.na(trt)) next
      if (any(is.na(c(n, yb, sdb, yf, sdf)))) {
        stop(sprintf("Missing baseline/follow-up values in row %d, arm %d", r, arm))
      }

      mean_change_raw <- yf - yb
      var_change <- sdb^2 + sdf^2 - 2 * r_assumed * sdb * sdf
      if (!is.finite(var_change) || var_change <= 0) {
        stop(sprintf("Non-positive change variance in row %d, arm %d", r, arm))
      }
      sd_change_raw <- sqrt(var_change)

      arms[[length(arms) + 1]] <- list(
        arm = arm, trt = trt, n = n, yb = yb, sdb = sdb, yf = yf, sdf = sdf,
        mean_change_raw = mean_change_raw, sd_change_raw = sd_change_raw
      )
    }

    arm_n <- vapply(arms, function(a) a$n, numeric(1))
    arm_sd_for_pool <- vapply(arms, function(a) max(a$sdb, a$sdf), numeric(1))
    pooled_sd <- pooled_sd_from_arms(n = arm_n, sd = arm_sd_for_pool)
    j <- hedges_j(sum(arm_n) - length(arm_n))

    for (a in arms) {
      row <- list(
        studyid = as.character(studyid),
        na = na,
        arm = a$arm,
        treatment = a$trt,
        n = a$n,
        mean_change = (a$mean_change_raw / pooled_sd) * j,
        sd_change = a$sd_change_raw / pooled_sd,
        source = "baseline_followup",
        y_baseline = a$yb,
        sd_baseline = a$sdb,
        y_followup = a$yf,
        sd_followup = a$sdf,
        r_used = r_assumed,
        responders = NA_integer_,
        p_response = NA_real_,
        q = NA_real_,
        cutoff = NA_real_
      )
      out <- append_row(out, row)
    }
  }

  out
}

parse_response_block <- function(df, start_row, end_row) {
  out <- empty_output()

  for (r in start_row:end_row) {
    studyid <- cell_value(df, r, 33)
    na <- as_int(cell_value(df, r, 1))
    q <- as_num(cell_value(df, r, 27))

    if (is.na(na) || is.na(studyid) || identical(studyid, "#")) next
    if (is.na(q)) stop(sprintf("Missing q in response row %d", r))

    sd_multiplier <- sqrt(1 + (1 - q) * (1 - q - 2 * ASSUMED_R_BF))
    if (!is.finite(sd_multiplier) || sd_multiplier <= 0) {
      stop(sprintf("Invalid responder SD multiplier in row %d", r))
    }

    arms <- list()
    for (arm in 1:5) {
      trt <- as_int(cell_value(df, r, 1 + arm))
      responders <- as_int(cell_value(df, r, 6 + arm))
      n <- as_int(cell_value(df, r, 11 + arm))
      ybr <- as_num(cell_value(df, r, 16 + arm))
      sdbr <- as_num(cell_value(df, r, 21 + arm))

      if (is.na(trt)) next
      if (any(is.na(c(responders, n, ybr, sdbr)))) {
        stop(sprintf("Missing responder values in row %d, arm %d", r, arm))
      }
      if (responders < 0 || responders > n) {
        stop(sprintf("Invalid responders count in row %d, arm %d", r, arm))
      }

      p_response <- responders / n
      z <- qnorm(clip_p(p_response))

      sd_rule <- sdbr
      cutoff <- -(q * ybr)
      mean_change_raw <- -(q * ybr + z * sd_rule * sd_multiplier)
      sd_change_raw <- sd_rule * sd_multiplier

      arms[[length(arms) + 1]] <- list(
        arm = arm, trt = trt, n = n, responders = responders, p_response = p_response,
        ybr = ybr, sdbr = sdbr, cutoff = cutoff,
        sd_rule = sd_rule, mean_change_raw = mean_change_raw, sd_change_raw = sd_change_raw
      )
    }

    arm_n <- vapply(arms, function(a) a$n, numeric(1))
    arm_sd_for_pool <- vapply(arms, function(a) a$sd_change_raw, numeric(1))
    pooled_sd <- pooled_sd_from_arms(n = arm_n, sd = arm_sd_for_pool)
    j <- hedges_j(sum(arm_n) - length(arm_n))

    for (a in arms) {
      row <- list(
        studyid = as.character(studyid),
        na = na,
        arm = a$arm,
        treatment = a$trt,
        n = a$n,
        mean_change = (a$mean_change_raw / pooled_sd) * j,
        sd_change = a$sd_change_raw / pooled_sd,
        source = "binary_response_converted",
        y_baseline = a$ybr,
        sd_baseline = a$sdbr,
        y_followup = NA_real_,
        sd_followup = NA_real_,
        r_used = NA_real_,
        responders = a$responders,
        p_response = a$p_response,
        q = q,
        cutoff = a$cutoff
      )
      out <- append_row(out, row)
    }
  }

  out
}

main <- function() {
  base_dir <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli"
  out_dir <- file.path(base_dir, "binfixed_class_ms")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[[1]] %in% c("-h", "--help")) {
    cat("Usage: Rscript build_mmc5_ms_smd_bias_adj_dataset.R [input_xlsx] [output_csv] [sheet]\n")
    cat("Defaults: input_xlsx='<base_dir>/article_supplements/mmc5.xlsx', output_csv='<base_dir>/binfixed_class_ms/combined_long_mean_change_dataset_ms_smd_bias_adj.csv', sheet='MS SMD bias-adj'\n")
    cat("Strict assumptions: rho=0.2, response conversion uses q, strict responder SD rule, pooled SD standardization, Hedges correction\n")
    return(invisible(NULL))
  }

  default_input_xlsx <- file.path(base_dir, "article_supplements", "mmc5.xlsx")
  default_output_csv <- file.path(out_dir, "combined_long_mean_change_dataset_ms_smd_bias_adj.csv")

  input_xlsx <- if (length(args) >= 1) args[[1]] else default_input_xlsx
  output_csv <- if (length(args) >= 2) args[[2]] else default_output_csv
  sheet <- if (length(args) >= 3) args[[3]] else "MS SMD bias-adj"

  if (!file.exists(input_xlsx)) {
    stop(sprintf("Input workbook not found: %s", input_xlsx))
  }

  raw <- readxl::read_excel(
    input_xlsx,
    sheet = sheet,
    col_names = FALSE,
    .name_repair = "minimal"
  )

  cfb_header <- find_header_row(raw, "yCFB[,1]")
  bf_header <- find_header_row(raw, "yB[,1]")
  resp_header <- find_header_row(raw, "r[,1]")

  cfb <- parse_cfb_block(raw, cfb_header + 1, bf_header - 1)
  bf <- parse_bf_block(raw, bf_header + 1, resp_header - 1, r_assumed = ASSUMED_R_BF)
  resp <- parse_response_block(raw, resp_header + 1, nrow(raw))

  out <- if (INCLUDE_RESPONDERS) rbind(cfb, bf, resp) else rbind(cfb, bf)
  out <- out[, OUTPUT_COLUMNS]

  utils::write.csv(out, file = output_csv, row.names = FALSE, na = "NA", quote = TRUE)
  cat(sprintf("Wrote %d rows to %s\n", nrow(out), output_csv))
}

main()
