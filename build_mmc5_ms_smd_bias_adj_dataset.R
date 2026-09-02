#!/usr/bin/env Rscript

# Build an arm-level SMD dataset from the 'MS SMD bias-adj' sheet in mmc5.xlsx.
# The output schema matches combined_long_mean_change_dataset_ms.csv columns.

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

pooled_sd_from_arms <- function(n, sd) {
  dfree <- sum(n) - length(n)
  if (dfree <= 0) stop("Invalid pooled SD degrees of freedom")
  pooled_var <- sum((n - 1) * (sd^2)) / dfree
  if (pooled_var <= 0) stop("Invalid pooled SD variance")
  sqrt(pooled_var)
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
      if (is.na(n) || is.na(y) || is.na(sd)) {
        stop(sprintf("Missing CFB values in row %d, arm %d", r, arm))
      }

      arms[[length(arms) + 1]] <- list(arm = arm, trt = trt, n = n, y = y, sd = sd)
    }

    pooled_sd <- pooled_sd_from_arms(
      n = vapply(arms, function(a) a$n, numeric(1)),
      sd = vapply(arms, function(a) a$sd, numeric(1))
    )

    for (a in arms) {
      row <- list(
        studyid = as.character(studyid),
        na = na,
        arm = a$arm,
        treatment = a$trt,
        n = a$n,
        mean_change = a$y / pooled_sd,
        sd_change = a$sd / pooled_sd,
        source = "smd",
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

parse_bf_block <- function(df, start_row, end_row, rho) {
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
      var_change <- (sdf^2) + (sdb^2) - (2 * rho * sdf * sdb)
      if (var_change <= 0) {
        stop(sprintf("Non-positive change variance in row %d, arm %d", r, arm))
      }
      sd_change_raw <- sqrt(var_change)

      arms[[length(arms) + 1]] <- list(
        arm = arm, trt = trt, n = n, yb = yb, sdb = sdb,
        yf = yf, sdf = sdf, mean_change_raw = mean_change_raw,
        sd_change_raw = sd_change_raw
      )
    }

    pooled_sd <- pooled_sd_from_arms(
      n = vapply(arms, function(a) a$n, numeric(1)),
      sd = vapply(arms, function(a) a$sdb, numeric(1))
    )

    for (a in arms) {
      row <- list(
        studyid = as.character(studyid),
        na = na,
        arm = a$arm,
        treatment = a$trt,
        n = a$n,
        mean_change = a$mean_change_raw / pooled_sd,
        sd_change = a$sd_change_raw / pooled_sd,
        source = "baseline_followup",
        y_baseline = a$yb,
        sd_baseline = a$sdb,
        y_followup = a$yf,
        sd_followup = a$sdf,
        r_used = rho,
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

parse_response_block <- function(df, start_row, end_row, rho, p_clip) {
  out <- empty_output()

  for (r in start_row:end_row) {
    studyid <- cell_value(df, r, 33)
    na <- as_int(cell_value(df, r, 1))
    q <- as_num(cell_value(df, r, 27))

    if (is.na(na) || is.na(studyid) || identical(studyid, "#")) next
    if (is.na(q)) stop(sprintf("Missing q in response row %d", r))

    adj <- sqrt(1 + (1 - q) * (1 - q - 2 * rho))
    if (adj <= 0) stop(sprintf("Invalid response adjustment in row %d", r))

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

      p <- responders / n
      p_for_z <- min(max(p, p_clip), 1 - p_clip)
      z <- qnorm(p_for_z)

      cutoff <- -(q * ybr)
      mean_change_raw <- -(q * ybr + (z * sdbr * adj))
      sd_change_raw <- sdbr * adj

      arms[[length(arms) + 1]] <- list(
        arm = arm, trt = trt, n = n, responders = responders, p = p,
        ybr = ybr, sdbr = sdbr, cutoff = cutoff,
        mean_change_raw = mean_change_raw, sd_change_raw = sd_change_raw
      )
    }

    pooled_sd <- pooled_sd_from_arms(
      n = vapply(arms, function(a) a$n, numeric(1)),
      sd = vapply(arms, function(a) a$sdbr, numeric(1))
    )

    for (a in arms) {
      row <- list(
        studyid = as.character(studyid),
        na = na,
        arm = a$arm,
        treatment = a$trt,
        n = a$n,
        mean_change = a$mean_change_raw / pooled_sd,
        sd_change = a$sd_change_raw / pooled_sd,
        source = "responders",
        y_baseline = a$ybr,
        sd_baseline = a$sdbr,
        y_followup = NA_real_,
        sd_followup = NA_real_,
        r_used = NA_real_,
        responders = a$responders,
        p_response = a$p,
        q = q,
        cutoff = a$cutoff
      )
      out <- append_row(out, row)
    }
  }

  out
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[[1]] %in% c("-h", "--help")) {
    cat("Usage: Rscript build_mmc5_ms_smd_bias_adj_dataset.R [input_xlsx] [output_csv] [sheet] [rho] [p_clip]\n")
    cat("Defaults: input_xlsx='mmc5.xlsx', output_csv='combined_long_mean_change_dataset_ms_smd_bias_adj.csv', sheet='MS SMD bias-adj', rho=0.5, p_clip=1e-6\n")
    return(invisible(NULL))
  }

  input_xlsx <- if (length(args) >= 1) args[[1]] else "mmc5.xlsx"
  output_csv <- if (length(args) >= 2) args[[2]] else "combined_long_mean_change_dataset_ms_smd_bias_adj.csv"
  sheet <- if (length(args) >= 3) args[[3]] else "MS SMD bias-adj"
  rho <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.5
  p_clip <- if (length(args) >= 5) as.numeric(args[[5]]) else 1e-6

  if (!file.exists(input_xlsx)) {
    stop(sprintf(
      "Input workbook not found: %s\nEither pass the path explicitly as the first argument or run from a directory containing mmc5.xlsx.",
      input_xlsx
    ))
  }
  if (is.na(rho)) stop("rho must be numeric")
  if (is.na(p_clip) || p_clip <= 0 || p_clip >= 0.5) stop("p_clip must be in (0, 0.5)")

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
  bf <- parse_bf_block(raw, bf_header + 1, resp_header - 1, rho = rho)
  resp <- parse_response_block(raw, resp_header + 1, nrow(raw), rho = rho, p_clip = p_clip)

  out <- rbind(cfb, bf, resp)
  out <- out[, OUTPUT_COLUMNS]

  utils::write.csv(out, file = output_csv, row.names = FALSE, na = "NA", quote = TRUE)
  cat(sprintf("Wrote %d rows to %s\n", nrow(out), output_csv))
}

main()
