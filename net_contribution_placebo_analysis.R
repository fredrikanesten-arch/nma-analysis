# =============================================================================
# Net Contribution Analysis: Non-Placebo vs Placebo Comparisons
# =============================================================================
# This script analyses the net contribution data to summarise which direct
# comparisons drive each non-placebo vs Placebo network estimate, and
# produces one stacked bar plot per comparison.
#
# Input files (in the working directory):
#   - net_contribution_long.csv
#   - net_contribution_top10_by_estimate.csv
#
# Output:
#   - Console: detailed text summary for each comparison
#   - PDF: one page per comparison showing the contribution breakdown
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

# ── 0. Configuration ──────────────────────────────────────────────────────────

TOP_N          <- 10      # Number of top contributors to show individually
OUTPUT_PDF     <- "net_contribution_placebo_plots.pdf"
PLACEBO_LABEL  <- "Placebo"

# ── 1. Load data ──────────────────────────────────────────────────────────────

nc_long <- read.csv("net_contribution_long.csv", stringsAsFactors = FALSE)

# Verify expected columns
stopifnot(all(c("direct_comparison", "network_estimate", "contribution") %in%
                names(nc_long)))

# ── 2. Filter: network estimates where exactly one side is "Placebo" ──────────

nc_long <- nc_long %>%
  mutate(
    treat1 = str_trim(sub(":.*", "", network_estimate)),
    treat2 = str_trim(sub(".*:", "", network_estimate))
  )

placebo_estimates <- nc_long %>%
  filter(treat1 == PLACEBO_LABEL | treat2 == PLACEBO_LABEL) %>%
  # Standardise so Placebo is always treat2 (non-placebo : Placebo)
  mutate(
    non_placebo = if_else(treat1 == PLACEBO_LABEL, treat2, treat1),
    network_estimate_std = paste0(non_placebo, " : ", PLACEBO_LABEL)
  )

comparisons <- sort(unique(placebo_estimates$non_placebo))

cat(sprintf("Found %d non-placebo vs Placebo network estimates.\n\n",
            length(comparisons)))

# ── 3. Helper: build contribution table for one estimate ───────────────────────

build_contrib_table <- function(non_placebo_name) {
  df <- placebo_estimates %>%
    filter(non_placebo == non_placebo_name) %>%
    arrange(desc(contribution))

  total  <- sum(df$contribution)
  df <- df %>%
    mutate(
      pct         = contribution / total * 100,
      cum_pct     = cumsum(pct),
      # Is this the "direct" comparison for the network estimate?
      is_direct   = (direct_comparison == network_estimate |
                     direct_comparison == paste0(non_placebo_name, ":", PLACEBO_LABEL) |
                     direct_comparison == paste0(PLACEBO_LABEL, ":", non_placebo_name)),
      rank_in_est = row_number()
    )
  list(df = df, total = total)
}

# ── 4. Text summary ───────────────────────────────────────────────────────────

cat("=======================================================================\n")
cat("DETAILED SUMMARY: Non-Placebo vs Placebo Net Contributions\n")
cat("=======================================================================\n\n")

for (comp in comparisons) {
  res    <- build_contrib_table(comp)
  df     <- res$df
  total  <- res$total
  n_rows <- nrow(df)

  # Identify the direct-evidence row (if any)
  direct_row <- df %>% filter(is_direct)
  has_direct <- nrow(direct_row) > 0

  cat(sprintf("-----------------------------------------------------------------------\n"))
  cat(sprintf("Comparison: %s vs %s\n", comp, PLACEBO_LABEL))
  cat(sprintf("  Total absolute net contribution: %.4f\n", total))
  cat(sprintf("  Number of contributing direct comparisons: %d\n", n_rows))

  if (has_direct) {
    cat(sprintf("  Direct evidence (%s:%s): %.4f (%.1f%% of total)\n",
                comp, PLACEBO_LABEL,
                direct_row$contribution[1], direct_row$pct[1]))
  } else {
    cat("  Direct evidence: NONE (estimate is entirely indirect)\n")
  }
  cat("\n")

  cat(sprintf("  Top %d contributing direct comparisons:\n", min(TOP_N, n_rows)))
  cat(sprintf("  %-55s %8s %7s\n", "Direct comparison", "Contrib.", "%"))
  cat(sprintf("  %-55s %8s %7s\n", strrep("-", 55), strrep("-", 8), strrep("-", 7)))

  top_df <- head(df, TOP_N)
  for (i in seq_len(nrow(top_df))) {
    flag <- if (top_df$is_direct[i]) " <-- direct" else ""
    cat(sprintf("  %-55s %8.4f %6.1f%%%s\n",
                top_df$direct_comparison[i],
                top_df$contribution[i],
                top_df$pct[i],
                flag))
  }

  # Cumulative coverage of top N
  cum_top <- sum(top_df$pct)
  cat(sprintf("\n  Top %d comparisons account for %.1f%% of total contribution.\n",
              min(TOP_N, n_rows), cum_top))

  # How many comparisons needed to reach 80 % / 90 %?
  for (threshold in c(80, 90)) {
    n_needed <- which(df$cum_pct >= threshold)[1]
    cat(sprintf("  %d%% of contribution reached after %d comparisons.\n",
                threshold, n_needed))
  }
  cat("\n")
}

# ── 5. Plots ──────────────────────────────────────────────────────────────────

# Colour palette: direct comparison gets a distinct colour, others use
# a qualitative palette, "Other" is always grey.

make_plot_data <- function(non_placebo_name) {
  res   <- build_contrib_table(non_placebo_name)
  df    <- res$df
  total <- res$total

  # Assign group labels: top-N shown individually, rest → "Other"
  df <- df %>%
    mutate(
      group = if_else(rank_in_est <= TOP_N,
                      direct_comparison,
                      "Other")
    )

  grouped <- df %>%
    group_by(group) %>%
    summarise(
      pct       = sum(pct),
      is_direct = any(is_direct),
      min_rank  = min(rank_in_est),
      .groups   = "drop"
    ) %>%
    arrange(min_rank) %>%
    # Put "Other" at the end
    mutate(order = if_else(group == "Other", Inf, as.numeric(min_rank))) %>%
    arrange(order)

  # Wrap long labels
  grouped <- grouped %>%
    mutate(group_label = str_wrap(group, width = 40))

  # Factor levels: top-N in descending contribution, Other last
  grouped$group_label <- factor(grouped$group_label,
                                levels = rev(grouped$group_label))

  list(data = grouped, total = total)
}

# Build a colour palette for up to TOP_N + 1 groups
# Direct comparison: dark red; top-N: viridis-like; Other: light grey
make_colours <- function(grouped_df) {
  n      <- nrow(grouped_df)
  labels <- as.character(levels(grouped_df$group_label))

  # Base palette (excluding Other and direct)
  base_pal <- hcl.colors(max(n, 3), palette = "Dark 2")

  cols <- setNames(base_pal[seq_len(n)], labels)

  # Override "Other"
  other_label <- as.character(grouped_df$group_label[grouped_df$group == "Other"])
  if (length(other_label) > 0)
    cols[other_label] <- "grey80"

  # Override direct comparison
  direct_label <- as.character(
    grouped_df$group_label[grouped_df$is_direct & grouped_df$group != "Other"])
  if (length(direct_label) > 0)
    cols[direct_label] <- "#c0392b"   # strong red

  cols
}

cat("Generating plots …\n")

pdf(OUTPUT_PDF, width = 11, height = 7)

for (comp in comparisons) {
  pd  <- make_plot_data(comp)
  gdf <- pd$data

  colours <- make_colours(gdf)

  # One horizontal stacked bar (proportion)
  p <- ggplot(gdf, aes(x = 1, y = pct, fill = group_label)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(
      aes(label = if_else(pct >= 1.5,
                          paste0(round(pct, 1), "%"),
                          "")),
      position = position_stack(vjust = 0.5),
      size = 3.2, colour = "white", fontface = "bold"
    ) +
    coord_flip() +
    scale_fill_manual(values = colours, name = "Contributing\ndirect comparison") +
    scale_y_continuous(labels = percent_format(scale = 1),
                       expand = expansion(mult = c(0, 0.01))) +
    labs(
      title    = sprintf("Net contribution to: %s vs %s", comp, PLACEBO_LABEL),
      subtitle = sprintf(
        "Top %d direct comparisons shown individually (direct evidence in red, if present).\nTotal absolute contribution = %.4f",
        TOP_N, pd$total),
      x = NULL,
      y = "% of total net contribution"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position  = "right",
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.45, "cm"),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(size = 9, colour = "grey40")
    ) +
    guides(fill = guide_legend(reverse = TRUE, ncol = 1))

  print(p)
}

dev.off()

cat(sprintf("\nDone. Plots written to '%s'.\n", OUTPUT_PDF))
