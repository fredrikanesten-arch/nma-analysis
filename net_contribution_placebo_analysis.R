# =============================================================================
# Net Contribution Analysis: ALL Non-Placebo vs Placebo Comparisons
# =============================================================================
# Analyses the net contribution data for all 49 treatments vs Placebo.
# For treatments with a direct network estimate (X:Placebo present in
# net_contribution_long.csv), all contributing comparisons are used.
# For indirect-only treatments, the pre-computed top-10 contributions from
# indirect_only_contribmatrix_top10.csv are used.
#
# Input files (in the working directory):
#   - net_contribution_long.csv
#   - indirect_only_contribmatrix_top10.csv
#   - bayes_relative_effects_vs_placebo_main_labeled.csv
#
# Output:
#   - Console: detailed text summary for every comparison (top 5 shown)
#   - PDF:     one stacked bar chart per comparison (top 5 individually)
# =============================================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

# ── 0. Configuration ──────────────────────────────────────────────────────────

TOP_N         <- 5        # contributors to show individually in plots/summary
OUTPUT_PDF    <- "net_contribution_placebo_plots.pdf"
PLACEBO_LABEL <- "Placebo"

# ── 1. Load data ──────────────────────────────────────────────────────────────

nc_long   <- read.csv("net_contribution_long.csv",
                      stringsAsFactors = FALSE)
nc_ind    <- read.csv("indirect_only_contribmatrix_top10.csv",
                      stringsAsFactors = FALSE)
bayes_eff <- read.csv("bayes_relative_effects_vs_placebo_main_labeled.csv",
                      sep = ";", stringsAsFactors = FALSE)

stopifnot(all(c("direct_comparison", "network_estimate", "contribution") %in%
                names(nc_long)))
stopifnot(all(c("network_estimate", "direct_comparison_raw",
                "contribution", "contribution_pct") %in% names(nc_ind)))
stopifnot(all(c("class", "mean") %in% names(bayes_eff)))

# ── 2. Extract all treatments vs Placebo from the DIRECT estimates ────────────

split_ne <- function(ne) {
  # network_estimate is "TreatA:TreatB"; split on first colon only
  list(
    t1 = str_trim(sub(":.*", "", ne)),
    t2 = str_trim(sub("^[^:]+:", "", ne))
  )
}

nc_long <- nc_long %>%
  mutate(
    treat1 = str_trim(sub(":.*", "", network_estimate)),
    treat2 = str_trim(sub("^[^:]+:", "", network_estimate))
  )

direct_placebo <- nc_long %>%
  filter(treat1 == PLACEBO_LABEL | treat2 == PLACEBO_LABEL) %>%
  mutate(
    non_placebo = if_else(treat1 == PLACEBO_LABEL, treat2, treat1),
    # flag whether this row IS the direct X:Placebo comparison
    is_own_direct = (
      direct_comparison == paste0(non_placebo, ":", PLACEBO_LABEL) |
      direct_comparison == paste0(PLACEBO_LABEL, ":", non_placebo)
    )
  )

direct_treatments <- unique(direct_placebo$non_placebo)

# ── 3. Extract indirect-only treatments from the pre-computed top-10 ──────────

nc_ind <- nc_ind %>%
  mutate(
    treat1 = str_trim(sub(":.*", "", network_estimate)),
    treat2 = str_trim(sub("^[^:]+:", "", network_estimate))
  )

indirect_placebo <- nc_ind %>%
  filter(treat1 == PLACEBO_LABEL | treat2 == PLACEBO_LABEL) %>%
  mutate(
    non_placebo = if_else(treat1 == PLACEBO_LABEL, treat2, treat1),
    direct_comparison = direct_comparison_raw,
    is_own_direct = FALSE       # by definition all indirect
  ) %>%
  select(non_placebo, network_estimate, direct_comparison,
         contribution, contribution_pct, is_own_direct)

indirect_treatments <- unique(indirect_placebo$non_placebo)

# ── 4. Build ordered list of all 49 treatments (by NMA effect size) ───────────

treat_order <- bayes_eff %>%
  arrange(desc(mean)) %>%
  pull(class)

all_treatments <- treat_order[treat_order %in%
                               c(direct_treatments, indirect_treatments)]

cat(sprintf(
  "Treatments with direct X:Placebo evidence: %d\n",
  length(direct_treatments)
))
cat(sprintf(
  "Treatments with indirect-only evidence:    %d\n",
  length(indirect_treatments)
))
cat(sprintf("Total comparisons to analyse:              %d\n\n",
            length(all_treatments)))

# ── 5. Helper: build contribution table for one treatment ─────────────────────

build_contrib_table <- function(treat_name) {
  is_indirect <- treat_name %in% indirect_treatments &&
                 !(treat_name %in% direct_treatments)

  if (!is_indirect) {
    # Full rows available
    df <- direct_placebo %>%
      filter(non_placebo == treat_name) %>%
      arrange(desc(contribution)) %>%
      mutate(
        total       = sum(contribution),
        pct         = contribution / total * 100,
        cum_pct     = cumsum(pct),
        rank_in_est = row_number()
      )
    list(
      df              = df,
      total           = df$total[1],
      evidence_type   = "direct + indirect",
      top10_only      = FALSE
    )
  } else {
    # Only top-10 rows available; pct is already computed
    df <- indirect_placebo %>%
      filter(non_placebo == treat_name) %>%
      arrange(desc(contribution)) %>%
      mutate(
        pct         = contribution_pct,
        rank_in_est = row_number()
      )
    list(
      df              = df,
      total           = NA_real_,   # total not available for indirect-only
      evidence_type   = "indirect only",
      top10_only      = TRUE
    )
  }
}

# ── 6. Text summary ───────────────────────────────────────────────────────────

cat(strrep("=", 72), "\n")
cat("SUMMARY: Top 5 contributors to each non-Placebo vs Placebo estimate\n")
cat("(ordered by NMA posterior mean effect size, largest first)\n")
cat(strrep("=", 72), "\n\n")

for (treat in all_treatments) {
  res   <- build_contrib_table(treat)
  df    <- res$df
  nma_mean <- bayes_eff$mean[bayes_eff$class == treat]

  cat(strrep("-", 72), "\n")
  cat(sprintf("%-52s NMA mean = %.3f\n", paste0(treat, " vs Placebo"),
              nma_mean))
  cat(sprintf("Evidence type: %s\n", res$evidence_type))

  if (!res$top10_only) {
    cat(sprintf("Total abs. net contribution: %.4f\n", res$total))
    has_direct <- any(df$is_own_direct)
    cat(sprintf("Own direct trial evidence:   %s\n",
                if (has_direct) {
                  direct_pct <- df$pct[df$is_own_direct][1]
                  sprintf("YES (%.1f%% of total)", direct_pct)
                } else "NO"))
  } else {
    cat(sprintf("(Only top-10 pre-computed; top 10 cover %.1f%% of total)\n",
                sum(df$pct)))
  }

  cat(sprintf("\nTop %d contributing direct comparisons:\n", min(TOP_N, nrow(df))))
  cat(sprintf("  %-55s %6s\n", "Direct comparison", "%"))
  cat(sprintf("  %-55s %6s\n", strrep("-", 55), strrep("-", 6)))
  for (i in seq_len(min(TOP_N, nrow(df)))) {
    flag <- if (df$is_own_direct[i]) " <-- DIRECT" else ""
    cat(sprintf("  %-55s %5.1f%%%s\n",
                df$direct_comparison[i], df$pct[i], flag))
  }
  cat("\n")
}

# ── 7. Plots ──────────────────────────────────────────────────────────────────

make_plot_data <- function(treat_name) {
  res <- build_contrib_table(treat_name)
  df  <- res$df

  # Assign labels: top-N individually, rest → "Other (remaining)"
  other_label <- "Other (remaining)"
  df <- df %>%
    mutate(
      group = if_else(rank_in_est <= TOP_N, direct_comparison, other_label)
    )

  grouped <- df %>%
    group_by(group) %>%
    summarise(
      pct         = sum(pct),
      is_own_direct = any(is_own_direct),
      min_rank    = min(rank_in_est),
      .groups     = "drop"
    ) %>%
    mutate(order = if_else(group == other_label, Inf, as.numeric(min_rank))) %>%
    arrange(order) %>%
    mutate(group_label = str_wrap(group, width = 38))

  # Preserve display order (largest first, Other last)
  grouped$group_label <- factor(grouped$group_label,
                                levels = rev(grouped$group_label))
  list(data = grouped, evidence_type = res$evidence_type,
       top10_only = res$top10_only, total = res$total)
}

make_colours <- function(grouped_df) {
  n      <- nrow(grouped_df)
  labels <- as.character(levels(grouped_df$group_label))
  base_pal <- hcl.colors(max(n, 3), palette = "Dark 2")
  cols <- setNames(base_pal[seq_len(n)], labels)

  # Other → grey
  other_lbl <- as.character(
    grouped_df$group_label[grouped_df$group == "Other (remaining)"])
  if (length(other_lbl) > 0) cols[other_lbl] <- "grey80"

  # Own direct → red
  direct_lbl <- as.character(
    grouped_df$group_label[grouped_df$is_own_direct &
                             grouped_df$group != "Other (remaining)"])
  if (length(direct_lbl) > 0) cols[direct_lbl] <- "#c0392b"

  cols
}

cat("Generating plots …\n")

pdf(OUTPUT_PDF, width = 12, height = 7)

for (treat in all_treatments) {
  nma_mean <- bayes_eff$mean[bayes_eff$class == treat]
  pd  <- make_plot_data(treat)
  gdf <- pd$data
  colours <- make_colours(gdf)

  total_lbl <- if (!pd$top10_only) {
    sprintf("Total abs. contribution = %.4f", pd$total)
  } else {
    sprintf("Top-10 shown (cover %.1f%% of total)", sum(gdf$pct[gdf$group != "Other (remaining)"]))
  }

  subtitle_txt <- sprintf(
    "Evidence type: %s  |  NMA posterior mean = %.3f  |  %s\n%s",
    pd$evidence_type, nma_mean, total_lbl,
    "Direct trial evidence highlighted in red (if present)."
  )

  p <- ggplot(gdf, aes(x = 1, y = pct, fill = group_label)) +
    geom_col(width = 0.55, colour = "white", linewidth = 0.25) +
    geom_text(
      aes(label = if_else(pct >= 2,
                          paste0(round(pct, 1), "%"), "")),
      position = position_stack(vjust = 0.5),
      size = 3, colour = "white", fontface = "bold"
    ) +
    coord_flip() +
    scale_fill_manual(values = colours,
                      name = "Contributing\ndirect comparison") +
    scale_y_continuous(labels = percent_format(scale = 1),
                       expand = expansion(mult = c(0, 0.01))) +
    labs(
      title    = sprintf("%s vs %s", treat, PLACEBO_LABEL),
      subtitle = subtitle_txt,
      x = NULL,
      y = "% of total net contribution"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position    = "right",
      legend.text        = element_text(size = 7),
      legend.key.size    = unit(0.4, "cm"),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 8.5, colour = "grey35")
    ) +
    guides(fill = guide_legend(reverse = TRUE, ncol = 1))

  print(p)
}

dev.off()

cat(sprintf("\nDone. %d plots written to '%s'.\n",
            length(all_treatments), OUTPUT_PDF))
