# =============================================================================
# SSRI Removal Impact Analysis: Effect on All Non-Placebo vs Placebo Estimates
# =============================================================================
# Uses the leave-one-treatment-out (LOTO) data to show how removing all SSRI
# studies from the network changes each non-Placebo vs Placebo effect estimate.
#
# Input files (in the working directory):
#   - loto_detail.csv               (per-treatment impact for each removed arm)
#   - loto_summary.csv              (summary stats per removed treatment)
#   - bayes_relative_effects_vs_placebo_main_labeled.csv  (for ordering)
#
# Output:
#   - Console: detailed text summary for each comparison
#   - PDF:     one forest-style dot-and-CI plot per comparison (baseline vs
#              post-SSRI-removal), plus one overview lollipop chart
# =============================================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

# ── 0. Configuration ──────────────────────────────────────────────────────────

REMOVED_TREATMENT <- "SSRIs"
OUTPUT_PDF        <- "ssri_removal_impact_plots.pdf"
PLACEBO_LABEL     <- "Placebo"

# ── 1. Load data ──────────────────────────────────────────────────────────────

# loto_detail uses semicolons and European decimal commas
read_loto <- function(path) {
  df <- read.csv(path, sep = ";", stringsAsFactors = FALSE, dec = ",")
  df
}

loto_det  <- read_loto("loto_detail.csv")
loto_sum  <- read_loto("loto_summary.csv")
bayes_eff <- read.csv("bayes_relative_effects_vs_placebo_main_labeled.csv",
                      sep = ";", stringsAsFactors = FALSE, dec = ".")

stopifnot(all(c("treat", "TE_base", "lower_base", "upper_base",
                "TE_drop", "lower_drop", "upper_drop",
                "removed_treatment", "delta_TE", "abs_delta_TE",
                "pct_shift_vs_base", "estimable", "ci_overlap") %in%
                names(loto_det)))

# ── 2. Filter to SSRI removal, all non-Placebo treatments ────────────────────

ssri_rows <- loto_det %>%
  filter(removed_treatment == REMOVED_TREATMENT,
         treat != PLACEBO_LABEL) %>%
  arrange(treat)

ssri_summary <- loto_sum %>%
  filter(removed_treatment == REMOVED_TREATMENT)

n_removed_studies   <- ssri_summary$n_removed_studies
n_remaining_studies <- ssri_summary$n_remaining_studies
n_estimable         <- ssri_summary$n_estimable_treatments

cat(strrep("=", 72), "\n")
cat(sprintf("SSRI REMOVAL IMPACT ANALYSIS\n"))
cat(sprintf("Removed treatment:  %s\n", REMOVED_TREATMENT))
cat(sprintf("Studies removed:    %s (%.1f%% of network)\n",
            n_removed_studies,
            (1 - as.numeric(sub(",", ".", ssri_summary$frac_remaining_studies))) * 100))
cat(sprintf("Studies remaining:  %s\n", n_remaining_studies))
cat(sprintf("Estimable contrasts after removal: %s\n", n_estimable))
cat(sprintf("Max |delta|: %s  |  Median |delta|: %s  |  P90 |delta|: %s\n",
            ssri_summary$max_abs_delta,
            ssri_summary$median_abs_delta,
            ssri_summary$p90_abs_delta))
cat(strrep("=", 72), "\n\n")

# ── 3. Order treatments by baseline effect size (largest first) ───────────────

treat_order_df <- bayes_eff %>%
  arrange(desc(mean)) %>%
  select(class, mean)

# Keep only treatments present in loto (some may become inestimable)
treat_order <- treat_order_df$class[
  treat_order_df$class %in% ssri_rows$treat
]

# Treatments that became inestimable after SSRI removal
inestimable <- loto_det %>%
  filter(removed_treatment == REMOVED_TREATMENT,
         treat != PLACEBO_LABEL,
         estimable == FALSE) %>%
  pull(treat)

# ── 4. Text summary ───────────────────────────────────────────────────────────

cat(strrep("=", 72), "\n")
cat("DETAILED SUMMARY: Impact of removing SSRI studies\n")
cat("on each non-Placebo vs Placebo effect estimate\n")
cat(sprintf("(ordered by baseline NMA effect size, largest first)\n"))
cat(strrep("=", 72), "\n\n")

if (length(inestimable) > 0) {
  cat(sprintf("NOTE: %d treatment(s) became inestimable after SSRI removal:\n",
              length(inestimable)))
  for (t in inestimable) cat(sprintf("  - %s\n", t))
  cat("\n")
}

for (treat in treat_order) {
  row <- ssri_rows %>% filter(treat == !!treat)
  if (nrow(row) == 0) next

  nma_base <- row$TE_base
  nma_drop <- row$TE_drop
  delta    <- row$delta_TE
  pct      <- row$pct_shift_vs_base * 100   # already a fraction in source
  ci_ok    <- row$ci_overlap

  direction <- if (delta > 0) "increased" else "decreased"
  ci_msg    <- if (ci_ok) "CIs overlap (estimate stable)" else
                           "CIs do NOT overlap (estimate unstable!)"

  cat(strrep("-", 72), "\n")
  cat(sprintf("%-55s  (base NMA mean = %.3f)\n", paste0(treat, " vs Placebo"), nma_base))
  cat(sprintf("  Baseline effect (full network):  %.3f  [%.3f, %.3f]\n",
              nma_base, row$lower_base, row$upper_base))
  cat(sprintf("  After SSRI removal:              %.3f  [%.3f, %.3f]\n",
              nma_drop, row$lower_drop, row$upper_drop))
  cat(sprintf("  Change (delta):                  %+.3f  (%+.1f%%)\n",
              delta, row$pct_shift_vs_base * 100))
  cat(sprintf("  Direction of shift:              Effect %s after removal\n", direction))
  cat(sprintf("  CI stability:                    %s\n", ci_msg))
  cat("\n")
}

# ── 5. Overview lollipop plot (all treatments, sorted by delta) ───────────────

plot_df <- ssri_rows %>%
  filter(treat %in% treat_order) %>%
  mutate(
    pct_shift_label = paste0(sprintf("%+.1f", pct_shift_vs_base * 100), "%"),
    unstable = !ci_overlap,
    treat_wrap = str_wrap(treat, width = 38)
  ) %>%
  arrange(delta_TE)

# Factor for y-axis order
plot_df$treat_wrap <- factor(plot_df$treat_wrap,
                             levels = unique(plot_df$treat_wrap))

overview <- ggplot(plot_df,
                   aes(x = delta_TE, y = treat_wrap,
                       colour = unstable, shape = unstable)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = delta_TE,
                   y = treat_wrap, yend = treat_wrap),
               colour = "grey70", linewidth = 0.4) +
  geom_point(size = 3) +
  geom_text(aes(label = pct_shift_label),
            hjust = ifelse(plot_df$delta_TE >= 0, -0.25, 1.25),
            size = 2.8, colour = "grey30") +
  scale_colour_manual(
    values = c("FALSE" = "#2c7bb6", "TRUE" = "#d7191c"),
    labels = c("FALSE" = "CIs overlap (stable)", "TRUE" = "CIs do NOT overlap"),
    name = "CI stability"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "CIs overlap (stable)", "TRUE" = "CIs do NOT overlap"),
    name = "CI stability"
  ) +
  labs(
    title    = sprintf("Effect of removing %s studies on all non-Placebo vs Placebo estimates",
                       REMOVED_TREATMENT),
    subtitle = sprintf(
      "%d studies removed (%.1f%% of network)  |  %s studies remaining  |  %d contrasts estimable",
      as.integer(n_removed_studies),
      (1 - as.numeric(sub(",", ".", ssri_summary$frac_remaining_studies))) * 100,
      n_remaining_studies,
      as.integer(n_estimable)
    ),
    x     = sprintf("Change in effect estimate (delta TE) after removing %s", REMOVED_TREATMENT),
    y     = NULL,
    caption = "Positive delta = estimate increased after removal; triangles = 95% CIs no longer overlap."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 8.5, colour = "grey40"),
    plot.caption     = element_text(size = 7.5, colour = "grey50"),
    axis.text.y      = element_text(size = 7.5),
    panel.grid.minor = element_blank()
  )

# ── 6. Per-treatment forest plots (baseline vs post-removal) ──────────────────

make_forest_plot <- function(treat_name, row) {
  nma_mean <- bayes_eff$mean[bayes_eff$class == treat_name]

  # Two rows: baseline and after removal
  forest_df <- data.frame(
    label  = factor(c("After SSRI removal", "Full network (baseline)"),
                    levels = c("After SSRI removal", "Full network (baseline)")),
    mean   = c(row$TE_drop,   row$TE_base),
    lower  = c(row$lower_drop, row$lower_base),
    upper  = c(row$upper_drop, row$upper_base),
    colour = c("drop", "base")
  )

  ci_msg <- if (row$ci_overlap)
    "95% CIs overlap — estimate remains stable"
  else
    "95% CIs do NOT overlap — estimate is UNSTABLE"

  subtitle_txt <- sprintf(
    "Baseline mean = %.3f  |  Post-removal mean = %.3f  |  Delta = %+.3f (%+.1f%%)\n%s",
    row$TE_base, row$TE_drop, row$delta_TE,
    row$pct_shift_vs_base * 100,
    ci_msg
  )

  ggplot(forest_df, aes(x = mean, y = label, colour = colour)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   height = 0.2, linewidth = 0.8) +
    geom_point(size = 4, shape = 18) +
    scale_colour_manual(
      values = c("base" = "#2c7bb6", "drop" = "#d7191c"),
      guide  = "none"
    ) +
    scale_x_continuous(
      # ensure both CIs visible
      limits = function(x) {
        pad <- diff(range(c(forest_df$lower, forest_df$upper))) * 0.15
        c(min(c(forest_df$lower, 0)) - pad,
          max(c(forest_df$upper, 0)) + pad)
      }
    ) +
    labs(
      title    = sprintf("%s vs %s", treat_name, PLACEBO_LABEL),
      subtitle = subtitle_txt,
      x        = "Effect estimate (vs Placebo)",
      y        = NULL,
      caption  = sprintf("Blue = full network; Red = after removing %s studies.",
                         REMOVED_TREATMENT)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey35"),
      plot.caption  = element_text(size = 8, colour = "grey50"),
      axis.text.y   = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
}

# ── 7. Write PDF ──────────────────────────────────────────────────────────────

cat("Generating plots …\n")

pdf(OUTPUT_PDF, width = 12, height = 8)

# Page 1: overview lollipop
print(overview)

# One page per treatment (ordered by baseline effect size)
for (treat in treat_order) {
  row <- ssri_rows %>% filter(treat == !!treat)
  if (nrow(row) == 0) next
  print(make_forest_plot(treat, row))
}

dev.off()

cat(sprintf("\nDone. %d plots written to '%s'.\n",
            1 + length(treat_order), OUTPUT_PDF))
