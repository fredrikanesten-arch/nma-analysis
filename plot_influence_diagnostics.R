#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(igraph)
  library(ggraph)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "."
output_dir <- if (length(args) >= 2) args[[2]] else file.path(data_dir, "diagnostic_plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_sc <- function(path, decimal_mark = ",") {
  read_delim(
    path,
    delim = ";",
    locale = locale(decimal_mark = decimal_mark),
    show_col_types = FALSE,
    trim_ws = TRUE
  )
}

make_key <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = ":"), paste(b, a, sep = ":"))
}

make_key_from_estimate <- function(x) {
  parts <- str_split_fixed(x, ":", 2)
  make_key(parts[, 1], parts[, 2])
}

is_pharmacological <- function(treatment) {
  str_detect(
    treatment,
    paste(
      c("^Any AD$", "^SSRIs$", "^SNRIs$", "^TCAs$", "^Trazodone$", "^Mirtazapine$"),
      collapse = "|"
    )
  )
}

controls <- c("Placebo", "No treatment", "Waitlist", "TAU", "Attention placebo", "Sham acupuncture")

loto_detail <- read_sc(file.path(data_dir, "loto_detail.csv"), decimal_mark = ",")
loto_summary <- read_sc(file.path(data_dir, "loto_summary.csv"), decimal_mark = ",")
loso_summary <- read_sc(file.path(data_dir, "loso_summary.csv"), decimal_mark = ",")
loo_fragility <- read_sc(file.path(data_dir, "loo_most_fragile_network_contrasts.csv"), decimal_mark = ".")
treatment_legend <- read_csv(file.path(data_dir, "treatment_legend_with_counts.csv"), show_col_types = FALSE)
collapsed <- read_csv(file.path(data_dir, "collapsed_class_strict_standardized_with_ids.csv"), show_col_types = FALSE)
edge_table <- read_csv(file.path(data_dir, "nma_edge_table.csv"), show_col_types = FALSE)

treatments <- treatment_legend$class

# 1) Heatmap of LOTO impact
heatmap_df <- loto_detail %>%
  filter(estimable) %>%
  group_by(removed_treatment, treat) %>%
  summarise(abs_delta = max(abs_delta_TE, na.rm = TRUE), .groups = "drop") %>%
  complete(removed_treatment = treatments, treat = treatments, fill = list(abs_delta = 0))

p1 <- ggplot(heatmap_df, aes(x = treat, y = removed_treatment, fill = abs_delta)) +
  geom_tile(color = "grey95", linewidth = 0.1) +
  scale_fill_gradient2(
    low = "white",
    mid = "#fcae91",
    high = "#99000d",
    midpoint = 0.25,
    limits = c(0, 0.5),
    oob = squish,
    name = "|delta_TE|"
  ) +
  labs(
    title = "LOTO impact heatmap",
    x = "Affected treatment",
    y = "Removed treatment"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  )

# Shared study-level metadata
study_sizes <- collapsed %>%
  group_by(studyid) %>%
  summarise(n_participants = sum(n, na.rm = TRUE), .groups = "drop")

study_categories <- collapsed %>%
  distinct(studyid, class) %>%
  mutate(
    active_non_control = !class %in% controls,
    pharm = is_pharmacological(class)
  ) %>%
  group_by(studyid) %>%
  summarise(
    study_category = case_when(
      any(active_non_control & !pharm) ~ "psychological/behavioural",
      any(active_non_control & pharm) ~ "pharmacological",
      TRUE ~ "control only"
    ),
    .groups = "drop"
  )

# 2) LOSO volcano scatter
volcano_df <- loso_summary %>%
  left_join(study_sizes, by = c("removed_study" = "studyid")) %>%
  mutate(log_n = log(pmax(n_participants, 1)))

top_outliers <- volcano_df %>%
  slice_max(order_by = max_abs_delta, n = 8)

p2 <- ggplot(volcano_df, aes(x = log_n, y = max_abs_delta)) +
  geom_point(alpha = 0.7, color = "#2b8cbe") +
  geom_text(
    data = top_outliers,
    aes(label = removed_study),
    size = 2.8,
    vjust = -0.6,
    check_overlap = TRUE
  ) +
  labs(
    title = "LOSO volcano plot",
    x = "Study size (log n participants)",
    y = "Max. absolute change in treatment effect"
  ) +
  theme_minimal(base_size = 11)

# 3) LOO fragility bar chart
loo_plot_df <- loo_fragility %>%
  separate(network_estimate, into = c("treat1", "treat2"), sep = ":", remove = FALSE) %>%
  mutate(
    contrast_category = if_else(
      is_pharmacological(treat1) & is_pharmacological(treat2),
      "pharmacological",
      "psychological"
    ),
    network_estimate = reorder(network_estimate, max_abs_change)
  )

p3 <- ggplot(loo_plot_df, aes(x = network_estimate, y = max_abs_change, fill = contrast_category)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(pharmacological = "#3182bd", psychological = "#f16913")) +
  labs(
    title = "LOO fragility by network contrast",
    x = "Network contrast",
    y = "Max. absolute change in treatment effect",
    fill = "Category"
  ) +
  theme_minimal(base_size = 11)

# 4) Influence network diagram
edge_weights <- loo_fragility %>%
  transmute(contrast_key = make_key_from_estimate(network_estimate), loo_weight = max_abs_change) %>%
  group_by(contrast_key) %>%
  summarise(loo_weight = max(loo_weight, na.rm = TRUE), .groups = "drop")

edges_plot <- edge_table %>%
  mutate(contrast_key = make_key(from, to)) %>%
  left_join(edge_weights, by = "contrast_key") %>%
  mutate(loo_weight = replace_na(loo_weight, 0))

node_influence <- loto_summary %>%
  filter(component_id == 1) %>%
  group_by(removed_treatment) %>%
  summarise(loto_max_delta = max(max_abs_delta, na.rm = TRUE), .groups = "drop")

node_tier <- treatment_legend %>%
  mutate(
    tier_num = ntile(n_studies, 3),
    evidence_tier = recode(tier_num, `1` = "low", `2` = "medium", `3` = "high")
  ) %>%
  select(class, evidence_tier)

nodes_plot <- tibble(name = sort(unique(c(edges_plot$from, edges_plot$to)))) %>%
  left_join(node_influence, by = c("name" = "removed_treatment")) %>%
  left_join(node_tier, by = c("name" = "class")) %>%
  mutate(
    loto_max_delta = replace_na(loto_max_delta, 0),
    evidence_tier = replace_na(evidence_tier, "low")
  )

graph_obj <- graph_from_data_frame(
  d = edges_plot %>% select(from, to, loo_weight),
  vertices = nodes_plot,
  directed = FALSE
)

set.seed(123)
p4 <- ggraph(graph_obj, layout = "fr") +
  geom_edge_link(aes(width = loo_weight), alpha = 0.35, colour = "grey45") +
  geom_node_point(aes(size = loto_max_delta, colour = evidence_tier), alpha = 0.9) +
  geom_node_text(aes(label = name), size = 2.4, repel = TRUE) +
  scale_edge_width(range = c(0.2, 2.5), guide = guide_legend(title = "LOO fragility")) +
  scale_size_continuous(range = c(2, 12), name = "LOTO max_Δ") +
  scale_colour_manual(values = c(low = "#9ecae1", medium = "#6baed6", high = "#2171b5")) +
  labs(title = "Influence network diagram", colour = "Evidence tier") +
  theme_void(base_size = 11)

# 5) LOSO ranked dot/caterpillar plot
loso_rank_df <- loso_summary %>%
  left_join(study_categories, by = c("removed_study" = "studyid")) %>%
  arrange(desc(max_abs_delta)) %>%
  mutate(rank = row_number())

p5 <- ggplot(loso_rank_df, aes(x = rank, y = max_abs_delta, colour = study_category)) +
  geom_segment(aes(xend = rank, y = 0, yend = max_abs_delta), alpha = 0.2) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = c(
    pharmacological = "#3182bd",
    `psychological/behavioural` = "#f16913",
    `control only` = "#636363"
  )) +
  labs(
    title = "LOSO ranked max_abs_delta (caterpillar)",
    x = "Study rank (1 = most influential)",
    y = "Max. absolute change in treatment effect",
    colour = "Study category"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(output_dir, "plot1_loto_heatmap.png"), p1, width = 14, height = 12, dpi = 300)
ggsave(file.path(output_dir, "plot2_loso_volcano.png"), p2, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "plot3_loo_fragility_bars.png"), p3, width = 11, height = 10, dpi = 300)
ggsave(file.path(output_dir, "plot4_influence_network.png"), p4, width = 13, height = 11, dpi = 300)
ggsave(file.path(output_dir, "plot5_loso_ranked_caterpillar.png"), p5, width = 12, height = 6, dpi = 300)

message("Saved all 5 plots to: ", normalizePath(output_dir))
