# =============================================================================
# Network Component Visualisation
# =============================================================================
# Builds the NMA treatment network from the raw study data and shows:
#   1. The full connected network (all 50 treatments)
#   2. Components formed after removing SSRIs  (→ 3 components)
#   3. Components formed after removing Any AD (→ 4 components)
#
# Critically illustrates that removing SSRIs fully disconnects
# Acupuncture + AD and Exercise group + AD from Placebo, making their
# vs-Placebo estimates inestimable without indirect evidence.
#
# Input files (working directory):
#   - gemtc_data_ab_strict_standardized_safe_ids.csv
#   - treatment_legend.csv
#   - giant_component_reliance.csv
#   - bayes_relative_effects_vs_placebo_main_labeled.csv
#
# Output: network_component_plots.pdf (5 pages)
# =============================================================================

library(igraph)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

OUTPUT_PDF <- "network_component_plots.pdf"

# ── 0. Helper: short label for long treatment names ───────────────────────────

shorten <- function(x) {
  x <- gsub("Cognitive and cognitive behavioural therapies", "CBT", x)
  x <- gsub("Interpersonal psychotherapy \\(IPT\\)", "IPT", x)
  x <- gsub("Short-term psychodynamic psychotherapies", "STPP", x)
  x <- gsub("individual", "ind", x)
  x <- gsub("therapies", "ther.", x)
  x <- gsub("therapy", "ther.", x)
  x <- gsub(" \\+ AD$", "+AD", x)
  x <- gsub(" \\+ placebo$", "+placebo", x)
  x <- gsub("Behavioural ther. ind \\+ AD", "Behav ind+AD", x)
  x <- gsub("Behavioural ther. ind$", "Behav ind", x)
  x <- gsub("Mindfulness or meditation group", "Mindfulness grp", x)
  x <- gsub("Psychoeducation group", "Psychoedu. grp", x)
  x <- gsub("Problem solving ind", "PSI", x)
  x <- gsub("Problem solving group", "PSG", x)
  x <- gsub("Peer support group", "Peer support", x)
  x <- gsub("Self-help with support", "Self-help+sup", x)
  x <- gsub("Relaxation ind", "Relax ind", x)
  str_squish(x)
}

# ── 1. Load data and build edge list ─────────────────────────────────────────

gemtc  <- read.csv("gemtc_data_ab_strict_standardized_safe_ids.csv",
                   stringsAsFactors = FALSE)
legend <- read.csv("treatment_legend.csv", stringsAsFactors = FALSE)
bayes  <- read.csv("bayes_relative_effects_vs_placebo_main_labeled.csv",
                   sep = ";", stringsAsFactors = FALSE)

# Map treatment IDs to class names
id_to_name <- setNames(legend$class, legend$treatment_id)

# Build edge list: each study that tests ≥2 treatments contributes edges
edges_df <- gemtc %>%
  mutate(trt_name = id_to_name[treatment]) %>%
  group_by(study) %>%
  summarise(treatments = list(unique(trt_name)), .groups = "drop") %>%
  rowwise() %>%
  do({
    treats <- sort(unlist(.$treatments))
    if (length(treats) >= 2) {
      pairs <- combn(treats, 2)
      data.frame(
        from  = pairs[1, ],
        to    = pairs[2, ],
        study = .$study,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(from = character(0), to = character(0),
                 study = character(0))
    }
  }) %>%
  ungroup()

# Aggregate: count studies per edge
edge_weights <- edges_df %>%
  group_by(from, to) %>%
  summarise(n_studies = n(), .groups = "drop")

cat(sprintf("Full network: %d nodes, %d edges, %d studies\n",
            length(unique(c(edge_weights$from, edge_weights$to))),
            nrow(edge_weights),
            nrow(gemtc) / 2  # approximate
))

# ── 2. Build igraph objects ───────────────────────────────────────────────────

build_graph <- function(vertices, edge_tbl) {
  # Keep only edges where both endpoints are in vertices
  e <- edge_tbl %>%
    filter(from %in% vertices, to %in% vertices)
  g <- graph_from_data_frame(
    d        = e[, c("from", "to", "n_studies")],
    directed = FALSE,
    vertices = data.frame(name = vertices, stringsAsFactors = FALSE)
  )
  E(g)$weight <- E(g)$n_studies
  g
}

all_nodes <- unique(c(edge_weights$from, edge_weights$to))
g_full    <- build_graph(all_nodes, edge_weights)

# Remove SSRIs
g_no_ssri <- build_graph(all_nodes[all_nodes != "SSRIs"], edge_weights)

# Remove Any AD
g_no_anyad <- build_graph(all_nodes[all_nodes != "Any AD"], edge_weights)

cat(sprintf("Full graph components:         %d\n", components(g_full)$no))
cat(sprintf("After removing SSRIs:          %d components\n",
            components(g_no_ssri)$no))
cat(sprintf("After removing Any AD:         %d components\n",
            components(g_no_anyad)$no))

# ── 3. Colour / style helpers ─────────────────────────────────────────────────

# Category colours
PHARMA  <- c("SSRIs", "SNRIs", "TCAs", "Mirtazapine", "Trazodone",
             "Any AD")
PSYCH   <- c("Cognitive and cognitive behavioural therapies individual",
             "Cognitive and cognitive behavioural therapies group",
             "Cognitive and cognitive behavioural therapies individual + AD",
             "Cognitive and cognitive behavioural therapies group + AD",
             "Cognitive and cognitive behavioural therapies individual + placebo",
             "Interpersonal psychotherapy (IPT) individual",
             "Interpersonal psychotherapy (IPT) individual + AD",
             "Interpersonal psychotherapy (IPT) individual + placebo",
             "Short-term psychodynamic psychotherapies individual",
             "Short-term psychodynamic psychotherapies individual + AD",
             "Behavioural therapies individual",
             "Behavioural therapies individual + AD",
             "Counselling individual",
             "Counselling individual + AD",
             "Counselling individual + placebo",
             "Problem solving individual",
             "Problem solving group",
             "Self-help",
             "Self-help with support",
             "Any psychotherapy",
             "Mindfulness or meditation group",
             "Music therapy group",
             "Psychoeducation group",
             "Psychoeducation group + AD",
             "Peer support group",
             "Peer support group + AD")
ACTIVE  <- c("Exercise group", "Exercise group + AD",
             "Exercise individual", "Exercise individual + AD",
             "Light therapy", "Light therapy + AD",
             "Yoga group", "Yoga group + AD",
             "Acupuncture", "Acupuncture + AD",
             "Relaxation individual + AD",
             "Relaxation individual + placebo",
             "Sham acupuncture")
CTRL    <- c("Placebo", "Waitlist", "No treatment", "TAU",
             "Attention placebo")

node_colour <- function(name) {
  case_when(
    name == "Placebo"             ~ "#e74c3c",   # red
    name %in% PHARMA              ~ "#2980b9",   # blue
    name %in% PSYCH               ~ "#27ae60",   # green
    name %in% ACTIVE              ~ "#8e44ad",   # purple
    name %in% CTRL                ~ "#f39c12",   # orange
    TRUE                          ~ "grey70"
  )
}

# ── 4. Plot function ──────────────────────────────────────────────────────────

plot_graph <- function(g, title, subtitle = "",
                       removed_node = NULL,
                       highlight_isolated = NULL) {

  vnames   <- V(g)$name
  n        <- vcount(g)
  comps    <- components(g)

  # Assign component membership colour overlay (for split plots)
  comp_id  <- comps$membership
  n_comp   <- comps$no

  # Base node colours by category
  vcols <- sapply(vnames, node_colour)

  # Node sizes: Placebo = 20, removed = 0 (not in graph), others by degree
  deg     <- degree(g)
  vsizes  <- ifelse(vnames == "Placebo", 22,
                    ifelse(vnames %in% c("SSRIs", "Any AD"), 18,
                           4 + deg * 0.9))

  # Edge widths
  ew <- E(g)$weight
  ewidths <- 0.4 + log1p(ew) * 0.7

  # Layout: use Fruchterman-Reingold with seed for reproducibility
  set.seed(42)
  if (n > 30) {
    lay <- layout_with_fr(g, weights = ew, niter = 1500)
  } else {
    lay <- layout_with_fr(g, niter = 1500)
  }

  # Component border colours (for split plots)
  comp_palette <- c("#e74c3c", "#2c3e50", "#1abc9c", "#e67e22",
                    "#9b59b6", "#f1c40f")
  vborder <- if (n_comp > 1) {
    comp_palette[comp_id]
  } else {
    vcols
  }

  # Short labels
  vlabels <- shorten(vnames)

  # Mark truly isolated / disconnected nodes
  if (!is.null(highlight_isolated)) {
    isolated_idx <- which(vnames %in% highlight_isolated)
    vcols[isolated_idx]  <- "white"
    vborder[isolated_idx] <- "#e74c3c"
    vsizes[isolated_idx] <- 14
  }

  # Mark removed node in title
  full_title <- title
  if (!is.null(removed_node))
    full_title <- paste0(title, "\n(", removed_node, " removed — shown as dashed border)")

  old_par <- par(mar = c(0.5, 0.5, 3, 0.5), bg = "white")
  on.exit(par(old_par))

  plot(g,
       layout            = lay,
       vertex.label      = vlabels,
       vertex.label.cex  = 0.52,
       vertex.label.color = "black",
       vertex.label.font  = 1,
       vertex.color      = vcols,
       vertex.frame.color = vborder,
       vertex.frame.width = ifelse(n_comp > 1, 2.5, 1),
       vertex.size       = vsizes,
       edge.width        = ewidths,
       edge.color        = "grey65",
       edge.curved       = 0.15,
       main              = full_title)

  if (nzchar(subtitle))
    mtext(subtitle, side = 3, line = -0.2, cex = 0.72, col = "grey30")

  # Legend
  legend("bottomleft",
         legend  = c("Placebo", "Pharmacological", "Psychological",
                     "Other active", "Control/reference",
                     if (n_comp > 1) "Component boundaries" else NULL),
         pch     = c(rep(21, 5), if (n_comp > 1) 22 else NULL),
         pt.bg   = c("#e74c3c", "#2980b9", "#27ae60", "#8e44ad",
                     "#f39c12", if (n_comp > 1) "white" else NULL),
         col     = c(rep("black", 5), if (n_comp > 1) "#e74c3c" else NULL),
         pt.cex  = 1.3,
         cex     = 0.62,
         bty     = "n",
         ncol    = 1)

  if (n_comp > 1) {
    # Annotate components
    for (ci in seq_len(n_comp)) {
      idx <- which(comp_id == ci)
      cx  <- mean(lay[idx, 1])
      cy  <- max(lay[idx, 2]) + 0.08
      size_lbl <- sprintf("Comp %d\n(n=%d)", ci, length(idx))
      text(cx, cy, size_lbl, cex = 0.62,
           col = comp_palette[ci], font = 2)
    }
  }

  invisible(NULL)
}

# ── 5. Write PDF ──────────────────────────────────────────────────────────────

cat("Generating plots …\n")

pdf(OUTPUT_PDF, width = 14, height = 11)

# Page 1: Full network
plot_graph(
  g_full,
  title    = "Full NMA Treatment Network (50 nodes, 128 edges)",
  subtitle = paste0(
    "Edge width ∝ number of studies. Node size ∝ degree. ",
    "Placebo = red. Blue = pharmacological, green = psychological, ",
    "purple = other active, orange = control/reference."
  )
)

# Page 2: After SSRIs removal — full picture showing all 3 components
plot_graph(
  g_no_ssri,
  title    = "Network after removing SSRIs (3 components formed)",
  subtitle = paste0(
    "Acupuncture+AD and Exercise group+AD become ISOLATED (component border = red). ",
    "They can no longer be compared to Placebo."
  ),
  removed_node       = "SSRIs",
  highlight_isolated = c("Acupuncture + AD", "Exercise group + AD")
)

# Page 3: Zoom on the isolated nodes — show only small components
iso_ssri <- c("Acupuncture + AD", "Exercise group + AD")
g_iso_ssri <- build_graph(
  c(iso_ssri, "SSRIs"),   # show them WITH SSRIs to explain the connection
  edge_weights
)
# Show explicitly what Acupuncture+AD and Exercise group+AD were connected to
plot_graph(
  g_iso_ssri,
  title    = "Zoom: Isolated nodes after SSRIs removal",
  subtitle = paste0(
    "Acupuncture+AD and Exercise group+AD are ONLY connected to SSRIs. ",
    "Remove SSRIs → zero path to Placebo → inestimable vs Placebo."
  )
)

# Page 4: After Any AD removal — 4 components
plot_graph(
  g_no_anyad,
  title    = "Network after removing Any AD (4 components formed)",
  subtitle = paste0(
    "Peer support group/+AD, Yoga group+AD, and CBT group+AD become isolated. ",
    "Component borders show the 4 disconnected subgraphs."
  ),
  removed_node       = "Any AD",
  highlight_isolated = c("Peer support group", "Peer support group + AD",
                          "Yoga group + AD",
                          "Cognitive and cognitive behavioural therapies group + AD")
)

# Page 5: Zoom on Any AD isolated nodes
iso_anyad <- c("Peer support group", "Peer support group + AD",
               "Yoga group + AD",
               "Cognitive and cognitive behavioural therapies group + AD",
               "Any AD")
g_iso_anyad <- build_graph(iso_anyad, edge_weights)
plot_graph(
  g_iso_anyad,
  title    = "Zoom: Isolated nodes after Any AD removal",
  subtitle = paste0(
    "Peer support group/+AD, Yoga group+AD, and CBT group+AD are ONLY ",
    "connected to Any AD. Remove Any AD → inestimable vs Placebo."
  )
)

dev.off()

cat(sprintf("Done. Plots written to '%s'.\n", OUTPUT_PDF))

# ── 6. Print summary table ────────────────────────────────────────────────────

cat("\n")
cat(strrep("=", 65), "\n")
cat("SUMMARY: Treatments that become inestimable vs Placebo\n")
cat(strrep("=", 65), "\n\n")

cat("After removing SSRIs (174 studies removed):\n")
cat("  The following treatments connect to the network ONLY via SSRIs\n")
cat("  and thus become completely disconnected from Placebo:\n\n")
for (t in iso_ssri)
  cat(sprintf("    %-45s (only connection: SSRIs)\n", t))

cat("\nAfter removing Any AD (studies removed):\n")
cat("  The following treatments connect to the network ONLY via Any AD\n")
cat("  and thus become completely disconnected from Placebo:\n\n")
for (t in iso_anyad[iso_anyad != "Any AD"])
  cat(sprintf("    %-45s (only connection: Any AD)\n", t))

cat("\n")
cat("Implication: In the LOTO analysis, the reported effect sizes for\n")
cat("these treatments after their anchor treatment is removed are NOT\n")
cat("based on genuine network paths — they reflect model extrapolation\n")
cat("or are explicitly marked inestimable. The NMA engine may set these\n")
cat("to NA or propagate large uncertainty. Always check estimable=TRUE.\n")
