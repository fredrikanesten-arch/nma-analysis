suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(netmeta)
})

read_csv_robust <- function(path) {
  dat_csv <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(dat_csv) && ncol(dat_csv) > 1) return(dat_csv)
  dat_scsv <- tryCatch(readr::read_delim(path, delim = ";", show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(dat_scsv) && ncol(dat_scsv) > 1) return(dat_scsv)
  stop("Could not parse CSV: ", path)
}

safe_char <- function(x) {
  out <- as.character(x)
  out[is.na(out)] <- ""
  out
}

infer_treatment_node <- function(dat) {
  label_candidates <- c("treatment_label", "trt_label", "treatment_name", "intervention", "trtname")
  label_col <- label_candidates[label_candidates %in% names(dat)][1]

  if (!is.na(label_col) && !is.null(label_col)) {
    node <- safe_char(dat[[label_col]])
    node[node == ""] <- NA_character_
  } else {
    node <- rep(NA_character_, nrow(dat))
  }

  if ("treatment" %in% names(dat)) {
    trt_num <- suppressWarnings(as.numeric(dat$treatment))
    trt_chr <- safe_char(dat$treatment)
    node_from_trt <- ifelse(!is.na(trt_num), paste0("trt_", sprintf("%03d", as.integer(trt_num))), trt_chr)
    node[is.na(node) | !nzchar(node)] <- node_from_trt[is.na(node) | !nzchar(node)]
  }

  if (!"class" %in% names(dat)) {
    stop("Input must include class or be mappable to class.")
  }

  fallback <- safe_char(dat$class)
  fallback[fallback == ""] <- NA_character_
  node[is.na(node) | !nzchar(node)] <- fallback[is.na(node) | !nzchar(node)]

  if (any(is.na(node) | !nzchar(node))) {
    bad <- which(is.na(node) | !nzchar(node))[1]
    stop("Could not infer treatment node on row ", bad)
  }

  node
}

build_analysis_node <- function(dat, mode, bridge_classes = c("Any AD", "Any psychotherapy")) {
  mode <- match.arg(mode, c("class", "treatment", "hybrid"))
  trt_node <- infer_treatment_node(dat)
  class_chr <- safe_char(dat$class)

  if (mode == "class") return(class_chr)
  if (mode == "treatment") return(trt_node)

  use_bridge <- class_chr %in% bridge_classes
  ifelse(use_bridge, class_chr, trt_node)
}

collapse_same_node <- function(dat, node_col = "analysis_node") {
  if (!node_col %in% names(dat)) stop("Missing node column: ", node_col)

  dat %>%
    group_by(studyid, .data[[node_col]]) %>%
    summarise(
      n = sum(n),
      mean_change = sum(mean_change * n) / sum(n),
      sd_change = {
        n_total <- sum(n)
        mu <- sum(mean_change * n) / n_total
        ss_within <- sum((n - 1) * (sd_change^2))
        ss_between <- sum(n * (mean_change - mu)^2)
        if (n_total <= 1) NA_real_ else sqrt((ss_within + ss_between) / (n_total - 1))
      },
      class = if (dplyr::n_distinct(class) == 1) dplyr::first(class) else "mixed",
      treatment_node = if (dplyr::n_distinct(treatment_node) == 1) dplyr::first(treatment_node) else "mixed",
      source = dplyr::first(source),
      control_class = dplyr::first(control_class),
      control_sd_ref = dplyr::first(control_sd_ref),
      control_sd_source = dplyr::first(control_sd_source),
      .groups = "drop"
    ) %>%
    rename(analysis_node = .data[[node_col]])
}

compute_components <- function(nodes, edge_df) {
  nodes <- unique(as.character(nodes))
  nodes <- nodes[nzchar(nodes)]
  if (length(nodes) == 0) {
    return(list(n_components = 0L, giant_component_size = 0L, component_df = data.frame(node = character(0), component = integer(0))))
  }

  adj <- stats::setNames(vector("list", length(nodes)), nodes)
  if (!is.null(edge_df) && nrow(edge_df) > 0) {
    for (i in seq_len(nrow(edge_df))) {
      a <- as.character(edge_df$treat1[i])
      b <- as.character(edge_df$treat2[i])
      if (!nzchar(a) || !nzchar(b) || is.na(a) || is.na(b) || a == b) next
      if (!a %in% nodes || !b %in% nodes) next
      adj[[a]] <- unique(c(adj[[a]], b))
      adj[[b]] <- unique(c(adj[[b]], a))
    }
  }

  visited <- stats::setNames(rep(FALSE, length(nodes)), nodes)
  comp_idx <- integer(0)
  comp_nodes <- character(0)
  comp_id <- 0L

  for (start in nodes) {
    if (visited[[start]]) next
    comp_id <- comp_id + 1L
    queue <- c(start)
    visited[[start]] <- TRUE

    while (length(queue) > 0) {
      v <- queue[[1]]
      queue <- queue[-1]
      comp_nodes <- c(comp_nodes, v)
      comp_idx <- c(comp_idx, comp_id)

      nbrs <- adj[[v]]
      if (length(nbrs) == 0) next
      for (u in nbrs) {
        if (!visited[[u]]) {
          visited[[u]] <- TRUE
          queue <- c(queue, u)
        }
      }
    }
  }

  comp_tbl <- data.frame(node = comp_nodes, component = comp_idx, stringsAsFactors = FALSE)
  sizes <- comp_tbl %>% count(component, name = "size")

  list(
    n_components = nrow(sizes),
    giant_component_size = if (nrow(sizes) == 0) 0L else max(sizes$size),
    component_df = comp_tbl
  )
}

make_nma_summary <- function(nma, collapsed_data, ref, digits = 3) {
  if (!ref %in% nma$trts) stop("Reference treatment not found: ", ref)

  te <- if (!is.null(nma$TE.random)) nma$TE.random else nma$TE.common
  lo <- if (!is.null(nma$lower.random)) nma$lower.random else nma$lower.common
  hi <- if (!is.null(nma$upper.random)) nma$upper.random else nma$upper.common
  se <- if (!is.null(nma$seTE.random)) nma$seTE.random else nma$seTE.common

  pick_vs_ref <- function(m) {
    if (ref %in% colnames(m)) return(data.frame(analysis_node = rownames(m), value = as.numeric(m[, ref]), row.names = NULL))
    if (ref %in% rownames(m)) return(data.frame(analysis_node = colnames(m), value = as.numeric(m[ref, ]), row.names = NULL))
    stop("Reference not found in matrix")
  }

  n_tbl <- collapsed_data %>%
    group_by(analysis_node) %>%
    summarise(
      N = sum(n),
      class_count = n_distinct(class),
      class = ifelse(class_count == 1, first(class), "mixed"),
      .groups = "drop"
    ) %>%
    select(-class_count)

  out <- pick_vs_ref(te) %>%
    rename(estimate = value) %>%
    left_join(rename(pick_vs_ref(lo), lower = value), by = "analysis_node") %>%
    left_join(rename(pick_vs_ref(hi), upper = value), by = "analysis_node") %>%
    left_join(rename(pick_vs_ref(se), se = value), by = "analysis_node") %>%
    left_join(n_tbl, by = "analysis_node") %>%
    mutate(
      estimate = ifelse(analysis_node == ref, 0, estimate),
      lower = ifelse(analysis_node == ref, 0, lower),
      upper = ifelse(analysis_node == ref, 0, upper),
      se = ifelse(analysis_node == ref, NA_real_, se),
      z = ifelse(!is.na(se) & se > 0, estimate / se, NA_real_),
      p_value = ifelse(!is.na(z), 2 * pnorm(-abs(z)), NA_real_)
    ) %>%
    arrange(estimate)

  fmt <- function(x) formatC(x, digits = digits, format = "f")
  out$effect_95CI <- paste0(fmt(out$estimate), " (", fmt(out$lower), " to ", fmt(out$upper), ")")
  out
}

summarize_by_class_posthoc <- function(node_summary) {
  node_summary %>%
    filter(!is.na(class), class != "mixed") %>%
    group_by(class) %>%
    summarise(
      n_nodes = n(),
      N = sum(N, na.rm = TRUE),
      estimate_weighted = weighted.mean(estimate, w = pmax(N, 1), na.rm = TRUE),
      lower_weighted = weighted.mean(lower, w = pmax(N, 1), na.rm = TRUE),
      upper_weighted = weighted.mean(upper, w = pmax(N, 1), na.rm = TRUE),
      estimate_min = min(estimate, na.rm = TRUE),
      estimate_max = max(estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(estimate_weighted)
}

get_largest_component_nodes <- function(comp_obj) {
  if (is.null(comp_obj$component_df) || nrow(comp_obj$component_df) == 0) return(character(0))
  sizes <- comp_obj$component_df %>% count(component, name = "size") %>% arrange(desc(size), component)
  keep_component <- sizes$component[[1]]
  comp_obj$component_df %>% filter(component == keep_component) %>% pull(node)
}

run_variant <- function(dat_std, variant_name, mode, reference_group, out_dir,
                        bridge_classes = c("Any AD", "Any psychotherapy"),
                        force_lcc = FALSE) {
  dat_mode <- dat_std %>%
    mutate(
      analysis_mode = mode,
      analysis_node = build_analysis_node(., mode = mode, bridge_classes = bridge_classes)
    )

  study_node_counts <- dat_mode %>% count(studyid, analysis_node, name = "arm_rows") %>% count(studyid, name = "n_nodes")
  dropped_single <- study_node_counts %>% filter(n_nodes < 2)

  loss_summary <- bind_rows(
    dat_std %>% mutate(reason = ifelse(drop_reason == "kept", "precollapse_kept", drop_reason)),
    dat_mode %>%
      semi_join(dropped_single, by = "studyid") %>%
      mutate(reason = "single_node_after_collapse")
  ) %>%
    group_by(reason) %>%
    summarise(
      n_rows = n(),
      n_studies = n_distinct(studyid),
      participants = sum(n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_studies), reason)

  collapsed <- collapse_same_node(dat_mode, node_col = "analysis_node")
  collapsed_pre_lcc <- collapsed %>% group_by(studyid) %>% filter(n_distinct(analysis_node) >= 2) %>% ungroup()

  edge_df <- collapsed_pre_lcc %>%
    group_by(studyid) %>%
    summarise(
      edges = list({
        nd <- unique(analysis_node)
        if (length(nd) < 2) return(data.frame(treat1 = character(0), treat2 = character(0)))
        cm <- t(combn(nd, 2))
        data.frame(treat1 = cm[, 1], treat2 = cm[, 2], stringsAsFactors = FALSE)
      }),
      .groups = "drop"
    ) %>%
    pull(edges) %>%
    bind_rows()

  component_obj <- compute_components(nodes = unique(collapsed_pre_lcc$analysis_node), edge_df = edge_df)

  collapsed_use <- collapsed_pre_lcc
  if (force_lcc && nrow(collapsed_pre_lcc) > 0) {
    keep_nodes <- get_largest_component_nodes(component_obj)
    collapsed_use <- collapsed_pre_lcc %>% filter(analysis_node %in% keep_nodes)
  }

  if (!reference_group %in% collapsed_use$analysis_node) {
    stop("Reference group '", reference_group, "' is not present in variant '", variant_name, "'.")
  }

  collapsed_use <- collapsed_use %>% group_by(studyid) %>% filter(n_distinct(analysis_node) >= 2) %>% ungroup()
  if (nrow(collapsed_use) == 0) stop("No usable data for variant '", variant_name, "'.")

  pw <- pairwise(
    treat = analysis_node,
    mean = mean_change,
    sd = sd_change,
    n = n,
    studlab = studyid,
    data = collapsed_use,
    sm = "MD"
  )

  nma <- netmeta(
    TE = pw$TE,
    seTE = pw$seTE,
    treat1 = pw$treat1,
    treat2 = pw$treat2,
    studlab = pw$studlab,
    data = pw,
    sm = "MD",
    random = TRUE,
    common = FALSE,
    reference.group = reference_group
  )

  node_summary <- make_nma_summary(nma = nma, collapsed_data = collapsed_use, ref = reference_group)
  class_summary <- summarize_by_class_posthoc(node_summary)

  graph_summary <- data.frame(
    variant = variant_name,
    mode = mode,
    bridge_classes = paste(bridge_classes, collapse = "|"),
    n_studies_input = dplyr::n_distinct(dat_std$studyid),
    n_studies_after_precollapse_filter = dplyr::n_distinct(collapsed_pre_lcc$studyid),
    n_studies_analyzed = dplyr::n_distinct(collapsed_use$studyid),
    n_nodes = dplyr::n_distinct(collapsed_use$analysis_node),
    n_components = component_obj$n_components,
    giant_component_size = component_obj$giant_component_size,
    force_lcc = force_lcc,
    stringsAsFactors = FALSE
  )

  prefix <- file.path(out_dir, paste0("analysis_", variant_name))

  write_csv(collapsed_use, paste0(prefix, "_collapsed.csv"))
  write_csv(loss_summary, paste0(prefix, "_loss_accounting.csv"))
  write_csv(component_obj$component_df, paste0(prefix, "_components.csv"))
  write_csv(graph_summary, paste0(prefix, "_graph_summary.csv"))
  write_csv(node_summary, paste0(prefix, "_nma_summary_nodes.csv"))
  write_csv(class_summary, paste0(prefix, "_nma_summary_classes_posthoc.csv"))

  ns <- tryCatch(netsplit(nma), error = function(e) NULL)
  if (!is.null(ns)) {
    ns_df <- tryCatch(as.data.frame(ns), error = function(e) NULL)
    if (!is.null(ns_df)) write_csv(ns_df, paste0(prefix, "_netsplit.csv"))
  }

  nc <- tryCatch(netcontrib(nma), error = function(e) NULL)
  if (!is.null(nc)) {
    nc_df <- tryCatch(as.data.frame(nc), error = function(e) NULL)
    if (!is.null(nc_df)) write_csv(nc_df, paste0(prefix, "_netcontrib.csv"))
  }

  dd <- tryCatch(decomp.design(nma), error = function(e) NULL)
  if (!is.null(dd)) {
    capture.output(print(dd), file = paste0(prefix, "_decomp_design.txt"))
  }

  saveRDS(
    list(
      variant = variant_name,
      mode = mode,
      bridge_classes = bridge_classes,
      loss_summary = loss_summary,
      components = component_obj,
      graph_summary = graph_summary,
      collapsed_use = collapsed_use,
      pw = pw,
      nma = nma,
      node_summary = node_summary,
      class_summary = class_summary,
      netsplit = ns,
      netcontrib = nc,
      decomp_design = dd
    ),
    paste0(prefix, "_objects.rds")
  )

  graph_summary
}

main <- function() {
  # -----------------------------
  # 0. Paths
  # -----------------------------
  base_dir <- "C:/Users/fredr/OneDrive/Desktop/nma_project/mavranezouli"
  in_data <- file.path(base_dir, "binfixed_class_ms", "combined_long_mean_change_dataset_ms_smd_bias_adj.csv")
  in_map <- file.path(base_dir, "clean_data", "trt_to_class_ms.csv")
  out_dir <- file.path(base_dir, "binfixed_class_ms")

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[[1]] %in% c("-h", "--help")) {
    cat("Usage: Rscript netmeta_class_ms_control_arm_baseline_sd.R [input_csv] [mapping_csv] [output_dir] [reference] [fallback] [mode] [bridge]")
    cat("\nmode: class | treatment | hybrid | multi (default: multi)")
    cat("\nbridge (for hybrid): any_ad_only | any_ad_any_psychotherapy (default: any_ad_any_psychotherapy)\n")
    cat("fallback: drop (default) or sd_change when control baseline SD is missing.\n")
    return(invisible(NULL))
  }

  in_data <- if (length(args) >= 1) args[[1]] else in_data
  in_map <- if (length(args) >= 2) args[[2]] else in_map
  out_dir <- if (length(args) >= 3) args[[3]] else out_dir
  reference_group <- if (length(args) >= 4) args[[4]] else "Placebo"
  fallback_mode <- if (length(args) >= 5) args[[5]] else "drop"
  run_mode <- if (length(args) >= 6) args[[6]] else "multi"
  bridge_key <- if (length(args) >= 7) args[[7]] else "any_ad_any_psychotherapy"

  if (!fallback_mode %in% c("drop", "sd_change")) {
    stop("Invalid fallback mode. Use 'drop' or 'sd_change'.")
  }

  if (!run_mode %in% c("class", "treatment", "hybrid", "multi")) {
    stop("Invalid mode. Use class, treatment, hybrid, or multi.")
  }

  bridge_sets <- list(
    any_ad_only = c("Any AD"),
    any_ad_any_psychotherapy = c("Any AD", "Any psychotherapy")
  )
  if (!bridge_key %in% names(bridge_sets)) {
    stop("Invalid bridge key. Use any_ad_only or any_ad_any_psychotherapy.")
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!file.exists(in_data)) stop("Input data not found: ", in_data)

  dat <- read_csv_robust(in_data)

  if (!"class" %in% names(dat)) {
    if (!file.exists(in_map)) stop("Mapping file required when class column is absent: ", in_map)
    trt_map <- read_csv_robust(in_map)
    req_map <- c("trtcode", "class")
    if (!all(req_map %in% names(trt_map))) stop("Mapping file missing columns: trtcode, class")

    trt_map_use <- trt_map %>%
      mutate(trtcode = suppressWarnings(as.numeric(trtcode)), class = as.character(class)) %>%
      filter(!is.na(trtcode), !is.na(class), nzchar(class)) %>%
      group_by(trtcode) %>%
      slice(1) %>%
      ungroup()

    dat <- dat %>%
      mutate(treatment = suppressWarnings(as.numeric(treatment))) %>%
      left_join(trt_map_use %>% rename(treatment = trtcode), by = "treatment")
  }

  req <- c("studyid", "class", "n", "mean_change", "sd_change")
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) stop("Input data missing required columns: ", paste(miss, collapse = ", "))

  if (!"sd_baseline" %in% names(dat)) dat$sd_baseline <- NA_real_
  if (!"source" %in% names(dat)) dat$source <- "unknown"

  dat <- dat %>%
    mutate(
      studyid = as.character(studyid),
      class = as.character(class),
      n = suppressWarnings(as.numeric(n)),
      mean_change = suppressWarnings(as.numeric(mean_change)),
      sd_change = suppressWarnings(as.numeric(sd_change)),
      sd_baseline = suppressWarnings(as.numeric(sd_baseline)),
      source = as.character(source)
    ) %>%
    filter(!is.na(studyid), !is.na(class), nzchar(class))

  dat$treatment_node <- infer_treatment_node(dat)

  control_priority <- c("Placebo", "Attention placebo", "No treatment", "Waitlist", "TAU")

  control_pick <- dat %>%
    filter(class %in% control_priority, !is.na(n), n > 1, !is.na(sd_change), sd_change > 0) %>%
    mutate(
      control_order = match(class, control_priority),
      usable_baseline = !is.na(sd_baseline) & sd_baseline > 0
    ) %>%
    arrange(studyid, control_order, desc(usable_baseline), desc(n)) %>%
    group_by(studyid) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      control_sd_ref = dplyr::case_when(
        usable_baseline ~ sd_baseline,
        fallback_mode == "sd_change" ~ sd_change,
        TRUE ~ NA_real_
      ),
      control_sd_source = dplyr::case_when(
        usable_baseline ~ "sd_baseline",
        fallback_mode == "sd_change" ~ "sd_change_fallback",
        TRUE ~ "missing"
      )
    ) %>%
    select(studyid, control_class = class, control_sd_ref, control_sd_source)

  dat_std <- dat %>%
    left_join(control_pick, by = "studyid") %>%
    mutate(
      drop_reason = case_when(
        is.na(control_class) ~ "no_control_arm_found",
        is.na(control_sd_ref) ~ "control_sd_not_available",
        is.na(n) | n <= 1 ~ "invalid_n",
        is.na(mean_change) ~ "missing_mean_change",
        is.na(sd_change) | sd_change <= 0 ~ "invalid_sd_change",
        TRUE ~ "kept"
      )
    )

  write_csv(
    dat_std %>% count(drop_reason, sort = TRUE),
    file.path(out_dir, "drop_reason_counts_control_arm_baseline_sd.csv")
  )
  write_csv(
    dat_std %>% distinct(studyid, control_class, control_sd_ref, control_sd_source) %>% arrange(studyid),
    file.path(out_dir, "control_arm_selection_control_arm_baseline_sd.csv")
  )

  dat_std <- dat_std %>%
    filter(drop_reason == "kept") %>%
    mutate(
      mean_change = mean_change / control_sd_ref,
      sd_change = sd_change / control_sd_ref
    )

  variants <- list()
  if (run_mode == "multi") {
    variants <- list(
      list(name = "primary_hybrid_any_ad_any_psychotherapy", mode = "hybrid", bridge = bridge_sets$any_ad_any_psychotherapy, lcc = FALSE),
      list(name = "sensitivity_class_only", mode = "class", bridge = bridge_sets$any_ad_any_psychotherapy, lcc = FALSE),
      list(name = "sensitivity_treatment_lcc", mode = "treatment", bridge = bridge_sets$any_ad_any_psychotherapy, lcc = TRUE),
      list(name = "sensitivity_hybrid_any_ad_only", mode = "hybrid", bridge = bridge_sets$any_ad_only, lcc = FALSE)
    )
  } else {
    variants <- list(list(
      name = paste0("single_", run_mode, "_", bridge_key),
      mode = run_mode,
      bridge = bridge_sets[[bridge_key]],
      lcc = FALSE
    ))
  }

  all_graph <- lapply(variants, function(v) {
    run_variant(
      dat_std = dat_std,
      variant_name = v$name,
      mode = v$mode,
      bridge_classes = v$bridge,
      force_lcc = v$lcc,
      reference_group = reference_group,
      out_dir = out_dir
    )
  }) %>% bind_rows()

  write_csv(all_graph, file.path(out_dir, "analysis_variants_graph_summary.csv"))

  saveRDS(
    list(
      run_mode = run_mode,
      bridge_key = bridge_key,
      reference_group = reference_group,
      fallback_mode = fallback_mode,
      variants = variants,
      graph_summary = all_graph
    ),
    file.path(out_dir, "analysis_run_manifest.rds")
  )

  cat("Done.\n")
  cat("Wrote outputs to: ", out_dir, "\n", sep = "")
}

main()
