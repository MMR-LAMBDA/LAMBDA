#!/usr/bin/env Rscript
#
#================================================================================
# Fitness Analysis Pipeline
#
# Purpose:
# An integrated R analysis pipeline. It takes variant count outputs from Enrich2
# and performs a complete analysis to generate final, WT-normalized fitness
# scores, volcano plots, heatmaps, and correlation plots.
#
# Workflow:
# 1. (Merge)     Merge individual sample count files from Enrich2 output.
# 2. (Filter)    Filter and pad counts against a master list (possible_combinations.txt).
# 3. (Combine)   Combine synonymous (_sy) and wild-type (_wt) counts.
# 4. (Calculate) Calculate fitness scores (logFC) using the limma-voom workflow.
# 5. (Plot)      Generate Volcano Plots.
# 6. (Plot)      Generate a Heatmap of significant variants and a logFC
#                correlation scatter plot.
#
# Dependencies:
#   - R Environment
#   - R Packages: dplyr, readr, tidyr, purrr, tibble, edgeR, limma,
#                 ggplot2, ggrepel, patchwork, pheatmap
#   - 'count_variants/tsv/' directory (output from Enrich2)
#   - 'possible_combinations.txt' file (master variant list)
#
# Produces:
#   - 'fitness_score_results/' (directory with .tsv result tables)
#   - 'plots_output/' (directory with diagnostic plots, e.g., voom trend)
#   - 'normalized_counts_output/' (directory with normalized count tables)
#   - 'combined_volcano_plots_R1_R2_vs_Input.pdf/png' (Volcano Plots)
#   - 'plots_output/heatmap_significant_mutations.pdf' (Heatmap)
#   - 'plots_output/correlation_logFC_R1_vs_R2.pdf' (Correlation Plot)
#================================================================================

cat("Master analysis pipeline started.\n")

# --- 0. ENVIRONMENT SETUP ---

# 0.1. Load Libraries
# On first run, ensure all packages are installed:
# install.packages(c("dplyr", "readr", "tidyr", "purrr", "tibble", "ggplot2", "ggrepel", "patchwork", "pheatmap"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("edgeR", "limma"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
})

# 0.2. Define Key Parameters & File Paths

# --- Analysis Parameters ---
WT_IDENTIFIER <- "_wt" # Wild-type identifier for relative fitness normalization

# --- Plotting Parameters ---
VOLCANO_LOG_FC_THRESHOLD <- 5.0
VOLCANO_ADJ_PVAL_THRESHOLD <- 0.05
VOLCANO_TOP_N_LABEL <- 40
VOLCANO_PLOT_TITLE_BASE <- "Difference between"

# --- Input Files/Dirs ---
ENRICH2_COUNTS_DIR <- file.path("count_variants", "tsv")
COMBINATIONS_FILE <- "possible_combinations.txt"

# --- Intermediate Files (will be created and overwritten) ---
STEP_03_OUTPUT <- "merged_synonymous_counts.tsv"
STEP_04_OUTPUT <- "filted_merged_all_combinations.tsv"
STEP_05_OUTPUT <- "final_merged_with_wt_sy_combined.tsv" # Input for Step 4

# --- Final Output Dirs/Files ---
FITNESS_RESULTS_DIR <- "fitness_score_results"
PLOTS_DIR <- "plots_output"
NORMALIZED_COUNTS_DIR <- "normalized_counts_output"

FITNESS_FILE_R1 <- file.path(FITNESS_RESULTS_DIR, "fitness_scores_Round1_vs_Input.tsv")
FITNESS_FILE_R2 <- file.path(FITNESS_RESULTS_DIR, "fitness_scores_Round2_vs_Input.tsv")
FITNESS_FILE_R2_R1 <- file.path(FITNESS_RESULTS_DIR, "fitness_scores_Round2_vs_Round1.tsv")

VOLCANO_PLOT_PDF <- "combined_volcano_plots_R1_R2_vs_Input.pdf"
VOLCANO_PLOT_PNG <- "combined_volcano_plots_R1_R2_vs_Input.png"
HEATMAP_PLOT_PDF <- file.path(PLOTS_DIR, "heatmap_significant_mutations.pdf")
CORRELATION_PLOT_PDF <- file.path(PLOTS_DIR, "correlation_logFC_R1_vs_R2.pdf")

# --- Create Output Dirs ---
dir.create(FITNESS_RESULTS_DIR, showWarnings = FALSE)
dir.create(PLOTS_DIR, showWarnings = FALSE)
dir.create(NORMALIZED_COUNTS_DIR, showWarnings = FALSE)


# --- 1. MERGE ENRICH2 COUNTS ---
#================================================================================
cat(paste("\n--- Step 1: Merging Enrich2 counts from", ENRICH2_COUNTS_DIR, "---\n"))

tryCatch({
  lib_folders <- list.dirs(path = ENRICH2_COUNTS_DIR, full.names = TRUE, recursive = FALSE)
  lib_folders <- lib_folders[grepl("_lib$", lib_folders)]
  
  if (length(lib_folders) == 0) {
    stop("No subfolders ending in _lib were found in ", ENRICH2_COUNTS_DIR)
  }
  
  all_data_list <- list()
  for (folder_path in lib_folders) {
    file_path <- file.path(folder_path, "main_synonymous_counts.tsv")
    if (file.exists(file_path)) {
      tryCatch({
        data_temp <- readr::read_tsv(
          file_path,
          skip = 1, # Skip header
          col_names = c("mutation_info", "count"),
          col_types = readr::cols(
            mutation_info = readr::col_character(),
            count = readr::col_integer()
          ),
          lazy = FALSE,
          progress = FALSE
        )
        folder_name <- basename(folder_path)
        new_col_name <- sub("_lib$", "", folder_name)
        data_temp <- data_temp %>% rename(!!new_col_name := count)
        all_data_list[[new_col_name]] <- data_temp
      }, error = function(e) {
        warning("Failed to read or process file: ", file_path, "\nError: ", e$message)
      })
    } else {
      warning("File not found: ", file_path)
    }
  }
  
  if (length(all_data_list) == 0) stop("No 'main_synonymous_counts.tsv' files were successfully processed.")
  
  if (length(all_data_list) == 1) {
    merged_data <- all_data_list[[1]]
  } else {
    valid_data_frames <- Filter(is.data.frame, all_data_list)
    if (length(valid_data_frames) < 1) stop("No valid data frames to merge.")
    merged_data <- valid_data_frames %>%
      purrr::reduce(function(df1, df2) full_join(df1, df2, by = "mutation_info"))
  }
  
  if (!is.null(merged_data) && nrow(merged_data) > 0) {
    count_cols <- setdiff(names(merged_data), "mutation_info")
    if(length(count_cols) > 0) {
      merged_data <- merged_data %>%
        mutate(across(all_of(count_cols), ~tidyr::replace_na(., 0)))
    }
    readr::write_tsv(merged_data, STEP_03_OUTPUT)
    cat(paste("Step 1 complete. Merged counts saved to:", STEP_03_OUTPUT, "\n"))
  } else {
    stop("No data available to merge after processing.")
  }
}, error = function(e) {
  stop(paste("ERROR in Step 1 (Merge Counts):", e$message))
})


# --- 2. FILTER AND PAD COUNTS ---
#================================================================================
cat(paste("\n--- Step 2: Filtering counts against", COMBINATIONS_FILE, "---\n"))

tryCatch({
  if (!file.exists(STEP_03_OUTPUT)) stop(paste("Input file not found:", STEP_03_OUTPUT))
  if (!file.exists(COMBINATIONS_FILE)) stop(paste("Input file not found:", COMBINATIONS_FILE))
  
  merged_counts_df <- readr::read_tsv(STEP_03_OUTPUT, show_col_types = FALSE)
  if (!"mutation_info" %in% names(merged_counts_df)) stop(paste("'", STEP_03_OUTPUT, "' is missing 'mutation_info' column."))
  
  possible_combinations_df <- readr::read_tsv(COMBINATIONS_FILE, show_col_types = FALSE)
  if (!"mutation_info" %in% names(possible_combinations_df)) stop(paste("'", COMBINATIONS_FILE, "' is missing 'mutation_info' column."))
  
  possible_combinations_df <- possible_combinations_df %>%
    select(mutation_info) %>%
    mutate(mutation_info = as.character(mutation_info)) %>%
    distinct(mutation_info)
  
  final_merged_df <- dplyr::left_join(possible_combinations_df, merged_counts_df, by = "mutation_info")
  
  count_columns <- setdiff(names(final_merged_df), "mutation_info")
  if (length(count_columns) > 0) {
    final_merged_df <- final_merged_df %>%
      mutate(across(all_of(count_columns), ~tidyr::replace_na(., 0)))
  }
  
  readr::write_tsv(final_merged_df, STEP_04_OUTPUT)
  cat(paste("Step 2 complete. Filtered/padded counts saved to:", STEP_04_OUTPUT, "\n"))
  
}, error = function(e) {
  stop(paste("ERROR in Step 2 (Filter Counts):", e$message))
})


# --- 3. COMBINE _WT AND _SY ROWS ---
#================================================================================
cat("\n--- Step 3: Combining _wt and _sy (synonymous) rows ---\n")

tryCatch({
  if (!file.exists(STEP_04_OUTPUT)) stop(paste("Input file not found:", STEP_04_OUTPUT))
  
  df <- readr::read_tsv(STEP_04_OUTPUT, show_col_types = FALSE)
  if (!"mutation_info" %in% names(df)) stop(paste("'", STEP_04_OUTPUT, "' is missing 'mutation_info' column."))
  
  wt_row_exists <- "_wt" %in% df$mutation_info
  sy_row_exists <- "_sy" %in% df$mutation_info
  
  if (wt_row_exists && sy_row_exists) {
    cat("Found '_wt' and '_sy' rows. Merging...\n")
    
    idx_wt <- which(df$mutation_info == "_wt")[1]
    idx_sy <- which(df$mutation_info == "_sy")[1]
    
    numeric_columns <- setdiff(names(df), "mutation_info")
    
    for (col_name in numeric_columns) {
      if (is.numeric(df[[col_name]])) {
        num_value_wt <- ifelse(is.na(df[[col_name]][idx_wt]), 0, as.numeric(df[[col_name]][idx_wt]))
        num_value_sy <- ifelse(is.na(df[[col_name]][idx_sy]), 0, as.numeric(df[[col_name]][idx_sy]))
        df[idx_wt, col_name] <- num_value_wt + num_value_sy
      }
    }
    
    df <- df[-idx_sy, ]
    cat("'_sy' row data merged into '_wt' row and '_sy' row removed.\n")
    
  } else {
    cat("'_wt' and/or '_sy' rows not found. Skipping row merge.\n")
  }
  
  readr::write_tsv(df, STEP_05_OUTPUT)
  cat(paste("Step 3 complete. Final analysis count matrix saved to:", STEP_05_OUTPUT, "\n"))
  
}, error = function(e) {
  stop(paste("ERROR in Step 3 (Combine WT/SY):", e$message))
})


# --- 4. CALCULATE FITNESS SCORES (limma-voom) ---
#================================================================================
cat(paste("\n--- Step 4: Calculating fitness scores from", STEP_05_OUTPUT, "---\n"))

tryCatch({
  if (!file.exists(STEP_05_OUTPUT)) stop(paste("Input file not found:", STEP_05_OUTPUT))
  
  counts_df <- read.delim(STEP_05_OUTPUT, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (any(duplicated(counts_df$mutation_info))) {
    warning("Duplicate 'mutation_info' identifiers found. Keeping first occurrence.")
    counts_df <- counts_df[!duplicated(counts_df$mutation_info), ]
  }
  
  counts_matrix <- counts_df %>%
    tibble::column_to_rownames(var = "mutation_info") %>%
    as.matrix()
  mode(counts_matrix) <- "numeric"
  counts_matrix[is.na(counts_matrix)] <- 0
  
  if (!WT_IDENTIFIER %in% rownames(counts_matrix)) {
    stop(paste("CRITICAL ERROR: Specified WT_IDENTIFIER '", WT_IDENTIFIER,
               "' not found in the rownames of the count matrix (", STEP_05_OUTPUT, ").",
               " Cannot proceed with relative fitness calculation.", sep = ""))
  } else {
    cat(paste("WT identifier '", WT_IDENTIFIER, "' found. Proceeding with analysis.\n", sep = ""))
  }
  
  sample_ids <- colnames(counts_matrix)
  sample_metadata <- data.frame(SampleID = sample_ids, stringsAsFactors = FALSE) %>%
    mutate(
      BaseName = gsub("_\\d$", "", SampleID),
      Replicate = as.integer(gsub(".*_", "", SampleID)),
      Condition = factor(BaseName, levels = c("plasmid", "sel_c1", "sel_c2")),
      SelectionRound = case_when(
        grepl("plasmid", Condition) ~ 0,
        grepl("sel_c1", Condition) ~ 1,
        grepl("sel_c2", Condition) ~ 2,
        TRUE ~ NA_integer_
      )
    ) %>%
    select(SampleID, Condition, Replicate, SelectionRound)
  
  if (any(is.na(sample_metadata$Condition))) {
    stop("Failed to parse sample metadata from column names. Check naming convention (e.g., 'plasmid_1', 'sel_c1_1').")
  }
  
  genes_info_df <- data.frame(mutation_info = rownames(counts_matrix), row.names = rownames(counts_matrix))
  
  dge_obj <- DGEList(counts = counts_matrix, samples = sample_metadata, genes = genes_info_df)
  
  design_matrix <- model.matrix(~0 + Condition, data = dge_obj$samples)
  colnames(design_matrix) <- levels(dge_obj$samples$Condition)
  
  contrasts_matrix <- makeContrasts(
    sel_c1_vs_plasmid = sel_c1 - plasmid,
    sel_c2_vs_plasmid = sel_c2 - plasmid,
    sel_c2_vs_sel_c1 = sel_c2 - sel_c1,
    levels = design_matrix
  )
  
  cat("Filtering low-expression variants...\n")
  keep <- filterByExpr(dge_obj, design = design_matrix)
  dge_obj_filtered <- dge_obj[keep, , keep.lib.sizes = FALSE]
  cat(paste(nrow(dge_obj_filtered), "variants kept out of", nrow(dge_obj), "\n"))
  
  cat("Calculating normalization factors (TMM)...\n")
  dge_obj_norm <- calcNormFactors(dge_obj_filtered, method = "TMM")
  
  cat("Running voom transformation...\n")
  voom_plot_file <- file.path(PLOTS_DIR, "voom_mean_variance_trend.pdf")
  pdf(voom_plot_file, width = 7, height = 6)
  voom_data <- voom(dge_obj_norm, design_matrix, plot = TRUE)
  invisible(dev.off())
  cat(paste("Voom mean-variance plot saved to:", voom_plot_file, "\n"))
  
  # --- Output normalized counts ---
  normalized_cpm_values <- cpm(dge_obj_norm, log = FALSE)
  normalized_cpm_df <- as.data.frame(normalized_cpm_values) %>%
    tibble::rownames_to_column(var = "mutation_info")
  
  log2_normalized_cpm_values <- voom_data$E
  log2_normalized_cpm_df <- as.data.frame(log2_normalized_cpm_values) %>%
    tibble::rownames_to_column(var = "mutation_info")
  
  write.table(normalized_cpm_df, file.path(NORMALIZED_COUNTS_DIR, "normalized_counts_per_sample_cpm.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(log2_normalized_cpm_df, file.path(NORMALIZED_COUNTS_DIR, "normalized_counts_per_sample_log2cpm.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("Normalized counts (CPM, log2-CPM) saved to:", NORMALIZED_COUNTS_DIR, "\n"))
  
  cat("Fitting linear models (limma)...\n")
  fit <- lmFit(voom_data, design_matrix)
  fit_contrasts <- contrasts.fit(fit, contrasts_matrix)
  fit_contrasts_bayes <- eBayes(fit_contrasts)
  
  cat("Extracting results and calculating relative fitness...\n")
  results_sel_c1_vs_plasmid <- topTable(fit_contrasts_bayes, coef = "sel_c1_vs_plasmid", number = Inf, sort.by = "P")
  results_sel_c2_vs_plasmid <- topTable(fit_contrasts_bayes, coef = "sel_c2_vs_plasmid", number = Inf, sort.by = "P")
  results_sel_c2_vs_sel_c1 <- topTable(fit_contrasts_bayes, coef = "sel_c2_vs_sel_c1", number = Inf, sort.by = "P")
  
  # --- Helper function for relative fitness ---
  calculate_relative_fitness <- function(results_table, wt_identifier_str) {
    if (!"mutation_info" %in% colnames(results_table)) {
      warning("'mutation_info' column not found.")
      return(results_table %>% mutate(Relative_logFC = NA_real_))
    }
    if (!wt_identifier_str %in% results_table$mutation_info) {
      warning(paste("WT identifier '", wt_identifier_str, "' not found in results table."))
      return(results_table %>% mutate(Relative_logFC = NA_real_))
    }
    lfc_wt_row <- results_table[results_table$mutation_info == wt_identifier_str, ]
    lfc_wt <- lfc_wt_row$logFC[1]
    if (is.na(lfc_wt)) {
      warning(paste("LFC for WT '", wt_identifier_str, "' is NA."))
      return(results_table %>% mutate(Relative_logFC = NA_real_))
    }
    results_table <- results_table %>%
      mutate(Relative_logFC = logFC - lfc_wt)
    return(results_table)
  }
  
  results_sel_c1_vs_plasmid_relative <- calculate_relative_fitness(results_sel_c1_vs_plasmid, WT_IDENTIFIER)
  results_sel_c2_vs_plasmid_relative <- calculate_relative_fitness(results_sel_c2_vs_plasmid, WT_IDENTIFIER)
  
  # --- Save final result tables ---
  write.table(results_sel_c1_vs_plasmid_relative, FITNESS_FILE_R1, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(results_sel_c2_vs_plasmid_relative, FITNESS_FILE_R2, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(results_sel_c2_vs_sel_c1, FITNESS_FILE_R2_R1, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat(paste("Step 4 complete. Fitness score tables saved to:", FITNESS_RESULTS_DIR, "\n"))
  
  # --- Generate R1 vs R2 fitness scatter plot (diagnostic) ---
  fitness_comparison_df <- inner_join(
    results_sel_c1_vs_plasmid_relative %>% select(mutation_info, Relative_Fitness_R1 = Relative_logFC, FDR_R1 = adj.P.Val),
    results_sel_c2_vs_plasmid_relative %>% select(mutation_info, Relative_Fitness_R2 = Relative_logFC, FDR_R2 = adj.P.Val),
    by = "mutation_info"
  ) %>% filter(!is.na(Relative_Fitness_R1) & !is.na(Relative_Fitness_R2))
  
  if(nrow(fitness_comparison_df) > 0) {
    fitness_scatter_plot <- ggplot(fitness_comparison_df, aes(x = Relative_Fitness_R1, y = Relative_Fitness_R2)) +
      geom_point(alpha = 0.5, aes(color = pmin(FDR_R1, FDR_R2))) +
      scale_color_viridis_c(name = "Min. FDR", option = "C") +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      geom_smooth(method = "lm", se = FALSE, color = "blue", formula = 'y ~ x') +
      labs(
        title = "Relative Fitness Score Comparison (WT-Normalized)",
        x = "Relative Fitness: Round 1 vs Input",
        y = "Relative Fitness: Round 2 vs Input"
      ) + theme_bw(base_size = 12)
    
    ggsave(file.path(PLOTS_DIR, "scatter_fitness_R1_vs_R2_relative.pdf"), fitness_scatter_plot, width = 8, height = 7)
  }
  
}, error = function(e) {
  stop(paste("ERROR in Step 4 (Calculate Fitness):", e$message))
})


# --- 5. GENERATE VOLCANO PLOTS ---
#================================================================================
cat("\n--- Step 5: Generating Volcano Plots ---\n")

tryCatch({
  # --- 5.1 Helper Function: Amino Acid Conversion ---
  convert_amino_acid_volcano <- function(seq) {
    seq <- gsub("Ala", "A", seq, ignore.case = TRUE); seq <- gsub("Arg", "R", seq, ignore.case = TRUE)
    seq <- gsub("Asn", "N", seq, ignore.case = TRUE); seq <- gsub("Asp", "D", seq, ignore.case = TRUE)
    seq <- gsub("Cys", "C", seq, ignore.case = TRUE); seq <- gsub("Gln", "Q", seq, ignore.case = TRUE)
    seq <- gsub("Glu", "E", seq, ignore.case = TRUE); seq <- gsub("Gly", "G", seq, ignore.case = TRUE)
    seq <- gsub("His", "H", seq, ignore.case = TRUE); seq <- gsub("Ile", "I", seq, ignore.case = TRUE)
    seq <- gsub("Leu", "L", seq, ignore.case = TRUE); seq <- gsub("Lys", "K", seq, ignore.case = TRUE)
    seq <- gsub("Met", "M", seq, ignore.case = TRUE); seq <- gsub("Phe", "F", seq, ignore.case = TRUE)
    seq <- gsub("Pro", "P", seq, ignore.case = TRUE); seq <- gsub("Ser", "S", seq, ignore.case = TRUE)
    seq <- gsub("Thr", "T", seq, ignore.case = TRUE); seq <- gsub("Trp", "W", seq, ignore.case = TRUE)
    seq <- gsub("Tyr", "Y", seq, ignore.case = TRUE); seq <- gsub("Val", "V", seq, ignore.case = TRUE)
    seq <- gsub("Ter", "*", seq, ignore.case = TRUE)
    return(seq)
  }
  
  # --- 5.2 Core Function: Generate Single Volcano Plot ---
  generate_volcano_plot <- function(input_file_path,
                                    plot_title_suffix = "",
                                    logFC_thresh = VOLCANO_LOG_FC_THRESHOLD,
                                    adj_pval_thresh = VOLCANO_ADJ_PVAL_THRESHOLD,
                                    top_n_genes = VOLCANO_TOP_N_LABEL) {
    
    cat(paste("  Generating volcano plot for:", basename(input_file_path), "\n"))
    if (!file.exists(input_file_path)) {
      cat(paste("  ERROR: Input file not found:", input_file_path, "\n")); return(NULL)
    }
    
    edgeRes <- read.delim(input_file_path, header = TRUE, stringsAsFactors = FALSE)
    if (!"mutation_info" %in% names(edgeRes)) {
      cat(paste("  ERROR:", basename(input_file_path), "is missing 'mutation_info' column.\n")); return(NULL)
    }
    
    edgeRes <- na.omit(edgeRes)
    if (nrow(edgeRes) == 0) {
      cat(paste("  ERROR: No data remaining after na.omit in", basename(input_file_path), "\n")); return(NULL)
    }
    
    original_row_names <- edgeRes$mutation_info
    new_row_names <- sapply(original_row_names, function(row_name) {
      row_name <- gsub("p\\.Phe", "", row_name, ignore.case = TRUE)
      row_name <- gsub("p\\.Glu", "", row_name, ignore.case = TRUE)
      row_name <- gsub("p\\.", "", row_name)
      row_name <- gsub(", ", "/", row_name)
      row_name <- convert_amino_acid_volcano(row_name)
      return(row_name)
    }, USE.NAMES = FALSE)
    edgeRes$gene_name <- new_row_names
    
    required_cols <- c("logFC", "adj.P.Val", "AveExpr")
    if (!all(required_cols %in% colnames(edgeRes))) {
      cat(paste("  ERROR: Dataframe missing required columns:", paste(setdiff(required_cols, colnames(edgeRes)), collapse = ", "), "\n")); return(NULL)
    }
    
    edgeRes_clean <- edgeRes %>%
      mutate(across(all_of(required_cols), ~suppressWarnings(as.numeric(as.character(.))))) %>%
      filter(is.finite(logFC) & is.finite(adj.P.Val) & is.finite(AveExpr) & adj.P.Val > 0)
    
    if (nrow(edgeRes_clean) == 0) {
      cat(paste("  ERROR: No finite data remaining for plotting in", basename(input_file_path), "\n")); return(NULL)
    }
    
    edgeRes_clean$neg_log10_adj_P_Val <- -log10(edgeRes_clean$adj.P.Val)
    current_plot_title <- paste(VOLCANO_PLOT_TITLE_BASE, plot_title_suffix)
    
    volcano_plot_obj <- ggplot(data = edgeRes_clean, aes(x = logFC, y = neg_log10_adj_P_Val)) +
      geom_point(aes(fill = AveExpr, color = AveExpr), shape = 21, size = 2, alpha = 0.7, stroke = 0) +
      
      # Changed "AveExpr" to "Avg. log2(CPM)" for clarity
      scale_fill_viridis_c(option = "viridis", name = "Avg. log2(CPM)", direction = -1, guide = "none") +
      scale_color_viridis_c(
        option = "viridis", name = "Avg. log2(CPM)", direction = -1,
        guide = guide_colorbar(
          barheight = unit(0.8, "npc"), barwidth = unit(1.2, "lines"),
          title.position = "right", title.theme = element_text(angle = 90, hjust = 0.5, vjust = 0.5), title.hjust = 0.5
        )
      ) +
      
      geom_hline(yintercept = -log10(adj_pval_thresh), linetype = "dashed", color = "grey50", linewidth = 0.7) +
      geom_vline(xintercept = c(-logFC_thresh, logFC_thresh), linetype = "dashed", color = "grey50", linewidth = 0.7) +
      labs(
        title = current_plot_title,
        x = expression("log"[2]*" Fold Change (logFC)"),
        y = expression("-log"[10]*" Adjusted P-value")
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = rel(1.2), face = "bold"),
        axis.title = element_text(size = rel(1.0), face = "bold"),
        axis.text = element_text(size = rel(0.9), color = "black"),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
      )
    
    edgeRes_clean$significance_group <- "Not Significant"
    edgeRes_clean$significance_group[which(edgeRes_clean$adj.P.Val < adj_pval_thresh & edgeRes_clean$logFC > logFC_thresh)] <- "Up"
    edgeRes_clean$significance_group[which(edgeRes_clean$adj.P.Val < adj_pval_thresh & edgeRes_clean$logFC < -logFC_thresh)] <- "Down"
    
    genes_to_label_df <- edgeRes_clean %>%
      filter(significance_group != "Not Significant") %>%
      arrange(adj.P.Val, desc(abs(logFC))) %>%
      head(top_n_genes)
    
    if (nrow(genes_to_label_df) > 0) {
      volcano_plot_obj <- volcano_plot_obj +
        geom_text_repel(
          data = genes_to_label_df, aes(label = gene_name), size = 4, color = "black",
          box.padding = unit(0.40, "lines"), point.padding = unit(0.35, "lines"),
          segment.color = 'grey50', segment.alpha = 0.7, max.overlaps = 10
        )
    }
    return(volcano_plot_obj)
  }
  
  # --- 5.3 Execute Plotting ---
  plot_r2 <- generate_volcano_plot(input_file_path = FITNESS_FILE_R2, plot_title_suffix = "(R2 vs Input)")
  plot_r1 <- generate_volcano_plot(input_file_path = FITNESS_FILE_R1, plot_title_suffix = "(R1 vs Input)")
  
  if (!is.null(plot_r1) && !is.null(plot_r2)) {
    combined_plot <- plot_r1 + plot_r2 + 
      plot_layout(ncol = 2) +
      plot_annotation(title = "Comparison of R1 vs Input & R2 vs Input",
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
    
    ggsave(VOLCANO_PLOT_PDF, plot = combined_plot, width = 16, height = 8, units = "in", dpi = 300)
    ggsave(VOLCANO_PLOT_PNG, plot = combined_plot, width = 16, height = 8, units = "in", dpi = 300)
    cat(paste("Step 5 complete. Combined volcano plots saved to:\n", VOLCANO_PLOT_PDF, "\n", VOLCANO_PLOT_PNG, "\n"))
    
  } else {
    cat("Could not generate one or both volcano plots. Skipping final plot combination.\n")
    if(!is.null(plot_r1)) ggsave(sub(".pdf", "_R1_only.pdf", VOLCANO_PLOT_PDF), plot = plot_r1, width = 8, height = 8)
    if(!is.null(plot_r2)) ggsave(sub(".pdf", "_R2_only.pdf", VOLCANO_PLOT_PDF), plot = plot_r2, width = 8, height = 8)
  }
  
}, error = function(e) {
  stop(paste("ERROR in Step 5 (Volcano Plots):", e$message))
})


# --- 6. GENERATE HEATMAP AND CORRELATION PLOT ---
#================================================================================
cat("\n--- Step 6: Generating Heatmap and Correlation Plot ---\n")

tryCatch({
  
  # --- 6.1 Read and prepare data ---
  cat("  Reading data for heatmap/correlation analysis...\n")
  fitness_r1 <- read_tsv(FITNESS_FILE_R1, show_col_types = FALSE)
  fitness_r2 <- read_tsv(FITNESS_FILE_R2, show_col_types = FALSE)
  norm_counts_file <- file.path(NORMALIZED_COUNTS_DIR, "normalized_counts_per_sample_log2cpm.tsv")
  norm_counts <- read_tsv(norm_counts_file, show_col_types = FALSE)
  
  cols_to_rename_r1 <- setdiff(colnames(fitness_r1), "mutation_info")
  colnames(fitness_r1)[colnames(fitness_r1) %in% cols_to_rename_r1] <- paste0(cols_to_rename_r1, "_R1")
  cols_to_rename_r2 <- setdiff(colnames(fitness_r2), "mutation_info")
  colnames(fitness_r2)[colnames(fitness_r2) %in% cols_to_rename_r2] <- paste0(cols_to_rename_r2, "_R2")
  
  merged_data <- norm_counts %>%
    full_join(fitness_r1, by = "mutation_info") %>%
    full_join(fitness_r2, by = "mutation_info")
  
  # --- 6.2 Helper Function: Amino Acid Conversion (for heatmap labels) ---
  convert_amino_acid_heatmap <- function(seq) {
    seq <- gsub("Ala", "A", seq); seq <- gsub("Arg", "R", seq)
    seq <- gsub("Asn", "N", seq); seq <- gsub("Asp", "D", seq)
    seq <- gsub("Cys", "C", seq); seq <- gsub("Gln", "Q", seq)
    seq <- gsub("Glu", "E", seq); seq <- gsub("Gly", "G", seq)
    seq <- gsub("His", "H", seq); seq <- gsub("Ile", "I", seq)
    seq <- gsub("Leu", "L", seq); seq <- gsub("Lys", "K", seq)
    seq <- gsub("Met", "M", seq); seq <- gsub("Phe", "F", seq)
    seq <- gsub("Pro", "P", seq); seq <- gsub("Ser", "S", seq)
    seq <- gsub("Thr", "T", seq); seq <- gsub("Trp", "W", seq)
    seq <- gsub("Tyr", "Y", seq); seq <- gsub("Val", "V", seq)
    seq <- gsub("Ter", "*", seq)
    return(seq)
  }
  
  current_mutation_info_values <- merged_data$mutation_info
  new_mutation_info_values <- sapply(current_mutation_info_values, function(row_name) {
    if (is.na(row_name)) return(NA)
    row_name <- gsub("p\\.Phe", "", row_name)
    row_name <- gsub("p\\.Glu", "", row_name)
    row_name <- gsub("p\\.", "", row_name)
    row_name <- gsub(", ", "/", row_name)
    row_name <- convert_amino_acid_heatmap(row_name)
    return(row_name)
  }, USE.NAMES = FALSE)
  
  merged_data$mutation_info <- new_mutation_info_values
  
  # --- 6.3 Generate Heatmap ---
  cat("  Generating heatmap...\n")
  logFC_threshold_R2 <- 4.0
  adj_p_value_threshold_R2 <- 0.02
  ave_expr_threshold_R2 <- 0
  
  required_filter_cols <- c("logFC_R2", "adj.P.Val_R2", "AveExpr_R2")
  if (!all(required_filter_cols %in% colnames(merged_data))) {
    stop("ERROR: Required columns for heatmap filtering are missing.")
  }
  
  significant_mutations_df <- merged_data %>%
    filter(
      !is.na(logFC_R2) & !is.na(adj.P.Val_R2) & !is.na(AveExpr_R2),
      abs(logFC_R2) >= logFC_threshold_R2,
      adj.P.Val_R2 < adj_p_value_threshold_R2,
      AveExpr_R2 > ave_expr_threshold_R2
    )
  
  if(nrow(significant_mutations_df) == 0) {
    warning("WARNING: No mutations passed the filtering criteria for the heatmap. Skipping heatmap generation.")
  } else {
    cat(sprintf("  Found %d significant mutations for heatmap.\n", nrow(significant_mutations_df)))
    
    log2cpm_columns <- c("plasmid_1", "plasmid_2", "plasmid_3",
                         "sel_c1_1", "sel_c1_2", "sel_c1_3",
                         "sel_c2_1", "sel_c2_2", "sel_c2_3")
    
    heatmap_matrix <- significant_mutations_df %>%
      select(mutation_info, all_of(log2cpm_columns)) %>%
      distinct(mutation_info, .keep_all = TRUE) %>%
      column_to_rownames(var = "mutation_info") %>%
      mutate_all(as.numeric) %>%
      as.matrix()
    
    if(nrow(heatmap_matrix) > 0 && ncol(heatmap_matrix) > 0) {
      annotation_col_df <- data.frame(
        SampleGroup = factor(rep(c("Input_Plasmid", "Round1_Selection", "Round2_Selection"), each = 3))
      )
      rownames(annotation_col_df) <- colnames(heatmap_matrix)
      
      # Display only Top 20 (or fewer)
      top_n_heatmap <- min(20, nrow(heatmap_matrix))
      heatmap_matrix_display <- heatmap_matrix[order(abs(significant_mutations_df$logFC_R2[match(rownames(heatmap_matrix), significant_mutations_df$mutation_info)]), decreasing=TRUE)[1:top_n_heatmap], , drop = FALSE]
      
      show_row_names_bool <- nrow(heatmap_matrix_display) <= 30
      row_font_size <- ifelse(nrow(heatmap_matrix_display) <= 50, 10, 6)
      
      pheatmap_plot_object <- pheatmap(
        heatmap_matrix_display,
        main = sprintf("Heatmap of log2(CPM) for Top %d Significant Mutations", top_n_heatmap),
        scale = "row",
        annotation_col = annotation_col_df,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        cutree_rows = 2,
        clustering_distance_rows = "correlation",
        show_rownames = show_row_names_bool,
        show_colnames = TRUE,
        fontsize_row = row_font_size,
        fontsize_col = 10,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        border_color = "grey60",
        gaps_col = c(3, 6),
        legend = TRUE,
        treeheight_row = 0,
        silent = TRUE # Prevent pheatmap from plotting to the default R device
      )
      
      # Save to file
      pdf(HEATMAP_PLOT_PDF, width = 8, height = 7)
      grid::grid.draw(pheatmap_plot_object$gtable)
      invisible(dev.off())
      
      cat(paste("  Heatmap saved to:", HEATMAP_PLOT_PDF, "\n"))
      
    } else {
      warning("WARNING: Heatmap matrix is empty. Skipping heatmap generation.")
    }
  }
  
  # --- 6.4 Generate Correlation Scatter Plot ---
  cat("  Generating correlation scatter plot...\n")
  sig_threshold <- 0.05
  
  required_cols_updated <- c("mutation_info", "logFC_R1", "logFC_R2", "adj.P.Val_R1", "adj.P.Val_R2")
  if (!all(required_cols_updated %in% colnames(merged_data))) {
    stop("ERROR: Required columns for correlation analysis are missing from merged_data.")
  }
  
  correlation_plot_data_pre_filter <- merged_data %>%
    select(all_of(required_cols_updated)) %>%
    filter(complete.cases(.[, c("logFC_R1", "logFC_R2", "adj.P.Val_R1", "adj.P.Val_R2")]))
  
  # *** Retain Option 1 only (Filter by R2 significance) ***
  correlation_plot_data <- correlation_plot_data_pre_filter %>%
    filter(adj.P.Val_R2 < sig_threshold)
  cat(sprintf("  Correlation filter: Using only Round 2 significant (adj.P.Val_R2 < %.2f) data points.\n", sig_threshold))
  
  min_data_points <- 3
  if(nrow(correlation_plot_data) < min_data_points) {
    warning("WARNING: Not enough significant data points for correlation analysis. Skipping correlation plot.")
  } else {
    cat(sprintf("  Found %d data points for correlation analysis.\n", nrow(correlation_plot_data)))
    
    correlation_test_result <- with(correlation_plot_data, cor.test(logFC_R1, logFC_R2, method = "pearson"))
    pearson_r <- correlation_test_result$estimate
    cor_p_value <- correlation_test_result$p.value
    
    linear_model <- lm(logFC_R2 ~ logFC_R1, data = correlation_plot_data)
    r_squared <- summary(linear_model)$r.squared
    slope <- coef(linear_model)[2]
    
    label_pearson <- sprintf("Pearson's r = %.2f", pearson_r)
    label_p_value <- ifelse(cor_p_value < 0.001, "p < 0.001", sprintf("p = %.3f", cor_p_value))
    label_r_sq_text <- sprintf("R^2 == %.2f", r_squared) # String for parse=TRUE
    label_slope <- sprintf("Slope = %.2f", slope)
    
    x_range <- range(correlation_plot_data$logFC_R1, na.rm = TRUE)
    y_range <- range(correlation_plot_data$logFC_R2, na.rm = TRUE)
    x_diff <- if(diff(x_range) == 0) 1 else diff(x_range)
    y_diff <- if(diff(y_range) == 0) 1 else diff(y_range)
    
    ann_x_pos <- x_range[1] + 0.05 * x_diff
    ann_y_pos <- c(y_range[2] - 0.05 * y_diff, y_range[2] - 0.12 * y_diff,
                   y_range[2] - 0.19 * y_diff, y_range[2] - 0.26 * y_diff)
    
    correlation_scatterplot <- ggplot(correlation_plot_data, aes(x = logFC_R1, y = logFC_R2)) +
      geom_point(alpha = 0.6, shape = 16, color = "grey20", size = 2) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#0072B2", fill = "#56B4E9", linewidth = 0.8) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick", linewidth = 0.7) +
      labs(
        x = expression("log"[2]*" Fold Change (Round 1 vs Input)"),
        y = expression("log"[2]*" Fold Change (Round 2 vs Input)"),
        title = "Round 1 vs Round 2 logFC Correlation (Filtered Data)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10, color = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
      ) +
      annotate("text", x = ann_x_pos, y = ann_y_pos[1], label = label_pearson, hjust = 0, size = 3.8) +
      annotate("text", x = ann_x_pos, y = ann_y_pos[2], label = label_p_value, hjust = 0, size = 3.8) +
      annotate("text", x = ann_x_pos, y = ann_y_pos[3], label = label_r_sq_text, hjust = 0, size = 3.8, parse = TRUE) +
      annotate("text", x = ann_x_pos, y = ann_y_pos[4], label = label_slope, hjust = 0, size = 3.8)
    
    ggsave(CORRELATION_PLOT_PDF, plot = correlation_scatterplot, width = 7, height = 6)
    cat(paste("  Correlation scatter plot saved to:", CORRELATION_PLOT_PDF, "\n"))
  }
  
  cat("Step 6 complete.\n")
  
}, error = function(e) {
  stop(paste("ERROR in Step 6 (Heatmap/Correlation):", e$message))
})


cat("\n================================================================================\n")
cat("Master analysis pipeline (6 steps) completed successfully.\n")
cat("================================================================================\n")