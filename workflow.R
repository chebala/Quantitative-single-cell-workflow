
# ===============================================================================
# 1. SETUP AND PACKAGE LOADING
# ===============================================================================

# Clear workspace
rm(list = ls())

# Load required packages
required_packages <- c(
  "tidyverse",    # Data manipulation and visualization
  "ggplot2",      # Advanced plotting
  "dplyr",        # Data manipulation
  "tidyr",        # Data reshaping
  "readr",        # CSV reading
  "readxl",       # Excel reading
  "ggridges",     # Ridge plots
  "pheatmap",     # Heatmaps
  "ggpubr",       # Publication-ready plots
  "broom",        # Tidy statistical output
  "scales",       # Plot scaling
  "RColorBrewer", # Color palettes
  "effectsize",   # effect size bcl2 distribution (ADDED MISSING COMMA)
  "limma" ,       # limma package
  "lmerTest"      # linear mixed effect model
)

# 1. Install standard CRAN packages first (excluding limma)
cran_packages <- setdiff(required_packages, "limma")
missing_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(missing_cran) > 0) {
  install.packages(missing_cran, dependencies = TRUE)
}

# 2. Install limma via BiocManager (Therough Bioconductor packages)
if (!"limma" %in% installed.packages()[,"Package"]) {
  if (!"BiocManager" %in% installed.packages()[,"Package"]) {
    install.packages("BiocManager", dependencies = TRUE)
  }
  BiocManager::install("limma", update = FALSE, ask = FALSE)
}

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Set random seed for reproducibility
set.seed(42)
# Load ggrepel
library(ggrepel)

# ===============================================================================
# 2. CONFIGURATION
# ===============================================================================

# Define paths (UPDATE THESE TO YOUR DATA LOCATIONS)
data_dir <- "E:/QmIF/qupath/export"  # Directory containing CSV files
viability_file <- "N:/Master thesis/Combinatorial treatment/synergy_analysis/data/combinatorial_treatment_tidy_data.xlsx"  # Path to viability Excel file
output_dir <- "E:/QmIF/output"  # Output directory for plots and tables

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE)

# Define cell lines (define your cell lines here)
cell_lines <- c("Cell_line_1", "Cell_line_2", "Cell_line_3", "Cell_line_4")

# Define treatment conditions (define your treatment conditions here) 
#I have 6 conditions including 3 single-agent treatment and two combinaotrial treatment

treatments <- c(
  "Untreated",  # Make sure to have a condition that can be used as a baseline to determine the treatment induced expression shifts
  "Condition_1", 
  "Condition_2",
  "Condition_3",
  "Condition_4",  
  "Condition_5"  
) # You can include more conditions but make sure to include them in the downstream analysis



# Define markers (mention then markers that you want to analyze)
markers <- c("Marker_1", "Marker_2", "Marker_3", "Marker_4")

# Define theme (I used a standardized theme across all of my plots)
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      # White background
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Remove grid
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      
      # Axes
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      
      # Legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      
      # Titles
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size + 1, hjust = 0.5),
      
      # Strip (facet labels)
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

# Set default theme
theme_set(theme_publication())

# Define color palette

treatment_colors <- c(
  "Untreated"       = "#1F77B4", 
  "Condition_1" = "#FF7F0E",
  "Condition_2" = "#2CA02C",
  "Condition_3" = "#D62728",
  "Condition_4"  = "#9467BD", 
  "Condition_5"  = "#8C564B"  
)

cell_line_colors <- c(
  "Cell_line_1" = "#E41A1C",
  "Cell_line_2" = "#377EB8",
  "Cell_line_3" = "#4DAF4A",
  "Cell_line_4" = "#984EA3"
)

# ===============================================================================
# 3. DATA IMPORT AND PREPROCESSING
# ===============================================================================

#' Import QuPath CSV files
#' 
#' This function reads all CSV files from the data directory and combines them
#' into a single dataframe with proper metadata.
#'
#' @param data_dir Directory containing QuPath CSV exports
#' @return Combined dataframe with all single-cell measurements
#' Import QuPath CSV files with Fuzzy Auto-Matching
# In each step, it is good to have a stop the function if the required file is not found as it might give false results in the downstream analysis
import_qupath_data <- function(data_dir) {
  cat("Importing QuPath data...\n")
  
  csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    stop("No CSV files found in the specified directory!")
  }
  
  cat(sprintf("Found %d CSV files\n", length(csv_files)))
  
  data_list <- lapply(csv_files, function(file) {
    cat(sprintf("Reading: %s\n", basename(file)))
    df <- read_csv(file, show_col_types = FALSE)
    
    fname <- tools::file_path_sans_ext(basename(file))
    df$SlideID <- fname
    
    parts <- unlist(strsplit(fname, "_"))
    raw_cell <- parts[1]
    
    # 1. Standardize Cell Line text variants safely
    # Here I have taken 4 cell lines and tried to standardize the names before using them in the pipelines
    df$CellLine <- case_when(
      raw_cell %in% c("Cellline1", "Cell_line-1") ~ "Cell_line_1",
      raw_cell %in% c("Cellline2", "Cell_line-2") ~ "Cell_line_2",
      raw_cell %in% c("Cellline3", "Cell_line-3") ~ "Cell_line_3",
      raw_cell %in% c("Cellline4", "Cell_line-4") ~ "Cell_line_4",
      TRUE ~ raw_cell
    )
    
    # 2. Extract and automatically correct Treatment naming gaps
    raw_treat <- paste(parts[2:(length(parts)-1)], collapse = "_")
    
    # Clean formatting anomalies (spaces, extra underscores, plus signs)
    raw_treat_clean <- gsub(" ", "_", raw_treat)
    raw_treat_clean <- gsub("\\+", "_", raw_treat_clean)
    raw_treat_clean <- gsub("__", "_", raw_treat_clean)
    # Here I have used some condition to find them in the df (This is just a example you have to specify the condition names by yourself
    #Tip: Try to use unique names as it is easier to use grepl function 
    # Match strings flexibly to lock down exact factor levels
    df$Treatment <- case_when(
      tolower(raw_treat_clean) == "untreated" ~ "Untreated",
      
      grepl("Condition", tolower(raw_treat_clean)) & 
        grepl("1", tolower(raw_treat_clean)) & 
        !grepl("_1", tolower(raw_treat_clean)) ~ "Condition_1",
      
      grepl("Condition", tolower(raw_treat_clean)) & 
        grepl("2", tolower(raw_treat_clean)) & 
        !grepl("_2", tolower(raw_treat_clean)) ~ "Condition_2",
      
      grepl("Condition", tolower(raw_treat_clean)) & 
        grepl("3", tolower(raw_treat_clean)) & 
        !grepl("_3", tolower(raw_treat_clean)) ~ "Condition_3",
      
      grepl("Cond", tolower(raw_treat_clean)) & grepl("Condition", tolower(raw_treat_clean)) & 
        grepl("4", tolower(raw_treat_clean)) & grepl("_4", tolower(raw_treat_clean)) ~ "Condition_4",
      
      grepl("Cond", tolower(raw_treat_clean)) & grepl("Condition", tolower(raw_treat_clean)) & 
        grepl("5", tolower(raw_treat_clean)) & grepl("_5", tolower(raw_treat_clean)) ~ "Condition_5",
      
      TRUE ~ raw_treat # Fallback if something doesn't match
    )
    
    return(df)
  })
  
  combined_data <- bind_rows(data_list)
  cat(sprintf("Total cells imported: %d\n", nrow(combined_data)))
  return(combined_data)
}

#' Clean and standardize column names
#' 
#' Standardizes QuPath column names to match expected format
#'
#' @param df Raw dataframe from QuPath
#' @return Dataframe with standardized column names
clean_column_names <- function(df) {
  cat("Standardizing column names...\n")
  
  df <- df %>%
    rename_with(~ gsub("Nucleus: ", "", .x)) %>%
    rename_with(~ gsub("Cell: ", "Cell_", .x)) %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    rename_with(~ gsub("Âµm", "um", .x)) %>%
    rename_with(~ gsub("\\^2", "2", .x))
  
  return(df)
}

#' Extract marker intensities
#' 
#' Extracts mean/median intensity values for each marker from QuPath data
#'
#' @param df Cleaned QuPath dataframe
#' @param markers Vector of marker names
#' @return Dataframe with extracted marker intensities
extract_marker_intensities <- function(df, markers) {
  cat("Extracting marker intensities...\n")
  # I have taken median instead of mean for the intensity value of each marker in specific compartments, you have to know your marker and you should know where your markers were expected to be found (nucleus/cytoplasm/cell membrane)
  # Expanded mapping list trying every possible QuPath capitalization and format
  marker_1_cols  <- c("Nucleus: Marker_1: Median", "Marker_1:_Median", "__Marker_1__Median", "Marker_1_Median", "Marker_1", "Nucleus: Marker_1 mean", "Marker_1_mean")
  marker_2_cols <- c("Nucleus: Marker_2: Median", "Nucleus_Marker_2:_Median", "Nucleus__Marker_2__Median", "Nucleus_Marker_2_Median", "Marker_2", "Nucleus: Marker_2 mean", "Marker_2_mean","Marker_2:_Median")
  marker_3_cols <- c("Cell: Marker_3: Median", "Cell_Marker_3:_Median", "Cell__Marker_3_Median", "Cell_Marker_3_Median", "Marker_3", "Cell: Marker_3 mean", "Marker_3_mean")
  marker_4_cols   <- c("Nucleus: Marker_4: Median", "Marker_4:_Median", "Nucleus__Marker_4__Median", "Marker_4_Median", "Marker_4", "Nucleus: Marker_4 mean", "Marker_4_mean")
  
  # Find which column actually exists
  actual_marker_1  <- intersect(marker_1_cols, colnames(df))
  actual_marker_2 <- intersect(marker_2_cols, colnames(df))
  actual_marker_3 <- intersect(marker_3_cols, colnames(df))
  actual_marker_4  <- intersect(marker_4_cols, colnames(df))
  
  # Build the renaming named vector dynamically
  rename_list <- c()
  if (length(actual_marker_1) > 0)  rename_list["Marker_1"]  <- actual_marker_1[1]
  if (length(actual_marker_2) > 0) rename_list["Marker_2"] <- actual_marker_2[1]
  if (length(actual_marker_3) > 0) rename_list["Marker_3"] <- actual_marker_3[1]
  if (length(actual_marker_4) > 0)  rename_list["Marker_4"]  <- actual_marker_4[1]
  
  if (length(rename_list) > 0) {
    df <- df %>% rename(any_of(rename_list))
  }
  
  # Final safety loop to check if they now exist
  for (marker in markers) {
    if (!marker %in% colnames(df)) {
      warning(sprintf("CRITICAL WARNING: Downstream script will fail. Could not find column for: %s", marker))
    }
  }
  
  return(df)
}

#' Filter low-quality cells
#' 
#' Removes cells based on quality criteria
#'
#' @param df Dataframe with single-cell data
#' @return Filtered dataframe
filter_low_quality_cells <- function(df) {
  cat("Filtering low-quality cells...\n")
  
  n_initial <- nrow(df)
  
  # Filter 1: Remove cells with missing marker values
  df <- df %>%
    filter(!is.na(Marker_1) & !is.na(Marker_2) & !is.na(Marker_3) & !is.na(Marker_4))
  
  # Filter 2: Remove cells with extreme area values (potential debris or clumps)
  if ("Nucleus_Area_um2" %in% colnames(df)) {
    df <- df %>%
      filter(Nucleus_Area_um2 > quantile(Nucleus_Area_um2, 0.01, na.rm = TRUE) &
               Nucleus_Area_um2 < quantile(Nucleus_Area_um2, 0.99, na.rm = TRUE))
  }
  
  n_final <- nrow(df)
  cat(sprintf("Removed %d cells (%.1f%%)\n", 
              n_initial - n_final, 
              100 * (n_initial - n_final) / n_initial))
  cat(sprintf("Retained %d cells\n", n_final))
  
  return(df)
}

# ===============================================================================
# 4b. CELL-LINE-ISOLATED PRINCIPAL COMPONENT ANALYSIS (QC STEP)
# ===============================================================================

run_cell_line_specific_pca <- function(df, markers, title_prefix, file_suffix, scale_flag = TRUE, max_cells_per_line = 15000) {
  cat(sprintf("\n[PCA ENGINE] Initializing Cell-Line-Isolated PCA for: %s State\n", title_prefix))
  
  missing_cols <- markers[!(markers %in% colnames(df))]
  if(length(missing_cols) > 0) {
    stop(paste("PCA Failed: Required columns missing from data matrix:", paste(missing_cols, collapse=", ")))
  }
  
  clean_df <- df %>% drop_na(all_of(markers))
  unique_lines <- unique(clean_df$CellLine)
  
  pca_models    <- list()
  combined_scores  <- list()
  combined_loadings <- list()
  variance_summaries <- list()
  
  for (line in unique_lines) {
    cat(sprintf(" Processing independent PCA matrix for Cell Line: %s...\n", line))
    
    line_df <- clean_df %>% filter(CellLine == line)
    
    if (nrow(line_df) > max_cells_per_line) {
      cat(sprintf("   Subsampling %d cells from %d total rows for %s\n", max_cells_per_line, nrow(line_df), line))
      line_df <- line_df %>% sample_n(max_cells_per_line)
    }
    
    pca_matrix <- line_df %>% select(all_of(markers)) %>% as.matrix()
    pca_res <- prcomp(pca_matrix, center = TRUE, scale. = scale_flag)
    pca_models[[line]] <- pca_res
    
    var_exp <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
    variance_summaries[[line]] <- var_exp
    
    scores_df <- as.data.frame(pca_res$x) %>%
      bind_cols(line_df %>% select(SlideID, CellLine, Treatment)) %>%
      mutate(PC1_Var = var_exp[1], PC2_Var = var_exp[2])
    combined_scores[[line]] <- scores_df
    
    load_df <- as.data.frame(pca_res$rotation[, 1:2]) %>%
      rownames_to_column("Marker") %>%
      mutate(
        CellLine = line,
        Standardized_Marker = gsub("log2_", "", Marker)
      )
    combined_loadings[[line]] <- load_df
    
    scores_df$Treatment <- factor(scores_df$Treatment, levels = treatments)
    
    # --- FIXED SPRINTF ARGS TO PREVENT VECTOR CRASH ---
    p_scatter <- ggplot(scores_df, aes(x = PC1, y = PC2, color = Treatment)) +
      geom_point(alpha = 0.25, size = 0.8) +
      stat_ellipse(level = 0.95, linewidth = 0.8) +
      scale_color_manual(values = treatment_colors) +
      labs(
        title = sprintf("%s PCA Landscape: %s", title_prefix, line),
        subtitle = "Independent Orthogonal Variance Projection with 95% Confidence Ellipses",
        x = sprintf("PC1 (%.1f%% Variance)", var_exp[1]),
        y = sprintf("PC2 (%.1f%% Variance)", var_exp[2]),
        color = "Treatment Condition"
      ) +
      theme_publication() +
      theme(legend.position = "bottom")
    
    clean_line_name <- gsub("-", "_", line)
    ggsave(
      filename = file.path(output_dir, sprintf("PCA_%s_scatter_%s.png", file_suffix, clean_line_name)),
      plot = p_scatter, width = 9, height = 7, dpi = 300
    )
    
    # ---------------------------------------------------------------------------
    # GENERATE PLOT PANEL 2: Isolated Biplot Arrow Map
    # ---------------------------------------------------------------------------
    load_df_scaled <- load_df %>%
      mutate(
        PC1_scaled = PC1 * max(abs(scores_df$PC1)) * 0.75,
        PC2_scaled = PC2 * max(abs(scores_df$PC2)) * 0.75
      )
    
    p_biplot <- ggplot() +
      geom_point(data = scores_df, aes(x = PC1, y = PC2, color = Treatment), alpha = 0.12, size = 0.7) +
      geom_segment(data = load_df_scaled, aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
                   arrow = arrow(length = unit(0.25, "cm"), type = "closed"), color = "black", linewidth = 1.1) +
      geom_text_repel(data = load_df_scaled, aes(x = PC1_scaled, y = PC2_scaled, label = Standardized_Marker),
                      size = 4.5, fontface = "bold", color = "black", box.padding = 0.4, max.overlaps = 30) +
      scale_color_manual(values = treatment_colors) +
      labs(
        title = sprintf("%s PCA Biplot: %s", title_prefix, line),
        subtitle = "Arrows project loading directions of markers across single cells",
        x = sprintf("PC1 (%.1f%% Variance)", var_exp[1]),
        y = sprintf("PC2 (%.1f%% Variance)", var_exp[2])
      ) +
      theme_publication() +
      theme(legend.position = "right")
    
    ggsave(
      filename = file.path(output_dir, sprintf("PCA_%s_biplot_%s.png", file_suffix, clean_line_name)),
      plot = p_biplot, width = 9, height = 7, dpi = 300
    )
  }
  
  final_scores   <- bind_rows(combined_scores)
  final_loadings <- bind_rows(combined_loadings)
  
  final_scores$Treatment <- factor(final_scores$Treatment, levels = treatments)
  final_scores$CellLine <- factor(final_scores$CellLine, levels = cell_lines)
  
  # ---------------------------------------------------------------------------
  # GENERATE PLOT PANEL 3: Standardized Faceted Global Overview
  # ---------------------------------------------------------------------------
  p_global_facet <- ggplot(final_scores, aes(x = PC1, y = PC2, color = Treatment)) +
    geom_point(alpha = 0.20, size = 0.7) +
    stat_ellipse(level = 0.95, linewidth = 0.7) +
    facet_wrap(~ CellLine, scales = "free", ncol = 2) +
    scale_color_manual(values = treatment_colors) +
    labs(
      title = sprintf("Global Comparative PCA View (%s Data)", title_prefix),
      subtitle = "Calculated independently within each phenotype background to protect treatment shifts",
      x = "Principal Component 1 (PC1 Axis)",
      y = "Principal Component 2 (PC2 Axis)"
    ) +
    theme_publication(base_size = 11) +
    theme(legend.position = "bottom")
  
  ggsave(
    filename = file.path(output_dir, sprintf("PCA_%s_global_faceted.png", file_suffix)),
    plot = p_global_facet, width = 13, height = 11, dpi = 300
  )
  
  return(list(
    models   = pca_models,
    scores   = final_scores,
    loadings = final_loadings,
    variance = variance_summaries
  ))
}

#=================================================================================
# 5 Plotting (Ridge)
#=================================================================================

#' Generate Distributional Ridge Plots
#' 
#' Generates 4 separate faceted plots (one per marker), stacking treatment conditions.
#'
#' @param df Transformed single-cell dataframe
#' @param markers Vector of target markers
#' @param output_dir Save path destination
generate_ridge_plots <- function(df, markers, output_dir) {
  cat("Generating single-cell distributional ridge plots...\n")
  
  # Ensure the factor strings match your treatment configuration layout perfectly
  df$Treatment <- factor(df$Treatment, levels = rev(treatments))
  
  for (marker in markers) {
    log_col <- paste0("log2_", marker)
    
    p <- ggplot(df, aes(x = .data[[log_col]], y = Treatment, fill = Treatment)) +
      geom_density_ridges(scale = 1.4, rel_min_height = 0.005, alpha = 0.75, bandwidth = NULL) +
      facet_wrap(~ CellLine, ncol = 2, scales = "free_x") +
      scale_fill_manual(values = treatment_colors) +
      labs(
        title = paste("Single-Cell Distribution profile:", marker),
        subtitle = "Faceted across Cell Lines (Unnormalized Shift Profiles)",
        x = "Log2 Intensity Value",
        y = "Treatment Protocol"
      ) +
      theme_publication() +
      theme(legend.position = "none") # Redundant due to categorical y-axis tracks
    
    # Save asset cleanly to drive
    file_out <- file.path(output_dir, paste0("ridge_plot_", marker, ".png"))
    ggsave(file_out, plot = p, width = 11, height = 8, dpi = 300)
    cat(sprintf("Saved ridge asset to: %s\n", basename(file_out)))
  }
}


#===============================================================================
# 6 Single-Cell Mixed effects model
#===============================================================================
#' Single-Cell Mixed Effects Statistical Engine with Robust Descriptive Effect Sizes
#' 
#' Natively handles individual cells as replicates while computing highly stable 
#' Cohen's d effect sizes directly from pooled single-cell variances to prevent optim crashes.
mixed_effect_model <- function(single_cell_df, markers, output_dir) {
  cat("Initializing Single-Cell Mixed-Effects Model Regression with Robust Effect Sizes...\n")
  library(lmerTest)
  
  clean_string_syntax <- function(x) {
    x <- gsub("-", "_", x)
    return(make.names(x))
  }
  
  local_cell_df <- single_cell_df %>%
    mutate(
      Treatment = clean_string_syntax(Treatment),
      CellLine  = clean_string_syntax(CellLine),
      SlideID   = clean_string_syntax(SlideID)
    )
  
  local_cell_df$Treatment <- factor(local_cell_df$Treatment, levels = clean_string_syntax(treatments))
  local_cell_df$CellLine  <- factor(local_cell_df$CellLine, levels = clean_string_syntax(cell_lines))
  
  log_markers <- paste0("log2_", markers)
  results_list <- list()
  
  # Loop through each marker at the pure single-cell level
  for (m in log_markers) {
    clean_marker_name <- gsub("log2_", "", m)
    cat(sprintf(" Fitting Single-Cell Mixed Model & Effect Sizes for: %s...\n", clean_marker_name))
    
    # 1. Fit the Linear Mixed Model
    model_formula <- as.formula(paste0(m, " ~ Treatment * CellLine + (1 | SlideID)"))
    lmm_fit <- lmer(model_formula, data = local_cell_df)
    
    # Extract standard regression coefficients and Satterthwaite p-values
    summary_stats <- as.data.frame(summary(lmm_fit)$coefficients) %>%
      rownames_to_column("Term") %>%
      mutate(Marker = clean_marker_name)
    
    # 2. BULLETPROOF DIRECT EFFECT SIZE CALCULATION (Bypasses t_to_d to prevent crashes)
    # Calculate global standard deviation for the baseline Untreated group
    untreated_sd <- local_cell_df %>%
      filter(Treatment == "Untreated") %>%
      pull(!!sym(m)) %>%
      sd(na.rm = TRUE)
    
    # If standard deviation is 0 or missing, default to global pool to protect against division by zero
    if(is.na(untreated_sd) || untreated_sd == 0) {
      untreated_sd <- sd(local_cell_df[[m]], na.rm = TRUE)
    }
    
    # Map raw Estimates directly to Cohen's d (d = mean_difference / pooled_sd)
    summary_stats <- summary_stats %>%
      mutate(
        Cohens_d = Estimate / untreated_sd,
        # Compute standard error of d based on sample size approximations
        d_SE = sqrt((1 / pmax(abs(Estimate), 1e-6)) + 
                      (Cohens_d^2 / (2 * nrow(local_cell_df)))))
    results_list[[m]] <- summary_stats
  }
  
  # 3. Consolidate single-cell regression outputs
  combined_stats <- bind_rows(results_list) %>%
    rename(
      Estimate  = `Estimate`,
      Std_Error = `Std. Error`,
      df_residual = `df`,
      t_value   = `t value`,
      P_Value   = `Pr(>|t|)`
    )
  
  # 4. Clean up and focus the output rows to show only your treatment adjustments
  final_mixed_results <- combined_stats %>%
    filter(grepl("Treatment", Term)) %>%
    mutate(
      adj_P_Val = p.adjust(P_Value, method = "BH"),
      Effect_Magnitude = case_when(
        abs(Cohens_d) < 0.2 ~ "Negligible",
        abs(Cohens_d) >= 0.2 & abs(Cohens_d) < 0.5 ~ "Small",
        abs(Cohens_d) >= 0.5 & abs(Cohens_d) < 0.8 ~ "Medium",
        abs(Cohens_d) >= 0.8 ~ "Large"
      ),
      Interpretation = case_when(
        !grepl(":", Term) ~ paste("Global drug shift inside reference cell line (Cell_line_1)"),
        grepl(":", Term)  ~ paste("Specific modification of drug response unique to this cell line background")
      )
    )
  
  # 5. Save the finalized single-cell replicate report to your folder
  write_csv(final_mixed_results, file.path(output_dir, "single_cell_mixed_model_effect_sizes.csv"))
  cat("Success! Mixed-effects models and Cohen's d effect sizes saved to output directory.\n")
  
  return(final_mixed_results)
}


#' Transform single-cell intensities
transform_single_cells <- function(df, markers) {
  cat("Applying log2(x + 1) transformation to single cells...\n")
  for (marker in markers) {
    df[[paste0("log2_", marker)]] <- log2(df[[marker]] + 1)
  }
  return(df)
}

#' Aggregate single cells to sample level
aggregate_to_samples <- function(df, markers) {
  cat("Aggregating single-cell data to sample-level medians...\n")
  sample_summary <- df %>%
    group_by(CellLine, Treatment, SlideID) %>%
    summarize(across(paste0("log2_", markers), ~ median(.x, na.rm = TRUE)), .groups = "drop")
  return(sample_summary)
}

#' Check Sample Normality Independently Within Each Cell Line
check_sample_normality <- function(sample_df, markers) {
  cat("\n=== SHAPIRO-WILK NORMALITY TEST (CELL-LINE ISOLATED) ===\n")
  
  unique_lines <- unique(sample_df$CellLine)
  log_markers <- paste0("log2_", markers)
  
  for (line in unique_lines) {
    cat(sprintf("\n--- Evaluating Normality for Cell Line: %s ---\n", line))
    
    # Isolate only the slides belonging to this specific cell line background
    line_data <- sample_df %>% filter(CellLine == line)
    
    for (marker in markers) {
      log_col <- paste0("log2_", marker)
      
      # Enforce a minimum sample threshold check to allow calculation
      if (nrow(line_data) >= 3 && length(unique(line_data[[log_col]])) >= 3) {
        sw_test <- shapiro.test(line_data[[log_col]])
        
        cat(sprintf("  Marker %-5s -> p-value: %.5f (%s)\n", 
                    marker, 
                    sw_test$p.value, 
                    if(sw_test$p.value > 0.05) "PASS (Normal)" else "FAIL (Skewed/Non-Normal)"))
      } else {
        cat(sprintf("  Marker %-5s -> Skipping (Not enough distinct slide replicates)\n", marker))
      }
    }
  }
  cat("========================================================\n\n")
}

#===========================================================================
#Import viability data##
#===========================================================================
#' Import viability data
#' 
#' Reads viability data from Excel file
#'
#' @param viability_file Path to Excel file
#' @return Dataframe with viability measurements
import_viability_data <- function(viability_file) {
  
  cat("Importing viability data...\n")
  
  # Read Excel file
  viability_data <- read_excel(viability_file)
  
  # Standardize column names by replacing spaces with underscores
  viability_data <- viability_data %>%
    rename_with(~ gsub(" ", "_", .x))
  
  cat(sprintf("Viability data: %d rows loaded\n", nrow(viability_data)))
  
  return(viability_data)
}
# ===============================================================================
# 7b. SINGLE-CELL NON-PARAMETRIC DISTRIBUTION ENGINE
# ===============================================================================

#' Analyze Single-Cell Population Shifts for All Markers
#' 
#' Loops through all markers and cell lines to execute Kruskal-Wallis 
#' and pairwise Wilcoxon tests with Rank-Biserial r effect sizes.
#'
#' @param single_cell_df Your transformed single-cell dataframe (results$single_cell_processed)
#' @param markers Vector of target markers (c("MYC", "PLK1", "BCL2", "P53"))
#' @param output_dir Destination directory for summary spreadsheets
analyze_all_marker_distributions <- function(single_cell_df, markers, output_dir) {
  cat("\n[DISTRIBUTION ENGINE] Running non-parametric single-cell shift metrics...\n")
  
  # Standardize level factors to enforce clean matching controls
  safe_treatments <- make.names(treatments)
  
  kw_master_list <- list()
  pairwise_master_list <- list()
  
  for (m in markers) {
    # Dynamically point to your log2-transformed column keys
    log_col <- paste0("log2_", m)
    
    cat(sprintf(" Processing single-cell population curves for: %s...\n", m))
    
    # 1. Compute Global Kruskal-Wallis per Cell Line for this specific marker
    kw_result <- single_cell_df %>%
      group_by(CellLine) %>%
      summarise(
        Marker = m,
        kruskal_H = kruskal.test(get(log_col) ~ Treatment)$statistic,
        p_value = kruskal.test(get(log_col) ~ Treatment)$p.value,
        df = kruskal.test(get(log_col) ~ Treatment)$parameter,
        .groups = "drop"
      )
    kw_master_list[[m]] <- kw_result
    
    # 2. Setup Pairwise Loop comparing every treated condition back to its own Untreated group
    groups_to_test <- single_cell_df %>%
      filter(Treatment != "Untreated") %>%
      distinct(CellLine, Treatment)
    
    results_list <- list()
    
    for(i in 1:nrow(groups_to_test)) {
      current_cell_line <- groups_to_test$CellLine[i]
      current_treatment <- groups_to_test$Treatment[i]
      
      tx_vector <- single_cell_df[[log_col]][
        single_cell_df$CellLine == current_cell_line & single_cell_df$Treatment == current_treatment
      ]
      
      ctrl_vector <- single_cell_df[[log_col]][
        single_cell_df$CellLine == current_cell_line & single_cell_df$Treatment == "Untreated"
      ]
      
      if(length(tx_vector) > 0 & length(ctrl_vector) > 0) {
        p_val <- wilcox.test(tx_vector, ctrl_vector)$p.value
        
        # Calculate Rank-Biserial r effect size via your imported effectsize library
        eff_r <- abs(effectsize::rank_biserial(tx_vector, ctrl_vector)$r)
        
        results_list[[i]] <- data.frame(
          Marker    = m,
          CellLine  = current_cell_line,
          Treatment = current_treatment,
          wilcox_p  = p_val,
          effect_size_r = eff_r
        )
      }
    }
    
    # Compile and apply false discovery adjustments to the tests for this marker
    pairwise_results <- bind_rows(results_list) %>%
      mutate(
        wilcox_p_adj = p.adjust(wilcox_p, method = "BH"),
        significance = case_when(
          wilcox_p_adj < 0.001 ~ "***",
          wilcox_p_adj < 0.01  ~ "**",
          wilcox_p_adj < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        effect_magnitude = case_when(
          effect_size_r >= 0.5 ~ "Large Biological Hit",
          effect_size_r >= 0.3 ~ "Moderate Shift",
          effect_size_r >= 0.1 ~ "Small/Weak Shift",
          TRUE                 ~ "Negligible (Sample Size Noise)"
        )
      )
    pairwise_master_list[[m]] <- pairwise_results
  }
  
  # 3. Consolidate everything into clean spreadsheets
  final_kw       <- bind_rows(kw_master_list)
  final_pairwise <- bind_rows(pairwise_master_list)
  
  write_csv(final_kw, file.path(output_dir, "single_cell_global_kruskal_wallis_results.csv"))
  write_csv(final_pairwise, file.path(output_dir, "single_cell_pairwise_wilcox_results.csv"))
  
  cat("Success! Single-cell population distribution datasets saved safely.\n")
  return(list(kruskal_wallis = final_kw, pairwise = final_pairwise))
}

# ===============================================================================
# 8. INTEGRATED HEATMAP ENGINE (SYNCHRONIZED & GAPPED)
# ===============================================================================

create_integrated_heatmap <- function(mixed_results, viability_data, output_dir) {
  cat("\n[HEATMAP ENGINE] Initializing Aligned Matrix Visualizer...\n")
  
  clean_text_format <- function(x) {
    x <- gsub("-", "_", x)
    x <- gsub(" ", "_", x)
    x <- gsub("\\+", "_", x)
    return(make.names(x))
  }
  
  # 1. Process viability sheet
  viability_summary <- viability_data %>%
    group_by(CellLine, Conditions) %>%
    summarise(mean_inhibition = mean(Inhibition, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Treatment = clean_text_format(Conditions),
      CellLine_Clean = case_when(
        CellLine %in% c("Cellline2", "Cell_line-2") ~ "Cell_line_2",
        CellLine %in% c("Cellline3", "Cell_line-3") ~ "Cell_line_3",
        CellLine %in% c("Cellline4", "Cell_line-4") ~ "Cell_line_4",
        TRUE ~ "Cell_line_1"
      )
    ) %>%
    select(CellLine = CellLine_Clean, Treatment, mean_inhibition)
  
  # 2. Reshape mixed model calculations 
  protein_wide <- mixed_results %>%
    mutate(
      CellLine = case_when(
        grepl("Cell_line_2", Term, ignore.case = TRUE) ~ "Cell_line_2",
        grepl("Cell_line_3", Term, ignore.case = TRUE) ~ "Cell_line_3",
        grepl("Cell_line_4", Term, ignore.case = TRUE) ~ "Cell_lin_4",
        TRUE ~ "Cell_line_1" 
      ),
      Clean_Treat = gsub("Treatment", "", Term),
      Treatment = gsub(":.*", "", Clean_Treat),
      Treatment = clean_text_format(Treatment)
    ) %>%
    select(CellLine, Treatment, Marker, Estimate) %>%
    pivot_wider(names_from = Marker, values_from = Estimate)
  
  # 3. Join data metrics
  heatmap_combined <- protein_wide %>%
    left_join(viability_summary, by = c("CellLine", "Treatment")) %>%
    mutate(
      Viability = 100 - coalesce(mean_inhibition, 0),
      Viability_Log2FC = log2(pmax(Viability, 0.1) / 100)
    ) %>%
    # FORCES rows to lock into your exact configured treatment level order sequence
    arrange(CellLine, factor(Treatment, levels = clean_text_format(treatments))) %>%
    mutate(RowName = paste(CellLine, Treatment, sep = "_"))
  
  # 4. Generate matrix core
  available_markers <- intersect(markers, colnames(heatmap_combined))
  heatmap_matrix <- heatmap_combined %>%
    select(RowName, all_of(available_markers)) %>%
    column_to_rownames("RowName") %>%
    as.matrix()
  
  heatmap_matrix[!is.finite(heatmap_matrix)] <- 0
  
  # 5. Synchronize annotation vector rows
  anno_data <- heatmap_combined %>%
    select(RowName, Viability_Log2FC) %>%
    column_to_rownames("RowName")
  anno_data$Viability_Log2FC[!is.finite(anno_data$Viability_Log2FC)] <- 0
  anno_data <- anno_data[rownames(heatmap_matrix), , drop = FALSE]
  
  anno_colors <- list(
    Viability_Log2FC = colorRampPalette(c("#4B0082", "#E6E6FA"))(100)
  )
  
  # 6. Render and save the high-resolution image asset
  png(filename = file.path(output_dir, "integrated_heatmap_grouped.png"), 
      width = 2600, height = 3000, res = 300)
  
  pheatmap(
    heatmap_matrix,
    annotation_row = anno_data,
    annotation_colors = anno_colors,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-0.5, 0.5, length.out = 101), 
    
    # === INCORPORATED ORDER FLAGS AS REQUESTED (CRASH PROTECTED) ===
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    gaps_row = which(diff(as.numeric(as.factor(heatmap_combined$CellLine))) != 0), 
    # ===============================================================
    
    main = "Protein (Mixed Model Estimate) & Viability Log2FC",
    border_color = "grey60",
    cellwidth = 45,
    cellheight = 18
  )
  
  dev.off()
  cat("Success! Cell line grouped heatmap saved cleanly with dividing lines.\n")
  return(heatmap_combined)
}




#' Automatically Deconvolute Interaction Terms into True Cell-Line Responses
#' 
#' @param mixed_results Dataframe generated by mixed_effect_model (results$final_mixed_results)
#' @param output_dir Folder path to save the simplified output spreadsheet
calculate_true_responses <- function(mixed_results, output_dir) {
  cat("\n[DECONVOLUTION ENGINE] Converting interaction terms to true cell-line responses...\n")
  
  # Ensure clean string matches
  safe_treatments <- make.names(treatments)
  safe_cell_lines <- make.names(cell_lines)
  # Reference line is always the first alphabetical item (scrubbed of hyphens)
  ref_cell_line <- safe_cell_lines[1] 
  alt_cell_lines <- setdiff(safe_cell_lines, ref_cell_line)
  
  # 1. Isolate the reference cell line (Z-138) baseline responses
  ref_data <- mixed_results %>%
    filter(!grepl(":", Term)) %>%
    mutate(
      CellLine = ref_cell_line,
      Treatment = gsub("Treatment", "", Term),
      True_Log2FC = Estimate,
      True_Cohens_d = Cohens_d
    ) %>%
    select(Marker, CellLine, Treatment, True_Log2FC, True_Cohens_d)
  
  # 2. Process alternative cell lines using the (Baseline + Modifier) formula
  alt_data_list <- list()
  
  for (cl in alt_cell_lines) {
    for (tr in setdiff(safe_treatments, "Untreated")) {
      
      # Target the specific modifier term row (e.g., "TreatmentCombo_IC50IC50:CellLineSU_DHL_4")
      modifier_term_pattern <- paste0("Treatment", tr, ":CellLine", cl)
      
      # Extract values for this specific marker combination loop
      for (m in unique(mixed_results$Marker)) {
        
        # Get baseline anchor value (Z-138)
        base_row <- ref_data %>% filter(Treatment == tr, Marker == m)
        base_est <- base_row$True_Log2FC
        base_d   <- base_row$True_Cohens_d
        
        # Get interaction modifier row
        mod_row <- mixed_results %>% filter(Term == modifier_term_pattern, Marker == m)
        
        if (nrow(mod_row) > 0) {
          # CORE FORMULA: Baseline + Modifier
          true_est <- base_est + mod_row$Estimate
          true_d   <- base_d + mod_row$Cohens_d
          
          alt_data_list[[paste(cl, tr, m, sep="_")]] <- data.frame(
            Marker = m,
            CellLine = cl,
            Treatment = tr,
            True_Log2FC = true_est,
            True_Cohens_d = true_d
          )
        }
      }
    }
  }
  
  # 3. Combine reference and alternative data blocks into one clean sheet
  deconvoluted_results <- bind_rows(ref_data, bind_rows(alt_data_list)) %>%
    mutate(
      # Re-evaluate descriptive magnitude tags based on the true localized Cohen's d values
      True_Effect_Magnitude = case_when(
        abs(True_Cohens_d) < 0.2 ~ "Negligible",
        abs(True_Cohens_d) >= 0.2 & abs(True_Cohens_d) < 0.5 ~ "Small",
        abs(True_Cohens_d) >= 0.5 & abs(True_Cohens_d) < 0.8 ~ "Medium",
        abs(True_Cohens_d) >= 0.8 ~ "Large"
      ),
      # Restore native clean cell line syntax strings for your final thesis table
      CellLine = gsub("_", "-", CellLine)
    ) %>%
    arrange(CellLine, factor(Treatment, levels = safe_treatments), Marker)
  
  # 4. Save clean data sheet directly to folder
  file_out <- file.path(output_dir, "true_cell_line_isolated_responses.csv")
  write_csv(deconvoluted_results, file_out)
  cat(sprintf("Success! Clean responses spreadsheet exported to: %s\n", basename(file_out)))
  
  return(deconvoluted_results)
}

# ===============================================================================
# 8. MAIN ANALYSIS PIPELINE
# ===============================================================================

#' Run complete analysis workflow
#'
#' @param data_dir Directory containing QuPath CSV files
#' @param viability_file Path to viability Excel file
#' @param output_dir Output directory
run_analysis <- function(data_dir, viability_file, output_dir) {
  
  cat("========================================\n")
  cat("mIF Single-Cell Analysis Pipeline\n")
  cat("========================================\n\n")
  
  # -----------------------------------------------------------------------------
  # STEP 1: Data Import and Preprocessing
  # -----------------------------------------------------------------------------
  cat("STEP 1: Data Import and Preprocessing\n")
  cat("========================================\n")
  
  # 1A. Import and basic column standardization
  raw_data     <- import_qupath_data(data_dir)
  cleaned_data <- clean_column_names(raw_data)
  
  # 1B. Capture Stage 0: Starting Cell Counts (Raw counts after column standardizing)
  stage0_counts <- cleaned_data %>%
    group_by(CellLine, Treatment, SlideID) %>%
    tally(name = "Stage0_Raw_Cells") %>%
    ungroup()
  
  # 1C. Filter Stage 1: Missing Marker Intensity Exclusions
  cat("\nRunning Filter Stage 1: Checking missing marker footprints...\n")
  extracted_data <- extract_marker_intensities(cleaned_data, markers)
  
  # Filter missing values manually here to capture the audit baseline checkpoint
  stage1_data <- extracted_data %>%
    filter(!is.na(Marker_1) & !is.na(Marker_2) & !is.na(Marker_3) & !is.na(Marker_4))
  
  stage1_counts <- stage1_data %>%
    group_by(CellLine, Treatment, SlideID) %>%
    tally(name = "Stage1_Complete_Markers") %>%
    ungroup()
  
  # 1D. Filter Stage 2: Extreme Morphological Debris Exclusions
  cat("\nRunning Filter Stage 2: Quality threshold filters...\n")
  
  # We reuse your structural logic directly to capture the exact final dataframe
  final_filtered_data <- filter_low_quality_cells(stage1_data)
  
  stage2_counts <- final_filtered_data %>%
    group_by(CellLine, Treatment, SlideID) %>%
    tally(name = "Stage2_Final_Retained") %>%
    ungroup()
  
  # 1E. Compile and Print/Save Filtration Audit Matrix
  filtration_audit <- stage0_counts %>%
    left_join(stage1_counts, by = c("CellLine", "Treatment", "SlideID")) %>%
    left_join(stage2_counts, by = c("CellLine", "Treatment", "SlideID")) %>%
    mutate(
      Total_Cells_Removed = Stage0_Raw_Cells - Stage2_Final_Retained,
      Pct_Cells_Retained  = round(100 * (Stage2_Final_Retained / Stage0_Raw_Cells), 1)
    )
  
  # Format factors for clean printing layout
  filtration_audit <- filtration_audit %>%
    arrange(CellLine, factor(Treatment, levels = treatments))
  
  # Print complete audit breakdown directly into the R Console
  cat("\n=== CELL FILTRATION AUDIT REPORT MATRIX ===\n")
  print(as.data.frame(filtration_audit), row.names = FALSE)
  cat("====================================================\n\n")
  
  # Save audit document to your designated output file directory
  audit_file_path <- file.path(output_dir, "cell_filtration_audit.csv")
  write_csv(filtration_audit, audit_file_path)
  cat(sprintf("Saved detailed filtration logging breakdown to: %s\n\n", basename(audit_file_path)))
  
  # 1F. Finish file intake with viability measurements
  viability_data <- import_viability_data(viability_file)
  
  # -----------------------------------------------------------------------------
  # STEP 2: Downstream Analysis Pipelines
  # -----------------------------------------------------------------------------
  cat("STEP 2: Downstream Analysis Pipelines\n")
  cat("========================================\n")
  
  # 2A. Apply Base-2 transformations smoothly to the passed cells
  single_cell_transformed <- transform_single_cells(final_filtered_data, markers)
  
  # 2B. Build the complete single-cell unnormalized distribution ridge series
  generate_ridge_plots(single_cell_transformed, markers, output_dir)
  
  # 2C. EXECUTE UPDATED PCA ENGINE
  log2_marker_columns <- paste0("log2_", markers)
  
  pca_outputs <- run_cell_line_specific_pca(
    df = single_cell_transformed,
    markers = log2_marker_columns,
    title_prefix = "Log2(X+1) Intensity",
    file_suffix = "log2_transformed",
    scale_flag = TRUE,
    max_cells_per_line = 15000
  )
  
  # 2D. Condense parameters down to sample median vectors for downstream Limma
  sample_level_data <- aggregate_to_samples(single_cell_transformed, markers)
  
  # 2E. Run sample normality diagnostic metrics checks
  check_sample_normality(sample_level_data, markers)
  
  # =======================================================

  
  # 2F. Calculate robust linear models using single cells as nested replicates
  mixed_effect_statistics <- mixed_effect_model(single_cell_transformed, markers, output_dir)
  
  # 2G. Build your cell line grouped integrated heatmap with dividing line blocks
  heatmap_package <- create_integrated_heatmap(mixed_effect_statistics, viability_data, output_dir)
  
  distribution_statistics <- analyze_all_marker_distributions(single_cell_transformed, markers, output_dir)
  # ===================================================================
  
  # Package and return data matrices cleanly to global workspace
  output_package <- list(
    single_cell_processed = single_cell_transformed,
    sample_medians        = sample_level_data,
    pca_results           = pca_outputs,
    final_mixed_results   = mixed_effect_statistics,
    distribution_stats    = distribution_statistics, 
    heatmap_data          = heatmap_package,
    filtration_report     = filtration_audit,
    viability             = viability_data
  )
}
# ===============================================================================
# 10. EXECUTE ANALYSIS
# ===============================================================================

# Run the complete analysis pipeline
# IMPORTANT: Update the paths in the CONFIGURATION section before running!

results <- run_analysis(data_dir, viability_file, output_dir)

# ===============================================================================
# END OF SCRIPT
# ===============================================================================
