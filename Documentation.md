*The workflow contains 10 stages*


*This workflow supports CSV files* 

## Stage 1: Setup and package loading

In this stage, all the required packages are loaded. 
Any packages that were not available will be installed automatically from CRAN repository.
limma package downloaded via Biomanager.

**Repository Setup & Package Dependency Management**

This initial pipeline module handles workspace standardization, automated dependency checking, multi-repository library installation (CRAN and Bioconductor), package loading, and pseudo-random seed synchronization for overall workflow reproducibility.

**Pipeline Setup Overview**

The package dependency infrastructure runs across three execution layers:
* **Workspace Environmental Initialization**: Clears out stale global memory constraints to guarantee isolated pipeline runs.
* **Dual-Channel Installer**: Automatically checks for missing packages and downloads them across Comprehensive R Archive Network (CRAN) and Bioconductor channels.
* **Namespace Registration**: Loads needed analysis frameworks into the working R session environment.

**Workspace Standardization**

To prevent data contamination from leftover data variables, the script triggers an automated workspace sweep:

```R
rm(list = ls())
```

This command flushes all previously saved dataframes, matrices, lists, and functions out of global memory. It guarantees that the pipeline runs inside an isolated environment, preventing silent naming collisions with variables left over from other scripts.

**Multi-Repository Dependency Management**

Multiplexed single-cell workflows require tools from distinct open-source software distributions. The engine coordinates a tailored dual-channel installation framework:

**1. Comprehensive R Archive Network (CRAN) Channel**
The pipeline screens your computer's local library path against a vector of core statistical and visualization dependencies:
* **Data Wrangling**: `tidyverse`, `dplyr`, `tidyr`, `readr`, `readxl`, `broom`
* **Plotting & Visual Aesthetics**: `ggplot2`, `ggridges`, `pheatmap`, `ggpubr`, `scales`, `RColorBrewer`, `ggrepel`
* **Modeling & Statistics**: `lmerTest`, `effectsize`

The script calculates exactly which CRAN packages are missing (`setdiff`), installs them along with their upstream dependencies (`dependencies = TRUE`), and skips already installed packages to save processing time.

**2. Bioconductor Channel**
The differential expression package **`limma`** is hosted exclusively through Bioconductor and cannot be installed via standard CRAN tools. 
* **The Fix**: The script identifies if `limma` is present. If missing, it installs the `BiocManager` interface from CRAN first, then automatically downloads `limma` without prompting the user (`update = FALSE, ask = FALSE`).


**Reproducibility & Numerical Integrity**

Single-cell algorithms—such as downsampling cell rows for Principal Component Analysis (PCA) or computing confidence ellipse borders—rely on pseudo-random number generators. 

```R
set.seed(42)
```

By locking in a static global random seed value (`42`), the script forces R to use identical numerical sequences across every execution loop. This step ensures that your downsampled PCA coordinates, biplot shapes, and calculations remain perfectly consistent and reproducible every time you or a collaborator runs the script.

## Stage 2: Configuration

Define the paths to CSV files that contains the measurements. Make sure to standardize the naming of files as the workflow is case sensitive
The output can be stored in a seperate directory (If the path does n't exist the workflow itself produce and save the output files.
With this setup all the input files and output files will be handled externally in either local drive or external hard drives. This allows to save R memory.
Here you can define the cell lines, treatment conditions, and markers.
All the plots were made by a standardized theme together with specific colors for treatment and cell lines.

**Pipeline Configuration & Visual Theming**

This pipeline module handles the central setup array for the workspace. It specifies storage file paths, establishes control thresholds, defines mapping vectors for experimental metadata, and sets up publication-ready visual formatting styles.

**Pipeline Configuration Overview**

The setup framework acts as a central control panel for the entire workflow:
* **Directory Provisioning**: Automates project folder generation to prevent downstream file export errors.
* **Factor Matrix Sequencing**: Defines the order of cell line and treatment vectors to ensure consistent layouts across tables and charts.
* **Visual Formatting Controls**: Installs a unified formatting template and standardizes theme elements across multiple visualization libraries.

 **Workspace Directory Mapping**

To set up the workspace, configure the following variables to point to your data files and storage locations:

* **`data_dir`**: The directory path containing your raw multiplexed single-cell segmentation tables exported from QuPath.
* **`viability_file`**: The file path pointing to your processed cellular synergy spreadsheet (`.xlsx`).
* **`output_dir`**: The destination directory where all summary charts, heatmaps, and audit statistics are saved.

*Note: The script runs an automated validation command (`dir.create(output_dir, showWarnings = FALSE)`) that checks your file system and automatically creates the target output folder if it does not exist.*

**Metadata Factors & Baseline Grounding**

The engine maps experimental variables into precise evaluation categories:

* **`cell_lines`**: Tracks the individual cellular backgrounds analyzed by the script (e.g., `Cell_line_1` through `Cell_line_4`).
* **`treatments`**: Sets up the mapping sequence for all six experimental arms, tracking single-agent variations alongside combinatorial treatments.
* **`markers`**: Defines the target multiplexed panel channels evaluated by downstream statistical components.

```text
CRITICAL DESIGN RULE:
The 'treatments' vector MUST always begin with your target control baseline ("Untreated"). 
The statistical engines use this baseline to calculate relative expression shifts, 
and the modeling blocks rely on it to derive custom Cohen's d effect sizes.
```

 **Graphic Styles & Publication Themes**

To prepare figures for publication, the script builds a custom layout engine (`theme_publication`) on top of the standard `ggplot2` framework. Running this engine clears default chart clutter and applies specific design parameters:

* **Canvas Isolation**: Applies a solid white workspace background (`element_rect`) and removes gridlines (`panel.grid.minor = element_blank()`) to keep plots clean.
* **Typography Control**: Enlarges axis titles (`face = "bold"`), centers chart headers (`hjust = 0.5`), and scales labels proportionally based on your chosen `base_size`.
* **Grid Formatting**: Adds clear borders and bold titles to multi-window facet labels (`strip.background`), making it easier to read charts across different cell line groups.

**Visual Identity Color Palettes**

The configuration sets up specific, high-contrast hex color strings to keep tracking labels consistent across independent plots, PCA layouts, and density profiles:

**1. Experimental Conditions** (`treatment_colors`)
Uses distinct colors to differentiate treatment setups, helping you visually separate control profiles from single-agent and combinatorial variations.
* `Untreated` (`#1F77B4` - Control Blue)
* `Condition_1` through `Condition_5` (High-contrast categorical shades)

**** 2. Biological Backgrounds**** (`cell_line_colors`)
Applies distinct color choices to cell lines to help you quickly identify individual biological backgrounds across multi-panel charts.

## Step 3: Data Import and Preprocessing

QuPath CSV files were combined and a master dataframe was created.
If there is no CSV file found in the specified path, the workflow will be halted.
The cell line names and treatment condition were standardized.
All the marker intensities werecaptured and checks whether the columns actually exists.

**Data Import & Preprocessing**

This pipeline module coordinates the file ingestion infrastructure. It automatically aggregates standalone QuPath multiplexed cell segmentation metrics, processes metadata tokens using fuzzy string-matching rules, standardizes encoding symbols, and maps subcellular marker expressions across cellular compartments.

**Pipeline Module Overview**
The data ingestion engine runs across three primary functional blocks:
* **`import_qupath_data`**: Discovers and scans raw `.csv` batch files from QuPath. It reads the source filename to extract `SlideID`, biological `CellLine`, and experimental `Treatment` properties on the fly.
* **`clean_column_names`**: Standardizes special characters, whitespace variations, and column prefix configurations exported by QuPath.
* **`extract_marker_intensities`**: Dynamically looks for and maps localized marker metrics (such as nuclear, cytoplasmic, or whole-cell segmentations) using an internal collection of known header variants.

**Fuzzy Auto-Matching & Extraction Rules**

To handle formatting differences or typo variations across multiple imaging runs, this module features built-in string cleaning rules:
** 1. Filename Processing Layout**
The script splits incoming filenames using an underscore delimiter (`_`). It assumes a standard naming layout: `[CellLine]_[Treatment_Conditions]_[Replicate_Details].csv`.
* **Cell Line Balancing**: Standardizes text variants into consistent factor labels (`Cell_line_1` to `Cell_line_4`).
* **Fuzzy Treatment Matching**: Uses regex substitution (`gsub`) to clean out space structures (` `), plus formatting (`+`), and double underscores (`__`). It then matches criteria keys using partial string parsing (`grepl`) to normalize categories (`Condition_1` through `Condition_5`).

**2. Targeted Subcellular Extraction**
Rather than relying on single, strict naming rules, the pipeline evaluates data columns across a predefined collection of potential QuPath column names to catch both Median and Mean measurements:
* **`Marker_1`, `Marker_2`, `Marker_4`**: Configured to look for nuclear expressions (`Nucleus: Marker: Median`, `Nucleus: Marker mean`).
* **`Marker_3`**: Configured to look for localized cytoplasmic or whole-cell expressions (`Cell: Marker_3: Median`, `Cell_Marker_3_Median`).

The engine runs a logical intersection check (`intersect`) to map the first valid matching column found directly to your downstream working variable name.

**Safety Intercepts & Failure Checks**

To prevent pipeline errors further down the script, this module includes immediate safety checks:

* **File Discovery Stop**: If the targeted folder directory is empty or missing, `import_qupath_data` immediately halts operations with a `stop()` command to prevent downstream function crashes.
* **Channel Warnings**: If an intersection search fails to map a target marker name, `extract_marker_intensities` logs a **CRITICAL WARNING** to your R console. This alerts you that a column name variant is missing from the search list, allowing you to update the naming definitions before running the rest of the workflow.

## Step 4a: Filter low-quality cells

Two stages of filtering

Filter 1: Remove cells with missing marker values instead it contains NA. 
Filter 2: Remove the potential debris (< 1) and clumps (>99) using the percentile of nucleus area (µm2). 
** Cell Quality Control & Filtering**

This pipeline module handles single-cell quality control (QC) exclusions. It programmatically screens out unsegmented debris, staining artifacts, and cell clumps using mathematical boundary thresholds to ensure only high-confidence cellular events reach downstream modeling blocks.

**Pipeline Module Overview**

The `filter_low_quality_cells` function cleans incoming datasets using a sequential, multi-stage gating strategy:
* **Completeness Enforcement**: Drops incomplete observations missing critical multiplexed channel readouts.
* **Morphological Trimming**: Truncates extreme structural outliers using automated percentile bounding.
* **Real-time Reporting**: Tracks exactly how many cell rows pass or fail the threshold benchmarks.


**Filtration Mechanics & Logic**

The function filters incoming single-cell tables using two specific data quality checkpoints:

**Checkpoint 1: Missing Footprint Removal**
Drops any cell record containing an `NA` missing entry across any primary target channel (`Marker_1` through `Marker_4`). This preserves downstream model matrix continuity.

**Checkpoint 2: Morphological Percentile Bounding**
If the column `Nucleus_Area_um2` exists in your dataset, the script builds data limits based on data distribution percentiles:

\[\text{Lower Bound} = \text{1st Percentile (0.01 Threshold Monitor)}\]
\[\text{Upper Bound} = \text{99th Percentile (0.99 Threshold Monitor)}\]

* **Debris Removal**: Cells falling below the 1st percentile are dropped as likely non-cellular staining debris or segmentation noise.
* **Clump Trimming**: Cells falling above the 99th percentile are dropped to exclude multi-cell aggregates, tissue folds, or optical fusion artifacts.


**Console Diagnostic Outputs**

When running, the module prints live processing updates to the R console. This lets you quickly monitor cell loss metrics:

```text
Filtering low-quality cells...
Removed 4321 cells (2.4%)
Retained 175679 cells
```

## Step 4b: Cell-line isolated principal component analysis (QC Step)

This is an additional QC step to know better about the data using the PCA plots.

**Cell-Line-Isolated Principal Component Analysis (QC Step)**

This pipeline module executes independent Principal Component Analysis (PCA) within each separate cell line background. Isolating the dimensional reduction step protects localized treatment variance from being masked by overwhelming baseline differences between distinct cell lines.

**Pipeline Module Overview**

The `run_cell_line_specific_pca` function processes single-cell matrices through an independent orthogonal variance loop:

* **Crash Protection**: Automatically validates input tables against targeted marker vectors, halting execution with clear messaging if columns are missing.
* **Smart Downsampling**: Includes a configurable cell ceiling parameter (`max_cells_per_line = 15000`) to subsample huge datasets. This prevents memory bottlenecks and graphical over-plotting while preserving overall population distributions.
* **Mathematical Centering**: Applies standard R `prcomp` algorithms to center and optionally scale (`scale_flag = TRUE`) log-transformed multiplexed channels.

**Environmental Dependencies**
The plotting routines look for specific variables and theme components built outside this module's footprint:

1. `treatments`: The standardized vector defining downstream factor level ordering.
2. `treatment_colors`: A named string array connecting hex color keys to each unique protocol condition.
3. `theme_publication()`: The publication-grade custom graphic layout engine.
4. `output_dir`: The global variable tracking where file assets are written to disk.


**Generated Visual Assets**

Running this module saves three distinct types of high-resolution visual plots (300 DPI) to your designated output folder:
** 1. Isolated PCA Scatter Landscapes** (`PCA_[suffix]_scatter_[cell_line].png`)
* **Visual Properties**: Renders individual single-cell data coordinates colored by treatment arm, utilizing an explicit 95% confidence data ellipse (`stat_ellipse`).
* **Sizing**: Fixed at 9" Width × 7" Height.

**2. Isolated Arrow Biplots** (`PCA_[suffix]_biplot_[cell_line].png`)
* **Visual Properties**: Layers vector loading segments directly over downscaled cell point coordinates. The direction and length of each arrow demonstrate how much a specific marker contributes to PC1 and PC2 variance.
* **Sizing**: Fixed at 9" Width × 7" Height.

**3. Global Comparative Faceted Panels** (`PCA_[suffix]_global_faceted.png`)
* **Visual Properties**: Merges all cell-line scores into a multi-column comparative layout grid (`facet_wrap`) with free coordinate scales. This lets you inspect overall treatment variance signatures in a single view.
* **Sizing**: Expanded template layout at 13" Width × 11" Height.



## Step 5: Ridge plot 

In this stage, high-resolution distributional ridge plots from transformed single-cell data were developed. It visualizes intesity shifts across different experiment treatments and cell lines.

The pipeline processes each marker sequentially and outputs standalone graphical files.
The images have a dimension of 11" width x 8" Height with 300 DPI resolution.

**Single-Cell Distributional Ridge Plots Visualization**

This pipeline module generates publication-quality, high-resolution density ridge plots. It stacks individual experimental treatment profiles to visually isolate unnormalized marker intensity shifts across different cell line backgrounds.


**Pipeline Module Overview**

The `generate_ridge_plots` function iterates through a vector of target markers to build standalone, multi-faceted visual arrays. 

* **Stacking Mechanics**: Layers density profiles along the y-axis to allow immediate comparison of marker distribution shapes between treatments.
* **Faceting Mechanics**: Splits data windows into isolated columns per `CellLine` with free x-axis scaling to accommodate native cell-line intensity baselines.


**Global Environmental Requirements**

To preserve pipeline automation, this function dynamically borrows three critical layout variables from your active runtime script workspace:

1. **`treatments`**: A character vector defining the definitive baseline order of your experimental arms. The script automatically reverses this vector (`rev(treatments)`) internally to ensure chronological top-to-bottom layering on the y-axis.
2. **`treatment_colors`**: A named character vector mapping explicit hex color codes to your specific treatment conditions.
3. **`theme_publication()`**: A custom, pre-configured `ggplot2` theme function that clears chart clutter and standardizes fonts.


**Output Asset Specifications**

Running this module outputs isolated image files directly into your specified destination directory:

* **File Name Scheme**: `ridge_plot_[MarkerName].png` (e.g., `ridge_plot_CD4.png`).
* **Resolution**: Fixed at 300 DPI (Print-Ready quality).
* **Dimensions**: 11 inches Width × 8 inches Height (Landscape layout).
* **Internal Performance Modifiers**: Sets `scale = 1.4` for optimal ridge overlap and drops tails via `rel_min_height = 0.005` to keep the canvas clean.



## Step 6: Single-Cell Mixed effects model

This module provides a robust, cell-isolated statistical framework for analyzing single-cell marker intensity shifts. It uses linear mixed-effects models (LMM) to account for random slide-to-slide variance, calculates bulletproof Cohen's d effect sizes, and evaluates sample normality.

**Pipeline Functions Overview**

The statistical pipeline runs across four sequentially linked operations:

1. **`transform_single_cells`**: Applies a deterministic `log2(x + 1)` transformation to protect against raw expression skewness.
2. **`mixed_effect_model`**: Fits a single-cell regression model (`Treatment * CellLine`) using `SlideID` as a random intercept. Computes direct, crash-proof Cohen's d magnitudes.
3. **`aggregate_to_samples`**: Compresses millions of single-cell events into robust sample-level slide medians.
4. **`check_sample_normality`**: Runs cell-line isolated Shapiro-Wilk normality profiles on the sample-level data to validate linear model assumptions.

 **Mathematical Logic & Effect Sizes**

To bypass the optimization crashes and division-by-zero errors common in standard R package effect-size calculators (like `effectsize::t_to_d`), this engine extracts standard deviation metrics directly from baseline controls:

$$\text{Cohen's } d = \frac{\text{Regression Coefficient Estimate}}{\text{Global Standard Deviation of "Untreated" Group}}$$

* If the "Untreated" group standard deviation is missing or evaluates to `0`, the engine defaults automatically to the **global pooled standard deviation** of that marker to protect pipeline continuity.
* **Effect Sizes** are thresholded via standard behavioral cutoffs: Negligible ($<0.2$), Small ($\ge0.2$), Medium ($\ge0.5$), and Large ($\ge0.8$).
* **P-Values** are adjusted for multiple testing across all markers using the **Benjamini-Hochberg (BH)** false discovery rate procedure.

Executing the primary model saves a focused comma-separated summary table directly to disk:

* **File Destination**: `[output_dir]/single_cell_mixed_model_effect_sizes.csv`
* **Column Metrics Provided**:
  * `Term`: The calculated baseline or interaction intercept (e.g., `Treatment_DrugA:CellLine_Line2`).
  * `Estimate`: Raw fixed-effect regression coefficient value.
  * `Cohens_d`: Direct effect sizing normalized to baseline group variance.
  * `P_Value`: Raw Satterthwaite analytical probability.
  * `adj_P_Val`: False-discovery optimized p-value (Benjamini-Hochberg corrected).
  * `Interpretation`: An automated text tag identifying if the result shows a **global drug shift** or a **cell-line specific modification**.
  
In this stage the viability data were read and standardized.

## Step 7: Non-parametric distribution ## 

**Non-Parametric Distribution Testing Engine**

This pipeline module provides an alternative statistical framework for handling non-normal single-cell populations. It uses distribution-free rank tests to analyze marker expression changes without assuming normal data curves.

**Pipeline Overview**

The `analyze_all_marker_distributions` function acts as a dual-stage non-parametric analytics engine:

1. **Global Variance Testing**: Computes a separate Kruskal-Wallis $H$ test for each cell line background to confirm significant population variations across all treatments.
2. **Targeted Post-Hoc Comparisons**: Extracts every treated condition and compares it back directly against its own native `Untreated` baseline via isolated pairwise Wilcoxon Rank-Sum (Mann-Whitney U) operations.

**Statistical Methods & Effect Sizes**

To balance out the high statistical power driven by massive single-cell event sizes (which can force tiny, meaningless variations to return extreme p-values), this engine couples p-values with strict effect size thresholds:

* **Rank-Biserial Correlation ($r$)**: Extracted using `effectsize::rank_biserial` to gauge the absolute degree of non-parametric group divergence.
* **Magnitude Thresholding**: Effect sizes are classified using customized single-cell biological impact tiers:
  * $r \ge 0.5$: **Large Biological Hit**
  * $r \ge 0.3$: **Moderate Shift**
  * $r \ge 0.1$: **Small/Weak Shift**
  * $r < 0.1$: **Negligible (Sample Size Noise)**
* **FDR Control**: Pairwise p-values are adjusted within each target marker profiling track using the **Benjamini-Hochberg (BH)** procedure. Significance stars are assigned directly based on adjusted parameters (`***` for $p < 0.001$, `**` for $p < 0.01$, `*` for $p < 0.05$, and `ns` for non-significant markers).

**Output File Specifications**

Running the distribution engine writes two distinct summary CSV spreadsheets to your output folder:

**1. Global Screening Report** (`single_cell_global_kruskal_wallis_results.csv`)
* Tracks overall variation signatures across all testing channels.
* Columns: `CellLine`, `Marker`, `kruskal_H` (test statistic), `p_value`, and `df` (degrees of freedom).

**2. Pairwise Contrast Report** (`single_cell_pairwise_wilcox_results.csv`)
* Tracks fine-grained treatment changes relative to localized controls.
* Columns:
  * `Marker` & `CellLine`: The exact target channel and biological background.
  * `Treatment`: The treated group contrasted against the baseline.
  * `wilcox_p`: The raw, unadjusted Mann-Whitney U test probability.
  * `wilcox_p_adj`: The final False Discovery Rate adjusted p-value.
  * `significance`: Readout flags mapping code markers (`***`, `**`, `*`, `ns`).
  * `effect_size_r`: Calculated absolute Rank-Biserial correlation coefficient value.
  * `effect_magnitude`: Explanatory label tracking single-cell signal shifts.

## Step 8: Integrated heatmap
**Integrated Heatmap Visualization & Interaction Deconvolution Engine**

This pipeline module handles advanced multi-omic data integration and downstream interpretation. It joins single-cell protein changes directly with cellular viability profiles on an aligned heatmap matrix, and features a deconvolution engine to untangle fixed-effect interaction terms.


**Pipeline Functions Overview**

This part of the analytical framework contains two core downstream processors:

1. **`create_integrated_heatmap`**: Reshapes single-cell regression estimates, fuses them with external bulk `Inhibition` percentage values, transforms viability values into a log2 fold-change format ($Log_2FC$), and renders a cell-line grouped matrix.
2. **`calculate_true_responses`**: Automates the statistical conversion of regression modifier terms into absolute, isolated fold changes and descriptive Cohen's d effect sizes for individual alternative cell lines.

**Interaction Deconvolution Mathematics**

Standard R regression summaries output coefficients relative to a reference intercept. Alternative cell line effects are reported strictly as comparative interaction modifiers (`Treatment:CellLine`). To establish the standalone biological response profile of alternative backgrounds, the deconvolution engine implements a recovery calculation:

$$\text{True } Log_2FC = \text{Baseline Estimate (Reference Cell Line)} + \text{Interaction Modifier Coefficient}$$

$$\text{True Cohen's } d = \text{Baseline Cohen's } d \text{ (Reference Cell Line)} + \text{Interaction Modifier Cohen's } d$$

Following deconvolution, Cohen's d magnitudes are re-evaluated across standard thresholds, and string hyphens are safely restored to cell-line identities for final publication tables.

**Output Asset Specifications**

Running this module saves two distinct assets to the designated destination folder:

 **1. High-Resolution Heatmap Plot** (`integrated_heatmap_grouped.png`)
* **Structural Behavior**: Order tracking is locked into your exact `treatments` factor configuration vector with hierarchical clustering disabled (`cluster_rows = FALSE`) to prevent alignment crashes.
* **Visual Annotations**: Row side-bars map cellular viability changes ($Log_2FC$) colored along an Indigo gradient (`#4B0082` to `#E6E6FA`).
* **Matrix Sizing**: Dimensions are fixed at 2600 × 3000 pixels at a print-ready 300 DPI resolution, segmented by cell line grid break gaps.

**2. Standalone Response Spreadsheet** (`true_cell_line_isolated_responses.csv`)
* **`CellLine` & `Treatment`**: Realigned experimental labels with standard hyphen formatting.
* **`True_Log2FC`**: The calculated real fold-change response.
* **`True_Cohens_d`**: Deconvoluted standard deviation magnitude tracking.
* **`True_Effect_Magnitude`**: Updated behavioral significance categorization tags based on true isolated values.


## Step 9: Main analysis pipeline ##
**Main Analysis Pipeline Orchestration Engine**

This is the primary orchestration script (`run_analysis`) for the multiplexed immunofluorescence (mIF) single-cell processing pipeline. It standardizes, filters, transforms, models, and visualizes multiplexed cell data inside a single command loop.

**Pipeline Architecture Overview**

The pipeline acts as a central coordinator that automatically handles data ingestion, quality control auditing, data transformations, dimensional reduction, and downstream statistical modeling. It enforces reproducible analytical rules across all patient slides, experimental conditions, and cell lines.

 **Data Flow & Execution Order**

Executing `run_analysis` triggers a two-tiered processing pipeline:

 **STEP 1: Data Import, Ingestion, and Quality Control**
1. **Raw Intake**: Consolidates cell-level data sheets from QuPath exports (`import_qupath_data`).
2. **Column Standardization**: Cleans raw syntax characters into compliant R names (`clean_column_names`).
3. **Stage 0 Logging**: Stores absolute starting raw counts.
4. **Filter Stage 1**: Flags and removes incomplete footprints missing marker intensity attributes.
5. **Filter Stage 2**: Filters out low-quality debris using cell size and shape rules (`filter_low_quality_cells`).
6. **Audit Compilation**: Outputs a cell-retention data sheet and simultaneously loads bulk cellular viability indexes (`import_viability_data`).

**STEP 2: Downstream Analysis Pipelines**
1. **Mathematical Transformation**: Applies a deterministic `log2(x + 1)` scaling formula across all targets.
2. **Ridge Profiling**: Generates single-cell population density curves (`generate_ridge_plots`).
3. **Dimensionality Reduction**: Downsamples and runs cell-line isolated Principal Component Analysis (`run_cell_line_specific_pca`).
4. **Aggregation & Normality**: Compresses cell records into slide medians, running Shapiro-Wilk checks.
5. **Parametric Modeling**: Fits linear mixed-effects model tracks using slides as random nested blocks.
6. **Multi-Omic Matrixing**: Pairs single-cell estimates with viability tables to render an aligned heatmap.
7. **Non-Parametric Contrast**: Runs Kruskal-Wallis screenings and pairwise Wilcoxon calculations.

 **Cell Filtration Audit Matrix**

During **Step 1E**, the pipeline compiles a strict data quality control (QC) table. This matrix calculates exactly how many single cells are excluded by your processing choices:

* **Console Logging**: A summary report automatically prints directly to the R console during processing for rapid pipeline feedback.
* **Audit File Summary (`cell_filtration_audit.csv`)**: A permanent, publication-ready table is saved to your disk containing these parameters:
  * `Stage0_Raw_Cells`: Total input cell count immediately after column standardization.
  * `Stage1_Complete_Markers`: Total cells remaining after dropping missing marker measurements.
  * `Stage2_Final_Retained`: Absolute count of high-quality cells used in modeling.
  * `Total_Cells_Removed`: Aggregate cell counts lost during filtering.
  * `Pct_Cells_Retained`: The precise percentage of biological data preserved ($100 \times \frac{\text{Stage2}}{\text{Stage0}}$).



**Master Output Package Architecture**

The function finishes by packing all data frames, matrices, and statistical indices into a named list object returned directly to your active global environment:

| List Key Element | Data Structure Type | Descriptional Use Case |
| :--- | :--- | :--- |
| `single_cell_processed` | Data Frame | Transformed, high-quality cell coordinates and marker metrics. |
| `sample_medians` | Data Frame | Slide-level median calculations used for assumption validation. |
| `pca_results` | Custom List | Eigenvalues, loadings, and coordinate embeddings from the PCA module. |
| `final_mixed_results` | Data Frame | Coefficients, Satterthwaite p-values, and adjusted Cohen's d values. |
| `distribution_stats` | Custom List | Parallel global Kruskal-Wallis and pairwise post-hoc Wilcoxon metrics. |
| `heatmap_data` | Wide Matrix | Fused data layout combining mixed-model results and viability profiles. |
| `filtration_report` | Data Frame | Cell conservation metrics saved to `cell_filtration_audit.csv`. |
| `viability` | Data Frame | Standardized percent growth inhibition metrics. |

## Step 10: Execute analysis
**Configuration & Pipeline Execution**

This final code block serves as the user-facing entry point of the pipeline script. It calls the master orchestration engine and saves all processed matrices and statistical arrays directly into your active R global environment.


**User Configuration Checklist**

Before executing the script, verify that your active R session contains the required input paths and parameters. The pipeline expects three variables to be pre-defined in a dedicated **CONFIGURATION** section at the top of your script:

* **`data_dir`**: A character string pointing to the local folder containing your raw QuPath segmentation `.csv` files.
* **`viability_file`**: A character string pointing directly to your bulk metabolic cell viability Excel document (`.xlsx`).
* **`output_dir`**: A character string specifying where all generated PNG charts, heatmaps, and audit spreadsheet files will be saved.

**Workspace Persistence**

When the script finishes running, it outputs a master variable named `results` into your active R environment. This variable holds the complete structured output package, allowing you to run ad-hoc queries, extract data subsets, or create custom charts without re-running the entire pipeline.

```R
# Quick Diagnostic Examples on Saved Workspace Objects:

# 1. View the cell-retention quality audit metrics
print(results\$filtration_report)

# 2. Extract top-performing drug treatments from the linear mixed model data
library(dplyr)
significant_hits <- results\$final_mixed_results %>%
  filter(adj_P_Val < 0.05 & abs(Cohens_d) >= 0.8)

# 3. Check the internal matrix coordinates used to render the heatmap
head(results\$heatmap_data)
```
