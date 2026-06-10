# Quantitative-single-cell-worflow-cytospin-based-
A base quantitative R based workflow for single cell analysis. The workflow is based on the intensity measurements captured by either olympus slide scanner. This workflow was developed with cytopsin/PFA fixed cell line samples (Untreated/Treated).

In this workflow, the images were processed in QuPath and segmented with Instanseg extension in QuPath. After manually validated the segmentation, the measurements were exported as detections in a CSV file. 

The workflow comprises two stage of filtering: 

Filter 1: Removes the cells (each row) which lack the intensity measurement for any of the 4(or more) markers. These cells usually contains "NA" instead of a value.

Filter2: As the cells were treated with inhibitors, the cells are likely smaller (pyknotic fragments) or large clumped debris which can be formed due to the cytospin method. To remove these artefacts, the cells were filtered based on nucleus area (µm2) [>1 percentile and <99 percentile].

As the raw intensities were usually skewed, the raw intesity values were log2(1+x) transformed. 

This normalization can be done with Z-scores [(x-median)/MAD] which I used for correlation analysis.

log2FC can be applied to determine the relative treated median intensity measurements when compared to the untreated median intensity, This log2FC can be used to make heatmaps. Here I used heatmap that integrates cell viability data of each condition(which I got through a cell viability assay). 

Using Z- scores, the expression distribution of each marker under different treatment condition can be visualized though ridge plots or violin plots with an integrated box plots.

**Notes and Limitations**

This workflow was developed using cytospin-prepared, PFA-fixed cell line samples and should be considered a baseline workflow. Users are encouraged to evaluate whether the filtering thresholds, normalization methods, and statistical approaches are appropriate for their own datasets before applying them.

Because analyses are performed at the single-cell level, pseudoreplication may influence statistical significance. Reporting effect sizes alongside p-values and incorporating biological replicates where possible is recommended.

**Feedback and Contributions**

Questions, suggestions, and bug reports are welcome through Issues.




<img width="748" height="373" alt="image" src="https://github.com/user-attachments/assets/a016797d-bd9b-4a17-abf7-edd5b5aaff3b" />
