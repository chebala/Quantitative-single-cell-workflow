The workflow contains 10 stages.

This workflow supports CSV files. 

**Stage 1: Setup and package loading**

In this stage, all the required packages are loaded. 
Any packages that were not available will be installed automatically from CRAN repository.
limma package downloaded via Biomanager.

**Stage 2: Configuration**

Define the paths to CSV files that contains the measurements. Make sure to standardize the naming of files as the workflow is case sensitive
The output can be stored in a seperate directory (If the path does n't exist the workflow itself produce and save the output files.
With this setup all the input files and output files will be handled externally in either local drive or external hard drives. This allows to save R memory.
Here you can define the cell lines, treatment conditions, and markers.
All the plots were made by a standardized theme together with specific colors for treatment and cell lines.

**Step 3: Data Import and Preprocessing**

QuPath CSV files were combined and a master dataframe was created.
If there is no CSV file found in the specified path, the workflow will be halted.
The cell line names and treatment condition were standardized.
All the marker intensities werecaptured and checks whether the columns actually exists.

**Step 4: Filter low-quality cells**

Two stages of filtering

Filter 1: Remove cells with missing marker values instead it contains NA. 
Filter 2: Remove the potential debris (< 1) and clumps (>99) using the percentile of nucleus area (µm2). 

Step 4b: Cell-line isolated principal component analysis (QC Step) 

This is an additional QC step to know better about the data using the PCA plots.

Step 5: Ridge plot 

