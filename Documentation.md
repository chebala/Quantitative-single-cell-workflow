The workflow contains 10 stages.

This workflow supports CSV files. 

**Stage 1: Setup and package loading**

In this stage, all the required packages are loaded. 
Any packages that were not available will be installed automatically from CRAN repository.
limma package downloaded via Biomanager.

**Stage 2: Configuration**

Define the paths to CSV files that contains the measurements. Make sure to standardize the naming of files as the workflow is case sensitive

The output can be stored in a seperate directory (If the path does n't exist the workflow itself produce and save the output files.
