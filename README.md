# Native Fold Delay (NFD)

This repository provides tools for calculating the **Native Fold Delay (NFD)** of any protein using an R notebook (`FoldDelay.Rmd`). Below is an overview of the contents and instructions for use.

In addition, it also contaisn two  datasets with NFD profiles for the *Saccharomyces cerevisiae* (Yeast) and *Escherichia coli* (E. coli) proteomes within **`Duran&Houben2025/`** folder.
Both datasets were used in the publication: "Native Fold Delay and its implications for co-translational chaperone binding and protein aggregation" (Duran & Houben, 2025).

## Run NFD tool

1. **Software Requirements**  
   - Install **R** and **RStudio** to execute the code in the R notebook.  
     Download the latest version here: [RStudio Desktop](https://posit.co/download/rstudio-desktop/).

2. **Instructions for Use**
   - Download and unzip the repository.
   - Open the **`NativeFoldDelay_tool/`** folder.
   - Place the PDB file of the protein to analyze in the `Structures/` folder (preferentially an AlphaFold protein structure).  
   - Open `FoldDelay.Rmd` in RStudio and follow the instructions within the notebook to perform your analysis.  
   - Results will be automatically saved in the `Output/` folder.
  
### Case Study
To help you get started, we provide a case study using the AlphaFold model of the *E. coli* protein **Peptidyl-prolyl cis-trans isomerase B**, included in this repository.

## Contact Information
For any questions or further information, please contact: [ramonenrique.duranromana@kuleuven.be](mailto:ramonenrique.duranromana@kuleuven.be)

