# Native Fold Delay (NFD)

This repository provides tools for calculating the **Native Fold Delay (NFD)** of any protein using an R notebook (`FoldDelay.Rmd`). Below is an overview of the contents and instructions for use.

In addition, it also contaisn two  datasets with NFD profiles for the *Saccharomyces cerevisiae* (Yeast) and *Escherichia coli* (E. coli) proteomes within **`Duran&Houben2025/`** folder.
Both datasets were used in the publication: "Native Fold Delay and its implications for co-translational chaperone binding and protein aggregation" (Duran & Houben, 2025).

**Note: For calculating the NFD of single proteins you can try our web server: https://folddelay.switchlab.org.** 

## Run NFD tool

1. **Software Requirements**  
   - Install **R** and **RStudio** to execute the code in the R notebook.  
     Download the latest version here: [RStudio Desktop](https://posit.co/download/rstudio-desktop/).

2. **Instructions for Use**
   - Download and unzip the repository.
   - Open the **`NativeFoldDelay_tool/`** folder.
   - Place the PDB file of the protein to analyze in the `Structures/` folder (preferentially an AlphaFold protein structure obtained from the AlphaFold Protein Structure Database https://alphafold.ebi.ac.uk).  
   - Open `FoldDelay.Rmd` in RStudio and follow the instructions within the notebook to perform your analysis.  
   - Results will be automatically saved in the `Output/` folder.
  
### Case Study
To help you get started, we provide a case study using the AlphaFold model of the *E. coli* protein **Peptidyl-prolyl cis-trans isomerase B**, included in this repository.

### Important
Currently, PAE filtering is supported only for AlphaFold structures downloaded from the AlphaFold Protein Structure Database (https://alphafold.ebi.ac.uk) using version 6, which is the current release.

## Contact Information
For any questions or further information, please contact: [ramonenrique.duranromana@kuleuven.be](mailto:ramonenrique.duranromana@kuleuven.be)

## References
If you use this repository, please cite:
Duran-Romana, R. et al. Native fold delay and its implications for co-translational chaperone binding and protein aggregation. Nature Communications 16, 1673 (2025). https://doi.org/10.1038/s41467-025-57033-z


