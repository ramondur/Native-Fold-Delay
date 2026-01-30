# Native Fold Delay (NFD)

This repository provides an R script to calculate the Native Fold Delay (NFD) of proteins at proteome scale. The workflow is designed for AlphaFold-predicted structures and therefore requires UniProt identifiers to automatically retrieve the correct models and metadata.

For users interested in single-protein analysis, including non-AlphaFold PDB structures, we also provide an interactive web server:

🔗 https://folddelay.switchlab.org

The web server supports NFD calculation as well as interactive analysis and visualization of the results.

---

## Repository contents

- `FoldDelay_script.R`  
  Main R script for large-scale NFD calculations using AlphaFold structures.

- `Duran&Houben2025/`  
  Precomputed NFD datasets for:
  - *Saccharomyces cerevisiae* (yeast)
  - *Escherichia coli*  

  These datasets were used in the publication:  
  **Duran & Houben (2025)**, *Native Fold Delay and its implications for co-translational chaperone binding and protein aggregation*.

---

## What the script does

When executed, the script performs the following steps:

1. **Structure retrieval**  
   Automatically downloads AlphaFold models and their corresponding Predicted Aligned Error (PAE) matrices from the AlphaFold Database for the provided UniProt IDs.

2. **Contact identification**  
   Identifies all residue-residue contacts within a specified distance cutoff (default: 6 Å), considering only N-terminal to C-terminal residue pairs.

3. **Native Fold Delay estimation**  
   For each contact, estimates the minimal time required for its formation by:
   - Computing the sequence separation between residues
   - Converting this separation into time using a global translation rate (default: 5 aa/s)

4. **Confidence filtering**  
   Contacts are filtered using PAE values to remove residue pairs with low positional confidence, reducing false-positive interactions.

5. **Domain annotation**  
   Protein domain boundaries are automatically retrieved from the Encyclopedia of Domains, allowing contacts to be classified as:
   - Intra-domain
   - Inter-domain

6. **Contact selection**  
   By default, only the most C-terminal interacting partner for each residue is retained. This behavior can be changed to keep all contacts (see command-line options).

---

## Running the script

###  Requirements

- **R** (≥ 4.0 recommended): https://cran.rstudio.com  
- Internet access for AlphaFold data retrieval

###  Input

A plain text file containing UniProt IDs, one per line:

```
P23869
Q9Y2X3
P0A7V8
```

### Basic Usage

```
Rscript FoldDelay_script.R --input uniprot_ids.txt
```

###  Command-Line Options

```
-i, --input CHARACTER
    Path to a text file containing UniProt identifiers.

-t, --transrate NUMBER
    Translation rate (aa/s). Default: 5

-d, --distance NUMBER
    Distance cutoff (Å). Default: 6

-k, --keep BOOLEAN
    TRUE  → keep all contacts
    FALSE → keep only most C-terminal partner per residue (default)
```

###  Examples

```
Rscript FoldDelay_script.R --input yeast_uniprot_ids.txt --transrate 20 --distance 6 --keep TRUE
```

### Output

The script produces a single CSV file, folddelay.csv, written to the current working directory. Each row represents a predicted residue-residue contact used in the Native Fold Delay calculation.
The output contains the following columns:

- **UniProt id:** UniProt accession identifier of the protein.
- **Index 1:** Index of the N-terminal residue involved in the contact
- **Aa 1:** Amino acid identity of the N-terminal residue
- **Index 2:** Index of the C-terminal residue involved in the contact
- **Aa 2:** Amino acid identity of the C-terminal residue
- **NFD (aa):** Sequence separation between the two residues (in amino acids)
- **NFD (s):** Estimated Native Fold Delay (NFD) time for the contact, based on the specified translation rate
- **Cath label:** CATH domain classification assigned to the residue
- **Domain ID:** Numeric identifier of the domain to which the residue belongs
- **Scope:** Domain interaction type:
  - `intra` — both residues belong to the same domain  
  - `inter` — residues belong to different domains
    
---

## Contact

For any questions or further information, please contact: [ramonenrique.duranromana@kuleuven.be](mailto:ramonenrique.duranromana@kuleuven.be)

---

## Citation

If you use this repository, please cite:
Duran-Romana, R. et al. Native fold delay and its implications for co-translational chaperone binding and protein aggregation. Nature Communications 16, 1673 (2025). https://doi.org/10.1038/s41467-025-57033-z
