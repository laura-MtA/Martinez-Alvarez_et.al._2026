# Figure 1 

This folder contains scripts used to generate the data and plots shown in **Figure 1** of the manuscript.  
The scripts are provided to support transparency and reproducibility of the analyses underlying each panel, rather than as a fully automated pipeline.

**Key input tables (provided):**
- `padloc_bac_noVSPRDMS.txt`
- `padloc_ar53_noVSPRDMS.txt`

### Panel 1A – Comparison of PADLOC and DefenseFinder predictions

Scripts associated with **Fig. 1A** quantify the overlap between defense-related proteins predicted by **PADLOC** and **DefenseFinder**.
**Key scripts:**
- `script_cdhitoutput.py` – parses CD-HIT `.clstr` files into a long-format table (cluster ID + protein ID)
- `script_compare.sh` – classifies clusters by tool membership and reports counts used for the Venn diagram

### Panel 1D – Distribution of defense systems per genome

Scripts associated with **Fig. 1D** generate histograms showing the distribution of the number of defense systems per genome, comparing Archaea and Bacteria.
**Key scripts:**
- `script_histogram.R` – reads the tables above and produces the histogram shown in Fig. 1D

### Panels 1E–F – Defense systems across archaeal phyla

Scripts associated with **Fig. 1E–F** summarise the average number of defense systems per genome and the average number of defense system families per genome across archaeal phyla
**Key scripts:**
- `script_arch_stats.R` - produces the bar plots shown in Fig.1E-F

For further details on methodology and filtering choices, please refer to the Methods section of the manuscript.
### R environment

Analyses and plotting scripts were written in **R** (v4.4.0) and developed using **RStudio** (2024.04.1 Build 748).
