Prediction_of_defense_systems/README.md
## Prediction of defense systems

This directory contains scripts associated with the prediction of microbial defense systems used in the analyses described in:

**Martínez-Alvarez et al. (2026)**  
*Diversity and Evolution of Archaeal Immune Strategies*

The scripts are provided as **examples of the analytical procedures** applied in this part of the study. They are **not intended to function as a complete or automated pipeline** and may require modification to run on other datasets or computing environments.

---

## Requirements

The analyses illustrated here rely on the following external tools (versions used in the study):

- **PADLOC** v2.0.0  
  https://github.com/padlocbio/padloc-db/tree/master  
  Payne,L.J., Hughes,T.C.D., Fineran,P.C. and Jackson,S.A. (2024) New antiviral defences are genetically embedded within prokaryotic immune systems. 10.1101/2024.01.29.577857.  
  Payne,L.J., Todeschini,T.C., Wu,Y., Perry,B.J., Ronson,C.W., Fineran,P.C., Nobrega,F.L. and Jackson,S.A. (2021) Identification and classification of antiviral defence systems in bacteria and archaea with PADLOC reveals new system types. Nucleic Acids Research, 49, 10868–10878.
- **DefenseFinder** v1.2.2  
  https://github.com/mdmparis/defense-finder  
  Tesson,F., Hervé,A., Mordret,E., Touchon,M., d’Humières,C., Cury,J. and Bernheim,A. (2022) Systematic and quantitative view of the antiviral arsenal of prokaryotes. Nat Commun, 13, 2561.  
- **CRISPRCasTyper (CCTyper)** v1.8.0  
  https://github.com/Russel88/CRISPRCasTyper  
  Russel,J., Pinilla-Redondo,R., Mayo-Muñoz,D., Shah,S.A. and Sørensen,S.J. (2020) CRISPRCasTyper: Automated Identification, Annotation, and Classification of CRISPR-Cas Loci. CRISPR J, 3, 462–469.
- **MMseqs2** v17.b804f  
  https://github.com/soedinglab/MMseqs2  
  Steinegger,M. and Söding,J. (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol, 35, 1026–1028.

---

## Repository contents

The files in this directory exemplify the custom scripts and computational steps used to identify and annotate defense systems.  
They are intended to support **transparency and reproducibility** rather than exact reproduction of the full workflow.

File paths, parameters and intermediate steps may be simplified for clarity.

---

## Notes on reproducibility

Exact reproduction of the results presented in the manuscript may additionally require:
- Specific database versions
- Manual filtering or curation steps
- Software dependencies described in the manuscript and supplementary materials
