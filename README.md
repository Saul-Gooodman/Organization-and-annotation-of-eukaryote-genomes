
⸻

#   Organization and Annotation of Eukaryote Genomes

##  Project Overview

This repository contains all scripts and analysis results generated for the “Organization and Annotation of Eukaryote Genomes (UE-SBL.30004)” course at the University of Bern.
The workflow covers genome annotation, transposable element (TE) analysis, gene prediction and refinement, functional annotation, and comparative genomics using OrthoFinder and GENESPACE.

⸻



##  Repository Structure

```
Organization-and-annotation-of-eukaryote-genomes/
├── scripts/                     # All scripts used for each analysis step
│   ├── 01_run_EDTA.sh
│   ├── 02_run_TEsorter.sh
│   ├── ...
│   └── 25_genespace_core_accessory_from_txt.sh
│
└── results/                     # Output results organized by week
    ├── 1. TE Annotation using EDTA
    ├── 2. Visualizing and comparing TE annotations
    ├── 3. Refining TE Classification with TEsorter
    ├── 4. Dynamics of Transposable Elements (TEs)
    ├── 5. Run MAKER with MPI
    ├── 6. Filtering and Refining Gene Annotations
    ├── 7. BUSCO
    ├── 8. Visualizing Gene annotation with Geneious
    ├── 9. Sequence homology to functionally validated proteins
    └── 10. Comparative Genomics with OrthoFinder and GENESPACE
```


⸻

##  Uploaded Content

Includes all scripts and corresponding result files for:
	•	TE annotation and classification
	•	Gene prediction and refinement
	•	Functional annotation
	•	BUSCO completeness assessment
	•	Functional validation using UniProt
	•	Comparative genomics (OrthoFinder + GENESPACE)

⸻

##  Software and Environment

All analyses were executed on the University of Bern bioinformatics cluster using Singularity containers.

Main tools used:
	•	EDTA, TEsorter, MAKER, BUSCO, AGAT, OrthoFinder, GENESPACE
	•	R 4.1.0, Python 3.9, DIAMOND 2.19, MCScanX

⸻

##  Author

Yuwei Liu（24-100-430）
Master student in Bioinformatics and Computational Biology
University of Bern


⸻

