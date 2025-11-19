# Organization and Annotation of the Kas-1 Genome

Author: **Yuwei Liu** (Student ID: 24-100-430)  
Course: *Organization and annotation of eukaryote genomes* 
Accession analysed: **Arabidopsis thaliana – Kas-1**

---

## 1. Project Overview

This repository contains all scripts and major output files generated for the analysis of the Arabidopsis thaliana accession **Kas-1**.

The work includes:

- Transposable element (TE) annotation with **EDTA**
- TE classification refinement with **TEsorter**
- TE abundance, density, and divergence (age) analyses
- Gene prediction and annotation using **MAKER**
- Gene-set completeness assessment with **BUSCO**
- Gene/TE visualization (Geneious / IGV)
- Comparative genomics with **OrthoFinder** and **GENESPACE**

The goal of the repository is to make the workflow transparent and reproducible and to store all major figures and outputs referenced in the final report.

---

## 2. Final Report

The final PDF report for the course is:

### **`Kas-1_Report_YuweiLiu.pdf`**

The report includes selected figures and tables with extended figure legends explaining:

1. What the figure or table shows  
2. How the reader should interpret it  
3. A short scientific interpretation

A link to this GitHub repository is provided inside the report.

---

## 3. Repository Structure

### **3.1 `results/` — Main Output Files**

The folder `results/` contains the key output files for each major step of the pipeline.  
Subdirectories follow the numbering of the course practicals:

- **1. TE Annotation using EDTA**  
  EDTA GFF outputs, TE summaries, and TE composition plots.

- **2. Visualizing and comparing TE annotations from EDTA**  
  Circos-style TE density plots and distribution summaries.

- **3. Refining TE Classification with TEsorter**  
  Clade-level classification (Copia/Gypsy), bar plots, and summary tables.

- **4. Dynamics of Transposable Elements (TEs)**  
  Divergence histograms and TE landscape plots (age structure).

- **5. Run MAKER with MPI**  
  MAKER GFF annotation files, transcript/protein FASTA files, and logs.

- **6. Filtering and Refining Gene Annotations**  
  Filtered gene set (`filtered.genes.renamed.gff3`) and summary statistics.

- **7. BUSCO**  
  BUSCO output tables and summary figures.

- **8. Visualizing Gene Annotation with Geneious**  
  IGV/Geneious snapshots showing gene models and nearby TEs.

- **9. Sequence homology to functionally validated proteins**  
  BLAST results against UniProt and TAIR10 (used to calculate gene hits).

- **10. Comparative Genomics with OrthoFinder and GENESPACE**  
  OrthoFinder orthogroups, GENESPACE synteny plots, and genome comparison outputs.

---

### **3.2 `scripts/` — Pipeline Scripts**

The `scripts/` folder contains all shell and R scripts used to run the complete workflow.  
Script numbers follow the practical session order.

---

## 4. Key Statistics for Kas-1

Key genome statistics used in the report:

- **Assembly size (quickmerge):** 137.2 Mb  
- **Contig number:** 35  
- **Contig N50:** 8.6 Mb  
- **TE content:** ~16% of the genome  
- **Filtered genes:** 30,313  
- **Genes with BLAST hits:**  
  - UniProt: 26,767  
  - TAIR10: 23,550  
- **Core orthogroups (TAIR10 + Kas-1):** 20,085  
- **Kas-1-specific orthogroups:** 5,687  

These statistics reflect the major features of the Kas-1 genome structure and annotation quality.

---
