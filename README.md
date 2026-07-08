# Multi-Omic Deconvolution of Phenotypic Plasticity & Quantum Evolution in Glioblastoma (GBM) 
note:- dont use this project cuz work is still needed and also this is so meesy 

This repository contains the complete codebase, manuscript drafts, analysis scripts, and quantum simulation resources for our multi-omic study investigating chemotherapy-induced proneural-to-mesenchymal transition (PMT) and subclonal evolution in glioblastoma.

---

## 1. Project Structure & File Guide

### 📄 Research Manuscripts & Findings
*   **[GBM_Research_Paper_NMohanaKrishna.md](GBM_Research_Paper_NMohanaKrishna.md)**: Main research manuscript covering the STAT3 PMT driver, the ATRX genetic fork, the IDH-mutant progression map, CD44/YKL-40 deconvolution, and proteomic subtype classification.
*   **[Unbiased_Driver_and_CpG_Switches_Paper.md](Unbiased_Driver_and_CpG_Switches_Paper.md)**: Dedicated manuscript detailing Stage 4 and Stage 5 results, including unbiased mutation scans, pathway burden tests, matched longitudinal promoter CpG methylation correlation, and independent CGGA cross-validation.
*   **[1st_findings.md](1st_findings.md)**, **[2nd_findings.md](2nd_findings.md)**, **[3rd_findings.md](3rd_findings.md)**, **[4th_findings.md](4th_findings.md)**, **[stage4_5_findings.md](stage4_5_findings.md)**: Progressive stage logs and scientific findings documents.

### 🧪 Analysis & Testing Scripts
*   **[stage1_unbiased_scan.py](stage1_unbiased_scan.py)**: Performs Mann-Whitney U tests on 3,034 mutated genes against YKL-40 and CD44 expression fold-changes with Benjamini-Hochberg FDR correction.
*   **[pathway_burden_scan.py](pathway_burden_scan.py)**: Groups somatic mutations into 1,084 Reactome pathways and correlates pathway-level mutation burden with expression fold-changes.
*   **[stage5_epigenome_expression_check.py](stage5_epigenome_expression_check.py)**: Correlates matched longitudinal promoter CpG methylation changes ($\Delta\text{Beta}$) with matched RNA expression changes ($\Delta\text{RNA}$) in $n=28$ patients.
*   **[cgga_cpg_cross_validation.py](cgga_cpg_cross_validation.py)**: Cross-validates baseline promoter CpG methylation-expression correlations in the independent CGGA IDH-A cohort.
*   **[test_hypothesis_1.py](test_hypothesis_1.py)** $\rightarrow$ **[test_hypothesis_4.py](test_hypothesis_4.py)**: Individual test scripts for the four longitudinal evolutionary hypotheses (Clonal Selection vs. Plasticity, Subtype Shuffling, ATRX Decoupling, and the MMR hypermutation reversal audit).
*   **[run_seurat_pipeline.R](run_seurat_pipeline.R)**: Seurat (v4) single-cell RNA-seq clustering and UMAP pipeline for 74,563 glioblastoma cells.

### ⚛️ Quantum Simulation (Q-TEM / Q-HTM)
*   **[qtem_simulator.py](qtem_simulator.py)**: Program to compile our 20-qubit Quantum Tumor Evolution Map circuit, run local Aer statevector simulations, and submit transpiled jobs directly to physical IBM Quantum Heron QPUs (`ibm_marrakesh`).
*   **[implementation_plan.md](implementation_plan.md)**: Design plan for the Quantum Hamilton Transition Map (Q-HTM), mapping our 30-qubit architecture to the five critical spatio-temporal and clonal gaps in GBM research.

---

## 2. Dataset Access & Directory Instructions

Due to file size limits, the raw datasets are stored locally or must be downloaded from their primary repositories:

### 🧬 Longitudinal GLASS Cohort (Western Population)
*   **WES/WGS Mutations**: `variants.passgeno.csv.gz` (contains somatic mutations across 271 longitudinal cases).
*   **RNA-seq Matrix**: `gene_tpm_matrix_all_samples.tsv` (case identifiers parsed via dot-separated values, e.g., `GLSS.19.0266.TP.01R...`).
*   **DNA Methylation**: Illumina 450K array beta values.
*   **Source**: Download via Synapse (ID: `syn17010685`) or the GLASS Consortium Portal (https://www.glass-consortium.org/).

### 🇨🇳 CGGA Cohort (Independent Validation & Single-Cell)
*   **Deep Multi-Omic Cohort**: 35 IDH-mutant, non-codel astrocytoma patients with matched WES, RNA-seq, EPIC DNA Methylation, Mass Spectrometry Proteomics, and Phosphoproteomics.
*   **Single-Cell RNA-seq**: 16 patients (74,563 single cells) stored as `scRNAseq.rdata` / `CGGA_IDH_A_scRNA_SeuratObj_20250915.rdata`.
*   **Source**: Download from the Chinese Glioma Genome Atlas database (http://www.cgga.org.cn/).

---

## 3. Core Scientific Discoveries

1.  **Dual-Compartment Mesenchymal Illusion**: We show that the mesenchymal transition in bulk datasets is a hybrid signal. **CD44** is predominantly expressed by tumor-associated macrophages (forming an immune shield), while **YKL-40 (CHI3L1)** is predominantly expressed by cancer cells (driving cell-intrinsic invasion).
2.  **The ATRX Decoupling Fork**: The ATRX mutation status does not directly alter YKL-40 baseline expression. Instead, it acts as a trans-regulatory genetic fork that decouples YKL-40 expression from macrophage-driven microenvironmental signaling.
3.  **Chemotherapy-Induced Subtype Instability**: paied Wilcoxon testing on matched patient lines showed that chemotherapy does not drive a linear progression to a stem state (PPR score decreased, $p=0.0031$), but rather triggers active, bidirectional epigenetic shuffling (subtype switching in 52.51% of cases).
4.  **FDR Unbiased Scans**: Scans across 3,034 genes and 1,084 Reactome pathways revealed that no single driver mutation directly commands YKL-40 or CD44 expression (all $q \approx 0.996$ after BH-FDR correction), demonstrating that resistance is polygenic/diffuse.
5.  **Validated CpG Promoter Switches**: Cross-validation in the independent CGGA cohort confirmed three promoter CpG switches (`cg15490070`, `cg03625911`, and `cg08530414`) that consistently correlate with expression across cohorts.
6.  **Quantum State Mapping**: Using our 20-qubit circuit (Q-TEM) on IBM Quantum hardware (`ibm_marrakesh`), we mapped tumor subclonal transitions into quantum phase space, proving that cancer cell state evolution can be represented as a physical trajectory in Hilbert space.
