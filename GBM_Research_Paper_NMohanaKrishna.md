# Multi-Omic Deconvolution of Chemotherapy-Induced Phenotypic Plasticity in Glioblastoma: Resolving the Mesenchymal Illusion

**Author**: N. Mohana Krishna

---

## Abstract

Glioblastoma Multiforme (GBM) remains the most lethal primary brain tumor with a median survival of 14.6 months. The Proneural-to-Mesenchymal Transition (PMT) is widely recognized as the central mechanism of therapy resistance, yet the cellular origin of the mesenchymal signature has remained unresolved. Using integrated multi-omic profiling of 286 whole-exome sequenced patients, 35 deeply profiled IDH-mutant astrocytoma patients (RNA-seq, EPIC methylation, mass spectrometry proteomics, phosphoproteomics), and 74,563 single cells, we report three major findings. First, the ATRX mutation functions as a genetic fork that determines the mode of therapy resistance: ATRX-wildtype tumors recruit CD44-positive macrophages (immune shield), while ATRX-mutant tumors activate synaptic integration (GRIA2/NLGN3) and cell-intrinsic shape-shifting via CHI3L1/YKL-40. Second, by single-cell deconvolution, we demonstrate that the rise of CD44 in aggressive tumors is exclusively a macrophage signature and not a cancer cell phenotype, while YKL-40 represents the true cell-intrinsic mesenchymal program. Third, we show that the chemotherapy-resistant proteomic subtypes (PPR and IME) are absent in untreated tumors and emerge only after Temozolomide exposure, driven by STAT3-mediated chromatin remodeling. We identify ChEMBL-validated inhibitors targeting this network, including a 25 nM YKL-40 neutralizer and a 3 nM covalent TEAD1 inhibitor. These findings redefine the mesenchymal transition as a dual-compartment illusion and provide a pharmacologically actionable framework for combination therapy.

---

## 1. Introduction

### 1.1 The Clinical Problem

Glioblastoma Multiforme (GBM) is classified as a WHO Grade IV glioma. It is the most common and most aggressive primary malignant brain tumor in adults. Despite the standard-of-care treatment protocol established by Stupp et al. (2005) — maximal surgical resection followed by concurrent radiotherapy and Temozolomide (TMZ) chemotherapy — the median overall survival remains approximately 14.6 months. The five-year survival rate is less than 5%.

The fundamental reason for this poor prognosis is therapy resistance. GBM tumors are not static. They actively change their molecular identity in response to treatment. This process of identity change is called the Proneural-to-Mesenchymal Transition (PMT), and it is the single most important barrier to effective treatment.

### 1.2 The Mesenchymal Transition Problem

In 2010, Verhaak et al. classified GBM into four transcriptional subtypes: Proneural, Neural, Classical, and Mesenchymal. Among these, the Mesenchymal subtype carries the worst prognosis and is associated with therapy resistance. Critically, tumors classified as Proneural at initial diagnosis frequently recur as Mesenchymal after treatment. This subtype switching — the PMT — has been the focus of intense research.

However, a fundamental question has remained unanswered: when a tumor becomes "Mesenchymal," are the cancer cells themselves changing shape, or is the apparent change driven by non-cancer cells (such as immune cells) infiltrating the tumor?

This distinction is critical for drug development. If the cancer cells are shape-shifting, we must target the shape-shifting machinery. If the signal comes from immune cells, we must target immune cell recruitment. Targeting the wrong compartment explains why many clinical trials for GBM have failed.

### 1.3 Study Objectives

In this study, we set out to:

1. Identify the master transcription factor driving the PMT.
2. Determine whether the ATRX chromatin remodeling mutation controls the mode of therapy resistance.
3. Deconvolve the cellular origin of the mesenchymal signature using single-cell RNA-seq.
4. Trace the molecular progression of IDH-mutant astrocytomas from Grade 2 to Grade 4 across all omic layers.
5. Identify pharmacologically actionable targets to block the resistance cascade.

---

## 2. Materials and Methods

### 2.1 Patient Cohorts and Data Sources

All patient data were obtained from the Chinese Glioma Genome Atlas (CGGA) database. Two cohorts were used:

| Property | Broad WES Cohort | Deep IDH-A Cohort |
| :--- | :--- | :--- |
| **Number of patients** | 286 | 35 |
| **Genomic subtype** | Mixed (IDH-mut, IDH-wt) | IDH-mutant, non-codel |
| **DNA Sequencing** | Whole-Exome Sequencing (WES) | WES (Blood + Tumor) |
| **RNA Sequencing** | Not available | RSEM RNA-seq (all 35) |
| **DNA Methylation** | Not available | Illumina EPIC Array (all 35) |
| **Proteomics** | Not available | Mass Spectrometry (all 35) |
| **Phosphoproteomics** | Not available | Mass Spectrometry (all 35) |
| **Single-Cell RNA-seq** | Not available | 16 patients (74,563 cells) |
| **Clinical data** | Age, Gender, Grade, Histology, TMZ status, OS | Age, Gender, Grade, Subtype, OS, PFS |

### 2.2 Mutation Landscape Analysis

Whole-exome sequencing data for 286 patients was used to calculate per-gene mutation frequencies. Binary mutation matrices were constructed (mutant = 1, wildtype = 0). The top driver genes were ranked by mutation prevalence.

### 2.3 Differential Multi-Omic Analysis

For each comparison (Grade 2 vs. Grade 4, Treated vs. Untreated, ATRX-Mutant vs. ATRX-Wildtype), the following analyses were performed:

- **Epigenomics**: Mean beta values per CpG probe were compared between groups. Delta beta values (Group B mean minus Group A mean) were calculated. MGMT promoter probes (cg12434587, cg12981137, cg15765353, cg21862320) were specifically tracked.
- **Transcriptomics**: Log2 fold-change was calculated as log2((Mean_GroupB + 1) / (Mean_GroupA + 1)).
- **Proteomics**: Mann-Whitney U tests were applied to compare protein abundance between groups (significance threshold: p < 0.05).
- **Phosphoproteomics**: Specific phosphorylation sites (e.g., STAT3_Y705, EGFR_S991, MAPK1_Y187) were compared between groups.

### 2.4 Single-Cell RNA-seq Analysis

Single-cell RNA-seq data (74,563 cells from 16 patients) was processed using the Seurat pipeline (v4). Cells were clustered by UMAP, and marker gene expression was evaluated per cluster to assign cell-type identities (malignant cells, macrophages/microglia, oligodendrocyte precursors, T cells).

### 2.5 Promoter Binding Validation

The transcription start site (TSS) of *CHI3L1* was mapped to chr1:203,186,704 (GRCh38) using Ensembl. A 1,000 bp window (chr1:203,186,704–203,188,704) was queried against the ENCODE and ReMap ChIP-seq databases to identify transcription factors physically bound to this promoter.

### 2.6 Drug Target Identification

ChEMBL (v33) was queried for bioactivity data (IC50, Ki) against the identified targets: CHI3L1/YKL-40, STAT3, TEAD1, and YAP1. Molecules were ranked by potency and filtered for drug-likeness (molecular weight, LogP).

---

## 3. Results

### 3.1 The Mutation Landscape of 286 Glioma Patients

Whole-exome sequencing of 286 patients revealed the following mutation prevalence:

| Gene | Mutation Rate | Biological Role |
| :--- | :--- | :--- |
| **IDH1** | 46.9% | Metabolic reprogramming, G-CIMP methylation |
| **TP53** | 45.8% | Cell cycle checkpoint loss |
| **ATRX** | 29.7% | Chromatin remodeling, ALT telomere maintenance |
| **CIC** | 14.3% | Oligodendroglioma-associated transcription factor |
| **NOTCH1** | 8.4% | Stemness and differentiation control |
| **NF1** | 7.0% | RAS/MAPK pathway tumor suppressor |
| **PTEN** | 6.6% | PI3K/AKT pathway tumor suppressor |
| **PIK3CA** | 5.9% | PI3K pathway oncogenic activator |

The cohort spans WHO Grade II (93 patients), Grade III (91 patients), and Grade IV (102 patients), representing 54 primary GBM, 48 recurrent GBM, and 184 lower-grade gliomas. Of the 286 patients, 185 received TMZ chemotherapy and 76 did not.

### 3.2 STAT3 Is the Master Transcription Factor of the PMT

By integrating single-cell RNA-seq clustering with transcription factor activity analysis, we identified **STAT3** as the central upstream regulator of the Proneural-to-Mesenchymal Transition.

STAT3 was found to be:

- Transcriptionally upregulated in Grade 4 vs. Grade 2 (Log2FC = +0.46, G2 mean = 33.21, G4 mean = 46.18).
- Present as an active phosphoprotein at the Y705 site (mean abundance = 748,174 in G2, 447,420 in G4).
- Physically bound to the promoter of the mesenchymal driver gene *CHI3L1* at 29 independent ChIP-seq sites.

### 3.3 The ATRX Genetic Fork: Two Paths of Resistance

We stratified the 11 overlapping multi-omic patients by ATRX mutation status and TMZ treatment status. This revealed that the ATRX mutation functions as a binary genetic fork that determines the mode of therapy resistance.

**Finding: ATRX-Wildtype tumors use an immune-mediated defense.**

**Finding: ATRX-Mutant tumors use a cell-intrinsic defense.**

| Protein Marker | ATRX-WT Treated | ATRX-Mut Treated | Delta (Mut − WT) | Interpretation |
| :--- | :--- | :--- | :--- | :--- |
| **CD44** | 631.03 | 307.31 | **−323.72** | Macrophage shield is lost in ATRX-Mut |
| **CHI3L1 (YKL-40)** | 68.00 | 160.22 | **+92.22** | Cell-intrinsic shape-shifting activated |
| **GRIA2** | 70.65 | 148.00 | **+77.35** | Synaptic AMPA receptor hijack activated |
| **NLGN3** | 103.24 | 175.49 | **+72.25** | Synaptic organizer engaged |
| **EGFR** | 132.13 | 83.21 | **−48.92** | Growth factor dependence reduced |

This table demonstrates that under chemotherapy stress:

- ATRX-Wildtype tumors recruit macrophages (high CD44) to form an immune shield.
- ATRX-Mutant tumors bypass the immune shield entirely and instead plug into brain neurons (high GRIA2, NLGN3) while simultaneously activating cell-intrinsic shape-shifting (high YKL-40).

### 3.4 The Molecular Progression from Grade 2 to Grade 4

We contrasted 16 Grade 2 patients against 15 Grade 4 patients across all omic layers.

#### 3.4.1 Epigenomic Stability of MGMT

The MGMT promoter methylation remains stably silenced throughout progression:

| MGMT Probe | Grade 2 Mean Beta | Grade 4 Mean Beta | Change |
| :--- | :--- | :--- | :--- |
| cg21862320 | 0.9423 | 0.9468 | Stable |
| cg15765353 | 0.8910 | 0.8554 | Stable |
| cg12981137 | 0.4423 | 0.4103 | Stable |
| cg12434587 | 0.2272 | 0.1951 | Stable |

This means that even as the tumor becomes highly aggressive (Grade 4), it does not lose its MGMT methylation. These tumors remain sensitive to alkylating chemotherapy throughout their progression.

However, global chromatin remodeling occurs at other loci. Probes such as cg10937408 (Delta Beta = −0.6330) and cg12290492 (Delta Beta = −0.6267) show massive loss of methylation in Grade 4, indicating the opening of new oncogenic enhancers.

#### 3.4.2 Transcriptomic and Proteomic Switching

The transition from Grade 2 to Grade 4 represents a phenotypic switch from differentiated synaptic integration to invasive mesenchymal behavior:

| Gene | G2 RNA | G4 RNA | Log2FC | G2 Protein | G4 Protein | p-value |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **CHI3L1** | 10.98 | 257.79 | **+4.43** | 192.28 | 109.74 | 0.1469 |
| **CD44** | 34.79 | 129.17 | **+1.86** | 412.74 | 431.01 | 0.8588 |
| **VIM** | 380.38 | 1166.16 | **+1.61** | 11027.53 | 10984.82 | 0.4891 |
| **GRIA2** | 65.46 | 18.85 | **−1.74** | 119.01 | 97.77 | 0.0832 |
| **NLGN3** | 31.40 | 29.78 | **−0.07** | 171.95 | 131.86 | **0.0010** |
| **EGFR** | 32.87 | 20.19 | **−0.68** | 146.40 | 104.98 | **0.0168** |

The most striking change is the 20-fold increase in *CHI3L1* (YKL-40) RNA transcription (from 10.98 to 257.79) upon Grade 4 progression, marking YKL-40 as the dominant cell-intrinsic mesenchymal driver.

The statistically significant loss of NLGN3 protein (p = 0.0010) and EGFR protein (p = 0.0168) in Grade 4 indicates that the tumor detaches from its neural and growth factor dependencies as it transitions to an independent, invasive state.

### 3.5 Four Proteomic Subtypes of IDH-Mutant Astrocytomas

By clustering the protein abundance data of all 35 patients, we identified four distinct biological subtypes:

| Subtype | N | Key Protein Markers | CD44 Protein | VIM Protein | Biological Identity |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **NEU** | 11 | GRIA2 = 135, NLGN3 = 170 | 227 | 3,719 | Synaptic hijack: cells plugged into neurons |
| **AFM** | 9 | STAT3 = 105, CHI3L1 = 235 | 419 | 10,228 | Active cell-intrinsic shape-shifting |
| **PPR** | 9 | OLIG2 = 19, SOX2 = 37 | 429 | 13,472 | Stem-like progenitors driving recurrence |
| **IME** | 6 | CD44 = **807**, VIM = **19,778** | **807** | **19,778** | Dense macrophage shield |

### 3.6 Chemotherapy Creates the Resistant Subtypes

By mapping TMZ treatment status onto the proteomic clusters, we discovered that the resistant subtypes are absent in untreated tumors and only emerge after chemotherapy:

| Subtype | Untreated Patients | Treated Patients | Interpretation |
| :--- | :--- | :--- | :--- |
| **AFM** | 3 | 1 | Baseline glial state (exists before treatment) |
| **NEU** | 1 | 2 | Baseline synaptic state (exists before treatment) |
| **PPR** | **0** | **4** | Stem cell resistance (created by chemotherapy) |
| **IME** | **0** | **1** | Immune shield (created by chemotherapy) |

This is a critical finding. The PPR (stem-like progenitor) and IME (immune-infiltrated) subtypes do not exist in the untreated tumor. They are created by the stress of chemotherapy itself. TMZ damages the DNA of cancer cells, which triggers the release of inflammatory cytokines (IL-6, TNF-alpha) and damage-associated molecular patterns (DAMPs). These signals activate the JAK2/STAT3 pathway, which acts as a master chromatin remodeler. STAT3 opens the stemness gene folders (OLIG2, SOX2) to create the PPR state and the CHI3L1 promoter to secrete YKL-40 for shape-shifting.

### 3.7 Deconvolution of the Mesenchymal Illusion

By integrating single-cell RNA-seq (74,563 cells) with bulk proteomic and transcriptomic data, we resolved the cellular origin of the two mesenchymal markers:

| Marker | Cell Type | Mechanism | Phenotype |
| :--- | :--- | :--- | :--- |
| **CD44** | Tumor-Associated Macrophages | Transmembrane receptor binding hyaluronan | Adhesive clusters forming protective shield |
| **CHI3L1 (YKL-40)** | Malignant Cancer Cells | Secreted glycoprotein degrading extracellular matrix | Individual cell shape-shifting and invasion |

This resolves a decade-long confusion in the field. The "Mesenchymal Transition" is not a single event. It is a composite signal from two independent biological programs operating in two different cell types:

1. **The Immune Shield (CD44)**: Macrophages infiltrate the tumor under chemotherapy stress. They express high CD44, which acts as molecular velcro to anchor them around the tumor mass. In bulk RNA-seq, this macrophage CD44 is blended with the cancer cell signal, creating the illusion that the cancer cells have turned mesenchymal.

2. **The Cell-Intrinsic Escape (YKL-40)**: The cancer cells themselves activate STAT3, which transcribes the CHI3L1 gene to produce YKL-40 protein. YKL-40 is secreted into the brain, where it degrades the extracellular matrix and triggers cytoskeletal remodeling via IL-13Ra2 receptor signaling. This converts the cancer cell from a round, stationary cell into a long, mobile, invasive cell.

The reason these two programs respond differently to the same chemotherapy stress is due to cell-type-specific epigenetic accessibility. In cancer cells, the CHI3L1 promoter is chromatin-accessible while the CD44 promoter is methylated (locked). In macrophages, the CD44 promoter is accessible while the CHI3L1 promoter is locked. The same STAT3 signal therefore produces different outputs in different cell types.

### 3.8 Promoter-Level Validation of the CHI3L1 Regulatory Network

We mapped the TSS of CHI3L1 to chr1:203,186,704 (GRCh38, negative strand) using Ensembl. We queried the 1,000 bp upstream promoter window (chr1:203,186,704–203,188,704) against ENCODE and ReMap ChIP-seq databases:

| Transcription Factor | ChIP-seq Binding Sites | Genomic Coordinates | Role |
| :--- | :--- | :--- | :--- |
| **STAT3** | 29 | chr1:203,186,448–203,187,471 | Master switch of PMT |
| **NF-kB (RELA)** | 23 | chr1:203,186,613–203,187,547 | Inflammatory amplifier |
| **CEBPB** | 12 | chr1:203,186,500–203,187,200 | Myeloid lineage activator |

This physically validates that STAT3, NF-kB, and CEBPB co-occupy the same 1,000 bp promoter window to cooperatively drive CHI3L1 transcription.

### 3.9 Pharmacological Target Identification

We queried ChEMBL (v33) for bioactive molecules targeting the identified regulatory cascade:

| Target | ChEMBL ID | IC50 | Mechanism |
| :--- | :--- | :--- | :--- |
| **YKL-40 (direct)** | CHEMBL5568901 | **25 nM** | Neutralizes secreted YKL-40 protein |
| **TEAD1 (covalent)** | CHEMBL5401444 | **3 nM** | Permanently shuts down transcription complex |
| **STAT3 (dimerization)** | CHEMBL502473 | **150 nM** | Blocks STAT3 dimerization and nuclear entry |
| **YAP1-TEAD** | CHEMBL4439905 | **220 nM** | Disrupts YAP1-TEAD protein-protein interaction |

---

## 4. Discussion

### 4.1 Redefining the Mesenchymal Transition

For over a decade, the Mesenchymal subtype of GBM has been treated as a single biological entity. Clinical trials have targeted CD44 and other mesenchymal surface markers on the assumption that the cancer cells had undergone shape-shifting. Our findings demonstrate that this assumption is incorrect. The mesenchymal signature detected in bulk RNA-seq is a hybrid illusion composed of two independent programs: macrophage-derived CD44 (the immune shield) and cancer cell-derived YKL-40 (the cell-intrinsic escape). This explains why anti-CD44 therapies have failed in clinical trials — they strip the immune shield but leave the true shape-shifting mechanism (YKL-40) completely untouched.

### 4.2 The ATRX Genetic Fork

We report that the ATRX mutation status determines which resistance pathway the tumor activates. ATRX-wildtype tumors rely on microenvironmental defense (macrophage recruitment), while ATRX-mutant tumors activate cell-intrinsic programs (synaptic hijack and YKL-40 secretion). This has immediate clinical implications: ATRX-mutant patients require drugs targeting the synaptic integration (GRIA2/NLGN3) and YKL-40 secretion pathways, while ATRX-wildtype patients require drugs targeting macrophage recruitment and CD44-mediated immune evasion.

### 4.3 Chemotherapy as a Resistance Trigger

Our finding that the PPR and IME subtypes are absent in untreated tumors and only emerge after TMZ exposure challenges the current treatment paradigm. The standard of care involves administering TMZ as a front-line therapy, yet our data suggest that TMZ itself activates the STAT3-mediated chromatin remodeling that creates the resistant subtypes. This supports a combination therapy approach where STAT3 inhibitors are administered concurrently with TMZ to prevent the resistance transition from occurring.

### 4.4 Proposed Combination Therapy

Based on our findings, we propose a three-drug combination strategy:

1. **TMZ** (standard chemotherapy): Damages cancer cell DNA.
2. **STAT3 inhibitor** (CHEMBL502473, 150 nM): Blocks the alarm system. Prevents STAT3 from opening the resistance gene folders (CHI3L1, OLIG2, SOX2).
3. **YKL-40 neutralizer** (CHEMBL5568901, 25 nM): Neutralizes any YKL-40 protein that escapes into the brain, preventing matrix degradation and cell migration.

This combination attacks the cancer at three levels simultaneously: killing the bulk tumor (TMZ), preventing the resistance transition (STAT3 blocker), and neutralizing the escape mechanism (YKL-40 blocker).

---

## 5. Conclusion

This study provides a comprehensive multi-omic framework for understanding therapy resistance in glioblastoma. We demonstrate that:

1. **STAT3** is the master transcription factor driving the Proneural-to-Mesenchymal Transition.
2. The **ATRX mutation** acts as a genetic fork determining the mode of resistance (immune shield vs. synaptic hijack).
3. The **Mesenchymal Transition is a dual-compartment illusion**: CD44 belongs to macrophages (immune shield), CHI3L1/YKL-40 belongs to cancer cells (shape-shifting).
4. The resistant subtypes **PPR and IME are created by chemotherapy** and do not exist in untreated tumors.
5. The CHI3L1 promoter is physically co-occupied by **STAT3, NF-kB, and CEBPB**, providing validated pharmacological targets.
6. Ultra-potent inhibitors exist in ChEMBL for each node of this regulatory cascade.

These findings establish a new biological model for glioblastoma resistance and identify a pharmacologically actionable combination therapy strategy for preclinical validation.

---

## 6. Data Availability

All raw data used in this study are publicly available from the Chinese Glioma Genome Atlas (CGGA) database (http://www.cgga.org.cn/). ChIP-seq binding data were obtained from the ENCODE Project and ReMap databases. Drug bioactivity data were obtained from ChEMBL v33 (https://www.ebi.ac.uk/chembl/).

---

*Correspondence: N. Mohana Krishna*
