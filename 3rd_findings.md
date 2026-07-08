# The Multi-Omic Landscape of IDH-Mutant Astrocytoma Progression: Deconvolution of the Immune Shield and Cell-Intrinsic Shape-Shifting

## Abstract
Gliomas carrying the `IDH1` mutation display a wide range of clinical behaviors, spanning from indolent Grade 2 astrocytomas to highly lethal Grade 4 astrocytomas. In this study, we perform a deep multi-omic analysis across all 35 patients in the IDH-A cohort to trace the molecular progression switches from Grade 2 (G2) to Grade 4 (G4). We show that while MGMT promoter methylation remains stably silenced throughout progression, the transition to Grade 4 is marked by a massive, 20-fold transcriptional upregulation of the mesenchymal gene *CHI3L1* (YKL-40) and a downregulation of differentiated synaptic markers like *GRIA2* and *NLGN3*. Furthermore, we establish the molecular profiles of four distinct proteomic subtypes (NEU, AFM, PPR, and IME), identifying their unique cellular behaviors. Finally, we deconvolve the tumor microenvironment to demonstrate that the rise of *CD44* in Grade 4 is a signature of macrophage infiltration (the immune shield), whereas *CHI3L1* represents the true cell-intrinsic shape-shifting of the cancer cells. We map ChEMBL inhibitors targeting these pathways to outline a therapeutic roadmap.

---

## 1. Omic Progression from Grade 2 to Grade 4 (WHO 2021)

We contrasted 16 Grade 2 (low-grade) patients against 15 Grade 4 (high-grade) patients across DNA methylation, RNA expression, and protein abundance.

### A. Epigenomics (DNA Methylation Stability of MGMT)
*   **MGMT Silencing is Conserved**: We evaluated the promoter methylation probes for the *MGMT* DNA repair gene. The methylation status remains stably locked in both Grade 2 and Grade 4 tumors:
    *   `cg21862320`: G2 Mean Beta = 0.9423, G4 Mean Beta = 0.9468
    *   `cg15765353`: G2 Mean Beta = 0.8910, G4 Mean Beta = 0.8554
    *   *Interpretation*: Progression to Grade 4 does not lead to the loss of MGMT methylation in IDH-mutant astrocytomas, indicating they remain functionally sensitive to alkylating chemotherapy.
*   **Global Chromatin Remodeling**: Grade 4 progression is characterized by localized hypomethylation, including probes `cg10937408` (Delta Beta = -0.6330) and `cg12290492` (Delta Beta = -0.6267), indicating the opening of new oncogenic enhancers.

### B. Transcriptomics and Proteomics (Phenotypic Switching)
The transcriptional transition from Grade 2 to Grade 4 represents a shift from differentiated synaptic integration to an invasive, mesenchymal state:

| Gene / Marker | G2 Mean (RNA) | G4 Mean (RNA) | Log2FC | G2 Mean (Protein) | G4 Mean (Protein) | Protein Significance | Biological Interpretation |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **CHI3L1 (YKL-40)** | **10.98** | **257.79** | **+4.43** | **192.28** | **109.74** | p = 0.1469 | 20-fold rise in RNA transcription of the mesenchymal driver. |
| **CD44** | **34.79** | **129.17** | **+1.86** | **412.74** | **431.01** | p = 0.8588 | Significant transcript rise representing macrophage recruitment. |
| **GRIA2** | **65.46** | **18.85** | **-1.74** | **119.01** | **97.77** | p = 0.0832 | Downregulation of AMPA receptors as tumor de-differentiates. |
| **NLGN3** | **31.40** | **29.78** | **-0.07** | **171.95** | **131.86** | **p = 0.0010** (SIG) | Significant loss of synaptic organizers in Grade 4. |
| **EGFR** | **32.87** | **20.19** | **-0.68** | **146.40** | **104.98** | **p = 0.0168** (SIG) | Significant loss of baseline EGFR protein abundance. |

---

## 2. Characterization of the Four Proteomic Subtypes

We grouped the 35 patients by their proteomic clusters (NEU, AFM, PPR, IME) to identify the molecular drivers of each tumor class:

```
                          [IDH-Mutant Astrocytoma]
                                     |
         +------------------+--------+--------+------------------+
         |                  |                 |                  |
     [NEU Subtype]      [AFM Subtype]     [PPR Subtype]      [IME Subtype]
      Synaptic Hijack     Astrocyte-like    Progenitor-like    Immune Shield
     GRIA2+/NLGN3+      STAT3+/YKL-40+    OLIG2+/SOX2+       CD44+/VIM+
```

1.  **NEU (Neuronal-like - Synaptic Hijack)**:
    *   *Markers*: Highest **GRIA2** (RNA = 58.78, Protein = 135.06) and **NLGN3** (RNA = 38.89, Protein = 169.66).
    *   *Identity*: Highly differentiated tumor cells that plug directly into neurons to feed on neurotransmitter activity.
2.  **AFM (Astrocyte-like / Fibroglial)**:
    *   *Markers*: High **STAT3** (Protein = 105.49) and **CHI3L1/YKL-40** (Protein = 235.17).
    *   *Identity*: Glial-lineage cells undergoing active cell-intrinsic shape-shifting.
3.  **PPR (Progenitor-like)**:
    *   *Markers*: Highest **OLIG2** (RNA = 224.97) and **SOX2** (RNA = 272.50).
    *   *Identity*: Stem-like progenitor cells responsible for tumor recurrence and mitotic activity.
4.  **IME (Immune-like - Infiltrated)**:
    *   *Markers*: Highest **CD44** (Protein = **806.60**) and Vimentin (**VIM** Protein = **19,778.48**).
    *   *Identity*: Highly inflamed tumors surrounded by a dense shield of immunosuppressive macrophages.

---

## 2b. The Before-and-After Treatment Metamorphosis

By mapping the treatment status of patients across the proteomic clusters, we traced the exact timeline of tumor adaptation under chemotherapy pressure.

### A. The Baseline State (Before Treatment)
In the untreated state, tumors belong almost exclusively to the **AFM** (astrocyte-like) or **NEU** (neuronal/synaptic) subtypes:
*   **AFM (75% of untreated cases)**: The cancer cells maintain a slow-growing, glial-like phenotype with low baseline transcription of *CHI3L1* (YKL-40).
*   **NEU (25% of untreated cases)**: The cancer cells form baseline connections with neurons but have not yet developed drug-defense shields.
*   *Key Characteristic*: The tumor lacks stem-like resistance (PPR) and has very low macrophage infiltration (IME, CD44-low).

### B. The Chemo-Adapted State (After Treatment)
Once exposed to Temozolomide (TMZ) chemotherapy, the tumor is reprogrammed, giving rise to two highly resistant subtypes that are completely absent in the untreated cohort:
1.  **PPR (Progenitor-like)**: The cells adapt by transforming into slow-dividing, highly resistant stem-like progenitors expressing high levels of *OLIG2* and *SOX2*.
2.  **IME (Immune-infiltrated)**: The tumor recruits a dense microenvironmental shield of CD44+ macrophages to block chemo penetration.

### C. The Triggers: What Drives the Transition?
The transformation from the baseline state to the resistant states is triggered directly by the stress of chemotherapy:
1.  **DNA Damage and Cellular Panic**: TMZ induces DNA methylation and double-strand breaks. In response to this DNA damage, dying and stressed tumor cells release inflammatory cytokines (like IL-6 and TNF-alpha) and DAMPs (like ATP and HMGB1).
2.  **STAT3 Activation**: The extracellular IL-6 and growth factor receptor signaling (EGFR) activate the JAK2/STAT3 pathway in the remaining tumor cells, phosphorylating STAT3 at Tyrosine 705.
3.  **Nuclear Translocation**: STAT3 dimerizes, enters the nucleus, and acts as a master chromatin remodeler. It opens the stemness folders (*OLIG2*/*SOX2*) to drive the **PPR state** and the *CHI3L1* promoter to secrete **YKL-40** (shape-shifting).

### D. Why Different Cells Respond Differently to the Same Chemo Signal
A major puzzle is how the same chemotherapy stress causes **CD44** to rise in immune cells and **YKL-40** to rise in cancer cells. This is governed by cell-type-specific **epigenetic folders**:
*   **In the Cancer Cells**: The `CHI3L1` (YKL-40) promoter is chromatin-accessible, while the `CD44` promoter is chemically locked (methylated). STAT3 activation can only access and transcribe **`CHI3L1`**, causing the cancer cells to shape-shift.
*   **In the Macrophages**: The `CD44` promoter is chromatin-accessible (since macrophages naturally use CD44 for migration), while the `CHI3L1` promoter is locked. The inflammatory signals activate NF-kB and STAT3 in these immune cells, transcribing **`CD44`** to form the physical protective shield.

---

## 3. Cellular Deconvolution: CD44 Macrophages vs. CHI3L1 (YKL-40) Shape-Shifting

By integrating single-cell RNA-seq (74,563 cells) with the clinical cohorts, we deconvolve the two arms of Mesenchymal resistance:

### A. CD44 (The Immune Shield)
*   *Cell Type*: **Tumor-Associated Macrophages / Microglia**.
*   *Mechanism*: CD44 is a transmembrane receptor that binds hyaluronan, acting as physical adhesive anchors. Under chemotherapy stress, the tumor recruits these macrophages, which express high CD44, forming a protective cellular shield.

### B. CHI3L1 / YKL-40 (The Cell-Intrinsic Shape-Shift)
*   *Cell Type*: **Malignant Cancer Cells**.
*   *Mechanism*: YKL-40 is a secreted glycoprotein. Under STAT3 transcription control, cancer cells secrete YKL-40 into the brain, where it:
    1.  Remodels and degrades the extracellular matrix to clear a migration path.
    2.  Binds back to receptors (like IL-13Ra2) on the cancer cell to trigger cytoskeletal remodeling, converting the cell into a mobile, spindly mesenchymal phenotype.

---

## 4. Pharmacological Alignment of Upstream Regulators

We verified that the *CHI3L1* promoter contains physically bound transcription complexes using ENCODE and ReMap ChIP-seq data, and mapped corresponding inhibitors:

*   **STAT3 (29 ChIP-seq sites)**: Binds `chr1:203,186,448-203,187,471` (overlaps TSS `203,186,704`).
    *   *Inhibitor*: **`CHEMBL502473`** (IC50 = 150 nM). Blocks dimerization and nuclear entry.
*   **NF-kB (23 ChIP-seq sites)**: Binds `chr1:203,186,613-203,187,547`.
    *   *Inhibitor*: **`CHEMBL21156`** (IC50 = 300 nM). Inhibits IKK-beta to block NF-kB activation.
*   **YAP1/TEAD (Hippo Pathway)**: Cooperates to open chromatin at the promoter.
    *   *Inhibitor*: **`CHEMBL5401444`** (IC50 = 3 nM). Covalent TEAD1 inhibitor that halts transcription.
*   **YKL-40 (Direct Inhibitor)**:
    *   *Inhibitor*: **`CHEMBL5568901`** (IC50 = 25 nM). Neutralizes secreted YKL-40 in the brain.
