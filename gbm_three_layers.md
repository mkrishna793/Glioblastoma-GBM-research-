# 🧬 Glioblastoma (GBM): The Three-Layer Research Framework

This document verifies and expands on the 3-layer biological framework of Glioblastoma (GBM) research. 

---

## 🔍 Scientific Verification Summary
**Verdict: The framework is biologically and conceptually CORRECT.** 
Your analogies (Genomics = Hardware, Epigenomics = Operating System/Software, Transcriptomics = Real-time Output) are highly accurate and serve as an excellent foundation for bioinformatics research. 

Below is the verified breakdown with additional technical details to ground the data analysis.

---

## 💾 Layer 1: Genomics (The Hardwired Corrupted Code)
* **Analogy**: Hardware. Permanent structural changes to the DNA sequence (A, T, C, G). Once mutated, these changes remain in the tumor cell lineage.

### Verified Elements:
1. **The Broken Brakes (Deletions)**: 
   * *The Data*: Deletion of the [PTEN](https://www.uniprot.org/uniprotkb/P60484/entry) gene (typically via loss of Chromosome 10) and mutations/deletions in [TP53](https://www.uniprot.org/uniprotkb/P04637/entry).
   * *Function*: PTEN normally acts as a brake on the PI3K cell-growth pathway. Without it, growth signals run unchecked. TP53 is the "guardian of the genome"—without it, damaged cells refuse to undergo apoptosis (suicide).
2. **The Glued Gas Pedal (Amplification)**:
   * *The Data*: Focal amplification of the [EGFR](https://www.uniprot.org/uniprotkb/P00533/entry) (Epidermal Growth Factor Receptor) gene, often accompanied by the `EGFRvIII` mutation (deletion of exons 2-7).
   * *Function*: Leads to massive overexpression of growth receptors on the cell surface, flooding the cell with growth signals even in the absence of external growth factors.
3. **The Core Dark Matter Hit (TERT Promoter Mutation)**:
   * *The Data*: Single-nucleotide mutations in the promoter region of the [TERT](https://www.uniprot.org/uniprotkb/O14746/entry) gene (most commonly `C228T` and `C250T`).
   * *Function*: Creates a new binding site for GABP/Ets transcription factors, causing the cell to over-express telomerase. This preserves telomeres (chromosome caps) and renders the tumor cells replicatively immortal.

---

## 💿 Layer 2: Epigenomics (The 3D Packaging Software)
* **Analogy**: The Operating System. It does not change the DNA letters, but controls which parts of the DNA blueprint are wound tightly (silenced) or uncoiled (accessible).

### Verified Elements:
1. **The Spool Network (Chromatin Packaging)**:
   * *The Data*: Histone modifications (like acetylation and methylation) and DNA methylation.
   * *Function*: DNA is wrapped around histone protein "spools." Adding acetyl tags (via Histone Acetyltransferases - HATs) loosens the packaging, allowing reading machinery to access genes (Active/Open Chromatin). Removing acetyl tags or adding methyl tags (via HDACs and Methyltransferases) locks the DNA tight, silencing the genes (Silent/Closed Chromatin).
2. **The MGMT Epigenetic Switch (Crucial GBM Specific)**:
   * *The Data*: Promoter methylation of the [MGMT](https://www.uniprot.org/uniprotkb/P16455/entry) gene.
   * *Function*: If the MGMT promoter is methylated (silenced), the tumor cell cannot repair the DNA damage caused by Temozolomide (TMZ) chemotherapy, making the patient much more responsive to treatment.
3. **Epigenetic Subtype Plasticity**:
   * *The Data*: Chromatin remodeling under treatment.
   * *Function*: Under the stress of radiation or chemotherapy, chromatin remodeling complexes shift the active DNA loops. Proneural chromatin regions are packed away, and Mesenchymal chromatin regions are uncoiled, resetting the cell's baseline state.

---

## 📠 Layer 3: Transcriptomics (The Real-Time Live Orders)
* **Analogy**: The Live Output Messages. It counts the temporary mRNA copies currently circulating in the cell.

### Verified Elements:
1. **The Expression Count**:
   * *The Data*: Measured in TPM (Transcripts Per Million) or FPKM from RNA-Seq.
   * *Function*: Reflects which genes are actively being printed into proteins at any given moment.
2. **The Proneural Subtype Signature**:
   * *The Data*: High expression of marker genes like [OLIG2](https://www.uniprot.org/uniprotkb/Q9H6Z4/entry), [SOX2](https://www.uniprot.org/uniprotkb/P48431/entry), and `ASCL1`.
   * *Function*: Directs cells to mimic rapid, early-stage neural development.
3. **The Mesenchymal Subtype Signature (The Post-Chemo Shift)**:
   * *The Data*: Spike in transcripts of mesenchymal markers like [CD44](https://www.uniprot.org/uniprotkb/P16070/entry), [VIM](https://www.uniprot.org/uniprotkb/P08670/entry) (Vimentin), and `CHI3L1` (YKL-40).
   * *Function*: Programmed for survival, migration, invasiveness, and inflammatory signaling, which are highly resistant to traditional therapy.
