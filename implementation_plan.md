# Implementation Plan: Quantum Hamilton Transition Map (Q-HTM)

This plan outlines the design and mathematical layout of the **Quantum Hamilton Transition Map (Q-HTM)**. The Q-HTM uses the unitary mechanics of real quantum computing to solve the five fundamental unknowns (temporal, spatial, clonal, multi-gene, and patient-specific gaps) in Glioblastoma Multiforme (GBM) evolution under Temozolomide (TMZ) chemotherapy.

---

## 1. How We Solve the 5 Gaps Using Quantum Mechanics

| Biological Unknown (The Gap) | Physical/Quantum Concept | Mathematical Mechanism |
| :--- | :--- | :--- |
| **1. The Temporal Gap** (When does YKL-40 start rising?) | **Continuous Unitary Evolution** | Trotter steps $\Delta t$ in $e^{-iHt}$ to interpolate exact intermediate states. |
| **2. The Spatial Gap** (Where in the tumor does it start?) | **Spatial Qubit Lattice** | Entangled neighbor-state registers representing 3D coordinates. |
| **3. The Clonal Gap** (Plasticity vs. Selection?) | **Superposition Wave Collapse** | Ratio of mutation qubits ($q_{\text{mut}}$) to expression qubits ($q_{\text{expr}}$). |
| **4. The Multi-Gene Gap** (What else is changing?) | **Multi-Qubit Encoding** | Expanded 30-qubit central dogma registers with QFT phase mapping. |
| **5. The Patient-Specific Gap** (Different paths, same end?) | **Degenerate State Paths** | Mapping different initial rotations to the same final $|11\rangle$ state. |

---

## 2. Mathematical Formalism of the Q-HTM

We define the state of the tumor as a wave function $|\Psi(t)\rangle$ on the IBM Heron QPU:

### 2.1 Qubit Registry Layout (30 Qubits)
*   **$q_0 \dots q_9$ (Mutational Status)**: Track somatic VAFs of top 10 driver mutations (e.g. *ATRX*, *TP53*, *IDH1*, *PTEN*, *EGFR*, *NF1*, *CIC*, *NOTCH1*, *RB1*, *PIK3CA*).
*   **$q_{10} \dots q_{19}$ (Epigenetic Switches)**: Track promoter CpG methylation states.
*   **$q_{20} \dots q_{29}$ (Multi-Gene Phenotype)**: Track expression of *CHI3L1*, *CD44*, *VIM*, *GRIA2*, *NLGN3*, and upstream regulators (STAT3, YAP1).

### 2.2 Solving the Gaps (Mathematical Formulation)

#### Gap 1 & 2: The Spatio-Temporal Transition Hamiltonian ($\hat{H}_{\text{transition}}$)
To track *when* and *where* resistance starts, we model the chemotherapy selective pressure as a time-dependent Hamiltonian:
$$\hat{H}_{\text{transition}} = \sum_{j} \omega_j \sigma_z^{(j)} + \sum_{j,k} J_{jk} \left( \sigma_x^{(j)}\sigma_x^{(k)} + \sigma_y^{(j)}\sigma_y^{(k)} \right)$$
*   **The Drift**: $\omega_j$ represents the internal epigenetic plasticity rate of gene $j$.
*   **The Spatial coupling**: $J_{jk}$ represents the physical interaction/diffusion rate of signal between cell coordinates $j$ and $k$.
By running Trotterized steps on the QPU, we measure the state at $t = \epsilon, 2\epsilon, 3\epsilon$, charting the exact curve of YKL-40 activation over time.

#### Gap 3: The Clonal vs. Plasticity Proof (Entanglement Entropy)
To discover if resistance is driven by selection (genetics) or plasticity (epigenetics), we measure the **Quantum Entanglement Entropy** ($S_E$) between the Mutation qubits ($A$) and the Expression qubits ($B$) after drug evolution:
$$S_E(\rho_A) = -\text{Tr}(\rho_A \log_2 \rho_A) \quad \text{where} \quad \rho_A = \text{Tr}_B(|\Psi\rangle\langle\Psi|)$$
*   **If $S_E \approx 0$**: The mutations and expression are unentangled. Resistance is driven by **pure selection** of pre-existing clones.
*   **If $S_E > 0$**: High entanglement. Resistance is driven by **active cell-intrinsic plasticity** induced by the microenvironment.

#### Gap 5: The Degeneracy Proof (Different Paths, Same Result)
Do different patients take different routes to get to high YKL-40?
*   We initialize different starting configurations representing different patient genotypes: $|\Psi_{P1}(0)\rangle$ and $|\Psi_{P2}(0)\rangle$.
*   We evolve them under the same $U_{\text{TMZ}}(t)$.
*   We measure the overlap (Fidelity) of the final states: $F = |\langle\Psi_{P1}(t) | \Psi_{P2}(t)\rangle|^2$.
*   If $F \approx 1$ but the intermediate paths diverged, we prove mathematically that the tumor has **convergent evolution paths** that funnel different patients into the same resistant state.

---

## 3. Verification & Execution Strategy

1.  **Calibration**: Load baseline mutation frequencies (VAFs) and methylation beta values from the primary (`TP`) samples of the GLASS cohort.
2.  **Execution on ibm_marrakesh**: Compile the Trotter steps for $\hat{H}_{\text{transition}}$ to measure intermediate state amplitudes.
3.  **Output**: Generate the **Tumor Evolution Map (Q-TEM)** showing:
    *   The transition curve (YKL-40 probability vs. time).
    *   The Entanglement Entropy score (defining the selection vs. plasticity ratio).
    *   The path fidelity score (proving convergence).
