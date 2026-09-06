# A Cluster-Specific First-principles Network Pharmacology Framework for Molecular-Level Mechanism Deduction

[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3100/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.64898%2F2026.08.01.742237-blue)](https://doi.org/10.64898/2026.08.01.742237)

**Independent Drug Discovery Lab (IDD Lab)** 
**Authors:** Thang T. Dang, Viet H. Pham, Ngan T. T. Nguyen, Phong X. Nguyen, Duc M. Trinh[cite: 1]

This repository contains the official code, methodology scripts, and validation data for the research article: *"A Cluster-Specific First-principles Network Pharmacology Framework for Molecular-Level Mechanism Deduction: Application to the HL-60-Selective Cytotoxicity of 3-Deoxycardiobutanolide"*[cite: 1].

## 📖 Abstract
Standard network pharmacology workflows relying on bulk pathway enrichment frequently produce broad, associative terms rather than molecular-resolution, testable mechanisms[cite: 1]. To address this, we introduce a network pharmacology framework designed to propose molecular-level mechanistic hypotheses using a cluster-specific protein-protein interaction (PPI) network expansion strategy and a first-principles deduction protocol[cite: 1]. 

By explicitly mapping the direct consequences of partial node inhibition—substrate accumulation, product depletion, and feedback disruption—before introducing cell-line-specific transcriptomic and dependency data, the architecture separates mechanistic reasoning from contextualization[cite: 1]. This repository hosts the full computational workflow, demonstrated on 3-deoxycardiobutanolide (Compound 2), to generate falsifiable, node-resolved hypotheses regarding its HL-60 leukemic selectivity and Bax/Bcl-2 anomalies[cite: 1].

## ⚙️ System Requirements & External Dependencies

The framework leverages several standalone software packages for physics-based simulations, docking, and network visualization that must be installed separately:

* **AutoDock Vina (v1.2.5)** - Utilized for the primary reverse docking library screen[cite: 1].
* **AutoDock4** - Required for orthogonal in-depth docking validation and Inter-Molecular Energy (IME) calculation[cite: 1].
* **GROMACS (v2026.1)** - Required for explicit-solvent molecular dynamics production and MM/PBSA thermodynamic estimation[cite: 1].
* **Cytoscape (with STRING app)** - Utilized for PPI network architecture clustering and visualization[cite: 1].
* **fpocket** - Utilized for binding site and pocket prediction[cite: 1].

## 📂 Repository Structure

The repository is organized into four primary modules corresponding to the methodological pipeline:

### `1_Target_Dataset_Construction/`
Contains the foundational datasets and generation scripts for the dual-strategy target library[cite: 1].
* **Group A (Exploratory):** Scripts and data utilizing ChEMBL Target Prediction, ChEMBL Multitask Neural Network (MTNN), and SwissTargetPrediction[cite: 1].
* **Group B (Investigative):** Scripts for extracting and processing Cancer Dependency Map (DepMap) DEMETER2 RNAi dependency scores[cite: 1].
* Includes target overlap analysis and final unified target list generation.

### `2_Reverse_Docking/`
Houses the execution scripts and summarized outputs for the large-scale docking protocol.
* **Ligand:** Contains the SMILES string and prepared structures for Compound 2[cite: 1].
* **Reverse_Docking_Procedure:** Execution sequences and usage guides for the Vina screening pipeline.
* **Reverse_Docking_Analysis:** Contains chain text files, ranked docking result summaries, and the final binding grid CSVs for Groups A and B[cite: 1].

### `3_General_Biological_Framework/`
Contains the analytical core of the cluster-specific deduction protocol.
* **DepMap Data:** Raw and extracted expression profiles $(Log_2(TPM+1))$ alongside processing scripts[cite: 1].
* **PPI_Network_Construction:** Contains independent network topology data, node lists, and enrichment outputs for Clusters 1 through 4, as well as the bulk enrichment ablation study[cite: 1].
* Includes high-resolution illustrations (`BGF_Lined.png`) and the original hand-drawn conceptualization (`GBF_Hand_Written.pdf`) of the General Biological Framework[cite: 1].

### `4_Hit_Validation/`
Contains rigorous post-generation filtering and physical scrutiny data[cite: 1].
* **In-Depth and Decoy Docking:** Summaries and interaction metrics for the orthogonal AutoDock4 validation against co-crystallized reference ligands and decoy proteins[cite: 1].
* **Molecular Dynamics Simulations:** Contains MDP parameter files, initial configurations, MM/PBSA energy decomposition summaries, and trajectory analyses (RMSD/RMSF)[cite: 1]. 

## 💾 Data Availability

To ensure reproducibility, solvent-stripped, compressed coordinates (`.xtc`) for all 100 ns molecular dynamics trajectories, along with initial topologies and MDP parameter files, are publicly available via Zenodo. 

* **CYP1A1 Trajectories:** [Insert Zenodo DOI]
* **RFC4 Trajectories:** [Insert Zenodo DOI]
* **TOP2A Trajectories:** [Insert Zenodo DOI]
* **NAMPT Trajectories:** [Insert Zenodo DOI]

*Note: Due to file size constraints exceeding public repository limits (>250 GB total), the raw, explicit-solvent trajectory files (`.trr`) and full reverse docking datasets (including all decoy structure files) are available from the corresponding author upon reasonable request.*

## 📜 Citation & Usage

If you utilize the code, conceptual framework, or raw simulation data provided in this repository, please cite the bioRxiv preprint alongside any specific Zenodo DOIs associated with the datasets.

```bibtex
@article{Dang2026NetworkPharm,
  title={A Cluster-Specific First-principles Network Pharmacology Framework for Molecular-Level Mechanism Deduction: Application to the HL-60-Selective Cytotoxicity of 3-Deoxycardiobutanolide},
  author={Dang, Thang T. and Pham, Viet H. and Nguyen, Ngan T. T. and Nguyen, Phong X. and Trinh, Duc M.},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.08.01.742237}
}
