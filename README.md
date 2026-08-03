# Morphology and Connectome Complexity Analysis

Code and analysis accompanying the paper:

**Uncovering the Basis of Human Connectome Complexity: The Role of Neuronal Morphology**

**Authors:** Natali Barros-Zulaica, Daniela Egas Santander, Lida Kanari, Ying Shi, Rodrigo Perin, Maurizio Pezzoli, Ruth Benavides-Piccione, Javier DeFelipe, Christiaan P. J. de Kock, Idan Segev, Henry Markram, and Michael W. Reimann.

---

## Overview

This repository contains the code used to analyze how neuronal morphology influences the complexity of predicted synaptic connectomes in human and rodent cortical microcircuits.

The study investigates how axonal and dendritic geometries constrain potential synaptic connectivity and how these geometric constraints contribute to non-random relevant connectome structures.

The repository includes:

- data preprocessing utilities,
- connectivity and complexity analyses,
- figure-generation notebooks,
- and auxiliary scripts used in the paper.

---

## Repository structure

```text
.
├── connectivity_matrix_analysis
│   ├── Fig1B_C_D_connectivity_statistics.ipynb
│   ├── Fig2B_C_E_ConnProb_of_pruned_networks.ipynb
│   ├── Fig2D_vertical_bias_of_connprob.ipynb
│   ├── Fig3A_SurvivalProb.ipynb
│   ├── Fig3B_GlobalConnProb.ipynb
│   ├── Fig3C_4D_SimplexCounts.ipynb
│   ├── Fig4A_S3_DistanceDependent_ConnProb.ipynb
│   ├── Fig4B_S5_S6A_reciprocal_overexpression_and_sampling.ipynb
│   ├── Fig4C_S6_motives_simplex.ipynb
│   ├── Fig5A_Conn_Prob_previous_touches.ipynb
│   ├── Fig6_left_middle_ConnProb_S8C.ipynb
│   ├── Fig6_right_ConnProb_ATC.ipynb
│   ├── FigS2_compute_T_E_counts.ipynb
│   ├── FigS2_plot_T_E_counts.ipynb
│   ├── FigS3_compute_DistanceDependent_ConnProb_Bins.ipynb
│   ├── FigS4_weighted_simplices_Plot_csv.ipynb
│   ├── FigS5_sampled_connectivity_for_em.ipynb
│   ├── FigS8_Analyze_connectivity_higher_order_sources.ipynb
│   ├── FigS8_Plot_connectivity_higher_order_sources_modified.ipynb
│   ├── connection_probability_results.csv
│   ├── connectivity_higher_order_effect_final.h5
│   ├── connectivity_higher_order_interactions_final.h5
│   ├── survival_probabilites.pkl
│   ├── simplex_counts
│   │   ├── hc_rd_EtoE_sc.pkl
│   │   ├── human_EtoE_sc.pkl
│   │   ├── human_full_sc.pkl
│   │   ├── rat_EtoE_sc.pkl
│   │   └── rat_full_sc.pkl
│   ├── simplices
│   │   ├── run_simplex_counts.py
│   │   ├── plot_simplex_counts.ipynb
│   │   ├── plot_simplex_counts_with_pvalues.ipynb
│   │   └── log_simplex_counts.txt
│   └── triads
│       └── *.npy
│
├── matrix_generation_extraction
│   ├── GetConnMatBIG_Human.py
│   ├── GetConnMatBIG_Rat.py
│   ├── GetConnMatSMALL.ipynb
│   ├── Prunning_example.ipynb
│   ├── perform_virtual_sampling_experiment_from_reference.ipynb
│   ├── prunning.py
│   ├── run_human.sh
│   └── run_rat.sh
│
├── morphology_data_analysis
│   └── Fig1E_5B_S9S10S11S12.ipynb
│
└── README.md
```

### connectivity_matrix_analysis

Contains the analysis notebooks and intermediate results used to reproduce the connectivity-related figures from the manuscript.

The notebooks correspond to the main and supplementary figures of the paper and analyze:

- connection probability statistics,
- distance-dependent connectivity,
- pruning effects,
- reciprocal connections,
- higher-order connectivity motifs,
- simplex counts and higher-order structures,
- comparisons between human and rat connectomes.

Generated intermediate results include connectivity matrices, simplex counts, survival probabilities, and higher-order connectivity measures.

### matrix_generation_extraction

Contains the scripts and notebooks required to generate and extract connectivity matrices used in the analysis.

This folder includes:

- generation of human and rat connectivity matrices:
- GetConnMatBIG_Human.py
- GetConnMatBIG_Rat.py
- pruning procedures:
- prunning.py
- Prunning_example.ipynb
- virtual sampling experiments:
- perform_virtual_sampling_experiment_from_reference.ipynb
- execution scripts:
- run_human.sh
- run_rat.sh

These scripts reproduce the generation of different prunings of the networks used for connectivity analyses.

### morphology_data_analysis

Contains the analysis notebook associated with neuronal morphology reconstruction data.

The notebook:

Fig1E_5B_S9S10S11S12.ipynb

reproduces the analyses related to morphological properties of reconstructed neurons and their relationship with connectome complexity.

---

## Data

The repository may contain only a subset of the data required to reproduce all analyses.
 Large datasets and reconstructed neuronal morphologies can be found in zenodo: https://zenodo.org/records/21216278

Expected organization:

```text
├── Human_connectivity_matrics.zip
├── Rat_connectivity_matrices.zip
└── Morphological_reconstructions.zip
└── README.txt
```

---

## Citation

If you use this code, please cite the paper:

```bibtex
@article{barroszulaica2026connectome,
  title={Uncovering the Basis of Human Connectome Complexity: The Role of Neuronal Morphology},
  author={Barros-Zulaica, Natali and Egas Santander, Daniela and Kanari, Lida and Shi, Ying and Perin, Rodrigo and Pezzoli, Maurizio and Benavides-Piccione, Ruth and DeFelipe, Javier and de Kock, Christiaan P. J. and Segev, Idan and Markram, Henry and Reimann, Michael W.},
  year={2026}
}
```

---

## License

Please add the appropriate license for this repository (CCBY- v0.4).

---

## Contact

**Natali Barros-Zulaica**

GitHub: https://github.com/NataliBZ