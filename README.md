# :wave: About
This repository contains the code and data for the paper: "Predicting ecological interactions across space through pairwise integration of latent network patterns".


# :page_facing_up: Paper and citing
Kesem Abramov, Barry Biton, Geut Galai, Rami Puzis, and Shai Pilosof. **Predicting ecological interactions across space through pairwise integration of latent network patterns**. EcoEvoRxiv? (preprint) 2025. [DOI:TBD](TBD).


# Abstract:
1. Ecological communities are complex and exhibit considerable spatial variability, presenting challenges in accurately understanding these systems. A primary obstacle in ecological research is the existence of 'missing links' between species: inevitable unobserved interactions that limit our comprehension of ecological networks and their response to change. While link prediction methods have been developed to address this challenge, most approaches overlook the intrinsic spatial variability of ecological systems.
2. We introduce a flexible, spatially explicit framework based on matrix decomposition that leverages latent structural patterns to predict missing interactions and their strength, without requiring species traits or environmental data. The framework integrates information from paired auxiliary and target networks (locations) using thresholded SVD for link prediction. We applied it to plant–pollinator networks across the Canary Islands, performing pairwise predictions between locations, comparing them to within-location predictions (as a control), and quantifying how spatial variability influences predictive performance.
3. Predictions revealed that latent network structure contains substantial predictive information, with F0.5 scores consistently exceeding a random baseline (mean F0.5 = 0.67 ± 0.02 SD), while being less sensitive to interaction strength. The method enabled identifying plausible gaps in the data and producing ecologically coherent predictions. Incorporating information from auxiliary locations enhanced predictive accuracy in certain cases, but success depended on spatial context: predictions were most reliable when derived from nearby, ecologically similar locations, and declined with increasing geographic and ecological distance, consistent with a distance-decay effect.
4. We conclude that the predictability of missing links is spatially variable, reflecting both network and species-level heterogeneity. These patterns provide insights into network structure and the ecological processes shaping it, complementing trait-based approaches. While network structure offers rich predictive information, spatial context is essential for applying it effectively: ignoring spatial variability can obscure ecological signals and inflate predictive error. Our framework is computationally efficient, transferable, and readily applicable to any system with spatial or temporal replication. It can be used for a variety of ecological contexts, including island systems, fragmented landscapes, and environmental gradients, making it a practical and scalable tool for advancing link prediction in ecology.

# :computer: Code:
Instructions for running the code and reproducing the results are in the repository Wiki under "Code".

# :file_folder: Folder breakdown:
Detailed in the repository Wiki (under "Directories").

# :file_cabinet: Data:
Detailed in the repository Wiki (under "Data").

Case studies: 
Trøjelsgaard, K., Jordano, P., Carstensen, D. W., & Olesen, J. M. (2015). Geographical variation in mutualistic networks: Similarity, turnover and partner fidelity. Proceedings. Biological sciences / The Royal Society, 282 (1802), 20142925. https://doi.org/10.1098/rspb.2014.2925

Krasnov, Boris R., Sonja Matthee, Marcela Lareschi, Natalia P. Korallo‐Vinarskaya, and Maxim V. Vinarski. (2010). Co‐occurrence of Ectoparasites on Rodent Hosts: Null Model Analyses of Data from Three Continents. Oikos 119 (1): 120–28. https://doi.org/10.1111/j.1600-0706.2009.17902.x

Distance between sites: distance_between_sites_canary.csv

The data are available in the repository set up in original publications:
main analysis:
[https://datadryad.org/dataset/doi:
10.5061/dryad.76173.](https://datadryad.org/dataset/doi:10.5061/dryad.76173)

Host-parasite network analysis: 
[https://datadryad.org/dataset/doi:10.5061/dryad.d3d36.](https://datadryad.org/dataset/doi:10.5061/dryad.d3d36)

# :bar_chart: Output:
Plots in PDF or PNG format. Files are saved to `results/paper_figs/` (main and scale analyses) or `results/hp_analysis/paper_figs/` (host-parasite analysis).

---

## Main analysis (`Abramov_et_al_spatial_prediction_analysis.R`)

### Main figures
| File | Figure |
|------|--------|
| `subset_heatmap_hist_v2.pdf` | Fig. 2  (composite: panels a,b from subset analysis, panels c, d from main)|
| `island_heatmap_f05.pdf` | Fig. 2c |
| `hist_f05a_legend_bottom.pdf` | Fig. 2d |
| `missing_interactions_degree2.pdf` | Fig. 3 (composite: panels a, b, c) |
| `map_missing_links_merged.pdf` | Fig. 3a (individual panel) |
| `pie_chart.pdf` | Fig. 3b (individual panel) |
| `degree_unobserved_links.pdf` | Fig. 3c (individual panel) |
| `isl_jaccard_distance.pdf` | Fig. 4 |

### Supplementary figures
| File | Figure |
|------|--------|
| `pr_roc.pdf` | Fig. S11 |
| `roc_curve.pdf` | Fig. S12 |
| `predicted_original.png` | Fig. S13 |
| `degree_occurrence.pdf` | Fig. S14 |
| `local_degree_predicted_links.pdf` | Fig. S15 |
| `degree_binning.pdf` | Fig. S16 |
| `plant_island_degree.pdf` | Fig. S17 |
| `netdensity_f05_nnse.pdf` | Fig. S18 |
| `netsize_f05_nnse.pdf` | Fig. S19 |

---

## Scale analysis (`site_scale_spatial_prediction_analysis.R`)

| File | Figure |
|------|--------|
| `nnse_f05_scales.pdf` | Fig. S1 |
| `site_heatmap_f05.pdf` | Fig. S2 |
| `hist_f05_site.pdf` | Fig. S3 |
| `jaccard_site_f05.pdf` | Fig. S4 |
| `cor_plot_site_dif_f05.pdf` | Fig. S5 |
| `site_netsize_f05_nnse.pdf` | Fig. S6 |
| `site_netdensity_f05_nnse.pdf` | Fig. S7 |

---

## k sensitivity analysis (`k_swap.R`)

| File | Figure |
|------|--------|
| `k_overall_sensitivity_v2.pdf` | Fig. S8 |

---

## Alternative link withholding strategies (`alternative_link_withholding_strategies.R`)

| File | Figure |
|------|--------|
| `fig_imbalance_pr_roc.pdf` | Fig. S9 |
| `island_heatmap_degree_holdout_combined.pdf` | Fig. S10 |

---

## Temporal host-parasite analysis (`temporal_analysis_clean_version.R`)

| File | Figure |
|------|--------|
| `pr_roc_host_parasite.pdf` | Fig. S20 |
| `map_missing_links_host_parasite.pdf` | Fig. S21 |
| `year_jaccard_distance.pdf` | Fig. S22 |

---

## Subset analysis (`subset_analysis_comparison.R`)

| File | Figure |
|------|--------|
| `island_subset_f05.pdf` | Fig. 2a |
| `island_subset_nnse.pdf` | Fig. 2b |
