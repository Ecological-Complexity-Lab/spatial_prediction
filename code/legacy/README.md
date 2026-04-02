# :wave: About
This repository contains the code and data for the paper: "Structure knows best: predicting ecological interactions across space through pairwise integration of latent network patterns".


# :page_facing_up: Paper and citing
Kesem Abramov, Barry Biton, Geut Galai, Rami Puzis, and Shai Pilosof. **Structure knows best: predicting ecological interactions across space through pairwise integration of latent network patterns**. EcoEvoRxiv? (preprint) 2025. [DOI:TBD](TBD).


# Abstract:
1. Ecological communities are complex and exhibit considerable spatial variability, presenting challenges in accurately understanding these systems. A primary obstacle in ecological research is the existence of ‘missing links’ between species: inevitable unobserved interactions that limit our comprehension of ecological networks and their response to change. While link prediction methods have been developed to address this challenge, most approaches overlook the intrinsic spatial variability of ecological systems.
2. We introduce a flexible, spatially explicit framework based on matrix decomposition that leverages latent structural patterns to predict missing interactions and their strength, without requiring species traits or environmental data. The framework integrates information from paired auxiliary and target networks (locations) using thresholded SVD for link prediction. We applied it to plant–pollinator networks across the Canary Islands, performing pairwise predictions between locations, comparing them to within-location predictions (as a control), and quantifying how spatial variability influences predictive performance.
3. Predictions revealed that latent network structure contains substantial predictive information, with F1 scores consistently exceeding a random baseline (mean F1 = 0.68 ± 0.04 SD), while being less sensitive to interaction strength. The method enabled identifying plausible gaps in the data and producing ecologically coherent predictions. Incorporating information from auxiliary locations enhanced predictive accuracy in certain cases, but success depended critically on spatial context: predictions were most reliable when derived from nearby, ecologically similar sites or islands with larger networks, and declined with increasing geographic and ecological distance, demonstrating a clear distance-decay effect.
4. We conclude that the predictability of missing links is spatially variable, reflecting both network and species-level heterogeneity. These patterns provide insights into network structure and the ecological processes shaping it, complementing trait-based approaches. While network structure offers rich predictive information, spatial context is essential for applying it effectively: ignoring spatial variability can obscure ecological signals and inflate predictive error. Our framework is computationally efficient, transferable, and readily applicable to any system with spatial or temporal replication. It can be used for a variety of across ecological contexts, including island systems, fragmented landscapes, and environmental gradients, making it a practical and scalable tool for advancing link prediction in ecology.

# :computer: Code:
Instructions for running the code and reproducing the results are in the repository Wiki under "Code".

# :file_folder: Folder breakdown:
Detailed in the repository Wiki (under "Directories").

# :file_cabinet: Data:
Detailed in the repository Wiki (under "Data").

Case study: Trøjelsgaard, K., Jordano, P., Carstensen, D. W., & Olesen, J. M. (2015). Geographical variation in mutualistic networks: Similarity, turnover and partner fidelity. Proceedings. Biological sciences / The Royal Society,282 (1802), 20142925. https://doi.org/10.1098/rspb.2014.2925

Distance between sites: distance_between_sites_canary.csv
The data are available in the repository set up in original publication: https://datadryad.org/dataset/doi:
10.5061/dryad.76173.

# :Output:
plots is pdf format.
Main code: 

- optimal_threshold (Fig. S8)
- predicted_original (Fig. S10)
- pr_roc (Fig. S9)
- hist_f1a_legend_bottom (Fig. 2d)
- netsize_f1_nnse (Fig. S13)
- netdensity_f1_nnse (Fig. 5)
- degree_unobserved_links (Fig. 3c)
- map_missing_links (Fig. 3a)
- map_missing_links_diags_offs (Fig. S11)
- pie_chart (Fig. 3b)
- island_heatmap_f1 (Fig. 2c)
- degree_occuurrence (Fig. S12)
- isl_jaccard_distance (Fig. 4)

Site scale analysis:
- hist_f1_site (Fig. S3)
- site_netsize_f1_nnse (Fig. S6)
- site_netdensity_f1_nnse (Fig. S7)
- site_heatmap_f1 (Fig. S2)
- jaccard_site_f1 (Fig. S4)
- cor_plot_site_dif_f1 (Fig. S5)
- nnse_f1_scales (Fig. S1)

Subset analysis:
- subset_analysis (Fig 2a,b)

