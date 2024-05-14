# PRAGMA_treatment-patterns

This repository builds on data prepared in another repository of the PRAGMA project:
https://github.com/jakobmanthey/PRAGMA_treatment-patterns

In this repository, data from patients with an incident alcohol use disorder diagnosis is analysed. At the core of this repository is the identification of treatment patterns reflecting seven different alcohol-specific treatment types through latent class analyses (LCA). The groups are subsequently verified and identified via visualizations of intervention compositions, frequencies and crossovers within classes. As a second step, a hierarchical cluster analysis is computed to evaluate robustness of the identified groups by the lca. 

The main analysis and all computations are performed within "AnalyzePragmaDataQuartersAUD_NoRelTime.qmd".
This quarto script relies on the helper / utility functions contained in "upsample_daily_utils.R". 

File paths are relational (here package), thus the local folder/project structure should equal the github repository. 

