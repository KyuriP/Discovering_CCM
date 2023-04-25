[![DOI](https://zenodo.org/badge/576782527.svg)](https://zenodo.org/badge/latestdoi/576782527)

***MSBBSS Thesis Project***
<div align="left">
<img src="https://github.com/KyuriP/Thesis_KP/blob/main/cyclelogo.png" width=10% height=10% align="left">
<h1> Discovering Cyclic Causal Models in Psychological Research </h1>
</div>
<br>

This repository serves as an archive for the thesis project ***Discovering Cyclic Causal Models in Psychological Research***. 
It contains the `R` code, data, figures, manuscript, supplementary material, and other relevant files related to the project.

## Structure
```
├── Thesis_KP.Rproj
├── simulation
│   ├── 01_simulation.R
│   ├── 02_simulation_results.R
│   ├── 03_sensitivity_analysis1.R
│   ├── 04_sensitivity_analysis2.R
│   ├── 05_sensitivity_analysis3.R
│   ├── 06_extra_analysis.R
│   └── data
├── empirical_example
│   ├── 01_empircal_analysis.R
│   ├── 02_sensitivity_analysis.R
│   └── data
├── utils
│    ├── CCD_fnc.R
│    ├── data_generating_fnc.R
│    ├── eval_metrics.R
│    ├── plot_fnc.R
│    ├── searchAM_KP_fnc.R
│    └── true_ancestral_fnc.R
├── supplementary_material
│   ├── Supplementary_material.html
│   ├── ...
│   └── Supplementary_material.qmd
├── figures
│   ├── Fig1_DAG_DCG_GGM.pdf
│   ├── ...
│   └── FigJ1_McNallydistribution.pdf
└── others
    ├── FETC-approval
    ├── Presentation
    ├── Proposal
    └── Research Report
```

## Descrition

| Folder                                     | Contents                                                         |
| :----------------------------------------- | :--------------------------------------------------------------- |
| [`simulation`](./simulation)       | It contains all `R` scripts required to run the simulation study.  <br> - `01_simulation.R`: simulate models, generate data, and run the algorithms. <br> - `02_simulation_results.R`: evaluate performance of algorithms and create figures. <br> - `03_sensitivity_analysis1.R`: run the secondary simulation with randomly sampled coefficients for the regression matrix B. <br> - `04_sensitivity_analysis2.R`: conduct the secondary simulation by randomly sampling coefficients with *positive* values for the regression matrix B. <br> - `05_sensitivity_analysis3.R`: run secondary analysis  with varying alpha levels. <br> - `06_extra_analysis.R`: investigate the unexpected patterns in the B5 dense conditions from the main simulation study.<br> - `data`: output objects from all simulations are stored as `Rdata` files due to considerable processing time involved.      |
| [`empirical_example`](./empirical_example) | It contains all R scripts required to run the empirical analysis. <br> - `01_empirical_analysis.R`: run empirical analysis (apply the CCD, FCI and CCI algorithms to an empirical data from [McNally et al. (2017)](https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900#supplementary-materials). <br> - `02_sensitivity_analysis.R`: perform the stability analysis on the empirical data (*Appendix I* in the paper). <br> - `data`: raw data that is publicly available [here](https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900#supplementary-materials) and output objects from the sensitivity analysis saved as `Rdata` files.               |
| [`utils`](./utils)               | It includes all the supporting functions required to run the `R` scripts mentioned above. |
| [`manuscript`](./manuscript)               | It contains all relevant files associated with the main manuscript. |
| [`supplementary_material`](./supplementary_material)           |It contains all files related to the supplementary material of the paper. The supplementary material can be accessed via this [link](https://kyurip.quarto.pub/discovering-cyclic-causal-models/).|
| [`figures`](./figures)               | It contains all the figures that are included in the paper.                    |
| [`others`](./others)               | It contains all types of subsidiary files  created during the course of the project (e.g., research proposal, presentation, interim report).                        |


## Ethics & Access
All data are simulated in this study except for the one used in the emprical analysis, which is publicly available on the [Psychological Medicine Journal](https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900\#supplementary-materials) webpage.
Ethical approval for this study has been granted by the Ethical Review Board of the Faculty of Social and
Behavioural Sciences [(FETC)](https://ferb.sites.uu.nl/) at Utrecht University, with reference numbers 22-1810 and 22-1825 (see the confirmation of approval: [fetcapproval1.pdf](https://github.com/KyuriP/Thesis_KP/blob/main/others/FETC-approval/fetcapproval1.pdf), [fetcapproval2.pdf](https://github.com/KyuriP/Thesis_KP/blob/main/others/FETC-approval/fetcapproval2.pdf)).

This repository, which serves as an archive of the study is accessible to the public on [GitHub](https://github.com/KyuriP/Thesis_KP) under the license type of `GNU General Public License v3.0`. 


## Contact
The public accessibility of the project repository is maintained by Kyuri Park.
For any inquiries or feedback about the project, please feel free to contact [Kyuri Park](https://kyurip.github.io/).
