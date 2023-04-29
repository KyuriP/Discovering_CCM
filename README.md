[![DOI](https://zenodo.org/badge/576782527.svg)](https://zenodo.org/badge/latestdoi/576782527) [![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/) [![Active](http://img.shields.io/badge/Status-Active-green.svg)](https://github.com/KyuriP/Thesis_KP)  

***MSBBSS Thesis Project***
<div align="left">
<img src="https://github.com/KyuriP/Thesis_KP/blob/main/cyclelogo.png" width=10% height=10% align="left">
<h1> Discovering Cyclic Causal Models in Psychological Research </h1>
</div>



This repository serves as a research archive for the thesis project ***Discovering Cyclic Causal Models in Psychological Research***. 
It contains the `R` code, data, figures, manuscript, supplementary material, and other subsidiary files related to the project.

<details>
           <summary><b><i>Abstract</i></b></summary>
           <p> 
Statistical network models have become popular tools for analyzing multivariate psychological data. 
In empirical practice, network parameters are often interpreted as reflecting causal relationships – an approach that can be characterized as a form of causal discovery. 
Recent research has shown that undirected network models are likely to perform poorly as causal discovery tools in the context of discovering acyclic causal structures, a task for which many alternative methods are available. 
However, acyclic causal models are likely unsuitable for many psychological phenomena, such as psychopathologies, which are often characterized by cycles or feedback loop relationships between symptoms. 
A number of cyclic causal discovery methods have been developed, largely in the computer science literature, but they are not as well studied or widely applied in empirical practice. 
In this paper, we provide an accessible introduction to the basics of cyclic causal discovery for empirical researchers. 
We examine three different cyclic causal discovery methods and investigate their performance in typical psychological research contexts by means of a simulation study. We also demonstrate the practical applicability of these methods using an empirical example and conclude the paper with a discussion of how the insights we gain from cyclic causal discovery relate to statistical network analysis.
           </p>
         </details>


## Structure
```
├── Thesis_KP.Rproj
├── manuscript
│   ├── thesis.tex
│   ├── ...
│   └── kyuriapa.bst
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
│   ├── 02_stability_analysis.R
│   └── data
├── utils
│    ├── CCD_fnc.R
│    ├── data_generating_fnc.R
│    ├── eval_metrics_fnc.R
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

## Description

| Folder                                     | Contents                                                         |
| :----------------------------------------- | :--------------------------------------------------------------- |
| [`manuscript`](https://github.com/KyuriP/CCM_Discovery/tree/main/%20%20%20%20%20manuscript)               | It contains all relevant files associated with the main manuscript. <br> • `thesis.pdf?`: the deliverable manuscript. <br> • `thesis_draft.tex`:  a plaintext file that contains the source code to generate `thesis.pdf??`.  <br> • `references.bib`: a bibliography file that contains all the references used in the manuscript. <br> • `kyuriapa.bst`: a BibTeX style file used to render [APA](https://apastyle.apa.org/style-grammar-guidelines/paper-format) formatting and style of references. |
| [`simulation`](https://github.com/KyuriP/CCM_Discovery/tree/main/%20%20%20%20simulation)       | It contains all `R` scripts required to run the simulation study.  <br> • `01_simulation.R`: simulate models, generate data, and run the algorithms. <br> • `02_simulation_results.R`: evaluate performance of algorithms and create figures. <br> • `03_sensitivity_analysis1.R`: run the secondary simulation with randomly sampled coefficients for the regression matrix $\mathbf{B}$. <br> • `04_sensitivity_analysis2.R`: conduct the secondary simulation by randomly sampling coefficients with *positive* values for the regression matrix $\mathbf{B}$. <br> • `05_sensitivity_analysis3.R`: run secondary analysis  with varying $\alpha$ levels. <br> • `06_extra_analysis.R`: investigate the unexpected patterns in the *5-variable dense* conditions from the main simulation study.<br> • `data`: output objects from all simulations are stored as `Rdata` files due to considerable processing time involved.      |
| [`empirical_example`](https://github.com/KyuriP/CCM_Discovery/tree/main/%20%20%20empirical_example) | It contains all R scripts required to run the empirical analysis. <br> • `01_empirical_analysis.R`: run empirical analysis (apply the CCD, FCI and CCI algorithms to an empirical data from [McNally et al. (2017)](https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900#supplementary-materials). <br> • `02_stability_analysis.R`: perform the stability analysis on the empirical data (*Appendix I* in the paper). <br> • `data`: raw data that is publicly available [here](https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900#supplementary-materials) and output objects from the sensitivity analysis saved as `Rdata` files.               |
| [`utils`](https://github.com/KyuriP/CCM_Discovery/tree/main/%20%20utils)               | It includes all the supporting functions required to run the other `R` scripts for conducting analyses. |
| [`supplementary_material`](./supplementary_material)           |It contains all files related to the supplementary material of the paper. <br> The supplementary material can be accessed via this [link](https://kyurip.quarto.pub/discovering-cyclic-causal-models/).|
| [`figures`](./figures)               | It contains all the figures that are presented in the paper.                    |
| [`others`](./others)               | It contains all types of subsidiary files  created during the course of the project (e.g., [research proposal](./others/Proposal), [presentation](./others/Presentation), [interim report](https://github.com/KyuriP/CCM_Discovery/tree/main/others/Research%20Report).                        |

## Reproducing study
If interested in reproducing the results of:
- simulation study, please refer to [`guide_simuation`](https://github.com/KyuriP/CCM_Discovery/blob/main/%20%20%20%20simulation/README.md).
- empirical anlaysis, please refer to [`guide_empirical`](https://github.com/KyuriP/CCM_Discovery/blob/main/%20%20%20empirical_example/README.md).


## Ethics & Access
All data are simulated in this study except for the one used in the emprical analysis, which is publicly available on the [Psychological Medicine Journal](https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900\#supplementary-materials) webpage.
Ethical approval for this study has been granted by the Ethical Review Board of the Faculty of Social and
Behavioural Sciences [(FETC)](https://ferb.sites.uu.nl/) at Utrecht University, with reference numbers 22-1810 and 22-1825 (see the confirmation of approval: [fetcapproval1.pdf](https://github.com/KyuriP/Thesis_KP/blob/main/others/FETC-approval/fetcapproval1.pdf), [fetcapproval2.pdf](https://github.com/KyuriP/Thesis_KP/blob/main/others/FETC-approval/fetcapproval2.pdf)).

This repository, which serves as an archive of the study is accessible to the public on [GitHub](https://github.com/KyuriP/Thesis_KP) under the license type of `GNU General Public License v3.0`. 


## Contact
The public accessibility of the project repository is maintained by Kyuri Park.
For any inquiries or feedback about the project, please feel free to contact [Kyuri Park](https://kyurip.github.io/).
