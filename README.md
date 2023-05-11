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
│   ├── manuscript_KP.tex
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
| [`manuscript`](./%20%20%20%20%20manuscript)               | It contains all relevant files associated with the main manuscript. <br> • `manuscript_KP.pdf`: the finalized manuscript. <br> • `manuscript_KP.tex`:  a plaintext file that contains the source code to generate the manuscript.  <br> • `references.bib`: a bibliography file that contains all the references used in the manuscript. <br> • `kyuriapa.bst`: a BibTeX style file used to render [APA](https://apastyle.apa.org/style-grammar-guidelines/paper-format) formatting and style of references.|
| [`simulation`](./%20%20%20%20simulation)       | It contains all `R` scripts required to run the simulation study.  <br> • `01_simulation.R`: simulate models, generate data, and run the algorithms. <br> • `02_simulation_results.R`: evaluate performance of algorithms and create figures. <br> • `03_sensitivity_analysis1.R`: run the secondary simulation with randomly sampled coefficients for the regression matrix $\mathbf{B}$. <br> • `04_sensitivity_analysis2.R`: conduct the secondary simulation by randomly sampling coefficients with *positive* values for the regression matrix $\mathbf{B}$. <br> • `05_sensitivity_analysis3.R`: run secondary analysis  with varying $\alpha$ levels. <br> • `06_extra_analysis.R`: investigate the unexpected patterns in the *5-variable dense* conditions from the main simulation study.<br> • `data`: results from all simulations are stored as `Rdata` files because of the considerable processing time required to run them.   |
| [`empirical_example`](./%20%20%20empirical_example) | It contains all R scripts required to run the empirical analysis. <br> • `01_empirical_analysis.R`: run the empirical analysis by applying the CCD, FCI and CCI algorithms to empirical data from [McNally et al. (2017)](https://www.cambridge.org/core/journals/psychological-medicine/article/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900#article). <br> • `02_stability_analysis.R`: perform the stability analysis on the empirical data (*Appendix I* in the paper). <br> • `data`: raw data that is publicly available [here](https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900#supplementary-materials) and output objects from the sensitivity analysis saved as `Rdata` files.               |
| [`utils`](./%20%20utils)               | It includes all the supporting functions required to run the other `R` scripts for conducting analyses. |
| [`supplementary_material`](./%20supplementary_material)           |It contains all files related to the supplementary material of the paper. <br> The supplementary material can be accessed via this [link](https://kyurip.quarto.pub/discovering-cyclic-causal-models/).|
| [`figures`](./figures)               | It contains all the figures that are presented in the paper.                    |
| [`others`](./others)               | It contains all types of subsidiary files  created during the course of the project (e.g., [research proposal](./others/Proposal), [presentation](./others/Presentation), [interim report](https://github.com/KyuriP/CCM_Discovery/tree/main/others/Research%20Report)).                        |

## Reproducing study
If interested in reproducing the results of:
- simulation study, please refer to [`guide_simuation`](https://github.com/KyuriP/CCM_Discovery/blob/main/%20%20%20%20simulation/README.md).
- empirical anlaysis, please refer to [`guide_empirical`](https://github.com/KyuriP/CCM_Discovery/blob/main/%20%20%20empirical_example/README.md).

## Session Information
<details>
<summary><b><i>R session & Pacakge information</i></b></summary>

           
```
─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.3 (2023-03-15)
 os       macOS Ventura 13.0
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Amsterdam
 date     2023-05-01
 rstudio  2022.07.2+576 Spotted Wakerobin (desktop)
 pandoc   2.19.2 @ /Applications/RStudio.app/Contents/MacOS/quarto/bin/tools/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────
 package      * version  date (UTC) lib source
 abind          1.4-5    2016-07-21 [1] CRAN (R 4.2.0)
 backports      1.4.1    2021-12-13 [1] CRAN (R 4.2.0)
 base64enc      0.1-3    2015-07-28 [1] CRAN (R 4.2.0)
 bdsmatrix      1.3-6    2022-06-03 [1] CRAN (R 4.2.0)
 BiocGenerics * 0.42.0   2022-04-26 [1] Bioconductor
 BiocManager  * 1.30.20  2023-02-24 [1] CRAN (R 4.2.0)
 broom          1.0.4    2023-03-11 [1] CRAN (R 4.2.0)
 cachem         1.0.7    2023-02-24 [1] CRAN (R 4.2.0)
 callr          3.7.3    2022-11-02 [1] CRAN (R 4.2.0)
 car            3.1-2    2023-03-30 [1] CRAN (R 4.2.3)
 carData        3.0-5    2022-01-06 [1] CRAN (R 4.2.0)
 CCI.KP       * 0.1.0    2023-01-21 [1] Github (KyuriP/CCI_KP@984bf12)
 checkmate      2.1.0    2022-04-21 [1] CRAN (R 4.2.0)
 cli            3.6.1    2023-03-23 [1] CRAN (R 4.2.0)
 clipr          0.8.0    2022-02-22 [1] CRAN (R 4.2.0)
 clue           0.3-64   2023-01-31 [1] CRAN (R 4.2.0)
 cluster        2.1.4    2022-08-22 [1] CRAN (R 4.2.3)
 codetools      0.2-19   2023-02-01 [1] CRAN (R 4.2.3)
 colorspace     2.1-0    2023-01-23 [1] CRAN (R 4.2.0)
 corpcor        1.6.10   2021-09-16 [1] CRAN (R 4.2.0)
 crayon         1.5.2    2022-09-29 [1] CRAN (R 4.2.0)
 curl           5.0.0    2023-01-12 [1] CRAN (R 4.2.0)
 data.table     1.14.8   2023-02-17 [1] CRAN (R 4.2.0)
 DEoptimR       1.0-11   2022-04-03 [1] CRAN (R 4.2.0)
 devtools     * 2.4.5    2022-10-11 [1] CRAN (R 4.2.0)
 digest         0.6.31   2022-12-11 [1] CRAN (R 4.2.0)
 DOT          * 0.1      2016-04-16 [1] CRAN (R 4.2.0)
 dplyr        * 1.1.2    2023-04-20 [1] CRAN (R 4.2.0)
 ellipsis       0.3.2    2021-04-29 [1] CRAN (R 4.2.0)
 evaluate       0.20     2023-01-17 [1] CRAN (R 4.2.0)
 fansi          1.0.4    2023-01-22 [1] CRAN (R 4.2.0)
 fastICA        1.2-3    2021-09-25 [1] CRAN (R 4.2.0)
 fastmap        1.1.1    2023-02-24 [1] CRAN (R 4.2.0)
 fdrtool        1.2.17   2021-11-13 [1] CRAN (R 4.2.0)
 foreign        0.8-84   2022-12-06 [1] CRAN (R 4.2.3)
 Formula        1.2-5    2023-02-24 [1] CRAN (R 4.2.0)
 fs             1.6.1    2023-02-06 [1] CRAN (R 4.2.0)
 furrr        * 0.3.1    2022-08-15 [1] CRAN (R 4.2.0)
 future       * 1.32.0   2023-03-07 [1] CRAN (R 4.2.0)
 generics       0.1.3    2022-07-05 [1] CRAN (R 4.2.0)
 ggh4x        * 0.2.4    2023-04-04 [1] CRAN (R 4.2.0)
 ggm            2.5      2020-02-16 [1] CRAN (R 4.2.0)
 ggplot2      * 3.4.2    2023-04-03 [1] CRAN (R 4.2.0)
 ggpubr       * 0.6.0    2023-02-10 [1] CRAN (R 4.2.0)
 ggsignif       0.6.4    2022-10-13 [1] CRAN (R 4.2.0)
 glasso         1.11     2019-10-01 [1] CRAN (R 4.2.0)
 globals        0.16.2   2022-11-21 [1] CRAN (R 4.2.0)
 glue           1.6.2    2022-02-24 [1] CRAN (R 4.2.0)
 graph        * 1.74.0   2022-04-26 [1] Bioconductor
 gridExtra      2.3      2017-09-09 [1] CRAN (R 4.2.0)
 gtable         0.3.3    2023-03-21 [1] CRAN (R 4.2.0)
 gtools         3.9.4    2022-11-27 [1] CRAN (R 4.2.0)
 Hmisc          5.0-1    2023-03-08 [1] CRAN (R 4.2.0)
 htmlTable      2.4.1    2022-07-07 [1] CRAN (R 4.2.0)
 htmltools      0.5.5    2023-03-23 [1] CRAN (R 4.2.3)
 htmlwidgets    1.6.2    2023-03-17 [1] CRAN (R 4.2.0)
 httpuv         1.6.9    2023-02-14 [1] CRAN (R 4.2.0)
 igraph         1.4.1    2023-02-24 [1] CRAN (R 4.2.0)
 jpeg           0.1-10   2022-11-29 [1] CRAN (R 4.2.0)
 jsonlite       1.8.4    2022-12-06 [1] CRAN (R 4.2.0)
 knitr          1.42     2023-01-25 [1] CRAN (R 4.2.0)
 later          1.3.0    2021-08-18 [1] CRAN (R 4.2.0)
 lattice        0.21-8   2023-04-05 [1] CRAN (R 4.2.0)
 lavaan         0.6-15   2023-03-14 [1] CRAN (R 4.2.0)
 lifecycle      1.0.3    2022-10-07 [1] CRAN (R 4.2.0)
 listenv        0.9.0    2022-12-16 [1] CRAN (R 4.2.0)
 magrittr     * 2.0.3    2022-03-30 [1] CRAN (R 4.2.0)
 MASS         * 7.3-58.3 2023-03-07 [1] CRAN (R 4.2.0)
 Matrix         1.5-4    2023-04-04 [1] CRAN (R 4.2.0)
 memoise        2.0.1    2021-11-26 [1] CRAN (R 4.2.0)
 mime           0.12     2021-09-28 [1] CRAN (R 4.2.0)
 miniUI         0.1.1.1  2018-05-18 [1] CRAN (R 4.2.0)
 mnormt         2.1.1    2022-09-26 [1] CRAN (R 4.2.0)
 munsell        0.5.0    2018-06-12 [1] CRAN (R 4.2.0)
 nlme           3.1-162  2023-01-31 [1] CRAN (R 4.2.3)
 nnet           7.3-18   2022-09-28 [1] CRAN (R 4.2.3)
 parallelly     1.35.0   2023-03-23 [1] CRAN (R 4.2.3)
 pbapply        1.7-0    2023-01-13 [1] CRAN (R 4.2.0)
 pbivnorm       0.6.0    2015-01-23 [1] CRAN (R 4.2.0)
 pcalg        * 2.7-8    2022-12-21 [1] CRAN (R 4.2.0)
 pillar         1.9.0    2023-03-22 [1] CRAN (R 4.2.0)
 pkgbuild       1.4.0    2022-11-27 [1] CRAN (R 4.2.0)
 pkgconfig      2.0.3    2019-09-22 [1] CRAN (R 4.2.0)
 pkgload        1.3.2    2022-11-16 [1] CRAN (R 4.2.0)
 plyr           1.8.8    2022-11-11 [1] CRAN (R 4.2.0)
 png            0.1-8    2022-11-29 [1] CRAN (R 4.2.0)
 ppcor          1.1      2015-12-03 [1] CRAN (R 4.2.0)
 prettyunits    1.1.1    2020-01-24 [1] CRAN (R 4.2.0)
 processx       3.8.0    2022-10-26 [1] CRAN (R 4.2.0)
 profvis        0.3.7    2020-11-02 [1] CRAN (R 4.2.0)
 promises       1.2.0.1  2021-02-11 [1] CRAN (R 4.2.0)
 ps             1.7.4    2023-04-02 [1] CRAN (R 4.2.0)
 psych          2.3.3    2023-03-18 [1] CRAN (R 4.2.0)
 purrr        * 1.0.1    2023-01-10 [1] CRAN (R 4.2.0)
 qgraph       * 1.9.4    2023-03-21 [1] CRAN (R 4.2.0)
 quadprog       1.5-8    2019-11-20 [1] CRAN (R 4.2.0)
 R6             2.5.1    2021-08-19 [1] CRAN (R 4.2.0)
 RBGL           1.72.0   2022-04-26 [1] Bioconductor
 rcausal      * 1.2.1    2022-09-06 [1] Github (bd2kccd/r-causal@cc74f8d)
 Rcpp           1.0.10   2023-01-22 [1] CRAN (R 4.2.0)
 remotes        2.4.2    2021-11-30 [1] CRAN (R 4.2.0)
 reshape2       1.4.4    2020-04-09 [1] CRAN (R 4.2.0)
 Rgraphviz    * 2.40.0   2022-04-26 [1] Bioconductor
 rJava        * 1.0-6    2021-12-10 [1] CRAN (R 4.2.0)
 rlang          1.1.0    2023-03-14 [1] CRAN (R 4.2.0)
 rmarkdown      2.21     2023-03-26 [1] CRAN (R 4.2.3)
 robustbase     0.95-1   2023-03-29 [1] CRAN (R 4.2.0)
 rpart          4.1.19   2022-10-21 [1] CRAN (R 4.2.3)
 rstatix        0.7.2    2023-02-01 [1] CRAN (R 4.2.0)
 rstudioapi     0.14     2022-08-22 [1] CRAN (R 4.2.0)
 scales         1.2.1    2022-08-20 [1] CRAN (R 4.2.0)
 sessioninfo    1.2.2    2021-12-06 [1] CRAN (R 4.2.0)
 sfsmisc        1.1-14   2022-11-24 [1] CRAN (R 4.2.0)
 shiny          1.7.4    2022-12-15 [1] CRAN (R 4.2.0)
 stringi        1.7.12   2023-01-11 [1] CRAN (R 4.2.0)
 stringr        1.5.0    2022-12-02 [1] CRAN (R 4.2.0)
 tibble         3.2.1    2023-03-20 [1] CRAN (R 4.2.0)
 tidyr          1.3.0    2023-01-24 [1] CRAN (R 4.2.0)
 tidyselect     1.2.0    2022-10-10 [1] CRAN (R 4.2.0)
 urlchecker     1.0.1    2021-11-30 [1] CRAN (R 4.2.0)
 usethis      * 2.1.6    2022-05-25 [1] CRAN (R 4.2.0)
 utf8           1.2.3    2023-01-31 [1] CRAN (R 4.2.0)
 V8             4.2.2    2022-11-03 [1] CRAN (R 4.2.0)
 vctrs          0.6.2    2023-04-19 [1] CRAN (R 4.2.0)
 withr          2.5.0    2022-03-03 [1] CRAN (R 4.2.0)
 xfun           0.38     2023-03-24 [1] CRAN (R 4.2.0)
 xtable         1.8-4    2019-04-21 [1] CRAN (R 4.2.0)

 [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library

──────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

</details>

## Ethics & Access
All data are simulated in this study except for the one used in the emprical analysis, which is publicly available on the [Psychological Medicine Journal](https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900\#supplementary-materials) webpage.
Ethical approval for this study has been granted by the Ethical Review Board of the Faculty of Social and
Behavioural Sciences [(FETC)](https://ferb.sites.uu.nl/) at Utrecht University, with reference numbers 22-1810 and 22-1825 (see the confirmation of approval: [fetcapproval1.pdf](https://github.com/KyuriP/Thesis_KP/blob/main/others/FETC-approval/fetcapproval1.pdf), [fetcapproval2.pdf](https://github.com/KyuriP/Thesis_KP/blob/main/others/FETC-approval/fetcapproval2.pdf)).

This repository, which serves as an archive of the study is accessible to the public on [GitHub](https://github.com/KyuriP/Thesis_KP) under the license type of `GNU General Public License v3.0`. 


## Contact
The public accessibility of the project repository is maintained by Kyuri Park.
For any inquiries or feedback about the project, please feel free to contact [Kyuri Park](https://kyurip.github.io/).
