## How to reproduce the results of the simulation study

>`01_simulation.R` and `02_simulation_result.R` concern the main simulation study.
The rest of files contain respective analyses that are conducted as sensitivity checks, except for `06_extra_analysis.R`, which investigates the unexpected patterns observed in the main simulation study.

|Executing order| File name          | Description of functionality |
|:---:|--------------------------------|----------------------|
|***1***| `01_simulation.R`              | It generates data from 8 different pre-specified models and runs three algorithms (CCD, FCI, and CCI) using the simulated data. <br> ***Note***: Due to the long running time, the code that executes algorithms is currently commented out. If interested, simply uncomment those lines and run the code.  |
|***2***| `02_simulation_results.R`      | It assesses the performance of each algorithm for every condition, arranges the performance metrics in a well-structured dataframe, and generates a figure of the performance for each evaluation metric.  |
|***3***| `03_sensitivity_analysis1.R`   | It conducts a secondary simulation study where parameters for B matrices are randomly sampled. The simulation is repeated to generate data, run the algorithms, evaluate their performance, and create figures of the performance metrics. |
|***4***| `04_sensitivity_analysis2.R`   | It conducts a secondary simulation study where parameters for B matrices are randomly sampled, but this time only from positive values. The simulation is repeated to generate data, run the algorithms, evaluate their performance, and create figures of the performance metrics. |
|***5***| `05_sensitivity_analysis3.R`   | It performs a secondary simulation study where alpha levels are adjusted based on the sample size. The simulation is repeated to generate data, run the algorithms (with varying alpha levels depending on N), evaluate their performance, and create figures of the performance metrics.| 
|***6***| `06_extra_analysis.R`          | It investigates the unexpected patterns found in the *5-variable dense conditions* by first examining the partial correlations and second check all the conditional independence test results. |
