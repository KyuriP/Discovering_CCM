## How to reproduce the results of the simulation study

`01_simulation.R` and `02_simulation_result.R` concern the main simulation study. <br>
The rest of files contain respective analyses that are conducted as sensitivity checks, except for `06_extra_analysis.R`, which investigates the unexpected patterns observed in the main simulation study.

|Executing order| File name          | Description of functionality |
|:---:|--------------------------------|----------------------|
|***1***| `01_simulation.R`              | Generates data from eight different pre-specified models and applies three algorithms (CCD, FCI, and CCI) to the simulated data. <br> ***Note***: Due to the long running time, the code that executes algorithms is currently commented out. If interested, simply uncomment those lines and run the code.  |
|***2***| `02_simulation_results.R`      | Assesses the performance of each algorithm for every condition, organizes the performance metrics into a structured dataframe, and generates a figure for each evaluation metric.  |
|***3***| `03_sensitivity_analysis1.R`   | Conducts a secondary simulation study by randomly sampling parameters for the weight matrix $\mathbf{B}$. The simulation is then repeated to generate data, apply the algorithms, evaluate their performance, and create figures representing the performance metrics. |
|***4***| `04_sensitivity_analysis2.R`   | Conducts a secondary simulation where the parameters for the weight matrix $\mathbf{B}$ are randomly sampled from *positive* values only. The simulation is repeated to generate data, run the algorithms, evaluate their performance, and create figures of the performance metrics. |
|***5***| `05_sensitivity_analysis3.R`   | Conducts a simulation study where the $\alpha$ level is adjusted based on the sample size ($N$). The simulation is repeated to generate data, run the algorithms with varying alpha levels depending on $N$, evaluate their performance, and generate figures of the performance metrics.| 
|***6***| `06_extra_analysis.R`          | It investigates the unexpected patterns found in the *5-variable dense conditions* by <br> *(i)* first, examining the partial correlations and <br> *(ii)* second, checking all the results of conditional independence tests. |

> ***Note.*** Each script header contains a more detailed description of the step-by-step procedure.
