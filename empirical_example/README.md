## How to reproduce the results of the empirical example


|Executing order| File name          | Description of functionality |
|:---:|--------------------------------|----------------------|
|***1***| `01_empircal_analysis.R`     | Contains the code to perform the empirical analysis (*Section 5* of the paper). <br> The stpes are as follows: <br> 1. Load the data from the `data` folder. <br> 2. Apply the algorithms to the data. <br> 3. Plot a statistical network model (Gaussian graphical model) and PAGs.   |
|***2***| `02_stability_analysis.R`  | Contains the code to perform the stability analysis on empirical data (*Appendix I* of the paper)  using the following steps: <br> 1. Generate 1000 data subsets. <br> 2. Apply the algorithms to each subset. <br> 3. Calculate the frequency of the estimated edge-endpoints. <br>4. Generate PAGs containing edge-endpoints that are consistently predicted (more than 70% of the time).|

> ***Note.*** Each script header contains a more detailed description of the step-by-step procedure.
