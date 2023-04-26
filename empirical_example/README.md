## How to reproduce the results of the empirical example


|Executing order| File name          | Description of functionality |
|:---:|--------------------------------|----------------------|
|***1***| `01_empircal_analysis.R`     | It contains the code to run empirical analysis. The stpes are as follows: <br> 1. It loads the data from the `data` folder. <br> 2. It applies the algorithms to the data <br> 3. It plots a statistical network model and PAGs.   |
|***2***| `02_sensitivity_analysis.R`  | It contains the code to perform the stability analysis on empirical data using the following steps: <br> 1. It generates 1000 data subsets. <br> 2. It applies the algorithms to each subset. <br> 3. It calculates the frequency of the estimated edge-endpoints. <br>4. It generates PAGs containing edge-endpoints that are consistently predicted (more than 70% of the time).|
