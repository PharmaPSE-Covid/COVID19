# Codes and data used for the paper

Codes and data used for creating figures are presented separately. In each folder, multiple sub-folders are uploaded to explain the steps of calculating the results. <br><br>
The meaning of each step is as followed:<br>

* Step 1: calculating *&beta;* from data of COVID-19
* Step 2: developing a random forest model that predicts *&beta;* from mobility, environmental, and PCR data
* Step 3: Recalculating the number of infectious from predicted *&beta;*
* Step 4: (only for sensitivity analysis) calculating the total and the peak number of infections for each scenario
* Step 5: (only for sensitivity analysis) performing global sensitivity analysis

Depending on the types of figures, some steps can be skipped. For example, scenario analysis (Fig. 3) doesn't predict *&beta;* from mobility etc but defines the values of *&beta;* for each scenario.
Due to the file size, some intemidiate data are missing (e.g., input of step 4). However, you can calculate it by yourself and compare your results with our final resutls (e.g., calculating the output of step 4 by obtaining the input of step 4 from step 3).
