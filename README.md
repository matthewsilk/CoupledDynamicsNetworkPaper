# CoupledDynamicsNetworkPaper

R code for the paper "Improving pandemic mitigation policies across communities through coupled dynamics of risk perception and infection".

R code to conduct the simulations and analyse the results of the paper:

- `FunctionsForHealthPaper.R` is the functions used for the network modelling
- `ScriptForHealthPaper2.R` runs the simulations
- `NetworksPaper1AnalysisCode_tidy.R` analyses the results to produce the figures and tables in the manuscript and supplementary


The folder `params/` contains parameters sts and networks configuration:
- Four `csv`s with the  parameter sets used to run the simulations and in the analysis.
- Nine `.RDS` in the subfolder `params/networks2/`, used to create the multiplex networks used in the simulations


You will first need to run the `ScriptForHealthPaper2.R`, which will population the folder `output/`, and once all simulations are done you will be able to run the script `NetworksPaper1AnalysisCode_tidy.R` 

```bash
Rscript ScriptForHealthPaper2.R 
Rscript NetworksPaper1AnalysisCode_tidy.R
```
