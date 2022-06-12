# ElbeMetalAnalysis
The codes for analyzing the trace metal pollution in the Elbe River

Bayesian network
It is a Bootstrapped tabu-learning Bayesian network for analyzing the relationship between multi-influential factors and the metal concentrations in rivers.
The program includes Bayesian network structure learning, arc confidence & strength analysis, model learning and parameter analysis, uncertainty analysis and outcome output
The files for the Bayesian network
- R code: BN_BOOT_ELBE3.R 
- Input data: BN_Elbe.xlsx

Bayesian multivariate receptor model
The program is used for the second step of the model development. The first step is done by the Positive Matrix Factorization model (designed by EPA). 
This program includes the optimization of identifiability condition in the prior factor profile estimated by the PMF model, Bayesian inference for solving chemical mass balance, and outcome output
The files for the Bayesian multivariate receptor model (This is an example for source apportionment of metals in the Elbe downstream)
- Matlab code: BMRM_D.m
- Input data: PMF_D.xls
