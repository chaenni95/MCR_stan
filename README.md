README: The MCR_stan repository contains the materials for fully Bayesian mixture cure rate model with diverse survival distributions and link functions. 


Code and Data for Simulation and Real Data Analysis

This folder contains R code and data used to replicate the results presented in the manuscript, including simulation studies and real data application on recidivism.

----------------------------------------
R Scripts Overview (Files 1–7):
----------------------------------------
1. generating_data.R  
   - Contains functions used to generate synthetic datasets for the simulation study.

2. summary_result.R  
   - Summarizes and extracts key metrics from the simulation outputs.

3_1.Weibull_model_assessment.R  
3_2.Lognormal_model_assessment.R  
3_3.Loglogistic_model_assessment.R  
   - Run model assessment simulations for each of the three AFT distributions.

4.linegraph_CI.R  
   - Function to draw line plots of the 95% coverage results across the 1000 replicates.

5.single_cox_snell.R  
   - Function to draw the Cox–Snell residual Kaplan–Meier plots for a single simulation replicate.

6.coverage.R  
   - Computes coverage probabilities for all model parameters.

7.augmented_cox_snell.R  
   - Function to generate the Cox–Snell residual plot using the full set of 1000 replicates.

----------------------------------------
Folders:
----------------------------------------
- sim500/  
  - Contains simulation outputs for sample size n = 500.

- sim1000/  
  - Contains simulation outputs for sample size n = 1000.

- data_rcode/  
  - Contains the real dataset used in the application section and the associated R scripts to run the twelve model variations.  
  - Includes Cox–Snell residual KM plots for each fitted model.

- stan_final/  
  - Contains Stan model files used to fit the Bayesian mixture cure rate models.

----------------------------------------
Supplementary Material:
----------------------------------------
- Supplementary_JDS.pdf  
  - Supplementary figures and tables referenced in the manuscript.

