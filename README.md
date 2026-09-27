# Serial Correlation Monte Carlo

R script for AEM 6850 (Empirical Methods, Cornell, Fall 2025) that uses a Monte Carlo simulation to compare naive OLS and Newey-West standard errors when the regressor and the error are serially correlated.

script_files:
- serial_correlation_sim.R : simulates y = x + e with AR(1) x and e for autocorrelations of 0, 0.25, 0.5 and 0.75 and sample sizes from 25 to 1,600, tests the true null with naive and Newey-West (4 lag) standard errors 1,000 times each, and plots the rejection rates
- Serial-Correlation-Monte-Carlo.Rproj : RStudio project, open this first so the script finds output_figure/

output_figure:
- Rejection_rates_Serial_Dependence.png : rejection rates for naive and Newey-West standard errors

Other files:
- readme.rtf : full readme with general, methodological and data-specific information
