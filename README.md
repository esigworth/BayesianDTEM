# Bayesian dose toxo-equivalence model (BayesianDTEM)

The R scripts contained in this repository contain the code necessary to reproduce all simulations from our manuscript, "Building a dose toxo-equivalence model from a Bayesian meta-analysis of published clinical trials." In each folder, there are four main script types as well as four .Rda dataframes containing the seeds and dose samples needed to exactly recreate our findings. The four script types are as follows:

- component_functions : set of functions used to generate data, fit models, and calculate equivalence relationships based on posterior samples
- full_run : set of functions used to complete one full replication of the simulation, incorporating all functions from component_functions
- example_execution : a script demonstrating the setting of modeling conditions and the execution of one replication, formatted to work well when running multiple models in parallel
- post_sim_calculations : function that takes in a list of simulation replications (each element of which is output from the example_execution script) and performs calculations to summarize the performance, including calculating MSE, coverage percentages, bias, efficiency, and the median estimated equivalence curves

In the correctly_specified folder, there is also a file example_figure_generation.R that shows how we created the summary figures for bias, efficiency, and equivalence presented in our manuscript, written for the particular case of a correctly specified model with additional covariates. 
