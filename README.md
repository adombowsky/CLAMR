# CLAMR
This is the Github Repo for Bayesian Learning of Clinically Meaningful Sepsis Phenotypes in Northern Tanzania. Here, you will find code to implement CLustering Around Meaningful Regions (CLAMR) and the simulation study discussed in Section 3.

* ```simulation_study_parallel.R``` implements the simulation study in Section 3.1.
* ```simulation_study2_parallel.R``` implements the simulation study in Section 3.2.
* ```simulation_study_well_specified.R``` implements the simulation study in Section D of the supplementary material.
* ```write_simulation_priors.R ``` creates the CLAMR hyperparameters used in the simulation studies by applying equation (6) to the known MRs.
* ```rfuncts``` contains all R functions, ```rcppfuncts``` contains all RCPP functions.
* ```rfuncts/marginal_cluster_list.R``` implements the Gibbs sampler for CLAMR.
* ```rfuncts/vanilla_gmm.R``` implements the Gibbs sampler for a Bayesian GMM.
