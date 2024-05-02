This is the Github Repo for Bayesian Learning of Clinically Meaningful Sepsis Phenotypes in Northern Tanzania. Here, you will find code to implement CLustering Around Meaningful Regions (CLAMR) and the simulation study discussed in Section 3.

* ```simulation_study_parallel.R``` computes the ARI and number of clusters for CLAMR and other clustering methods on repeated replications of synthetic data. Note that this document uses parallel computing for the different replications.
* ```write_simulation_priors.R ``` creates the CLAMR hyperparameters used in the simulation study by applying equation (6) to the known MRs.
* ```rfuncts``` contains all R functions, ```rcppfuncts``` contains all RCPP functions.
* ```rfuncts/marginal_cluster_list.R``` implements the Gibbs sampler for CLAMR.
* ```rfuncts/vanilla_gmm.R``` implements the Gibbs sampler for a Bayesian GMM.
