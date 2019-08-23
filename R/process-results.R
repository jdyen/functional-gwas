# add post-processing of results

### load fitted
draws <- readRDS("outputs/sgwas_draws.rds")

# want samples in one matrix
samples <- do.call(rbind, draws)

### convergence checks

# plot chains
bayesplot::mcmc_trace(draws, pars = "beta_dominant[1,1]")

### posterior summaries


### check (co)variances and associated parameters


### calculate model fit


### identify included/excluded SNPs (balance swamping vs masking)
additive <- samples[, grep("beta_additive", colnames(samples))]
dominant <- samples[, grep("beta_dominant", colnames(samples))]

delta_seq <- seq(0.1, 0.7, length = 20)
ncoef_additive <- ncoef_dominant <- matrix(NA, nrow = nrow(additive), ncol = length(delta_seq))
for (i in seq_along(delta_seq)) {
  
  delta <- delta_seq[i]
  
  ncoef_additive[, i] <- apply(additive, 1, sequential_m2)
  ncoef_dominant[, i] <- apply(dominant, 1, sequential_m2)
  
}





### create plots (or outputs to be used to create plots)
