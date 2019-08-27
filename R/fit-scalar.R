# gwas analysis with scalar phenotype

# where is everything?
setwd("~/Dropbox/research/functional-gwas/")

# we want to use the greta package to fit models
library(greta)

# need some helper functions
source("R/utils.R")

# load some data setso
growth_stats <- read.csv("data/compiled/gp_fgwas_scalar_pheno.csv", stringsAsFactors = FALSE)
snp_data <- read.csv("data/compiled/gp_fgwas_snps.csv", stringsAsFactors = FALSE)
covariate_data <- read.csv("data/compiled/gp_fgwas_covariates.csv", stringsAsFactors = FALSE)
genetic_data <- read.csv("data/compiled/gp_fgwas_genetics.csv", stringsAsFactors = FALSE)

# subset SNP data while testing
snp_data <- snp_data[, seq_len(20)]

# snp_data as a matrix will make distance calcs faster
snp_data <- as.matrix(snp_data)

# impute missing SNPs using a nearest-neighbour approach
snp_data <- impute_from_nearest(snp_data)

# indices for basis expansion
n_obs <- nrow(growth_stats)
n_snp <- ncol(snp_data)

# define prior for overall intercept
alpha <- normal(0, 10)

# prepare covariate data
x_numeric <- as.matrix(covariate_data[, c(4, 6:8)])
x_numeric <- impute_from_nearest(x_numeric)
x_design <- model.matrix( ~ -1 + oto_age + latitude + longitude + total_length_mm,
                          data = as.data.frame(x_numeric))
n_covar <- ncol(x_design)

# define priors for standard covariates
beta_covariate <- normal(0, 1, dim = n_covar)

# define hyperpriors for shrinkage regression prior
lambda_additive <- gamma(0.001, 0.001)
tau_additive <- gamma(1, 1 / (lambda_additive), dim = n_snp)
lambda_dominant <- gamma(0.001, 0.001)
tau_dominant <- gamma(1, 1 / (lambda_dominant), dim = n_snp)
sigma_main <- normal(0, 1, truncation = c(0, Inf))

# define shrinkage prior on additive SNP effects
sigma_additive <- sigma_main * tau_additive
beta_additive <- normal(0, sigma_additive)

# repeat for dominant SNP effects
sigma_dominant <- sigma_main * tau_dominant
beta_dominant <- normal(0, sigma_dominant)

# define priors for random effects
individual <- rebase_index(covariate_data$id)
pop <- rebase_index(covariate_data$pop)
sigma_ind <- normal(0, 1, truncation = c(0, Inf))
sigma_pop <- normal(0, 1, truncation = c(0, Inf))
gamma_ind <- normal(0, sigma_ind, dim = max(individual))
gamma_pop <- normal(0, sigma_pop, dim = max(pop))

# create some SNP data sets that capture our key processes
snps_additive <- snp_data - 1
snps_dominant <- ifelse(snp_data == 1, 1, 0)

# define linear predictor
mu <- alpha + 
  x_design %*% beta_covariate +
  snps_additive %*% beta_additive + 
  snps_dominant %*% beta_dominant +
  gamma_ind[individual] + gamma_pop[pop]

# define final variance prior
sigma_total <- normal(0, 1, truncation = c(0, Inf))

# define likelihood
growth_vec <- growth_stats$ave_growth
distribution(growth_vec) <- normal(mu, sigma_total)

# define greta model
mod <- model(beta_additive, beta_dominant)

# give mcmc settings
chains <- 4
n_samples <- 1000
warmup <- n_samples

# set some random initial values
inits <- lapply(seq_len(chains),
                function(i) initials(beta_additive = rnorm(n_snp),
                                     beta_dominant = rnorm(n_snp)))

# sample from model
draws <- mcmc(mod,
              sampler = hmc(Lmin = 5, Lmax = 100, epsilon = 0.1, diag_sd = 1),
              n_samples = n_samples,
              warmup = warmup,
              chains = chains,
              initial_values = inits)

# save fitted
saveRDS(draws, file = paste0("outputs/sgwas_draws_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"))
