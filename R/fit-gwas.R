## SPLIT THIS INTO utils.R and fit-gwas.R
## ADD a univariate version, coupled with the growth metrics script.

# where is everything?
setwd("~/Dropbox/research/functional-gwas/")

# we want to use the greta package to fit models
library(greta)

# need some helper functions
source("R/utils.R")

# load some data setso
growth_data <- read.csv("data/compiled/gp_fgwas_pheno.csv", stringsAsFactors = FALSE)
growth_stats <- read.csv("data/compiled/gp_fgwas_scalar_pheno.csv", stringsAsFactors = FALSE)
snp_data <- read.csv("data/compiled/gp_fgwas_snps.csv", stringsAsFactors = FALSE)
covariate_data <- read.csv("data/compiled/gp_fgwas_covariates.csv", stringsAsFactors = FALSE)
genetic_data <- read.csv("data/compiled/gp_fgwas_genetics.csv", stringsAsFactors = FALSE)

# subset SNP data while testing
snp_data <- snp_data[, seq_len(100)]

# snp_data as a matrix will make distance calcs faster
snp_data <- as.matrix(snp_data)

# impute missing SNPs using a nearest-neighbour approach
snp_data <- impute_from_nearest(snp_data)

# indices for basis expansion
n_obs <- nrow(growth_data)
n_snp <- ncol(snp_data)
n_time <- ncol(growth_data)

# set legendre degree (order = nu - 1 = degree + 1)
order <- 4
degree <- order - 1

# define prior for overall intercept
alpha <- normal(0, 10, dim = order)

# prepare covariate data
x_numeric <- as.matrix(covariate_data[, c(4, 6:8)])
x_numeric <- impute_from_nearest(x_numeric)
x_design <- model.matrix( ~ -1 + oto_age + latitude + longitude + total_length_mm,
                          data = as.data.frame(x_numeric))
n_covar <- ncol(x_design)

# define priors for standard covariates
beta_covariate <- normal(0, 1, dim = c(order, n_covar))

# define hyperpriors for shrinkage regression prior
lambda_additive <- gamma(0.001, 0.001)
tau_additive <- gamma(1, 1 / (lambda_additive), dim = n_snp)
lambda_dominant <- gamma(0.001, 0.001)
tau_dominant <- gamma(1, 1 / (lambda_dominant), dim = n_snp)
sigma_main <- normal(0, 1, truncation = c(0, Inf))

# define shrinkage prior on additive SNP effects
sigma_additive <- sigma_main * rep(tau_additive, times = order)
dim(sigma_additive) <- c(n_snp, order)
beta_additive <- normal(0, sigma_additive)

# repeat for dominant SNP effects
sigma_dominant <- sigma_main * rep(tau_dominant, times = order)
dim(sigma_dominant) <- c(n_snp, order)
beta_dominant <- normal(0, sigma_dominant)

# define priors for random effects
individual <- rebase_index(covariate_data$id)
pop <- rebase_index(covariate_data$pop)
sigma_ind <- normal(0, 1, truncation = c(0, Inf))
sigma_pop <- normal(0, 1, truncation = c(0, Inf))
gamma_ind <- normal(0, sigma_ind, dim = c(max(individual), n_time))
gamma_pop <- normal(0, sigma_pop, dim = c(max(pop), n_time))

# pull out the number of timepoints for which we need covariance estimates
num_obs <- apply(growth_data, 1, function(x) max(which(!is.na(x))))
t_vec <- seq_len(n_time)

# create a multinormal prior for the residual error with AR1 covariance across times
rho <- uniform(0, 1)
sigma_eps <- normal(0, 1, truncation = c(0, Inf))
eps_ar1 <- ar1(t_vec, rho, sigma_eps)

# define legendre polynomials
legendre_basis <- legendre_poly(n = degree, x = seq_len(n_time))
legendre_basis_prime <- t(legendre_basis)

# create some SNP data sets that capture our key processes
snps_additive <- snp_data - 1
snps_dominant <- ifelse(snp_data == 1, 1, 0)

# define linear predictor
mu <- x_design %*% beta_covariate %*% legendre_basis +
  snps_additive %*% beta_additive %*% legendre_basis + 
  snps_dominant %*% beta_dominant %*% legendre_basis +
  gamma_ind[individual, ] + gamma_pop[pop, ]
mu <- sweep(mu, 2,  legendre_basis_prime %*% alpha, "+")
mu <- sweep(mu, 2, eps_ar1, "+")

# define final variance prior
sigma_total <- normal(0, 1, truncation = c(0, Inf))

# define likelihood
growth_vec <- unlist(growth_data)
mu_vec <- c(mu)[!is.na(growth_vec)]
growth_vec <- growth_vec[!is.na(growth_vec)]
distribution(growth_vec) <- normal(mu_vec, sigma_total)

# define greta model
mod <- model(beta_additive, beta_dominant)

# give mcmc settings
chains <- 4
n_samples <- 1000
warmup <- n_samples

# sample from the model
inits <- lapply(seq_len(chains),
                function(i) initials(beta_additive = matrix(rnorm(order * n_snp),
                                                            nrow = n_snp, ncol = order),
                                     beta_dominant = matrix(rnorm(order * n_snp),
                                                            nrow = n_snp, ncol = order)))
draws <- mcmc(mod,
              n_samples = n_samples,
              warmup = warmup,
              chains = chains,
              initial_values = inits)
