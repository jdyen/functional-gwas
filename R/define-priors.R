library(greta)

# load some data sets
growth_data <- read.csv("data/compiled/gp_fgwas_pheno.csv", stringsAsFactors = FALSE)
snp_data <- read.csv("data/compiled/gp_fgwas_snps.csv", stringsAsFactors = FALSE)
covariate_data <- read.csv("data/compiled/gp_fgwas_covariates.csv", stringsAsFactors = FALSE)
genetic_data <- read.csv("data/compiled/gp_fgwas_genetics.csv", stringsAsFactors = FALSE)

# subset SNP data while testing
# snp_data <- snp_data[, seq_len(50)]

# snp_data as a matrix will make distance calcs faster
snp_data <- as.matrix(snp_data)

# define a distance function to calculate Euclidean distance between individuals
sq_sum <- function(x) sqrt(sum(x ^ 2, na.rm = TRUE))
dist_fn <- function(target, remainder) {
  diffs <- sweep(remainder, 2, target, "-")
  apply(diffs, 1, sq_sum)
}


### PULL THIS OUT and make it a function, re-use for covariate data

# impute missing SNPs using a nearest-neighbour approach
for (i in seq_len(nrow(snp_data))) {
  snp_sub <- snp_data[i, ]
  snp_remain <- snp_data[-i, ]
  missing <- is.na(snp_sub)
  if (any(missing)) {
    if (sum(missing) > 1)
      complete <- apply(snp_remain[, missing], 1, function(x) !any(is.na(x)))
    else
      complete <- !is.na(snp_remain[, missing])
    snp_remain <- snp_remain[complete, ]
    dists <- dist_fn(snp_sub, snp_remain)
    nearest_neighbours <- dists == min(dists)
    if (sum(nearest_neighbours) > 1) {
      if (sum(missing) > 1)
        snp_data[i, missing] <- apply(snp_remain[nearest_neighbours, missing], 2, median)
      else
        snp_data[i, missing] <- median(snp_remain[nearest_neighbours, missing])
    } else {
      snp_data[i, missing] <- snp_remain[nearest_neighbours, missing]
    }
  }
}

# set up an AR1 covariance with parameter rho
abs_diff <- function(x, y) abs(x - y)
ar1 <- function(t, rho) {
  t_diff <- outer(t, t, FUN = abs_diff)
  rho ^ (t_diff)
}

# setup a legendre polynomial function
legendre_list <- list(function(x) rep(1, length(x)),
                      function(x) x,
                      function(x) -0.5 + 1.5 * (x ^ 2),
                      function(x) -1.5 * x + 2.5 * (x ^ 3),
                      function(x) 0.375 - 3.75 * (x ^ 2) + 4.375 * (x ^ 4),
                      function(x) 1.875 * x - 8.75 * (x ^ 3) + 7.875 * (x ^ 5))
legendre_poly <- function(n, x) {
  sapply(seq_len(n + 1), function(i, x) legendre_list[[i]](x), x)
}

# indices for basis expansion
n_obs <- nrow(growth_data)
n_snp <- ncol(snp_data)
n_time <- ncol(growth_data)

# set legendre degree (order = nu - 1 = degree + 1)
order <- 4
nu <- order + 1
degree <- order - 1

# define prior for overall intercept
alpha <- normal(0, 10, dim = order)

# prepare covariate data
## IMPUTE FIRST OR AVOID na.rm
x_design <- model.matrix( ~ pop + latitude + longitude + total_length_mm,
                          data = covariate_data)
n_covar <- ncol(x_design)

# define priors for standard covariates
beta_covariate <- normal(0, 10, dim = c(order, n_covar))

# define hyperpriors for shrinkage regression prior
lambda2_additive <- gamma(0.01, 0.01)  # actually lambda2 / 2, factored out 2 from numerator below
tau2_additive <- gamma((nu + 1) / 2, 1 / (nu * lambda2_additive), dim = n_snp)
lambda2_dominant <- gamma(0.01, 0.01)
tau2_dominant <- gamma((nu + 1) / 2, 1 / (nu * lambda2_dominant), dim = n_snp)

# define shrinkage prior on additive SNP effects
sigma_main <- normal(0, 10, truncation = c(0, Inf))
sigma_additive <- sigma_main * rep(sqrt(tau2_additive), each = order)
beta_additive <- normal(0, sigma_additive)
dim(beta_additive) <- c(order, n_snp)

# repeat for dominant SNP effects
sigma_dominant <- sigma_main * rep(sqrt(tau2_dominant), each = order)
beta_dominant <- normal(0, sigma_dominant)
dim(beta_dominant) <- c(order, n_snp)

# define priors for random effects
gamma_year
gamma_pop

# pull out the number of timepoints for which we need covariance estimates
num_obs <- apply(growth_data, 1, function(x) max(which(!is.na(x))))
t_vec <- seq_len(n_time)
rho <- uniform(0, 1)

# create an AR1 covariance matrix for the max number of observations
covar_full <- ar1(t_vec, rho)

# pull out matrix subsets for each number of observations
#  - could flatten or run this in-line
covar_obs <- list()
unique_num_obs <- sort(unique(num_obs))
for (i in seq_along(unique_num_obs)) {
  idx <- seq_len(unique_num_obs[i])
  covar_obs[[i]] <- covar_full[idx, idx]
}
eps_mvn <- zeros(n_obs, n_time)
eps_mvn[num_obs == 1, ] <- normal(0, covar_obs[[1]])
for (i in seq_along(unique_num_obs)[-1]) {
  n_obs <- unique_num_obs[i]
  idx <- num_obs == n_obs
  eps_mvn[idx, seq_len(n_obs)] <- multivariate_normal(zeros(1, n_obs), covar_obs[[i]])
}

# define legendre polynomials
legendre_basis <- legendre_poly(n = degree, x = seq_len(n_time))

# test out linear predictors
mu <- t(legendre_basis %*% alpha) +
  x_design %*% t(legendre_basis %*% beta_covariate) +
  snp_data %*% t(legendre_basis %*% beta_additive) + 
  snp_data %*% t(legendre_basis %*% beta_dominant) +
  gamma_year + gamma_pop + eps_mvn


# the following gives a list of covariances -- can we flatten it?
# covar_obs[num_obs]


