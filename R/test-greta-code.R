# trial variable selection without discrete sampler
#
# non-functional version

## FUNCTIONAL MODEL:
# is just MVN prior on regression parameters with hierarchical shrinkage prior
#   - need to think through coding of SNPs (they separate additive and dominant effets but it seems odd)
#   - use Stan versions of spike-and-slab approximations

# load packages
library(greta)

# settings
n <- 100
npred <- 50
prop_zero <- 0.9

# set seed
set.seed(2018-07-24)

# simulate some data
x <- matrix(rnorm(n * npred), ncol = npred)
alpha <- rnorm(1)
beta <- rnorm(npred)
beta[sample(seq_len(npred), size = floor(prop_zero * npred), replace = FALSE)] <- 0
y <- alpha + x %*% beta + rnorm(n)

# model priors
alpha_est <- greta::normal(mean = 0.0, sd = 1.0, dim = 1)
beta_est <- greta::normal(mean = 0.0, sd = 1.0, dim = npred)
beta_inc <- greta::normal(mean = 0.0, sd = 1.0, dim = npred)

# define linear predictor
beta_indicator <- (1 + sign(beta_inc)) / 2
mu <- alpha_est + x %*% (beta_est * beta_indicator)

# set likelihood
sigma_main <- greta::uniform(min = 0.0, max = 10.0, dim = 1)
distribution(y) <- greta::normal(mean = mu, sd = sigma_main)

# define model
mod <- greta::model(mu, alpha_est, beta_est, beta_inc, beta_indicator, sigma_main)

# sample from model
samples <- greta::mcmc(mod, n_samples = 100000, warmup = 10000)

# summarise fitted model
mu_fitted <- samples[[1]][, grep('mu', colnames(samples[[1]]))]
beta_fitted <- samples[[1]][, grep('beta_est', colnames(samples[[1]]))]
beta_prob <- samples[[1]][, grep('beta_indicator', colnames(samples[[1]]))]
beta_prob <- apply(beta_prob, 2, mean)


# functional version
# settings
n <- 500
npred <- 400
nbin <- 6
nbasis <- 5
ndegree <- 2
prop_zero <- 0.95

# set seed
set.seed(2018-07-24)

# simulate some data
basis_fn <- splines::bs(seq_len(nbin), df = nbasis, degree = ndegree)
x <- matrix(rnorm(n * npred), ncol = npred)
alpha <- basis_fn %*% rnorm(nbasis)
beta <- basis_fn %*% matrix(rnorm((npred * nbasis)), nrow = nbasis)
beta[, sample(seq_len(npred), size = floor(prop_zero * npred), replace = FALSE)] <- 0
y <- sweep(x %*% t(beta), 2, alpha, '+') +
  matrix(rnorm(n * nbin), ncol = nbin)

# model priors
nbasisfit <- 5
ndegreefit <- 2
alpha_est <- greta::normal(mean = 0.0, sd = 1.0, dim = nbasisfit)
beta_est <- greta::normal(mean = 0.0, sd = 1.0, dim = c(nbasisfit, npred))
beta_inc <- greta::normal(mean = 0.0, sd = 1.0, dim = npred)
basis_fit <- splines::bs(seq_len(nbin), df = nbasisfit, degree = ndegreefit)

# define linear predictor
beta_indicator <- (1 + sign(beta_inc)) / 2
mu <- sweep(x %*% (sweep(t(basis_fit %*% beta_est),
                         1, beta_indicator, '*')),
            2, (basis_fit %*% alpha_est), '+')

# set likelihood
sigma_main <- greta::uniform(min = 0.0, max = 10.0, dim = 1)
distribution(y) <- greta::normal(mean = mu, sd = sigma_main)

# define model
mod <- greta::model(mu, alpha_est, beta_est, beta_inc, beta_indicator, sigma_main)

# sample from model
samples <- greta::mcmc(mod, n_samples = 50000, warmup = 40000)

# summarise fitted model
mu_fitted <- samples[[1]][, grep('mu', colnames(samples[[1]]))]
beta_fitted <- samples[[1]][, grep('beta_est', colnames(samples[[1]]))]
beta_prob <- samples[[1]][, grep('beta_indicator', colnames(samples[[1]]))]
beta_prob <- apply(beta_prob, 2, mean)
