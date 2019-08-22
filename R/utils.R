# helper functions for fGWAS analysis

# calculate differences but include first value as well
diff_with_first <- function(x) {
  if (sum(!is.na(x)) > 1)
    out <- c(x[1], diff(x))
  else
    out <- x[1]
  out
}

# rebase a vector intended as an index to start at 1 and count up with no gaps
rebase_index <- function(x)
  as.integer(as.factor(x))

# define a distance function to calculate Euclidean distance between individuals
sq_sum <- function(x) sqrt(sum(x ^ 2, na.rm = TRUE))
dist_fn <- function(target, remainder) {
  diffs <- sweep(remainder, 2, target, "-")
  apply(diffs, 1, sq_sum)
}

# internal function used by impute_from_nearest
impute_single <- function(i, data) {
  target <- data[i, ]
  remainder <- data[-i, ]
  missing <- is.na(target)
  if (any(missing)) {
    if (sum(missing) > 1)
      complete <- apply(remainder[, missing], 1, function(x) !any(is.na(x)))
    else
      complete <- !is.na(remainder[, missing])
    remainder <- remainder[complete, ]
    dists <- dist_fn(target, remainder)
    nearest_neighbours <- dists == min(dists)
    if (sum(nearest_neighbours) > 1) {
      if (sum(missing) > 1)
        data[i, missing] <- apply(remainder[nearest_neighbours, missing], 2, median)
      else
        data[i, missing] <- median(remainder[nearest_neighbours, missing])
    } else {
      data[i, missing] <- remainder[nearest_neighbours, missing]
    }
  }
  data[i, ]
}

# impute missing data from nearest neighbour
impute_from_nearest <- function(data) {
  t(sapply(seq_len(nrow(data)), impute_single, data))
}

# set up an AR1 covariance with parameters rho and sigma
ar1 <- function(t, rho, sigma) {
  
  # matrix of time modifiers
  t_mat <- outer(t, t, FUN = "-")
  t_mat <- pmax(t_mat, 0)
  
  # which elements to include (don't want upper triangular ones)
  mask <- lower.tri(t_mat, diag = TRUE)
  
  # matrix of rho ^ n contributions
  rho_mat <- (rho ^ t_mat) * mask
  
  # add overall variance and epsilons
  innovations <- normal(0, 1, dim = length(t)) * sigma
  
  # return ar1-ified values
  rho_mat %*% innovations
  
}

# setup a legendre polynomial function
legendre_list <- list(function(x) rep(1, length(x)),
                      function(x) x,
                      function(x) -0.5 + 1.5 * (x ^ 2),
                      function(x) -1.5 * x + 2.5 * (x ^ 3),
                      function(x) 0.375 - 3.75 * (x ^ 2) + 4.375 * (x ^ 4),
                      function(x) 1.875 * x - 8.75 * (x ^ 3) + 7.875 * (x ^ 5))
legendre_poly <- function(n, x) {
  t(sapply(seq_len(n + 1), function(i, x) legendre_list[[i]](x), x))
}
