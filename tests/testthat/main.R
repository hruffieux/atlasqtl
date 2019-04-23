rm(list = ls())

set.seed(123)

############################
## simulate basic dataset ##
############################

n <- 100; p <- 75; q <- 20; p0 <- 10

# candidate predictors (subject to selection)
X_act <- matrix(rbinom(n * p0, size = 2, p = 0.2), nrow = n)
X_inact <- matrix(rbinom(n * (p - p0), size = 2, p = 0.2), nrow = n)
X <- cbind(X_act, X_inact)[, sample(p)]

beta <-  matrix(rnorm(p0 * q), nrow = p0)

# Gaussian outcomes
Y <- matrix(rnorm(n * q, mean = X_act %*% beta, sd = 1), nrow = n)

# remove constant variables (needed for checking dimension consistency)
X <- scale(X)
rm_cst <- function(mat_sc) mat_sc[, !is.nan(colSums(mat_sc))]
rm_coll <- function(mat_sc) mat_sc[, !duplicated(mat_sc, MARGIN = 2)]

X <- rm_cst(X)
X <- rm_coll(X)

p <- ncol(X)


########################
## atlasqtl inference ##
########################

p0_av <- c(5, 25)

# Continuous outcomes, no covariates
#
vb <- atlasqtl(Y = Y, X = X, p0_av = p0_av)

