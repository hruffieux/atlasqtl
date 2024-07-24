library(atlasqtl)
library(ggplot2)
library(dplyr)
seed <- 123; set.seed(seed)

###################
## Simulate data ##
###################

# Example with small problem sizes:
#
n <- 100; p <- 500; p_act <- 10; q <- 100; q_act <- 100

# Candidate predictors (subject to selection)
#
# Here example with common genetic variants under Hardy-Weinberg equilibrium
#
X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.25), nrow = n)
X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.25), nrow = n)

# shuffle indices 
shuff_x_ind <- sample(p)
shuff_y_ind <- sample(q)

X <- cbind(X_act, X_inact)[, shuff_x_ind]

# Association pattern and effect sizes
#
pat <- matrix(FALSE, ncol = q, nrow = p)
bool_x <- shuff_x_ind <= p_act
bool_y <- shuff_y_ind <= q_act

pat_act <- beta_act <- matrix(0, nrow = p_act, ncol = q_act)
pat_act[sample(p_act * q_act, floor(p_act * q_act / 5))] <- 1
beta_act[as.logical(pat_act)] <-  rnorm(sum(pat_act))

pat[bool_x, bool_y] <- pat_act

# Gaussian responses
#
Y_act <- matrix(rnorm(n * q_act, mean = X_act %*% beta_act), nrow = n)
Y_inact <- matrix(rnorm(n * (q - q_act)), nrow = n)

Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]

########################
## Infer associations ##
########################

# Expectation and variance for the prior number of predictors associated with
# each response
#
p0 <- c(mean(colSums(pat)), 10)
mu_t <- 1; v_t <- 4 

devtools::load_all()


system.time(res_atlas <- atlasqtl(as.matrix(Y), as.matrix(X),
                                  p0 = c(mu_t, v_t),
                                  user_seed = 1, maxit= 50000,
                                  thinned_elbo_eval = F,
                                  batch = "y",
                                  tol = 0.1,
                                  # anneal_tol = c(1, 0.1),
                                  anneal_tol = NULL,
                                  # anneal = c(1, 2, 10),
                                  anneal = NULL,
                                  burn_in = 20,
                                  epsilon_lb = c(2, 1.5, 0.25),
                                  epsilon_it = c(0.1, 50, 0), #k, x_0, m
                                  partial_elbo = F,
                                  eval_perform = T))

########################
## Algorithm evaluation ##
########################

perform_df = res_atlas$perform_df

# Check how the ELBO changes over iterations 
perform_df %>% 
  ggplot(aes(x = iter, y = ELBO, color = annealing)) + geom_point()

# Check how the difference of ELBO changes over iterations 
perform_df %>% 
  ggplot(aes(x = iter, y = ELBO_diff, color = annealing)) + geom_point()

# Check how e changes over iterations
perform_df %>% 
  ggplot(aes(x = iter, y = e, color = annealing)) + geom_point()

perform_df %>% 
  ggplot(aes(x = iter, y = e_lb, color = annealing)) + geom_point()

perform_df %>% 
  ggplot(aes(x = iter, y = e_it, color = annealing)) + geom_point()


# lognormal_cdf <- function(x, mu, sigma, m) {
#   return(m + (1 - m) * pnorm((log(x) - mu) / sigma))
# }
# 
# perform_df %>% mutate(
#    new_e = lognormal_cdf(ELBO_diff, mu =1, sigma = 1, m = 0.1)
# ) %>% ggplot(aes(x = iter, y =  new_e)) + geom_point()
# 
# lognormal_cdf(perform_df$ELBO_diff,mu=2, sigma = 1, m = 0.25)
# 
# perform_df %>% mutate(
#   new_e = lognormal_cdf(iter, mu = 1, sigma = 0.5)
# ) %>% ggplot(aes(x = iter, y =  new_e, color = subsample)) + geom_point()


#evaluate AUROC and AUPRC
roc <- PRROC::roc.curve(scores.class0 = as.vector(res_atlas$gam_vb), weights.class0 = as.vector(as.matrix(pat)), curve =TRUE)
AUROC = roc$auc
AUROC 

pr <- PRROC::pr.curve(scores.class0 = as.vector(res_atlas$gam_vb), weights.class0 = as.vector(as.matrix(pat)), curve =TRUE)
AUPRC = pr$auc.integral
AUPRC


# system.time(res_atlas <- atlasqtl(as.matrix(Y), as.matrix(X),
#                                   p0 = c(mu_t, v_t),
#                                   user_seed = 1, maxit= 50000,
#                                   tol=0.1))


