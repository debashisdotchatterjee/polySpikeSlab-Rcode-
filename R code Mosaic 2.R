# Install qtl package if not already installed
if (!requireNamespace("qtl", quietly = TRUE)) {
  install.packages("qtl")
}

# Load the qtl package
library(qtl)

# Load the hyper dataset
data(hyper)

# Summarize the dataset to understand its structure
summary(hyper)

# Calculate genotype probabilities
hyper <- calc.genoprob(hyper, step = 1)

# Perform initial QTL analysis using Haley-Knott regression
scanone_result <- scanone(hyper, method = "hk")

# Plot the scanone result
plot(scanone_result, main = "Initial QTL Analysis")

# Define a function for Bayesian fine-mapping using a spike-and-slab prior
bayesian_finemap <- function(genoprobs, phenotype, prior_prob, slab_sd) {
  library(MCMCpack)
  
  # Define the prior for the beta coefficients (spike-and-slab)
  spike_and_slab_prior <- function(beta, prior_prob, slab_sd) {
    spike <- dnorm(beta, 0, 1e-6) * prior_prob
    slab <- dnorm(beta, 0, slab_sd) * (1 - prior_prob)
    return(spike + slab)
  }
  
  # Set up the MCMC
  n_iter <- 10000
  n_markers <- ncol(genoprobs)
  beta_samples <- matrix(0, n_iter, n_markers)
  beta_current <- rep(0, n_markers)
  
  for (i in 1:n_iter) {
    for (j in 1:n_markers) {
      # Sample beta from the posterior distribution
      beta_proposal <- rnorm(1, beta_current[j], 0.1)
      
      # Calculate the log-likelihoods
      log_lik_current <- sum(dnorm(phenotype, genoprobs %*% beta_current, sd(phenotype), log = TRUE))
      log_lik_proposal <- sum(dnorm(phenotype, genoprobs %*% replace(beta_current, j, beta_proposal), sd(phenotype), log = TRUE))
      
      # Calculate the log-priors
      log_prior_current <- log(spike_and_slab_prior(beta_current[j], prior_prob, slab_sd))
      log_prior_proposal <- log(spike_and_slab_prior(beta_proposal, prior_prob, slab_sd))
      
      # Check for finite values
      if (is.finite(log_lik_current) && is.finite(log_lik_proposal) &&
          is.finite(log_prior_current) && is.finite(log_prior_proposal)) {
        log_acceptance_ratio <- (log_lik_proposal + log_prior_proposal) - (log_lik_current + log_prior_current)
        
        if (log(runif(1)) < log_acceptance_ratio) {
          beta_current[j] <- beta_proposal
        }
      }
    }
    beta_samples[i, ] <- beta_current
  }
  
  return(beta_samples)
}

# Extract genotype probabilities and phenotype data
genoprobs <- pull.geno(hyper)
phenotype <- hyper$pheno[, 1]

# Run the Bayesian fine-mapping
prior_prob <- 0.01
slab_sd <- 0.1
beta_samples <- bayesian_finemap(genoprobs, phenotype, prior_prob, slab_sd)

# Summarize the posterior samples
beta_means <- apply(beta_samples, 2, mean)
beta_sds <- apply(beta_samples, 2, sd)

# Plot the posterior means
plot(beta_means, type = "h", main = "Posterior Means of Beta Coefficients", xlab = "Marker Index", ylab = "Posterior Mean")
