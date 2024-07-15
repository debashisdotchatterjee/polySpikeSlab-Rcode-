# Load necessary libraries
if (!requireNamespace("qtl", quietly = TRUE)) {
  install.packages("qtl")
}
if (!requireNamespace("MCMCpack", quietly = TRUE)) {
  install.packages("MCMCpack")
}
if (!requireNamespace("coda", quietly = TRUE)) {
  install.packages("coda")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("missForest", quietly = TRUE)) {
  install.packages("missForest")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
if (!requireNamespace("susieR", quietly = TRUE)) {
  install.packages("susieR")
}

library(qtl)
library(MCMCpack)
library(coda)
library(ggplot2)
library(missForest)
library(tidyr)
library(susieR)

# Load the hyper dataset
data(hyper)

# Calculate genotype probabilities
hyper <- calc.genoprob(hyper, step = 1)

# Extract genotype probabilities and phenotypes
genoprobs <- pull.geno(hyper)
phenotype <- hyper$pheno[,1]

# Convert genoprobs to double precision
genoprobs_double <- apply(genoprobs, 2, as.double)

# Impute missing values if any (using missForest for demonstration purposes)
if (sum(is.na(genoprobs_double)) > 0) {
  genoprobs_double <- missForest(genoprobs_double)$ximp
}

# Bayesian fine-mapping function (modified for real data)
bayesian_finemap <- function(genoprobs, phenotype, prior_prob, slab_sd, n_iter = 10000) {
  spike_and_slab_prior <- function(beta, prior_prob, slab_sd) {
    spike <- dnorm(beta, 0, 1e-6) * prior_prob
    slab <- dnorm(beta, 0, slab_sd) * (1 - prior_prob)
    return(spike + slab)
  }
  
  n_markers <- ncol(genoprobs)
  beta_samples <- matrix(0, n_iter, n_markers)
  beta_current <- rnorm(n_markers, 0, 0.01)
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  
  for (i in 1:n_iter) {
    for (j in 1:n_markers) {
      beta_proposal <- rnorm(1, beta_current[j], 0.2)
      log_lik_current <- sum(dnorm(phenotype, genoprobs %*% beta_current, sd(phenotype), log = TRUE))
      log_lik_proposal <- sum(dnorm(phenotype, genoprobs %*% replace(beta_current, j, beta_proposal), sd(phenotype), log = TRUE))
      log_prior_current <- log(spike_and_slab_prior(beta_current[j], prior_prob, slab_sd))
      log_prior_proposal <- log(spike_and_slab_prior(beta_proposal, prior_prob, slab_sd))
      
      if (is.finite(log_lik_current) && is.finite(log_lik_proposal) && is.finite(log_prior_current) && is.finite(log_prior_proposal)) {
        log_acceptance_ratio <- (log_lik_proposal + log_prior_proposal) - (log_lik_current + log_prior_current)
        if (is.finite(log_acceptance_ratio) && log(runif(1)) < log_acceptance_ratio) {
          beta_current[j] <- beta_proposal
        }
      }
      
      if (is.na(beta_current[j])) {
        stop(paste("NA detected in beta_current at iteration", i, "marker", j))
      }
    }
    beta_samples[i, ] <- beta_current
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(beta_samples)
}

# Run Bayesian fine-mapping on real data
n_iter <- 10000
prior_prob <- 0.1
slab_sd <- 1.0
beta_samples <- bayesian_finemap(genoprobs_double, phenotype, prior_prob, slab_sd, n_iter)
beta_mcmc <- mcmc(beta_samples)

# Summarize the posterior samples
beta_means <- apply(beta_samples, 2, mean)
beta_sds <- apply(beta_samples, 2, sd)
beta_df <- data.frame(Marker = 1:length(beta_means), Mean = beta_means, SD = beta_sds)

# Run SuSiE fine-mapping on the same dataset
X <- as.matrix(genoprobs_double)
y <- as.numeric(phenotype)
susie_fit <- susie(X, y)
susie_posterior_means <- susie_fit$beta

# Ensure the length matches by trimming or padding (if necessary)
if (length(susie_posterior_means) > length(beta_means)) {
  susie_posterior_means <- susie_posterior_means[1:length(beta_means)]
} else if (length(susie_posterior_means) < length(beta_means)) {
  susie_posterior_means <- c(susie_posterior_means, rep(0, length(beta_means) - length(susie_posterior_means)))
}

# Create a data frame for plotting
comparison_df <- data.frame(
  Marker = 1:length(beta_means),
  Bayesian_Means = beta_means,
  SuSiE_Means = susie_posterior_means
)

# Plot the comparison
png(filename = "plots/bayesian_vs_susie_means.png")
ggplot(comparison_df, aes(x = Marker)) +
  geom_line(aes(y = Bayesian_Means, color = "Bayesian Fine-Mapping"), size = 1) +
  geom_line(aes(y = SuSiE_Means, color = "SuSiE"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "Method", values = c("Bayesian Fine-Mapping" = "red", "SuSiE" = "blue")) +
  labs(title = "Comparison of Bayesian Fine-Mapping and SuSiE Results",
       x = "Marker Index",
       y = "Posterior Mean Beta Coefficient") +
  theme_minimal()
dev.off()


############################

# Initial QTL analysis using Haley-Knott regression
hyper <- calc.genoprob(hyper, step = 1)
out_hk <- scanone(hyper, method = "hk")

# Extracting QTL results
qtl_results <- out_hk$lod
marker_names <- rownames(out_hk)

# Bayesian fine-mapping posterior means (already calculated)
bayesian_means <- beta_means

# Creating a data frame for comparison
model_comparison_df <- data.frame(
  Marker = 1:length(marker_names),
  QTL_Results = qtl_results,
  Bayesian_Means = bayesian_means
)

# Plotting the comparison
png(filename = "plots/model_comparison.png")
ggplot(model_comparison_df, aes(x = Marker)) +
  geom_line(aes(y = QTL_Results, color = "Initial QTL Scan"), size = 1, color = "black") +
  geom_line(aes(y = Bayesian_Means, color = "Bayesian Fine-Mapping"), size = 1, color = "red") +
  labs(title = "Model Comparison: QTL Scan vs Bayesian Fine-Mapping",
       x = "Marker Index",
       y = "LOD / Posterior Mean Beta Coefficient") +
  theme_minimal() +
  scale_color_manual(name = "Method", values = c("Initial QTL Scan" = "black", "Bayesian Fine-Mapping" = "red")) +
  theme(legend.position = "top")
dev.off()


