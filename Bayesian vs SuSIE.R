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

# Initial QTL analysis using Haley-Knott regression
out_hk <- scanone(hyper, method = "hk")

# Extract LOD scores and ensure they match the marker names
lod_scores <- out_hk$lod[, 1]
marker_names_lod <- rownames(out_hk$lod)

# Ensure genoprobs_double matches the length of lod_scores
aligned_genoprobs <- genoprobs_double[, colnames(genoprobs_double) %in% marker_names_lod]

# Reorder aligned_genoprobs columns to match marker_names_lod
aligned_genoprobs <- aligned_genoprobs[, match(marker_names_lod, colnames(aligned_genoprobs))]

# Fit models using the estimated beta coefficients
qtl_fitted <- aligned_genoprobs %*% lod_scores
bayes_fitted <- aligned_genoprobs %*% beta_means

# Check lengths of qtl_fitted and bayes_fitted
cat("Length of qtl_fitted: ", length(qtl_fitted), "\n")
cat("Length of bayes_fitted: ", length(bayes_fitted), "\n")

# Create a data frame for model comparison
comparison_df <- data.frame(
  Observed = phenotype[1:min(length(qtl_fitted), length(bayes_fitted))],  # Ensure matching lengths
  QTL_Fitted = qtl_fitted[1:min(length(qtl_fitted), length(bayes_fitted))],
  Bayes_Fitted = bayes_fitted[1:min(length(qtl_fitted), length(bayes_fitted))]
)

# Model Comparison Table
model_comparison_table <- data.frame(
  Model = c("QTL Scan", "Bayesian Fine-Mapping"),
  RMSE = c(sqrt(mean((comparison_df$Observed - comparison_df$QTL_Fitted)^2)),
           sqrt(mean((comparison_df$Observed - comparison_df$Bayes_Fitted)^2))),
  R2 = c(cor(comparison_df$Observed, comparison_df$QTL_Fitted)^2,
         cor(comparison_df$Observed, comparison_df$Bayes_Fitted)^2)
)

print(model_comparison_table)

# Save the table
write.csv(model_comparison_table, file = "plots/model_comparison_table.csv", row.names = FALSE)

# Plot the comparison
png(filename = "plots/model_comparison.png")
ggplot(comparison_df, aes(x = Observed)) +
  geom_point(aes(y = QTL_Fitted), color = "black") +
  geom_point(aes(y = Bayes_Fitted), color = "red") +
  labs(title = "Model Comparison: Observed vs Fitted Values",
       x = "Observed Phenotype",
       y = "Fitted Values") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(name = "Method", values = c("Initial QTL Scan" = "black", "Bayesian Fine-Mapping" = "red"))
dev.off()
