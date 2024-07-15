library(MCMCpack)
library(ggplot2)
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


############################


# Define a function for Bayesian fine-mapping using a spike-and-slab prior
bayesian_finemap <- function(genoprobs, phenotype, prior_prob, slab_sd) {
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
  beta_current <- rnorm(n_markers, 0, 0.01)  # Initialize with small random values
  
  # Progress indicator setup
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  
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
      
      # Check for finite values and prevent NAs
      if (is.finite(log_lik_current) && is.finite(log_lik_proposal) &&
          is.finite(log_prior_current) && is.finite(log_prior_proposal)) {
        log_acceptance_ratio <- (log_lik_proposal + log_prior_proposal) - (log_lik_current + log_prior_current)
        
        if (is.finite(log_acceptance_ratio) && log(runif(1)) < log_acceptance_ratio) {
          beta_current[j] <- beta_proposal
        }
      }
      
      # Debugging: Check for NAs in beta_current
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

# Plot the posterior means with enhanced visualization
beta_df <- data.frame(Marker = 1:length(beta_means), Mean = beta_means, SD = beta_sds)

ggplot(beta_df, aes(x = Marker, y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "red") +
  labs(title = "Posterior Means of Beta Coefficients",
       x = "Marker Index", y = "Posterior Mean") +
  theme_minimal()
##############################

# Create a directory to save the plots
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Raw data plot: Phenotype distribution
png(filename = "plots/phenotype_distribution.png")
ggplot(data.frame(phenotype), aes(x = phenotype)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(title = "Phenotype Distribution", x = "Phenotype", y = "Frequency") +
  theme_minimal()
dev.off()

# Posterior density plots using ggplot2
posterior_df <- as.data.frame(beta_samples)
posterior_df_long <- tidyr::gather(posterior_df, key = "Marker", value = "Beta")

png(filename = "plots/posterior_density_plots.png")
ggplot(posterior_df_long, aes(x = Beta, fill = Marker)) +
  geom_density(alpha = 0.5) +
  labs(title = "Posterior Density Plots", x = "Beta", y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# Trace plots for MCMC diagnostics using ggplot2
traceplot_df <- data.frame(Iteration = 1:n_iter, posterior_df)

traceplot_df_long <- tidyr::gather(traceplot_df, key = "Marker", value = "Beta", -Iteration)

png(filename = "plots/trace_plots.png")
ggplot(traceplot_df_long, aes(x = Iteration, y = Beta, color = Marker)) +
  geom_line(alpha = 0.5) +
  labs(title = "Trace Plots", x = "Iteration", y = "Beta") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# Model comparison plot: Comparing the posterior means with initial scanone results
png(filename = "plots/model_comparison.png")
plot(scanone_result, main = "Initial QTL Scan with Bayesian Fine-Mapping Overlay")
lines(beta_means, col = "red", lwd = 2)
legend("topright", legend = c("Scanone", "Bayesian Fine-Mapping"), col = c("black", "red"), lwd = 2)
dev.off()

# Posterior means with error bars
png(filename = "plots/posterior_means_with_error_bars.png")
ggplot(beta_df, aes(x = Marker, y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "red") +
  labs(title = "Posterior Means of Beta Coefficients with Error Bars",
       x = "Marker Index", y = "Posterior Mean") +
  theme_minimal()
dev.off()

# Pairwise scatter plots for beta coefficients
png(filename = "plots/pairwise_scatter_plots.png")
pairs(as.data.frame(beta_samples[, 1:10]), main = "Pairwise Scatter Plots of Beta Coefficients")
dev.off()

# Summary statistics of posterior samples
summary_beta_samples <- summary(beta_mcmc)
write.csv(summary_beta_samples$statistics, file = "plots/posterior_summary_statistics.csv")

# Plotting MCMC diagnostics using the coda package
png(filename = "plots/mcmc_diagnostics.png")
par(mfrow = c(2, 2))
autocorr.plot(beta_mcmc)
geweke.plot(beta_mcmc)
heidel.diag(beta_mcmc)
par(mfrow = c(1, 1))
dev.off()

# Autocorrelation plots for MCMC diagnostics
png(filename = "plots/autocorrelation_plots.png")
autocorr.plot(beta_mcmc, lag.max = 50, auto.layout = TRUE)
dev.off()

# Cross-correlation between the posterior samples of different markers
png(filename = "plots/cross_correlation_plots.png")
crosscorr.plot(beta_mcmc, lag.max = 50)
dev.off()

# Effective sample size for MCMC diagnostics
ess <- effectiveSize(beta_mcmc)
write.csv(ess, file = "plots/effective_sample_size.csv")

# Gelman-Rubin diagnostic for convergence
gelman_diag <- gelman.diag(beta_mcmc)
write.csv(gelman_diag$psrf, file = "plots/gelman_rubin_diagnostic.csv")


###################

# Install and load the missForest package if not already installed
if (!requireNamespace("missForest", quietly = TRUE)) {
  install.packages("missForest")
}
library(missForest)

# Impute missing values using missForest
set.seed(123)  # Set seed for reproducibility
imputed_data <- missForest(genoprobs_double)

# Extract the imputed matrix
genoprobs_imputed <- imputed_data$ximp

# Check for missing values after imputation
cat("Number of missing values in genoprobs_imputed after missForest imputation: ", sum(is.na(genoprobs_imputed)), "\n")

# Prepare the data for SuSiE
X <- as.matrix(genoprobs_imputed)
y <- as.numeric(phenotype)

# Run SuSiE
susie_fit <- susie(X, y)

# Extract posterior means from SuSiE results
susie_posterior_means <- susie_fit$mu

# Plot SuSiE results for comparison
png(filename = "plots/susie_results.png")
plot(susie_posterior_means, type = "h", main = "SuSiE Posterior Means", 
     xlab = "Marker Index", ylab = "Posterior Mean", col = "blue", lwd = 2)
dev.off()

# Benchmarking against Bayesian fine-mapping results
# Combine the results for comparison
comparison_df <- data.frame(
  Marker = 1:length(beta_means),
  Bayesian_Means = beta_means,
  SuSiE_Means = susie_posterior_means
)

# Plot the comparison
png(filename = "plots/comparison_bayesian_susie.png")
ggplot(comparison_df, aes(x = Marker)) +
  geom_line(aes(y = Bayesian_Means, color = "Bayesian"), size = 1) +
  geom_line(aes(y = SuSiE_Means, color = "SuSiE"), size = 1, linetype = "dashed") +
  labs(title = "Comparison of Bayesian Fine-Mapping and SuSiE Results",
       x = "Marker Index", y = "Posterior Mean") +
  scale_color_manual(values = c("Bayesian" = "red", "SuSiE" = "blue")) +
  theme_minimal()
dev.off()

# Save comparison results to a CSV file
write.csv(comparison_df, file = "plots/comparison_results.csv")

# Additional metrics for comparison
# Calculate correlation between the methods
correlation <- cor(comparison_df$Bayesian_Means, comparison_df$SuSiE_Means)
cat("Correlation between Bayesian Fine-Mapping and SuSiE results: ", correlation, "\n")

# Save correlation result
write.table(correlation, file = "plots/correlation.txt", row.names = FALSE, col.names = "Correlation")

# Generate summary statistics
summary_stats <- data.frame(
  Method = c("Bayesian", "SuSiE"),
  Mean = c(mean(beta_means), mean(susie_posterior_means)),
  SD = c(sd(beta_means), sd(susie_posterior_means)),
  Min = c(min(beta_means), min(susie_posterior_means)),
  Max = c(max(beta_means), max(susie_posterior_means))
)

# Save summary statistics to a CSV file
write.csv(summary_stats, file = "plots/summary_statistics.csv")
