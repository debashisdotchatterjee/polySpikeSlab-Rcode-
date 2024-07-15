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

library(qtl)
library(MCMCpack)
library(coda)
library(ggplot2)
library(missForest)
library(tidyr)

# Simulate data
set.seed(123)  # For reproducibility

# Parameters for simulation
n_individuals <- 100  # Number of individuals
n_markers <- 50  # Number of markers
true_beta <- rnorm(n_markers, 0, 1)  # True beta coefficients
true_beta[1:5] <- c(2, -2, 1.5, -1.5, 1)  # Strong effects for first 5 markers

# Simulate genotype matrix (X) with minor allele frequencies around 0.3
genoprobs <- matrix(rbinom(n_individuals * n_markers, 2, 0.3), nrow = n_individuals, ncol = n_markers)

# Simulate phenotype (y) with some added noise
phenotype <- genoprobs %*% true_beta + rnorm(n_individuals, 0, 1)

# Bayesian fine-mapping function (modified for simulated data)
bayesian_finemap <- function(genoprobs, phenotype, prior_prob, slab_sd, n_iter = 10000) {
  # Define the prior for the beta coefficients (spike-and-slab)
  spike_and_slab_prior <- function(beta, prior_prob, slab_sd) {
    spike <- dnorm(beta, 0, 1e-6) * prior_prob
    slab <- dnorm(beta, 0, slab_sd) * (1 - prior_prob)
    return(spike + slab)
  }
  
  # Set up the MCMC
  n_markers <- ncol(genoprobs)
  beta_samples <- matrix(0, n_iter, n_markers)
  beta_current <- rnorm(n_markers, 0, 0.01)  # Initialize with small random values
  
  # Progress indicator setup
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  
  for (i in 1:n_iter) {
    for (j in 1:n_markers) {
      # Sample beta from the posterior distribution
      beta_proposal <- rnorm(1, beta_current[j], 0.2)  # Increase proposal variance
      
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

# Run Bayesian fine-mapping on simulated data
n_iter <- 10000  # Define the number of iterations
prior_prob <- 0.1  # Increase prior probability for stronger effect
slab_sd <- 1.0  # Increase slab standard deviation for broader prior
beta_samples <- bayesian_finemap(genoprobs, phenotype, prior_prob, slab_sd, n_iter = n_iter)

# Convert to mcmc object for diagnostics
beta_mcmc <- mcmc(beta_samples)

# Summarize the posterior samples
beta_means <- apply(beta_samples, 2, mean)
beta_sds <- apply(beta_samples, 2, sd)
beta_df <- data.frame(Marker = 1:length(beta_means), Mean = beta_means, SD = beta_sds)

# Create a directory to save the plots
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Raw data plot: Phenotype distribution
png(filename = "plots/simulated_phenotype_distribution.png")
ggplot(data.frame(phenotype), aes(x = phenotype)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(title = "Simulated Phenotype Distribution", x = "Phenotype", y = "Frequency") +
  theme_minimal()
dev.off()

# Posterior density plots with true beta lines
posterior_df <- as.data.frame(beta_samples)
posterior_df_long <- tidyr::gather(posterior_df, key = "Marker", value = "Beta")
true_beta_df <- data.frame(Marker = as.factor(rep(1:n_markers, each = n_iter)), Beta = rep(true_beta, each = n_iter))

png(filename = "plots/simulated_posterior_density_plots.png")
ggplot(posterior_df_long, aes(x = Beta, fill = Marker, color = Marker)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = true_beta_df, aes(xintercept = Beta, color = Marker), linetype = "dashed") +
  labs(title = "Posterior Density Plots (Simulated Data) with True Betas", x = "Beta", y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# Trace plots for MCMC diagnostics using ggplot2
traceplot_df <- data.frame(Iteration = 1:nrow(beta_samples), posterior_df)
traceplot_df_long <- tidyr::gather(traceplot_df, key = "Marker", value = "Beta", -Iteration)

png(filename = "plots/simulated_trace_plots.png")
ggplot(traceplot_df_long, aes(x = Iteration, y = Beta, color = Marker)) +
  geom_line(alpha = 0.5) +
  labs(title = "Trace Plots (Simulated Data)", x = "Iteration", y = "Beta") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# Plot true beta coefficients vs. posterior means
png(filename = "plots/true_vs_posterior_means.png")
ggplot(data.frame(True = true_beta, Posterior = beta_means), aes(x = True, y = Posterior)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "True vs Posterior Beta Coefficients", x = "True Beta", y = "Posterior Mean Beta") +
  theme_minimal()
dev.off()

# Posterior means with error bars
png(filename = "plots/simulated_posterior_means_with_error_bars.png")
ggplot(beta_df, aes(x = Marker, y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "red") +
  labs(title = "Posterior Means of Beta Coefficients with Error Bars (Simulated Data)",
       x = "Marker Index", y = "Posterior Mean") +
  theme_minimal()
dev.off()

# Pairwise scatter plots for beta coefficients
png(filename = "plots/simulated_pairwise_scatter_plots.png")
pairs(as.data.frame(beta_samples[, 1:10]), main = "Pairwise Scatter Plots of Beta Coefficients (Simulated Data)")
dev.off()

# Summary statistics of posterior samples
summary_beta_samples <- summary(beta_mcmc)
write.csv(summary_beta_samples$statistics, file = "plots/simulated_posterior_summary_statistics.csv")

# Summary statistics of posterior samples
summary_beta_samples <- summary(beta_mcmc)
write.csv(summary_beta_samples$statistics, file = "plots/simulated_posterior_summary_statistics.csv")

# Plotting autocorrelation diagnostics using the coda package
png(filename = "plots/simulated_autocorrelation_plot.png")
autocorr.plot(beta_mcmc, lag.max = 50, auto.layout = FALSE)
dev.off()

# Plotting Geweke diagnostic
png(filename = "plots/simulated_geweke_plot.png")
geweke.plot(beta_mcmc)
dev.off()

# Plotting Heidelberger-Welch diagnostic
png(filename = "plots/simulated_heidel_plot.png")
heidel.diag(beta_mcmc)
dev.off()
##########################

# Ensure the SuSiE library is installed
if (!requireNamespace("susieR", quietly = TRUE)) {
  install.packages("susieR")
}

# Load the SuSiE library
library(susieR)

# Prepare the data for SuSiE
X <- as.matrix(genoprobs_double)
y <- as.numeric(phenotype)

# Run SuSiE
susie_fit <- susie(X, y)

# Extract posterior means from SuSiE
susie_posterior_means <- susie_fit$beta

# Check the dimensions and length consistency
cat("Length of beta_means: ", length(beta_means), "\n")
cat("Length of susie_posterior_means: ", length(susie_posterior_means), "\n")

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

#############################

# Ensure the SuSiE library is installed
if (!requireNamespace("susieR", quietly = TRUE)) {
  install.packages("susieR")
}

# Load the SuSiE library
library(susieR)

# Prepare the data for SuSiE
X <- as.matrix(genoprobs_double)
y <- as.numeric(phenotype)

# Run SuSiE
susie_fit <- susie(X, y)

# Extract posterior means from SuSiE
susie_posterior_means <- susie_fit$beta

# Check the dimensions and length consistency
cat("Length of beta_means: ", length(beta_means), "\n")
cat("Length of susie_posterior_means: ", length(susie_posterior_means), "\n")

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


