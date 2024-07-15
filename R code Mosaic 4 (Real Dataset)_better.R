# Load necessary libraries
library(qtl)
library(MCMCpack)
library(coda)
library(ggplot2)
library(missForest)
library(tidyr)

# Load the hyper dataset
data(hyper)

# Calculate genotype probabilities
hyper <- calc.genoprob(hyper, step = 1)

# Initial QTL analysis
qtl_result <- scanone(hyper, method = "hk")

# Plot QTL results
plot(qtl_result)

# Extract genotype probabilities and phenotypes
genoprobs <- pull.geno(hyper)
phenotype <- hyper$pheno[,1]

# Convert genoprobs to double precision
genoprobs_double <- apply(genoprobs, 2, as.double)

# Impute missing values
if (sum(is.na(genoprobs_double)) > 0) {
  genoprobs_double <- missForest(genoprobs_double)$ximp
}

# Bayesian fine-mapping function
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

# Run Bayesian fine-mapping
n_iter <- 10000
prior_prob <- 0.1
slab_sd <- 1.0
beta_samples <- bayesian_finemap(genoprobs_double, phenotype, prior_prob, slab_sd, n_iter)
beta_mcmc <- mcmc(beta_samples)

####################


# Create a directory to save the plots if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Raw data plot: Phenotype distribution
png(filename = "plots/hyper_phenotype_distribution.png")
ggplot(data.frame(phenotype), aes(x = phenotype)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(title = "Phenotype Distribution in the hyper Dataset", x = "Phenotype", y = "Frequency") +
  theme_minimal()
dev.off()

# Posterior density plots for beta coefficients
posterior_df <- as.data.frame(beta_samples)
posterior_df_long <- tidyr::gather(posterior_df, key = "Marker", value = "Beta")

png(filename = "plots/hyper_posterior_density_plots.png")
ggplot(posterior_df_long, aes(x = Beta, fill = Marker, color = Marker)) +
  geom_density(alpha = 0.5) +
  labs(title = "Posterior Density Plots of Beta Coefficients", x = "Beta", y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# Trace plots for MCMC diagnostics using ggplot2
traceplot_df <- data.frame(Iteration = 1:nrow(beta_samples), posterior_df)
traceplot_df_long <- tidyr::gather(traceplot_df, key = "Marker", value = "Beta", -Iteration)

png(filename = "plots/hyper_trace_plots.png")
ggplot(traceplot_df_long, aes(x = Iteration, y = Beta, color = Marker)) +
  geom_line(alpha = 0.5) +
  labs(title = "Trace Plots of MCMC Samples", x = "Iteration", y = "Beta") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# Posterior means with error bars
beta_means <- apply(beta_samples, 2, mean)
beta_sds <- apply(beta_samples, 2, sd)
beta_df <- data.frame(Marker = 1:length(beta_means), Mean = beta_means, SD = beta_sds)

png(filename = "plots/hyper_posterior_means_with_error_bars.png")
ggplot(beta_df, aes(x = Marker, y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "red") +
  labs(title = "Posterior Means of Beta Coefficients with Error Bars", x = "Marker Index", y = "Posterior Mean") +
  theme_minimal()
dev.off()

# Pairwise scatter plots for beta coefficients
png(filename = "plots/hyper_pairwise_scatter_plots.png")
pairs(as.data.frame(beta_samples[, 1:10]), main = "Pairwise Scatter Plots of Beta Coefficients")
dev.off()

# Plotting autocorrelation diagnostics using the coda package
png(filename = "plots/hyper_autocorrelation_plot.png")
autocorr.plot(beta_mcmc, lag.max = 50, auto.layout = TRUE)
dev.off()

# Plotting Geweke diagnostic
png(filename = "plots/hyper_geweke_plot.png")
geweke.plot(beta_mcmc)
dev.off()

# Plotting Heidelberger-Welch diagnostic
png(filename = "plots/hyper_heidel_plot.png")
heidel.diag(beta_mcmc)
dev.off()

