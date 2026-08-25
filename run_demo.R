##################################################
## Attributing uncertatiny in Aotearoa/NZ agricultural GHG emissions simulations
## Matthieu Vignes, Aug 2026
##################################################

# This script runs the full preliminary demo: DAG diagram, three-way method comparison (Monte 
# Carlo / analytical / Bayesian), and Sobol-based uncertainty attribution

# Output figures in figures/*.png, a printed summary: results/summary.txt and a per 
# variable uncertatiny contribution measure in results/sobol.csv 

# See the README for more details

suppressWarnings(suppressMessages({
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) "."
  )
}))
if (is.null(script_dir) || is.na(script_dir) || script_dir == "") script_dir <- "."
setwd(script_dir)

for (f in c("dag_model.R", "monte_carlo.R", "analytical.R", "bayesian.R",
            "attribution.R", "plotting.R")) {
  source(file.path("R", f))
}

dir.create("figures")
dir.create("results")

params <- default_params()

cat(strrep("=", 70), "\n")
cat("1. Monte Carlo propagation (correct vs naive-independent)\n")
cat(strrep("=", 70), "\n")
mc_correct <- run_correct(params)
mc_naive <- run_naive(params)
s_correct <- summarise_mc(mc_correct$national_total)
s_naive <- summarise_mc(mc_naive$national_total)
cat(sprintf("  Correlation-aware MC : mean=%.0f  sd=%.0f  CV=%.1f%%  95%%CI=(%.0f, %.0f)\n",
            s_correct$mean, s_correct$sd, s_correct$cv_pct, s_correct$ci95_lo, s_correct$ci95_hi))
cat(sprintf("  Naive independent MC : mean=%.0f  sd=%.0f  CV=%.1f%%  95%%CI=(%.0f, %.0f)\n",
            s_naive$mean, s_naive$sd, s_naive$cv_pct, s_naive$ci95_lo, s_naive$ci95_hi))
understatement <- 100 * (1 - s_naive$sd / s_correct$sd)
cat(sprintf("  -> Ignoring the shared roots understates the total-uncertainty sd by ~%.0f%% in this toy example.\n",
            understatement))

cat("\n", strrep("=", 70), "\n", sep = "")
cat("2. Analytical (delta-method / Fenton-Wilkinson) propagation\n")
cat(strrep("=", 70), "\n")
an_correct <- national_total_moments(params, shared = TRUE)
an_naive <- national_total_moments(params, shared = FALSE)
cat(sprintf("  Correlation-aware analytical : mean=%.0f  sd=%.0f  CV=%.1f%%  95%%CI=(%.0f, %.0f)\n",
            an_correct$mean, an_correct$sd, an_correct$cv_pct, an_correct$ci95_lo, an_correct$ci95_hi))
cat(sprintf("  Naive independent analytical : mean=%.0f  sd=%.0f  CV=%.1f%%  95%%CI=(%.0f, %.0f)\n",
            an_naive$mean, an_naive$sd, an_naive$cv_pct, an_naive$ci95_lo, an_naive$ci95_hi))
cat(sprintf("  -> Analytical sd matches Monte Carlo sd to within %.1f%% (exact first/second moments; CI uses a log-normal moment-matching approximation).\n",
            100 * abs(an_correct$sd - s_correct$sd) / s_correct$sd))

cat("\n", strrep("=", 70), "\n", sep = "")
cat("3. Bayesian inference over the graph (posterior predictive)\n")
cat(strrep("=", 70), "\n")
pq_obs <- calibration_observation("pasture_quality", value = 0.03, obs_sd = 0.025)
me_obs <- calibration_observation("methodology_error", value = -0.02, obs_sd = 0.03)
bayes_out <- run_posterior_predictive(params, pq_obs, me_obs)
s_bayes <- summarise_mc(bayes_out$national_total)
post <- bayes_out$posterior
cat(sprintf("  Synthetic calibration data -> posterior pasture_quality: mean=%.3f, sd=%.3f (prior sd=%.3f)\n",
            post$pasture_quality[["mean"]], post$pasture_quality[["sd"]], params$sd_pasture_quality))
cat(sprintf("  Synthetic calibration data -> posterior methodology_error: mean=%.3f, sd=%.3f (prior sd=%.3f)\n",
            post$methodology_error[["mean"]], post$methodology_error[["sd"]], params$sd_methodology_error))
cat(sprintf("  Posterior-predictive national total: mean=%.0f  sd=%.0f  CV=%.1f%%  95%%CI=(%.0f, %.0f)\n",
            s_bayes$mean, s_bayes$sd, s_bayes$cv_pct, s_bayes$ci95_lo, s_bayes$ci95_hi))
tightening <- 100 * (1 - s_bayes$sd / s_correct$sd)
cat(sprintf("  -> Conditioning on (synthetic) calibration data narrows the total-uncertainty sd by ~%.0f%% relative to the prior forward propagation.\n",
            tightening))

cat("\n", strrep("=", 70), "\n", sep = "")
cat("4. Uncertainty attribution (Sobol indices)\n")
cat(strrep("=", 70), "\n")
sobol <- sobol_indices(params)
for (nm in ROOT_NAMES) {
  cat(sprintf("  %-24s  S=%.3f   ST=%.3f\n", nm, sobol[[nm]][["S"]], sobol[[nm]][["ST"]]))
}
grouped <- grouped_first_order(sobol)
cat("  Grouped (MPI-brief-style categories):\n")
for (nm in names(sort(grouped, decreasing = TRUE))) {
  cat(sprintf("    %-38s %.3f\n", nm, grouped[[nm]]))
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("5. Figures\n")
cat(strrep("=", 70), "\n")
plot_dag(file.path("figures", "dag_structure.png"))
plot_method_comparison(mc_correct, mc_naive, an_correct, an_naive, bayes_out,
                        file.path("figures", "method_comparison.png"))
plot_attribution(grouped, file.path("figures", "uncertainty_attribution.png"))
cat("  Saved figures to figures/\n")

# Writes results
sink(file.path("results", "summary.txt"))
cat("GHGIR toy demo -- results summary (synthetic data; see README.md)\n\n")
cat("Monte Carlo (correlation-aware):\n"); str(s_correct)
cat("\nMonte Carlo (naive/independent):\n"); str(s_naive)
cat("\nAnalytical (correlation-aware):\n"); str(an_correct)
cat("\nAnalytical (naive/independent):\n"); str(an_naive)
cat("\nBayesian posterior predictive:\n"); str(s_bayes)
cat(sprintf("\nsd understatement if naive: %.1f%%\n", understatement))
cat(sprintf("sd tightening from calibration: %.1f%%\n", tightening))
sink()

sobol_df <- data.frame(
  root = ROOT_NAMES,
  S = vapply(ROOT_NAMES, function(nm) sobol[[nm]][["S"]], numeric(1)),
  ST = vapply(ROOT_NAMES, function(nm) sobol[[nm]][["ST"]], numeric(1))
)
write.csv(sobol_df, file.path("results", "sobol_indices.csv"), row.names = FALSE)
