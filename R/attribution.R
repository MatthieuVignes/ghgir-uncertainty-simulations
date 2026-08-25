##################################################
## Variance-based (Sobol) uncertainty attribution on the DAG
## Matthieu Vignes, Aug 2026
##################################################

# MPI's brief (obj 3) asks for a ranked, defensible attribution of uncertainty to activity data, 
# emission factors, and model parameters. We critique "plain Monte Carlo" for only ranking inputs 
# by their own marginal variance; not the same as ranking them by their influence on the output's 
# variance, once shared/correlated structure is accounted for.
# This script estimates first-order and total-order Sobol sensitivity indices for the 8 
# independent root factors of the toy DAG (see ROOT_NAMES) with respect to the national total, 
# using the standard Saltelli/Jansen pick-freeze estimator. The 8 roots are independent by 
# construction (correlation between categories arises from two of them being *shared* by both 
# branches, not from pairwise dependence among the 8). A Sobol decomposition over these roots is 
# well defined and gives a genuine, defensible ranking rather than a marginal-variance heuristic.

sobol_indices <- function(params, n = 20000, seed = 7) {
  # First- (S) and total-order (ST) Sobol indices of each root factor with respect to the national 
  # total, via the Saltelli/Jansen pick-freeze estimator
  set.seed(seed)
  k <- length(ROOT_NAMES)

  A <- matrix(rnorm(n * k), nrow = n, ncol = k, dimnames = list(NULL, ROOT_NAMES))
  B <- matrix(rnorm(n * k), nrow = n, ncol = k, dimnames = list(NULL, ROOT_NAMES))

  fA <- forward(A, params)$national_total
  fB <- forward(B, params)$national_total

  f_all <- c(fA, fB)
  var_total <- var(f_all)
  mean_all <- mean(f_all)

  results <- list()
  for (nm in ROOT_NAMES) {
    AB_i <- A
    AB_i[, nm] <- B[, nm]
    fABi <- forward(AB_i, params)$national_total

    ## Jansen (1999) estimators -- more stable than the original Sobol ones
    s_first <- mean(fB * (fABi - fA)) / var_total
    s_total <- mean((fA - fABi)^2) / (2 * var_total)

    results[[nm]] <- c(S = s_first, ST = s_total)
  }
  attr(results, "meta") <- list(variance = var_total, mean = mean_all, n = n)
  results
}

# Grouping of roots into the 3 MPI-brief categories (activity data/emission factors & methodology/model 
# parameters), used only for presentation. The above estimator works at the individual-root level. This 
# grouping is applied afterwards in the summary bar chart
CATEGORY_GROUPS <- list(
  "Activity data (animal numbers)" = c("log_n_dairy", "log_n_beef"),
  "Activity data (feed/DMI)" = c("eps_dmi_dairy", "eps_dmi_beef", "pasture_quality"),
  "Emission factor / methodology" = c("eps_ch4_dairy", "eps_ch4_beef", "methodology_error")
)

grouped_first_order <- function(sobol_result) {
  # Sum first-order indices within each MPI-brief-style category, for the bar chart
  vapply(CATEGORY_GROUPS, function(names) {
    sum(vapply(names, function(nm) sobol_result[[nm]][["S"]], numeric(1)))
  }, numeric(1))
}
