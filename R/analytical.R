##################################################
## Analytical (closed-form) error propagation through the toy DAG.
## Matthieu Vignes, Aug 2026
##################################################

# Every node up to the category-emissions stage is a *linear* function of Gaussian 
# roots (log-space). Hence each category's emissions is log-normal. We can get its 
# mean/variance in closed form (the delta-method / "analytical error propagation" 
# MPI's brief asks to compare against Monte Carlo)

# Final step an approximation: national total = dairy emissions + beef emissions, 
# which is a sum of two correlated log-normals, not log-normal in general. Hence we:
# 1. get the mean and variance of the national total exactly (Var(X+Y) = Var(X) + 
#    Var(Y) + 2 Cov(X, Y); cov of joint log-normal variables has closed form)
# 2. approximate the *distribution* of the total (for interval estimates) by 
#    moment-matching a log-normal to that mean/variance (Fenton-Wilkinson approx)
# Analytical method's key trade-off: exact first- and second-moment matching, but an 
# approximate tail/CI, whereas Monte Carlo gets exact shape (but has sampling error)

category_log_moments <- function(params, category, share_pasture_quality, share_methodology_error) {
  # Linear coefficients of log(category emissions) on the 8 standardised root factors z 
  # each root contributes coefficient * sd_root * z_root -> exact mean and variance for the log of that 
  # category's emissions.
  sd <- root_sd_vector(params)
  coef <- setNames(rep(0, length(ROOT_NAMES)), ROOT_NAMES)

  if (category == "dairy") {
    a0 <- params$a0_dairy; a1 <- params$a1_dairy
    b0 <- params$b0_dairy; b1 <- params$b1_dairy; b2 <- params$b2_dairy
    mean_log_n <- params$mean_log_n_dairy
    dmi_key <- "eps_dmi_dairy"; ch4_key <- "eps_ch4_dairy"; n_key <- "log_n_dairy"
  } else {
    a0 <- params$a0_beef; a1 <- params$a1_beef
    b0 <- params$b0_beef; b1 <- params$b1_beef; b2 <- params$b2_beef
    mean_log_n <- params$mean_log_n_beef
    dmi_key <- "eps_dmi_beef"; ch4_key <- "eps_ch4_beef"; n_key <- "log_n_beef"
  }
  pq_root <- "pasture_quality"
  me_root <- "methodology_error"

  const <- mean_log_n + b0 + b1 * a0

  coef[[pq_root]] <- coef[[pq_root]] + if (share_pasture_quality) b1 * a1 else 0
  coef[[dmi_key]] <- coef[[dmi_key]] + b1
  coef[[me_root]] <- coef[[me_root]] + if (share_methodology_error) b2 else 0
  coef[[ch4_key]] <- coef[[ch4_key]] + 1.0
  coef[[n_key]]   <- coef[[n_key]]   + 1.0

  # If not shared, category still has its own independent draw of the same marginal variance. We add 
  # an extra category-specific adequate variance component with the same coefficient*sd magnitude the 
  # shared root would have had.
  extra_var <- 0
  if (!share_pasture_quality) extra_var <- extra_var + (b1 * a1 * params$sd_pasture_quality)^2
  if (!share_methodology_error) extra_var <- extra_var + (b2 * params$sd_methodology_error)^2

  weighted <- coef * sd
  var <- sum(weighted^2) + extra_var
  mean <- const

  list(mean = mean, var = var, coef = coef, sd = sd)
}

national_total_moments <- function(params, shared) {
  # Exact mean and variance (natural scale, tonnes CH4/yr) of national total, plus a Fenton-Wilkinson 
  # log-normal approximation for the CI.
  const <- params$days_per_year * 1e-3   # kg -> t and per-day -> per-year (outside logs)

  d <- category_log_moments(params, "dairy", shared, shared)
  b <- category_log_moments(params, "beef", shared, shared)

  # covariance of the 2 category log-emissions: only shared roots contribute; independ roots have 0 cross term
  cov_log <- if (shared) sum(d$coef * b$coef * d$sd^2) else 0

  mean_d <- const * exp(d$mean + d$var / 2)
  mean_b <- const * exp(b$mean + b$var / 2)
  var_nat_d <- (exp(d$var) - 1) * exp(2 * d$mean + d$var) * const^2
  var_nat_b <- (exp(b$var) - 1) * exp(2 * b$mean + b$var) * const^2
  cov_nat <- (exp(cov_log) - 1) * exp(d$mean + b$mean + (d$var + b$var) / 2) * const^2

  total_mean <- mean_d + mean_b
  total_var <- var_nat_d + var_nat_b + 2 * cov_nat

  # Fenton-Wilkinson: match a log-normal to (total_mean, total_var)
  fw_var_log <- log(1 + total_var / total_mean^2)
  fw_mu_log <- log(total_mean) - fw_var_log / 2
  ci <- exp(qnorm(c(0.025, 0.975), mean = fw_mu_log, sd = sqrt(fw_var_log)))

  list(
    mean = total_mean,
    sd = sqrt(total_var),
    cv_pct = 100 * sqrt(total_var) / total_mean,
    ci95_lo = ci[1],
    ci95_hi = ci[2],
    cov_between_categories = cov_nat
  )
}
