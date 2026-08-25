##################################################
## Toy causal-graphical model of an agricultural GHG inventory chain.
## Matthieu Vignes, Aug 2026
##################################################

# This is a synthetic, illustrative model, which mimics a plausible qualitative 
# structure: pasture quality -> dry matter intake (DMI) -> enteric methane
#                -> emission factor -> national aggregation, for 2 livestock 
# categories (dairy, beef) that:
# (a) draw on a shared pasture-quality signal, and
# (b) share a common "methodology" error term, because both categories' emission 
#     factors are derived from the same underlying Tier 2 equation form.
# No real MPI / AIM data, parameters, or outputs are used anywhere in this work. 
# Parameter values below are illustrative and were chosen to produce a plausible
# demonstration; they must not be read as estimates of actual NZ agricultural 
# emissions or uncertainty.
# We show a methodological point: how correlated/shared uncertainty changes the total 
# uncertainty estimate and the attribution ranking relative to treating every input as 
# independent on a (simple enough to inspect fully) structure

# Root ("exogenous") variables of the DAG.
# All modelled as a Gaussian in log-space (or, for the 2 activity-data nodes, directly 
# as a log-normal count). All 8 roots are independent of one another a priori, hence we  
# double this vector as the factor list for the Sobol-style attribution analysis

ROOT_NAMES <- c(
  "pasture_quality", # shared regional/national pasture ME signal
  "methodology_error", # shared Tier-2 equation-form uncertainty
  "eps_dmi_dairy", # dairy-specific DMI residual
  "eps_dmi_beef", # beef-specific DMI residual
  "eps_ch4_dairy", # dairy-specific enteric-CH4 residual
  "eps_ch4_beef", # beef-specific enteric-CH4 residual
  "log_n_dairy", # dairy animal-numbers activity data (log head)
  "log_n_beef" # beef animal-numbers activity data (log head)
)

default_params <- function() {
  list(# root sd (log-scale unless noted)
    sd_pasture_quality = 0.06,
    sd_methodology_error = 0.08,
    sd_eps_dmi_dairy = 0.05,
    sd_eps_dmi_beef = 0.07,
    sd_eps_ch4_dairy = 0.04,
    sd_eps_ch4_beef = 0.05,
    sd_log_n_dairy = 0.03,
    sd_log_n_beef = 0.05,
    # DMI regression on pasture quality: log(DMI) = a0 + a1*PQ + eps
    a0_dairy = log(18.0), # kg DM/animal/day
    a1_dairy = 0.9,
    a0_beef = log(9.0),
    a1_beef = 0.6,
    # CH4-per-animal regression: log(CH4) = b0 + b1*log(DMI) + b2*ME_err + eps
    # CH4 (kg/day) = exp(b0) * DMI^b1, i.e. exp(b0) plays the role of a (simplified) methane 
    # conversion rate-like scale factor
    b0_dairy = log(0.33/18.0),
    b1_dairy = 1.0,
    b2_dairy = 1.0,
    b0_beef = log(0.15/9.0),
    b1_beef = 1.0,
    b2_beef = 1.0,
    # activity data: mean log animal numbers
    mean_log_n_dairy = log(6.2e6), # national dairy herd scale
    mean_log_n_beef = log(3.7e6), # national beef herd scale
    # days/year multiplier to annualise per-day CH4
    days_per_year = 365.0
  )
}

ROOT_SD_KEY <- c(
  pasture_quality = "sd_pasture_quality",
  methodology_error = "sd_methodology_error",
  eps_dmi_dairy = "sd_eps_dmi_dairy",
  eps_dmi_beef = "sd_eps_dmi_beef",
  eps_ch4_dairy = "sd_eps_ch4_dairy",
  eps_ch4_beef = "sd_eps_ch4_beef",
  log_n_dairy = "sd_log_n_dairy",
  log_n_beef = "sd_log_n_beef"
)

root_sd_vector <- function(params) {
  vapply(ROOT_NAMES, function(nm) params[[ROOT_SD_KEY[[nm]]]], numeric(1))
}

# Runs `expr` under a temporary, fixed RNG seed without changing the caller 
# global RNG state
with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv))
  }
  set.seed(seed)
  force(expr)
}

draw_z <- function(n, seed = NULL) {
  # Draws n standard-normal samples for the 8 root factors
  if (!is.null(seed)) set.seed(seed)
  matrix(rnorm(n * length(ROOT_NAMES)), nrow = n, ncol = length(ROOT_NAMES),
         dimnames = list(NULL, ROOT_NAMES))
}

forward <- function(z, params,
                     share_pasture_quality = TRUE,
                     share_methodology_error = TRUE,
                     root_mean_overrides = NULL,
                     root_sd_overrides = NULL) {
  # Deterministic forward pass through the DAG given standardised root draws
  # z: n x 8 matrix of standard-normal draws for the 8 root factors (columns named as in 
  #  ROOT_NAMES). Scaling by the root sd done inside this function.
  # share_pasture_quality, share_methodology_error: if FALSE, the dairy and beef branches 
  #  each get an indep draw for that root instead of sharing one. This reproduces the 
  #  "naive" plain-Monte-Carlo we criticise and compare to the correct, correlation-corrected model
  # root_mean_overrides, root_sd_overrides: optional named numeric vectors to give specific roots 
  #  a non-default mean/sd (e.g. a Bayesian posterior mean/sd instead of the prior mean 0/prior sd)
  # Returns a list of numerics (each length n): all intermediate and final node values. National 
  #  totals are in tonnes CH4/year
  n <- nrow(z)
  sd <- root_sd_vector(params)
  mean <- setNames(rep(0, length(ROOT_NAMES)), ROOT_NAMES)
  if (!is.null(root_sd_overrides)) {
    for (nm in names(root_sd_overrides)) sd[[nm]] <- root_sd_overrides[[nm]]
  }
  if (!is.null(root_mean_overrides)) {
    for (nm in names(root_mean_overrides)) mean[[nm]] <- root_mean_overrides[[nm]]
  }
  roots <- lapply(ROOT_NAMES, function(nm) mean[[nm]] + z[, nm] * sd[[nm]])
  names(roots) <- ROOT_NAMES

  if (share_pasture_quality) {
    pq_dairy <- roots$pasture_quality
    pq_beef <- roots$pasture_quality
  } else {
    pq_dairy <- roots$pasture_quality
    pq_beef <- with_seed(0, rnorm(n, 0, params$sd_pasture_quality))
  }

  if (share_methodology_error) {
    me_dairy <- roots$methodology_error
    me_beef <- roots$methodology_error
  } else {
    me_dairy <- roots$methodology_error
    me_beef <- with_seed(1, rnorm(n, 0, params$sd_methodology_error))
  }

  # DMI
  log_dmi_dairy <- params$a0_dairy + params$a1_dairy * pq_dairy + roots$eps_dmi_dairy
  log_dmi_beef  <- params$a0_beef  + params$a1_beef  * pq_beef  + roots$eps_dmi_beef

  # enteric CH4 per animal per day
  log_ch4_dairy <- params$b0_dairy + params$b1_dairy * log_dmi_dairy +
    params$b2_dairy * me_dairy + roots$eps_ch4_dairy
  log_ch4_beef <- params$b0_beef + params$b1_beef * log_dmi_beef +
    params$b2_beef * me_beef + roots$eps_ch4_beef

  # activity data: animal numbers
  log_n_dairy <- params$mean_log_n_dairy + roots$log_n_dairy
  log_n_beef  <- params$mean_log_n_beef  + roots$log_n_beef

  # category emissions, tonnes CH4/year
  kg_to_t <- 1e-3
  emis_dairy <- exp(log_n_dairy + log_ch4_dairy) * params$days_per_year * kg_to_t
  emis_beef  <- exp(log_n_beef  + log_ch4_beef)  * params$days_per_year * kg_to_t

  national_total <- emis_dairy + emis_beef

  list(
    pasture_quality_dairy = pq_dairy,
    pasture_quality_beef = pq_beef,
    methodology_error_dairy = me_dairy,
    methodology_error_beef = me_beef,
    log_dmi_dairy = log_dmi_dairy,
    log_dmi_beef = log_dmi_beef,
    log_ch4_dairy = log_ch4_dairy,
    log_ch4_beef = log_ch4_beef,
    emissions_dairy = emis_dairy,
    emissions_beef = emis_beef,
    national_total = national_total
  )
}
