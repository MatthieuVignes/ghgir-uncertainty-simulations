##################################################
## Bayesian inference over the causal graph
## Matthieu Vignes, Aug 2026
##################################################

# Monte Carlo and the analytical method both do forward uncertainty propagation by
# pushing prior uncertainty on the roots through the DAG, We report the resulting 
# output uncertainty. Neither method uses data to narrow that uncertainty
# A proper Bayesian treatment of the same graph instead conditions the two shared 
# latent roots (pasture_quality and methodology_error) on synthetic "calibration" 
# observations (e.g. an independent pasture-quality survey, and a small animal-trial 
# comparison of measured vs modelled enteric methane), using conjugate Gaussian 
# updating, and then forward-propagates the resulting (narrower) posterior instead of 
# the prior. This demonstrates a qualitative benefit we claim for a graphical-model 
# framework: data collected anywhere on the graph can tighten the national-total estimate,
# not just the node it was measured at, because pasture_quality and methodology_error are 
# shared parents of both categories.
# All "calibration observations" below are synthetic placeholders, not real measurements, 
# and exist only to demonstrate the mechanism.

calibration_observation <- function(root_name, value, obs_sd) {
  list(root_name = root_name, value = value, obs_sd = obs_sd)
}

conjugate_update <- function(prior_mean, prior_sd, obs) {
  # Standard 1D Gaussian conjugate update: N(prior) x N(likelihood) -> N(posterior)
  prior_var <- prior_sd^2
  obs_var <- obs$obs_sd^2
  post_var <- 1 / (1 / prior_var + 1 / obs_var)
  post_mean <- post_var * (prior_mean / prior_var + obs$value / obs_var)
  c(mean = post_mean, sd = sqrt(post_var))
}

posterior_over_shared_roots <- function(params, pasture_quality_obs, methodology_error_obs) {
  post <- list()
  post$pasture_quality <- if (!is.null(pasture_quality_obs)) {
    conjugate_update(0, params$sd_pasture_quality, pasture_quality_obs)
  } else {
    c(mean = 0, sd = params$sd_pasture_quality)
  }
  post$methodology_error <- if (!is.null(methodology_error_obs)) {
    conjugate_update(0, params$sd_methodology_error, methodology_error_obs)
  } else {
    c(mean = 0, sd = params$sd_methodology_error)
  }
  post
}

run_posterior_predictive <- function(params, pasture_quality_obs, methodology_error_obs,
                                      n = 50000, seed = 42) {
  # Forward-propagate the DAG using the posterior (data-informed) distribution for the 2 
  # shared roots, and the prior for everything else
  post <- posterior_over_shared_roots(params, pasture_quality_obs, methodology_error_obs)
  mean_overrides <- c(pasture_quality = post$pasture_quality[["mean"]],
                       methodology_error = post$methodology_error[["mean"]])
  sd_overrides <- c(pasture_quality = post$pasture_quality[["sd"]],
                     methodology_error = post$methodology_error[["sd"]])

  z <- draw_z(n, seed = seed)
  out <- forward(z, params,
                 share_pasture_quality = TRUE, share_methodology_error = TRUE,
                 root_mean_overrides = mean_overrides, root_sd_overrides = sd_overrides)
  out$posterior <- post
  out
}
