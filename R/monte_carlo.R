##################################################
## Monte Carlo uncertainty propagation through the toy DAG.
## Matthieu Vignes, Aug 2026
##################################################

# We provide 2 variants:
# * run_correct(), which samples the DAG's roots as specified (pasture quality and 
#   methodology error are shared between the dairy and beef branches), i.e. the 
#   correlation/structure-aware Monte Carlo the proposal argues for
# * run_naive(), which samples every input independently, incl. giving dairy and 
#   beef their own independent draws of pasture quality and methodology error. This 
#   is our "plain Monte Carlo" that implicitly assumes independence and tends to 
#   understate total uncertainty
# Both use the same marginal distributions and the same DAG equations; the only 
# difference is whether the shared roots are actually shared

run_correct <- function(params, n = 50000, seed = 42) {
  z <- draw_z(n, seed = seed)
  forward(z, params, share_pasture_quality = TRUE, share_methodology_error = TRUE)
}

run_naive <- function(params, n = 50000, seed = 42) {
  z <- draw_z(n, seed = seed)
  forward(z, params, share_pasture_quality = FALSE, share_methodology_error = FALSE)
}

summarise_mc <- function(national_total) {
  # Mean, sd, coefficient of variation, and a 95% interval, in tonnes CH4/yr
  m <- mean(national_total)
  s <- sd(national_total)
  ci <- quantile(national_total, c(0.025, 0.975))
  list(
    mean = m,
    sd = s,
    cv_pct = 100 * s / m,
    ci95_lo = unname(ci[1]),
    ci95_hi = unname(ci[2])
  )
}
