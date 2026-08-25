##################################################
## Plotting functions for the demo figures
## Matthieu Vignes, Aug 2026
##################################################

GREEN <- "#1F5C3D"
RED <- "#B23A48"
BLUE <- "#3A6EA5"

plot_method_comparison <- function(mc_correct, mc_naive, analytical_correct,
                                    analytical_naive, bayesian_out, path) {
  png(path, width = 1100, height = 700, res = 130)
  on.exit(dev.off())

  x_all <- c(mc_correct$national_total, mc_naive$national_total, bayesian_out$national_total)
  xlim <- range(x_all) * c(0.9, 1.05)

  h1 <- hist(mc_correct$national_total, breaks = 80, plot = FALSE)
  h2 <- hist(mc_naive$national_total, breaks = 80, plot = FALSE)
  h3 <- hist(bayesian_out$national_total, breaks = 80, plot = FALSE)
  ylim <- c(0, max(h1$density, h2$density, h3$density) * 1.15)

  plot(h1, freq = FALSE, col = adjustcolor(GREEN, alpha.f = 0.45), border = NA,
       xlim = xlim, ylim = ylim, main = "Comparing uncertainty-propagation methods on the same toy DAG",
       xlab = "National total, illustrative t CH4/yr (synthetic units)", ylab = "Density")
  plot(h2, freq = FALSE, col = adjustcolor(RED, alpha.f = 0.45), border = NA, add = TRUE)
  plot(h3, freq = FALSE, col = adjustcolor(BLUE, alpha.f = 0.45), border = NA, add = TRUE)

  x <- seq(xlim[1], xlim[2], length.out = 400)
  add_lognormal_curve <- function(res, colour) {
    var_log <- log(1 + res$sd^2 / res$mean^2)
    mu_log <- log(res$mean) - var_log / 2
    dens <- dlnorm(x, meanlog = mu_log, sdlog = sqrt(var_log))
    lines(x, dens, col = colour, lwd = 2, lty = 2)
  }
  add_lognormal_curve(analytical_correct, GREEN)
  add_lognormal_curve(analytical_naive, RED)

  legend("topright", bty = "n", cex = 0.75,
         legend = c("Monte Carlo (correlation-aware)", "Monte Carlo (naive, independent)",
                    "Bayesian posterior predictive",
                    "Analytical / Fenton-Wilkinson (correlation-aware)",
                    "Analytical / Fenton-Wilkinson (naive)"),
         fill = c(adjustcolor(GREEN, 0.45), adjustcolor(RED, 0.45), adjustcolor(BLUE, 0.45), NA, NA),
         border = c(NA, NA, NA, NA, NA),
         lty = c(NA, NA, NA, 2, 2), lwd = c(NA, NA, NA, 2, 2),
         col = c(NA, NA, NA, GREEN, RED),
         merge = TRUE)
}

plot_attribution <- function(grouped_indices, path) {
  png(path, width = 900, height = 550, res = 130)
  on.exit(dev.off())

  ord <- order(grouped_indices)
  vals <- grouped_indices[ord]
  names_wrapped <- names(vals)

  par(mar = c(4.5, 13, 3.5, 2))
  bp <- barplot(vals, horiz = TRUE, col = GREEN, border = NA, las = 1,
                names.arg = names_wrapped, xlab = "First-order Sobol index (share of output variance)",
                main = "Uncertainty attribution: ranked contributors\nto national-total variance",
                xlim = c(0, max(vals) * 1.15))
}

plot_dag <- function(path) {
  png(path, width = 1300, height = 700, res = 130)
  on.exit(dev.off())

  nodes <- list(
    pq   = list(label = "Pasture quality\n(shared)", x = 0.5, y = 3.0, shared = TRUE),
    me   = list(label = "Methodology error\n(shared)", x = 0.5, y = 1.0, shared = TRUE),
    nd   = list(label = "Animal numbers\ndairy", x = 2.0, y = 5.6, shared = FALSE),
    nb   = list(label = "Animal numbers\nbeef", x = 2.0, y = -0.6, shared = FALSE),
    dmid = list(label = "DMI dairy", x = 2.0, y = 4.0, shared = FALSE),
    dmib = list(label = "DMI beef", x = 2.0, y = 2.0, shared = FALSE),
    chd  = list(label = "CH4/animal\ndairy", x = 3.8, y = 4.0, shared = FALSE),
    chb  = list(label = "CH4/animal\nbeef", x = 3.8, y = 1.0, shared = FALSE),
    ed   = list(label = "Emissions\ndairy", x = 5.4, y = 4.0, shared = FALSE),
    eb   = list(label = "Emissions\nbeef", x = 5.4, y = 1.0, shared = FALSE),
    nat  = list(label = "National\ntotal", x = 7.0, y = 2.5, shared = FALSE)
  )
  edges <- list(
    c("pq", "dmid"), c("pq", "dmib"),
    c("dmid", "chd"), c("dmib", "chb"),
    c("me", "chd"), c("me", "chb"),
    c("nd", "ed"), c("chd", "ed"),
    c("nb", "eb"), c("chb", "eb"),
    c("ed", "nat"), c("eb", "nat")
  )

  par(mar = c(1, 1, 3, 1))
  plot(NA, xlim = c(-0.5, 8), ylim = c(-1.5, 6.5), axes = FALSE, xlab = "", ylab = "",
       main = "Toy inventory DAG: red nodes are shared across the dairy/beef branches\n(illustrative structure only; see README)",
       cex.main = 0.95)

  r <- 0.55  # node "radius" in plot units, for edge trimming
  for (e in edges) {
    n1 <- nodes[[e[1]]]; n2 <- nodes[[e[2]]]
    dx <- n2$x - n1$x; dy <- n2$y - n1$y
    dist <- sqrt(dx^2 + dy^2)
    ux <- dx / dist; uy <- dy / dist
    x0 <- n1$x + ux * r; y0 <- n1$y + uy * r
    x1 <- n2$x - ux * r; y1 <- n2$y - uy * r
    arrows(x0, y0, x1, y1, length = 0.1, col = "#555555", lwd = 1.4)
  }

  for (nm in names(nodes)) {
    nd <- nodes[[nm]]
    col <- if (nd$shared) RED else GREEN
    symbols(nd$x, nd$y, circles = r * 0.9, inches = FALSE, add = TRUE,
            bg = col, fg = "white", lwd = 1.5)
    text(nd$x, nd$y, labels = nd$label, col = "white", cex = 0.62, font = 2)
  }
}

