############################################################
## Power analysis -- "Are Altruists Parochial?"
## Stage 1 Registered Report, Nature Human Behaviour
##
## 1. Power of the confirmatory test (OLS of beta on theta)
## 2. Power of the equivalence test (TOST) under b = 0
############################################################

out_dir <- here::here("power_analysis", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

library(tidyverse)

alpha        <- 0.025
n_sims       <- 10000
target_power <- 0.95

beta_mean <- 0.30
beta_sd   <- 0.291
sd_theta  <- 0.30

n_grid <- seq(100, 2000, by = 100)
b_grid <- seq(0, 0.20, by = 0.01)

set.seed(123)

# ---- utilities ----

a_from_sd_theta <- function(sd_theta) {
  (1 / sd_theta^2 - 1) / 2
}

rtheta_beta <- function(n, sd_theta) {
  a <- a_from_sd_theta(sd_theta)
  2 * rbeta(n, a, a) - 1
}

phi_from_mean_sd <- function(mu, sd) {
  mu * (1 - mu) / sd^2 - 1
}

phi_beta <- phi_from_mean_sd(beta_mean, beta_sd)

clamp01 <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

simulate_main <- function(n, b) {
  theta <- rtheta_beta(n, sd_theta)
  mu    <- clamp01(beta_mean + b * theta)
  beta  <- rbeta(n, mu * phi_beta, (1 - mu) * phi_beta)
  tibble(beta, theta)
}

one_p_main <- function(n, b) {
  df <- simulate_main(n, b)
  coef(summary(lm(beta ~ theta, data = df)))["theta", "Pr(>|t|)"]
}

power_main <- function(n, b) {
  mean(replicate(n_sims, one_p_main(n, b)) < alpha)
}

############################################################
## 1. Power of the confirmatory test
############################################################

grid <- expand_grid(n = n_grid, b = b_grid)

message("Running ", nrow(grid), " (n,b) combinations x ", n_sims, " sims each ...")

results <- grid %>%
  mutate(
    power = map2_dbl(n, b, ~ {
      message("main: n=", .x, " b=", .y)
      power_main(.x, .y)
    })
  )

write_csv(results, file.path(out_dir, "power_curves_main_sdtheta_0.30.csv"))

# ---- MDE at 95% power for the two analysed sample sizes ----

mde <- results %>%
  filter(n %in% c(1000, 2000)) %>%
  group_by(n) %>%
  summarise(
    b = approx(x = power, y = b, xout = target_power, ties = "ordered")$y,
    .groups = "drop"
  ) %>%
  mutate(
    label = c("Within-country\n(N = 1,000)", "Cross-country\n(N = 2,000)"),
    hjust = c(-0.1, 1),
    vjust = c(-0.4, -1)
  )

print(mde)

# ---- heatmap power surface (Figure 3) ----

p_heat <- ggplot(results, aes(x = n, y = b, fill = power)) +
  geom_tile(color = "white", linewidth = 0) +
  geom_vline(xintercept = seq(100, 2000, 100),
             color = "white", alpha = 0.15, linewidth = 0.2) +
  geom_hline(yintercept = seq(0, 0.20, 0.02),
             color = "white", alpha = 0.15, linewidth = 0.2) +
  geom_contour(aes(z = power), breaks = target_power,
               color = "red", linewidth = 1.2) +
  annotate("text", x = 700, y = 0.19, label = "95% power", color = "red") +
  geom_segment(data = mde, aes(x = n, xend = n, y = -Inf, yend = b),
               inherit.aes = FALSE, linetype = "dashed", color = "black") +
  geom_segment(data = mde, aes(x = -Inf, xend = n, y = b, yend = b),
               inherit.aes = FALSE, linetype = "dashed", color = "black") +
  geom_point(data = mde, aes(x = n, y = b),
             inherit.aes = FALSE, size = 3, color = "black") +
  geom_text(data = mde,
            aes(x = n, y = b, label = label, hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, size = 3.5, lineheight = 0.9, color = "black") +
  scale_fill_viridis_c(name = "Power", limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 2000, by = 200)) +
  scale_y_continuous(breaks = seq(0, 0.20, by = 0.02)) +
  labs(x = "Sample size (N)", y = "Effect size (b)") +
  theme_minimal(base_size = 14)

print(p_heat)

ggsave(file.path(out_dir, "p_heat_main.png"), p_heat, width = 6, height = 6, dpi = 300)

############################################################
## 2. Power of the equivalence test (TOST) under b = 0
## Each one-sided test at alpha = 0.025.
############################################################

set.seed(456)

n_tost     <- c(1000, 2000)                 # within-country / cross-country
delta_grid <- seq(0.001, 0.20, by = 0.001)  # candidate bounds (raw scale)

fit_null_once <- function(n) {
  df  <- simulate_main(n, b = 0)
  fit <- summary(lm(beta ~ theta, data = df))
  tibble(
    b_hat = fit$coefficients["theta", "Estimate"],
    se    = fit$coefficients["theta", "Std. Error"],
    dof   = fit$df[2]
  )
}

message("TOST: simulating under b = 0 ...")

null_fits <- map_dfr(n_tost, function(n) {
  message("TOST: n = ", n)
  map_dfr(seq_len(n_sims), ~ fit_null_once(n)) %>%
    mutate(n = n)
})

tost_power <- expand_grid(n = n_tost, delta = delta_grid) %>%
  mutate(
    power = map2_dbl(n, delta, function(nn, d) {
      f     <- filter(null_fits, n == nn)
      tcrit <- qt(1 - alpha, f$dof)
      mean((f$b_hat + d) / f$se > tcrit &
             (f$b_hat - d) / f$se < -tcrit)
    })
  )

write_csv(tost_power, file.path(out_dir, "tost_power_by_delta.csv"))

# ---- smallest bound achieving >= 95% power (raw and SD units) ----

sd_ratio <- sd_theta / beta_sd

delta_95 <- tost_power %>%
  group_by(n) %>%
  filter(power >= target_power) %>%
  slice_min(delta, n = 1) %>%
  ungroup() %>%
  mutate(delta_std = delta * sd_ratio)

print(delta_95)

# ---- TOST power curves ----

p_tost <- ggplot(tost_power, aes(x = delta, y = power,
                                 color = factor(n), group = n)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = target_power, linetype = "dashed", color = "grey40") +
  geom_vline(data = delta_95, aes(xintercept = delta, color = factor(n)),
             linetype = "dotted", linewidth = 0.8, show.legend = FALSE) +
  geom_point(data = delta_95, aes(x = delta, y = target_power),
             size = 3, show.legend = FALSE) +
  geom_text(data = delta_95,
            aes(x = delta, y = 0.05, label = sprintf("Delta = %.3f", delta)),
            hjust = -0.1, show.legend = FALSE, size = 3.5) +
  scale_color_viridis_d(name = "Sample size (N)", end = 0.7) +
  scale_x_continuous(breaks = seq(0, 0.20, by = 0.02)) +
  labs(
    x = expression("Equivalence bound " * Delta * " (raw scale)"),
    y = "Power to declare equivalence (true b = 0)"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

print(p_tost)

ggsave(file.path(out_dir, "p_tost_power.png"), p_tost, width = 8, height = 5, dpi = 300)
