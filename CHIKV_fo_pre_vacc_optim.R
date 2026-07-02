library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(splines)
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# ------------------------------------------------------------
# 1. SEIR function (unchanged in structure, takes any base_beta)
# ------------------------------------------------------------
seir_baseline <- function(
    T_weeks,        # number of weeks to simulate
    A,              # number of age groups
    N,              # vector length A: population by age group
    R_init_prop,    # vector length A: proportion already immune by age
    I0,             # vector length A: initial infections by age
    base_beta,      # vector length T_weeks: weekly transmission rate
    sigma,          # 1 / latent period (per week)
    gamma,          # recovery rate (per week)
    rho,            # reporting rate
    prop_symp = 0.5242478,
    sub_steps = 7
) {
  pmax0 <- function(x) pmax(0, x)
  N_total <- sum(N)
  dt      <- 1 / sub_steps   # length of each sub-step in "weeks"
  
  # Storage at WEEKLY resolution
  S <- E <- I <- R <- matrix(0, nrow = A, ncol = T_weeks)
  new_infections <- new_symptomatic <- matrix(0, nrow = A, ncol = T_weeks)
  
  # Working state at fine resolution
  S_now <- N - I0 - R_init_prop * N
  E_now <- rep(0, A)
  I_now <- I0
  R_now <- R_init_prop * N
  
  # Record t = 1
  S[, 1] <- S_now;  E[, 1] <- E_now;  I[, 1] <- I_now;  R[, 1] <- R_now
  
  for (t in 2:T_weeks) {
    
    # Track new infections accumulated within this week
    new_I_week <- rep(0, A)
    
    # beta is constant across the sub-steps within a week
    beta_t <- base_beta[t - 1]
    
    for (k in 1:sub_steps) {
      foi <- beta_t * sum(I_now) / N_total
      
      new_E <- foi   * S_now * dt
      new_I <- sigma * E_now * dt
      new_R <- gamma * I_now * dt
      
      S_now <- pmax0(S_now - new_E)
      E_now <- pmax0(E_now + new_E - new_I)
      I_now <- pmax0(I_now + new_I - new_R)
      R_now <- pmax0(R_now + new_R)
      
      new_I_week <- new_I_week + new_I
    }
    
    S[, t] <- S_now;  E[, t] <- E_now;  I[, t] <- I_now;  R[, t] <- R_now
    new_infections[, t]  <- new_I_week
    new_symptomatic[, t] <- prop_symp * new_I_week
  }
  
  list(
    S = S, E = E, I = I, R = R,
    new_infections  = new_infections,
    new_symptomatic = new_symptomatic,
    new_reported    = rho * new_symptomatic   # reported = rho * symptomatic
  )
}

# ------------------------------------------------------------
# 2. Spline helper: turn a vector of spline coefs into beta_t
# ------------------------------------------------------------
# We use a natural cubic spline (ns) with knots evenly spaced across weeks.
# Number of degrees of freedom (df) = number of spline coefficients.
# More df = more flexible curve; fewer = smoother.
# Try df = 6 to start; can experiment with 5, 7, 8.

make_beta_t <- function(coefs, T_weeks, df) {
  # Build spline basis: matrix with T_weeks rows, df columns
  weeks <- 1:T_weeks
  basis <- ns(weeks, df = df, intercept = TRUE)
  # Linear combination of basis functions, exponentiated to ensure beta > 0
  log_beta <- as.numeric(basis %*% coefs)
  exp(log_beta)
}

# ------------------------------------------------------------
# 3. Load Fortaleza-specific data (age structure, observed cases)
# ------------------------------------------------------------
municipality_name <- "fortaleza"
state_name        <- "Ceará"

normalize <- function(x) tolower(gsub("[áàâã]", "a",
                                      gsub("[éèê]",   "e",
                                           gsub("[íì]",    "i",
                                                gsub("[óòôõ]", "o",
                                                     gsub("[úù]",   "u",
                                                          gsub("ç",      "c", x)))))))

age_df <- read_excel("population.xlsx", sheet = "prop_immune")
age_df <- age_df[normalize(age_df$municipality) == normalize(municipality_name), ]
stopifnot(nrow(age_df) > 0)
age_df$age_group <- factor(age_df$age_group, levels = age_df$age_group)

N           <- age_df$pop_num
R_init_prop <- age_df$prop_ever_infected
A           <- nrow(age_df)

# Observed weekly cases
raw <- read_excel("weekly_case.xlsx", sheet = "weekly_all")
fortaleza_obs <- raw |>
  filter(Code == 230440 & Year == 2022) |>
  pivot_longer(cols = starts_with("Week"),
               names_to = "week", values_to = "cases") |>
  mutate(week  = as.integer(sub("Week ", "", week)),
         cases = ifelse(is.na(cases), 0, cases)) |>
  arrange(week)

observed_cases <- fortaleza_obs$cases
T_weeks <- length(observed_cases)
stopifnot(T_weeks == 52, sum(observed_cases) == 29660)

# Fixed parameters
sigma     <- 1 / 0.60
gamma     <- 0.54
prop_symp <- 0.5242478

# Initial infections (back-calculated, frozen)
week_1_cases <- 10
rho_for_seed <- 0.10
I0_total <- week_1_cases / rho_for_seed / prop_symp
susceptible_pop <- N * (1 - R_init_prop)
I0 <- round(I0_total * susceptible_pop / sum(susceptible_pop))

cat("Loaded ", T_weeks, " weeks of observed data\n", sep = "")
cat("Observed total cases:", sum(observed_cases), "\n")
cat("I0 (seeded):", sum(I0), "\n")

# ------------------------------------------------------------
# 4. Negative log-likelihood function for optim()
# ------------------------------------------------------------
# Parameters to estimate:
#   - spline_coefs: df-length vector (controls shape of beta_t)
#   - log_rho:     log of reporting rate (so optim works on unconstrained scale)
#   - log_theta:   log of overdispersion (likewise)
#
# We work on log-scale for rho and theta to keep them positive and let optim
# search freely without bounds.

df_spline <- 6   # number of spline coefficients; experiment with 5-8

neg_log_lik <- function(params, 
                        observed, 
                        T_weeks, 
                        df, 
                        A, 
                        N, 
                        R_init_prop,
                        I0, 
                        sigma, 
                        gamma, 
                        prop_symp) {
  
  # Unpack parameter vector
  spline_coefs <- params[1:df]
  log_rho      <- params[df + 1]
  log_theta    <- params[df + 2]
  
  rho   <- exp(log_rho)
  theta <- exp(log_theta)
  
  # Build beta_t from spline coefs
  beta_t <- make_beta_t(spline_coefs, T_weeks, df)
  
  # Reject pathological beta values
  if (any(!is.finite(beta_t)) || any(beta_t > 20) || any(beta_t < 1e-6)) {
    return(1e10)  # huge penalty
  }
  
  # Run SEIR
  out <- tryCatch(
    seir_baseline(
      T_weeks = T_weeks, 
      A = A, 
      N = N, 
      R_init_prop = R_init_prop,
      I0 = I0, 
      base_beta = beta_t,
      sigma = sigma, 
      gamma = gamma, 
      rho = rho, 
      prop_symp = prop_symp
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(1e10)
  
  predicted <- colSums(out$new_reported)
  
  # Guard against zeros (negbin requires mean > 0)
  predicted <- pmax(predicted, 1e-6)
  
  # Negative binomial log-likelihood
  # dnbinom in R: mean = mu, size = theta
  ll <- sum(dnbinom(observed, mu = predicted, size = theta, log = TRUE))
  
  if (!is.finite(ll)) return(1e10)
  -ll
}

# ------------------------------------------------------------
# 5. Run optim()
# ------------------------------------------------------------
# Starting values: fit spline to log of Hyolim's Ceará beta as initial guess
beta_all <- read_excel("beta_hyolim.xlsx", sheet = "beta")
beta_state <- beta_all[normalize(beta_all$region) == normalize(state_name), ]
beta_state$week <- as.integer(sub("beta_week", "", beta_state$parameter))
beta_state <- beta_state[order(beta_state$week), ]
hyolim_beta <- beta_state$median

# Initial spline coefs: regress log(hyolim_beta) on the spline basis
weeks <- 1:T_weeks
basis <- ns(weeks, df = df_spline, intercept = TRUE)
init_spline_coefs <- coef(lm(log(hyolim_beta) ~ basis - 1))

init_params <- c(
  init_spline_coefs,
  log(0.25),   # log_rho starting at rho = 0.25
  log(10)      # log_theta starting at theta = 10
)

cat("\nStarting optim()...\n")
fit <- optim(
  par = init_params,
  fn  = neg_log_lik,
  observed = observed_cases,
  T_weeks  = T_weeks, 
  df = df_spline,
  A = A, 
  N = N, 
  R_init_prop = R_init_prop,
  I0 = I0, 
  sigma = sigma, 
  gamma = gamma, 
  prop_symp = prop_symp,
  method  = "BFGS",
  hessian = TRUE,
  control = list(maxit = 1000, trace = 1, REPORT = 10)
)

# ------------------------------------------------------------
# 6. Extract results
# ------------------------------------------------------------
best_coefs <- fit$par[1:df_spline]
best_rho   <- exp(fit$par[df_spline + 1])
best_theta <- exp(fit$par[df_spline + 2])

best_beta_t <- make_beta_t(best_coefs, T_weeks, df_spline)

cat("\n--- Fit results ---\n")
cat("Convergence code (0 = success):", fit$convergence, "\n")
cat("Final negative log-likelihood:", fit$value, "\n")
cat("Best rho:  ", round(best_rho, 4), "\n")
cat("Best theta:", round(best_theta, 2), "\n")
cat("beta_t range: [", round(min(best_beta_t), 3), ", ",
    round(max(best_beta_t), 3), "]\n", sep = "")

# Run SEIR with best-fit parameters for plotting
out_best <- seir_baseline(
  T_weeks = T_weeks, 
  A = A, 
  N = N, 
  R_init_prop = R_init_prop,
  I0 = I0, 
  base_beta = best_beta_t,
  sigma = sigma, 
  gamma = gamma, 
  rho = best_rho
)

predicted_cases <- colSums(out_best$new_reported)

cat("Predicted total reported cases:", round(sum(predicted_cases)), "\n")
cat("Observed total reported cases: ", sum(observed_cases), "\n")
cat("Predicted peak week:", which.max(predicted_cases),
    "(", round(max(predicted_cases)), "cases)\n")
cat("Observed peak week: ", which.max(observed_cases),
    "(", max(observed_cases), "cases)\n")

# ------------------------------------------------------------
# 7. Plots
# ------------------------------------------------------------
# (a) Fitted beta_t vs Hyolim's Ceará beta_t
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(weeks, hyolim_beta, type = "l", lwd = 2, col = "grey50", lty = 2,
     xlab = "Week", ylab = expression(beta[t]),
     main = "Fitted Fortaleza beta_t vs Hyolim's Ceará beta_t",
     ylim = range(c(hyolim_beta, best_beta_t)))
lines(weeks, best_beta_t, lwd = 2, col = "#d6604d")
legend("topright", legend = c("Ceará (Hyolim)", "Fortaleza (MLE)"),
       col = c("grey50", "#d6604d"), lty = c(2, 1), lwd = 2, bty = "n")

# (b) Predicted vs observed reported cases
plot(weeks, observed_cases, type = "p", pch = 16, cex = 0.6, col = "black",
     xlab = "Week", ylab = "Reported CHIKV cases",
     main = "Fortaleza 2022: predicted vs observed reported cases")
lines(weeks, predicted_cases, lwd = 2, col = "#d6604d")
legend("topright", legend = c("Observed (SINAN)", "Predicted (MLE fit)"),
       col = c("black", "#d6604d"), lty = c(NA, 1), pch = c(16, NA),
       lwd = c(NA, 2), bty = "n")
par(mfrow = c(1, 1))

# ------------------------------------------------------------
# 8. Approximate confidence intervals from the Hessian
# ------------------------------------------------------------
# Standard errors come from inv(Hessian) under the asymptotic MLE result.
# Only meaningful if optim converged cleanly and Hessian is positive-definite.
if (fit$convergence == 0) {
  vcov_mat <- tryCatch(solve(fit$hessian), error = function(e) NULL)
  if (!is.null(vcov_mat)) {
    se <- sqrt(diag(vcov_mat))
    cat("\n--- Approximate 95% CIs ---\n")
    cat("rho:   [",
        round(exp(fit$par[df_spline + 1] - 1.96 * se[df_spline + 1]), 4), ", ",
        round(exp(fit$par[df_spline + 1] + 1.96 * se[df_spline + 1]), 4), "]\n",
        sep = "")
    cat("theta: [",
        round(exp(fit$par[df_spline + 2] - 1.96 * se[df_spline + 2]), 2), ", ",
        round(exp(fit$par[df_spline + 2] + 1.96 * se[df_spline + 2]), 2), "]\n",
        sep = "")
  } else {
    cat("\nHessian could not be inverted; no CIs reported.\n")
  }
}

# ------------------------------------------------------------
# 9. Uncertainty bands via parametric bootstrap from the MLE
# ------------------------------------------------------------
library(MASS)  # for mvrnorm

set.seed(123)  # reproducibility
n_samples <- 500   # number of bootstrap samples (more = smoother bands, slower)

# Sample from multivariate normal centred at MLE
param_samples <- mvrnorm(n = n_samples,
                         mu = fit$par,
                         Sigma = vcov_mat)

# For each sample, run the SEIR and store predicted weekly cases + beta_t
pred_matrix <- matrix(NA, nrow = n_samples, ncol = T_weeks)
beta_matrix <- matrix(NA, nrow = n_samples, ncol = T_weeks)

cat("\nRunning", n_samples, "bootstrap samples...\n")
for (i in 1:n_samples) {
  spline_coefs_i <- param_samples[i, 1:df_spline]
  rho_i          <- exp(param_samples[i, df_spline + 1])
  beta_t_i       <- make_beta_t(spline_coefs_i, T_weeks, df_spline)
  
  # Skip if beta_t is pathological
  if (any(!is.finite(beta_t_i)) || any(beta_t_i > 20) || any(beta_t_i < 1e-6)) next
  
  out_i <- tryCatch(
    seir_baseline(
      T_weeks = T_weeks, 
      A = A, 
      N = N, 
      R_init_prop = R_init_prop,
      I0 = I0, 
      base_beta = beta_t_i,
      sigma = sigma, 
      gamma = gamma, 
      rho = rho_i
    ),
    error = function(e) NULL
  )
  if (is.null(out_i)) next
  
  pred_matrix[i, ] <- colSums(out_i$new_reported)
  beta_matrix[i, ] <- beta_t_i
  if (i %% 50 == 0) cat("  ", i, "/", n_samples, "\n")
}

# Compute pointwise quantiles (2.5%, 50%, 97.5%) at each week
pred_lower  <- apply(pred_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
pred_median <- apply(pred_matrix, 2, quantile, probs = 0.500, na.rm = TRUE)
pred_upper  <- apply(pred_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)

beta_lower  <- apply(beta_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
beta_median <- apply(beta_matrix, 2, quantile, probs = 0.500, na.rm = TRUE)
beta_upper  <- apply(beta_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)

# ------------------------------------------------------------
# 10. Plot with uncertainty bands (ggplot2 version)
# ------------------------------------------------------------
plot_df <- data.frame(
  week     = weeks,
  observed = observed_cases,
  pred_med = pred_median,
  pred_lo  = pred_lower,
  pred_hi  = pred_upper
)

ggplot(plot_df, aes(x = week)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = "Predicted (95% band)"),
              alpha = 0.25) +
  geom_line(aes(y = pred_med, colour = "Predicted (MLE fit)"), linewidth = 1) +
  geom_point(aes(y = observed, colour = "Observed (SINAN)"), size = 1.5) +
  scale_colour_manual(name = NULL,
                      values = c("Observed (SINAN)"  = "black",
                                 "Predicted (MLE fit)" = "#f4a582")) +
  scale_fill_manual(name = NULL,
                    values = c("Predicted (95% band)" = "#f4a582")) +
  scale_x_continuous(breaks = seq(0, 52, 10)) +
  labs(x = "Week",
       y = "Reported CHIKV cases",
       title = "Fortaleza 2022: predicted vs observed reported cases") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 11),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# Same idea for beta_t
beta_df <- data.frame(
  week     = weeks,
  beta_med = beta_median,
  beta_lo  = beta_lower,
  beta_hi  = beta_upper
)

ggplot(beta_df, aes(x = week)) +
  geom_ribbon(aes(ymin = beta_lo, ymax = beta_hi),
              fill = "#a8d1e7", alpha = 0.25) +
  geom_line(aes(y = beta_med), colour = "#a8d1e7", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 52, 10)) +
  labs(x = "Epidemiological week (2022)",
       y = expression(beta[t]),
       title = "Fitted weekly transmission rate (Fortaleza) with 95% band") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", colour = "grey80"))