library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(splines)
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# Shared helpers incl. the canonical Caldas case loader load_caldas_age_cases().
if (!exists("load_caldas_age_cases")) source("ca_common.R")

# ------------------------------------------------------------
# 1. SEIR simulator: weekly time-stepped, age-structured.
#    Takes a weekly transmission-rate vector (base_beta) and returns weekly
#    compartments + incidence. Seeds E0/I0 so it can start mid-epidemic.
# ------------------------------------------------------------
seir_baseline <- function(
    T_weeks,        # number of weeks to simulate
    A,              # number of age groups
    N,              # vector length A: population by age group
    R_init_prop,    # vector length A: proportion already immune by age
    I0,             # vector length A: initial infectious by age (a stock)
    base_beta,      # vector length T_weeks: weekly transmission rate
    sigma,          # 1 / latent period (per week)
    gamma,          # recovery rate (per week)
    rho,            # reporting rate
    prop_symp = 0.5242478,
    sub_steps = 7,
    E0 = rep(0, A)  # initial exposed by age; >0 seeds an already-running epidemic
) {
  pmax0 <- function(x) pmax(0, x)
  N_total <- sum(N)
  dt      <- 1 / sub_steps   # length of each sub-step in "weeks"
  
  # Storage at weekly resolution
  S <- E <- I <- R <- matrix(0, nrow = A, ncol = T_weeks)
  new_infections <- new_symptomatic <- matrix(0, nrow = A, ncol = T_weeks)
  
  # Working state at the start of week 1. E0 may be > 0 because the outbreak was
  # ongoing: reported cases are driven by the E->I flow (sigma * E),
  # so a non-zero E0 makes week 1 incidence start at the
  # observed level instead of spinning up from zero.
  S_now <- pmax0(N - I0 - E0 - R_init_prop * N)
  E_now <- E0
  I_now <- I0
  R_now <- R_init_prop * N
  
  # Loop from week 1 so incidence is recorded for every week (including week 1).
  for (t in 1:T_weeks) {
    
    # Track new infections accumulated within this week
    new_I_week <- rep(0, A)
    
    # beta is constant across the sub-steps within a week
    beta_t <- base_beta[t]
    
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
    
    S[, t] <- S_now;  
    E[, t] <- E_now;  
    I[, t] <- I_now;  
    R[, t] <- R_now
    new_infections[, t]  <- new_I_week
    new_symptomatic[, t] <- prop_symp * new_I_week
  }
  
  list(
    S = S, 
    E = E, 
    I = I, 
    R = R,
    new_infections  = new_infections,
    new_symptomatic = new_symptomatic,
    new_reported    = rho * new_symptomatic   # reported = rho * symptomatic
  )
}

# ------------------------------------------------------------
# 2. Spline helper: spline coefficients -> weekly beta_t
# ------------------------------------------------------------
# beta_t = exp(natural-spline basis %*% coefs), so beta stays positive.
# The spline is fitted over the ACTIVE window (weeks 1 to active_weeks); beyond it
# beta is HELD CONSTANT at its last active value. This makes the data-poor tail flat
# by construction (no spline-boundary curl-up). If active_weeks == T_weeks this
# reduces to a plain full-window spline.
# Uses globals basis_full (active_weeks x df), active_weeks, T_weeks.
make_beta_t <- function(coefs) {
  beta_active <- as.numeric(exp(basis_full %*% coefs))         # length active_weeks, > 0
  beta_full   <- rep(beta_active[active_weeks], T_weeks)       # flat-held tail
  beta_full[seq_len(active_weeks)] <- beta_active
  beta_full
}

# ------------------------------------------------------------
# 3. Load specific data (age structure, observed cases)
# ------------------------------------------------------------
municipality_name <- "Caldas Novas"
state_name        <- "Goiás"

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

# Population by age group is recorded for 2022 in the spreadsheet. The Caldas Novas
# outbreak we model runs over 2025-2026, so we grow the 2022 age-group counts to
# 2025 using the town-level exponential growth rate:
#   r = ln(pop_2025 / pop_2022) / (years between)  [per year]
#   pop_2025 = pop_2022 * exp(r * 3)               [3 years on from 2022]
# The growth rate is applied uniformly to every age group (counts scale, the
# ever-infected proportions are unchanged).
pop_2022_total <- 98622
pop_2025_total <- 106820
growth_r       <- log(pop_2025_total / pop_2022_total) / 3   # per-year growth rate
N           <- age_df$pop_num * exp(growth_r * 3)            # 2022 -> 2025 counts
R_init_prop <- age_df$prop_ever_infected
A           <- nrow(age_df)
cat(sprintf("Population growth rate r = %.5f /yr; total pop 2022 = %d -> 2025 = %d\n",
            growth_r, round(sum(age_df$pop_num)), round(sum(N))))

# Overall proportion of the population already immune at the outbreak start
# (population-weighted average of the per-age-group ever-infected proportions).
pct_immune <- 100 * sum(R_init_prop * N) / sum(N)
cat(sprintf("Total population already immune at start: %.1f%%\n", pct_immune))

# Observed weekly cases: age-stratified SINAN download ("ca_combined" sheet), loaded
# via the shared loader in ca_common.R (single source of truth, also used by
# weekly_age_stratified.R). Window 2025-W23 -> 2026-W22 = 53 weeks (2025 has an epi
# Semana 53), one more than the older plain weekly_all series (weekly_case.R).
ca_cases       <- load_caldas_age_cases()
caldas_obs     <- ca_cases$caldas_obs
observed_cases <- ca_cases$observed_cases
ca_age         <- ca_cases$ca_age          # age x week (for the age-stratified burden)
age_totals     <- ca_cases$age_totals      # observed cases by age group
obs_band_prop  <- ca_cases$obs_band_prop   # observed case share across the 9 CFR bands
T_weeks <- length(observed_cases)
stopifnot(T_weeks == 53, sum(observed_cases) == 8209)

# Fixed parameters
sigma     <- 1 / 0.60
gamma     <- 0.54
prop_symp <- 0.5242478
# Reporting rate is only weakly identified (the NLL is nearly flat in rho because
# beta rescales to compensate), so we fix it on epidemiological grounds rather than
# estimate it. rho = 0.45 implies ~45% attack rate among susceptibles, typical for
# a first CHIKV wave. The fit is insensitive to this within the feasible range
# (for Caldas Novas any rho from ~0.25 up is feasible; lower rho => higher attack rate).
rho_fixed <- 0.25

# Initial conditions (back-calculated from observed week-1 reported cases, frozen).
# We capture the outbreak FROM ITS START (2025-W23, a low baseline before take-off),
# so the epidemic is NOT already running: we seed a small infectious stock I0 from the
# week-1 reported count and start with no exposed (E0 = 0), letting the SEIR spin the
# epidemic up naturally.
#   week-1 incidence  inc1 = reported / (rho * prop_symp)
#   I0 = inc1 / gamma   (infectious stock consistent with that incidence)
week_1_cases <- observed_cases[1]
inc1     <- week_1_cases / rho_fixed / prop_symp
I0_total <- inc1 / gamma
susceptible_pop <- N * (1 - R_init_prop)
I0 <- round(I0_total * susceptible_pop / sum(susceptible_pop))
E0 <- rep(0, A)

cat("Loaded ", T_weeks, " weeks of observed data\n", sep = "")
cat("Observed total cases:", sum(observed_cases), "\n")
cat("I0 (seeded infectious):", sum(I0), " E0 (seeded exposed):", sum(E0), "\n")

# ------------------------------------------------------------
# 4. Objective for optim(): negative log-posterior
# ------------------------------------------------------------
# Estimated parameters (rho is fixed, not estimated):
#   - spline_coefs : df-length vector controlling the shape of beta_t
#   - log_theta    : log(theta - 1), kept on an unconstrained scale so theta > 1
# Objective = negative-binomial log-likelihood + Lognormal log-prior on beta_t.

# Spline df is chosen by a sweep below (section 5). Set df_choice to a number (e.g. 6)
# to FORCE a particular df; leave it NA to auto-pick the lowest-BIC model.
# NOTE: lowest BIC is only a default -- always sanity-check the chosen model for
# MEANINGFULNESS (plausible R0 range and attack rate) before trusting it.
df_choice <- 5

# Active window: fit a time-varying beta only up to this week_index, then hold beta
# constant (flat) afterwards, which prevents a spline curl-up where beta is not identified.
# For Caldas Novas the outbreak is in the BACK half of the window, but from ~2026-W19 the
# susceptible pool is depleted (<20%): incidence there is insensitive to beta, so a
# full-window spline curls UPWARD -- a spurious late rise (R0 -> 3.24 at the very last
# week) as the model raises beta to squeeze the declining tail out of a depleted pool.
# Holding beta flat from index 50 (2026-W19) puts R0_max back at the true 2025-W50 peak
# (~3.1) and matches the observed total. (Set active_weeks <- T_weeks for a full-window spline.)
#   2025-W23 = index 1 ; 2026-W09 (case peak) = index 40 ; 2026-W19 = index 50 ; 2026-W22 = index 53.
active_weeks <- 50

#     Lognormal prior (regularisation) on weekly beta_t, in the spirit of Hyolim's
#     beta_t ~ Lognormal. We add it as a MAP penalty so the fit stays MLE/optim-based.
#     log(beta_t) ~ Normal(prior_logmean, prior_logsd).
#     SCALE NOTE: this is our homogeneous-mixing beta (= R0 * gamma), NOT Hyolim's
#     contact-matrix-scaled beta. Her prior median exp(-1)=0.37 is on her scale; on
#     ours the epidemiologically-neutral baseline is beta ~ gamma = 0.54 (R0 ~ 1), so
#     we centre the prior there. In the data-rich peak the likelihood dominates and
#     beta rises well above this; in the data-poor tail beta reverts toward ~0.54
#     (a low, flat off-season level) instead of curling up.
#     Smaller prior_logsd = stronger regularisation. Set prior_logsd = Inf to disable.
prior_logmean <- log(0.54)   # prior median beta = 0.54 (R0 ~ 1 on our scale)
prior_logsd   <- 0.70        # loosened from 0.40: at rho = 0.25 the outbreak implies a
# ~76% attack rate (R0 up to ~3.3), which the tighter 0.40 prior
# suppressed (the fit undershot the total by ~5%). Caldas has data
# across the WHOLE window (active_weeks = T_weeks), so heavy beta
# regularisation is not needed here; 0.70 keeps R0_max ~3.2 and the
# predicted total within ~1% of observed at both rho = 0.25 and 0.45.

# Peak-emphasis weighting (to track the trend / peak, as in Hyolim's fitted curves
# which sit slightly ABOVE the observed points at the peak). A negative-binomial fit
# alone tends to UNDERSHOOT the peak: at high counts the NB variance is large, so
# missing the peak is "cheap" and the smooth spline saves effort by sitting under it.
# We up-weight each week's log-likelihood by w_t = 1 + peak_emphasis * (obs_t / max obs),
# so the high-incidence weeks dominate and the curve is pulled up to match (and lightly
# overshoot) the peak. peak_emphasis = 0 recovers the plain (unweighted) likelihood.
peak_emphasis <- 10

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
                        prop_symp, 
                        rho, 
                        E0) {
  
  spline_coefs <- params[1:df]
  log_theta    <- params[df + 1]
  
  theta <- 1 + exp(log_theta)   # floored at 1, can't collapse to ~0
  
  beta_t <- make_beta_t(spline_coefs)
  
  if (any(!is.finite(beta_t)) || any(beta_t > 20) || any(beta_t < 1e-6)) {
    return(1e10)
  }
  
  out <- tryCatch(
    seir_baseline(T_weeks = T_weeks, 
                  A = A, 
                  N = N, 
                  R_init_prop = R_init_prop,
                  I0 = I0, 
                  base_beta = beta_t,
                  sigma = sigma, 
                  gamma = gamma, 
                  rho = rho, 
                  prop_symp = prop_symp,
                  E0 = E0),
    error = function(e) NULL
  )
  if (is.null(out)) return(1e10)
  
  predicted <- pmax(colSums(out$new_reported), 1e-6)
  
  # Peak-emphasis weights: up-weight high-incidence weeks so the fit tracks the peak
  # (Hyolim-style overshoot) instead of the NB-driven undershoot. peak_emphasis = 0
  # gives the standard unweighted likelihood.
  wts <- 1 + peak_emphasis * (observed / max(observed))
  ll  <- sum(wts * dnbinom(observed, mu = predicted, size = theta, log = TRUE))
  if (!is.finite(ll)) return(1e10)
  
  # MAP penalty: Lognormal prior on weekly beta_t, applied over all T_weeks.
  # It mainly regularises the long low-incidence tail (2024-W20 onward), where the
  # data barely inform beta, pulling it toward the prior median.
  # prior_logsd = Inf disables the prior (recovers the plain MLE).
  log_prior <- if (is.finite(prior_logsd)) {
    sum(dnorm(log(beta_t), mean = prior_logmean, sd = prior_logsd, log = TRUE))
  } else 0
  
  -(ll + log_prior)   # negative log-posterior; optim now finds the penalised (MAP) estimate
}

# ------------------------------------------------------------
# 5. Fit the model, sweeping over spline df to compare models
# ------------------------------------------------------------
weeks    <- 1:T_weeks
peak_idx <- which.max(observed_cases)              # case peak (index 39, 2026-W09)

# Spline basis is built over the ACTIVE window only (1:active_weeks); make_beta_t holds
# beta flat beyond it. Helper keeps fit_for_df and diagnose_fit consistent.
build_basis <- function(df) ns(seq_len(active_weeks), df = df, intercept = TRUE)

# Fit the model for a given spline df, using 3 starting beta shapes (multistart) and
# keeping the best. Sets the GLOBAL basis_full so make_beta_t() uses the matching
# basis. No Hyolim reference is used: the case curve rises from the 2025-W23 seed to
# a peak then decays, so a beta that starts moderate, rises, then falls is reasonable.
fit_for_df <- function(df) {
  basis_full <<- build_basis(df)   # global, matched to this df, over active window
  coefs_for_logbeta <- function(log_beta_vec) coef(lm(log_beta_vec ~ basis_full - 1))
  start_betas <- list(   # length active_weeks (peak_idx falls inside the active window)
    flat      = rep(1.2, active_weeks),
    humped    = c(seq(1.2, 2.2, length.out = peak_idx),
                  seq(2.2, 0.5, length.out = active_weeks - peak_idx)),
    declining = seq(2.2, 0.5, length.out = active_weeks)
  )
  best <- NULL
  for (nm in names(start_betas)) {
    init_params <- c(coefs_for_logbeta(log(start_betas[[nm]])), log(50))
    f <- tryCatch(
      optim(par = init_params, fn = neg_log_lik,
            observed = observed_cases, T_weeks = T_weeks, df = df,
            A = A, N = N, R_init_prop = R_init_prop, I0 = I0,
            sigma = sigma, gamma = gamma, prop_symp = prop_symp,
            rho = rho_fixed, E0 = E0,
            method = "BFGS", hessian = TRUE, control = list(maxit = 1000)),
      error = function(e) NULL
    )
    if (!is.null(f) && (is.null(best) || f$value < best$value)) best <- f
  }
  best
}

# Diagnostics for a fitted model: likelihood-only logLik, AIC, BIC, the implied
# R0 = beta/gamma range, predicted total, peak week, and the attack rate. These let
# us judge MEANINGFULNESS (sensible R0 & attack rate), not just the BIC score.
diagnose_fit <- function(f, df) {
  basis_full <<- build_basis(df)
  coefs  <- f$par[1:df]
  theta  <- 1 + exp(f$par[df + 1])
  beta_t <- make_beta_t(coefs)
  out    <- seir_baseline(T_weeks = T_weeks, A = A, N = N, R_init_prop = R_init_prop,
                          I0 = I0, base_beta = beta_t, sigma = sigma, gamma = gamma,
                          rho = rho_fixed, E0 = E0)
  pred   <- pmax(colSums(out$new_reported), 1e-6)
  ll     <- sum(dnbinom(observed_cases, mu = pred, size = theta, log = TRUE))
  k      <- df + 1                                 # spline coefs + log_theta (rho fixed)
  R0_t   <- beta_t / gamma                         # basic reproduction number over time
  immune_end <- 100 * (1 - sum(out$S[, T_weeks]) / sum(N))             # ever infected by end
  attack     <- 100 * sum(out$new_infections) / sum(N * (1 - R_init_prop)) # of initially susceptible
  data.frame(df = df, k = k, logLik = ll,
             AIC = -2 * ll + 2 * k, BIC = -2 * ll + k * log(T_weeks),
             R0_min = min(R0_t), R0_max = max(R0_t),
             pred_total = round(sum(colSums(out$new_reported))),
             peak_wk = which.max(pred),
             immune_end_pct = immune_end, attack_pct = attack)
}

# # ---- Sweep over candidate df values and compare ----
# df_grid <- 5:8
# fits_by_df <- lapply(df_grid, fit_for_df)
# names(fits_by_df) <- df_grid
# sweep <- do.call(rbind, Map(function(f, d) if (!is.null(f)) diagnose_fit(f, d),
#                             fits_by_df, df_grid))
# cat("\n--- Spline df comparison (lower BIC = better fit/complexity tradeoff) ---\n")
# cat("    Sanity-check R0 (~2-4 plausible for a first CHIKV wave) and attack rate\n")
# cat("    BEFORE trusting the min-BIC pick.\n")
# print(sweep, row.names = FALSE, digits = 4)
# 
# # ---- Choose the model to carry forward ----
df_spline <- if (is.na(df_choice)) sweep$df[which.min(sweep$BIC)] else df_choice
cat(sprintf("\nChosen spline df: %d (%s)\n", df_spline,
            if (is.na(df_choice)) "auto: lowest BIC" else "manually set via df_choice"))
fit <- fit_for_df(df_spline)   # also resets basis_full to match df_spline

# ------------------------------------------------------------
# 6. Extract results
# ------------------------------------------------------------
best_coefs  <- fit$par[1:df_spline]
best_rho    <- rho_fixed
best_theta  <- 1 + exp(fit$par[df_spline + 1])
best_beta_t <- make_beta_t(best_coefs)

# ---- Export the fitted weekly transmission rate beta_t for downstream models ----
# (e.g. MAYV_ca_pre_vacc.R: replaces the hard-coded season_mean1). The MAYV window is
# aligned to this one (2025-W23 -> 2026-W22), so index t maps 1:1; MAYV uses a
# mean-normalised seasonal SHAPE (it rescales by R0 * gamma), so we also export that.
beta_season_shape <- best_beta_t / mean(best_beta_t)   # mean = 1: drop-in for season_mean1
saveRDS(best_beta_t,       "caldas_beta_fitted.rds")   # absolute fitted beta_t (length T_weeks)
saveRDS(beta_season_shape, "caldas_beta_season.rds")   # same curve, normalised to mean 1
write.csv(data.frame(week_index   = weeks,
                     week_label    = caldas_obs$week_label,
                     beta_t        = best_beta_t,
                     beta_shape    = beta_season_shape,   # mean-1 normalised
                     R0            = best_beta_t / gamma),
          "caldas_beta_fitted.csv", row.names = FALSE)
cat("Exported fitted beta_t to caldas_beta_fitted.{rds,csv} and caldas_beta_season.rds\n")

cat("\n--- Fit results ---\n")
cat("Convergence code (0 = success):", fit$convergence, "\n")
cat("Final objective (neg log-posterior):", fit$value, "\n")
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
  rho = best_rho,
  E0 = E0
)

predicted_cases <- colSums(out_best$new_reported)

cat("Predicted total reported cases:", round(sum(predicted_cases)), "\n")
cat("Observed total reported cases: ", sum(observed_cases), "\n")
cat("Predicted peak week:", which.max(predicted_cases),
    "(", round(max(predicted_cases)), "cases)\n")
cat("Observed peak week: ", which.max(observed_cases),
    "(", max(observed_cases), "cases)\n")

# ------------------------------------------------------------
# 6b. Goodness of fit: log-likelihood, AIC, BIC
# ------------------------------------------------------------
# BIC = -2 * logLik + k * ln(n)
#   logLik : the maximised LOG-LIKELIHOOD only (the negbin data fit) -- NOT the
#            log-posterior. We strip the Lognormal prior penalty back out so this is
#            a likelihood-based score comparable across models.
#   k      : number of FREELY ESTIMATED parameters = df_spline spline coefs + 1
#            (log_theta). rho is fixed, so it does not count.
#   n      : number of data points = T_weeks (52 weekly case counts).
# AIC = -2 * logLik + 2 * k is also reported. Lower BIC/AIC = better, but only as a
# RELATIVE score between competing models on the SAME data (e.g. df = 5 vs 6 vs 7).
loglik_best <- sum(dnbinom(observed_cases, mu = pmax(predicted_cases, 1e-6),
                           size = best_theta, log = TRUE))
k_params <- df_spline + 1     # spline coefs + log_theta (rho fixed)
n_obs    <- T_weeks
AIC_val  <- -2 * loglik_best + 2 * k_params
BIC_val  <- -2 * loglik_best + k_params * log(n_obs)

cat("\n--- Goodness of fit ---\n")
cat(sprintf("Log-likelihood: %.2f\n", loglik_best))
cat(sprintf("Parameters (k): %d   Observations (n): %d\n", k_params, n_obs))
cat(sprintf("AIC: %.2f\n", AIC_val))
cat(sprintf("BIC: %.2f\n", BIC_val))

# Meaningfulness checks for the chosen model: implied R0 = beta/gamma over time, and
# the attack rate. These must be epidemiologically plausible, not just BIC-optimal.
R0_t          <- best_beta_t / gamma
immune_end    <- 100 * (1 - sum(out_best$S[, T_weeks]) / sum(N))   # ever infected by end
attack_suscep <- 100 * sum(out_best$new_infections) / sum(susceptible_pop)
cat("\n--- Meaningfulness checks ---\n")
cat(sprintf("Implied R0 = beta/gamma range: [%.2f, %.2f]\n", min(R0_t), max(R0_t)))
cat(sprintf("Immune at start: %.1f%%  ->  ever infected by end: %.1f%%\n",
            pct_immune, immune_end))
cat(sprintf("Outbreak attack rate (of initially susceptible): %.1f%%\n", attack_suscep))

# ------------------------------------------------------------
# 7. Plots
# ------------------------------------------------------------
# Shared x-axis ticks for the 2025-W23 -> 2026-W22 window: week 30, 40, 50 (2025),
# the 2026 boundary, then week 10, 20 (2026).
x_ticks <- caldas_obs |>
  filter((Year == 2025 & week %in% c(30, 40, 50)) |
           (Year == 2026 & week %in% c(10, 20)))
year_break <- mean(c(max(caldas_obs$week_index[caldas_obs$Year == 2025]),
                     min(caldas_obs$week_index[caldas_obs$Year == 2026])))

# (a) Fitted beta_t (base-R quick look)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(weeks, best_beta_t, type = "l", lwd = 2, col = "#d6604d",
     xlab = "Week", ylab = expression(beta[t]),
     main = "Fitted Caldas Novas beta_t (2025-W23 to 2026-W22)",
     xaxt = "n")
axis(1, at = x_ticks$week_index, labels = x_ticks$week)
abline(v = year_break, lty = 2, col = "grey50")

# (b) Predicted vs observed reported cases
plot(weeks, observed_cases, type = "p", pch = 16, cex = 0.6, col = "black",
     xlab = "Week", ylab = "Reported CHIKV cases",
     main = "Caldas Novas: predicted vs observed reported cases", xaxt = "n")
axis(1, at = x_ticks$week_index, labels = x_ticks$week)
lines(weeks, predicted_cases, lwd = 2, col = "#d6604d")
legend("topright", legend = c("Observed (SINAN)", "Predicted (MLE fit)"),
       col = c("black", "#d6604d"), lty = c(NA, 1), pch = c(16, NA),
       lwd = c(NA, 2), bty = "n")
par(mfrow = c(1, 1))

# ------------------------------------------------------------
# 8. Approximate confidence intervals from the Hessian
# ------------------------------------------------------------
# Standard errors come from inv(Hessian), the large-sample (Laplace) approximation.
# Only meaningful if optim converged cleanly and the Hessian is positive-definite.
if (fit$convergence == 0) {
  vcov_mat <- tryCatch(solve(fit$hessian), error = function(e) NULL)
  if (!is.null(vcov_mat)) {
    se <- sqrt(diag(vcov_mat))
    cat("\n--- Approximate 95% CIs ---\n")
    cat("rho:   fixed at", rho_fixed, "(only weakly identified; set externally)\n")
    cat("theta: [",
        round(1 + exp(fit$par[df_spline+1] - 1.96*se[df_spline+1]), 2), ", ",
        round(1 + exp(fit$par[df_spline+1] + 1.96*se[df_spline+1]), 2), "]\n", sep="")
  } else {
    cat("\nHessian could not be inverted; no CIs reported.\n")
  }
}

# ------------------------------------------------------------
# 9. Uncertainty bands: normal approximation to the (MAP) posterior
# ------------------------------------------------------------
library(MASS)  # for mvrnorm

# Draw parameters from a multivariate normal centred at the estimate with
# covariance = inv(Hessian) (the Laplace / large-sample approximation). We sample
# only the spline coefficients: rho is fixed and theta drives observation noise,
# not the mean curve, so beta_t and predicted cases are deterministic in the coefs.
idx      <- seq_len(df_spline)
vcov_sub <- vcov_mat[idx, idx, drop = FALSE]
set.seed(123)        # reproducibility
n_samples <- 500     # more samples = smoother bands, slower

param_samples <- mvrnorm(n = n_samples, mu = fit$par[idx], Sigma = vcov_sub)

# For each sample, run the SEIR and store predicted reported cases, TRUE infections + beta_t
pred_matrix <- matrix(NA, nrow = n_samples, ncol = T_weeks)
inf_matrix  <- matrix(NA, nrow = n_samples, ncol = T_weeks)   # all infections (true epidemic)
beta_matrix <- matrix(NA, nrow = n_samples, ncol = T_weeks)

cat("\nRunning", n_samples, "bootstrap samples...\n")
for (i in 1:n_samples) {
  spline_coefs_i <- param_samples[i, seq_len(df_spline)]
  beta_t_i       <- make_beta_t(spline_coefs_i)
  
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
      rho = rho_fixed,
      E0 = E0
    ),
    error = function(e) NULL
  )
  if (is.null(out_i)) next
  
  pred_matrix[i, ] <- colSums(out_i$new_reported)
  inf_matrix[i, ]  <- colSums(out_i$new_infections)   # all infections (reported + unreported)
  beta_matrix[i, ] <- beta_t_i
  if (i %% 50 == 0) cat("  ", i, "/", n_samples, "\n")
}

# Compute pointwise quantiles (2.5%, 50%, 97.5%) at each week
pred_lower  <- apply(pred_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
pred_median <- apply(pred_matrix, 2, quantile, probs = 0.500, na.rm = TRUE)
pred_upper  <- apply(pred_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)

inf_lower   <- apply(inf_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
inf_median  <- apply(inf_matrix, 2, quantile, probs = 0.500, na.rm = TRUE)
inf_upper   <- apply(inf_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)

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
  pred_hi  = pred_upper,
  inf_med  = inf_median,
  inf_lo   = inf_lower,
  inf_hi   = inf_upper
)

ggplot(plot_df, aes(x = week)) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  annotate("text", x = year_break, y = 0, label = "2026", angle = 90,
           vjust = -0.5, hjust = -11.5, fontface = "bold", size = 3.5) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = "Predicted (95% band)"),
              alpha = 0.25) +
  geom_line(aes(y = pred_med, colour = "Predicted"), linewidth = 1) +
  geom_point(aes(y = observed, colour = "Observed"), size = 1.5) +
  scale_colour_manual(name = NULL,
                      values = c("Observed"  = "black",
                                 "Predicted" = "#f4a582")) +
  scale_fill_manual(name = NULL,
                    values = c("Predicted (95% band)" = "#f4a582")) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week",
       y = "Reported CHIKV cases",
       title = "CHIKV cases in Caldas Novas (2025-2026)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 11),
    legend.position        = "inside", 
    legend.position.inside = c(0.32, 0.98), 
    legend.justification   = c(1, 1),  
    legend.background      = element_rect(fill = scales::alpha("white", 0.6), colour = NA),
    legend.spacing.y       = unit(0, "pt"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# ------------------------------------------------------------
# 10b. Estimated TRUE infections vs reported cases
# ------------------------------------------------------------
# new_infections = ALL infections (symptomatic + asymptomatic, reported + unreported).
# reported = rho * prop_symp * infections, so the true epidemic is ~1/(rho*prop_symp)
# (~4x at rho = 0.45) larger than the reported counts -- the "iceberg" below surveillance.
ggplot(plot_df, aes(x = week)) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  annotate("text", x = year_break + 0.6, y = max(plot_df$inf_hi) * 0.96,
           label = "2026", hjust = 0, fontface = "bold", size = 3.4, colour = "grey40") +
  geom_ribbon(aes(ymin = inf_lo, ymax = inf_hi, fill = "True infections (95% band)"),
              alpha = 0.25) +
  geom_line(aes(y = inf_med, colour = "True infections (estimated)"), linewidth = 1) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi, fill = "Reported (95% band)"),
              alpha = 0.30) +
  geom_line(aes(y = pred_med, colour = "Reported (predicted)"), linewidth = 1) +
  geom_point(aes(y = observed, colour = "Reported (observed)"), size = 1.3) +
  scale_colour_manual(name = NULL,
                      values = c("True infections (estimated)" = "#3182bd",
                                 "Reported (predicted)"        = "#d6604d",
                                 "Reported (observed)"         = "black")) +
  scale_fill_manual(name = NULL, guide = "none",   # bands match their line colours; no separate legend
                    values = c("True infections (95% band)" = "#a8d1e7",
                               "Reported (95% band)"        = "#f4a582")) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week",
       y = "Weekly cases",
       title = "Estimated true vs. reported cases in Caldas Novas (2025-2026)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    legend.position        = "inside",
    legend.position.inside = c(0.35, 0.98),
    legend.justification   = c(1, 1),
    legend.background      = element_rect(fill = scales::alpha("white", 0.6), colour = NA),
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
  # 95% band + median over the full window. The tail is weakly identified (few cases,
  # susceptibles depleted) and is pulled toward the prior median there.
  geom_ribbon(aes(ymin = beta_lo, ymax = beta_hi), fill = "#a8d1e7", alpha = 0.35) +
  geom_line(aes(y = beta_med), colour = "#3182bd", linewidth = 1) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  annotate("text", x = year_break, y = min(beta_df$beta_lo), label = "2026",
           angle = 90, vjust = -0.5, hjust = -11.5, fontface = "bold", size = 3.5) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  labs(x = "Week",
       y = expression(beta[t]),
       title = "Fitted weekly transmission rate for Caldas Novas (2025-2026)") +
  theme_bw(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5),
        axis.title       = element_text(size = 12),
        axis.text        = element_text(size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
  )