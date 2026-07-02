library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
#  MAYV hypothetical-outbreak scenario for Sete Lagoas
#  Forward simulation (NO fitting): R0 is ASSUMED, not estimated,
#  because there is no observed MAYV outbreak to fit to.
#  -> the spline + negative-binomial + optim() machinery is gone.
# ============================================================

# ------------------------------------------------------------
# 1. SEIR function
# ------------------------------------------------------------
seir_baseline_MAYV <- function(
    T_weeks,        # number of weeks to simulate
    A,              # number of age groups
    N,              # vector length A: population by age group
    R_init_prop,    # vector length A: proportion already immune by age
    I0,             # vector length A: initial infections by age
    base_beta,      # vector length T_weeks: weekly transmission rate
    sigma,          # 1 / latent period (per week)
    gamma,          # recovery rate (per week)
    rho,            # reporting rate
    prop_symp,
    sub_steps = 7
) {
  pmax0 <- function(x) pmax(0, x)
  N_total <- sum(N)
  dt      <- 1 / sub_steps   # length of each sub-step in "weeks"
  
  S <- E <- I <- R <- matrix(0, nrow = A, ncol = T_weeks)
  new_infections <- new_symptomatic <- matrix(0, nrow = A, ncol = T_weeks)
  
  S_now <- N - I0 - R_init_prop * N
  E_now <- rep(0, A)
  I_now <- I0
  R_now <- R_init_prop * N
  
  S[, 1] <- S_now;  E[, 1] <- E_now;  I[, 1] <- I_now;  R[, 1] <- R_now
  
  for (t in 2:T_weeks) {
    new_I_week <- rep(0, A)
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
    new_reported    = rho * new_symptomatic
  )
}

# ------------------------------------------------------------
# 2. Natural-history parameters (MAYV)  -- day-based, converted to WEEKS
# ------------------------------------------------------------
infectious_days <- 6
gen_time_days   <- 15.2
latent_days     <- gen_time_days - infectious_days   # effective latent period

sigma <- 7 / latent_days        # E -> I rate per week  (= 0.76)
gamma <- 7 / infectious_days    # I -> R rate per week  (= 1.17)

prop_symp = 0.5242478
rho       <- 0.10               # reporting rate

# NOTE on the latent period: 9.2 d is NOT the human intrinsic incubation (3 d).
# A human-only SEIR has no mosquito compartments, so the E compartment must absorb
# the whole vector loop (extrinsic incubation etc.). Mean generation time of an
# SEIR ~ (1/sigma) + (1/gamma) = latent + infectious, so to reproduce Caicedo's
# 15.2-d generation time we set latent = 15.2 - infectious = 9.2 d.
# R0 and final size depend only on gamma (the infectious period), so this choice
# changes only the PEAK TIMING, not the cases-averted total.

# ------------------------------------------------------------
# 3. Load Sete Lagoas population; set MAYV baseline immunity
#    Population & age structure come from the sheet. We DO NOT use its
#    prop_ever_infected: the only MAYV serosurvey (Abad-Franch, Amazonian) is
#    flat with age with a high intercept -- the signature of alphavirus assay
#    cross-reactivity, not age-accumulating MAYV infection -- so its catalytic
#    FOI (0.0225) is not transferable to urban Sete Lagoas. Primary assumption:
#    urban Sete Lagoas is essentially MAYV-naive.
# ------------------------------------------------------------
age_df <- read_excel("population_MAYV.xlsx", sheet = "prop_immune")
age_df <- age_df[tolower(age_df$municipality) == "sete lagoas", ]
stopifnot(nrow(age_df) > 0)
age_df$age_group <- factor(age_df$age_group, levels = age_df$age_group)

N <- age_df$pop_num
A <- nrow(age_df)

# Naive-population baseline immunity. Primary = 0; sensitivity band: 0.00, 0.05, 0.10.
baseline_immunity <- 0.0
R_init_prop <- rep(baseline_immunity, A)
# To use age-specific values from the sheet instead, set:
#  R_init_prop <- age_df$prop_ever_infected

T_weeks <- 260    # 5-yr horizon: R0 = 1.1 is a slow burn and needs a long window to complete

# ------------------------------------------------------------
# 4. Seed the epidemic  (a SMALL number of introductions)
# ------------------------------------------------------------
# With a naive population, R_eff = R0 * S(0) ~ R0 > 1, so the epidemic is
# super-critical and the seed affects only TIMING, not final size. A small,
# realistic seed (a few introductions) is therefore fine. Avoid Mark's 1%
# proportion (~24,000 people) -- it is unrealistic as an introduction and would
# only matter near threshold, which we are no longer in.
I0_total <- 10                              
susceptible_pop <- N * (1 - R_init_prop)
I0 <- round(I0_total * susceptible_pop / sum(susceptible_pop))

# ------------------------------------------------------------
# 5. THRESHOLD PRE-CHECK  (do this BEFORE simulating)
#    R_eff = R0 * susceptible fraction. If < 1, no self-sustaining outbreak.
# ------------------------------------------------------------
S0_frac <- sum(susceptible_pop) / sum(N)            # mean susceptible fraction
R0_grid <- c(1.1, 1.2, 1.3)                         # near-threshold (outside Amazon)

diag_tbl <- data.frame(
  R0           = R0_grid,
  R_eff        = round(R0_grid * S0_frac, 3),
  takeoff_prob = round(ifelse(R0_grid * S0_frac > 1,
                              1 - 1 / (R0_grid * S0_frac), 0), 3)
)
cat("Mean baseline immunity:", round(1 - S0_frac, 3),
    " | susceptible fraction S(0):", round(S0_frac, 3), "\n")
cat("(takeoff_prob = chance a single introduction sustains; ignores heterogeneity,\n",
    " which would lower it further.)\n", sep = "")
print(diag_tbl)

# ------------------------------------------------------------
# 6. Run a scenario for any R0  (constant beta = R0 * gamma)
# ------------------------------------------------------------
run_scenario <- function(R0) {
  base_beta <- rep(R0 * gamma, T_weeks)
  seir_baseline_MAYV(
    T_weeks = T_weeks, A = A, N = N, R_init_prop = R_init_prop,
    I0 = I0, base_beta = base_beta,
    sigma = sigma, gamma = gamma, rho = rho, prop_symp = prop_symp
  )
}

# Analytical final size: fraction of initial susceptibles ever infected.
# Solves a = 1 - exp(-R_eff * a), R_eff = R0 * S(0). Horizon- and seed-independent,
# so this -- not the simulated sum -- is the authoritative attack rate (the
# simulated total can under-count if a slow R0 = 1.1 burn hasn't finished by T_weeks).
final_size <- function(R0, s0) {
  Reff <- R0 * s0
  if (Reff <= 1) return(0)
  uniroot(function(a) 1 - exp(-Reff * a) - a,
          lower = 1e-6, upper = 1 - 1e-9)$root
}

# Helper: largest reported count in any contiguous w-week window (the "outbreak year").
# Uses a cumulative sum so it is fast; falls back to the whole series if shorter than w.
window_max <- function(x, w = 52) {
  if (length(x) < w) return(sum(x))
  cs <- cumsum(c(0, x))
  max(cs[(w + 1):length(cs)] - cs[1:(length(cs) - w)])
}

summary_df <- do.call(rbind, lapply(R0_grid, function(R0) {
  out      <- run_scenario(R0)
  wk_rep   <- colSums(out$new_reported)      # reported cases per week
  wk_inf   <- colSums(out$new_infections)    # TOTAL new infections per week (all ages)
  sim_inf  <- sum(out$new_infections)        # full-epidemic total infections
  a_final  <- final_size(R0, S0_frac)
  peak_yr  <- window_max(wk_rep, 52)
  data.frame(
    R0               = R0,
    R_eff            = round(R0 * S0_frac, 3),
    attack_rate_susc = round(a_final, 4),
    total_reported   = round(a_final * sum(susceptible_pop) * prop_symp * rho),
    reported_yr1     = round(sum(wk_rep[1:52])),                 # reported, weeks 1-52
    infections_yr1   = round(sum(wk_inf[1:52])),                 # TOTAL infections, weeks 1-52  <- what you wanted
    symptomatic_yr1  = round(sum(wk_inf[1:52]) * prop_symp),     # symptomatic, weeks 1-52
    peak_52wk        = round(peak_yr),
    peak_inc_per100k = round(peak_yr / sum(N) * 1e5, 1),
    peak_week        = which.max(wk_rep),
    sim_infections   = round(sim_inf)
  )
}))
cat("\n--- Scenario summary ---\n")
cat("Columns over the FULL epidemic: attack_rate_susc, total_reported (analytical).\n",
    "Columns on a one-year basis (for comparison with CHIKV's single 52-wk season):\n",
    "  reported_yr1  = reported in weeks 1-52 from introduction (slow R0 = few cases here)\n",
    "  peak_52wk     = reported in the busiest contiguous 52-week window (the 'outbreak year')\n",
    "  peak_inc_per100k = that busiest-year count per 100,000 population\n",
    "CHIKV 2024 Sete Lagoas for reference: 21,272 reported in one season.\n",
    "(If sim_infections << total_reported/(prop_symp*rho), the slow burn had not finished by T_weeks.)\n",
    sep = "")
print(summary_df)

# ------------------------------------------------------------
# 7. Plot reported epidemic curves across R0 scenarios
# ------------------------------------------------------------
curve_df <- do.call(rbind, lapply(R0_grid, function(R0) {
  data.frame(week     = 1:T_weeks,
             reported = colSums(run_scenario(R0)$new_reported),
             R0       = factor(R0))
}))

ggplot(curve_df, aes(week, reported, colour = R0)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(0, T_weeks, 52)) +
  labs(x = "Week", y = "Predicted reported MAYV cases",
       colour = expression(R[0]),
       title = "Hypothetical MAYV outbreak in Sete Lagoas") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", colour = "grey80"))

# ------------------------------------------------------------
# 8. (Optional) Propagate R0 uncertainty over the 1.1-1.3 range
#    Scenario-equivalent of your CHIKV bootstrap: sample the uncertain INPUT
#    (R0), run the model, summarise the spread -- NOT a fitted posterior.
# ------------------------------------------------------------
set.seed(123)
n_draw   <- 500
R0_draws <- runif(n_draw, 1.1, 1.3)         # flat over Caicedo's outside-Amazon range
pred_mat <- matrix(NA, n_draw, T_weeks)
for (i in seq_len(n_draw)) {
  pred_mat[i, ] <- colSums(run_scenario(R0_draws[i])$new_reported)
}
band <- data.frame(
  week = 1:T_weeks,
  lo   = apply(pred_mat, 2, quantile, 0.025),
  med  = apply(pred_mat, 2, quantile, 0.500),
  hi   = apply(pred_mat, 2, quantile, 0.975)
)
ggplot(band, aes(week)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#a8d1e7", alpha = 0.35) +
  geom_line(aes(y = med), colour = "#2c7fb8", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, T_weeks, 52)) +
  labs(x = "Week", y = "Reported MAYV cases (modelled)",
       title = "MAYV scenario, R0 ~ Uniform(1.1, 1.3): 95% band") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", colour = "grey80"))