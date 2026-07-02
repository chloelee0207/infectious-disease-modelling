# ============================================================
# Burden table via Monte Carlo propagation
# Combines: (a) calibration uncertainty in beta_t / rho  -> Hessian mvrnorm
#           (b) input-parameter uncertainty             -> Beta draws (your xlsx)
# Outputs: infections / symptomatic / deaths AVERTED, with 95% UI,
#          for disease-blocking and disease+infection-blocking.
#
# Assumes the following already exist from your fit script:
#   fit, vcov_mat, df_spline, make_beta_t(), seir_baseline(),
#   T_weeks, A, N, R_init_prop, I0, sigma, gamma
# and that you have a way to run a vaccinated scenario (your seirv code).
# ============================================================

library(MASS)   # mvrnorm
set.seed(2026)

n_draws <- 1000

# ------------------------------------------------------------
# 1. Input-parameter distributions (from disease_progression.xlsx)
#    Each is Beta(alpha, beta) as you extracted from Hyolim's Table S4.
# ------------------------------------------------------------
# prop_symp: probability symptomatic among infections (Overall)
ps_a   <- 35.84;  ps_b   <- 32.56

# hospitalisation among symptomatic (single, not age-specific)
hosp_a <- 59;     hosp_b <- 1415.2

# case fatality, decadal bands: [0,10) [10,20) ... [80,90)  -> length 9
cfr_hosp_a <- c( 26.7,  133.0,  407.9,  493.2,  789.5,  944.9, 1805.2, 4466.9, 3801.2)
cfr_hosp_b <- c(1455.7,8443.6,22412.5,25471.5,33197.9,34239.3,37879.3,36999.1,16439.2)
cfr_nonh_a <- c(  9.7,    6.1,    2.6,    9.7,    3.4,   16.2,   19.9,   27.5,  39.98)
cfr_nonh_b <- c(30950.3,54677.4,70752.9,93187.2,76516.7,79639.98,54639.2,28387.1,9914.8)

# ------------------------------------------------------------
# 2. Crosswalk: map each of your 12 model age groups to a decadal CFR band (1..9)
#    Your 12 groups (from population.xlsx, prop_immune sheet):
#      1  <1      -> band 1  [0,10)
#      2  1-9     -> band 1  [0,10)
#      3  10-17   -> band 2  [10,20)
#      4  18-19   -> band 2  [10,20)   (decadal banding puts 18-19 in [10,20))
#      5  20-29   -> band 3  [20,30)
#      6  30-39   -> band 4  [30,40)
#      7  40-49   -> band 5  [40,50)
#      8  50-59   -> band 6  [50,60)
#      9  60-69   -> band 7  [60,70)
#      10 70-79   -> band 8  [70,80)
#      11 80-89   -> band 9  [80,90)
#      12 90-100  -> band 9  [80,90)   (no band beyond 80-89; use oldest available)
age_to_band <- c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9)
stopifnot(!any(is.na(age_to_band)), length(age_to_band) == A)

# ------------------------------------------------------------
# 3. Burden extractor: turn a scenario output into inf / symp / deaths
#    Expects out$new_infections and out$new_symptomatic as A x T_weeks matrices.
#    cfr_vec is the age-specific death-per-symptomatic-case rate (length A).
# ------------------------------------------------------------
burden <- function(out, cfr_vec) {
  symp_age <- rowSums(out$new_symptomatic)          # symptomatic by age
  list(
    inf   = sum(out$new_infections),
    symp  = sum(symp_age),
    death = sum(symp_age * cfr_vec)
  )
}

# ------------------------------------------------------------
# 4. Scenario runner. Baseline uses your seir_baseline() directly.
#    For the vaccinated arms, wrap YOUR seirv_vaccinated code in run_vacc().
#    It must return an object with $new_infections and $new_symptomatic.
# ------------------------------------------------------------
run_vacc <- function(beta_t, prop_symp, VE_inf, VE_block) {
  # >>> PLUG IN your seirv_vaccinated() call here, passing beta_t, prop_symp,
  #     VE_inf, VE_block, and your fixed rollout settings (2-wk delay,
  #     10%/wk, Fortaleza dose share). Return its output list. <<<
  stop("run_vacc(): insert your seirv_vaccinated call")
}

# ------------------------------------------------------------
# 5. Draw calibration parameters jointly (carries beta_t<->rho correlation)
# ------------------------------------------------------------
param_samples <- mvrnorm(n = n_draws, mu = fit$par, Sigma = vcov_mat)

# storage
res <- data.frame(
  base_inf   = NA_real_, base_symp  = NA_real_, base_death = NA_real_,
  inf_av_both  = NA_real_,                                  # disblock averts 0 infections
  symp_av_disb = NA_real_, symp_av_both = NA_real_,
  death_av_disb= NA_real_, death_av_both= NA_real_
)[rep(1, n_draws), ]

# ------------------------------------------------------------
# 6. Monte Carlo loop
# ------------------------------------------------------------
for (i in 1:n_draws) {

  # --- (a) transmission from the Hessian draw ---
  coefs_i  <- param_samples[i, 1:df_spline]
  beta_t_i <- make_beta_t(coefs_i, T_weeks, df_spline)
  if (any(!is.finite(beta_t_i)) || any(beta_t_i > 20) || any(beta_t_i < 1e-6)) next

  # --- (b) input-parameter draws ---
  ps_i    <- rbeta(1, ps_a, ps_b)
  hosp_i  <- rbeta(1, hosp_a, hosp_b)
  cfr_h_i <- rbeta(9, cfr_hosp_a, cfr_hosp_b)
  cfr_n_i <- rbeta(9, cfr_nonh_a, cfr_nonh_b)
  cfr_band <- hosp_i * cfr_h_i + (1 - hosp_i) * cfr_n_i   # death per symptomatic, by band
  cfr_vec  <- cfr_band[age_to_band]                       # expand to A age groups

  # --- (c) run the three scenarios with this beta_t and prop_symp ---
  out_base <- seir_baseline(
    T_weeks = T_weeks, A = A, N = N, R_init_prop = R_init_prop, I0 = I0,
    base_beta = beta_t_i, sigma = sigma, gamma = gamma,
    rho = 1, prop_symp = ps_i                # rho irrelevant for burden; set 1
  )
  out_disb <- run_vacc(beta_t_i, ps_i, VE_inf = 0,     VE_block = 0.989)
  out_both <- run_vacc(beta_t_i, ps_i, VE_inf = 0.989, VE_block = 0.989)

  # --- (d) burden + averted ---
  b  <- burden(out_base, cfr_vec)
  d  <- burden(out_disb, cfr_vec)
  bo <- burden(out_both, cfr_vec)

  res$base_inf[i]    <- b$inf
  res$base_symp[i]   <- b$symp
  res$base_death[i]  <- b$death
  res$inf_av_both[i]   <- b$inf   - bo$inf       # disease-blocking infections averted = 0
  res$symp_av_disb[i]  <- b$symp  - d$symp
  res$symp_av_both[i]  <- b$symp  - bo$symp
  res$death_av_disb[i] <- b$death - d$death
  res$death_av_both[i] <- b$death - bo$death

  if (i %% 100 == 0) cat("  ", i, "/", n_draws, "\n")
}

# ------------------------------------------------------------
# 7. Summarise: median (2.5%-97.5%)
# ------------------------------------------------------------
fmt <- function(x, d = 0) {
  q <- quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE)
  sprintf("%s (%s-%s)",
          formatC(round(q[1], d), big.mark = ",", format = "f", digits = d),
          formatC(round(q[2], d), big.mark = ",", format = "f", digits = d),
          formatC(round(q[3], d), big.mark = ",", format = "f", digits = d))
}
pct <- function(av, base) {
  r <- 100 * av / base
  q <- quantile(r, c(0.5, 0.025, 0.975), na.rm = TRUE)
  sprintf("%.1f%% (%.1f-%.1f%%)", q[1], q[2], q[3])
}

burden_table <- data.frame(
  Outcome = c("Infections averted", "Symptomatic cases averted", "Deaths averted",
              "% reduction (symptomatic)"),
  `Disease-blocking only` = c(
    "-",                                   # infections averted = 0 by construction
    fmt(res$symp_av_disb),
    fmt(res$death_av_disb, 1),
    pct(res$symp_av_disb, res$base_symp)
  ),
  `Disease + infection blocking` = c(
    fmt(res$inf_av_both),
    fmt(res$symp_av_both),
    fmt(res$death_av_both, 1),
    pct(res$symp_av_both, res$base_symp)
  ),
  check.names = FALSE
)

cat("\n--- Fortaleza burden table (median, 95% UI) ---\n")
print(burden_table, row.names = FALSE)

# baseline burden (for the "before vaccination" sentence)
cat("\nBaseline (no vaccine):\n")
cat("  Infections:  ", fmt(res$base_inf),   "\n")
cat("  Symptomatic: ", fmt(res$base_symp),  "\n")
cat("  Deaths:      ", fmt(res$base_death, 1), "\n")
