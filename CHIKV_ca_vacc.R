library(dplyr)
library(tidyr)
library(ggplot2)
library(splines)

# Shared helpers (week_to_index, fmtq, burden, load_burden_params) live here,
# deduplicated with MAYV_ca_vacc.R.
source("ca_common.R")

# ============================================================
# Caldas Novas (Goiás) CHIKV vaccination scenarios (IXCHIQ)
# ------------------------------------------------------------
# Two efficacy arms per scenario, following Hyolim's approach (as for Sete Lagoas):
#   - Disease-blocking ONLY (VE_inf = 0, VE_block = eff): conservative floor;
#     averts symptomatic/hospitalisations/deaths but NOT infections.
#   - Disease + infection blocking (VE_inf = eff, VE_block = eff): optimistic
#     ceiling; susceptibles move S->V so infections are averted too.
#
# Depends on the fitted objects from CHIKV_ca_pre_vacc_optim.R:
#   best_beta_t, best_rho, N, R_init_prop, A, age_df, I0, E0,
#   sigma, gamma, prop_symp, T_weeks, observed_cases, caldas_obs, weeks
# Sourced once if not already in memory (the fit + bootstrap is slow).
# ============================================================
if (!exists("best_beta_t") || !exists("caldas_obs")) {
  message("Fitted objects not found; sourcing CHIKV_ca_pre_vacc_optim.R ...")
  source("CHIKV_ca_pre_vacc_optim.R")
}

# ------------------------------------------------------------
# 1. SEIRV function with vaccination
# ------------------------------------------------------------
# S -> V transition: a fraction VE_inf of vaccinated susceptibles becomes fully
# immune `immun_delay` weeks after their dose (infection blocking). VE_block
# additionally reduces the symptomatic fraction in the covered population
# (disease blocking).
seirv_vaccinated <- function(
    T_weeks, A, N, R_init_prop, I0, E0,
    base_beta, sigma, gamma, rho,
    target_age,            # length-A binary vector; 1 = eligible
    total_coverage,        # final coverage of target pop
    weekly_delivery_speed, # fraction of supply delivered per week
    delay,                 # week_index at which vaccination begins
    VE_inf  = 0.989,       # infection-blocking efficacy
    VE_block = 0,          # disease-blocking efficacy
    immun_delay = 2,       # weeks from dose to immunity
    prop_symp = 0.5242478,
    sub_steps = 7
) {

  pmax0 <- function(x) pmax(0, x)
  N_total <- sum(N)
  dt      <- 1 / sub_steps

  S <- E <- I <- R <- V <- matrix(0, nrow = A, ncol = T_weeks)
  V_covered <- matrix(0, nrow = A, ncol = T_weeks)
  vacc_delayed <- matrix(0, nrow = A, ncol = T_weeks)
  coverage_frac <- matrix(0, nrow = A, ncol = T_weeks)
  new_infections <- new_symptomatic <- matrix(0, nrow = A, ncol = T_weeks)

  # Vaccine supply machinery
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  total_supply      <- target_pop * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed

  total_avail_age <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  unvaccinated   <- N

  # Working state at the start of week 1 (Caldas Novas seeds I0 only; E0 = 0).
  S_now <- pmax0(N - I0 - E0 - R_init_prop * N)
  E_now <- E0
  I_now <- I0
  R_now <- R_init_prop * N
  V_now <- rep(0, A)

  for (t in 1:T_weeks) {

    # ---- (a) People vaccinated `immun_delay` weeks ago become immune now
    prev_V_covered <- if (t > 1) V_covered[, t - 1] else rep(0, A)
    if (t - immun_delay >= 1) {
      effective_dose <- vacc_delayed[, t - immun_delay]
      immunized      <- round(VE_inf * effective_dose)   # only VE_inf fraction become immune
      V_covered[, t] <- prev_V_covered + effective_dose  # "covered" = received dose (for disease blocking)
    } else {
      immunized      <- rep(0, A)
      V_covered[, t] <- prev_V_covered
    }

    # Move immunized people S -> V (only when VE_inf > 0)
    S_now <- pmax0(S_now - immunized)
    V_now <- V_now + immunized

    coverage_frac[, t] <- V_covered[, t] / N

    # ---- (b) Allocate this week's doses (if past the start week `delay`)
    if (t >= delay && target_pop > 0) {
      rem <- weekly_dose_total
      for (a in target_idx) {
        alloc <- min(
          ceiling(weekly_dose_total * (N[a] / target_pop)),
          rem,
          unvaccinated[a],
          total_avail_age[a] - total_used_age[a]
        )
        if (alloc > 0) {
          prop_S    <- if (N[a] > 0) S_now[a] / N[a] else 0
          vacc_to_S <- round(alloc * prop_S)  # only S can be successfully vacc
          vacc_delayed[a, t] <- vacc_to_S
          total_used_age[a]  <- total_used_age[a] + alloc
          unvaccinated[a]    <- unvaccinated[a] - alloc
          rem                <- rem - alloc
        }
      }
    }

    # ---- (c) SEIR dynamics with sub-steps
    new_I_week <- rep(0, A)
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

    S[, t] <- S_now;  E[, t] <- E_now;  I[, t] <- I_now
    R[, t] <- R_now;  V[, t] <- V_now

    new_infections[, t]  <- new_I_week
    # Symptomatic cases reduced by disease-blocking efficacy in the covered pop
    new_symptomatic[, t] <- prop_symp * new_I_week *
      (1 - VE_block * coverage_frac[, t])
  }

  list(
    S = S, E = E, I = I, R = R, V = V,
    V_covered = V_covered, coverage_frac = coverage_frac,
    vacc_delayed = vacc_delayed, total_supply = total_supply,
    total_used_age = total_used_age,
    new_infections = new_infections, new_symptomatic = new_symptomatic,
    new_reported = rho * new_symptomatic
  )
}

# ============================================================
# 2. Vaccine programme parameters (common to all scenarios)
# ============================================================
# Eligible population: adults 18-59 = age groups 4..8.
age_group  <- as.character(age_df$age_group)
target_age <- rep(0, A)
target_age[c(4, 5, 6, 7, 8)] <- 1
target_pop <- sum(N[target_age == 1])

total_coverage        <- 0.30   # 30% of the eligible 18-59 population
weekly_delivery_speed <- 0.10   # 10% of supply / week => 10-week full rollout
immun_delay           <- 2      # 2 weeks from dose to immunity (IXCHIQ 14 days)
ixchiq_efficacy       <- 0.989  # IXCHIQ efficacy (Hyolim)

cat("Eligible 18-59 population:", round(target_pop), "\n")
cat("Target doses (30% coverage):", round(target_pop * total_coverage), "\n")
cat("Weekly delivery:", round(target_pop * total_coverage * weekly_delivery_speed),
    "doses/week over ~", round(1 / weekly_delivery_speed), "weeks\n")

# ------------------------------------------------------------
# 2b. Scenario start weeks (week_index within the 2025-W24 -> 2026-W22 window)
# ------------------------------------------------------------
# Calendar -> index map for this fit window:
#   2025-W24 = index 1  => 2025 weeks: index = week - 23
#   2026-W01 = index 31 => 2026 weeks: index = 30 + week
#   peak (2026-W09) = index 39
# week_to_index() is defined in ca_common.R.
start_s1 <- week_to_index(2026, 16)   # IXCHIQ real rollout, 18 Apr 2026 (2026-W16) -> 46
start_s2 <- week_to_index(2026, 1)    # start of 2026 (2026-W01)                     -> 31
start_s3 <- week_to_index(2025, 26)   # middle of 2025, before the outbreak          -> 3
cat(sprintf("Scenario start indices: S1=%d (2026-W16), S2=%d (2026-W01), S3=%d (2025-W26)\n",
            start_s1, start_s2, start_s3))
stopifnot(all(c(start_s1, start_s2, start_s3) >= 1),
          all(c(start_s1, start_s2, start_s3) <= T_weeks))

# ============================================================
# 3. Burden parameters (read from disease_progression.xlsx / Hyolim Table S4)
# ============================================================
# The spreadsheet holds the Beta(alpha, beta) hyperparameters (columns "Value 1" =
# alpha, "Value 2" = beta) for each disease-progression probability, by Group.
# We read alpha/beta here; the means feed the point estimates and the alpha/beta
# feed the Monte Carlo draws (section 7).
# Read the Beta(alpha, beta) hyperparameters and unpack the severity parameters
# (ps_*, hosp_*, cfr_*, age_to_band, cfr_vec) into the global environment. The
# loader and burden() live in ca_common.R (shared with MAYV_ca_vacc.R).
invisible(list2env(load_burden_params(A), globalenv()))

# Age re-weighting for DEATHS: correct the model's infection age split (population-
# structure driven) to the OBSERVED case age split from ca_combined. out_best is the
# fitted no-vaccine run (from CHIKV_ca_pre_vacc_optim.R); obs_band_prop comes from the
# shared loader. age_weight is applied by burden() (ca_common.R) to all scenarios, so
# the baseline reproduces the observed age distribution AND the vaccine's age-targeting
# is preserved. Set age_weight <- 1 to disable and recover the population-based deaths.
age_weight <- compute_age_weight(rowSums(out_best$new_infections), obs_band_prop, age_to_band)
cat(sprintf("Age re-weighting for deaths: w ranges [%.2f, %.2f] (1 = no correction)\n",
            min(age_weight), max(age_weight)))

# ============================================================
# 4. Run the no-vaccine baseline and 3 timings x 2 efficacy arms
# ============================================================
run_scenario <- function(coverage, start_week, VE_inf, VE_block,
                         base_beta = best_beta_t, prop_symp_use = prop_symp) {
  seirv_vaccinated(
    T_weeks = T_weeks, A = A, N = N, R_init_prop = R_init_prop,
    I0 = I0, E0 = E0, base_beta = base_beta,
    sigma = sigma, gamma = gamma, rho = best_rho,
    target_age = target_age, total_coverage = coverage,
    weekly_delivery_speed = weekly_delivery_speed, delay = start_week,
    VE_inf = VE_inf, VE_block = VE_block, immun_delay = immun_delay,
    prop_symp = prop_symp_use
  )
}

# No-vaccine counterfactual (coverage 0 -> reproduces the fitted epidemic)
out_baseline <- run_scenario(0, start_s3, VE_inf = 0, VE_block = 0)

timings <- list(
  "actual rollout" = start_s1,
  "start of 2026"  = start_s2,
  "pre-outbreak"   = start_s3
)
arms <- list(
  "Disease-blocking"             = c(VE_inf = 0,               VE_block = ixchiq_efficacy),
  "Disease + infection blocking" = c(VE_inf = ixchiq_efficacy, VE_block = ixchiq_efficacy)
)

scenarios <- list("No vaccine (baseline)" = out_baseline)
scen_meta <- data.frame(name = "No vaccine (baseline)", timing = "No vaccine",
                        arm = "No vaccine", start = NA_integer_)
for (tn in names(timings)) {
  for (an in names(arms)) {
    nm <- paste0(tn, " | ", an)
    scenarios[[nm]] <- run_scenario(total_coverage, timings[[tn]],
                                    arms[[an]]["VE_inf"], arms[[an]]["VE_block"])
    scen_meta <- rbind(scen_meta,
                       data.frame(name = nm, timing = tn, arm = an, start = timings[[tn]]))
  }
}

# ============================================================
# 5. Burden comparison
# ============================================================
burden_mat <- t(sapply(scenarios, burden))          # rows = scenario, cols = outcome
base_b     <- burden_mat["No vaccine (baseline)", ]

# Averted vs the no-vaccine baseline (drop the baseline row)
vac_rows    <- rownames(burden_mat) != "No vaccine (baseline)"
averted     <- sweep(-burden_mat[vac_rows, , drop = FALSE], 2, -base_b)   # base - scenario
pct_averted <- 100 * sweep(averted, 2, base_b, "/")

summary_tbl <- data.frame(
  scen_meta[match(rownames(averted), scen_meta$name), c("timing", "arm")],
  inf_averted    = round(averted[, "infections"]),
  symp_averted   = round(averted[, "symptomatic"]),
  hosp_averted   = round(averted[, "hospitalisations"], 1),
  deaths_averted = round(averted[, "deaths"], 2),
  pct_inf        = round(pct_averted[, "infections"], 1),
  pct_symp       = round(pct_averted[, "symptomatic"], 1),
  row.names = NULL
)

cat("\n=== Baseline burden (no vaccine) ===\n")
print(round(base_b, 1))
cat("\n=== Total burden by scenario ===\n")
print(round(burden_mat, 1))
cat("\n=== Averted vs no-vaccine baseline ===\n")
print(summary_tbl, row.names = FALSE)

# S3 (pre-outbreak) is the maximum-achievable ceiling; show each scenario's
# infections averted as a share of S3's, for the infection-blocking arm.
# Derive the name from `timings` (keyed by the pre-outbreak start week) so it
# survives any relabeling of the scenario names.
s3_timing <- names(timings)[which(unlist(timings) == start_s3)]
s3_both   <- paste0(s3_timing, " | Disease + infection blocking")
s3_inf    <- averted[s3_both, "infections"]
cat(sprintf("\nS3 (pre-outbreak, infection-blocking) averts %d infections = the ceiling.\n",
            round(s3_inf)))
both_rows <- summary_tbl$arm == "Disease + infection blocking"
cat("Infection-blocking arm, infections averted as % of S3 ceiling:\n")
print(data.frame(timing = summary_tbl$timing[both_rows],
                 pct_of_S3 = round(100 * summary_tbl$inf_averted[both_rows] / round(s3_inf), 1)),
      row.names = FALSE)

write.csv(summary_tbl, "caldas_vacc_averted.csv", row.names = FALSE)
cat("\nWrote caldas_vacc_averted.csv\n")

# ============================================================
# 6. Plots  ->  moved to CHIKV_outputs.R (presentation layer)
# ============================================================
# All figures (weekly infections, averted-burden bars, MC averted-burden bars, and
# the symptomatic epidemic-curve ribbons) are now produced by CHIKV_outputs.R, which
# consumes this engine's objects. This script computes scenarios + Monte Carlo and
# prints/writes the numeric tables only; source CHIKV_outputs.R afterwards for figures.

# ============================================================
# 7. Full Monte Carlo: averted burden with 95% uncertainty intervals
# ============================================================
# Propagates two sources of uncertainty, paired per draw so baseline and each
# scenario share the same draw (carries correlations through to the averted):
#   (a) calibration uncertainty in beta_t -> mvrnorm on the fitted covariance
#   (b) input-parameter uncertainty in prop_symp, hosp rate, age-specific CFR
#       -> Beta draws (Hyolim Table S4). Vaccine efficacy is fixed at IXCHIQ's.
library(MASS)

n_draws <- 1000
set.seed(2026)
param_samples <- mvrnorm(n = n_draws, mu = fit$par, Sigma = vcov_mat)

vac_names  <- setdiff(names(scenarios), "No vaccine (baseline)")
outcomes   <- c("infections", "symptomatic", "hospitalisations", "deaths")
base_draws <- matrix(NA_real_, n_draws, 4, dimnames = list(NULL, outcomes))
av_draws   <- setNames(lapply(vac_names, function(x)
                  matrix(NA_real_, n_draws, 4, dimnames = list(NULL, outcomes))), vac_names)

# Per-draw WEEKLY symptomatic trajectory (summed over age) for every scenario,
# for the epidemic-curve ribbon plot (95% UI) in CHIKV_outputs.R.
wk_symp <- setNames(lapply(names(scenarios), function(x)
                  matrix(NA_real_, n_draws, T_weeks)), names(scenarios))


cat("\nRunning", n_draws, "Monte Carlo draws (7 SEIR runs each)...\n")
for (i in 1:n_draws) {
  coefs_i  <- param_samples[i, 1:df_spline]
  beta_t_i <- make_beta_t(coefs_i)
  if (any(!is.finite(beta_t_i)) || any(beta_t_i > 20) || any(beta_t_i < 1e-6)) next

  # input-parameter draws
  ps_i      <- rbeta(1, ps_a, ps_b)
  hosp_i    <- rbeta(1, hosp_a, hosp_b)
  cfr_h_i   <- rbeta(9, cfr_hosp_a, cfr_hosp_b)
  cfr_n_i   <- rbeta(9, cfr_nonh_a, cfr_nonh_b)
  cfr_vec_i <- (hosp_i * cfr_h_i + (1 - hosp_i) * cfr_n_i)[age_to_band]

  out_b <- run_scenario(0, start_s3, 0, 0, base_beta = beta_t_i, prop_symp_use = ps_i)
  bb    <- burden(out_b, hosp_i, cfr_vec_i)
  base_draws[i, ] <- bb
  wk_symp[["No vaccine (baseline)"]][i, ] <- colSums(out_b$new_symptomatic)

  for (tn in names(timings)) for (an in names(arms)) {
    nm  <- paste0(tn, " | ", an)
    out <- run_scenario(total_coverage, timings[[tn]], arms[[an]]["VE_inf"],
                        arms[[an]]["VE_block"], base_beta = beta_t_i, prop_symp_use = ps_i)
    av_draws[[nm]][i, ] <- bb - burden(out, hosp_i, cfr_vec_i)
    wk_symp[[nm]][i, ]  <- colSums(out$new_symptomatic)
  }
  if (i %% 100 == 0) cat("  ", i, "/", n_draws, "\n")
}

# ---- Summarise: median (2.5% - 97.5%); fmtq() is defined in ca_common.R ----

# Reporting-rate draws for the TRUE/iceberg scaling. Centred on rho_fixed (best_rho)
# for consistency: we use rho = 0.40 for Caldas Novas (see CHIKV_ca_pre_vacc_optim.R),
# so we sample Beta(32, 48) (mean 0.40, 95% ~0.29-0.51) -- the same concentration (80)
# as Hyolim's national Beta(20, 60) [mean 0.25], re-centred to our municipality rate.
# Drawn as an independent vector so the model draws above -- and thus the anchored
# REPORTED column -- are unchanged; rho is independent of them.
set.seed(2027)
stopifnot(abs(best_rho - 0.40) < 1e-9)   # keep Beta mean == best_rho
rho_draws <- rbeta(n_draws, 32, 48)

# ---- Baseline (no-vaccine) burden: TRUE vs REPORTED, all four outcomes ----
# The reporting rate rho maps true -> reported (new_reported = rho * new_symptomatic),
# so its uncertainty belongs on the TRUE/iceberg side, NOT the reported side:
#   REPORTED = best_rho * base_draws -- the surveillance-visible burden, anchored by
#              calibration to the observed case total (rho-independent by construction,
#              so reported deaths stay put; its width is model-fit + severity only).
#   TRUE     = REPORTED / rho, with rho ~ Beta(20, 60) drawn per iteration (rho_draws).
#              This is the epidemiologically meaningful iceberg: a less-certain reporting
#              rate makes the true epidemic more uncertain (lower rho -> bigger iceberg).
# NB Hyolim's supplement specifies rho ~ Beta(20, 60) (mean 0.25, 95% 0.162-0.350). Its
# quantile UI is a little wider than the paper's headline "25% (20.1-32.5)" (that figure
# is ~Beta(50,141)); we use the stated generative Beta(20,60) for method consistency.
cat("\n=== Baseline burden, no vaccine: TRUE vs REPORTED (median, 95% UI) ===\n")
cat(sprintf("    Reported = best_rho x base (best_rho = %.2f); True = Reported / rho, rho ~ Beta(20,60).\n",
            best_rho))
cat(sprintf("    Observed reported cases = %s.\n",
            formatC(sum(observed_cases), big.mark = ",", format = "d")))
cat(sprintf("    %-16s  %-26s  %-26s\n", "Outcome", "True", "Reported"))
base_tbl <- data.frame(outcome = character(), true = character(),
                       reported = character(), stringsAsFactors = FALSE)
for (o in outcomes) {
  d      <- if (o == "deaths") 1 else 0
  rep_v  <- best_rho * base_draws[, o]           # anchored, rho-independent
  true_v <- rep_v / rho_draws                     # iceberg, carries reporting-rate uncertainty
  cat(sprintf("    %-16s  %-26s  %-26s\n", o, fmtq(true_v, d), fmtq(rep_v, d)))
  base_tbl <- rbind(base_tbl, data.frame(outcome = o,
                                         true = fmtq(true_v, d), reported = fmtq(rep_v, d)))
}
write.csv(base_tbl, "caldas_baseline_burden.csv", row.names = FALSE)
cat("Wrote caldas_baseline_burden.csv\n")

# Averted burden is a TRUE-scale iceberg quantity (real infections/deaths prevented),
# so it scales with the reporting rate the same way the totals do: multiply the
# fixed-rho averted draws by best_rho / rho_draws to carry the rho ~ Beta(20,60)
# uncertainty (so baseline_true - scenario_true = averted_true reconciles per draw).
# NB the % reduction (summary_tbl / epicurve annotations) is rho-INVARIANT -- the
# 1/rho factor cancels in averted/baseline -- so those are left unscaled.
f_rho    <- best_rho / rho_draws
av_true  <- function(nm, o) av_draws[[nm]][, o] * f_rho
mc_tbl <- data.frame(
  timing = scen_meta$timing[match(vac_names, scen_meta$name)],
  arm    = scen_meta$arm[match(vac_names, scen_meta$name)],
  Infections       = sapply(vac_names, function(nm) fmtq(av_true(nm, "infections"))),
  Symptomatic      = sapply(vac_names, function(nm) fmtq(av_true(nm, "symptomatic"))),
  Hospitalisations = sapply(vac_names, function(nm) fmtq(av_true(nm, "hospitalisations"), 1)),
  Deaths           = sapply(vac_names, function(nm) fmtq(av_true(nm, "deaths"), 2)),
  row.names = NULL, check.names = FALSE
)
cat("\n=== Averted vs no-vaccine baseline, TRUE scale (median, 95% UI) ===\n")
cat("    (outcomes: infections, symptomatic, hospitalisations, deaths; rho ~ Beta(20,60))\n")
old_width <- getOption("width"); options(width = 220)   # avoid wrapping the Deaths column
print(mc_tbl, row.names = FALSE)
options(width = old_width)
write.csv(mc_tbl, "caldas_vacc_averted_mc.csv", row.names = FALSE)
cat("\nWrote caldas_vacc_averted_mc.csv\n")

# ---- Averted-burden bar chart (95% UI) is drawn in CHIKV_outputs.R from av_draws. ----
cat("\nEngine done. Source CHIKV_outputs.R for figures + the Excel workbook.\n")
