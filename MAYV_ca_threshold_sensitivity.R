# ============================================================
# MAYV_ca_threshold_sensitivity.R
# ------------------------------------------------------------
# Two structural sensitivity analyses for the MAYV model, both propagated over the
# SAME 1000-draw ensemble the main results use (gamma, sigma, rho, prop_symp and prior
# immunity per draw, paired with freshly drawn vaccine parameters), so the 95% UIs are
# comparable with the headline table.
#
# The vaccine draws are read from MAYV_ca_engine_results_high.rds rather than redrawn, so
# the R0 = 2.5 / seed-week-1 cell reproduces the headline result exactly.
#
#   1. TRANSMISSION POTENTIAL. Baseline and averted symptomatic cases across fixed R0,
#      reported alongside the reproduction numbers the model actually runs at. R0 here
#      is the WET-SEASON PEAK (the envelope is peak-normalised), so the operative
#      quantities are the effective reproduction numbers: mean R_eff over the 52 weeks,
#      peak R_eff, and the number of weeks with R_eff > 1. Mean R0(t) is deliberately
#      omitted -- it is just R_eff without the susceptible factor, and R_eff is what
#      determines whether transmission grows.
#
#   2. SEEDING WEEK. The index case is introduced at an arbitrary week, since there is
#      no observed MAYV introduction to calibrate to. Absolute burden is highly
#      sensitive to that choice; the % averted is not. Seed SIZE is held at one
#      infectious person throughout -- larger seeds are not a plausible base case.
#
# Run after: MAYV_ca_lhs.R (writes MAYV_ca_lhs_ensemble_high.rds)
# Output:    MAYV_ca_threshold_sensitivity.xlsx
# ============================================================
DEFS_ONLY <- TRUE
suppressMessages({library(writexl); source("ca_common.R"); source("MAYV_ca_engine.R")})

E <- readRDS("MAYV_ca_lhs_ensemble_high.rds")
N <- E$N; A <- E$A; T_weeks <- E$T_weeks; ND <- length(E$gamma)
season <- { v <- as.numeric(readRDS("caldas_hybrid_season.rds")); v / max(v) }   # peak-normalised
target_age <- rep(0, A); target_age[4:8] <- 1
start_pre <- 17L; immun_delay <- 2                        # 2025-W40, matching CHIKV

# Vaccine parameters are READ FROM THE ENGINE, not redrawn. Redrawing from the same
# distributions gives a different Monte Carlo realisation, which would make the R0 = 2.5
# / seed-week-1 row disagree with the headline results by a case or two -- a discrepancy
# a reader would reasonably query. Reusing the engine's draws makes that row reproduce
# the main analysis exactly, so the table is a strict extension of it.
Mres <- readRDS("MAYV_ca_engine_results_high.rds")
stopifnot(!is.null(Mres$del_d))          # re-run MAYV_ca_engine.R if this fails
cov_d <- Mres$cov_d; veb_d <- Mres$veb_d; del_d <- Mres$del_d; dly_d <- Mres$delay_d
stopifnot(length(cov_d) == ND, length(veb_d) == ND, length(del_d) == ND, length(dly_d) == ND)

wk_lab <- function(i) { o <- 23 + i
  if (o <= 53) sprintf("2025-W%02d", o) else sprintf("2026-W%02d", o - 53) }
q3  <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
fmt <- function(x, dp = 0) { q <- q3(x); sprintf("%s (%s-%s)",
        formatC(round(q[1], dp), big.mark=",", format="f", digits=dp),
        formatC(round(q[2], dp), big.mark=",", format="f", digits=dp),
        formatC(round(q[3], dp), big.mark=",", format="f", digits=dp)) }

# one draw: baseline and vaccinated symptomatic totals
one <- function(i, R0, seed_week) {
  imm <- E$immune_frac[i]; Rimm <- rep(imm, A); sus <- N * (1 - imm)
  bt  <- R0 * E$gamma[i] * season
  sim <- function(cov, vb) seirv_vaccinated_MAYV(
    T_weeks, A, N, Rimm, E$I0_total * sus/sum(sus), bt, E$sigma[i], E$gamma[i], E$rho[i],
    target_age, cov, del_d[i], start_pre + dly_d[i], VE_inf = 0, VE_block = vb,
    immun_delay = immun_delay, prop_symp = E$prop_symp[i], E0 = E$E0, seed_week = seed_week)
  b <- sim(0, 0); v <- sim(cov_d[i], veb_d[i]); psi <- E$prop_symp[i]
  c(base = psi * sum(b$new_infections),
    vacc = sum(colSums(psi * v$new_infections * (1 - veb_d[i] * v$coverage_frac))))
}
sweep_cfg <- function(R0, seed_week) {
  m <- vapply(seq_len(ND), function(i) one(i, R0, seed_week), numeric(2))
  list(base = m["base", ], averted = m["base", ] - m["vacc", ],
       pct = 100 * (m["base", ] - m["vacc", ]) / m["base", ])
}

# ---- 1. transmission potential ----------------------------------------------
R0_GRID  <- c(1.2, 2.0, 2.4, 2.5, 2.7, 3.0)
R0_LABEL <- c("1.2 (incidental transmission scenario)", "2.0", "2.4",
              "2.5 (future urban-adapted transmission scenario)", "2.7", "3.0")
# R0 is fixed, but R_eff = R0 x season(t) x S/N inherits the sampled prior immunity, so
# it DOES carry a 95% UI. Reported at the start of the epidemic; depletion lowers it
# further as the outbreak runs, though at MAYV's attack rates that shift is negligible.
SN_d <- 1 - E$immune_frac; imm_med <- median(E$immune_frac); SN <- median(SN_d)
reff <- function(x, d = 3) sprintf("%.*f (%.*f-%.*f)", d, median(x), d,
                                   quantile(x, .025), d, quantile(x, .975))
cat(sprintf("Propagating %d draws over %d R0 values...\n", ND, length(R0_GRID)))
tp <- do.call(rbind, lapply(seq_along(R0_GRID), function(k) {
  R0 <- R0_GRID[k]; r <- sweep_cfg(R0, E$seed_week)
  cat(sprintf("  R0 = %.1f done\n", R0))
  data.frame(fixed_R0 = R0_LABEL[k],
             mean_R_eff = reff(R0 * mean(season) * SN_d),
             peak_R_eff = reff(R0 * max(season)  * SN_d),
             weeks_R_eff_above_1 = { w <- sapply(SN_d, function(z) sum(R0 * season * z > 1))
                                     sprintf("%d (%d-%d)", median(w), quantile(w, .025),
                                             quantile(w, .975)) },
             baseline_symptomatic = fmt(r$base),
             averted_symptomatic  = fmt(r$averted),
             pct_averted          = fmt(r$pct, 1),
             stringsAsFactors = FALSE) }))

# ---- 2. seeding week --------------------------------------------------------
SEED_GRID <- c(1, 9, 17, 25)
cat("Propagating seeding weeks...\n")
sd <- do.call(rbind, lapply(SEED_GRID, function(sw) {
  r <- sweep_cfg(2.5, sw); cat(sprintf("  seed week %d done\n", sw))
  data.frame(seed_week_index = sw, seed_week = wk_lab(sw),
             season_at_seeding = round(season[sw], 3),
             R_eff_at_seeding  = round(2.5 * season[sw] * SN, 3),
             baseline_symptomatic = fmt(r$base),
             averted_symptomatic  = fmt(r$averted),
             pct_averted          = fmt(r$pct, 2),
             stringsAsFactors = FALSE) }))

notes <- data.frame(item = c(
  "Scope", "R0 convention", "Why R_eff and not mean R0(t)", "Susceptible factor",
  "R_eff uncertainty", "Comparison with published R0",
  "Seed size", "Seeding finding", "Draws", "Vaccine"),
  detail = c(
  "Structural sensitivity for MAYV. Both sheets use the same 1000-draw ensemble as the main results.",
  "R0 is the WET-SEASON PEAK: the seasonal envelope is normalised to a peak of 1, so R0(t) = R0 x season(t) <= R0 at all times.",
  "Mean R0(t) is mean R_eff without the susceptible factor. R_eff is what determines whether transmission grows, so only R_eff is tabulated.",
  sprintf("R_eff = R0 x season(t) x S/N, evaluated at the START of the epidemic. Prior immunity is %.1f%% (95%% UI 3.0-18.0%%), so S/N = %.3f (0.820-0.970). Depletion lowers R_eff further as the epidemic runs, though negligibly at MAYV attack rates.", 100*imm_med, SN),
  "R0 is fixed per scenario and has no interval by construction, but R_eff does: its 95% UI comes entirely from the sampled prior immunity, the only stochastic term in R0 x season(t) x S/N.",
  sprintf("Caicedo et al. 2021 derive R0 from age-stratified seroprevalence using catalytic models (P(a) = 1-exp(-lambda*a)), so their 2.1-2.9 for the Amazon is an ENDEMIC-AVERAGE reproduction number, not a seasonal peak. Taking 2.50 as the PEAK is therefore conservative: it implies a season-averaged R0 of %.2f, and reproducing 2.50 as an annual mean would need a peak of %.2f. The peak is instead comparable to their outbreak-derived estimate, 2.2 (95%% CrI 0.8-4.8) from the 1954-55 Santa Cruz epidemic.", 2.5*mean(season), 2.5/mean(season)),
  "Held at ONE infectious person in every row. Larger seeds were used only as a regime diagnostic and are not a plausible base case.",
  "Absolute burden is highly sensitive to seeding week; % averted is nearly invariant, because the vaccine's effect depends on the timing overlap between the epidemic and the coverage curve, and the seasonal envelope pins the peak regardless of seeding.",
  sprintf("%d. Transmission draws (gamma, sigma, rho, prop_symp, prior immunity) come from the ensemble and vaccine draws from the engine, so the R0 = 2.5 / seed-week-1 cell reproduces the headline result exactly.", ND),
  "Pre-outbreak campaign at 2025-W40, coverage of eligible 18-59 Beta(30%, 20-40%), disease-blocking efficacy Beta(50%, 25-75%), VE_inf = 0."),
  stringsAsFactors = FALSE)


# ---- 3. the seasonal envelope itself -----------------------------------------
# Documents season(t) week by week, and the arithmetic behind mean(season) = 0.5422.
# caldas_hybrid_season.rds is stored MEAN-1 (mean = 1 by construction). The engine
# peak-normalises it, season = v / max(v), so
#     mean(season) = mean(v) / max(v) = 1 / max(v) = 1 / 1.8443 = 0.5422
# i.e. the 0.5422 is not a separate estimate -- it is the reciprocal of the mean-1
# envelope's peak, and follows entirely from the peak-normalisation convention.
v_mean1 <- as.numeric(readRDS("caldas_hybrid_season.rds"))
SPLICE  <- 40L                                   # 2026-W10; weeks 1-40 CHIKV beta, 41-52 rainfall
env <- data.frame(
  week_index      = seq_along(v_mean1),
  epi_week        = vapply(seq_along(v_mean1), wk_lab, character(1)),
  source          = ifelse(seq_along(v_mean1) <= SPLICE,
                           "fitted CHIKV beta_t", "CHIRPS rainfall climatology"),
  envelope_mean1  = round(v_mean1, 4),
  season_peaknorm = round(v_mean1 / max(v_mean1), 4),
  R0_t_at_2.5     = round(2.5 * v_mean1 / max(v_mean1), 3),
  R_eff_t_at_2.5  = round(2.5 * v_mean1 / max(v_mean1) * SN, 3),
  stringsAsFactors = FALSE)
env_summary <- data.frame(quantity = c(
  "mean of the stored mean-1 envelope", "max of the stored mean-1 envelope",
  "mean(season) after peak-normalisation", "derivation",
  "min(season)", "weeks with season > 0.5", "peak week",
  "source split", "why the tail is rainfall"),
  value = c(
  sprintf("%.6f (1 by construction)", mean(v_mean1)),
  sprintf("%.6f", max(v_mean1)),
  sprintf("%.4f", mean(v_mean1)/max(v_mean1)),
  sprintf("mean(season) = mean(v)/max(v) = 1/%.4f = %.4f -- a consequence of peak-normalising, not a separate estimate", max(v_mean1), 1/max(v_mean1)),
  sprintf("%.4f", min(v_mean1)/max(v_mean1)),
  sprintf("%d of %d", sum(v_mean1/max(v_mean1) > 0.5), length(v_mean1)),
  sprintf("week %d (%s)", which.max(v_mean1), wk_lab(which.max(v_mean1))),
  sprintf("weeks 1-%d (%.0f%%) fitted CHIKV beta_t; weeks %d-%d (%.0f%%) rainfall climatology, rescaled to join continuously at 2026-W10",
          SPLICE, 100*SPLICE/52, SPLICE+1, length(v_mean1), 100*(52-SPLICE)/52),
  "The fitted CHIKV beta_t tail is confounded with susceptible depletion, not short of data: 39.5% of observed CHIKV cases fall at or after the join. Depletion alone explains the case decline, so beta_t stays flat (0.638 -> 0.630) and carries no seasonal signal. Rainfall is exogenous to the epidemic."),
  stringsAsFactors = FALSE)

write_xlsx(list(transmission_potential = tp, seeding_week = sd,
                envelope = env, envelope_summary = env_summary, notes = notes),
           "MAYV_ca_threshold_sensitivity.xlsx")
cat("\nWrote MAYV_ca_threshold_sensitivity.xlsx\n")
print(tp[, c("fixed_R0","mean_R_eff","peak_R_eff","weeks_R_eff_above_1","baseline_symptomatic","averted_symptomatic")], row.names = FALSE)
cat("\n"); print(sd[, c("seed_week","R_eff_at_seeding","baseline_symptomatic","averted_symptomatic","pct_averted")], row.names = FALSE)
