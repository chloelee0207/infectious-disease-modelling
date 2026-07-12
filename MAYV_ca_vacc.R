library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Shared helpers (week_to_index, fmtq, burden, load_burden_params) live here,
# deduplicated with CHIKV_ca_vacc.R.
source("ca_common.R")

# ============================================================
# Caldas Novas (Goias) MAYV HYPOTHETICAL vaccination -- pre-outbreak campaign
# ------------------------------------------------------------
# Combines two earlier pieces:
#   * MAYV_ca_pre_vacc.R  -- the forward (un-fitted) MAYV outbreak: R0 is ASSUMED
#     (Caicedo outside-Amazon 1.1-1.3), MAYV rides the FITTED CHIKV seasonal
#     envelope, and a single index case is seeded at the wet-season onset.
#   * CHIKV_ca_vacc.R     -- the vaccination machinery (S->V rollout, two efficacy
#     arms, averted-burden tables/plots, paired Monte Carlo).
#
# DESIGN: a PRE-OUTBREAK campaign (vaccinate before the wet-season seed) crossed
# over two axes:
#   * OUTBREAK SIZE (R0), fixed:  Current 1.2 (Caicedo) and Future 3.51 (Dodero-Rojas
#     urban Ae. aegypti upper bound).
#   * VACCINE: DISEASE-BLOCKING ONLY (VE_inf = 0 in EVERY scenario; VE_block = eff).
#     Such a vaccine reduces symptoms/hospitalisations/deaths but NEVER blocks
#     infection, so total infections are identical to no-vaccine. Central efficacy
#     VE_block = 50%; 25% and 75% are the efficacy sensitivity levels.
# MAYV_outputs.R turns the engine objects into three plots ((a) total burden
# [no vaccine vs 50%], (b) averted [50%], (c) averted efficacy sensitivity
# [25/50/75%]) and the Excel workbooks. R0 1.1-1.3 is a separate sensitivity.
#
# IMPORTANT CAVEATS
#   1. There is NO licensed Mayaro vaccine. Efficacy is a HYPOTHETICAL disease-
#      blocking value, exposed as `mayv_vacc_efficacy` (central) / `mayv_vacc_eff_sens`.
#   2. There is NO MAYV-specific severity distribution. To show hospitalisations
#      and deaths averted we BORROW the CHIKV disease-progression probabilities
#      (Hyolim Table S4, disease_progression.xlsx) as an alphavirus proxy. MAYV is
#      generally milder than CHIKV, so these hospitalisation/death figures are an
#      UPPER BOUND; read them as "CHIKV-equivalent severity", not MAYV truth.
#
# Depends on objects built by MAYV_ca_pre_vacc.R (sourced once if absent):
#   seir_baseline_MAYV, sigma, gamma, prop_symp, rho, N, A, age_df,
#   R_init_prop, susceptible_pop, S0_frac, T_weeks, season, I0, E0, seed_week,
#   R0_central, R0_sens, wk_num, x_ticks, year_break, wet_start, wet_end,
#   season_layers
# ============================================================
if (!exists("season") || !exists("seir_baseline_MAYV") || !exists("seed_week")) {
  message("MAYV pre-vacc objects not found; sourcing MAYV_ca_pre_vacc.R ...")
  source("MAYV_ca_pre_vacc.R")
}

# ------------------------------------------------------------
# 1. SEIRV function: MAYV SEIR + vaccination + delayed seed
# ------------------------------------------------------------
# This is seirv_vaccinated() (CHIKV) with two MAYV-specific features carried over
# from seir_baseline_MAYV(): the index case I0 is introduced at `seed_week` (not at
# t = 1), and the population is otherwise infection-free before then. With
# coverage = 0 this reproduces the pre-vacc baseline EXACTLY.
seirv_vaccinated_MAYV <- function(
    T_weeks, A, N, R_init_prop, I0,
    base_beta, sigma, gamma, rho,
    target_age,            # length-A binary vector; 1 = eligible
    total_coverage,        # final coverage of target pop
    weekly_delivery_speed, # fraction of supply delivered per week
    delay,                 # week_index at which vaccination begins
    VE_inf  = 0.989,       # infection-blocking efficacy
    VE_block = 0,          # disease-blocking efficacy
    immun_delay = 2,       # weeks from dose to immunity
    prop_symp = 0.5242478,
    sub_steps = 7,
    E0 = rep(0, A),        # initial exposed by age
    seed_week = 1          # week at which the index case(s) I0 are introduced
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

  # Working state at the start of week 1. The index case is NOT present yet (it is
  # introduced at seed_week below), so the town is fully susceptible minus baseline
  # immunity. E0 is normally 0 for MAYV but is accepted for generality.
  S_now <- pmax0(N - E0 - R_init_prop * N)
  E_now <- E0
  I_now <- rep(0, A)
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

    # ---- (seed) Introduce the index case(s) at the wet-season onset
    if (t == seed_week) {
      I_now <- I_now + I0
      S_now <- pmax0(S_now - I0)
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
# 2. Vaccine programme parameters
# ============================================================
# Eligible population: adults 18-59 = age groups 4..8 (same as the CHIKV model).
age_group  <- as.character(age_df$age_group)
target_age <- rep(0, A)
target_age[c(4, 5, 6, 7, 8)] <- 1
target_pop <- sum(N[target_age == 1])

total_coverage        <- 0.30   # 30% of the eligible 18-59 population
weekly_delivery_speed <- 0.10   # 10% of supply / week => ~10-week full rollout
immun_delay           <- 2      # 2 weeks from dose to immunity (IXCHIQ-like)

# HYPOTHETICAL MAYV vaccine efficacy. No licensed MAYV product exists. The vaccine
# is modelled as DISEASE-BLOCKING ONLY in every scenario: it reduces the symptomatic
# fraction (VE_block) but never blocks infection (VE_inf = 0), so it never alters
# transmission or the total number of infections. The central disease-blocking
# efficacy is 50%; 25% and 75% are used in the efficacy sensitivity analysis.
mayv_vacc_efficacy    <- 0.50            # central disease-blocking VE_block
mayv_vacc_eff_sens    <- c(0.25, 0.50, 0.75)   # efficacy sensitivity levels

# Pre-outbreak campaign: begin well before the wet-season seed (2025-W26 = index 3,
# vs seed_week = 26), so the immunised fraction is in place before any local
# transmission. This is the maximum-achievable ("ceiling") timing.
# week_to_index() is defined in ca_common.R.
start_pre <- week_to_index(2025, 26)   # 2025-W26 -> index 3

cat("Eligible 18-59 population:", round(target_pop), "\n")
cat("Target doses (30% coverage):", round(target_pop * total_coverage), "\n")
cat(sprintf("MAYV vaccine: DISEASE-BLOCKING ONLY, central VE_block = %.0f%% (sensitivity %s%%)\n",
            100 * mayv_vacc_efficacy, paste(100 * mayv_vacc_eff_sens, collapse = "/")))
cat(sprintf("Pre-outbreak campaign starts at index %d (2025-W26); wet-season seed at %d\n",
            start_pre, seed_week))
stopifnot(start_pre >= 1, start_pre <= T_weeks, start_pre < seed_week)

# ============================================================
# 3. Burden parameters (BORROWED CHIKV severity -- alphavirus proxy)
# ============================================================
# No MAYV-specific severity data exists, so hospitalisations/deaths use the CHIKV
# disease-progression Beta hyperparameters (Hyolim Table S4): means feed the point
# estimates, the alpha/beta feed the Monte Carlo draws (section 7). These are a
# CHIKV-equivalent UPPER BOUND on MAYV severity.
# Read the Beta(alpha, beta) hyperparameters and unpack the severity parameters
# (ps_*, hosp_*, cfr_*, age_to_band, cfr_vec) into the global environment. The
# loader and burden() live in ca_common.R (shared with CHIKV_ca_vacc.R). These are
# the BORROWED CHIKV severity values -- a CHIKV-equivalent upper bound for MAYV.
invisible(list2env(load_burden_params(A), globalenv()))

# ============================================================
# 4. Scenario runner (DISEASE-BLOCKING ONLY, VE_inf = 0 in every scenario)
# ============================================================
# R0 enters through the seasonal beta exactly as in the pre-vacc model:
#   base_beta = R0 * gamma * season   (season has mean 1 over the window).
# season_use / prop_symp_use let the Monte Carlo swap in draws.
run_scenario <- function(R0, coverage, start_week, VE_inf, VE_block,
                         season_use = season, prop_symp_use = prop_symp) {
  base_beta <- R0 * gamma * season_use
  seirv_vaccinated_MAYV(
    T_weeks = T_weeks, A = A, N = N, R_init_prop = R_init_prop, I0 = I0,
    base_beta = base_beta, sigma = sigma, gamma = gamma, rho = rho,
    target_age = target_age, total_coverage = coverage,
    weekly_delivery_speed = weekly_delivery_speed, delay = start_week,
    VE_inf = VE_inf, VE_block = VE_block, immun_delay = immun_delay,
    prop_symp = prop_symp_use, E0 = E0, seed_week = seed_week
  )
}

# ============================================================
# 5. Scenario grid: outbreak size x disease-blocking efficacy
# ============================================================
# Layer 1 -- OUTBREAK SIZE (R0), FIXED at each layer value (the 1.1-1.3 range is the
#   separate R0 sensitivity in section 8, NOT part of the 95% UI):
#   * Current R0 = 1.2  (Caicedo 2021, outside Amazon)
#   * Future  R0 = 3.51 (Dodero-Rojas 2020, urban Ae. aegypti UPPER bound)
# Layer 2 -- VACCINE: DISEASE-BLOCKING ONLY (VE_inf = 0). Central VE_block = 50%;
#   25% and 75% are the efficacy sensitivity levels. Infection-blocking is 0% in
#   EVERY scenario, so infections never differ from the no-vaccine baseline.
#
# This file is the ENGINE: it builds the scenarios + Monte-Carlo objects in memory
# (perm, perm_mc, av_mc, pe, tbl_total, tbl_averted, sens, ...). All FIGURES and
# EXCEL writing live in MAYV_outputs.R -- run that after this script.

size_R0 <- c("Current outbreak (R0 = 1.2)" = R0_central,
             "Future outbreak (R0 = 3.51)" = 3.51)
eff_lvl <- setNames(mayv_vacc_eff_sens, paste0(100 * mayv_vacc_eff_sens, "% efficacy"))  # 25/50/75
eff_central_lab <- paste0(100 * mayv_vacc_efficacy, "% efficacy")                        # "50% efficacy"

# Vaccine designs: a no-vaccine reference + one disease-blocking arm per efficacy.
mech <- data.frame(eff = "No vaccine", VE_block = 0, coverage = 0, stringsAsFactors = FALSE)
for (el in names(eff_lvl))
  mech <- rbind(mech, data.frame(eff = el, VE_block = eff_lvl[[el]], coverage = total_coverage))
perm <- do.call(rbind, lapply(names(size_R0), function(s) cbind(size = s, mech, stringsAsFactors = FALSE)))
perm$key <- paste(perm$size, perm$eff, sep = " | ")
np <- nrow(perm)   # 2 sizes x (1 no-vaccine + 3 efficacies) = 8

# One scenario: disease-blocking only (VE_inf = 0). season_use / prop_symp_use let
# the Monte Carlo swap in draws.
run_perm <- function(R0, coverage, VE_block, season_use = season, prop_symp_use = prop_symp)
  run_scenario(R0, coverage, start_pre, VE_inf = 0, VE_block = VE_block,
               season_use = season_use, prop_symp_use = prop_symp_use)

outcomes     <- c("infections", "symptomatic", "hospitalisations", "deaths")
outcome_labs <- c("Infections", "Symptomatic", "Hospitalisations", "Deaths")

# ---- Point estimates (R0 fixed at each layer value, central inputs) ----
pe <- t(sapply(seq_len(np), function(j)
  burden(run_perm(size_R0[[perm$size[j]]], perm$coverage[j], perm$VE_block[j]))))
rownames(pe) <- perm$key

# ---- Monte Carlo 95% UI: R0 FIXED; propagate seasonal shape + severity (prop_symp,
#      hosp rate, CFR), paired per draw across all scenarios. ----
shape_samples <- as.matrix(read.csv("chikv_ca_beta_shape_samples.csv"))
n_post  <- nrow(shape_samples)
n_draws <- 1000
set.seed(2027)
shape_draws <- sample.int(n_post, n_draws, replace = TRUE)
perm_mc <- lapply(seq_len(np), function(x) matrix(NA_real_, n_draws, 4, dimnames = list(NULL, outcomes)))

cat("\nRunning", n_draws, "MC draws for", np, "scenarios...\n")
for (i in 1:n_draws) {
  season_i  <- shape_samples[shape_draws[i], ]
  ps_i      <- rbeta(1, ps_a, ps_b)
  hosp_i    <- rbeta(1, hosp_a, hosp_b)
  cfr_h_i   <- rbeta(9, cfr_hosp_a, cfr_hosp_b)
  cfr_n_i   <- rbeta(9, cfr_nonh_a, cfr_nonh_b)
  cfr_vec_i <- (hosp_i * cfr_h_i + (1 - hosp_i) * cfr_n_i)[age_to_band]
  for (j in seq_len(np)) {
    out <- run_perm(size_R0[[perm$size[j]]], perm$coverage[j], perm$VE_block[j],
                    season_use = season_i, prop_symp_use = ps_i)
    perm_mc[[j]][i, ] <- burden(out, hosp_i, cfr_vec_i)
  }
  if (i %% 100 == 0) cat("  ", i, "/", n_draws, "\n")
}

# ============================================================
# 6. Tables: total burden + averted, all disease-blocking scenarios
# ============================================================
# Averted = same-R0 no-vaccine baseline - scenario, PAIRED per draw.
base_idx <- sapply(seq_len(np), function(j) which(perm$size == perm$size[j] & perm$eff == "No vaccine"))
vac_rows <- which(perm$eff != "No vaccine")
av_pe <- t(sapply(vac_rows, function(j) pe[base_idx[j], ] - pe[j, ]))
av_mc <- lapply(vac_rows, function(j) perm_mc[[base_idx[j]]] - perm_mc[[j]])

dgt <- c(infections = 0, symptomatic = 0, hospitalisations = 1, deaths = 2)
burden_tbl <- function(rows, mats) data.frame(
  `Outbreak size`             = perm$size[rows],
  `Vaccine (disease-blocking)`= perm$eff[rows],
  Infections       = sapply(seq_along(mats), function(k) fmtq(mats[[k]][, "infections"],       dgt[["infections"]])),
  Symptomatic      = sapply(seq_along(mats), function(k) fmtq(mats[[k]][, "symptomatic"],      dgt[["symptomatic"]])),
  Hospitalisations = sapply(seq_along(mats), function(k) fmtq(mats[[k]][, "hospitalisations"], dgt[["hospitalisations"]])),
  Deaths           = sapply(seq_along(mats), function(k) fmtq(mats[[k]][, "deaths"],           dgt[["deaths"]])),
  check.names = FALSE, row.names = NULL)

tbl_total   <- burden_tbl(seq_len(np), perm_mc)
tbl_averted <- burden_tbl(vac_rows,   av_mc)

old_width <- getOption("width"); options(width = 240)
cat("\n=== TOTAL burden (median, 95% UI) ===\n");                 print(tbl_total,   row.names = FALSE)
cat("\n=== Burden AVERTED vs no vaccine (median, 95% UI) ===\n"); print(tbl_averted, row.names = FALSE)
options(width = old_width)

# ============================================================
# 7. R0 SENSITIVITY (Caicedo 1.1 - 1.3) at the central 50% disease-blocking vaccine
# ============================================================
# Separate from the 95% UI: deterministic point estimates showing how the
# Current-outbreak results move with the assumed R0. The Future outbreak (R0 = 3.51)
# is a single Dodero-Rojas upper-bound value with no range, so it is not swept here.
sens_R0 <- c(1.1, 1.2, 1.3)
rnd <- function(v) mapply(function(x, d) round(x, d), v[outcomes], dgt[outcomes])
sens <- do.call(rbind, lapply(sens_R0, function(r0) {
  b0 <- burden(run_perm(r0, 0, 0))                               # no vaccine
  bv <- burden(run_perm(r0, total_coverage, mayv_vacc_efficacy)) # 50% disease-blocking
  data.frame(R0 = r0, Outcome = outcome_labs,
             `No vaccine`     = rnd(b0),
             `50% disease-blocking` = rnd(bv),
             Averted          = rnd(b0 - bv),
             `% averted`      = round(100 * (b0 - bv)[outcomes] / b0[outcomes], 1),
             check.names = FALSE, row.names = NULL)
}))
cat("\n=== R0 SENSITIVITY: burden at R0 = 1.1 / 1.2 / 1.3, 50% disease-blocking (point estimates) ===\n")
print(sens, row.names = FALSE)

cat("\nEngine ready. Run MAYV_outputs.R for figures + Excel workbooks.\n")
