# ============================================================
# Caldas Novas MAYV -- uncertainty propagation for the forward outbreak
# (self-contained).
#
# WHAT THIS IS. Parameter PRIORS are propagated through a forward SEIR: each Latin
# hypercube draw re-runs the simulation, so the reported / infection / attack-rate
# bands carry natural-history, reporting, R0 and prior-immunity uncertainty. There is
# NO fitting -- Caldas Novas has no observed MAYV outbreak -- so these are
# PRIOR-PREDICTIVE bands, not posteriors. The seasonal transmission SHAPE is a HYBRID
# (caldas_hybrid_season.rds, built by MAYV_build_hybrid_envelope.R): the fitted CHIKV
# beta_t for the rise/peak + a climatological dry-season tail.
# Window: 52 weeks, 2025-W24 -> 2026-W22.
#
# SAMPLED INPUTS (priors, from model_calibration.xlsx MAYV rows unless noted):
#     gamma   ~ Normal(rate)         recovery rate  (7 d central, 5-10 d range)
#     sigma   = 1 / Normal(period)   intrinsic incubation 3.0 d, 95% CrI 2.2-3.8 (Caicedo 2021)
#     rho     ~ Beta(20, 60)         reporting rate (mean 0.25)      [hardcoded]
#     prop_symp ~ Beta(35.84, 32.56) symptomatic fraction (med 0.524)[hardcoded]
#     immune  ~ Lognormal            prior-immune FRACTION, FLAT across ages;
#                                    Lima 2021 Central-West 8% (95% CI 3-18%)
#
# HOW IMMUNITY IS APPLIED (Route B: flat seroprevalence, NOT the CHIKV catalytic
# FOI*age model). Each draw takes a single fraction p and sets R_init_prop=(p,...,p),
# i.e. immune_a = p * N_a. p ~ Lognormal calibrated to Lima et al. 2021's pooled
# CENTRAL-WEST exposure rate (8%, 95% CI 3-18%, I2=98%): the lognormal reproduces that
# CI and its right skew. FLAT is used deliberately -- although MAYV exposure is really
# occupational (working-age adults enter the forest), transmission here is HOMOGENEOUS
# (foi = beta*sum(I)/N_total), so only the OVERALL p affects total infections/attack
# rate; re-slicing a fixed p across ages changes only the AGE DISTRIBUTION of infections,
# which matters solely for age-stratified burden (add an exposure-weight vector there).
# p sets S(0), so it widens the ATTACK-RATE / total-infection band.
#
# R0 = wet-season PEAK R_eff, and it is FIXED PER SCENARIO -- NOT sampled. base_beta =
# R0 * gamma * season, so R_eff(t) = R0*season(t)*S/N is INDEPENDENT of gamma
# (gamma/sigma move only PEAK TIMING & HEIGHT; size is driven by R0 & immunity).
# SEASONAL-PEAK scaling (envelope rescaled so max = 1) makes the cited R0 the wet-season
# PEAK R_eff -- the honest, load-bearing quantity: it avoids the mean-1 framing that hides
# a higher true peak. NB R0(t) = R0*season(t) is just beta_t/gamma with an imposed seasonal
# shape (standard seasonal forcing), NOT a new formula.
#
# WHY R0 IS FIXED, NOT SAMPLED. The published MAYV R0 figures are POINT ESTIMATES FROM
# DIFFERENT PLACES AND DECADES, each with its own interval -- Caicedo et al. 2021 give
# 1.11 (French Guiana 1960), 3.47 (Brazil 1966), 2.1-2.9 (Amazon basin), 1.1-1.3 (outside
# the Amazon); Dodero-Rojas et al. 2020 give 1.18-3.51 as LOWER/UPPER LIMITS. The spread
# BETWEEN them is between-setting heterogeneity (vectors, land use, host contact), not
# uncertainty about one municipality. Sampling R0 across that span therefore does not
# produce a meaningful 95% UI for Caldas Novas: because outbreak size is a steep convex
# function of R0, the draws split into "fizzle" and "explode" regimes, the output becomes
# bimodal, and the reported central estimate detaches from the cited central R0. (Under the
# old truncated-Lognormal[1.18,3.51] prior the stated median was 2.03, which alone yields
# ~25 infections, yet the reported conditional median was ~6,650 -- i.e. the number
# described R0 ~= 2.83, the take-off tail.) R0 is a SCENARIO-DEFINING variable here, so we
# fix it and propagate the parameters that ARE genuine single-setting uncertainty.
#
# THE TWO SCENARIOS (R0_SCENARIO, peak R_eff):
#   low  = 1.20  Caicedo et al. 2021 outside-Amazon-basin (1.1-1.3) -- current ecology.
#                Read as a PEAK (not annual-mean) value: this is the conservative reading,
#                since an annual-mean 1.2 would imply a peak of ~3+ and a far larger
#                outbreak. A single introduction does not self-sustain at this R0.
#   high = 2.04  Urban Aedes-borne transmission. This is the GEOMETRIC MEAN of Dodero-Rojas
#                et al. 2020's [1.18, 3.51] -- the only MAYV-SPECIFIC R0 range published --
#                i.e. the central value of that range: sqrt(1.18 * 3.51) = 2.035. Geometric
#                rather than arithmetic because R0 lives on a ratio scale.
#                Preferred over the lower bound of Caicedo's Amazon-basin range (2.1-2.9):
#                Caicedo derive R0 from age-stratified seroprevalence via catalytic models,
#                so 2.1-2.9 is an ENDEMIC-AVERAGE reproduction number, not a seasonal peak.
#                Importing one of its bounds as a PEAK R0 would silently mix two different
#                quantities; the geometric centre of an explicitly-stated R0 range does not.
#                Context: the peak R0 fitted for CHIKV in THIS municipality is 2.67
#                (95% UI 2.39-3.13), and the Belterra 1977-78 growth rate implies ~1.9-2.5,
#                so 2.04 is conservative against both. It is a COUNTERFACTUAL: it assumes
#                MAYV acquires Ae. aegypti transmission competence it is not currently
#                known to possess, and is not a prediction.
#                NB R0 and the latent period are partly INTERCHANGEABLE in setting outbreak
#                size, because both control how many generations fit inside the seasonal
#                window before the dry-season tail collapses transmission. Shortening the
#                latent period from 12 d to 3 d multiplies outbreak size ~50x at fixed R0,
#                which is why moving to Caicedo's 3.0 d required re-anchoring R0 downward.
# BECAUSE R0 IS FIXED, outbreak size is unimodal within a scenario and there is no
# fizzle/take-off split, so the engine reports over ALL draws (no conditioning).
# There is still NO MAYV outbreak to fit, so these remain PRIOR-PREDICTIVE bands.
#
# WINDOW: 2025-W24 -> 2026-W22 (52 epi weeks). Extended past the CHIKV fit window (which
# ends 2026-W22) so the dry-season tail sits inside the window and every outbreak resolves.
# 2025 carries an epi-week 53, so 2025-W24..W53 = 30 wks and 2026-W01..W22 = 22 wks.
#
# FIXED (NOT sampled):
#   * seasonal envelope = HYBRID (caldas_hybrid_season.rds, built by
#     MAYV_build_hybrid_envelope.R): fitted CHIKV beta_t for the data-constrained rise
#     and peak (2025-W24..2026-W10) spliced onto the CHIRPS climatological dry season for
#     the tail (2026-W10..W22). CHIKV and MAYV share the Aedes vector, so the fitted CHIKV
#     signal is the best proxy for the rise; the dry-season tail (which the beta lacks --
#     its post-outbreak tail is an unconstrained spline artefact) is what lets the outbreak
#     collapse. A DETERMINISTIC covariate (fixed weekly shape), so it does NOT enter the LHS.
#   * single index case seeded at the WINDOW OPEN (2025-W24, seed_week = 1), matching the
#     CHIKV engine (infection present from t = 1) so the MAYV outbreak rises with CHIKV.
# ============================================================
setwd("/Users/chloelee/Documents/R/summer_project")
suppressMessages({library(readxl); library(dplyr); library(tidyr); library(ggplot2)})

# ------------------------------------------------------------
# 1. SEIR forward simulator (seeded at seed_week)
# ------------------------------------------------------------
seir_baseline_MAYV <- function(
    T_weeks, A, N, R_init_prop, I0, base_beta, sigma, gamma, rho, prop_symp,
    sub_steps = 7, E0 = rep(0, A), seed_week = 1) {
  pmax0 <- function(x) pmax(0, x); N_total <- sum(N); dt <- 1/sub_steps
  S <- E <- I <- R <- matrix(0, A, T_weeks)
  new_infections <- new_symptomatic <- matrix(0, A, T_weeks)
  S_now <- pmax0(N - E0 - R_init_prop*N); E_now <- E0; I_now <- rep(0, A); R_now <- R_init_prop*N
  for (t in 1:T_weeks) {
    new_I_week <- rep(0, A); beta_t <- base_beta[t]
    if (t == seed_week) { I_now <- I_now + I0; S_now <- pmax0(S_now - I0) }
    for (k in 1:sub_steps) {
      foi <- beta_t * sum(I_now)/N_total
      new_E <- foi*S_now*dt; new_I <- sigma*E_now*dt; new_R <- gamma*I_now*dt
      S_now <- pmax0(S_now - new_E); E_now <- pmax0(E_now + new_E - new_I)
      I_now <- pmax0(I_now + new_I - new_R); R_now <- pmax0(R_now + new_R)
      new_I_week <- new_I_week + new_I
    }
    S[,t]<-S_now; E[,t]<-E_now; I[,t]<-I_now; R[,t]<-R_now
    new_infections[,t] <- new_I_week; new_symptomatic[,t] <- prop_symp*new_I_week
  }
  list(new_infections=new_infections, new_symptomatic=new_symptomatic,
       new_reported=rho*new_symptomatic)
}

# ------------------------------------------------------------
# 2. Population (Caldas Novas), grown 2022 -> 2025
#    (immunity is now a SAMPLED input, so susceptible_pop / I0 are computed per draw)
# ------------------------------------------------------------
age_df <- read_excel("population.xlsx", sheet = "prop_immune")
age_df <- as.data.frame(age_df[tolower(age_df$municipality) == "caldas novas", ])
stopifnot(nrow(age_df) > 0)
A <- nrow(age_df)
pop_2022_total <- 98622; pop_2025_total <- 106820
growth_r <- log(pop_2025_total/pop_2022_total)/3
N        <- age_df$pop_num * exp(growth_r*3)

T_weeks <- 52   # 2025-W24 -> 2026-W22, aligned with the CHIKV chain

# ------------------------------------------------------------
# 3. HYBRID seasonal envelope: CHIKV beta (rise/peak) + climatological dry-season tail
# ------------------------------------------------------------
# caldas_hybrid_season.rds (built by MAYV_build_hybrid_envelope.R): 52-week mean-1 envelope, 2025-W24 ->
# 2026-W22. It uses the FITTED CHIKV beta_t where it is data-constrained (the rise and
# peak, 2025-W24 .. 2026-W10) and splices on the CHIRPS climatological DRY SEASON for the
# tail (2026-W10 .. W22). Rationale: the CHIKV beta's post-outbreak tail is a spline
# artefact (flat ~0.9, unconstrained -- no CHIKV cases there), which lacks a real dry
# season; against MAYV's high R0 that flat tail keeps R_eff ~ 1 and the outbreak never
# resolves inside the window. The rainfall dry season (deep trough, min ~0.04) crashes
# R_eff after the peak so the outbreak collapses -- for the same reason CHIKV's own
# outbreak does. So each source is used only where it is informative.
season_mean1 <- readRDS("caldas_hybrid_season.rds")
stopifnot(length(season_mean1) == T_weeks, abs(mean(season_mean1) - 1) < 1e-6)
# R0 INTERPRETATION. r0_is_peak = TRUE: R0 is the SEASONAL-PEAK R_eff (rescale envelope
# so max = 1; yearly-avg R0 = R0*mean = ~0.58*R0). FALSE: R0 is the ANNUAL-MEAN (envelope
# left at mean 1; wet peak reaches R0*max ~ 1.7x). Seed/wet-band logic below always uses
# the mean-1 envelope so it is unaffected by this choice. NB the CHIKV beta_t envelope is
# FLATTER than the rainfall one it replaces (max 1.73 vs 2.44, min 0.58 vs 0.04).
r0_is_peak <- TRUE
season <- if (r0_is_peak) season_mean1 / max(season_mean1) else season_mean1

# ------------------------------------------------------------
# 4. Seed (single introduction at the WINDOW OPEN, matching the CHIKV chain)
# ------------------------------------------------------------
# The CHIKV engine has infection present from t = 1 (2025-W24). To make the MAYV
# outbreak RISE AT THE SAME TIME as CHIKV (both driven by the shared beta_t envelope),
# we introduce the index case at the window open too. The old "wet-season onset"
# rule (first week beta > its mean, ~2025-W43) was a leftover from the CHIRPS-rainfall
# design; against the flatter CHIKV beta it fired ~11 weeks after the envelope starts
# rising and, being R0-independent, seeded high-R0 outbreaks ~19 weeks too late so they
# ran off the end of the 52-week window. Seeding at t=1 fixes both: the beta trough
# over 2025-W24..W35 keeps early growth slow (so a W40 vaccination is still pre-surge),
# then the outbreak surges as beta peaks (~W50) and resolves inside the window.
I0_total  <- if (exists("MAYV_I0"))   MAYV_I0   else 1     # seed size (infectious persons); override for experiments
E0        <- rep(0, A)
seed_week <- if (exists("MAYV_SEED")) MAYV_SEED else 1     # seed week index; 1 = window open (2025-W24), 17 = 2025-W40

# ------------------------------------------------------------
# 5. Priors / samplers
# ------------------------------------------------------------
cal <- as.data.frame(read_excel("model_calibration.xlsx", sheet = 1))
cal <- cal[cal$Group == "MAYV" & !is.na(cal$Median), ]
row_for <- function(k) cal[grepl(k, cal$Parameter, ignore.case = TRUE), ][1, ]
sd_of   <- function(r) (r[["95% UI upper"]] - r[["95% UI lower"]]) / (2*1.96)

gr <- row_for("gamma"); sr <- row_for("sigma")
g_m <- gr$Median; g_sd <- sd_of(gr)            # gamma is a RATE (per week)
p_m <- sr$Median; p_sd <- sd_of(sr)            # sigma stored as a PERIOD (weeks) -> invert

# LATENT-PERIOD SCENARIO. Base case = the workbook value: Caicedo et al. 2021's intrinsic
# incubation period, 3.0 d (95% CrI 2.4-4.1) -- LOGNORMAL, fitted so the CrI endpoints are
# the 2.5th/97.5th percentiles. Lognormal suits a duration bounded below by zero and is
# right-skewed. NB the fitted MEDIAN is 3.14 d, not 3.0: no two-parameter lognormal can
# reproduce a mean of 3.0 AND the CrI 2.4-4.1, because that interval is not symmetric about
# 3.0 (its midpoint is 3.25). We preserve the INTERVAL, which carries the uncertainty, and
# accept a 4.6% shift in the central value; forcing the median to 3.0 would instead shrink
# the CrI to 2.30-3.92. Using Caicedo for BOTH the latent period and R0 keeps the natural
# history and the reproduction number internally consistent (they were estimated together).
#   NB Caicedo also report a distribution SD of 0.3 d (95% CrI 0.0-2.1). That is BETWEEN-
#   INDIVIDUAL variation in incubation time, which the SEIR's exponential E->I waiting time
#   already represents. The LHS samples uncertainty in the MEAN, so the CrI on the mean is
#   the correct input -- same convention as every other row in the workbook.
# Sensitivity: MAYV_LATENT <- "long12d" restores the previous 12-day assumption. That is a
# STRUCTURAL sensitivity, not a wider prior: it changes the generation time from ~1.43 to
# ~2.71 wk, so far fewer generations fit inside the supercritical season before the dry-
# season tail collapses transmission, and the outbreak is ~50x smaller at the same R0.
# Latent period and R0 are therefore partly interchangeable in setting outbreak size.
if (!exists("MAYV_LATENT")) MAYV_LATENT <- "caicedo3d"    # "caicedo3d" | "long12d"
stopifnot(MAYV_LATENT %in% c("caicedo3d", "long12d"))
# Lognormal on the PERIOD, fitted to the workbook's 95% UI endpoints (weeks).
lat_lo <- sr$`95% UI lower`; lat_hi <- sr$`95% UI upper`
if (MAYV_LATENT == "long12d") { lat_lo <- 11.31/7; lat_hi <- 12.69/7 }  # pre-2026-07 12 d assumption
lat_ml <- (log(lat_lo) + log(lat_hi)) / 2
lat_sl <- (log(lat_hi) - log(lat_lo)) / (2 * 1.96)
p_m <- exp(lat_ml)                              # median period, weeks
cat(sprintf("Latent period '%s': lognormal(meanlog %.6f, sdlog %.6f) -> median %.4f wk (%.2f d), 95%% UI %.2f-%.2f d\n",
            MAYV_LATENT, lat_ml, lat_sl, p_m, 7*p_m, 7*lat_lo, 7*lat_hi))
rab <- c(a = 20, b = 60)                        # rho ~ Beta(20,60), mean 0.25
# SYMPTOMATIC FRACTION -- MAYV-specific, replacing the borrowed CHIKV Beta(35.84, 32.56).
# Bounded by the only two usable MAYV denominators:
#   lower 48/71 = 0.6761  Belterra 1978 (71 seropositive, 23 of them asymptomatic)
#   upper 33/36 = 0.9167  Azevedo narrative review
# Lognormal fitted so those are the 2.5th/97.5th percentiles: right-skewed, so most mass
# sits toward the LOWER end. That is the conservative direction, because most MAYV case
# series ascertain only febrile presentations and so over-count the symptomatic fraction.
# TRUNCATED at 1 -- a proportion cannot exceed it, and the untruncated fit puts ~0.1% of
# draws above 1 (about 1 draw in 1000).
ps_ml  <- (log(48/71) + log(33/36)) / 2
ps_sl  <- (log(33/36) - log(48/71)) / (2 * 1.96)
ps_cap <- plnorm(1, ps_ml, ps_sl)               # CDF mass below 1, for the truncation

# R0 = wet-season PEAK R_eff (r0_is_peak = TRUE above), FIXED PER SCENARIO (see the
# header for why it is not sampled). Two scenarios only; each is one transmission regime.
if (!exists("R0_SCENARIO")) R0_SCENARIO <- "high"   # "high" = 2.1-2.9 | "low" = 1.1-1.3
# R0 is SAMPLED from the scenario's range, lognormal with the endpoints as the 2.5th/97.5th
# percentiles. BOTH ranges are Caicedo et al. 2021, so the two scenarios are one quantity
# measured in two ecological zones rather than a mix of sources -- no geometric mean taken
# across a modelling paper's limits, and no single point estimate asserted.
R0_RANGE <- list(low  = c(1.1, 1.3),    # outside the Amazon basin -- current ecology
                 high = c(2.1, 2.9))    # Amazon basin, applied to Goias as a PEAK R0
stopifnot(R0_SCENARIO %in% names(R0_RANGE))
R0_LO <- R0_RANGE[[R0_SCENARIO]][1]; R0_HI <- R0_RANGE[[R0_SCENARIO]][2]
R0_meanlog <- (log(R0_LO) + log(R0_HI)) / 2
R0_sdlog   <- (log(R0_HI) - log(R0_LO)) / (2 * 1.96)
R0_median  <- exp(R0_meanlog)                   # geometric centre of the range
R0_VALUE   <- R0_median                         # name kept for downstream compatibility

# Prior immunity (FLAT fraction), Lima et al. 2021 pooled CENTRAL-WEST exposure rate
# 8% (95% CI 3-18%, I2=98%). Right-skewed -> Lognormal calibrated to the CI endpoints:
# it reproduces [0.03,0.18] exactly (median 7.3%, mean ~8%). Swap this block for a
# Beta/other if the evidence changes.
imm_lo95 <- 0.03; imm_hi95 <- 0.18
imm_meanlog <- (log(imm_lo95) + log(imm_hi95)) / 2         # median = exp(meanlog) ~ 0.073
imm_sdlog   <- (log(imm_hi95) - log(imm_lo95)) / (2*1.96)
imm_base    <- 0.08                                        # Lima central estimate, for the baseline run

cat(sprintf("SAMPLED priors:\n  gamma     ~ N(%.3f, %.4f) /wk\n  latent    ~ logN(%.4f, %.4f) wk  -> median %.2f d, 95%% %.2f-%.2f d\n  rho       ~ Beta(%d, %d)\n  prop_symp ~ logN(%.4f, %.4f) trunc at 1 -> median %.3f, 95%% %.3f-%.3f\n  immune    ~ logN(med %.3f, 95%% [%.2f, %.2f])\n",
            g_m, g_sd,
            lat_ml, lat_sl, 7*exp(lat_ml), 7*qlnorm(.025, lat_ml, lat_sl), 7*qlnorm(.975, lat_ml, lat_sl),
            rab["a"], rab["b"],
            ps_ml, ps_sl, exp(ps_ml), qlnorm(.025, ps_ml, ps_sl), qlnorm(.975, ps_ml, ps_sl),
            exp(imm_meanlog), imm_lo95, imm_hi95))
cat(sprintf("R0 ~ lognormal on [%.1f, %.1f] (%s scenario, %s): median %.3f\n",
            R0_LO, R0_HI, R0_SCENARIO, if (r0_is_peak) "seasonal-peak R_eff" else "annual-mean", R0_median))

# ------------------------------------------------------------
# Seasonality summary. R0 above is the WET-SEASON PEAK, because the envelope is
# peak-normalised (r0_is_peak). The mean of season(t) is well below 1, so the
# year-average transmission intensity is much lower than the headline R0 -- report both,
# or a reader will read "R0 = 2.04" as an unforced R0 and expect a large attack rate.
# The model therefore runs close to the epidemic threshold for most of the year, which
# is why burden is a steep function of R0 rather than a saturating one.
# ------------------------------------------------------------
season_mean <- mean(season); imm_med <- exp(imm_meanlog)   # lognormal median
cat(sprintf("Seasonality: season(t) peak %.3f, MEAN %.4f, min %.3f; %d of %d weeks above 0.5\n",
            max(season), season_mean, min(season), sum(season > 0.5), length(season)))
cat(sprintf("Implied MEAN R0(t) = %.3f; mean R_eff at S/N = %.3f is %.3f (peak R_eff %.3f); %d of %d weeks with R_eff > 1\n",
            R0_VALUE*season_mean, 1-imm_med, R0_VALUE*season_mean*(1-imm_med),
            R0_VALUE*max(season)*(1-imm_med),
            sum(R0_VALUE*season*(1-imm_med) > 1), length(season)))

# ------------------------------------------------------------
# 6. One forward run -> weekly reported / infections + summary scalars
#    imm = flat prior-immune fraction; sets S(0), susceptible_pop and the seed split.
# ------------------------------------------------------------
run_draw <- function(R0, g, s, r, ps, imm) {
  Rimm <- rep(imm, A)                            # FLAT immune fraction across ages
  sus  <- N * (1 - Rimm)                         # susceptibles by age
  I0i  <- I0_total * sus / sum(sus)              # single seed split by susceptibility
  out  <- seir_baseline_MAYV(T_weeks, A, N, Rimm, I0i, R0*g*season,
                             sigma = s, gamma = g, rho = r, prop_symp = ps,
                             E0 = E0, seed_week = seed_week)
  rep_wk <- colSums(out$new_reported); inf_wk <- colSums(out$new_infections)
  list(rep = rep_wk, inf = inf_wk,
       tot_rep = sum(rep_wk), tot_inf = sum(inf_wk),
       peak_rep_wk = which.max(rep_wk), peak_rep = max(rep_wk),
       attack = 100 * sum(inf_wk) / sum(sus), immune = 100 * imm)
}

# Baseline at median inputs (dashed reference on the plots), at the scenario's MEDIAN R0
PSYMP_MED <- exp(ps_ml)                          # median symptomatic fraction (~0.787)
base <- run_draw(R0_median, g_m, 1/p_m, rab["a"]/sum(rab), PSYMP_MED, imm_base)
cat(sprintf("Baseline R0=%.2f (%s scenario, median inputs, immune %.1f%%): total infections %.0f | total reported %.1f | peak reported wk %d (%.1f) | attack %.2f%%\n",
            R0_median, R0_SCENARIO, base$immune, base$tot_inf, base$tot_rep, base$peak_rep_wk, base$peak_rep, base$attack))
# Reference: the OTHER scenario's R0 under the SAME scaling, for context in the log.
R0_other <- exp(mean(log(R0_RANGE[[setdiff(names(R0_RANGE), R0_SCENARIO)]])))
base_other <- run_draw(R0_other, g_m, 1/p_m, rab["a"]/sum(rab), PSYMP_MED, imm_base)
cat(sprintf("Reference R0=%.2f (other scenario, same scaling): total infections %.0f | total reported %.1f | attack %.3f%%\n",
            R0_other, base_other$tot_inf, base_other$tot_rep, base_other$attack))

# ------------------------------------------------------------
# 7. Latin Hypercube (5 SAMPLED inputs), forward-simulate each.
#    R0 is FIXED per scenario (header), so it is NOT an LHS dimension: the propagated
#    bands below carry natural-history / reporting / symptomatic-fraction / prior-immunity
#    uncertainty ONLY. Every draw is the same transmission regime, so the outbreak-size
#    distribution is unimodal and needs no take-off conditioning.
# ------------------------------------------------------------
set.seed(2024); n <- 1000
lhs_col <- function(n) (sample.int(n) - runif(n)) / n
U   <- sapply(1:6, function(j) lhs_col(n))
gam <- qnorm(U[,1], g_m, g_sd)                  # rate
sig <- 1/qlnorm(U[,2], lat_ml, lat_sl)          # lognormal PERIOD -> rate
rho <- qbeta(U[,3], rab["a"], rab["b"])
psy <- qlnorm(U[,4]*ps_cap, ps_ml, ps_sl)       # lognormal, truncated at 1
imm <- qlnorm(U[,5], imm_meanlog, imm_sdlog)    # flat prior-immune fraction (Lima 2021)
R0v <- qlnorm(U[,6], R0_meanlog, R0_sdlog)      # SAMPLED from the scenario range

rep_mat <- inf_mat <- matrix(NA_real_, n, T_weeks)
tot_rep <- tot_inf <- peak_rep_wk <- peak_rep <- attack <- R0peak <- immune <- rep(NA_real_, n)
cat("Forward-simulating", n, "LHS draws...\n")
for (i in 1:n) {
  d <- tryCatch(run_draw(R0v[i], gam[i], sig[i], rho[i], psy[i], imm[i]), error = function(e) NULL)
  if (is.null(d)) next
  rep_mat[i,] <- d$rep; inf_mat[i,] <- d$inf
  tot_rep[i] <- d$tot_rep; tot_inf[i] <- d$tot_inf
  peak_rep_wk[i] <- d$peak_rep_wk; peak_rep[i] <- d$peak_rep; attack[i] <- d$attack
  immune[i] <- d$immune; R0peak[i] <- R0v[i] * max(season)
  if (i %% 100 == 0) cat("  ", i, "/", n, "\n")
}
# Every forward draw is feasible (no data-consistency filter: there is no observed
# MAYV outbreak to compare against). Keep only numerically finite draws.
ok <- which(is.finite(tot_inf) & is.finite(tot_rep))
cat(sprintf("Kept %d / %d finite draws.\n", length(ok), n))

# ------------------------------------------------------------
# 8. Summaries (baseline point vs propagated median [95% UI])
# ------------------------------------------------------------
q3   <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
band <- function(M) apply(M[ok, , drop = FALSE], 2, quantile, c(.025, .5, .975), na.rm = TRUE)
cat("\n=========== PROPAGATED MAYV RESULTS (", length(ok), " draws, R0 SAMPLED ",
    sprintf("%.1f-%.1f", R0_LO, R0_HI), ") ===========\n", sep = "")
cat("Bands below carry gamma / sigma / rho / prop_symp / prior-immunity AND R0 uncertainty.\n")
cat(sprintf("Prior immunity:      propagated %.1f%% [%.1f%%, %.1f%%]\n",
            q3(immune[ok])[1], q3(immune[ok])[2], q3(immune[ok])[3]))
cat(sprintf("R0 at seasonal peak: SAMPLED %.2f [%.2f, %.2f] (lognormal on the %s range)\n",
            median(R0v[ok]), quantile(R0v[ok],.025), quantile(R0v[ok],.975), R0_SCENARIO))
cat(sprintf("Attack rate (of susceptibles): baseline %.2f%% -> propagated %.2f%% [%.2f%%, %.2f%%]\n",
            base$attack, q3(attack[ok])[1], q3(attack[ok])[2], q3(attack[ok])[3]))
cat(sprintf("Total infections:  baseline %.0f -> propagated %.0f [%.0f, %.0f]\n",
            base$tot_inf, q3(tot_inf[ok])[1], q3(tot_inf[ok])[2], q3(tot_inf[ok])[3]))
cat(sprintf("Total reported:    baseline %.1f -> propagated %.1f [%.1f, %.1f]\n",
            base$tot_rep, q3(tot_rep[ok])[1], q3(tot_rep[ok])[2], q3(tot_rep[ok])[3]))
cat(sprintf("Busiest reported week (count): baseline %.1f -> propagated %.1f [%.1f, %.1f]\n",
            base$peak_rep, q3(peak_rep[ok])[1], q3(peak_rep[ok])[2], q3(peak_rep[ok])[3]))

# ------------------------------------------------------------
# 9. Plots: propagated 95% bands (reported + true infections), with season shading
# ------------------------------------------------------------
wk_num  <- function(idx) ifelse(idx <= 30, idx + 23, idx - 30)   # 2025-W24..W53 | 2026-W01..W22
tick_idx <- c(7, 17, 27, 40, 50)                      # 2025-W30/40/50, 2026-W10/20
x_ticks  <- data.frame(week_index = tick_idx, week = wk_num(tick_idx))
year_break <- 30.5                                    # 2025-W53 (idx 30) | 2026-W01 (idx 31)
wet_start <- which(season_mean1 >= 1)[1]              # above-mean (wet) span (mean-1 envelope)
wet_end   <- tail(which(season_mean1 >= 1), 1)
season_layers <- list(
  annotate("rect", xmin = 0.5,       xmax = wet_start,   ymin = -Inf, ymax = Inf, fill = "#f4ead7", alpha = 0.55),
  annotate("rect", xmin = wet_start, xmax = wet_end,     ymin = -Inf, ymax = Inf, fill = "#cfe6f2", alpha = 0.55),
  annotate("rect", xmin = wet_end,   xmax = T_weeks+0.5, ymin = -Inf, ymax = Inf, fill = "#f4ead7", alpha = 0.55))

weeks <- 1:T_weeks
rb <- band(rep_mat); ib <- band(inf_mat)

p_rep <- ggplot(data.frame(week = weeks, lo = rb[1,], med = rb[2,], hi = rb[3,], base = base$rep)) +
  season_layers +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(week, ymin = lo, ymax = hi), fill = "grey40", alpha = 0.30) +
  geom_line(aes(week, med), colour = "grey25", linewidth = 1) +
  geom_line(aes(week, base), colour = "#d6604d", linewidth = 0.8, linetype = "dashed") +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  labs(x = "Week", y = "Predicted reported MAYV cases",
       title = "Hypothetical MAYV outbreak: propagated 95% band",
       subtitle = sprintf("Median (solid) + 95%% band from gamma/sigma/rho/prop_symp/immunity/R0 priors (R0 ~ %.1f-%.1f); baseline at median inputs (dashed)", R0_LO, R0_HI)) +
  theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5),
                       plot.subtitle = element_text(hjust = 0.5, size = 9),
                       panel.grid.minor = element_blank())
ggsave("MAYV_ca_lhs_reported.png", p_rep, width = 8, height = 4.5, dpi = 120)

p_inf <- ggplot(data.frame(week = weeks, lo = ib[1,], med = ib[2,], hi = ib[3,])) +
  season_layers +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(week, ymin = lo, ymax = hi), fill = "#a8d1e7", alpha = 0.5) +
  geom_line(aes(week, med), colour = "#3182bd", linewidth = 1) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  labs(x = "Week", y = "True MAYV infections (all)",
       title = "Hypothetical MAYV outbreak: true infections, propagated 95% band") +
  theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5),
                       panel.grid.minor = element_blank())
ggsave("MAYV_ca_lhs_infections.png", p_inf, width = 8, height = 4.5, dpi = 120)

# ------------------------------------------------------------
# 10. Save per-draw table + ensemble (mirrors CHIKV_ca_lhs_ensemble.rds shape so a
#     future MAYV vaccine/engine variant can iterate over these draws directly)
# ------------------------------------------------------------
# This script does not source ca_common.R, so DRAWS_DIR is defined here if absent.
if (!exists("DRAWS_DIR")) DRAWS_DIR <- "lhs_draws"
dir.create(DRAWS_DIR, showWarnings = FALSE)
write.csv(data.frame(draw = 1:n, R0 = R0v, gamma = gam, sigma = sig, rho = rho, prop_symp = psy,
                     immune_frac = imm, total_infections = tot_inf, total_reported = tot_rep,
                     peak_reported_wk = peak_rep_wk, attack_pct = attack,
                     finite = (seq_len(n) %in% ok)),
          file.path(DRAWS_DIR, "MAYV_ca_lhs_draws.csv"), row.names = FALSE)

mayv_lhs_ensemble <- list(
  rep = rep_mat[ok, , drop = FALSE], inf = inf_mat[ok, , drop = FALSE],
  R0 = R0v[ok], gamma = gam[ok], sigma = sig[ok], rho = rho[ok], prop_symp = psy[ok],
  immune_frac = imm[ok],                          # per-draw FLAT immune fraction (Rimm = rep(imm, A))
  base_R0 = R0_median, R0_scenario = R0_SCENARIO, r0_is_peak = r0_is_peak,
  # R0 is now SAMPLED. R0_fixed is kept (= the scenario median) so downstream labels and
  # the OWSA base case keep working; R0_sampled/R0_lo/R0_hi describe the actual prior.
  R0_fixed = R0_median, R0_sampled = TRUE, R0_lo = R0_LO, R0_hi = R0_HI,
  R0_meanlog = R0_meanlog, R0_sdlog = R0_sdlog,
  lat_meanlog = lat_ml, lat_sdlog = lat_sl, MAYV_LATENT = MAYV_LATENT,
  ps_meanlog = ps_ml, ps_sdlog = ps_sl,
  base_gamma = g_m, base_sigma = 1/p_m, base_prop_symp = PSYMP_MED,
  base_immune = imm_base,
  season = season, N = N, A = A, age_df = age_df,
  I0_total = I0_total, E0 = E0, seed_week = seed_week,
  T_weeks = T_weeks, weeks = weeks, x_ticks = x_ticks, year_break = year_break)
saveRDS(mayv_lhs_ensemble, "MAYV_ca_lhs_ensemble.rds")
# Scenario-tagged copy, so both R0 scenarios can coexist on disk for comparison.
saveRDS(mayv_lhs_ensemble, sprintf("MAYV_ca_lhs_ensemble_%s.rds", R0_SCENARIO))

cat(sprintf("\nSaved MAYV_ca_lhs_reported.png, MAYV_ca_lhs_infections.png, MAYV_ca_lhs_draws.csv,\n  MAYV_ca_lhs_ensemble.rds and MAYV_ca_lhs_ensemble_%s.rds (R0 = %.2f)\n",
            R0_SCENARIO, R0_VALUE))
