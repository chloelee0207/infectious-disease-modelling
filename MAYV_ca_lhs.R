# ============================================================
# Caldas Novas MAYV -- uncertainty propagation for the forward outbreak
# (self-contained; the standalone LHS twin of MAYV_ca_pre_vacc.R).
#
# WHAT THIS IS. MAYV_ca_pre_vacc.R runs the outbreak at FIXED natural-history /
# reporting values and builds its band only from the borrowed CHIKV seasonal-shape
# posterior. This script instead PROPAGATES the parameter PRIORS through the same
# forward SEIR: for each Latin-Hypercube draw it re-runs the simulation, so the
# reported / infection / attack-rate bands carry the natural-history + reporting
# + R0 + prior-immunity uncertainty. There is NO fitting -- there is no observed
# MAYV outbreak in Caldas Novas to fit, so these are PRIOR-predictive bands, not
# posteriors. This is MAYV's OWN pipeline; the CHIKV LHS/engine are left untouched.
#
# SAMPLED INPUTS (priors, from model_calibration.xlsx MAYV rows unless noted):
#     gamma   ~ Normal(rate)         recovery rate  (7 d central, 5-10 d range)
#     sigma   = 1 / Normal(period)   intrinsic incubation ~12 d (tight)
#     rho     ~ Beta(20, 60)         reporting rate (mean 0.25)      [hardcoded]
#     prop_symp ~ Beta(35.84, 32.56) symptomatic fraction (med 0.524)[hardcoded]
#     R0      ~ Lognormal (SAMPLED), the wet-season PEAK R_eff, from a scenario prior:
#             low = Caicedo 1.1-1.3 (~no outbreak); high = urban 1.18-3.51 (med ~2.03)
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
# R0 = wet-season PEAK, SAMPLED from a scenario prior. base_beta = R0 * gamma * season,
# so R_eff(t) = R0*season(t)*S/N is INDEPENDENT of gamma (gamma/sigma move only PEAK
# TIMING & HEIGHT; size is driven by R0 & immunity). SEASONAL-PEAK scaling (envelope
# rescaled so max = 1) makes the cited R0 the wet-season PEAK R_eff -- the honest,
# load-bearing quantity: it avoids the mean-1 framing that hides a higher true peak. The
# shape is the FITTED CHIKV beta_t envelope (same vector, same town), which is easier to
# justify than a rainfall proxy; the cost is that the peakedness is inherited from the
# CHIKV fit. R0 is drawn per LHS row from a Lognormal (R0_SCENARIO): "low" (Caicedo outside-
# Amazon 1.1-1.3) barely clears 1 at the peak -> no self-sustaining outbreak (deterministic
# gives a tiny outbreak; stochastically ~extinction); "high" (urban-adapted 1.18-3.51,
# med ~2.03) gives an outbreak with a genuine R0 band. NB R0(t) = R0*season(t) is just
# beta_t/gamma with an imposed seasonal shape (standard seasonal forcing), NOT a new
# formula. There is NO MAYV outbreak to fit, so this is PRIOR sampling, not fitting.
#
# WINDOW: 2025-W24 -> 2026-W22 (52 epi weeks), set by the CHIKV fit window. 2025 carries
# an epi-week 53, so 2025-W24..W53 = 30 wks and 2026-W01..W22 = 22 wks.
#
# FIXED (NOT sampled):
#   * seasonal envelope = FITTED CHIKV beta_t shape (caldas_beta_season.rds, written by
#     the Caldas Novas CHIKV fit: best_beta_t / mean(best_beta_t), mean-1, 52 weeks,
#     2025-W24..2026-W22). CHIKV and MAYV are assumed to share the same Aedes vector
#     season, so the fitted CHIKV transmission signal is the best available EMPIRICAL
#     proxy for MAYV seasonality -- more easily justified than a CHIRPS rainfall shape,
#     which is only a driver of the vector rather than a transmission measurement.
#     It is a DETERMINISTIC covariate (a fixed weekly shape), so it does NOT enter the LHS.
#   * single index case seeded at the wet-season onset (first above-mean week = index 19,
#     2025-W42).
# ============================================================
setwd("/Users/chloelee/Documents/R/summer_project")
suppressMessages({library(readxl); library(dplyr); library(tidyr); library(ggplot2)})

# ------------------------------------------------------------
# 1. SEIR forward simulator (identical to MAYV_ca_pre_vacc.R: seeded at seed_week)
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

T_weeks <- 52

# ------------------------------------------------------------
# 3. Seasonal envelope from the FITTED CHIKV beta_t (shared Aedes vector season)
# ------------------------------------------------------------
# caldas_beta_season.rds: mean-1 weekly transmission envelope for 2025-W24 -> 2026-W22,
# taken from the Caldas Novas CHIKV fit (best_beta_t / mean(best_beta_t)). Using the
# fitted CHIKV shape rather than a CHIRPS rainfall proxy is the more easily justified
# choice: it IS an observed transmission signal for the same vector in the same town.
season_mean1 <- readRDS("caldas_beta_season.rds")
stopifnot(length(season_mean1) == T_weeks, abs(mean(season_mean1) - 1) < 1e-6)
# R0 INTERPRETATION. r0_is_peak = TRUE: R0 is the SEASONAL-PEAK R_eff (rescale envelope
# so max = 1; yearly-avg R0 = R0*mean = ~0.58*R0). FALSE: R0 is the ANNUAL-MEAN (envelope
# left at mean 1; wet peak reaches R0*max ~ 1.7x). Seed/wet-band logic below always uses
# the mean-1 envelope so it is unaffected by this choice. NB the CHIKV beta_t envelope is
# FLATTER than the rainfall one it replaces (max 1.73 vs 2.44, min 0.58 vs 0.04).
r0_is_peak <- TRUE
season <- if (r0_is_peak) season_mean1 / max(season_mean1) else season_mean1

# ------------------------------------------------------------
# 4. Seed (single introduction at the wet-season onset)
# ------------------------------------------------------------
I0_total  <- 1
E0        <- rep(0, A)
seed_week <- which(season_mean1 >= 1)[1]   # single introduction at wet-season onset (first above-mean week)

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
rab <- c(a = 20, b = 60)                        # rho ~ Beta(20,60), mean 0.25
ps_a <- 35.84; ps_b <- 32.56                    # prop_symp ~ Beta, median 0.524

# R0 = wet-season PEAK R_eff (r0_is_peak = TRUE above), SAMPLED from a scenario range.
# The endpoints are LOWER/UPPER LIMITS (hard bounds), NOT a 95% UI: high = Dodero-Rojas
# 2020 estimated MAYV R0 range [1.18, 3.51]; low = Caicedo outside-Amazon [1.10, 1.30].
# We keep a right-skewed Lognormal SHAPE (median = geometric centre) but TRUNCATE it at
# the two limits (below), so every draw lands strictly within [lo, hi] -- no tail past
# the upper limit that would blow up the outbreak-size band.
if (!exists("R0_SCENARIO")) R0_SCENARIO <- "high"   # "high" Dodero-Rojas | "low" Caicedo (~no outbreak); overridable via a pre-set var
R0_priors <- list(
  low  = c(lo = 1.10, hi = 1.30),               # Caicedo outside-Amazon limits -> median ~1.20
  high = c(lo = 1.18, hi = 3.51))               # Dodero-Rojas 2020 MAYV R0 limits -> median ~2.03
r0p        <- R0_priors[[R0_SCENARIO]]
R0_lo      <- unname(r0p["lo"]); R0_hi <- unname(r0p["hi"])
R0_meanlog <- (log(R0_lo) + log(R0_hi)) / 2     # median at the geometric centre of the range
R0_sdlog   <- (log(R0_hi) - log(R0_lo)) / (2 * 1.96)
R0_median  <- unname(exp(R0_meanlog))

# Prior immunity (FLAT fraction), Lima et al. 2021 pooled CENTRAL-WEST exposure rate
# 8% (95% CI 3-18%, I2=98%). Right-skewed -> Lognormal calibrated to the CI endpoints:
# it reproduces [0.03,0.18] exactly (median 7.3%, mean ~8%). Swap this block for a
# Beta/other if the evidence changes.
imm_lo95 <- 0.03; imm_hi95 <- 0.18
imm_meanlog <- (log(imm_lo95) + log(imm_hi95)) / 2         # median = exp(meanlog) ~ 0.073
imm_sdlog   <- (log(imm_hi95) - log(imm_lo95)) / (2*1.96)
imm_base    <- 0.08                                        # Lima central estimate, for the baseline run

cat(sprintf("Priors: gamma~N(%.3f,%.4f)  latent~N(%.3f,%.4f)wk (sigma=1/latent)  rho~Beta(%d,%d)  prop_symp~Beta(%.2f,%.2f)  R0~truncLogN[%s] med %.2f range[%.2f,%.2f] (%s)  immune~logN(med %.3f, 95%%[%.2f,%.2f])\n",
            g_m, g_sd, p_m, p_sd, rab["a"], rab["b"], ps_a, ps_b,
            R0_SCENARIO, R0_median, R0_lo, R0_hi, if (r0_is_peak) "seasonal-peak" else "annual-mean",
            exp(imm_meanlog), imm_lo95, imm_hi95))

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

# Baseline at median inputs (dashed reference on the plots), scenario R0 median
base <- run_draw(R0_median, g_m, 1/p_m, rab["a"]/sum(rab), 0.5242478, imm_base)
cat(sprintf("Baseline R0=%.2f (%s median inputs, immune %.1f%%): total infections %.0f | total reported %.1f | peak reported wk %d (%.1f) | attack %.2f%%\n",
            R0_median, R0_SCENARIO, base$immune, base$tot_inf, base$tot_rep, base$peak_rep_wk, base$peak_rep, base$attack))
# Reference: Caicedo outside-Amazon PEAK R0=1.2 under the SAME scaling -> no self-sustaining outbreak
base_low <- run_draw(1.20, g_m, 1/p_m, rab["a"]/sum(rab), 0.5242478, imm_base)
cat(sprintf("Reference R0=1.20 (Caicedo peak, same scaling): total infections %.0f | total reported %.1f | attack %.3f%%  <- deterministic; stochastically ~extinction\n",
            base_low$tot_inf, base_low$tot_rep, base_low$attack))

# ------------------------------------------------------------
# 7. Latin Hypercube (6 inputs), forward-simulate each
# ------------------------------------------------------------
set.seed(2024); n <- 1000
lhs_col <- function(n) (sample.int(n) - runif(n)) / n
U   <- sapply(1:6, function(j) lhs_col(n))
gam <- qnorm(U[,1], g_m, g_sd)                  # rate
sig <- 1/qnorm(U[,2], p_m, p_sd)                # period -> rate
rho <- qbeta(U[,3], rab["a"], rab["b"])
psy <- qbeta(U[,4], ps_a, ps_b)
# TRUNCATED Lognormal on [R0_lo, R0_hi] (Dodero's hard limits): map the LHS uniform onto
# the [lo, hi] CDF span so every peak-R0 draw lands strictly inside the range (no tails).
R0v <- qlnorm(plnorm(R0_lo, R0_meanlog, R0_sdlog) +
              U[,5] * (plnorm(R0_hi, R0_meanlog, R0_sdlog) - plnorm(R0_lo, R0_meanlog, R0_sdlog)),
              R0_meanlog, R0_sdlog)
imm <- qlnorm(U[,6], imm_meanlog, imm_sdlog)    # flat prior-immune fraction (Lima 2021)

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
cat("\n=========== PROPAGATED MAYV RESULTS (", length(ok), " draws) ===========\n", sep = "")
cat(sprintf("Prior immunity:      propagated %.1f%% [%.1f%%, %.1f%%]\n",
            q3(immune[ok])[1], q3(immune[ok])[2], q3(immune[ok])[3]))
cat(sprintf("R0 at seasonal peak: baseline %.2f -> propagated %.2f [%.2f, %.2f]\n",
            R0_median*max(season), q3(R0peak[ok])[1], q3(R0peak[ok])[2], q3(R0peak[ok])[3]))
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
       subtitle = "Median (solid) + 95% band from gamma/sigma/rho/prop_symp/R0/immunity priors; baseline at median inputs (dashed)") +
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
write.csv(data.frame(draw = 1:n, R0 = R0v, gamma = gam, sigma = sig, rho = rho, prop_symp = psy,
                     immune_frac = imm, total_infections = tot_inf, total_reported = tot_rep,
                     peak_reported_wk = peak_rep_wk, attack_pct = attack,
                     finite = (seq_len(n) %in% ok)),
          "MAYV_ca_lhs_draws.csv", row.names = FALSE)

mayv_lhs_ensemble <- list(
  rep = rep_mat[ok, , drop = FALSE], inf = inf_mat[ok, , drop = FALSE],
  R0 = R0v[ok], gamma = gam[ok], sigma = sig[ok], rho = rho[ok], prop_symp = psy[ok],
  immune_frac = imm[ok],                          # per-draw FLAT immune fraction (Rimm = rep(imm, A))
  base_R0 = R0_median, R0_scenario = R0_SCENARIO, r0_is_peak = r0_is_peak,
  base_gamma = g_m, base_sigma = 1/p_m, base_prop_symp = 0.5242478,
  base_immune = imm_base,
  season = season, N = N, A = A, age_df = age_df,
  I0_total = I0_total, E0 = E0, seed_week = seed_week,
  T_weeks = T_weeks, weeks = weeks, x_ticks = x_ticks, year_break = year_break)
saveRDS(mayv_lhs_ensemble, "MAYV_ca_lhs_ensemble.rds")

cat("\nSaved MAYV_ca_lhs_reported.png, MAYV_ca_lhs_infections.png, MAYV_ca_lhs_draws.csv, MAYV_ca_lhs_ensemble.rds\n")
