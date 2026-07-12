library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
#  MAYV hypothetical single-season outbreak for Caldas Novas, Goias
#  Forward simulation (NO fitting): R0 is ASSUMED (Caicedo's outside-Amazon
#  1.1-1.3), not estimated -- there is no observed MAYV outbreak here to fit.
#
#  Apples-to-apples with the CHIKV model: MAYV rides the SAME vector-seasonality
#  envelope as the fitted Caldas Novas CHIKV beta_t (a single wet-season hump),
#  rescaled to MAYV's R0. The two diseases then differ only in R0 and natural
#  history. Because R0 ~ 1.1-1.3 is interpreted as the SEASONAL-PEAK R0 and the
#  off-season trough is sub-threshold, the outbreak is small and self-limiting:
#  it builds, peaks in the wet season, and resolves within the 52-week window.
# ============================================================

# ------------------------------------------------------------
# 1. SEIR function (time-varying weekly beta; can seed E0/I0)
# ------------------------------------------------------------
seir_baseline_MAYV <- function(
    T_weeks,        # number of weeks to simulate
    A,              # number of age groups
    N,              # vector length A: population by age group
    R_init_prop,    # vector length A: proportion already immune by age
    I0,             # vector length A: initial infectious by age (a stock)
    base_beta,      # vector length T_weeks: weekly transmission rate
    sigma,          # 1 / latent period (per week)
    gamma,          # recovery rate (per week)
    rho,            # reporting rate
    prop_symp,
    sub_steps = 7,
    E0 = rep(0, A),   # initial exposed by age
    seed_week = 1     # week at which the index case(s) I0 are introduced
) {
  pmax0 <- function(x) pmax(0, x)
  N_total <- sum(N)
  dt      <- 1 / sub_steps

  S <- E <- I <- R <- matrix(0, nrow = A, ncol = T_weeks)
  new_infections <- new_symptomatic <- matrix(0, nrow = A, ncol = T_weeks)

  # The index case is NOT present at t = 1; it is introduced at seed_week (below),
  # so before that the town is fully susceptible (minus baseline immunity).
  S_now <- pmax0(N - E0 - R_init_prop * N)
  E_now <- E0
  I_now <- rep(0, A)
  R_now <- R_init_prop * N

  for (t in 1:T_weeks) {
    new_I_week <- rep(0, A)
    beta_t <- base_beta[t]

    # Introduce the seed (one infectious person, split across ages) at seed_week.
    if (t == seed_week) {
      I_now <- I_now + I0
      S_now <- pmax0(S_now - I0)
    }

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
infectious_days <- 6      # viraemic period
gen_time_days   <- 15.2   # MAYV human-to-human generation time (Caicedo)
latent_days     <- gen_time_days - infectious_days   # = 9.2 d effective latent

sigma <- 7 / latent_days        # E -> I rate per week  (= 0.76)
gamma <- 7 / infectious_days    # I -> R rate per week  (= 1.17)

prop_symp <- 0.5242478
rho       <- 0.10               # MAYV reporting rate (severe under-ascertainment)

# NOTE on the latent period: 9.2 d is NOT the human intrinsic incubation (3 d).
# A human-only SEIR has no mosquito compartments, so E must absorb the whole vector
# loop. Mean generation time ~ (1/sigma)+(1/gamma) = latent+infectious, so to
# reproduce Caicedo's 15.2-d generation time we set latent = 15.2 - 6 = 9.2 d.
# R0 and final size depend only on gamma; this choice changes only PEAK TIMING.

# ------------------------------------------------------------
# 3. Load Caldas Novas population; set MAYV baseline immunity
#    PRIMARY ASSUMPTION: the population is NAIVE to Mayaro (baseline_immunity = 0).
#    There is no recorded MAYV circulation in Caldas Novas, so a fully susceptible
#    town is the natural starting point. We then VARY the immune fraction in
#    sensitivity analysis using Lima et al. (2021)'s pooled MAYV seroprevalence:
#    8% (95% CI 3-18%; note I2 = 98%, i.e. very heterogeneous across studies).
#    Immunity is applied as a flat fraction across all age groups (the pooled
#    estimate is not age-resolved); we no longer use an age-specific catalytic FOI.
# ------------------------------------------------------------
age_df <- read_excel("population.xlsx", sheet = "prop_immune")
age_df <- age_df[tolower(age_df$municipality) == "caldas novas", ]
stopifnot(nrow(age_df) > 0)
age_df$age_group <- factor(age_df$age_group, levels = age_df$age_group)
A <- nrow(age_df)

# Grow 2022 age-group counts to 2025 (outbreak year) with the town growth rate,
# mirroring the CHIKV Caldas Novas model.
pop_2022_total <- 98622
pop_2025_total <- 106820
growth_r       <- log(pop_2025_total / pop_2022_total) / 3
N              <- age_df$pop_num * exp(growth_r * 3)

# --- MAYV baseline immunity ----------------------------------------------------
# Primary = 0 (naive). Sensitivity band (Lima et al. 2021): 0.08 central, 0.03 low,
# 0.18 high. Set baseline_immunity to one of these. NB at R0 <= 1.3 the outbreak
# stays super-critical across this whole band (max 18% immune -> S(0) = 0.82).
baseline_immunity <- 0.0
R_init_prop <- rep(baseline_immunity, A)

pct_immune <- 100 * baseline_immunity
cat(sprintf("Baseline MAYV immunity = %.0f%% (naive primary; vary 3/8/18%% per Lima 2021)\n",
            pct_immune))

susceptible_pop <- N * (1 - R_init_prop)
S0_frac <- sum(susceptible_pop) / sum(N)

T_weeks <- 52    # single-season window 2025-W24 -> 2026-W22 = 52 weeks (2025 has an
                 # epi-week 53, so 2025 = W24..W53 = 30 wks), aligned 1:1 with the
                 # 52-week CHIKV Caldas fit.

# ------------------------------------------------------------
# 4. Seasonal transmission envelope (the "vector seasonality"), taken from the
#    FITTED Caldas Novas CHIKV beta_t shape so MAYV and CHIKV share one season.
# ------------------------------------------------------------
# PROVENANCE: the mean-normalised fitted Caldas Novas CHIKV beta_t shape, written by
# CHIKV_ca_pre_vacc_optim.R as caldas_beta_season.rds (= best_beta_t / mean(best_beta_t),
# mean = 1, full 52-week window 2025-W24 -> 2026-W22). This loads it dynamically so the
# seasonal shape stays in sync with the current CHIKV fit instead of being hard-coded.
# REQUIRES running CHIKV_ca_pre_vacc_optim.R first to (re)generate the file.
season_mean1 <- readRDS("caldas_beta_season.rds")
stopifnot(length(season_mean1) == T_weeks, abs(mean(season_mean1) - 1) < 1e-6)

# R0 SCALING CHOICE -------------------------------------------------------------
# r0_is_peak = FALSE: R0 (1.1-1.3) is the ANNUAL-MEAN R0 (envelope mean = 1), so
#                     the wet-season PEAK instantaneous R0 rises to ~1.79x R0
#                     (= 2.0-2.3) while the dry-season trough (~0.48x) is
#                     sub-threshold. This is the DEFAULT: it keeps the yearly
#                     average at Caicedo's outside-Amazon value yet concentrates
#                     transmission into the wet season, yielding a small
#                     self-limiting outbreak that peaks in March and resolves
#                     within the year (busiest week ~0.6-15 reported across the
#                     R0 range, bracketing Goiania's ~10/wk reference).
# r0_is_peak = TRUE : R0 is instead the SEASONAL-PEAK R0 (envelope peak = 1), so
#                     R_eff only marginally clears 1 for a few weeks. With MAYV's
#                     short infectious period this CANNOT self-sustain from a
#                     single introduction -- you get sporadic spillover, not an
#                     outbreak. A more conservative reading; needs continuous
#                     importation to show any wet-season clustering.
r0_is_peak <- FALSE
season <- if (r0_is_peak) season_mean1 / max(season_mean1) else season_mean1

# ------------------------------------------------------------
# 5. Seed the epidemic (a SINGLE introduced case at the START of the wet season)
# ------------------------------------------------------------
# MAYV reaches the town through one introduction. We introduce it at the ONSET of the
# wet season (week_index 26 ~ 2025-W49, early Dec), NOT in mid-2025: a single case
# dropped into the dry season would just decay before transmission can support it, and
# realistically introductions that ignite an outbreak arrive when conditions already
# favour transmission. This replaces the earlier ongoing-importation scenario, which
# would have required modelling human movement (out of scope). We split the one case
# across age groups in proportion to susceptibility (a fractional seed is fine in a
# deterministic model and avoids an arbitrary age pick).
I0_total  <- 1
I0        <- I0_total * susceptible_pop / sum(susceptible_pop)
E0        <- rep(0, A)
seed_week <- 26        # introduce the index case at the wet-season onset (~2025-W49)

# ------------------------------------------------------------
# 6. Scenario runner: weekly beta = R0 * gamma * seasonal envelope
# ------------------------------------------------------------
# R0 = 1.2 is the CENTRAL estimate; 1.1 and 1.3 are SENSITIVITY scenarios (low/high
# transmissibility), NOT the bounds of the 95% CI. The CI is built separately around
# R0 = 1.2 by propagating the fitted CHIKV seasonal-shape posterior (section 10).
R0_central <- 1.2
R0_sens    <- c(1.1, 1.3)
R0_grid    <- c(R0_central, R0_sens)   # used only for the summary / sensitivity tables

run_scenario <- function(R0) {
  base_beta <- R0 * gamma * season
  seir_baseline_MAYV(
    T_weeks = T_weeks, A = A, N = N, R_init_prop = R_init_prop,
    I0 = I0, base_beta = base_beta, sigma = sigma, gamma = gamma,
    rho = rho, prop_symp = prop_symp, E0 = E0, seed_week = seed_week
  )
}

# ------------------------------------------------------------
# 7. Threshold pre-check: peak-season R_eff (constant-R0 final-size theory does
#    NOT apply under seasonal forcing, so we judge take-off by the peak R_eff and
#    read the outbreak SIZE off the simulation).
# ------------------------------------------------------------
peak_R0_of <- function(R0) R0 * max(season) / 1            # max R0 over the season
diag_tbl <- data.frame(
  R0            = R0_grid,
  peak_R0       = round(sapply(R0_grid, peak_R0_of), 3),
  peak_R_eff    = round(sapply(R0_grid, function(R0) peak_R0_of(R0) * S0_frac), 3),
  trough_R_eff  = round(sapply(R0_grid, function(R0) R0 * min(season) * S0_frac), 3)
)
cat(sprintf("Susceptible fraction S(0): %.3f  (immunity %.1f%%)\n", S0_frac, pct_immune))
cat("Outbreak needs peak_R_eff > 1; trough_R_eff < 1 makes it self-limiting.\n")
print(diag_tbl)

# ------------------------------------------------------------
# 8. Scenario summary over the season
# ------------------------------------------------------------
summary_df <- do.call(rbind, lapply(R0_grid, function(R0) {
  out    <- run_scenario(R0)
  wk_rep <- colSums(out$new_reported)
  wk_inf <- colSums(out$new_infections)
  data.frame(
    R0                = R0,
    peak_R_eff        = round(peak_R0_of(R0) * S0_frac, 3),
    total_infections  = round(sum(wk_inf)),
    total_symptomatic = round(sum(wk_inf) * prop_symp),
    total_reported    = round(sum(wk_rep)),
    peak_reported_wk  = round(max(wk_rep), 1),       # busiest week (reported)
    peak_week         = which.max(wk_rep),
    attack_rate_susc  = round(sum(wk_inf) / sum(susceptible_pop), 4)
  )
}))
cat("\n--- MAYV single-season scenario summary (Caldas Novas) ---\n")
cat("peak_reported_wk = highest weekly reported count (compare with Goiania 2018:\n",
    "  busiest MAYV week ~10 reported cases).\n", sep = "")
print(summary_df)

# ------------------------------------------------------------
# 9. Axis labels + wet/dry season shading helpers for the plots
# ------------------------------------------------------------
# Map a within-window index (1 = 2025-W24) to the plain epidemiological week number.
# 2025 has an epi-week 53, so index 1-30 = 2025 (W24-W53) and index 31-52 = 2026 (W01-W22).
wk_num <- function(idx) ifelse(idx <= 30, idx + 23, idx - 30)   # 2025: +23 ; 2026: -30
tick_idx <- c(7, 17, 27, 40, 50)                      # 2025-W30/40/50, 2026-W10/20
# Fixed 2025|2026 boundary sits between index 30 (2025-W53) and index 31 (2026-W01).
# (No observed-data lookup here -- this is a forward MAYV sim, so caldas_obs does not exist.)
year_break <- 30.5
# x-axis tick table: labels are the plain epi-week NUMBERS (30,40,50,10,20), not dates.
x_ticks <- data.frame(week_index = tick_idx, week = wk_num(tick_idx))

# Wet vs dry season bands. The window starts 2025-W24 (early June). Using an
# approximate epi-week -> month mapping, the WET season (Dec-Apr) spans week-index
# ~26-48; the rest (Jun-Nov 2025 and May-Jun 2026) is DRY (May-Nov).
wet_start <- 26   # ~2025-W49 = early Dec 2025
wet_end   <- 48   # ~2026-W18 = early May 2026

# Reusable layers: pale wet/dry rectangles behind the data + season labels pinned to
# the BOTTOM of the panel. annotate("rect") uses fixed fills (no fill scale), so it
# does not clash with a scenario fill/colour scale in the comparison plot.
season_layers <- list(
  annotate("rect", xmin = 0.5,       xmax = wet_start,     ymin = -Inf, ymax = Inf,
           fill = "#f4ead7", alpha = 0.55),
  annotate("rect", xmin = wet_start, xmax = wet_end,       ymin = -Inf, ymax = Inf,
           fill = "#cfe6f2", alpha = 0.55),
  annotate("rect", xmin = wet_end,   xmax = T_weeks + 0.5, ymin = -Inf, ymax = Inf,
           fill = "#f4ead7", alpha = 0.55),
  annotate("text", x = (0.5 + wet_start)/2,  y = -Inf, label = "Dry season (May-Nov)",
           vjust = -0.8, size = 3, colour = "grey35"),
  annotate("text", x = (wet_start + wet_end)/2, y = -Inf, label = "Wet season (Dec-Apr)",
           vjust = -0.8, size = 3, colour = "grey35"),
  annotate("text", x = (wet_end + T_weeks)/2, y = -Inf, label = "Dry",
           vjust = -0.8, size = 3, colour = "grey35")
)

# ------------------------------------------------------------
# 10. 95% CI around R0 = 1.2 from the CHIKV beta_t POSTERIOR (NOT from 1.1/1.3)
# ------------------------------------------------------------
# R0 is held at the central 1.2. The uncertainty is propagated from the fitted CHIKV
# seasonal shape: each posterior draw of beta_t (saved by CHIKV_ca_pre_vacc_optim.R,
# normalised to mean 1) is rescaled to R0*gamma and simulated. The 2.5/50/97.5
# percentiles across draws give the band. This is a PARAMETER CI (uncertainty in the
# borrowed seasonality); it is narrow because the CHIKV beta_t is well estimated. For
# a wider PREDICTION interval one would additionally add negative-binomial reporting
# noise around each weekly mean.
shape_samples <- as.matrix(read.csv("chikv_ca_beta_shape_samples.csv"))
if (r0_is_peak) shape_samples <- shape_samples / apply(shape_samples, 1, max)  # match scaling
n_post <- nrow(shape_samples)

# Central (point-estimate) weekly reported curve for a given R0 (single seed at
# seed_week; no importation).
central_curve <- function(R0) {
  colSums(seir_baseline_MAYV(
    T_weeks, A, N, R_init_prop, I0, R0 * gamma * season, sigma, gamma, rho, prop_symp,
    E0 = E0, seed_week = seed_week)$new_reported)
}

# 95% CI band: run every posterior shape draw at fixed R0, quantile across draws.
ci_band <- function(R0) {
  M <- matrix(NA_real_, n_post, T_weeks)
  for (i in seq_len(n_post)) {
    M[i, ] <- colSums(seir_baseline_MAYV(
      T_weeks, A, N, R_init_prop, I0, R0 * gamma * shape_samples[i, ],
      sigma, gamma, rho, prop_symp, E0 = E0, seed_week = seed_week)$new_reported)
  }
  data.frame(week = 1:T_weeks,
             lo  = apply(M, 2, quantile, 0.025),
             hi  = apply(M, 2, quantile, 0.975))
}

# ============================================================
#  GRAPH: R0 = 1.2 central + 95% CI (single introduction at the wet-season onset)
# ============================================================
g1 <- ci_band(R0_central)
g1$central <- central_curve(R0_central)

p_mayv_mc <- ggplot(g1, aes(week)) +
  season_layers +
  annotate("text", x = year_break, y = 0, label = "2026", angle = 90,
           vjust = -0.5, hjust = -11.5, fontface = "bold", size = 3.5) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  # geom_vline(xintercept = seed_week, linetype = "dotted", colour = "grey45") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey40", alpha = 0.30) +
  geom_line(aes(y = central), colour = "grey40", linewidth = 1) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  labs(x = "Week", y = "Predicted reported MAYV cases",
       title = expression(bold("Hypothetical MAYV outbreak"))) +
       # subtitle = sprintf("Central + 95%% CI from CHIKV seasonal-shape posterior; %.0f%% immunity; one introduction at wet-season onset",
       #                    pct_immune)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", colour = "grey85"))
print(p_mayv_mc)
ggsave("mayv_ca_vacc_infections.png", p_mayv_mc, width = 8, height = 4.5, dpi = 120)