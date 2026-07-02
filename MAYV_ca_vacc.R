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
# SINGLE SCENARIO: a PRE-OUTBREAK campaign (vaccinate before the wet-season seed).
# This is the only timing of interest for MAYV -- a reactive roll-out during a
# small, short, self-limiting outbreak averts almost nothing -- so the earlier
# "actual rollout" and "start of 2026" timings are dropped. Two efficacy arms,
# as for CHIKV:
#   - Disease-blocking ONLY (VE_inf = 0, VE_block = eff): conservative floor;
#     averts symptomatic/hospitalisations/deaths but NOT infections.
#   - Disease + infection blocking (VE_inf = eff, VE_block = eff): optimistic
#     ceiling; susceptibles move S->V so infections are averted too.
#
# IMPORTANT CAVEATS
#   1. There is NO licensed Mayaro vaccine. Efficacy is a HYPOTHETICAL: we borrow
#      IXCHIQ's value as an "alphavirus-vaccine analogue" (and/or assumed IXCHIQ
#      cross-protection). It is exposed as one lever (`mayv_vacc_efficacy`).
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

# HYPOTHETICAL MAYV vaccine efficacy. No licensed MAYV product exists; we borrow
# IXCHIQ's 98.9% as an upper-bound "alphavirus analogue" / assumed-cross-protection
# value. This is the headline assumption -- lower it to model a weaker product.
mayv_vacc_efficacy    <- 0.989

# Pre-outbreak campaign: begin well before the wet-season seed (2025-W26 = index 4,
# vs seed_week = 27), so the immunised fraction is in place before any local
# transmission. This is the maximum-achievable ("ceiling") timing.
# week_to_index() is defined in ca_common.R.
start_pre <- week_to_index(2025, 26)   # 2025-W26 -> index 4

cat("Eligible 18-59 population:", round(target_pop), "\n")
cat("Target doses (30% coverage):", round(target_pop * total_coverage), "\n")
cat(sprintf("HYPOTHETICAL MAYV vaccine efficacy: %.1f%% (no licensed product)\n",
            100 * mayv_vacc_efficacy))
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
# 4. Scenario runner + run baseline and the pre-outbreak campaign (2 arms)
# ============================================================
# R0 enters through the seasonal beta exactly as in the pre-vacc model:
#   base_beta = R0 * gamma * season   (season has mean 1 over the window).
# season_use / prop_symp_use let the Monte Carlo (section 7) swap in draws.
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

# No-vaccine counterfactual at the central R0 (reproduces the pre-vacc outbreak)
out_baseline <- run_scenario(R0_central, 0, start_pre, VE_inf = 0, VE_block = 0)

# Single timing (pre-outbreak); two efficacy arms, disease-blocking listed FIRST.
timings <- list("pre-outbreak" = start_pre)
arms <- list(
  "Disease-blocking"             = c(VE_inf = 0,                  VE_block = mayv_vacc_efficacy),
  "Disease + infection blocking" = c(VE_inf = mayv_vacc_efficacy, VE_block = mayv_vacc_efficacy)
)

scenarios <- list("No vaccine (baseline)" = out_baseline)
scen_meta <- data.frame(name = "No vaccine (baseline)", timing = "No vaccine",
                        arm = "No vaccine", start = NA_integer_)
for (tn in names(timings)) {
  for (an in names(arms)) {
    nm <- paste0(tn, " | ", an)
    scenarios[[nm]] <- run_scenario(R0_central, total_coverage, timings[[tn]],
                                    arms[[an]]["VE_inf"], arms[[an]]["VE_block"])
    scen_meta <- rbind(scen_meta,
                       data.frame(name = nm, timing = tn, arm = an, start = timings[[tn]]))
  }
}

# ============================================================
# 5. Burden comparison (point estimates at R0 = R0_central)
# ============================================================
burden_mat <- t(sapply(scenarios, burden))          # rows = scenario, cols = outcome
base_b     <- burden_mat["No vaccine (baseline)", ]

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

cat(sprintf("\n=== Baseline burden (no vaccine, R0 = %.1f) ===\n", R0_central))
print(round(base_b, 1))
cat("\n=== Total burden by scenario ===\n")
print(round(burden_mat, 1))
cat("\n=== Averted vs no-vaccine baseline (pre-outbreak campaign) ===\n")
print(summary_tbl, row.names = FALSE)

write.csv(summary_tbl, "caldas_mayv_vacc_averted.csv", row.names = FALSE)
cat("\nWrote caldas_mayv_vacc_averted.csv\n")

# ============================================================
# 6. Plots
# ============================================================
# (a) Weekly TRUE infections: baseline + BOTH efficacy arms of the pre-outbreak
#     campaign. NB the disease-blocking arm (VE_inf = 0) does not move anyone S->V,
#     so its infection curve sits exactly on the baseline -- drawn dashed so it is
#     still visible on top of the baseline line.
pre_lbl <- "Vaccine"
inf_keys <- c("No vaccine (baseline)",
              "pre-outbreak | Disease-blocking",
              "pre-outbreak | Disease + infection blocking")
inf_lbls <- c("No vaccine (baseline)",
              paste0(pre_lbl, ": disease-blocking"),
              paste0(pre_lbl, ": disease + infection blocking"))
inf_df <- do.call(rbind, Map(function(nm, lbl) {
  data.frame(week = 1:T_weeks, infections = colSums(scenarios[[nm]]$new_infections), scenario = lbl)
}, inf_keys, inf_lbls))
inf_df$scenario <- factor(inf_df$scenario, levels = inf_lbls)

p_inf <- ggplot(inf_df, aes(week, infections, colour = scenario, linetype = scenario)) +
  season_layers +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  # geom_vline(xintercept = seed_week, linetype = "dotted", colour = "grey45") +
  annotate("text", x = year_break, y = 0, label = "2026", angle = 90,
           vjust = -0.5, hjust = -11.5, fontface = "bold", size = 3.5) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = setNames(c("grey40", "#9ecae1", "#08519c"), inf_lbls),
                      name = NULL) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "solid"), inf_lbls),
                        name = NULL) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week", y = "Predicted true MAYV cases",
       title = "Hypothetical MAYV cases") +
  guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside", legend.position.inside = c(0.07, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p_inf
ggsave("ca_mayv_vacc_infections.png", p_inf, width = 8, height = 4.5, dpi = 120)

# (b) Averted burden: facet by outcome, x = timing (pre-outbreak), fill = arm
bar_df <- summary_tbl |>
  dplyr::select(timing, arm, Infections = inf_averted, Symptomatic = symp_averted,
                Hospitalisations = hosp_averted, Deaths = deaths_averted) |>
  pivot_longer(c(Infections, Symptomatic, Hospitalisations, Deaths),
               names_to = "outcome", values_to = "averted")
bar_df$timing  <- factor(bar_df$timing, levels = names(timings))
bar_df$arm     <- factor(bar_df$arm, levels = names(arms))   # disease-blocking first
bar_df$outcome <- factor(bar_df$outcome,
                         levels = c("Infections", "Symptomatic", "Hospitalisations", "Deaths"))

# p_bar <- ggplot(bar_df, aes(timing, averted, fill = arm)) +
#   geom_col(position = position_dodge(width = 0.8), width = 0.7) +
#   facet_wrap(~ outcome, scales = "free_y", nrow = 1) +
#   scale_fill_manual(values = c("Disease-blocking" = "#9ecae1",
#                                "Disease + infection blocking" = "#08519c"),
#                     name = NULL) +
#   scale_y_continuous(labels = scales::comma) +
#   labs(x = NULL, y = "Averted (vs no vaccine)",
#        title = "Caldas Novas MAYV cases - burden averted") +
#   theme_bw(base_size = 11) +
#   theme(plot.title = element_text(face = "bold", hjust = 0.5),
#         axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
#         legend.position = "bottom", panel.grid.minor = element_blank())
# p_bar
# ggsave("ca_mayv_vacc_averted.png", p_bar, width = 10, height = 4.2, dpi = 120)

cat("\nSaved plot: ca_mayv_vacc_infections.png\n")

# ============================================================
# 7. Monte Carlo: averted burden with 95% uncertainty intervals
# ============================================================
# MAYV has no fitted likelihood, so the 95% UI propagates input-parameter
# uncertainty, all PAIRED per draw across baseline and each scenario (so the shared
# draw carries through to the averted). R0 is held FIXED at the central value
# (R0_central = 1.2); the 1.1-1.3 range is reserved for the separate R0 sensitivity
# analysis (section 10), NOT used to build this CI:
#   (a) seasonal shape  -- a random draw from the CHIKV beta_t posterior
#                          (chikv_ca_beta_shape_samples.csv, mean 1)
#   (b) prop_symp, hosp rate, age-specific CFR -> Beta draws (Hyolim Table S4).
# Vaccine efficacy is fixed at mayv_vacc_efficacy. (As in the grid, total infections
# are ~pinned by the fixed R0, so their UI is near-degenerate.)
shape_samples <- as.matrix(read.csv("chikv_ca_beta_shape_samples.csv"))
n_post <- nrow(shape_samples)

n_draws <- 1000
set.seed(2026)
shape_draws <- sample.int(n_post, n_draws, replace = TRUE)

vac_names  <- setdiff(names(scenarios), "No vaccine (baseline)")
outcomes   <- c("infections", "symptomatic", "hospitalisations", "deaths")
base_draws <- matrix(NA_real_, n_draws, 4, dimnames = list(NULL, outcomes))
av_draws   <- setNames(lapply(vac_names, function(x)
                  matrix(NA_real_, n_draws, 4, dimnames = list(NULL, outcomes))), vac_names)

cat("\nRunning", n_draws, "Monte Carlo draws (3 SEIR runs each)...\n")
for (i in 1:n_draws) {
  season_i <- shape_samples[shape_draws[i], ]

  # input-parameter draws (borrowed CHIKV severity)
  ps_i      <- rbeta(1, ps_a, ps_b)
  hosp_i    <- rbeta(1, hosp_a, hosp_b)
  cfr_h_i   <- rbeta(9, cfr_hosp_a, cfr_hosp_b)
  cfr_n_i   <- rbeta(9, cfr_nonh_a, cfr_nonh_b)
  cfr_vec_i <- (hosp_i * cfr_h_i + (1 - hosp_i) * cfr_n_i)[age_to_band]

  out_b <- run_scenario(R0_central, 0, start_pre, 0, 0, season_use = season_i, prop_symp_use = ps_i)
  bb    <- burden(out_b, hosp_i, cfr_vec_i)
  base_draws[i, ] <- bb

  for (tn in names(timings)) for (an in names(arms)) {
    nm  <- paste0(tn, " | ", an)
    out <- run_scenario(R0_central, total_coverage, timings[[tn]], arms[[an]]["VE_inf"],
                        arms[[an]]["VE_block"], season_use = season_i, prop_symp_use = ps_i)
    av_draws[[nm]][i, ] <- bb - burden(out, hosp_i, cfr_vec_i)
  }
  if (i %% 100 == 0) cat("  ", i, "/", n_draws, "\n")
}

# ---- Summarise: median (2.5% - 97.5%); fmtq() is defined in ca_common.R ----

cat("\n=== Baseline burden, no vaccine (median, 95% UI) ===\n")
for (o in outcomes)
  cat(sprintf("  %-16s %s\n", o, fmtq(base_draws[, o], if (o == "deaths") 1 else 0)))

mc_tbl <- data.frame(
  timing = scen_meta$timing[match(vac_names, scen_meta$name)],
  arm    = scen_meta$arm[match(vac_names, scen_meta$name)],
  Infections       = sapply(vac_names, function(nm) fmtq(av_draws[[nm]][, "infections"])),
  Symptomatic      = sapply(vac_names, function(nm) fmtq(av_draws[[nm]][, "symptomatic"])),
  Hospitalisations = sapply(vac_names, function(nm) fmtq(av_draws[[nm]][, "hospitalisations"], 1)),
  Deaths           = sapply(vac_names, function(nm) fmtq(av_draws[[nm]][, "deaths"], 2)),
  row.names = NULL, check.names = FALSE
)
cat("\n=== Averted vs no-vaccine baseline (median, 95% UI) ===\n")
cat("    (outcomes: infections, symptomatic, hospitalisations, deaths)\n")
old_width <- getOption("width"); options(width = 220)   # avoid wrapping the Deaths column
print(mc_tbl, row.names = FALSE)
options(width = old_width)
write.csv(mc_tbl, "caldas_mayv_vacc_averted_mc.csv", row.names = FALSE)
cat("\nWrote caldas_mayv_vacc_averted_mc.csv\n")

# ---- Plot: averted burden with 95% UI error bars ----
mc_long <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- av_draws[[nm]]
  data.frame(timing  = scen_meta$timing[scen_meta$name == nm],
             arm     = scen_meta$arm[scen_meta$name == nm],
             outcome = outcomes,
             med = apply(m, 2, median, na.rm = TRUE),
             lo  = apply(m, 2, quantile, .025, na.rm = TRUE),
             hi  = apply(m, 2, quantile, .975, na.rm = TRUE),
             row.names = NULL)
}))
mc_long$timing  <- factor(mc_long$timing, levels = names(timings))
mc_long$arm     <- factor(mc_long$arm, levels = names(arms))   # disease-blocking first
mc_long$outcome <- factor(mc_long$outcome,
                          levels = c("infections", "symptomatic", "hospitalisations", "deaths"),
                          labels = c("Infections", "Symptomatic", "Hospitalisations", "Deaths"))

p_mc <- ggplot(mc_long, aes(timing, med, fill = arm)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8),
                width = 0.25, linewidth = 0.4) +
  facet_wrap(~ outcome, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Disease-blocking" = "#9ecae1",
                               "Disease + infection blocking" = "#08519c"), name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Burden averted, median + 95% UI",
       title = "Hypothetical MAYV cases - burden averted") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom", panel.grid.minor = element_blank())
p_mc
ggsave("ca_mayv_vacc_averted_mc.png", p_mc, width = 10, height = 4.2, dpi = 120)
cat("Saved plot: ca_mayv_vacc_averted_mc.png\n")

# ============================================================
# 8. TWO-LAYER comparison: outbreak size x (efficacy level x vaccine mechanism)
# ============================================================
# Layer 1 -- OUTBREAK SIZE (R0):
#   * Current  R0 = 1.2  : outside the Amazon basin (Caicedo et al. 2021, 1.1-1.3).
#   * Future   R0 = 3.51 : urban transmission adapted to Ae. aegypti, the UPPER
#                          bound of Dodero-Rojas et al. 2020 (1.18-3.51).
#   Both keep the SAME convention as the rest of the script: R0 is the annual-mean
#   scaling of the (mean-1) seasonal envelope, so base_beta = R0 * gamma * season.
# Layer 2 -- VACCINE: an efficacy LEVEL crossed with a MECHANISM:
#   * No vaccine                              (eff = 0; single reference bar)
#   * 50% efficacy   x {disease-blocking, disease + infection blocking}
#   * 98.9% efficacy x {disease-blocking, disease + infection blocking}
#   The two MECHANISMS (the two basic vaccine modes being compared):
#     - Disease-blocking ONLY       : VE_inf = 0,   VE_block = eff  (no infection block)
#     - Disease + infection blocking: VE_inf = eff, VE_block = eff
# => per R0 size: 1 + 2 x 2 = 5 scenarios; x 2 sizes = 10 scenarios.
# Reported as TOTAL burden (infections, symptomatic, hospitalisations, deaths),
# pre-outbreak campaign timing throughout. (Severity is the borrowed CHIKV proxy.)
library(writexl)

size_R0 <- c("Current outbreak (R0 = 1.2)"  = R0_central,
             "Future outbreak (R0 = 3.51)"  = 3.51)   # R0 FIXED at each layer value
eff_lvl <- c("50% efficacy" = 0.50, "98.9% efficacy" = mayv_vacc_efficacy)
# NOTE on uncertainty: R0 is held FIXED at the layer value (1.2 / 3.51). The 95% UI
# is NOT built from the 1.1-1.3 R0 range -- that range is reserved for the separate
# R0 sensitivity analysis (section 10). The UI instead propagates input-parameter
# uncertainty: the borrowed CHIKV seasonal-shape posterior + the severity Betas
# (prop_symp, hosp rate, age-CFR). Consequence: total INFECTIONS are ~pinned by the
# fixed R0 (final-size theory), so their UI is near-degenerate; the genuine UI lives
# in symptomatic / hospitalisations / deaths via the severity parameters.

# Vaccine design rows: a no-vaccine reference + (efficacy level x mechanism).
mech <- data.frame(eff = "No vaccine", arm = "No vaccine",
                   VE_inf = 0, VE_block = 0, coverage = 0, stringsAsFactors = FALSE)
for (el in names(eff_lvl)) {
  e <- eff_lvl[[el]]
  mech <- rbind(mech,
    data.frame(eff = el, arm = "Disease-blocking",
               VE_inf = 0, VE_block = e, coverage = total_coverage),
    data.frame(eff = el, arm = "Disease + infection blocking",
               VE_inf = e, VE_block = e, coverage = total_coverage))
}
# Cross every vaccine design with every outbreak size (size varies slowest).
perm <- do.call(rbind, lapply(names(size_R0), function(s)
  cbind(size = s, mech, stringsAsFactors = FALSE)))
perm$key <- paste(perm$size, perm$eff, perm$arm, sep = " | ")
np <- nrow(perm)   # 10

# One scenario at a given R0 (coverage 0 = no vaccine). season_use / prop_symp_use
# let the Monte Carlo swap in draws.
run_perm <- function(R0, coverage, VE_inf, VE_block,
                     season_use = season, prop_symp_use = prop_symp) {
  run_scenario(R0, coverage, start_pre, VE_inf = VE_inf, VE_block = VE_block,
               season_use = season_use, prop_symp_use = prop_symp_use)
}

# ---- Point estimates (R0 fixed at each layer value, central inputs) ----
pe <- t(sapply(seq_len(np), function(j)
  burden(run_perm(size_R0[[perm$size[j]]], perm$coverage[j], perm$VE_inf[j], perm$VE_block[j]))))
rownames(pe) <- perm$key

cat("\n=== TOTAL burden by permutation (point estimate) ===\n")
print(round(pe, 2))

# ---- Monte Carlo 95% UI: R0 FIXED per layer; propagate seasonal shape + severity
#      (prop_symp, hosp rate, CFR), paired per draw across all scenarios. ----
set.seed(2027)
shape_draws2 <- sample.int(n_post, n_draws, replace = TRUE)
perm_mc <- lapply(seq_len(np), function(x)
  matrix(NA_real_, n_draws, 4, dimnames = list(NULL, outcomes)))

cat("\nRunning", n_draws, "MC draws for", np, "permutations (", np, "SEIR runs each)...\n")
for (i in 1:n_draws) {
  season_i  <- shape_samples[shape_draws2[i], ]
  ps_i      <- rbeta(1, ps_a, ps_b)
  hosp_i    <- rbeta(1, hosp_a, hosp_b)
  cfr_h_i   <- rbeta(9, cfr_hosp_a, cfr_hosp_b)
  cfr_n_i   <- rbeta(9, cfr_nonh_a, cfr_nonh_b)
  cfr_vec_i <- (hosp_i * cfr_h_i + (1 - hosp_i) * cfr_n_i)[age_to_band]

  for (j in seq_len(np)) {
    out <- run_perm(size_R0[[perm$size[j]]], perm$coverage[j], perm$VE_inf[j],
                    perm$VE_block[j], season_use = season_i, prop_symp_use = ps_i)
    perm_mc[[j]][i, ] <- burden(out, hosp_i, cfr_vec_i)
  }
  if (i %% 100 == 0) cat("  ", i, "/", n_draws, "\n")
}

# ---- Tables (point estimate + median 95% UI) and Excel workbook ----
mk_tbl <- function(value_fun) data.frame(
  `Outbreak size`     = perm$size,
  `Vaccine efficacy`  = perm$eff,
  `Vaccine mechanism` = perm$arm,
  Infections       = value_fun("infections"),
  Symptomatic      = value_fun("symptomatic"),
  Hospitalisations = value_fun("hospitalisations"),
  Deaths           = value_fun("deaths"),
  check.names = FALSE, row.names = NULL
)
dgt <- c(infections = 0, symptomatic = 0, hospitalisations = 1, deaths = 2)
tbl_pe <- mk_tbl(function(o) round(pe[, o], dgt[[o]]))
tbl_ui <- mk_tbl(function(o) sapply(seq_len(np), function(j) fmtq(perm_mc[[j]][, o], dgt[[o]])))

cat("\n=== TOTAL burden by permutation (median, 95% UI) ===\n")
old_width <- getOption("width"); options(width = 240)
print(tbl_ui, row.names = FALSE)
options(width = old_width)

write_xlsx(list(point_estimate = tbl_pe, median_95UI = tbl_ui),
           "caldas_mayv_vacc_6perm_burden.xlsx")
cat("\nWrote caldas_mayv_vacc_6perm_burden.xlsx (sheets: point_estimate, median_95UI)\n")

# ---- Plot: total burden. x = efficacy level, fill = mechanism (two arms dodged);
#      no-vaccine is a single reference bar. One panel per (size x outcome). ----
library(ggh4x)   # facet_grid2: grouped size headers + independent per-panel y-axes
eff_x_levels  <- c("No vaccine", names(eff_lvl))
arm_levels    <- c("No vaccine", "Disease-blocking", "Disease + infection blocking")
arm_cols      <- c("No vaccine" = "grey60", "Disease-blocking" = "#9ecae1",
                   "Disease + infection blocking" = "#08519c")

mc_long <- function(rows, mats) do.call(rbind, lapply(seq_along(rows), function(k) {
  m <- mats[[k]]; j <- rows[k]
  data.frame(size = perm$size[j], eff = perm$eff[j], arm = perm$arm[j], outcome = outcomes,
             med = apply(m, 2, median, na.rm = TRUE),
             lo  = apply(m, 2, quantile, .025, na.rm = TRUE),
             hi  = apply(m, 2, quantile, .975, na.rm = TRUE), row.names = NULL)
}))
set_factors <- function(d, eff_levs) {
  d$size    <- factor(d$size, levels = names(size_R0))
  d$eff     <- factor(d$eff,  levels = eff_levs)
  d$arm     <- factor(d$arm,  levels = arm_levels)
  d$outcome <- factor(d$outcome,
                      levels = c("infections", "symptomatic", "hospitalisations", "deaths"),
                      labels = c("Infections", "Symptomatic", "Hospitalisations", "Deaths"))
  d
}

b6 <- set_factors(mc_long(seq_len(np), perm_mc), eff_x_levels)
# Complete the efficacy x mechanism grid with NA rows so regular position_dodge
# reserves the SAME 3 mechanism slots at every efficacy level. This keeps the bars
# and their error bars aligned -- position_dodge2 mis-places the error bars when the
# single "No vaccine" bar sits beside the 2-mechanism efficacy groups. NB merge on
# the already-labelled factors, then only RE-ORDER (do not relabel) the factors.
b6 <- merge(expand.grid(size = levels(b6$size), outcome = levels(b6$outcome),
                        eff = levels(b6$eff), arm = levels(b6$arm), stringsAsFactors = FALSE),
            b6, all.x = TRUE)
b6$size    <- factor(b6$size,    levels = names(size_R0))
b6$eff     <- factor(b6$eff,     levels = eff_x_levels)
b6$arm     <- factor(b6$arm,     levels = arm_levels)
b6$outcome <- factor(b6$outcome, levels = c("Infections", "Symptomatic",
                                            "Hospitalisations", "Deaths"))
dodge <- position_dodge(width = 0.8)

p6 <- ggplot(b6, aes(eff, med, fill = arm)) +
  geom_col(position = dodge, width = 0.75) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = 0.25, linewidth = 0.4) +
  ggh4x::facet_grid2(size ~ outcome, scales = "free_y", independent = "y") +
  scale_fill_manual(values = arm_cols, name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Total burden, median + 95% UI",
       title = "Hypothetical MAYV cases - total burden") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom", panel.grid.minor = element_blank())
p6
ggsave("ca_mayv_vacc_6perm_burden.png", p6, width = 11, height = 5.5, dpi = 120)
cat("Saved plot: ca_mayv_vacc_6perm_burden.png\n")

# ============================================================
# 9. Burden AVERTED (vs the same-R0 no-vaccine baseline)
# ============================================================
# No vaccine averts 0 by definition, so the averted view drops it: 2 sizes x 2
# efficacy levels x 2 mechanisms = 8 bars. Averted = (no-vaccine) - (scenario),
# PAIRED per draw (baseline and scenario share the same season/severity draw),
# reusing perm_mc from section 8.
base_idx <- sapply(seq_len(np), function(j)
  which(perm$size == perm$size[j] & perm$arm == "No vaccine"))
vac_rows <- which(perm$arm != "No vaccine")

av_pe <- t(sapply(vac_rows, function(j) pe[base_idx[j], ] - pe[j, ]))
rownames(av_pe) <- perm$key[vac_rows]
av_mc <- lapply(vac_rows, function(j) perm_mc[[base_idx[j]]] - perm_mc[[j]])

mk_av_tbl <- function(value_fun) data.frame(
  `Outbreak size`     = perm$size[vac_rows],
  `Vaccine efficacy`  = perm$eff[vac_rows],
  `Vaccine mechanism` = perm$arm[vac_rows],
  Infections       = value_fun("infections"),
  Symptomatic      = value_fun("symptomatic"),
  Hospitalisations = value_fun("hospitalisations"),
  Deaths           = value_fun("deaths"),
  check.names = FALSE, row.names = NULL
)
av_pe_tbl <- mk_av_tbl(function(o) round(av_pe[, o], dgt[[o]]))
av_ui_tbl <- mk_av_tbl(function(o) sapply(seq_along(vac_rows),
                                          function(k) fmtq(av_mc[[k]][, o], dgt[[o]])))

cat("\n=== Burden AVERTED vs no vaccine (median, 95% UI) ===\n")
old_width <- getOption("width"); options(width = 240)
print(av_ui_tbl, row.names = FALSE)
options(width = old_width)

write_xlsx(list(point_estimate = av_pe_tbl, median_95UI = av_ui_tbl),
           "caldas_mayv_vacc_6perm_averted.xlsx")
cat("\nWrote caldas_mayv_vacc_6perm_averted.xlsx (sheets: point_estimate, median_95UI)\n")

# ---- Plot: burden averted, same layout (x = efficacy level, fill = mechanism) ----
a6 <- set_factors(mc_long(vac_rows, av_mc), names(eff_lvl))
a6$arm <- droplevels(a6$arm)   # only the two vaccine mechanisms here (no "No vaccine")

p6_av <- ggplot(a6, aes(eff, med, fill = arm)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8),
                width = 0.25, linewidth = 0.4) +
  ggh4x::facet_grid2(size ~ outcome, scales = "free_y", independent = "y") +
  scale_fill_manual(values = arm_cols, name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Burden averted vs no vaccine, median + 95% UI",
       title = "Hypothetical MAYV cases - burden averted") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom", panel.grid.minor = element_blank())
p6_av
ggsave("ca_mayv_vacc_6perm_averted.png", p6_av, width = 11, height = 5.5, dpi = 120)
cat("Saved plot: ca_mayv_vacc_6perm_averted.png\n")

# ============================================================
# 10. R0 SENSITIVITY analysis (Caicedo 1.1 - 1.3) -- NOT the 95% CI
# ============================================================
# Separate from the uncertainty intervals above: this varies R0 across Caicedo's
# outside-Amazon range (1.1, 1.2, 1.3) to show how sensitive the Current-outbreak
# results are to the assumed R0. DETERMINISTIC point estimates (central severity
# inputs, central seasonal shape) -- no Monte Carlo, no CI. The Future outbreak
# (R0 = 3.51) is a single value (Dodero-Rojas upper bound) with no range, so it is
# not part of this sweep.
sens_R0 <- c(1.1, 1.2, 1.3)
sens <- do.call(rbind, lapply(sens_R0, function(r0) {
  b <- t(sapply(seq_len(nrow(mech)), function(m)
    burden(run_perm(r0, mech$coverage[m], mech$VE_inf[m], mech$VE_block[m]))))
  data.frame(R0 = r0,
             `Vaccine efficacy`  = mech$eff,
             `Vaccine mechanism` = mech$arm,
             Infections       = round(b[, "infections"]),
             Symptomatic      = round(b[, "symptomatic"]),
             Hospitalisations = round(b[, "hospitalisations"], 1),
             Deaths           = round(b[, "deaths"], 2),
             check.names = FALSE, row.names = NULL)
}))
cat("\n=== R0 SENSITIVITY: total burden at R0 = 1.1 / 1.2 / 1.3 (point estimates) ===\n")
old_width <- getOption("width"); options(width = 240)
print(sens, row.names = FALSE)
options(width = old_width)

write_xlsx(list(R0_sensitivity = sens), "caldas_mayv_vacc_R0_sensitivity.xlsx")
cat("\nWrote caldas_mayv_vacc_R0_sensitivity.xlsx (sheet: R0_sensitivity)\n")
