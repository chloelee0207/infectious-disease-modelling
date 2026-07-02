library(dplyr)
library(tidyr)
library(ggplot2)
library(splines)

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
# 2b. Scenario start weeks (week_index within the 2025-W23 -> 2026-W22 window)
# ------------------------------------------------------------
# Calendar -> index map for this fit window:
#   2025-W23 = index 1  => 2025 weeks: index = week - 22
#   2026-W01 = index 31 => 2026 weeks: index = 30 + week
#   peak (2026-W09) = index 39
week_to_index <- function(year, week) ifelse(year == 2025, week - 22L, 30L + week)

start_s1 <- week_to_index(2026, 16)   # IXCHIQ real rollout, 18 Apr 2026 (2026-W16) -> 46
start_s2 <- week_to_index(2026, 1)    # start of 2026 (2026-W01)                     -> 31
start_s3 <- week_to_index(2025, 26)   # middle of 2025, before the outbreak          -> 4
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
library(readxl)
dp <- read_excel("disease_progression.xlsx", sheet = "disease_progression")
names(dp)[1:10] <- c("parameter", "group", "median", "ui_lo", "ui_hi",
                     "dist", "p1", "alpha", "p2", "beta")
stopifnot(all(dp$dist == "Beta"))

# Pull alpha/beta for rows matching an exact Parameter (+ optional Group regex),
# ordered by the lower age bound when the Group is an "Age [lo, hi)" band.
get_beta_ab <- function(param, group_regex = NULL) {
  d <- dp[dp$parameter == param, ]
  if (!is.null(group_regex)) d <- d[grepl(group_regex, d$group), ]
  lo <- suppressWarnings(as.numeric(sub(".*\\[\\s*([0-9]+).*", "\\1", d$group)))
  if (!all(is.na(lo))) d <- d[order(lo), ]
  list(a = d$alpha, b = d$beta)
}

# prop_symp (symptomatic among infections), "Overall"
ps <- get_beta_ab("Probability of symptomatic cases among infections", "^Overall$")
ps_a <- ps$a; ps_b <- ps$b

# hospitalisation among symptomatic (single, not age-specific)
hp <- get_beta_ab("Probability of hospitalisation among symptomatic cases")
hosp_a <- hp$a; hosp_b <- hp$b
hosp_rate <- hosp_a / (hosp_a + hosp_b)             # mean (~0.040)

# case fatality by decadal band [0,10)..[80,90) (length 9), hospitalised vs not
ch <- get_beta_ab("Probability of death among hospitalised cases")
cn <- get_beta_ab("Probability of death among non-hospitalised cases")
cfr_hosp_a <- ch$a; cfr_hosp_b <- ch$b
cfr_nonh_a <- cn$a; cfr_nonh_b <- cn$b
stopifnot(length(ps_a) == 1, length(hosp_a) == 1,
          length(cfr_hosp_a) == 9, length(cfr_nonh_a) == 9)

cfr_h_mean <- cfr_hosp_a / (cfr_hosp_a + cfr_hosp_b)
cfr_n_mean <- cfr_nonh_a / (cfr_nonh_a + cfr_nonh_b)
cfr_band   <- hosp_rate * cfr_h_mean + (1 - hosp_rate) * cfr_n_mean  # death per symptomatic, by band

age_to_band <- c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9)   # 12 model groups -> 9 decadal bands
stopifnot(length(age_to_band) == A)
cfr_vec <- cfr_band[age_to_band]                       # death per symptomatic case, length A

# Burden extractor: infections, symptomatic, hospitalisations, deaths.
# hr (hosp rate) and cv (CFR-by-age) default to the means but accept per-draw values.
burden <- function(out, hr = hosp_rate, cv = cfr_vec) {
  symp_age <- rowSums(out$new_symptomatic)
  symp     <- sum(symp_age)
  c(infections       = sum(out$new_infections),
    symptomatic      = symp,
    hospitalisations = symp * hr,
    deaths           = sum(symp_age * cv))
}

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
# 6. Plots
# ============================================================
x_ticks <- caldas_obs |>
  filter((Year == 2025 & week %in% c(30, 40, 50)) |
           (Year == 2026 & week %in% c(10, 20)))
year_break <- mean(c(max(caldas_obs$week_index[caldas_obs$Year == 2025]),
                     min(caldas_obs$week_index[caldas_obs$Year == 2026])))

timing_cols <- c("actual rollout" = "#f4a582",
                 "start of 2026"  = "#2166ac",
                 "pre-outbreak"   = "#b2182b")

# (a) Weekly TRUE infections: baseline + the 3 timings (infection-blocking arm).
#     Disease-blocking-only is omitted here because it does not change infections.
inf_keys <- c("No vaccine (baseline)",
              paste0(names(timings), " | Disease + infection blocking"))
inf_df <- do.call(rbind, lapply(inf_keys, function(nm) {
  lbl <- if (nm == "No vaccine (baseline)") nm else scen_meta$timing[scen_meta$name == nm]
  data.frame(week = 1:T_weeks, infections = colSums(scenarios[[nm]]$new_infections), scenario = lbl)
}))
inf_df$scenario <- factor(inf_df$scenario,
                          levels = c("No vaccine (baseline)", names(timings)))

p_inf <- ggplot(inf_df, aes(week, infections, colour = scenario)) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey60") +
  # geom_vline(xintercept = unlist(timings), linetype = "dotted",
  #            colour = timing_cols[names(timings)], linewidth = 0.5) +
  annotate("text", x = year_break, y = 0, label = "2026", angle = 90,
           vjust = -0.5, hjust = -11.5, fontface = "bold", size = 3.5) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("No vaccine (baseline)" = "grey40", timing_cols), name = NULL) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week", y = "Weekly infections (true)",
       title = "Caldas Novas CHIKV cases - disease + infection blocking") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p_inf
ggsave("ca_vacc_infections.png", p_inf, width = 8, height = 4.5, dpi = 120)

# (b) Averted burden: facet by outcome, x = timing, fill = efficacy arm
bar_df <- summary_tbl |>
  dplyr::select(timing, arm, Infections = inf_averted, Symptomatic = symp_averted,
                Hospitalisations = hosp_averted, Deaths = deaths_averted) |>
  pivot_longer(c(Infections, Symptomatic, Hospitalisations, Deaths),
               names_to = "outcome", values_to = "averted")
bar_df$timing  <- factor(bar_df$timing, levels = names(timings))
bar_df$arm     <- factor(bar_df$arm, levels = names(arms))   # disease-blocking first
bar_df$outcome <- factor(bar_df$outcome,
                         levels = c("Infections", "Symptomatic", "Hospitalisations", "Deaths"))

p_bar <- ggplot(bar_df, aes(timing, averted, fill = arm)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ outcome, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Disease-blocking" = "#9ecae1",
                               "Disease + infection blocking" = "#08519c"),
                    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Averted (vs no vaccine)",
       title = "Caldas Novas CHIKV cases - burden averted") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom", panel.grid.minor = element_blank())
p_bar
ggsave("ca_vacc_averted.png", p_bar, width = 10, height = 4.2, dpi = 120)

cat("\nSaved plots: ca_vacc_infections.png, ca_vacc_averted.png\n")

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

  for (tn in names(timings)) for (an in names(arms)) {
    nm  <- paste0(tn, " | ", an)
    out <- run_scenario(total_coverage, timings[[tn]], arms[[an]]["VE_inf"],
                        arms[[an]]["VE_block"], base_beta = beta_t_i, prop_symp_use = ps_i)
    av_draws[[nm]][i, ] <- bb - burden(out, hosp_i, cfr_vec_i)
  }
  if (i %% 100 == 0) cat("  ", i, "/", n_draws, "\n")
}

# ---- Summarise: median (2.5% - 97.5%) ----
fmtq <- function(v, d = 0) {
  q <- quantile(v, c(.5, .025, .975), na.rm = TRUE)
  f <- function(x) formatC(round(x, d), big.mark = ",", format = "f", digits = d)
  sprintf("%s (%s - %s)", f(q[1]), f(q[2]), f(q[3]))
}

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
write.csv(mc_tbl, "caldas_vacc_averted_mc.csv", row.names = FALSE)
cat("\nWrote caldas_vacc_averted_mc.csv\n")

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
       title = "Caldas Novas CHIKV cases - burden averted") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom", panel.grid.minor = element_blank())
p_mc
ggsave("ca_vacc_averted_mc.png", p_mc, width = 10, height = 4.2, dpi = 120)
cat("Saved plot: ca_vacc_averted_mc.png\n")
