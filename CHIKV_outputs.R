# ============================================================
# Burden table via Monte Carlo propagation
# Combines: (a) calibration uncertainty in beta_t / rho  -> Hessian mvrnorm
#           (b) input-parameter uncertainty             -> Beta draws (your xlsx)
# Outputs: infections / symptomatic / deaths AVERTED, with 95% UI,
#          for disease-blocking and disease+infection-blocking.
#
# Assumes the following already exist from your fit script:
#   fit, vcov_mat, df_spline, make_beta_t(), seir_baseline(),
#   T_weeks, A, N, R_init_prop, I0, E0, sigma, gamma, best_rho
# and the vaccine scenario objects from CHIKV_se_vacc.R:
#   seirv_vaccinated(), target_age, weekly_delivery_speed, delay, immun_delay
# NOTE: seirv_vaccinated() now seeds E0 (exposed) so week-1 incidence matches the
# fit; E0 (the frozen back-calculated seed) must be in memory and is passed below.
# ============================================================

library(MASS)    # mvrnorm
library(readxl)  # read_excel
set.seed(2026)

n_draws <- 1000

# ------------------------------------------------------------
# 1. Input-parameter distributions — read straight from the Excel by name.
#    Edit the spreadsheet, not this script.
# ------------------------------------------------------------
dp <- read_excel("disease_progression.xlsx", sheet = "disease_progression")

# columns of interest (backticks because of spaces): `Value 1` = alpha, `Value 2` = beta
A_COL <- "Value 1"; B_COL <- "Value 2"

# verify a row's alpha/beta against its stated median (catches typos at point of use)
.check_row <- function(d) {
  implied <- d[[A_COL]] / (d[[A_COL]] + d[[B_COL]])
  bad <- abs(implied - d$Median) > 0.02
  if (any(bad)) warning("alpha/beta disagree with median for: ",
                        paste(d$Parameter[bad], d$Group[bad], sep = " / ", collapse = "; "))
}

# one (alpha, beta) for a parameter, optionally narrowed to a Group
get_beta <- function(param, group = NULL) {
  d <- dp[dp$Parameter == param, , drop = FALSE]
  if (!is.null(group)) d <- d[d$Group == group, , drop = FALSE]
  stopifnot(nrow(d) == 1)
  .check_row(d)
  c(alpha = d[[A_COL]], beta = d[[B_COL]])
}

# age-banded (alpha, beta) vectors, ordered youngest -> oldest by the [lower, ...) label
get_beta_age <- function(param) {
  d <- dp[dp$Parameter == param, , drop = FALSE]
  lower <- as.numeric(sub(".*\\[\\s*([0-9]+).*", "\\1", d$Group))  # "Age [60, 70)" -> 60
  d <- d[order(lower), , drop = FALSE]
  .check_row(d)
  list(alpha = d[[A_COL]], beta = d[[B_COL]])
}

# prop_symp: "Overall" matches the 0.5242478 you calibrated with.
# (Americas-specific Beta is also in the sheet if you prefer region-matching,
#  but keep it consistent with whatever prop_symp your beta_t fit assumed.)
ps     <- get_beta("Probability of symptomatic cases among infections", "Overall")
ps_a   <- ps["alpha"];  ps_b <- ps["beta"]

hosp   <- get_beta("Probability of hospitalisation among symptomatic cases")
hosp_a <- hosp["alpha"]; hosp_b <- hosp["beta"]

ch <- get_beta_age("Probability of death among hospitalised cases")
cfr_hosp_a <- ch$alpha;  cfr_hosp_b <- ch$beta            # length 9
cn <- get_beta_age("Probability of death among non-hospitalised cases")
cfr_nonh_a <- cn$alpha;  cfr_nonh_b <- cn$beta            # length 9

stopifnot(length(cfr_hosp_a) == 9, length(cfr_nonh_a) == 9)

# ------------------------------------------------------------
# 1b. Vaccine efficacy and coverage (sampled per draw, see loop)
# ------------------------------------------------------------
# Moment-match a Beta from a mean + 95% CI — the SAME method that produced the
# alpha/beta values in your spreadsheet.
beta_from_ci <- function(mean, lower, upper) {
  v <- ((upper - lower) / (2 * 1.96))^2          # variance implied by the CI
  s <- mean * (1 - mean) / v - 1                 # alpha + beta
  c(alpha = mean * s, beta = (1 - mean) * s)
}

# Ixchiq efficacy: 98.9% (95% CI 96.7-99.8)  -> Beta(~171, ~1.9)
ve   <- beta_from_ci(0.989, 0.967, 0.998)
ve_a <- ve["alpha"]; ve_b <- ve["beta"]

# Coverage: 20-40% of the 18-59 eligible population. Brazil's IXCHIQ rollout has
# extended beyond the original 10 municipalities and the doses to Sete Lagoas have
# not been publicly announced, so we model explicit target-coverage levels rather
# than back-calculating from a dose budget. We run the Monte Carlo at each FIXED
# coverage level (low and high bound), so each burden table's 95% UI reflects
# calibration + input-parameter uncertainty ONLY, at that coverage (Hyolim's
# approach, where coverage was fixed at 50%). This keeps the two coverage
# scenarios cleanly separated instead of blending them into one CI.
coverage_levels <- c(0.20, 0.40)

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
# 4. Reuse your existing functions — do NOT redefine them here.
#    Keep seir_baseline(), seirv_vaccinated() and make_beta_t() in ONE place.
#    If they live in their own definitions file, load it once:
# source("seir_functions.R")
#    (If you've already run those definitions in this R session, skip the
#     source() — they're already in memory.)
#
#    NOTE: seirv_vaccinated() must accept prop_symp as an argument and use it
#    for new_symptomatic (as seir_baseline does). If yours hardcodes prop_symp
#    internally, add it as an argument so the sampled value flows through.
# ------------------------------------------------------------

# ------------------------------------------------------------
# 5-6. Monte Carlo engine: run n_draws at a fixed coverage level.
#      Returns scalar burden draws (res) plus weekly symptomatic trajectories
#      (summed over age) for the three scenarios, for the ribbon plot.
#      set.seed() is reset inside, so every coverage level sees the same
#      calibration/parameter draws -> the 20% and 40% results are paired and
#      differ only because of coverage.
# ------------------------------------------------------------
run_mc <- function(cov_fixed, seed = 2026) {
  set.seed(seed)
  param_samples <- mvrnorm(n = n_draws, mu = fit$par, Sigma = vcov_mat)

  res <- data.frame(
    base_inf   = NA_real_,
    base_symp  = NA_real_,
    base_hosp  = NA_real_,
    base_death = NA_real_,
    inf_av_both  = NA_real_,                                # disblock averts 0 infections
    symp_av_disb = NA_real_,
    symp_av_both = NA_real_,
    hosp_av_disb = NA_real_,
    hosp_av_both = NA_real_,
    death_av_disb= NA_real_,
    death_av_both= NA_real_
  )[rep(1, n_draws), ]

  # Weekly symptomatic (summed over age), one row per draw, for the UI ribbons.
  wk_base <- wk_disb <- wk_both <- matrix(NA_real_, nrow = n_draws, ncol = T_weeks)

  for (i in 1:n_draws) {

    # --- (a) transmission from the Hessian draw ---
    coefs_i  <- param_samples[i, 1:df_spline]
    beta_t_i <- make_beta_t(coefs_i)
    if (any(!is.finite(beta_t_i)) || any(beta_t_i > 20) || any(beta_t_i < 1e-6)) next

    # --- (b) input-parameter draws ---
    ps_i    <- rbeta(1, ps_a, ps_b)
    hosp_i  <- rbeta(1, hosp_a, hosp_b)
    cfr_h_i <- rbeta(9, cfr_hosp_a, cfr_hosp_b)
    cfr_n_i <- rbeta(9, cfr_nonh_a, cfr_nonh_b)
    cfr_band <- hosp_i * cfr_h_i + (1 - hosp_i) * cfr_n_i  # death per symptomatic, by band
    cfr_vec  <- cfr_band[age_to_band]                      # expand to A age groups

    ve_i <- rbeta(1, ve_a, ve_b)        # vaccine efficacy this draw
    # NB: no coverage draw here - coverage is FIXED at cov_fixed for this run.

    # --- (c) run the three scenarios with this beta_t and prop_symp ---
    #     rho is irrelevant here (burden uses infections/symptomatic, not reported)
    #     but the function requires it, so pass best_rho.
    out_base <- seirv_vaccinated(
      T_weeks = T_weeks, A = A, N = N, R_init_prop = R_init_prop, I0 = I0, E0 = E0,
      base_beta = beta_t_i, sigma = sigma, gamma = gamma, rho = best_rho,
      prop_symp = ps_i,                              # <- sampled this draw
      target_age = target_age, total_coverage = 0,   # zero coverage = no vaccine
      weekly_delivery_speed = weekly_delivery_speed, delay = delay,
      VE_inf = 0, VE_block = 0, immun_delay = immun_delay
    )
    out_disb <- seirv_vaccinated(
      T_weeks = T_weeks, A = A, N = N, R_init_prop = R_init_prop, I0 = I0, E0 = E0,
      base_beta = beta_t_i, sigma = sigma, gamma = gamma, rho = best_rho,
      prop_symp = ps_i,
      target_age = target_age, total_coverage = cov_fixed,
      weekly_delivery_speed = weekly_delivery_speed, delay = delay,
      VE_inf = 0, VE_block = ve_i, immun_delay = immun_delay
    )
    out_both <- seirv_vaccinated(
      T_weeks = T_weeks, A = A, N = N, R_init_prop = R_init_prop, I0 = I0, E0 = E0,
      base_beta = beta_t_i, sigma = sigma, gamma = gamma, rho = best_rho,
      prop_symp = ps_i,
      target_age = target_age, total_coverage = cov_fixed,
      weekly_delivery_speed = weekly_delivery_speed, delay = delay,
      VE_inf = ve_i, VE_block = ve_i, immun_delay = immun_delay
    )

    # --- (d) burden + averted ---
    b  <- burden(out_base, cfr_vec)
    d  <- burden(out_disb, cfr_vec)
    bo <- burden(out_both, cfr_vec)

    # Hospitalisations = symptomatic * P(hospitalisation | symptomatic). hosp_i is
    # not age-specific, so it scales each scenario's symptomatic total directly.
    base_hosp <- b$symp  * hosp_i
    disb_hosp <- d$symp  * hosp_i
    both_hosp <- bo$symp * hosp_i

    res$base_inf[i]    <- b$inf
    res$base_symp[i]   <- b$symp
    res$base_hosp[i]   <- base_hosp
    res$base_death[i]  <- b$death
    res$inf_av_both[i]   <- b$inf   - bo$inf     # disease-blocking infections averted = 0
    res$symp_av_disb[i]  <- b$symp  - d$symp
    res$symp_av_both[i]  <- b$symp  - bo$symp
    res$hosp_av_disb[i]  <- base_hosp - disb_hosp
    res$hosp_av_both[i]  <- base_hosp - both_hosp
    res$death_av_disb[i] <- b$death - d$death
    res$death_av_both[i] <- b$death - bo$death

    # weekly symptomatic trajectories (summed over age) for the ribbons
    wk_base[i, ] <- colSums(out_base$new_symptomatic)
    wk_disb[i, ] <- colSums(out_disb$new_symptomatic)
    wk_both[i, ] <- colSums(out_both$new_symptomatic)

    if (i %% 100 == 0)
      cat("  coverage ", sprintf("%.0f%%", 100 * cov_fixed), ": ", i, "/", n_draws, "\n", sep = "")
  }

  list(res = res, wk_base = wk_base, wk_disb = wk_disb, wk_both = wk_both)
}

# Run both coverage bounds. Paired draws (same seed) -> directly comparable.
mc <- lapply(coverage_levels, run_mc)
names(mc) <- paste0(round(100 * coverage_levels), "%")

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

make_burden_table <- function(res) {
  data.frame(
    Outcome = c("Infections averted", "Symptomatic cases averted",
                "Hospitalisations averted", "Deaths averted",
                "% reduction (symptomatic)"),
    `Disease-blocking only` = c(
      "-",                                   # infections averted = 0 by construction
      fmt(res$symp_av_disb),
      fmt(res$hosp_av_disb),
      fmt(res$death_av_disb, 1),
      pct(res$symp_av_disb, res$base_symp)
    ),
    `Disease + infection blocking` = c(
      fmt(res$inf_av_both),
      fmt(res$symp_av_both),
      fmt(res$hosp_av_both),
      fmt(res$death_av_both, 1),
      pct(res$symp_av_both, res$base_symp)
    ),
    check.names = FALSE
  )
}

# One table per coverage level (baseline is identical across levels - same draws,
# coverage does not affect the no-vaccine scenario).
burden_tables <- lapply(mc, function(m) make_burden_table(m$res))
for (cl in names(mc)) {
  res <- mc[[cl]]$res
  cat("\n--- Sete Lagoas burden table @ ", cl,
      " coverage (median, 95% UI) ---\n", sep = "")
  print(burden_tables[[cl]], row.names = FALSE)
  cat("\nBaseline (no vaccine):\n")
  cat("  Infections:      ", fmt(res$base_inf),   "\n")
  cat("  Symptomatic:     ", fmt(res$base_symp),  "\n")
  cat("  Hospitalisations:", fmt(res$base_hosp),  "\n")
  cat("  Deaths:          ", fmt(res$base_death, 1), "\n")
}

# ------------------------------------------------------------
# 8. Diagnostics, per coverage level. Same format as the burden table:
#    median first, then (2.5%-97.5%), thousands-separated.
# ------------------------------------------------------------
print_diagnostics <- function(res, label) {
  # Recover scenario TOTALS as baseline - averted (the loop only stored averted).
  # Disease-blocking has VE_inf = 0, so it blocks no infections -> its infection
  # total equals baseline.
  disb_inf   <- res$base_inf                         # = baseline (no infections averted)
  disb_symp  <- res$base_symp  - res$symp_av_disb
  disb_hosp  <- res$base_hosp  - res$hosp_av_disb
  disb_death <- res$base_death - res$death_av_disb
  both_inf   <- res$base_inf   - res$inf_av_both
  both_symp  <- res$base_symp  - res$symp_av_both
  both_hosp  <- res$base_hosp  - res$hosp_av_both
  both_death <- res$base_death - res$death_av_both

  # one line: "<label>  median (2.5%-97.5%)" using the burden-table formatters
  row  <- function(name, x, d = 0) cat(sprintf("  %-30s %s\n", name, fmt(x, d)))
  rowp <- function(name, av, base) cat(sprintf("  %-30s %s\n", name, pct(av, base)))

  cat("\n========= RAW NUMBERS @ ", label, " (median, 95% UI) =========\n", sep = "")
  cat("Baseline:\n")
  row("infections",       res$base_inf)
  row("symptomatic",      res$base_symp)
  row("hospitalisations", res$base_hosp)
  row("deaths",           res$base_death, 1)

  cat("\nDisease-blocking only:\n")
  row("infections total (= baseline)", disb_inf)
  row("symptomatic total",            disb_symp)
  row("hospitalisations total",       disb_hosp)
  row("deaths total",                 disb_death, 1)
  row("symptomatic averted",          res$symp_av_disb)
  row("hospitalisations averted",     res$hosp_av_disb)
  row("deaths averted",               res$death_av_disb, 1)
  rowp("% reduction (symptomatic)",   res$symp_av_disb,  res$base_symp)
  rowp("% reduction (deaths)",        res$death_av_disb, res$base_death)

  cat("\nDisease + infection blocking:\n")
  row("infections total",         both_inf)
  row("symptomatic total",        both_symp)
  row("hospitalisations total",   both_hosp)
  row("deaths total",             both_death, 1)
  row("infections averted",       res$inf_av_both)
  row("symptomatic averted",      res$symp_av_both)
  row("hospitalisations averted", res$hosp_av_both)
  row("deaths averted",           res$death_av_both, 1)
  rowp("% reduction (symptomatic)", res$symp_av_both,  res$base_symp)
  rowp("% reduction (deaths)",      res$death_av_both, res$base_death)
  cat("\n")
}

for (cl in names(mc)) print_diagnostics(mc[[cl]]$res, paste0(cl, " coverage"))

# Scale check: reported = rho * symptomatic, so this should match the fitted
# predicted reported total from CHIKV_se_pre_vacc_optim.R (sum of predicted_cases).
sum(colSums(out_best$new_reported))
# Or equivalently
best_rho * sum(out_best$new_symptomatic)

# ------------------------------------------------------------
# 9. Epidemic-curve plot with 95% UI ribbons (cf. Hyolim Fig 1A).
#    One SEPARATE plot per coverage level (18-59 / Ixchiq).
#    y = predicted symptomatic cases; black x = reported (symptomatic) cases.
# ------------------------------------------------------------
library(ggplot2)

# weekly 2.5 / 50 / 97.5% across draws for one scenario matrix (n_draws x T_weeks)
ribbon_df <- function(wk, scenario) {
  qs <- apply(wk, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  data.frame(week = 1:T_weeks, lo = qs[1, ], med = qs[2, ], hi = qs[3, ],
             scenario = scenario)
}

scen_cols <- c("No vaccination"               = "grey60",
               "Disease-blocking only"        = "#4393c3",
               "Disease + infection blocking" = "#d6604d")

# Build one ribbon plot for a single coverage level.
make_ribbon_plot <- function(cl) {
  m <- mc[[cl]]
  plot_df <- rbind(
    ribbon_df(m$wk_base, "No vaccination"),
    ribbon_df(m$wk_disb, "Disease-blocking only"),
    ribbon_df(m$wk_both, "Disease + infection blocking")
  )
  plot_df$scenario <- factor(plot_df$scenario,
    levels = c("No vaccination", "Disease-blocking only", "Disease + infection blocking"))

  lab <- paste0(
    "% reduction (symptomatic):\n",
    "Disease-blocking only: ",        pct(m$res$symp_av_disb, m$res$base_symp), "\n",
    "Disease + infection blocking: ", pct(m$res$symp_av_both, m$res$base_symp))

  # Vaccine deployment (rollout) window: doses delivered from `delay` over
  # 1/weekly_delivery_speed weeks (e.g. weeks 2-12 with delay=2, speed=0.10).
  roll_start <- delay
  roll_end   <- delay + 1 / weekly_delivery_speed

  ggplot(plot_df, aes(week, med, colour = scenario, fill = scenario)) +
    annotate("rect", xmin = roll_start, xmax = roll_end, ymin = -Inf, ymax = Inf,
             fill = "#3a7d3a", alpha = 0.10) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 0.8) +
    # geom_point(data = data.frame(week = 1:T_weeks, med = observed_cases),
    #            aes(week, med), colour = "black", shape = 4, size = 1.1,
    #            inherit.aes = FALSE) +
    annotate("text", x = Inf, y = Inf, label = lab,
             hjust = 1.02, vjust = 1.2, size = 2.8, lineheight = 0.95) +
    scale_colour_manual(values = scen_cols, aesthetics = c("colour", "fill")) +
    scale_x_continuous(breaks = seq(0, 52, 10)) +
    labs(x = "Week", y = "Predicted symptomatic CHIKV cases", colour = NULL, fill = NULL,
         title = paste0("Sete Lagoas 2024: symptomatic cases at ", cl,
                        " vaccine coverage"),
         caption = sprintf("Shaded band = vaccine deployment window (weeks %g-%g)",
                           roll_start, roll_end)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

# One plot object per coverage level, printed separately (each on its own device/tab).
ribbon_plots <- lapply(names(mc), make_ribbon_plot)
names(ribbon_plots) <- names(mc)
for (cl in names(ribbon_plots)) print(ribbon_plots[[cl]])

# Access individually if needed, e.g.:  ribbon_plots[["20%"]]  /  ribbon_plots[["40%"]]