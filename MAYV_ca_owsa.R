# ============================================================
# MAYV_ca_owsa.R -- Caldas Novas MAYV one-way sensitivity analysis.
#
# DETERMINISTIC, at the scenario's FIXED peak R0. MAYV is not fitted to data, so no
# re-fitting is needed: every run is a forward simulation at chosen parameter values.
# No take-off conditioning: R0 is fixed, so there is no fizzle/take-off split.
#
# BASE CASE. Peak R0 is taken from the engine's fixed scenario value (M$R0_fixed --
# low = 1.20 Caicedo outside-Amazon, high = 2.50 sustained urban Aedes transmission).
# R0 is NOT a tornado row: it is scenario-defining, not a single-setting uncertain
# parameter, and outbreak size is a steep convex function of it, so one bar would swamp
# every other parameter and would misrepresent between-setting heterogeneity as
# uncertainty. Instead R0 gets its own DEDICATED RESPONSE CURVE (section 5b): a fine
# sweep over R0_SWEEP, reported as a table + figure. That is the honest way to show a
# threshold relationship -- a curve, not an interval.
#
# The tornado therefore answers: GIVEN a MAYV outbreak at the scenario's R0, what drives
# the vaccine's impact.
#
# The vaccine is DISEASE-BLOCKING ONLY (Ixchiq cross-protection against MAYV is
# hypothetical), so infections are identical across arms and only symptomatic cases,
# hospitalisations and DALYs move.
#
# The inherited seasonal envelope (hybrid: CHIKV beta rise + climatological dry-season
# tail) is reported as a structural sensitivity against its two pure components, not as
# a tornado row.
#
# Run order:  MAYV_ca_lhs.R -> MAYV_ca_engine.R (once) -> this file
# ============================================================
suppressMessages({library(dplyr); library(ggplot2); library(writexl)})

DEFS_ONLY <- TRUE
source("MAYV_ca_engine.R")        # seirv_vaccinated_MAYV, outcome_one, severity/DALY params
E <- readRDS("MAYV_ca_lhs_ensemble.rds")
M <- readRDS("MAYV_ca_engine_results.rds")

# ------------------------------------------------------------
# 1. Central values and one-way bounds
# ------------------------------------------------------------
ob <- M$outbreak                                   # ALL draws (fixed R0 -> no conditioning)
R0_BASE <- unname(M$R0_fixed)                      # the scenario's FIXED peak R0
stopifnot(is.finite(R0_BASE))
# Dedicated R0 response curve (section 5b). Step 0.1 is deliberate: in the steep region a
# 0.1 increment roughly doubles outbreak size, so a coarser grid would hide the threshold.
R0_SWEEP <- seq(1.1, 3.6, by = 0.1)                # 26 points spanning both scenarios

BASE <- list(R0 = R0_BASE, imm = 0.0735, rho = 0.25, ve = 0.50,
             cov = 0.30, deliv = 0.10, delay = 2, immun = 2, env = "hybrid")

# TWO inputs are deliberately NOT tornado rows, because neither is a parameter with an
# ordered uncertainty range; both are structural choices reported in the
# structural_sensitivity sheet instead:
#   R0       -- defines the SCENARIO (whether an outbreak occurs and how large). It is
#               swept as a response curve in section 5b instead of getting a tornado bar.
#   envelope -- a model-structure choice between the pure CHIKV beta shape, the pure
#               rainfall shape and the hybrid. "Lower" and "upper" are meaningless for
#               a categorical choice, and the rainfall arm additionally peaks before the
#               campaign completes, so a single bar would conflate size with timing.
# The reporting rate is also excluded, but for a different reason: MAYV is not fitted,
# so rho never enters the transmission model. The epidemic size comes from R0, prior
# immunity and the seed; rho is applied afterwards only to convert TRUE cases into
# REPORTED ones. It therefore has exactly zero effect on true cases averted.
BOUNDS <- list(
  imm   = c(0.03, 0.18),                               # Lima 2021 Central-West 95% CI
  ve    = c(0.25, 0.75),                               # hypothetical cross-protection
  cov   = c(0.20, 0.40),
  deliv = c(0.09, 0.11),
  delay = c(1, 3),
  immun = c(1, 3))

PAR_LAB <- c(R0 = "Peak R0", imm = "Rate of exposure",
             rho = "Reporting rate", ve = "Vaccine efficacy",
             cov = "Vaccine coverage", deliv = "Weekly delivery speed",
             delay = "Delay in deployment", immun = "Time to immunity",
             env = "Seasonal envelope (inherited)")

GAMMA <- E$base_gamma; SIGMA <- E$base_sigma; PSYMP <- 0.5242478
EVAL_WIN <- 1:T_weeks

# Seasonal envelopes: hybrid (base) vs its two pure components, each peak-normalised
env_of <- function(tag) {
  v <- switch(tag,
    hybrid = readRDS("caldas_hybrid_season.rds"),
    beta   = readRDS("caldas_beta_season.rds"),
    rain   = { r <- readRDS("caldas_rain_season.rds"); r / mean(r) })
  v <- as.numeric(v); stopifnot(length(v) == T_weeks)
  v / max(v)                                       # r0_is_peak = TRUE
}

# ------------------------------------------------------------
# 2. One scenario: baseline + disease-blocking arm, deterministic
# ------------------------------------------------------------
run_scenario <- function(p) {
  season <- env_of(p$env)
  Rimm   <- rep(p$imm, A)                          # flat seroprevalence (Route B)
  sus    <- N * (1 - p$imm)
  I0i    <- E$I0_total * sus / sum(sus)
  beta_t <- p$R0 * GAMMA * season
  st     <- min(E$seed_week + 0, T_weeks)
  sim <- function(cov, vb)
    seirv_vaccinated_MAYV(T_weeks, A, N, Rimm, I0i, beta_t, SIGMA, GAMMA, p$rho,
      target_age, cov, p$deliv, start_pre + p$delay, VE_inf = 0, VE_block = vb,
      immun_delay = p$immun, prop_symp = PSYMP, E0 = E$E0, seed_week = E$seed_week)
  b <- sim(0, 0); v <- sim(p$cov, p$ve)
  psi <- PSYMP
  symp_b <- sum(colSums(psi * b$new_infections)[EVAL_WIN])
  symp_v <- sum(colSums(psi * v$new_infections * (1 - p$ve * v$coverage_frac))[EVAL_WIN])
  inf    <- sum(b$new_infections[, EVAL_WIN])
  pool   <- sum(N * (1 - p$imm))
  list(symptomatic_base = symp_b, averted = symp_b - symp_v,
       attack = 100 * inf / pool,
       hosp_averted = (symp_b - symp_v) * hosp_rate)
}

# ------------------------------------------------------------
# 3. Base case + each parameter at both bounds
# ------------------------------------------------------------
cat(sprintf("MAYV OWSA: base peak R0 = %.2f (FIXED, '%s' scenario; %d draws, no conditioning)\n",
            R0_BASE, M$R0_scenario, length(ob)))
base_run <- run_scenario(BASE)
cat(sprintf("  base: attack %.1f%% | symptomatic %.0f | averted %.0f\n",
            base_run$attack, base_run$symptomatic_base, base_run$averted))

rows <- list()
for (nm in names(BOUNDS)) {
  for (side in c("lower", "upper")) {
    p <- BASE; p[[nm]] <- BOUNDS[[nm]][if (side == "lower") 1 else 2]
    r <- run_scenario(p)
    rows[[length(rows)+1]] <- data.frame(
      parameter = unname(PAR_LAB[nm]), par_key = nm, bound = side,
      value = as.character(p[[nm]]),
      symptomatic_base = r$symptomatic_base, averted = r$averted,
      hosp_averted = r$hosp_averted, attack_pct = r$attack,
      pct_reduction = 100 * r$averted / r$symptomatic_base, row.names = NULL)
  }
  cat("  ", nm, " done\n", sep = "")
}
owsa <- do.call(rbind, rows)

# ---- structural sensitivities (NOT tornado rows): R0 scenario + seasonal envelope ----
struct <- rbind(
  # Both fixed R0 scenarios, for context (the full curve is in the r0_response sheet).
  do.call(rbind, lapply(c(1.20, 2.50), function(x) {
    r <- run_scenario(modifyList(BASE, list(R0 = x)))
    data.frame(input = "Peak R0 (scenario, FIXED)", setting = sprintf("%.2f", x),
               symptomatic_base = r$symptomatic_base, averted = r$averted,
               attack_pct = r$attack, pct_reduction = 100*r$averted/r$symptomatic_base) })),
  do.call(rbind, lapply(c("beta", "hybrid", "rain"), function(x) {
    r <- run_scenario(modifyList(BASE, list(env = x)))
    data.frame(input = "Seasonal envelope (structure)",
               setting = c(beta = "pure CHIKV beta", hybrid = "hybrid (base case)",
                           rain = "pure rainfall")[[x]],
               symptomatic_base = r$symptomatic_base, averted = r$averted,
               attack_pct = r$attack, pct_reduction = 100*r$averted/r$symptomatic_base) })))
rownames(struct) <- NULL

# ------------------------------------------------------------
# 4. Tornado
# ------------------------------------------------------------
b0 <- base_run$averted
sw <- owsa |> group_by(parameter) |> summarise(swing = diff(range(averted)), .groups = "drop")
owsa$parameter <- factor(owsa$parameter, levels = sw$parameter[order(sw$swing)])
p_t <- ggplot(owsa, aes(y = parameter)) +
  geom_vline(xintercept = b0, linetype = "dashed", colour = "grey45") +
  geom_segment(aes(x = b0, xend = averted, yend = parameter, colour = bound),
               linewidth = 6, alpha = .85) +
  scale_colour_manual(values = c(lower = "#d6604d", upper = "#4393c3"),
                      labels = c(lower = "Lower", upper = "Upper"), name = "Bound") +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Symptomatic cases averted", y = NULL) +
  theme_bw(11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("MAYV_ca_owsa_symptomatic.png", p_t, width = 8.5, height = 4.6, dpi = 130)
print(p_t)

# ------------------------------------------------------------
# 4b. R0 RESPONSE CURVE (replaces giving R0 a tornado bar).
# Outbreak size is a steep CONVEX function of peak R0: over R0_SWEEP it spans several
# orders of magnitude, so a single "lower/upper" bar would be meaningless and would
# swamp every vaccine parameter. A curve is the honest presentation. Deterministic at
# median inputs, matching the rest of this script.
# Note what the curve shows: the vaccine's RELATIVE effect (% symptomatic reduced) is
# almost R0-INVARIANT -- disease-blocking scales symptomatic cases by (1 - VE*coverage)
# whatever the epidemic size -- while the ABSOLUTE cases averted track outbreak size.
# ------------------------------------------------------------
r0_response <- do.call(rbind, lapply(R0_SWEEP, function(x) {
  r <- run_scenario(modifyList(BASE, list(R0 = x)))
  data.frame(R0 = x, symptomatic_base = r$symptomatic_base, averted = r$averted,
             attack_pct = r$attack,
             pct_reduction = 100 * r$averted / max(r$symptomatic_base, 1e-9),
             row.names = NULL) }))
cat("\n=== R0 response curve (deterministic, median inputs) ===\n")
print(transform(r0_response, symptomatic_base = round(symptomatic_base),
                averted = round(averted), attack_pct = round(attack_pct, 3),
                pct_reduction = round(pct_reduction, 1)), row.names = FALSE)

scen_pts <- subset(r0_response, abs(R0 - 1.20) < 1e-9 | abs(R0 - 2.50) < 1e-9)
p_r0 <- ggplot(r0_response, aes(R0)) +
  geom_line(aes(y = pmax(symptomatic_base, 0.5), colour = "Baseline symptomatic"), linewidth = 1) +
  geom_line(aes(y = pmax(averted, 0.5),          colour = "Averted by vaccine"),   linewidth = 1) +
  geom_point(data = scen_pts, aes(y = pmax(symptomatic_base, 0.5)), size = 2.6, colour = "#c0392b") +
  geom_vline(xintercept = c(1.20, 2.50), linetype = "dotted", colour = "grey40") +
  annotate("text", x = 1.20, y = Inf, label = "low\n1.20",  vjust = 1.3, size = 3, colour = "grey30") +
  annotate("text", x = 2.50, y = Inf, label = "high\n2.50", vjust = 1.3, size = 3, colour = "grey30") +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = c("Baseline symptomatic" = "#c0392b",
                                 "Averted by vaccine"   = "#2c7fb8"), name = NULL) +
  labs(x = "Wet-season peak R0 (fixed)", y = "Symptomatic cases (log scale)",
       title = "MAYV: outbreak size and vaccine impact vs fixed peak R0",
       subtitle = sprintf("Deterministic, median inputs; %.0f%% coverage, %.0f%% disease-blocking VE. Dotted lines = the two reported scenarios.",
                          100*BASE$cov, 100*BASE$ve)) +
  theme_bw(11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("MAYV_ca_r0_response.png", p_r0, width = 8.5, height = 5, dpi = 130)
cat("Saved MAYV_ca_r0_response.png\n")

# ------------------------------------------------------------
# 4c. PROBABILISTIC R0 SWEEP -- median + 95% UI at each fixed peak R0.
# The section-4b curve is deterministic (median inputs). This one re-runs the FULL
# ensemble at every R0 on the grid, so each row is a proper median [95% UI] carrying the
# same uncertainty the engine propagates: gamma, sigma, rho, prop_symp and prior immunity
# from the LHS ensemble, plus the engine's OWN sampled coverage / VE / deployment delay
# (reused draw-for-draw from MAYV_ca_engine_results.rds). Because those are the identical
# draws, the R0 = 2.50 row reproduces the engine's headline numbers.
#
# Weekly delivery speed is held at its median rather than resampled: it is the weakest
# input in the tornado (swing ~0.07 symptomatic cases) and the engine does not store its
# per-draw vector. Everything else is paired exactly.
#
# Cost: length(R0_SWEEP) x N_DRAWS SEIRV runs (~26,000, about 2 minutes). Set
# R0_PSA_RUN <- FALSE before sourcing to skip.
# ------------------------------------------------------------
if (!exists("R0_PSA_RUN")) R0_PSA_RUN <- TRUE
r0_response_psa <- NULL
if (R0_PSA_RUN) {
  nd      <- length(E$gamma)
  cov_v   <- as.numeric(M$cov_d)      # engine's sampled coverage of eligible 18-59
  ve_v    <- as.numeric(M$veb_d)      # engine's sampled disease-blocking VE
  delay_v <- as.numeric(M$delay_d)    # engine's sampled deployment delay (wk)
  deliv_c <- M$deliv_median           # weekly delivery speed, held at median (see above)
  stopifnot(length(cov_v) == nd, length(ve_v) == nd, length(delay_v) == nd)
  season_psa <- env_of(BASE$env)
  q3 <- function(x) unname(quantile(x, c(.5, .025, .975), na.rm = TRUE))

  cat(sprintf("\n4c. Probabilistic R0 sweep: %d R0 values x %d draws = %d runs...\n",
              length(R0_SWEEP), nd, length(R0_SWEEP)*nd))
  t_start <- Sys.time()
  r0_response_psa <- do.call(rbind, lapply(seq_along(R0_SWEEP), function(k) {
    x <- R0_SWEEP[k]
    sb <- av <- at <- numeric(nd)
    for (i in seq_len(nd)) {
      Rimm <- rep(E$immune_frac[i], A)
      sus  <- N * (1 - E$immune_frac[i])
      I0i  <- E$I0_total * sus / sum(sus)
      # VE_inf = 0 -> infections are vaccine-invariant, so ONE run yields both arms.
      r <- seirv_vaccinated_MAYV(T_weeks, A, N, Rimm, I0i, x * E$gamma[i] * season_psa,
             E$sigma[i], E$gamma[i], E$rho[i], target_age, cov_v[i], deliv_c,
             start_pre + delay_v[i], VE_inf = 0, VE_block = ve_v[i],
             immun_delay = BASE$immun, prop_symp = E$prop_symp[i],
             E0 = E$E0, seed_week = E$seed_week)
      ninf <- r$new_infections; psi <- E$prop_symp[i]
      sb[i] <- sum(colSums(psi * ninf)[EVAL_WIN])
      av[i] <- sb[i] - sum(colSums(psi * ninf * (1 - ve_v[i] * r$coverage_frac))[EVAL_WIN])
      at[i] <- 100 * sum(ninf[, EVAL_WIN]) / sum(sus)
    }
    pr <- 100 * av / pmax(sb, 1e-12)
    qs <- q3(sb); qa <- q3(av); qp <- q3(pr); qt <- q3(at)
    cat(sprintf("   R0 %.1f done (%2d/%2d)  symptomatic %.0f [%.0f, %.0f]\n",
                x, k, length(R0_SWEEP), qs[1], qs[2], qs[3]))
    data.frame(R0 = x,
               symptomatic_med = qs[1], symptomatic_lo = qs[2], symptomatic_hi = qs[3],
               averted_med = qa[1], averted_lo = qa[2], averted_hi = qa[3],
               pct_reduction_med = qp[1], pct_reduction_lo = qp[2], pct_reduction_hi = qp[3],
               attack_med = qt[1], attack_lo = qt[2], attack_hi = qt[3],
               row.names = NULL) }))
  cat(sprintf("   sweep finished in %.1f min\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

  # formatted "median (95% UI)" table -- the presentation version
  fmt3 <- function(m, lo, hi, d = 0)
    sprintf(paste0("%.", d, "f (%.", d, "f - %.", d, "f)"), m, lo, hi)
  r0_table <- data.frame(
    `Peak R0`                       = sprintf("%.1f", r0_response_psa$R0),
    `Attack rate (%)`               = fmt3(r0_response_psa$attack_med, r0_response_psa$attack_lo,
                                           r0_response_psa$attack_hi, 2),
    `Baseline symptomatic cases`    = fmt3(r0_response_psa$symptomatic_med, r0_response_psa$symptomatic_lo,
                                           r0_response_psa$symptomatic_hi, 0),
    `Symptomatic cases averted`     = fmt3(r0_response_psa$averted_med, r0_response_psa$averted_lo,
                                           r0_response_psa$averted_hi, 0),
    `Symptomatic cases reduced (%)` = fmt3(r0_response_psa$pct_reduction_med, r0_response_psa$pct_reduction_lo,
                                           r0_response_psa$pct_reduction_hi, 1),
    check.names = FALSE, stringsAsFactors = FALSE)
  write.csv(r0_table, "MAYV_ca_r0_response_table.csv", row.names = FALSE)
  cat("\n=== R0 response TABLE: median (95% UI) over ", nd, " draws ===\n", sep = "")
  print(r0_table, row.names = FALSE)

  # figure with 95% UI ribbons (log scale; the two reported scenarios marked)
  p_r0p <- ggplot(r0_response_psa, aes(R0)) +
    geom_ribbon(aes(ymin = pmax(symptomatic_lo, .5), ymax = pmax(symptomatic_hi, .5),
                    fill = "Baseline symptomatic"), alpha = .22) +
    geom_ribbon(aes(ymin = pmax(averted_lo, .5), ymax = pmax(averted_hi, .5),
                    fill = "Averted by vaccine"), alpha = .22) +
    geom_line(aes(y = pmax(symptomatic_med, .5), colour = "Baseline symptomatic"), linewidth = 1) +
    geom_line(aes(y = pmax(averted_med, .5),     colour = "Averted by vaccine"),   linewidth = 1) +
    geom_vline(xintercept = c(1.20, 2.50), linetype = "dotted", colour = "grey40") +
    annotate("text", x = 1.20, y = Inf, label = "low\n1.20",  vjust = 1.3, size = 3, colour = "grey30") +
    annotate("text", x = 2.50, y = Inf, label = "high\n2.50", vjust = 1.3, size = 3, colour = "grey30") +
    scale_y_log10(labels = scales::comma) +
    scale_colour_manual(values = c("Baseline symptomatic" = "#c0392b",
                                   "Averted by vaccine"   = "#2c7fb8"), name = NULL,
                        aesthetics = c("colour", "fill")) +
    labs(x = "Wet-season peak R0 (fixed)", y = "Symptomatic cases (log scale)",
         title = "MAYV: outbreak size and vaccine impact vs fixed peak R0",
         subtitle = sprintf(paste("Median and 95%% UI over %d draws: gamma, sigma, rho, prop_symp, prior immunity,",
                                  "coverage, VE, deployment delay.\nR0 is FIXED at each point, so its span is NOT inside these bands."), nd)) +
    theme_bw(11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                         plot.subtitle = element_text(size = 8.5))
  ggsave("MAYV_ca_r0_response_psa.png", p_r0p, width = 9, height = 5.2, dpi = 130)
  cat("Saved MAYV_ca_r0_response_psa.png and MAYV_ca_r0_response_table.csv\n")
}

# ------------------------------------------------------------
# 5. Export
# ------------------------------------------------------------
notes <- data.frame(item = c("Analysis", "Base case R0", "R0 response curve", "Vaccine",
                             "Seasonal envelope", "Fixed (not varied)", "Window"),
  detail = c(
  "Deterministic one-way sensitivity; MAYV is not fitted, so every run is a forward simulation.",
  sprintf("Peak R0 %.2f, FIXED by the '%s' scenario (low = 1.20 Caicedo et al. 2021 outside-Amazon-basin; high = 2.50 sustained urban Aedes transmission, conservative vs the 2.67 peak R0 fitted for CHIKV in this municipality). Not sampled: published MAYV R0 figures are point estimates from different settings and decades.",
          R0_BASE, M$R0_scenario),
  sprintf("R0 is NOT a tornado row. It is scenario-defining, and outbreak size is a steep CONVEX function of it (%.0f -> %.0f symptomatic across R0 %.1f-%.1f), so one bar would swamp every other parameter. It is swept as a response curve instead: see the r0_response sheet and MAYV_ca_r0_response.png. No take-off conditioning is applied anywhere (R0 fixed -> one transmission regime, unimodal outbreak size).",
          min(r0_response$symptomatic_base), max(r0_response$symptomatic_base),
          min(R0_SWEEP), max(R0_SWEEP)),
  "Disease-blocking only (VE_inf = 0): infections are identical across arms; only symptomatic cases and downstream outcomes move.",
  "Inherited from the CHIKV fit (hybrid = CHIKV beta rise + climatological dry-season tail). Reported as a structural sensitivity. NB the pure-rainfall envelope peaks at week 16, before the campaign completes (~week 29), so its low averted burden reflects timing as well as epidemic size.",
  "gamma, sigma, prop_symp (as in the CHIKV OWSA). Also the reporting rate: MAYV is not fitted, so rho does not enter the transmission model at all -- it converts true cases to reported cases after the simulation and has exactly zero effect on true burden averted. (In the CHIKV OWSA rho ranks third, because there it works backwards from the observed 8,204 cases to infer the size of the true epidemic.)",
  sprintf("52 weeks, 2025-W24 -> 2026-W22 (indices %d-%d).", min(EVAL_WIN), max(EVAL_WIN))),
  stringsAsFactors = FALSE)
write_xlsx(c(list(notes = notes,
                base_case = data.frame(symptomatic_base = base_run$symptomatic_base,
                                       averted = base_run$averted, attack_pct = base_run$attack,
                                       pct_reduction = 100*base_run$averted/base_run$symptomatic_base),
                owsa = owsa, structural_sensitivity = struct,
                r0_response = r0_response),
             if (!is.null(r0_response_psa))
               list(r0_response_psa = r0_response_psa, r0_response_table = r0_table)),
           "MAYV_ca_owsa.xlsx")
saveRDS(list(owsa = owsa, struct = struct, base = base_run, BASE = BASE, BOUNDS = BOUNDS,
             R0_BASE = R0_BASE, R0_SWEEP = R0_SWEEP, r0_response = r0_response,
             r0_response_psa = r0_response_psa,
             R0_scenario = M$R0_scenario), "MAYV_ca_owsa.rds")

cat("\n=== Swing in symptomatic cases averted ===\n")
print(owsa |> group_by(parameter) |>
        summarise(low = min(averted), high = max(averted), swing = high - low, .groups = "drop") |>
        arrange(desc(swing)) |> as.data.frame(), row.names = FALSE, digits = 5)
cat("\n=== Structural sensitivities (reported separately, not tornado rows) ===\n")
print(struct, row.names = FALSE, digits = 5)
cat("\nWrote MAYV_ca_owsa.xlsx, MAYV_ca_owsa.rds, MAYV_ca_owsa_symptomatic.png\n")

# ------------------------------------------------------------
# 6. Combined figure: CHIKV (panel A) over MAYV (panel B).
# MAYV has only the disease-blocking arm, so panel B is a single facet occupying the
# left-hand column; patchwork's grid layout keeps it aligned with CHIKV's
# disease-blocking facet above, with the right-hand cell left empty.
# ------------------------------------------------------------
if (!file.exists("CHIKV_ca_owsa.rds")) {
  cat("Skipped the combined tornado (CHIKV_ca_owsa.rds not found -- run CHIKV_ca_owsa.R).\n")
} else {
  suppressMessages(library(patchwork))
  ARMS <- c("Disease-blocking", "Disease + infection blocking")
  CO <- readRDS("CHIKV_ca_owsa.rds")

  # ---- panel A: CHIKV, both arms
  ca <- CO$owsa; ca$val <- ca$symptomatic
  bA <- setNames(CO$base$symptomatic, CO$base$arm)
  ca$base <- bA[ca$arm]; ca$arm <- factor(ca$arm, levels = ARMS)
  swA <- ca |> group_by(parameter) |> summarise(s = max(abs(val - base)), .groups = "drop")
  ca$parameter <- factor(ca$parameter, levels = swA$parameter[order(swA$s)])
  blA <- data.frame(arm = factor(names(bA), levels = ARMS), base = as.numeric(bA))

  # ---- panel B: MAYV, disease-blocking only (single facet, no empty panel)
  mb <- owsa; mb$val <- mb$averted; mb$base <- base_run$averted
  mb$arm <- factor(ARMS[1], levels = ARMS[1])
  swB <- mb |> group_by(parameter) |> summarise(s = max(abs(val - base)), .groups = "drop")
  mb$parameter <- factor(mb$parameter, levels = swB$parameter[order(swB$s)])
  blB <- data.frame(arm = factor(ARMS[1], levels = ARMS[1]), base = base_run$averted)

  tor <- function(d, bl, ttl, legend) {
    ggplot(d, aes(y = parameter)) +
      geom_vline(data = bl, aes(xintercept = base), linetype = "dashed", colour = "grey45") +
      geom_segment(aes(x = base, xend = val, yend = parameter, colour = bound),
                   linewidth = 5.5, alpha = .85) +
      facet_wrap(~ arm, nrow = 1, scales = "free_x") +
      scale_colour_manual(values = c(lower = "#d6604d", upper = "#4393c3"),
                          labels = c(lower = "Lower", upper = "Upper"), name = "Bound") +
      scale_x_continuous(labels = scales::comma) +
      labs(title = ttl, x = "Symptomatic cases averted", y = NULL) +
      theme_bw(11) +
      theme(text = element_text(size = 14), legend.position = if (legend) "bottom" else "none",
            plot.title = element_text(face = "bold", size = 12),
            strip.text = element_text(face = "bold", size = 10),
            panel.grid.minor = element_blank())
  }
  pA <- tor(ca, blA, "A   Chikungunya", legend = FALSE)
  pB <- tor(mb, blB, "B   Mayaro", legend = TRUE)

  # A spans both columns; B occupies the left column only, aligned beneath it.
  design <- "AA\nB#"
  g <- pA + pB + plot_layout(design = design, heights = c(1, 1.05))
  ggsave("CHIKV_MAYV_owsa.png", g, width = 11, height = 8.4, dpi = 130)
  cat("Saved CHIKV_MAYV_owsa.png (panel A = CHIKV, panel B = MAYV).\n")
}

# ------------------------------------------------------------
# 7. Scenario sweep: vaccine efficacy x coverage.
# For a DISEASE-BLOCKING vaccine, infections are vaccine-invariant, so the sweep needs
# only one SEIR run per coverage level: symptomatic_vaccinated = prop_symp * infections
# * (1 - VE * coverage_frac(t)). The grid is then exact, not interpolated.
# Coverage is of the ELIGIBLE 18-59 group, who are 61.5% of the population, so the
# ceiling on % averted is VE x 0.615 less a small timing discount.
# R0-INVARIANCE: because disease-blocking scales symptomatic cases by (1 - VE*coverage_frac)
# regardless of epidemic size, this grid is essentially the SAME under either fixed-R0
# scenario, provided the outbreak peaks after the campaign completes (it does: peak ~wk 44
# vs campaign done ~wk 29). So one heatmap serves both scenarios -- it is a property of the
# vaccine and the age structure, not of transmission intensity.
# ------------------------------------------------------------
VE_GRID  <- seq(0.1, 1.0, by = 0.1)
COV_GRID <- seq(0.1, 1.0, by = 0.1)

sweep_season <- env_of(BASE$env)
sweep_Rimm   <- rep(BASE$imm, A)
sweep_sus    <- N * (1 - BASE$imm)
sweep_I0     <- E$I0_total * sweep_sus / sum(sweep_sus)
sweep_beta   <- BASE$R0 * GAMMA * sweep_season
sweep_sim <- function(cov) seirv_vaccinated_MAYV(T_weeks, A, N, sweep_Rimm, sweep_I0,
  sweep_beta, SIGMA, GAMMA, BASE$rho, target_age, cov, BASE$deliv,
  start_pre + BASE$delay, VE_inf = 0, VE_block = 0, immun_delay = BASE$immun,
  prop_symp = PSYMP, E0 = E$E0, seed_week = E$seed_week)

sw_base  <- sweep_sim(0)
symp_tot <- sum(colSums(PSYMP * sw_base$new_infections))
cf_list  <- lapply(COV_GRID, function(cv) sweep_sim(cv)$coverage_frac)

grid_df <- do.call(rbind, lapply(seq_along(COV_GRID), function(j) {
  cf <- cf_list[[j]]
  do.call(rbind, lapply(VE_GRID, function(ve) {
    symp_v <- sum(colSums(PSYMP * sw_base$new_infections * (1 - ve * cf)))
    data.frame(coverage = COV_GRID[j], ve = ve,
               averted = symp_tot - symp_v,
               pct = 100 * (symp_tot - symp_v) / symp_tot) })) }))

p_hm <- ggplot(grid_df, aes(factor(100*coverage), factor(100*ve), fill = pct)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.1f", pct)), size = 4.0, colour = "grey15") +
  scale_fill_distiller(palette = "Spectral", name = "% of total\nsymptomatic\ncases averted") +
  labs(x = "Vaccine coverage (% of eligible population aged 18-59)",
       y = "Vaccine efficacy (%)",
       title = sprintf("MAYV: %% of symptomatic cases averted, by VE x coverage (fixed peak R0 = %.2f)", BASE$R0),
       subtitle = paste("Disease-blocking only, so the grid is EXACT, not interpolated. Near-invariant to R0:",
                        "the reduction is ~VE x coverage x (share of cases in the eligible 18-59 group).")) +
  theme_bw(11) + theme(panel.grid = element_blank(),
                       plot.title = element_text(face = "bold", size = 11),
                       plot.subtitle = element_text(size = 8.5))
ggsave("MAYV_ca_ve_coverage_heatmap.png", p_hm, width = 8.5, height = 6, dpi = 130)
print(p_hm)

write_xlsx(list(ve_coverage_grid = grid_df), "MAYV_ca_ve_coverage_sweep.xlsx")
cat(sprintf("\nVE x coverage sweep: %d cells | base case (VE %.0f%%, cov %.0f%%) = %.2f%% averted | max (VE 100%%, cov 100%%) = %.2f%%\n",
            nrow(grid_df), 100*BASE$ve, 100*BASE$cov,
            grid_df$pct[grid_df$ve == 0.5 & grid_df$coverage == 0.3],
            max(grid_df$pct)))
cat("Saved MAYV_ca_ve_coverage_heatmap.png and MAYV_ca_ve_coverage_sweep.xlsx\n")
