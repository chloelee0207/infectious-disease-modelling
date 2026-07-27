# ============================================================
# MAYV_ca_owsa.R -- Caldas Novas MAYV one-way sensitivity analysis.
#
# DETERMINISTIC, and CONDITIONAL ON AN OUTBREAK. MAYV is not fitted to data, so no
# re-fitting is needed: every run is a forward simulation at chosen parameter values.
#
# BASE CASE. A single introduction only self-sustains when peak R0 clears ~2.4, and
# just 19% of prior draws do, so the base case fixes R0 at the OUTBREAK-CONDITIONAL
# median. R0 is treated as a SCENARIO dimension and is not varied in the tornado: it
# sets whether an outbreak happens and how large, and within even its conditional range
# it moves averted burden ~100-fold, which would swamp every other parameter. Its
# effect is reported separately in the notes sheet.
#
# The tornado therefore answers: GIVEN a MAYV outbreak of the central size, what drives
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
ob <- M$outbreak                                   # draws that took off
R0_BASE <- unname(median(E$R0[ob]))                # outbreak-conditional median

BASE <- list(R0 = R0_BASE, imm = 0.0735, rho = 0.25, ve = 0.50,
             cov = 0.30, deliv = 0.10, delay = 2, immun = 2, env = "hybrid")

# TWO inputs are deliberately NOT tornado rows, because neither is a parameter with an
# ordered uncertainty range; both are structural choices reported in the
# structural_sensitivity sheet instead:
#   R0       -- defines the SCENARIO (whether an outbreak occurs and how large). Even
#               within its outbreak-conditional range it moves averted burden ~100-fold.
#   envelope -- a model-structure choice between the pure CHIKV beta shape, the pure
#               rainfall shape and the hybrid. "Lower" and "upper" are meaningless for
#               a categorical choice, and the rainfall arm additionally peaks before the
#               campaign completes, so a single bar would conflate size with timing.
# The reporting rate is also excluded, but for a different reason: MAYV is not fitted,
# so rho never enters the transmission model. The epidemic size comes from R0, prior
# immunity and the seed; rho is applied afterwards only to convert TRUE cases into
# REPORTED ones. It therefore has exactly zero effect on true cases averted.
R0_COND <- unname(quantile(E$R0[ob], c(.025, .975)))
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
cat(sprintf("MAYV OWSA: base peak R0 = %.2f (outbreak-conditional median of %d draws)\n",
            R0_BASE, length(ob)))
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
  do.call(rbind, lapply(c(R0_COND[1], R0_BASE, R0_COND[2]), function(x) {
    r <- run_scenario(modifyList(BASE, list(R0 = x)))
    data.frame(input = "Peak R0 (scenario)", setting = sprintf("%.2f", x),
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
# 5. Export
# ------------------------------------------------------------
notes <- data.frame(item = c("Analysis", "Base case R0", "R0 as a scenario", "Vaccine",
                             "Seasonal envelope", "Fixed (not varied)", "Window"),
  detail = c(
  "Deterministic one-way sensitivity; MAYV is not fitted, so every run is a forward simulation.",
  sprintf("Peak R0 %.2f = median among the %d/%d prior draws that took off (attack > %.1f%%). At the unconditional prior median (2.03) a single introduction does not self-sustain.",
          R0_BASE, length(ob), M$N_DRAWS, M$OUTBREAK_ATTACK_THRESH),
  sprintf("Results are conditional on an outbreak occurring (%.1f%% of prior draws). Peak R0 and the seasonal envelope are reported in the structural_sensitivity sheet, not as tornado rows: neither is a parameter with an ordered uncertainty range, and each moves the result more than any vaccine parameter. Below peak R0 ~%.2f a single introduction does not self-sustain at all.",
          100*M$p_outbreak, min(E$R0[ob])),
  "Disease-blocking only (VE_inf = 0): infections are identical across arms; only symptomatic cases and downstream outcomes move.",
  "Inherited from the CHIKV fit (hybrid = CHIKV beta rise + climatological dry-season tail). Reported as a structural sensitivity. NB the pure-rainfall envelope peaks at week 16, before the campaign completes (~week 29), so its low averted burden reflects timing as well as epidemic size.",
  "gamma, sigma, prop_symp (as in the CHIKV OWSA). Also the reporting rate: MAYV is not fitted, so rho does not enter the transmission model at all -- it converts true cases to reported cases after the simulation and has exactly zero effect on true burden averted. (In the CHIKV OWSA rho ranks third, because there it works backwards from the observed 8,204 cases to infer the size of the true epidemic.)",
  sprintf("52 weeks, 2025-W24 -> 2026-W22 (indices %d-%d).", min(EVAL_WIN), max(EVAL_WIN))),
  stringsAsFactors = FALSE)
write_xlsx(list(notes = notes,
                base_case = data.frame(symptomatic_base = base_run$symptomatic_base,
                                       averted = base_run$averted, attack_pct = base_run$attack,
                                       pct_reduction = 100*base_run$averted/base_run$symptomatic_base),
                owsa = owsa, structural_sensitivity = struct), "MAYV_ca_owsa.xlsx")
saveRDS(list(owsa = owsa, struct = struct, base = base_run, BASE = BASE, BOUNDS = BOUNDS,
             R0_BASE = R0_BASE, p_outbreak = M$p_outbreak), "MAYV_ca_owsa.rds")

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
      theme(legend.position = if (legend) "bottom" else "none",
            plot.title = element_text(face = "bold", size = 11),
            strip.text = element_text(face = "bold", size = 9),
            panel.grid.minor = element_blank())
  }
  pA <- tor(ca, blA, "A   Chikungunya", legend = FALSE)
  pB <- tor(mb, blB, "B   Mayaro (conditional on outbreak)", legend = TRUE)

  # A spans both columns; B occupies the left column only, aligned beneath it.
  design <- "AA\nB#"
  g <- pA + pB + plot_layout(design = design, heights = c(1, 1.05))
  ggsave("CHIKV_MAYV_owsa.png", g, width = 11, height = 8.4, dpi = 130)
  cat("Saved CHIKV_MAYV_owsa.png (panel A = CHIKV, panel B = MAYV).\n")
}
