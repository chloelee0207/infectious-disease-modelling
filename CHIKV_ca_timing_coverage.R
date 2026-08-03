# ============================================================
# CHIKV_ca_timing_coverage.R
# ------------------------------------------------------------
# Scenario grid: WHEN the campaign starts x HOW MUCH coverage it reaches, against the
# share of symptomatic chikungunya cases averted. The point is the timing axis. Ixchiq
# was intended as a prophylactic, pre-outbreak intervention; in Caldas Novas it was
# announced on 18 April 2026 (2026-W16), well after the modelled peak. This figure shows
# what that delay costs, at every coverage level.
#
# Deterministic at the central parameter set (the same one CHIKV_ca_owsa.R calls BASE),
# not a Monte Carlo: 52 start weeks x 10 coverage levels x 2 arms is a scenario surface,
# and re-propagating every cell would obscure the pattern without changing its shape.
# Uncertainty for the two REPORTED timings is already carried by the main results.
#
# The epidemic fit is done ONCE -- FOI and rho are held at their central values, so the
# baseline epidemic is identical in every cell and only the campaign moves.
#
# Output: CHIKV_ca_timing_coverage.png / .pdf and CHIKV_ca_timing_coverage.xlsx
# ============================================================
DEFS_ONLY <- TRUE
suppressMessages({library(ggplot2); library(dplyr); library(writexl)})
source("CHIKV_ca_lhs.R")          # refit(), data, fixed choices, exposure_age, E0, gen_start
source("ca_common.R")             # seirv_vaccinated

GAMMA <- 0.54; SIGMA <- 1/0.60; PSYMP <- prop_symp
FOI <- 0.008; RHO <- 0.25; VE <- 263/266; DELIV <- 0.10; IMMUN <- 2
target_age <- rep(0, A); target_age[c(4,5,6,7,8)] <- 1        # eligible 18-59
EVAL_WIN <- seq_len(T_weeks)

# the two timings that matter, as week indices into the 52-week window
idx_of <- function(y, w) caldas_obs$week_index[caldas_obs$Year == y & caldas_obs$week == w]
WK_PRE <- idx_of(2025, 40)    # intended: pre-outbreak prophylactic campaign
WK_ACT <- idx_of(2026, 16)    # actual: announced 18 April 2026, after the peak

# ---- one fit, reused by every cell -------------------------------------------------
fit  <- refit(FOI, GAMMA, SIGMA, RHO, PSYMP, gen_start)
if (is.null(fit)) stop("re-fit failed at the central FOI/rho")
Rimm  <- 1 - exp(-FOI * exposure_age)
sfrac <- (N*(1-Rimm))/sum(N*(1-Rimm))
I0i   <- round(((week_1_cases/RHO/PSYMP)/GAMMA) * sfrac)

# Total symptomatic cases. The engine's age reweighting rescales to preserve the total,
# so for a TOTAL (rather than an age split) the raw sum is the same number.
symp_total <- function(out) sum(out$new_symptomatic[, EVAL_WIN])
sim <- function(cov, start, vi, vb)
  seirv_vaccinated(T_weeks, A, N, Rimm, I0i, E0, fit$beta, SIGMA, GAMMA, RHO,
                   target_age, cov, DELIV, start, vi, vb, IMMUN, prop_symp = PSYMP)
BASE_SYMP <- symp_total(sim(0, WK_PRE, 0, 0))     # campaign week is irrelevant at cov = 0
cat(sprintf("Baseline symptomatic (central inputs): %.0f\n", BASE_SYMP))

# ---- the grid ----------------------------------------------------------------------
WEEKS <- seq_len(T_weeks)
COVS  <- seq(0.10, 1.00, by = 0.10)
ARMS  <- c("Disease-blocking", "Disease + infection blocking")
cat(sprintf("Sweeping %d start weeks x %d coverage levels x %d arms...\n",
            length(WEEKS), length(COVS), length(ARMS)))
grid <- do.call(rbind, lapply(WEEKS, function(w) do.call(rbind, lapply(COVS, function(cv) {
  db <- symp_total(sim(cv, w, 0,  VE))
  di <- symp_total(sim(cv, w, VE, VE))
  data.frame(week = w, coverage = cv,
             arm = factor(ARMS, levels = ARMS),
             pct_averted = 100 * (BASE_SYMP - c(db, di)) / BASE_SYMP,
             stringsAsFactors = FALSE) }))))
grid$epi_week <- sprintf("%d-W%02d", caldas_obs$Year[grid$week], caldas_obs$week[grid$week])

# ---- the fold-loss between the intended and the actual timing ----------------------
fold <- grid |>
  filter(week %in% c(WK_PRE, WK_ACT)) |>
  mutate(timing = ifelse(week == WK_PRE, "intended_pre_outbreak", "actual_2026W16")) |>
  select(arm, coverage, timing, pct_averted) |>
  tidyr::pivot_wider(names_from = timing, values_from = pct_averted) |>
  mutate(fold_loss = intended_pre_outbreak / actual_2026W16) |>
  as.data.frame()
cat("\nIntended (2025-W40) vs actual (2026-W16) announcement:\n")
print(fold |> mutate(across(where(is.numeric), \(x) round(x, 2))), row.names = FALSE)

# ---- figure ------------------------------------------------------------------------
lab_at <- c(1, 10, 20, 30, 40, 52)
p <- ggplot(grid, aes(100*coverage, week, fill = pct_averted)) +
  geom_raster(interpolate = TRUE) +
  geom_hline(yintercept = WK_PRE, colour = "grey15", linetype = "22", linewidth = .5) +
  geom_hline(yintercept = WK_ACT, colour = "grey15", linetype = "solid", linewidth = .5) +
  annotate("text", x = 12, y = WK_PRE + 2.2, hjust = 0, size = 3.1, colour = "grey15",
           label = "Intended: pre-outbreak, 2025-W40") +
  annotate("text", x = 12, y = WK_ACT + 2.2, hjust = 0, size = 3.1, colour = "grey15",
           label = "Actual: announced 18 Apr 2026, 2026-W16") +
  facet_wrap(~ arm, nrow = 1) +
  scale_fill_gradientn(colours = c("#3b7fb6","#7fcdbb","#d9f0a3","#fee391","#fc8d59","#d73027"),
                       name = "% of symptomatic\ncases averted") +
  scale_x_continuous(breaks = seq(10, 100, 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = lab_at, labels = sprintf("%d\n(%s)", lab_at,
                       sprintf("%d-W%02d", caldas_obs$Year[lab_at], caldas_obs$week[lab_at])),
                     expand = c(0, 0)) +
  labs(x = "Vaccine coverage (% of eligible population aged 18-59)",
       y = "Start of vaccination campaign (week index)") +
  theme_bw(12) +
  theme(strip.text = element_text(face = "bold", size = 11),
        panel.grid = element_blank(), legend.position = "right")
ggsave("CHIKV_ca_timing_coverage.png", p, width = 11, height = 5.6, dpi = 150)
ggsave("CHIKV_ca_timing_coverage.pdf", p, width = 11, height = 5.6)

notes <- data.frame(item = c("Design", "Why deterministic", "Timings", "Fit", "Coverage",
                             "Outcome", "Reading the surface"),
  detail = c(
  sprintf("%d start weeks (the whole 52-week window) x %d coverage levels (10-100%%) x 2 vaccine arms.", length(WEEKS), length(COVS)),
  "A scenario surface at the central parameter set, not a Monte Carlo. Uncertainty for the two reported timings is carried by the main results; re-propagating 1,040 cells would obscure the pattern without changing its shape.",
  sprintf("Intended prophylactic campaign 2025-W40 (week %d); actual announcement 18 April 2026 = 2026-W16 (week %d), after the modelled peak.", WK_PRE, WK_ACT),
  sprintf("FOI = %.3f and rho = %.2f held at their central values, so the baseline epidemic (%.0f symptomatic cases) is identical in every cell and only the campaign moves.", FOI, RHO, BASE_SYMP),
  "Coverage is of the ELIGIBLE 18-59 group, who are 61.5% of the population -- not of the whole population.",
  "Share of TOTAL symptomatic cases averted over the 52-week window.",
  "Impact is governed by timing far more than by coverage: moving down a column (earlier start) gains more than moving right along a row (more doses)."),
  stringsAsFactors = FALSE)
write_xlsx(list(grid = grid, intended_vs_actual = fold, notes = notes),
           "CHIKV_ca_timing_coverage.xlsx")
cat("\nSaved CHIKV_ca_timing_coverage.png / .pdf / .xlsx\n")
