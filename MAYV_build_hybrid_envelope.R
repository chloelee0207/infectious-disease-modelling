# ============================================================
# MAYV_build_hybrid_envelope.R
# ------------------------------------------------------------
# Builds caldas_hybrid_season.rds -- the MAYV transmission seasonality envelope used
# by MAYV_ca_lhs.R. It splices two sources, using each only where it is informative:
#
#   * RISE + PEAK (2025-W24 .. 2026-W10): the FITTED Caldas Novas CHIKV beta_t
#     (caldas_beta_season.rds) -- an empirical transmission signal for the shared
#     Aedes vector, well constrained by CHIKV case data.
#   * DRY-SEASON TAIL (2026-W10 .. 2026-W22): the CHIRPS climatological rainfall
#     shape (caldas_rain_season.rds, from rainfall_season.R), scaled to join
#     continuously at 2026-W10. Rainfall is exogenous to the epidemic, so unlike the
#     fitted beta_t it gives a seasonal shape uncontaminated by susceptible depletion.
#
# WHY. NOT for lack of data: 39.5% of the observed CHIKV cases (3,237) fall at or after
# the 2026-W10 join, and only 2 of the 52 weeks have none. The reason is CONFOUNDING with
# susceptible depletion. Over the tail the susceptible pool drops to 45% of its starting
# value by 2026-W10 and 21% by 2026-W22, which on its own explains the 98% fall in
# reported cases -- so the likelihood applies no pressure for beta_t to decline, and it
# stays flat (0.638 -> 0.630, -1.2%). The tail therefore carries depletion, not
# seasonality, and cannot be read as a seasonal signal.
#
# Transplanting it to MAYV would be wrong twice over: the MAYV population is ~93%
# susceptible, so there is no depletion to damp transmission, and a flat beta against
# MAYV's higher R0 holds R_eff ~ 1 so the outbreak never resolves inside the window.
# The rainfall dry season (dropping to ~0.12 by 2026-W22) crashes R_eff after the peak,
# so the outbreak collapses -- as CHIKV's own (lower-R0) outbreak does. Result: a
# 52-week envelope, 2025-W24 -> 2026-W22 (CHIKV-aligned window), inside which the
# outbreak peaks and ~98% resolves (only the deep dry-season trough past W22 is cut off).
#
# Ordinal-week bookkeeping (handles the 2025-W53 leap week):
#   2025 W24..W53 = 24..53 ; 2026 W01..W38 = 54..91
#   caldas_beta_season.rds: idx i -> ordinal 23+i  (W24=1 .. W22=52)  covers 24..75
#   caldas_rain_season.rds: idx i -> ordinal 39+i  (W40=1 .. W38=52)  covers 40..91
#
# Output: caldas_hybrid_season.rds (52-wk, mean-1) + MAYV_ca_hybrid_envelope.png
#         and a 600-dpi copy, MAYV_ca_hybrid_envelope_600dpi.png, for the manuscript.
# Re-run whenever caldas_beta_season.rds or caldas_rain_season.rds changes.
# ============================================================
suppressMessages(library(ggplot2))

beta <- readRDS("caldas_beta_season.rds")   # 2025-W24 .. 2026-W22
rain <- readRDS("caldas_rain_season.rds")   # 2025-W40 .. 2026-W38
b_at <- function(o) beta[o - 23]            # valid ordinal 24..75
r_at <- function(o) rain[o - 39]            # valid ordinal 40..91

TRANS <- 63                                 # join at 2026-W10 (beta data-constrained up to here)
ords  <- 24:75                              # hybrid window 2025-W24 -> 2026-W22 (52 weeks, CHIKV-aligned)
scale <- b_at(TRANS) / r_at(TRANS)          # match level across the join -> continuous
hybrid <- sapply(ords, function(o) if (o <= TRANS) b_at(o) else r_at(o) * scale)
hybrid <- hybrid / mean(hybrid)             # store mean-1 (the engine peak-scales anyway)

stopifnot(length(hybrid) == 52, abs(mean(hybrid) - 1) < 1e-9)
saveRDS(hybrid, "caldas_hybrid_season.rds")
cat(sprintf("Wrote caldas_hybrid_season.rds: %d wks, mean=%.3f, max=%.2f, min-after-peak=%.3f\n",
            length(hybrid), mean(hybrid), max(hybrid), min(hybrid[which.max(hybrid):length(hybrid)])))

# --- verification plot: the three envelopes on a shared ordinal-week axis ---
wk_lab <- function(o) ifelse(o <= 53, o, o - 53)
# Series labels. The colours are PINNED with scale_colour_manual rather than left to the
# default hue palette: that palette assigns colours in factor order, which is alphabetical
# for a character column, and "Average rainfall..." would sort first and steal the red from
# the CHIKV series. Explicit hex values are the same three defaults, so the figure keeps the
# red/green/blue reading it has always had whatever the labels say.
LAB_BETA <- "CHIKV fitted weekly transmission rates \u03b2(t)"
LAB_HYB  <- "Hybrid seasonal envelope"
LAB_RAIN <- "Average rainfall in Caldas Novas"
ENV_COL  <- c("#F8766D", "#00BA38", "#619CFF")            # ggplot hue_pal()(3): red, green, blue
names(ENV_COL) <- c(LAB_BETA, LAB_HYB, LAB_RAIN)

df <- rbind(
  data.frame(ord = 24:75, val = beta,   env = LAB_BETA),
  data.frame(ord = 40:91, val = rain,   env = LAB_RAIN),
  data.frame(ord = ords,  val = hybrid, env = LAB_HYB))
df$env <- factor(df$env, levels = names(ENV_COL))          # legend order, not alphabetical
p <- ggplot(df, aes(ord, val, colour = env)) + geom_line(linewidth = 1) +
  geom_vline(xintercept = TRANS, linetype = "dashed", colour = "grey50") +
  annotate("text", x = TRANS, y = Inf, label = "join\n2026-W10", vjust = 1.2, size = 3) +
  scale_x_continuous(breaks = c(24, 40, 50, 54, 63, 75, 91),
                     labels = wk_lab(c(24, 40, 50, 54, 63, 75, 91))) +
  scale_colour_manual(values = ENV_COL, name = NULL) +
  labs(x = "Week", y = "transmission envelope (mean-1)", colour = NULL) +
  theme_bw(11) + theme(legend.position = "bottom")
ggsave("MAYV_ca_hybrid_envelope.png", p, width = 9, height = 4.6, dpi = 120)
# 600-dpi companion for the manuscript: identical physical size and layout, so the two
# are interchangeable in a document and only the pixel density differs.
ggsave("MAYV_ca_hybrid_envelope_600dpi.png", p, width = 9, height = 4.6, dpi = 600)
cat("Wrote MAYV_ca_hybrid_envelope.png and MAYV_ca_hybrid_envelope_600dpi.png\n")
