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
#     continuously at 2026-W10.
#
# WHY. The CHIKV beta_t's post-outbreak tail is a spline artefact -- flat (~0.9) and
# unconstrained (no CHIKV cases there), with no real dry season. Against MAYV's high
# R0 that flat tail keeps R_eff ~ 1, so the outbreak never resolves inside the window.
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
# Output: caldas_hybrid_season.rds (68-wk, mean-1) + MAYV_ca_hybrid_envelope.png.
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
df <- rbind(
  data.frame(ord = 24:75, val = beta,   env = "CHIKV beta (as-is)"),
  data.frame(ord = 40:91, val = rain,   env = "Rainfall (as-is)"),
  data.frame(ord = ords,  val = hybrid, env = "HYBRID (beta rise + dry tail)"))
p <- ggplot(df, aes(ord, val, colour = env)) + geom_line(linewidth = 1) +
  geom_vline(xintercept = TRANS, linetype = "dashed", colour = "grey50") +
  annotate("text", x = TRANS, y = Inf, label = "join\n2026-W10", vjust = 1.2, size = 3) +
  scale_x_continuous(breaks = c(24, 40, 50, 54, 63, 75, 91),
                     labels = wk_lab(c(24, 40, 50, 54, 63, 75, 91))) +
  labs(x = "epi-week", y = "transmission envelope (mean-1)", colour = NULL,
       title = "MAYV seasonal envelope: hybrid = CHIKV beta (rise/peak) + climatological dry-season tail") +
  theme_bw(11) + theme(legend.position = "bottom")
ggsave("MAYV_ca_hybrid_envelope.png", p, width = 9, height = 4.6, dpi = 120)
cat("Wrote MAYV_ca_hybrid_envelope.png\n")
