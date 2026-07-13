# ============================================================
# CHIKV_ca_nnv.R -- Caldas Novas CHIKV Number Needed to Vaccinate: PRESENTATION ONLY.
# ------------------------------------------------------------
# Reads CHIKV_ca_engine_results.rds (the unified single-LHS engine) and draws the NNV
# figure + writes the workbook. It does NOT run any SEIR/MC of its own -- NNV is one
# of the per-draw outcomes the engine already produced:
#     nnv_i(outcome) = doses_i / (baseline_i - scenario_i)
# at that draw's SAMPLED coverage (~Beta(30%,20-40)), so the NNV 95% UI is the same
# uncertainty propagation as the burden/DALY numbers. NNV < 1 is possible (one dose
# averts >1 case via herd protection).
#
# Run order:  source("CHIKV_ca_engine.R")   # produces the RDS (once)
#             source("CHIKV_ca_nnv.R")       # this file
# ============================================================
library(dplyr); library(ggplot2); library(writexl)
if (!exists("fmtq")) source("ca_common.R")

if (!file.exists("CHIKV_ca_engine_results.rds"))
  stop("CHIKV_ca_engine_results.rds not found -- run CHIKV_ca_engine.R first.")
G <- readRDS("CHIKV_ca_engine_results.rds")
nnv <- G$nnv; vac_names <- G$vac_names; timings <- G$timings; arm_names <- G$arm_names
NNV_OUT <- G$NNV_OUT; N_DRAWS <- G$N_DRAWS; cov_d <- G$cov_d

# NNV figure/table are reported for the PRE-OUTBREAK rollout (per request).
FOCUS_TIMING <- "pre-outbreak"
focus_names  <- grep(FOCUS_TIMING, vac_names, value = TRUE)

out_labs <- c(symptomatic = "Symptomatic case", hospitalisations = "Hospitalisation",
              deaths = "Death", daly = "DALY")

# ------------------------------------------------------------
# 1. Assemble the plotting frame from the per-draw NNV (median + 95% UI)
# ------------------------------------------------------------
plt <- do.call(rbind, lapply(focus_names, function(nm) {
  m <- nnv[[nm]]
  do.call(rbind, lapply(NNV_OUT, function(o) {
    q <- quantile(m[, o], c(.5, .025, .975), na.rm = TRUE)
    data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm), outcome = o,
               med = q[1], lo = q[2], hi = q[3], row.names = NULL)
  }))
}))
plt$outcome <- factor(out_labs[plt$outcome], levels = unname(out_labs))
plt$arm     <- factor(plt$arm, levels = arm_names)

# ------------------------------------------------------------
# 2. Figure: NNV to avert each outcome, pre-outbreak, by protection arm
# ------------------------------------------------------------
arm_cols <- c("Disease-blocking" = "#4393c3", "Disease + infection blocking" = "#d6604d")
p_nnv <- ggplot(plt, aes(arm, med, fill = arm)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45", linewidth = .3) +
  geom_col(width = .65, colour = "grey30", linewidth = .2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .2, linewidth = .35) +
  facet_wrap(~ outcome, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = arm_cols, name = NULL) +
  scale_y_log10(breaks = scales::breaks_log(n = 6),
                labels = scales::label_number(drop0trailing = TRUE)) +
  labs(x = NULL, y = "NNV to avert one outcome",
       title = "Number Needed to Vaccinate (NNV)") +
       # subtitle = sprintf("Pre-outbreak rollout, 52-wk window; sampled coverage (median %.0f%%); median + 95%% UI over %d draws",
       #                    100*median(cov_d), N_DRAWS),
       # caption = "Dashed line = NNV 1 (1 dose averts 1 case); NNV < 1 means each dose averts >1 case via herd protection") +
  theme_bw(11) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.grid.minor = element_blank())
print(p_nnv); ggsave("CHIKV_ca_nnv.png", p_nnv, width = 9, height = 4.6, dpi = 120)

# ------------------------------------------------------------
# 3. Workbook: NNV (all timings x arms) median (95% UI), from the same per-draw object
# ------------------------------------------------------------
nnv_tbl <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- nnv[[nm]]
  vals <- setNames(as.list(vapply(NNV_OUT, function(o) fmtq(m[, o], 1), character(1))),
                   unname(out_labs[NNV_OUT]))
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm), vals,
             row.names = NULL, check.names = FALSE)
}))

notes <- data.frame(
  parameter = c("Location / pathogen", "Source", "NNV definition", "Numerator",
                "Denominator", "Coverage", "Evaluation window", "Uncertainty"),
  value = c("Caldas Novas (Goias), CHIKV",
            "CHIKV_ca_engine_results.rds (unified single-LHS engine)",
            "NNV = doses / burden averted vs no-vaccine baseline, per draw (population-level, incl. herd effects)",
            "target population x coverage (eligible 18-59 pop x that draw's sampled coverage)",
            "baseline burden - scenario burden, per outcome",
            sprintf("sampled ~ Beta(30%%, 20-40); median %.1f%%", 100*median(cov_d)),
            "52 epi-weeks 2025-W40 -> 2026-W38",
            sprintf("%d-draw single LHS; NNV is a per-draw ratio so it shares the burden/DALY propagation", N_DRAWS)),
  stringsAsFactors = FALSE)
write_xlsx(list(notes = notes, nnv = nnv_tbl), "CHIKV_ca_nnv_outputs.xlsx")
cat("Wrote CHIKV_ca_nnv.png and CHIKV_ca_nnv_outputs.xlsx\n")
