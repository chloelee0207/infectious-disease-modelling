# ============================================================
# CHIKV_ca_daly.R -- Caldas Novas CHIKV DALYs: PRESENTATION ONLY.
# ------------------------------------------------------------
# Reads CHIKV_ca_engine_results.rds (the unified single-LHS engine) and draws the
# DALY figures + writes the workbook. It does NOT run any SEIR/MC of its own -- YLD
# (acute + chronic), YLL and DALYs are per-draw outcomes the engine already produced,
# consistent draw-for-draw with the burden totals (the deaths behind YLL ARE the
# burden-table deaths). DALY = YLD + YLL, undiscounted (Hyolim Kang CHIK_VIM).
#
# Run order:  source("CHIKV_ca_engine.R")   # produces the RDS (once)
#             source("CHIKV_ca_daly.R")      # this file
# ============================================================
library(dplyr); library(ggplot2); library(writexl)
if (!exists("fmtq")) source("ca_common.R")

if (!file.exists("CHIKV_ca_engine_results.rds"))
  stop("CHIKV_ca_engine_results.rds not found -- run CHIKV_ca_engine.R first.")
G <- readRDS("CHIKV_ca_engine_results.rds")
per_draw <- G$per_draw; averted <- G$averted
scen_names <- G$scen_names; vac_names <- G$vac_names
timings <- G$timings; arm_names <- G$arm_names; N_DRAWS <- G$N_DRAWS; PHASE_MODE <- G$PHASE_MODE
base <- per_draw[["No vaccine (baseline)"]]

# ------------------------------------------------------------
# 1. Console summary + workbook
# ------------------------------------------------------------
cat("\n=== Baseline DALYs, no vaccine (median, 95% UI) ===\n")
for (m in c("yld","yll","daly","deaths"))
  cat(sprintf("  %-8s %s\n", m, fmtq(base[, m], if (m=="deaths") 1 else 0)))

daly_by_scenario <- do.call(rbind, lapply(scen_names, function(nm) {
  m <- per_draw[[nm]]; b <- nm == "No vaccine (baseline)"
  data.frame(timing = if (b) "No vaccine" else sub(" \\|.*","",nm),
             arm    = if (b) "No vaccine" else sub(".*\\| ","",nm),
             YLD = fmtq(m[,"yld"]), YLL = fmtq(m[,"yll"]), DALY = fmtq(m[,"daly"]),
             deaths = fmtq(m[,"deaths"], 1),
             YLD_acute = fmtq(m[,"yld_acute"]), YLD_chronic = fmtq(m[,"yld_chronic"]),
             row.names = NULL, check.names = FALSE)
}))
daly_averted <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- averted[[nm]]
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
             DALY_averted = fmtq(m[,"daly"]),
             pct_DALY = sprintf("%.1f%%", 100*median(m[,"daly"]/base[,"daly"], na.rm=TRUE)),
             row.names = NULL, check.names = FALSE)
}))
cat("\n=== DALYs averted vs no-vaccine baseline (median, 95% UI) ===\n")
ow <- getOption("width"); options(width=250); print(daly_averted, row.names=FALSE); options(width=ow)

notes <- data.frame(
  parameter = c("Location / pathogen", "Source", "DALY framework", "Discounting",
                "YLD phases", "YLL", "Evaluation window", "Uncertainty"),
  value = c("Caldas Novas (Goias), CHIKV",
            "CHIKV_ca_engine_results.rds (unified single-LHS engine)",
            "DALY = YLD + YLL (Hyolim Kang CHIK_VIM structure)",
            "none (0%, no age-weighting)",
            sprintf("%s: acute (hospitalised->severe DW, non-hosp->mild/mod) + chronic (not recovered by 90d)", PHASE_MODE),
            "deaths (by decadal band) x remaining life-years at age of death",
            "52 epi-weeks 2025-W40 -> 2026-W38",
            sprintf("%d-draw single LHS shared with the burden/NNV outcomes (consistent draw-for-draw)", N_DRAWS)),
  stringsAsFactors = FALSE)
write_xlsx(list(notes = notes, daly_by_scenario = daly_by_scenario, daly_averted = daly_averted),
           "CHIKV_ca_daly_outputs.xlsx")
cat("\nWrote CHIKV_ca_daly_outputs.xlsx\n")

# ------------------------------------------------------------
# 2. Figure (a): baseline DALY composition (median)
# ------------------------------------------------------------
comp <- data.frame(component = c("YLD (acute)","YLD (chronic)","YLL"),
                   value = c(median(base[,"yld_acute"]), median(base[,"yld_chronic"]), median(base[,"yll"])))
comp$component <- factor(comp$component, levels = comp$component)
p_comp <- ggplot(comp, aes("Baseline", value, fill = component)) +
  geom_col(width = .55) +
  scale_fill_manual(values = c("#c6dbef","#2171b5","#d6604d"), name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "DALYs (median)",
       title = "Caldas Novas CHIKV: baseline DALY composition (no vaccine)",
       caption = sprintf("%s YLD; YLL undiscounted; %d-draw unified LHS", PHASE_MODE, N_DRAWS)) +
  theme_bw(11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
print(p_comp); ggsave("CHIKV_ca_daly_composition.png", p_comp, width = 5.5, height = 5, dpi = 120)

# ------------------------------------------------------------
# 3. Figure (b): DALYs averted by timing x arm (median + 95% UI)
# ------------------------------------------------------------
av_long <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- averted[[nm]]
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
             med = median(m[,"daly"], na.rm=TRUE),
             lo = quantile(m[,"daly"],.025,na.rm=TRUE), hi = quantile(m[,"daly"],.975,na.rm=TRUE),
             row.names = NULL)
}))
av_long$timing <- factor(av_long$timing, levels = names(timings))
av_long$arm    <- factor(av_long$arm, levels = arm_names)
p_av <- ggplot(av_long, aes(timing, med, fill = arm)) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.8), width = .25, linewidth = .4) +
  scale_fill_manual(values = c("Disease-blocking"="#9ecae1","Disease + infection blocking"="#08519c"), name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "DALYs averted (vs no vaccine), median + 95% UI",
       title = "Caldas Novas CHIKV: DALYs averted by vaccination timing") +
  theme_bw(11) + theme(plot.title = element_text(face="bold", hjust=.5),
                       axis.text.x = element_text(angle = 20, hjust = 1),
                       legend.position = "bottom", panel.grid.minor = element_blank())
print(p_av); ggsave("CHIKV_ca_daly_averted.png", p_av, width = 8, height = 4.5, dpi = 120)
cat("Saved CHIKV_ca_daly_composition.png and CHIKV_ca_daly_averted.png\n")
