# ============================================================
# MAYV_outputs.R -- Caldas Novas MAYV presentation & export layer
# ------------------------------------------------------------
# The MODEL / SCENARIO / Monte-Carlo work lives in MAYV_ca_vacc.R (the "engine").
# This script is PURELY presentation: it takes the engine's in-memory result
# objects and (1) draws the figures and (2) writes the Excel workbooks.
#
# Run order in one R session:
#     source("ca_common.R")          # shared helpers (fmtq, ...)
#     source("MAYV_ca_pre_vacc.R")   # forward MAYV model -> season, N, rho, ... (slow-ish)
#     source("MAYV_ca_vacc.R")       # scenarios + Monte Carlo (the engine)
#     source("MAYV_outputs.R")       # <- this file: figures + Excel
#
# Vaccine design (from the engine): DISEASE-BLOCKING ONLY (VE_inf = 0), central
# VE_block = 50% (25/75% sensitivity), two fixed outbreak sizes (R0 = 1.2 / 3.51).
# ============================================================
library(ggplot2)
library(ggh4x)
library(writexl)
if (file.exists("ca_common.R")) source("ca_common.R")   # fmtq(), ...

# ------------------------------------------------------------
# 0. Guard: the engine must have been run first
# ------------------------------------------------------------
.needed <- c("perm", "size_R0", "eff_lvl", "eff_central_lab", "perm_mc", "av_mc",
             "vac_rows", "tbl_total", "tbl_averted", "sens", "outcomes",
             "outcome_labs", "dgt", "rho", "total_coverage", "mayv_vacc_efficacy",
             "T_weeks", "start_pre", "seed_week", "R0_central")
.missing <- .needed[!vapply(.needed, exists, logical(1))]
if (length(.missing)) {
  stop("MAYV_outputs.R needs objects from the engine. Run MAYV_ca_vacc.R first.\n",
       "  Missing: ", paste(.missing, collapse = ", "))
}

# ============================================================
# 1. Figures (grid: outbreak-size rows x outcome columns, independent y-axes)
# ============================================================
# Each panel has one bar per x category (no dodging), so bars and error bars align.
mc_long <- function(mats, rows) do.call(rbind, lapply(seq_along(mats), function(k) {
  m <- mats[[k]]; j <- rows[k]
  data.frame(size = perm$size[j], x = perm$eff[j], outcome = outcomes,
             med = apply(m, 2, median,   na.rm = TRUE),
             lo  = apply(m, 2, quantile, .025, na.rm = TRUE),
             hi  = apply(m, 2, quantile, .975, na.rm = TRUE), row.names = NULL)
}))
fmt_axes <- function(d, x_levels) {
  d$size    <- factor(d$size, levels = names(size_R0))
  d$x       <- factor(d$x,    levels = x_levels)
  d$outcome <- factor(d$outcome, levels = outcomes, labels = outcome_labs)
  d
}
grid_plot <- function(d, fill_vals, ytitle, title) {
  ggplot(d, aes(x, med, fill = x)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.4) +
    ggh4x::facet_grid2(size ~ outcome, scales = "free_y", independent = "y") +
    scale_fill_manual(values = fill_vals, name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = ytitle, title = title) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          legend.position = "none", panel.grid.minor = element_blank())
}

# (a) TOTAL burden: no vaccine + central 50% disease-blocking
sel_tot <- which(perm$eff %in% c("No vaccine", eff_central_lab))
p_total <- grid_plot(
  fmt_axes(mc_long(perm_mc[sel_tot], sel_tot), c("No vaccine", eff_central_lab)),
  setNames(c("grey60", "#2166ac"), c("No vaccine", eff_central_lab)),
  "Total burden, median + 95% UI",
  "Hypothetical MAYV cases - total burden: no vaccine vs 50% disease-blocking")
print(p_total)
ggsave("ca_mayv_vacc_burden.png", p_total, width = 10, height = 5.5, dpi = 120)

# (b) AVERTED: central 50% disease-blocking only
av_is_50 <- which(perm$eff[vac_rows] == eff_central_lab)
p_av <- grid_plot(
  fmt_axes(mc_long(av_mc[av_is_50], vac_rows[av_is_50]), eff_central_lab),
  setNames("#08519c", eff_central_lab),
  "Burden averted vs no vaccine, median + 95% UI",
  "Hypothetical MAYV cases - burden averted by 50% disease-blocking vaccine")
print(p_av)
ggsave("ca_mayv_vacc_averted.png", p_av, width = 10, height = 5.5, dpi = 120)

# (c) AVERTED efficacy sensitivity: 25% / 50% / 75% disease-blocking
p_av_sens <- grid_plot(
  fmt_axes(mc_long(av_mc, vac_rows), names(eff_lvl)),
  setNames(c("#9ecae1", "#4292c6", "#08519c"), names(eff_lvl)),
  "Burden averted vs no vaccine, median + 95% UI",
  "Hypothetical MAYV cases - burden averted by disease-blocking efficacy (25/50/75%)")
print(p_av_sens)
ggsave("ca_mayv_vacc_averted_sensitivity.png", p_av_sens, width = 10, height = 5.5, dpi = 120)

# ============================================================
# 2. Excel: burden (total + averted) and R0 sensitivity (moved from the engine)
# ============================================================
write_xlsx(list(total_burden = tbl_total, averted = tbl_averted),
           "caldas_mayv_vacc_burden.xlsx")
write_xlsx(list(R0_sensitivity = sens), "caldas_mayv_vacc_R0_sensitivity.xlsx")

# ============================================================
# 3. Excel: baseline / vaccinated / averted, each as TRUE and REPORTED (95% UI)
# ============================================================
# TRUE  = the modelled burden (infections, symptomatic, hospitalisations, deaths).
# REPORTED = rho x TRUE, i.e. the surveillance-visible fraction at the fixed MAYV
# reporting rate rho = 0.10, applied uniformly to every outcome (a simplification --
# hospitalisations/deaths are in reality better ascertained than cases). rho is
# fixed, so TRUE and REPORTED carry the SAME parameter-propagation 95% UI, differing
# only by the constant factor rho.
tr_tbl <- function(rows, mats) do.call(rbind, lapply(seq_along(rows), function(k) {
  m <- mats[[k]]; j <- rows[k]
  data.frame(
    `Outbreak size` = perm$size[j],
    Vaccine         = perm$eff[j],
    Outcome         = outcome_labs,
    true            = vapply(outcomes, function(o) fmtq(m[, o],       dgt[[o]]), character(1)),
    reported        = vapply(outcomes, function(o) fmtq(rho * m[, o], dgt[[o]]), character(1)),
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
}))

base_rows  <- which(perm$eff == "No vaccine")              # perm rows (Current, Future)
vac50_rows <- which(perm$eff == eff_central_lab)           # perm rows (Current, Future)
av50_k     <- which(perm$eff[vac_rows] == eff_central_lab) # indices into av_mc / vac_rows

baseline_true_reported   <- tr_tbl(base_rows,        perm_mc[base_rows])
vaccinated_true_reported <- tr_tbl(vac50_rows,       perm_mc[vac50_rows])
averted_true_reported    <- tr_tbl(vac_rows[av50_k], av_mc[av50_k])
# Make the averted value columns self-describing.
colnames(averted_true_reported)[4:5] <- c("averted (true)", "averted (reported)")

# Live run parameters / assumptions.
future_R0 <- unname(size_R0[["Future outbreak (R0 = 3.51)"]])
notes <- data.frame(
  parameter = c("Outbreak window",
                "Outbreak size R0 (Current / Future)",
                "Vaccine mechanism",
                "Infection-blocking efficacy (VE_inf)",
                "Central disease-blocking efficacy (VE_block)",
                "Efficacy sensitivity levels",
                "Target coverage of adults 18-59",
                "Reporting rate (rho)",
                "REPORTED definition",
                "Severity (hosp, CFR) source",
                "95% UI source",
                "Pre-outbreak campaign start (week_index)",
                "Wet-season seed (week_index)"),
  value = c(sprintf("2025-W23 -> 2026-W22 (%d weeks)", T_weeks),
            sprintf("%.2f / %.2f", R0_central, future_R0),
            "Disease-blocking only",
            "0% (VE_inf = 0 in every scenario -> infections never averted)",
            sprintf("%.0f%%", 100 * mayv_vacc_efficacy),
            paste0(paste(100 * eff_lvl, collapse = "/"), "%"),
            sprintf("%.0f%%", 100 * total_coverage),
            sprintf("%.2f", rho),
            "REPORTED = rho x TRUE, applied uniformly to all outcomes (simplification)",
            "Borrowed CHIKV disease-progression (Hyolim Table S4) - CHIKV-equivalent UPPER bound",
            "Parameter propagation (CHIKV seasonal-shape posterior + severity Betas); R0 fixed",
            as.character(start_pre),
            as.character(seed_week)),
  stringsAsFactors = FALSE)

sheets <- list(
  notes                    = notes,
  baseline_true_reported   = baseline_true_reported,     # No vaccine, both R0
  vaccinated_true_reported = vaccinated_true_reported,   # 50% disease-blocking, both R0
  averted_true_reported    = averted_true_reported       # 50% disease-blocking averted, both R0
)
write_xlsx(sheets, "caldas_mayv_vacc_outputs.xlsx")

cat("\n=== Baseline (no vaccine) TRUE & REPORTED (median, 95% UI) ===\n")
print(baseline_true_reported, row.names = FALSE)
cat("\n=== Vaccinated 50% disease-blocking TRUE & REPORTED (median, 95% UI) ===\n")
print(vaccinated_true_reported, row.names = FALSE)
cat("\n=== Averted (50% disease-blocking) TRUE & REPORTED (median, 95% UI) ===\n")
print(averted_true_reported, row.names = FALSE)

cat("\nWrote caldas_mayv_vacc_outputs.xlsx (sheets:",
    paste(names(sheets), collapse = ", "), ")\n")
cat("Wrote caldas_mayv_vacc_burden.xlsx, caldas_mayv_vacc_R0_sensitivity.xlsx\n")
cat("Saved figures: ca_mayv_vacc_burden.png, ca_mayv_vacc_averted.png,",
    "ca_mayv_vacc_averted_sensitivity.png\n")
