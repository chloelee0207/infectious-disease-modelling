# ============================================================
# MAYV_ca_outputs.R -- Caldas Novas MAYV presentation & export layer.
# ------------------------------------------------------------
# Reads MAYV_ca_engine_results.rds and writes THREE workbooks that mirror the CHIKV
# outputs (CHIKV_ca_outputs.R) tab-for-tab, so the two diseases are directly
# comparable:
#   MAYV_ca_vacc_outputs.xlsx : notes, baseline_true_reported, vaccinated_true_reported,
#                               averted_MC_95UI, averted_per_100k_doses, scenario_totals
#   MAYV_ca_daly_outputs.xlsx : daly_by_scenario, daly_averted
#   MAYV_ca_nnv_outputs.xlsx  : nnv
# (The direct-medical cost workbook, MAYV_ca_costs.xlsx, is produced separately by
#  MAYV_ca_costs.R and is unchanged.)
#
# TWO IMPROVEMENTS over the CHIKV workbook, per request:
#   * averted_MC_95UI `pct_symp` now carries its 95% UI (the CHIKV excel had only the
#     point estimate; the UI previously lived only in the epidemic-curve figure).
#   * daly_averted `pct_DALY` likewise carries its 95% UI.
#
# MAYV-SPECIFIC FRAMING (differs from CHIKV, which is fitted to a real outbreak):
#   * NO TAKE-OFF CONDITIONING. R0 is FIXED per scenario in MAYV_ca_lhs.R (low = 1.20,
#     high = 2.50), so every draw is the same transmission regime and outbreak size is
#     unimodal. Every figure is summarised over ALL draws (G$outbreak is now every row).
#     The 95% UIs carry natural-history / reporting / symptomatic-fraction / prior-immunity
#     / vaccine / severity-DALY uncertainty, NOT the between-setting R0 span.
#   * ONE vaccine scenario: pre-outbreak, DISEASE-BLOCKING ONLY (VE_inf = 0), so
#     infections are never averted (Infections = 0, and Infection NNV = NA).
#   * DEATHS = 0: no confirmed MAYV-attributable death, so CFR = 0 in the engine ->
#     deaths and YLL are zero and DALY = YLD.
#   * Severity/DALY parameters are BORROWED from CHIKV (upper bound), and the seasonal
#     envelope is the hybrid CHIKV-beta + dry-season envelope (2025-W24 -> 2026-W22).
#
# Run order: source("MAYV_ca_lhs.R"); source("MAYV_ca_engine.R"); source(this)
# ============================================================
library(writexl)
if (!exists("fmtq")) source("ca_common.R")   # fmtq(v, d) -> "median (lo - hi)"

if (!file.exists("MAYV_ca_engine_results.rds"))
  stop("MAYV_ca_engine_results.rds not found -- run MAYV_ca_engine.R first.")
G <- readRDS("MAYV_ca_engine_results.rds")

ok        <- G$outbreak                                  # ALL draws (fixed R0 -> no conditioning)
rho_draw  <- G$rho_draw
scen_names <- G$scen_names
vac_names  <- setdiff(scen_names, "No vaccine (baseline)")
base_pd   <- G$per_draw[["No vaccine (baseline)"]]
outcomes  <- c("infections","symptomatic","hospitalisations","deaths")  # burden subset

lab_timing <- function(nm) if (nm == "No vaccine (baseline)") "No vaccine" else sub(" \\|.*","",nm)
lab_arm    <- function(nm) if (nm == "No vaccine (baseline)") "No vaccine" else sub(".*\\| ","",nm)
d_death    <- function(o) if (o == "deaths") 1 else 0

# % reduction (median + 95% UI), matching the figure format "8.1% (4.1 - 13.8%)".
fmtpct <- function(ratio) {
  q <- quantile(ratio, c(.5, .025, .975), na.rm = TRUE)
  sprintf("%.1f%% (%.1f - %.1f%%)", 100*q[1], 100*q[2], 100*q[3])
}

# averted (baseline - scenario), per draw, over ALL draws, burden outcomes.
av <- setNames(lapply(vac_names, function(nm)
  base_pd[ok, outcomes, drop = FALSE] - G$per_draw[[nm]][ok, outcomes, drop = FALSE]), vac_names)
base_symp_ok <- base_pd[ok, "symptomatic"]

# ============================================================
# WORKBOOK 1 -- MAYV_ca_vacc_outputs.xlsx
# ============================================================
# baseline true-vs-reported: reported = rho x true (per-draw rho), over all draws.
base_tbl <- do.call(rbind, lapply(outcomes, function(o) {
  tv <- base_pd[ok, o]; d <- d_death(o)
  data.frame(outcome = o, true = fmtq(tv, d), reported = fmtq(rho_draw[ok]*tv, d),
             row.names = NULL)
}))

# per-scenario TOTALS, true & reported (all scenarios, incl. baseline row).
vtr <- do.call(rbind, lapply(scen_names, function(nm) {
  m <- G$per_draw[[nm]]
  data.frame(timing = lab_timing(nm), arm = lab_arm(nm), outcome = outcomes,
             true     = sapply(outcomes, function(o) fmtq(m[ok,o],             d_death(o))),
             reported = sapply(outcomes, function(o) fmtq(rho_draw[ok]*m[ok,o], d_death(o))),
             row.names = NULL)
}))

# averted, median + 95% UI, per vaccine scenario. pct_symp NOW CARRIES ITS 95% UI.
mc_tbl <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- av[[nm]]
  data.frame(timing = lab_timing(nm), arm = lab_arm(nm),
             Infections = fmtq(m[,"infections"]), Symptomatic = fmtq(m[,"symptomatic"]),
             Hospitalisations = fmtq(m[,"hospitalisations"], 1), Deaths = fmtq(m[,"deaths"], 2),
             pct_symp = fmtpct(m[,"symptomatic"] / base_symp_ok),
             row.names = NULL)
}))

# outcomes averted per 100,000 doses (scale-free). Infections averted only by an
# infection-blocking arm, so disease-blocking shows NA there.
per100k_outcomes <- c(infections="Infections", symptomatic="Symptomatic",
                      hospitalisations="Hospitalisations", deaths="Deaths", daly="DALYs")
mc_per100k <- do.call(rbind, lapply(vac_names, function(nm) {
  doses <- G$per_draw[[nm]][ok, "doses"]
  b <- base_pd[ok, , drop = FALSE]; s <- G$per_draw[[nm]][ok, , drop = FALSE]
  cells <- lapply(names(per100k_outcomes), function(o) {
    r <- 1e5 * (b[,o] - s[,o]) / doses; r[!is.finite(r) | r < 0] <- NA; fmtq(r, 0)
  })
  setNames(data.frame(timing = lab_timing(nm), arm = lab_arm(nm), cells,
                      row.names = NULL, check.names = FALSE),
           c("timing","arm", unname(per100k_outcomes)))
}))

scenario_totals <- data.frame(
  scenario = scen_names,
  do.call(rbind, lapply(scen_names, function(nm)
    round(apply(G$per_draw[[nm]][ok, outcomes, drop=FALSE], 2, median, na.rm=TRUE), 1))),
  row.names = NULL, check.names = FALSE)

notes <- data.frame(
  parameter = c("Disease", "Seasonal envelope", "Evaluation window", "R0 interpretation",
                "R0 scenario (FIXED)", "Draws exceeding legacy 1% filter", "Conditioning",
                "Draws summarised (n)", "Vaccine mechanism", "Coverage of 18-59 (median)",
                "VE disease-blocking (median)", "Deaths / CFR", "Severity + DALY source",
                "Reporting rate rho", "REPORTED definition", "Uncertainty"),
  value = c("Mayaro virus (MAYV), Caldas Novas, hypothetical outbreak.",
            "Hybrid: Caldas CHIKV beta_t for the rise/peak + CHIRPS climatological dry-season tail (2026-W10 join), mean-1.",
            sprintf("2025-W24 -> 2026-W22 (weeks %d-%d).", min(G$EVAL_WIN), max(G$EVAL_WIN)),
            "R0 = wet-season PEAK R_eff (envelope rescaled so max = 1).",
            sprintf("'%s': R0 FIXED at %.2f, not sampled. Published MAYV R0 figures are point estimates from different settings/decades (Caicedo 2021: 1.1-1.3 outside Amazon, 2.1-2.9 Amazon; Dodero-Rojas 2020 limits 1.18-3.51), i.e. between-setting heterogeneity rather than uncertainty about one municipality. R0 is varied ACROSS scenarios instead.", G$R0_scenario, G$R0_fixed),
            sprintf("%.1f%% of %d draws exceed attack > %.1f%% of susceptibles. DIAGNOSTIC ONLY -- not used to filter: at fixed R0 the outbreak-size distribution is unimodal, so this threshold would bisect a single continuous distribution.",
                    100*G$frac_over_legacy_thresh, G$N_DRAWS, G$OUTBREAK_ATTACK_THRESH),
            sprintf("NONE. All %d draws are summarised (fixed R0 -> one transmission regime).", length(ok)),
            as.character(length(ok)),
            "Disease-blocking ONLY (VE_inf = 0), pre-outbreak campaign -> infections never averted.",
            sprintf("%.0f%%", 100*median(G$cov_d)),
            sprintf("%.0f%% (sensitivity 25-75%%)", 100*median(G$veb_d)),
            "Zero: no confirmed MAYV-attributable death -> CFR = 0, so deaths & YLL = 0, DALY = YLD.",
            "BORROWED CHIKV (Hyolim Table S4) -- CHIKV-equivalent UPPER bound, not measured MAYV.",
            "Per-draw (Beta, median ~0.25).",
            "REPORTED = rho x TRUE per draw. Severe outcomes (hosp) are usually better ascertained, so their REPORTED values are conservative lower bounds.",
            "Latin-hypercube over gamma/sigma/rho/prop_symp/prior-immunity + vaccine + severity/DALY, propagated jointly. R0 is FIXED (scenario-defining), so its between-setting span is NOT inside these UIs. Prior immunity dominates the spread (r ~ -0.86 with log symptomatic)."),
  stringsAsFactors = FALSE)

# doses actually delivered to the eligible 18-59, and dose wastage. Delivered PRE-
# OUTBREAK, so these are R0-independent (and, as everywhere now, over ALL draws).
# wastage = doses administered to already-immune eligible people (cannot benefit)
#         = 1 - on-target (reached susceptibles) / administered.
doses_wastage <- data.frame(
  timing = lab_timing(G$vac_name), arm = lab_arm(G$vac_name),
  `eligible population (18-59)` = round(G$target_pop_elig),
  `coverage of eligible`        = fmtpct(G$cov_d),
  `doses delivered`             = fmtq(G$doses_deliv, 0),
  `doses on-target (to susceptibles)` = fmtq(G$doses_ontarget, 0),
  `wastage %`                   = fmtpct(1 - G$doses_ontarget / G$doses_deliv),
  check.names = FALSE, row.names = NULL)

write_xlsx(list(notes = notes, baseline_true_reported = base_tbl,
                vaccinated_true_reported = vtr, averted_MC_95UI = mc_tbl,
                averted_per_100k_doses = mc_per100k, doses_wastage = doses_wastage,
                scenario_totals = scenario_totals),
           "MAYV_ca_vacc_outputs.xlsx")

# ============================================================
# WORKBOOK 2 -- MAYV_ca_daly_outputs.xlsx
# ============================================================
daly_by_scenario <- do.call(rbind, lapply(scen_names, function(nm) {
  m <- G$per_draw[[nm]]
  data.frame(timing = lab_timing(nm), arm = lab_arm(nm),
             YLD = fmtq(m[ok,"yld"]), YLL = fmtq(m[ok,"yll"]), DALY = fmtq(m[ok,"daly"]),
             YLD_acute = fmtq(m[ok,"yld_acute"]), YLD_subacute = fmtq(m[ok,"yld_subacute"]),
             YLD_chronic = fmtq(m[ok,"yld_chronic"]), row.names = NULL, check.names = FALSE)
}))
base_daly_ok <- base_pd[ok, "daly"]
daly_averted <- do.call(rbind, lapply(vac_names, function(nm) {
  ad <- base_daly_ok - G$per_draw[[nm]][ok, "daly"]
  data.frame(timing = lab_timing(nm), arm = lab_arm(nm),
             DALY_averted = fmtq(ad),
             pct_DALY = fmtpct(ad / base_daly_ok),          # now with 95% UI
             row.names = NULL, check.names = FALSE)
}))
write_xlsx(list(daly_by_scenario = daly_by_scenario, daly_averted = daly_averted),
           "MAYV_ca_daly_outputs.xlsx")

# ============================================================
# WORKBOOK 3 -- MAYV_ca_nnv_outputs.xlsx
# ============================================================
# NNV = doses / burden averted, per draw (all draws). Same outcome set
# as the CHIKV NNV tab. Infection & Death are NA here (disease-blocking averts no
# infections; MAYV has zero deaths).
nnv_outcomes <- c(infections = "Infection", symptomatic = "Symptomatic case",
                  hospitalisations = "Hospitalisation", deaths = "Death", daly = "DALY")
nnv_tbl <- do.call(rbind, lapply(vac_names, function(nm) {
  doses <- G$per_draw[[nm]][ok, "doses"]
  cells <- lapply(names(nnv_outcomes), function(o) {
    avo <- base_pd[ok, o] - G$per_draw[[nm]][ok, o]
    n   <- doses / avo; n[!is.finite(n) | n < 0] <- NA
    fmtq(n, 0)
  })
  setNames(data.frame(timing = lab_timing(nm), arm = lab_arm(nm), cells,
                      row.names = NULL, check.names = FALSE),
           c("timing","arm", unname(nnv_outcomes)))
}))
write_xlsx(list(nnv = nnv_tbl), "MAYV_ca_nnv_outputs.xlsx")

# ------------------------------------------------------------
cat("Wrote MAYV_ca_vacc_outputs.xlsx (notes, baseline_true_reported, vaccinated_true_reported,\n",
    "     averted_MC_95UI, averted_per_100k_doses, scenario_totals)\n", sep = "")
cat("Wrote MAYV_ca_daly_outputs.xlsx (daly_by_scenario, daly_averted)\n")
cat("Wrote MAYV_ca_nnv_outputs.xlsx  (nnv)\n\n")
cat(sprintf("Fixed R0 = %.2f, all %d draws (no conditioning).  pct symptomatic reduced: %s\n",
            G$R0_fixed, G$N_DRAWS, mc_tbl$pct_symp[1]))
print(mc_tbl, row.names = FALSE)
