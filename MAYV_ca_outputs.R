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
#     high = 2.1-2.9), sampled per draw. Every figure is summarised over ALL draws
#     (G$outbreak is every row); no take-off conditioning is applied.
#     The 95% UIs carry natural-history / reporting / symptomatic-fraction / prior-immunity
#     / vaccine / severity-DALY uncertainty, NOT the between-setting R0 span.
#   * ONE vaccine scenario: pre-outbreak, DISEASE-BLOCKING ONLY (VE_inf = 0), so
#     infections are never averted (Infections = 0, and Infection NNV = NA).
#   * DEATHS = 0: no confirmed MAYV-attributable death, so CFR = 0 in the engine ->
#     deaths and YLL are zero and DALY = YLD.
#   * Severity/DALY parameters are BORROWED from CHIKV (upper bound), and the seasonal
#     envelope is the hybrid CHIKV-beta + dry-season envelope (2025-W24 -> 2026-W22).
#
# Run order: source("MAYV_ca_lhs.R"); source("MAYV_ca_engine.R");
#            source("MAYV_ca_costs.R"); source(this)
#            (the cost layer reads only the engine .rds, and the residual-burden
#             figure needs its per-draw direct medical costs.)
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
            sprintf("'%s': R0 SAMPLED from %.1f-%.1f (lognormal, endpoints as 2.5th/97.5th percentiles; median %.2f). Both scenario ranges are Caicedo et al. 2021 -- low = outside the Amazon basin, high = Amazon basin applied to Goias as a PEAK R0 -- so this is within-source uncertainty, not a mix across settings.", G$R0_scenario, G$R0_lo, G$R0_hi, G$R0_fixed),
            sprintf("%.1f%% of %d draws exceed attack > %.1f%% of susceptibles. DIAGNOSTIC ONLY -- not used to filter: at fixed R0 the outbreak-size distribution is unimodal, so this threshold would bisect a single continuous distribution.",
                    100*G$frac_over_legacy_thresh, G$N_DRAWS, G$OUTBREAK_ATTACK_THRESH),
            sprintf("NONE. All %d draws are summarised (fixed R0 -> one transmission regime).", length(ok)),
            as.character(length(ok)),
            "Disease-blocking ONLY (VE_inf = 0), pre-outbreak campaign -> infections never averted.",
            sprintf("%.0f%%", 100*median(G$cov_d)),
            sprintf("%.1f%% mean (Kostecki et al. 2026: 4 of 12 CHIKV patients cross-neutralised MAYV; Beta(4,8), 95%% UI 10.9-61.0%%)", 100*mean(G$veb_d)),
            "Zero: no confirmed MAYV-attributable death -> CFR = 0, so deaths & YLL = 0, DALY = YLD.",
            "BORROWED CHIKV (Hyolim Table S4) -- CHIKV-equivalent UPPER bound, not measured MAYV.",
            "Per-draw (Beta, median ~0.25).",
            "REPORTED = rho x TRUE per draw. Severe outcomes (hosp) are usually better ascertained, so their REPORTED values are conservative lower bounds.",
            "Latin-hypercube over gamma/sigma/rho/prop_symp/prior-immunity/R0 + vaccine + severity/DALY, propagated jointly. R0 IS sampled within the scenario range, so its span IS inside these UIs and dominates them -- outbreak size is a steep convex function of R0, so the intervals are wide and right-skewed. Read the deciles, not just the median."),
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
# Residual burden under PRE-OUTBREAK vaccination, as a % of no vaccination.
# Two bars per panel: unvaccinated counterfactual vs vaccinated. Ratios are formed
# per draw (scenario / baseline) and then summarised, so the UI keeps the
# baseline-scenario pairing. Styled to match CHIKV_ca_residual_burden.png so the two
# can be merged into one figure.
#
# Deaths are omitted: MAYV_ZERO_DEATHS = TRUE sets the MAYV CFR to 0, so baseline
# deaths are identically 0 and the ratio is undefined. DALY = YLD for the same reason.
#
# Healthcare cost is TOTAL direct medical cost -- inpatient plus acute, sub-acute and
# chronic outpatient -- matching the "costs averted" definition in Table 2. It is read
# from the cost layer because the outpatient components need the recovery funnel, which
# only MAYV_ca_costs.R evaluates.
#
# NOTE on interpretation: every cost component is linear in case counts that all descend
# from total symptomatic through rates shared by both arms within a draw (a SINGLE
# all-ages hosp_rate, the recovery funnel, and the sampled unit costs). Those factors
# cancel in the ratio, so this PERCENTAGE equals the symptomatic-case percentage to ~0.1
# pp. It is the absolute R$ averted, not the %, that this panel adds over a case count.
# Deaths do differ (90% vs 84%) because cfr_vec IS age-specific, so the vaccine's
# 18-59 targeting shifts the age mix of fatal cases.
# ------------------------------------------------------------
library(ggplot2)
stopifnot(max(base_pd[ok, "deaths"]) == 0)        # guard the "no deaths panel" claim
stopifnot(length(vac_names) == 1)                 # single pre-outbreak disease-blocking arm

pct_of_base <- function(num, den) { r <- 100 * num / den; r[!is.finite(r)] <- NA; r }
res_rows <- list()
add_res <- function(outcome, scen, v) res_rows[[length(res_rows)+1]] <<- data.frame(
  outcome = outcome, scenario = scen,
  med = median(v, na.rm = TRUE), lo = quantile(v, .025, na.rm = TRUE),
  hi = quantile(v, .975, na.rm = TRUE), row.names = NULL)

if (!file.exists("MAYV_ca_costs.rds"))
  stop("MAYV_ca_costs.rds not found -- run MAYV_ca_costs.R before this script ",
       "(the residual-burden figure needs total direct medical cost).")
mayv_cost <- readRDS("MAYV_ca_costs.rds")
# GUARD: the cost layer is written by MAYV_ca_costs.R, which must run BEFORE this
# script (see the run order above). Without this check, running them out of order --
# or after switching R0_SCENARIO -- silently pairs one scenario's costs with another's
# burden, which shows up as a ~0% healthcare-cost reduction instead of ~8%.
if (!identical(mayv_cost$R0_scenario, G$R0_scenario))
  stop("MAYV_ca_costs.rds is from R0 scenario '", mayv_cost$R0_scenario %||% "<untagged>",
       "' but the engine results are '", G$R0_scenario,
       "'. Re-run MAYV_ca_costs.R for this scenario BEFORE MAYV_ca_outputs.R.")
cost_pd <- mayv_cost$cost_pd
stopifnot(nrow(cost_pd[[1]]) == nrow(base_pd))    # cost draws align with engine draws

# outcome -> (source, column). "engine" = per-draw burden; "cost" = per-draw cost layer.
resid_outcomes <- list(
  list(lab = "Cumulative DALYs", src = "engine", col = "daly"),
  list(lab = "Healthcare cost",  src = "cost",   col = "total_direct_medical"))
nm <- vac_names[1]
for (o in resid_outcomes) {
  v <- if (o$src == "engine") pct_of_base(G$per_draw[[nm]][ok, o$col], base_pd[ok, o$col])
       else pct_of_base(cost_pd[[nm]][ok, o$col],
                        cost_pd[["No vaccine (baseline)"]][ok, o$col])
  add_res(o$lab, "No vaccination", 100)
  add_res(o$lab, "Vaccination",    v)
}
resid <- do.call(rbind, res_rows)
resid$scenario <- factor(resid$scenario, levels = c("No vaccination", "Vaccination"))
resid$outcome  <- factor(resid$outcome, levels = sapply(resid_outcomes, `[[`, "lab"))
resid$arm      <- factor("Disease-blocking")      # MAYV has no infection-blocking arm

p_resid <- ggplot(resid, aes(scenario, med, fill = scenario)) +
  geom_col(width = .6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, linewidth = .35) +
  facet_grid(outcome ~ arm) +
  scale_fill_manual(values = c("No vaccination" = "grey60", "Vaccination" = "#4e79a7"),
                    name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, NA)) +
  labs(x = NULL, y = "Cumulative burden (% of no vaccination)") +
  theme_bw(11) +
  theme(axis.text.x = element_text(size = 9), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank())
print(p_resid)
ggsave("MAYV_ca_residual_burden.png", p_resid, width = 3.6, height = 5.4, dpi = 130)
# % reduction from no vaccination, the complement of the residual columns. The interval
# BOUNDS SWAP: a draw with a high residual burden is a draw with a small reduction, so
# red_lo is 100 - hi and red_hi is 100 - lo. Taking 100 - lo as the lower bound would
# report the interval backwards.
resid$red_med <- 100 - resid$med
resid$red_lo  <- 100 - resid$hi
resid$red_hi  <- 100 - resid$lo
write_xlsx(list(residual_burden_pct = resid), "MAYV_ca_residual_burden.xlsx")
saveRDS(resid, "MAYV_ca_residual_burden.rds")     # for the merged CHIKV|MAYV figure
cat("Saved MAYV_ca_residual_burden.png and .xlsx (burden as % of no vaccination;\n",
    "     deaths panel omitted -- MAYV CFR is fixed at 0).\n", sep = "")

# ------------------------------------------------------------
cat("Wrote MAYV_ca_vacc_outputs.xlsx (notes, baseline_true_reported, vaccinated_true_reported,\n",
    "     averted_MC_95UI, averted_per_100k_doses, scenario_totals)\n", sep = "")
cat("Wrote MAYV_ca_daly_outputs.xlsx (daly_by_scenario, daly_averted)\n")
cat("Wrote MAYV_ca_nnv_outputs.xlsx  (nnv)\n\n")
cat(sprintf("Fixed R0 = %.2f, all %d draws (no conditioning).  pct symptomatic reduced: %s\n",
            G$R0_fixed, G$N_DRAWS, mc_tbl$pct_symp[1]))
print(mc_tbl, row.names = FALSE)
