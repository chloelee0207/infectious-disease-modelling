# ============================================================
# MAYV_ca_burden_outputs.R -- TRANSMISSION-ONLY burden workbook for Caldas Novas MAYV.
# ------------------------------------------------------------
# A deliberately reduced companion to MAYV_ca_outputs.R. It reports ONLY the outcomes that
# fall out of the transmission model itself:
#     infections, symptomatic cases, hospitalisations, deaths
# and drops everything downstream of the borrowed CHIKV disease-progression parameters:
# DALYs, YLD/YLL, disability weights, illness durations, remaining life-years, the recovery
# funnel, and all healthcare costs. There is no MAYV-specific severity or cost-of-illness
# evidence, so those quantities are the CHIKV progression applied to a MAYV case count
# rather than independent MAYV estimates.
#
# WHY THE INTERVALS DO NOT NARROW. Dropping those parameters does NOT tighten the 95% UIs
# on infections or symptomatic cases, because they never entered them in the first place:
#   * infections  comes straight out of the SEIR (R0, gamma, sigma, prior immunity, seed)
#   * symptomatic = infections x prop_symp, both transmission-layer quantities
#   * reported    = rho x symptomatic, applied after the simulation
#   * the DALY/cost parameters enter outcome_one() only through yld_* and yll, which are
#     OUTPUTS -- nothing feeds back into the epidemic
# The width comes from the transmission layer alone. Against the current draws, log
# infections correlates +0.74 with sampled R0 and -0.52 with prior immunity; prop_symp and
# rho are ~0. Because outbreak size is a steep convex function of R0, sampling it across a
# range as wide as 2.1-2.9 produces a heavily right-skewed distribution -- hence the ~990x
# ratio between the 97.5th and 2.5th percentiles. That is a property of the R0 prior, not
# of the burden layer, and no reduction in downstream parameters can shrink it.
#
# TWO PARAMETERS ARE STILL NEEDED and are retained:
#   hosp_d  hospitalisation rate among symptomatic -- MAYV-specific, 5% (4-6%)
#   cfr     case fatality -- but MAYV_ZERO_DEATHS sets it to 0, so deaths are identically
#           zero and the deaths rows are reported as an explicit zero rather than omitted.
#
# Run after: MAYV_ca_engine.R  (does NOT need MAYV_ca_costs.R)
# Output:    MAYV_ca_burden_outputs.xlsx
# ============================================================
suppressMessages({library(writexl)})
if (!exists("fmtq")) source("ca_common.R")

if (!file.exists("MAYV_ca_engine_results.rds"))
  stop("MAYV_ca_engine_results.rds not found -- run MAYV_ca_engine.R first.")
G  <- readRDS("MAYV_ca_engine_results.rds")
ok <- G$outbreak                      # ALL draws (no take-off conditioning)
ND <- G$N_DRAWS
rho_draw <- G$rho_draw
base_pd  <- G$per_draw[["No vaccine (baseline)"]]
vac_pd   <- G$per_draw[[G$vac_name]]
doses    <- G$target_pop_elig * as.numeric(G$cov_d)

# Transmission-derived outcomes only. `reported` is carried as a separate COLUMN
# (rho x true) rather than as its own row, matching MAYV_ca_outputs.R.
OUT  <- c("infections", "symptomatic", "hospitalisations", "deaths")
LAB  <- c(infections = "Infections", symptomatic = "Symptomatic cases",
          hospitalisations = "Hospitalisations", deaths = "Deaths")
dp   <- function(o) if (o == "deaths") 2 else 0     # deaths are 0 for MAYV; keep decimals

# ---- 1. baseline, true vs reported ---------------------------------------
baseline <- do.call(rbind, lapply(OUT, function(o) {
  tv <- base_pd[ok, o]
  data.frame(outcome = LAB[[o]], true = fmtq(tv, dp(o)),
             reported = fmtq(rho_draw[ok] * tv, dp(o)), row.names = NULL) }))

# ---- 2. vaccinated arm, true vs reported ---------------------------------
vaccinated <- do.call(rbind, lapply(OUT, function(o) {
  tv <- vac_pd[ok, o]
  data.frame(scenario = G$vac_name, outcome = LAB[[o]], true = fmtq(tv, dp(o)),
             reported = fmtq(rho_draw[ok] * tv, dp(o)), row.names = NULL) }))

# ---- 3. averted, with the % reduction ------------------------------------
pctf <- function(a, b) { r <- 100 * a / b; r[!is.finite(r)] <- NA
  q <- quantile(r, c(.5, .025, .975), na.rm = TRUE)
  sprintf("%.1f%% (%.1f - %.1f%%)", q[1], q[2], q[3]) }
averted <- do.call(rbind, lapply(OUT, function(o) {
  av <- base_pd[ok, o] - vac_pd[ok, o]
  data.frame(outcome = LAB[[o]], averted = fmtq(av, dp(o)),
             averted_reported = fmtq(rho_draw[ok] * av, dp(o)),
             pct_reduced = pctf(av, base_pd[ok, o]), row.names = NULL) }))

# ---- 4. averted per 100,000 doses ----------------------------------------
# Scale-free: doses are the shared input across pathogens, so this is comparable with the
# CHIKV figures without dividing one model by the other.
per100k <- do.call(rbind, lapply(OUT, function(o) {
  r <- 1e5 * (base_pd[ok, o] - vac_pd[ok, o]) / doses[ok]
  r[!is.finite(r) | r < 0] <- NA
  data.frame(outcome = LAB[[o]], averted_per_100k_doses = fmtq(r, dp(o)), row.names = NULL) }))

# ---- 5. doses and wastage -------------------------------------------------
wastage <- 1 - G$doses_ontarget[ok] / G$doses_deliv[ok]
doses_tbl <- data.frame(
  quantity = c("Doses administered (eligible 18-59)", "Doses reaching susceptibles",
               "Dose wastage (%)"),
  value = c(fmtq(G$doses_deliv[ok], 0), fmtq(G$doses_ontarget[ok], 0),
            sprintf("%.1f%% (%.1f - %.1f%%)", 100*median(wastage),
                    100*quantile(wastage, .025), 100*quantile(wastage, .975))),
  row.names = NULL)

# ---- 6. notes -------------------------------------------------------------
notes <- data.frame(
  item = c("Scope", "Excluded", "Why the intervals do not narrow",
           "What drives the width", "Deaths", "Hospitalisation rate",
           "R0", "Reported", "Draws", "Conditioning"),
  detail = c(
  "Transmission-derived burden only: infections, symptomatic cases, hospitalisations, deaths. Baseline, vaccinated, averted, averted per 100,000 doses.",
  "DALYs, YLD/YLL, disability weights, illness durations, remaining life-years, the recovery funnel and ALL healthcare costs. Those depend on CHIKV disease-progression parameters borrowed wholesale; there is no MAYV-specific severity or cost-of-illness evidence, so they are the CHIKV progression applied to a MAYV case count, not independent MAYV estimates.",
  "Removing those parameters does NOT tighten these intervals, because they never entered them. Infections come straight out of the SEIR; symptomatic = infections x prop_symp; reported = rho x symptomatic. The DALY/cost parameters enter only through YLD and YLL, which are outputs -- nothing feeds back into the epidemic.",
  sprintf("The transmission layer alone. Correlation of log infections with sampled R0 is about +0.74 and with prior immunity about -0.52; prop_symp and rho are ~0. Outbreak size is a steep convex function of R0, so sampling it across the scenario range gives a heavily right-skewed distribution -- the 97.5th/2.5th percentile ratio for infections is about %.0fx.", quantile(base_pd[ok,"infections"], .975)/max(quantile(base_pd[ok,"infections"], .025), 1e-9)),
  "Identically zero. No MAYV-attributable death has been confirmed, so the engine sets CFR = 0 (MAYV_ZERO_DEATHS). The deaths rows are kept as an explicit zero rather than dropped.",
  "5% (4-6%), MAYV-specific -- a separate row in disease_progression.xlsx, not the CHIKV 4% (3-5%).",
  sprintf("SAMPLED from the '%s' scenario range %.1f-%.1f (lognormal, endpoints as 2.5th/97.5th percentiles).", G$R0_scenario, G$R0_lo, G$R0_hi),
  "REPORTED = rho x TRUE, applied per draw after the simulation. rho never enters transmission.",
  sprintf("%d Latin hypercube draws, propagated jointly with the vaccine parameters.", ND),
  "None. Every draw is reported; no take-off filter."),
  stringsAsFactors = FALSE)

write_xlsx(list(notes = notes,
                baseline_true_reported = baseline,
                vaccinated_true_reported = vaccinated,
                averted_95UI = averted,
                averted_per_100k_doses = per100k,
                doses_wastage = doses_tbl),
           "MAYV_ca_burden_outputs.xlsx")
cat("Wrote MAYV_ca_burden_outputs.xlsx (transmission-only: no DALYs, no costs)\n")
cat(sprintf("  scenario '%s', R0 sampled %.1f-%.1f, %d draws\n\n",
            G$R0_scenario, G$R0_lo, G$R0_hi, ND))
print(baseline, row.names = FALSE)
cat("\n"); print(averted, row.names = FALSE)
cat("\n"); print(per100k, row.names = FALSE)
