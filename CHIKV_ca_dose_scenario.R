# ============================================================
# CHIKV_ca_dose_scenario.R
# ------------------------------------------------------------
# Compares the modelled campaign against the DOSE ALLOCATION actually promised to
# Caldas Novas: 4,670 doses, i.e. 7.11% of the eligible 18-59 group (4.37% of the total
# population), versus the ~30% of eligible the main analysis samples.
#
# Coverage is FIXED in the dose-constrained run rather than sampled, because a promised
# allocation is a known quantity. Everything else -- transmission, reporting, severity,
# DALY and the remaining vaccine parameters -- stays sampled, so the 95% UI reflects
# epidemiological and disease uncertainty at that allocation instead of mixing in
# uncertainty about how many doses a programme might choose to buy.
#
# Inputs (produce them first):
#   CHIKV_ca_engine.R                                  -> CHIKV_ca_engine_results.rds
#   CHIKV_DOSES <- 4670; source("CHIKV_ca_engine.R")   -> CHIKV_ca_engine_results_doses4670.rds
# Output: CHIKV_ca_dose_scenario.xlsx
# ============================================================
suppressMessages({library(writexl)})

DOSES <- 4670
f_main <- "CHIKV_ca_engine_results.rds"
f_dose <- sprintf("CHIKV_ca_engine_results_doses%d.rds", DOSES)
for (f in c(f_main, f_dose))
  if (!file.exists(f)) stop("missing ", f, " -- see the header for how to produce it.")
A <- readRDS(f_main); B <- readRDS(f_dose)

q3  <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
fmt <- function(x, d = 0) { q <- q3(x); sprintf("%s (%s-%s)",
        formatC(round(q[1], d), big.mark = ",", format = "f", digits = d),
        formatC(round(q[2], d), big.mark = ",", format = "f", digits = d),
        formatC(round(q[3], d), big.mark = ",", format = "f", digits = d)) }
pfmt <- function(x) { q <- q3(x); d <- if (q[1] < 1) 3 else 1
  sprintf("%.*f%% (%.*f-%.*f%%)", d, q[1], d, q[2], d, q[3]) }

row_for <- function(G, lab, tn, an) {
  b  <- G$per_draw[["No vaccine (baseline)"]]
  v  <- G$per_draw[[paste0(tn, " | ", an)]]
  data.frame(coverage_scenario = lab, timing = tn, arm = an,
             doses            = fmt(v[, "doses"]),
             symptomatic_averted = fmt(b[, "symptomatic"] - v[, "symptomatic"]),
             pct_symptomatic_averted = pfmt(100*(b[, "symptomatic"] - v[, "symptomatic"])/b[, "symptomatic"]),
             hospitalisations_averted = fmt(b[, "hospitalisations"] - v[, "hospitalisations"], 1),
             deaths_averted   = fmt(b[, "deaths"] - v[, "deaths"], 2),
             DALYs_averted    = fmt(b[, "daly"] - v[, "daly"], 1),
             averted_per_100k_doses = fmt(1e5*(b[, "symptomatic"] - v[, "symptomatic"])/v[, "doses"]),
             stringsAsFactors = FALSE)
}
TIMINGS <- names(A$timings); ARMS <- A$arm_names
tbl <- do.call(rbind, c(
  lapply(TIMINGS, function(t) do.call(rbind, lapply(ARMS, function(a)
    row_for(A, sprintf("sampled ~30%% of eligible (%s doses)",
                       format(round(median(A$per_draw[[paste0(t," | ",a)]][,"doses"])), big.mark=",")), t, a)))),
  lapply(TIMINGS, function(t) do.call(rbind, lapply(ARMS, function(a)
    row_for(B, sprintf("FIXED %s doses (7.11%% of eligible)", format(DOSES, big.mark=",")), t, a))))))

# ratio of the two, per draw, so the pairing is kept
ratio <- do.call(rbind, lapply(TIMINGS, function(t) do.call(rbind, lapply(ARMS, function(a) {
  ba <- A$per_draw[["No vaccine (baseline)"]][, "symptomatic"]
  bb <- B$per_draw[["No vaccine (baseline)"]][, "symptomatic"]
  pa <- 100*(ba - A$per_draw[[paste0(t," | ",a)]][, "symptomatic"])/ba
  pb <- 100*(bb - B$per_draw[[paste0(t," | ",a)]][, "symptomatic"])/bb
  data.frame(timing = t, arm = a, pct_sampled_coverage = pfmt(pa),
             pct_fixed_4670_doses = pfmt(pb), fold_lower = fmt(pa/pb, 2),
             stringsAsFactors = FALSE) }))))

notes <- data.frame(item = c("Question", "Coverage", "What is still sampled", "Doses",
                             "Per 100,000 doses", "Baseline"),
  detail = c(
  sprintf("What the %s doses promised to Caldas Novas would have achieved, against the ~30%% of eligible coverage the main analysis samples.", format(DOSES, big.mark=",")),
  sprintf("%s doses = %.2f%% of the eligible 18-59 group (%.2f%% of the total population of 106,820). Coverage is FIXED at that value, not sampled, because a promised allocation is a known quantity.", format(DOSES, big.mark=","), 100*DOSES/A$target_pop_elig, 100*DOSES/106820),
  "Transmission (beta, FOI, gamma, sigma, reporting rate, symptomatic fraction), vaccine efficacy, delivery speed, deployment delay, severity, CFR and all DALY parameters. The interval therefore reflects epidemiological and disease uncertainty at a fixed allocation.",
  "Doses ADMINISTERED, which is below the allocation where the rollout does not complete inside the window.",
  "Scale-free comparison: with fewer doses the per-dose return is higher, because early doses avert more than later ones once the susceptible pool starts to deplete.",
  "Each arm is compared with its own no-vaccination baseline from the same run, so the percentages are internally consistent."),
  stringsAsFactors = FALSE)

write_xlsx(list(by_scenario = tbl, sampled_vs_fixed = ratio, notes = notes),
           "CHIKV_ca_dose_scenario.xlsx")
cat("Wrote CHIKV_ca_dose_scenario.xlsx\n")
print(ratio, row.names = FALSE)
