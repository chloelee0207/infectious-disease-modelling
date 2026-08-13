# ============================================================
# CHIKV_ca_coverage_timing.R
# ------------------------------------------------------------
# Symptomatic cases averted PER 100,000 DOSES, by campaign timing x vaccine coverage,
# with 95% UIs propagated over the 1,000-draw ensemble.
#
# Coverage is FIXED in each run rather than sampled: the table is indexed BY coverage, so
# sampling it would put the row label and the parameter in disagreement. Everything else
# stays sampled, so each cell's interval reflects epidemiological, vaccine and disease
# uncertainty at that coverage, which is what makes rows comparable.
#
# Per-dose return FALLS as coverage rises: the marginal dose at high coverage protects
# someone the epidemic might not have reached anyway. Absolute impact still rises. The
# two columns therefore move in opposite directions, which is the point of the table.
#
# Inputs: one engine run per coverage level, produced by
#   for (cv in COVS) { CHIKV_FIXED_COV <- cv; source("CHIKV_ca_engine.R") }
# Output: CHIKV_ca_coverage_timing.xlsx
# ============================================================
suppressMessages(library(writexl))

files <- sort(list.files(pattern = "^CHIKV_ca_engine_results_cov[0-9]+\\.rds$"))
if (!length(files)) stop("no fixed-coverage engine results found -- see the header.")

q3   <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
fmt  <- function(x, d = 0) { q <- q3(x); sprintf("%s (%s-%s)",
         formatC(round(q[1], d), big.mark = ",", format = "f", digits = d),
         formatC(round(q[2], d), big.mark = ",", format = "f", digits = d),
         formatC(round(q[3], d), big.mark = ",", format = "f", digits = d)) }
pfmt <- function(x) { q <- q3(x); d <- if (q[1] < 1) 3 else 1
  sprintf("%.*f%% (%.*f-%.*f%%)", d, q[1], d, q[2], d, q[3]) }

rows <- list()
for (f in files) {
  E <- readRDS(f); cov <- E$cov_d[1]
  b <- E$per_draw[["No vaccine (baseline)"]][, "symptomatic"]
  for (tn in names(E$timings)) for (an in E$arm_names) {
    v  <- E$per_draw[[paste0(tn, " | ", an)]]
    av <- b - v[, "symptomatic"]
    rows[[length(rows) + 1]] <- data.frame(
      coverage_pct_eligible = sprintf("%.2f%%", 100*cov),
      doses_median          = formatC(round(median(v[, "doses"])), big.mark = ",", format = "d"),
      timing = tn, arm = an,
      symptomatic_averted            = fmt(av),
      pct_symptomatic_averted        = pfmt(100*av/b),
      averted_per_100k_doses         = fmt(1e5*av/v[, "doses"]),
      DALYs_averted_per_100k_doses   = fmt(1e5*(E$per_draw[["No vaccine (baseline)"]][, "daly"] -
                                                v[, "daly"])/v[, "doses"], 1),
      cov_num = cov, stringsAsFactors = FALSE)
  }
}
tbl <- do.call(rbind, rows)
tbl <- tbl[order(tbl$cov_num, tbl$timing, tbl$arm), ]
tbl$cov_num <- NULL

notes <- data.frame(item = c("Unit", "Coverage", "Why coverage is fixed", "Doses",
                             "Direction of the two columns", "Timings"),
  detail = c(
  "Symptomatic cases (and DALYs) averted per 100,000 doses ADMINISTERED, median (95% UI), formed per draw.",
  sprintf("Coverage is of the ELIGIBLE 18-59 group (65,658 people), not of the total population (106,820). Levels run: %s.",
          paste(sprintf("%.2f%%", 100*sort(unique(sapply(files, function(f) readRDS(f)$cov_d[1])))), collapse = ", ")),
  "Fixed, not sampled. The table is indexed by coverage, so sampling it would put the row label and the parameter in disagreement. Every other input stays sampled, so each interval reflects epidemiological, vaccine and disease uncertainty at that coverage.",
  "Doses ADMINISTERED, which can fall below coverage x eligible where the rollout does not finish inside the 52-week window. Roughly 10% of doses reach people already immune and avert nothing.",
  "Per-dose return FALLS as coverage rises while absolute impact RISES: the marginal dose at high coverage protects someone the epidemic might not have reached. Read the two columns together.",
  "The campaign week is the DECISION week; dosing begins after the sampled 1-3 week deployment delay."),
  stringsAsFactors = FALSE)

write_xlsx(list(per_100k_doses = tbl, notes = notes), "CHIKV_ca_coverage_timing.xlsx")
cat("Wrote CHIKV_ca_coverage_timing.xlsx --", nrow(tbl), "rows\n")
