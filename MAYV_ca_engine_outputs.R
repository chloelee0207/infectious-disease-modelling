# ============================================================
# MAYV_ca_engine_outputs.R -- comprehensive Excel export for the MAYV engine.
# ------------------------------------------------------------
# Reads MAYV_ca_engine_results.rds and writes caldas_mayv_engine_outputs.xlsx with,
# for BOTH R0 treatments and TRUE + REPORTED versions of every outcome:
#   * sampled-R0  (truncated Lognormal on Dodero limits) -> burden CONDITIONAL on the
#                 outbreak taking off, + disease-blocking vaccine averted/NNV
#   * fixed-R0    (R0_FIX, representative outbreak; every draw is an outbreak) -> burden
#                 + vaccine averted/NNV
# Definitions (stated on the summary sheet):
#   TRUE     = full modelled epidemic quantity.
#   REPORTED = TRUE x rho (per-draw case reporting rate) = surveillance-observed. For
#              symptomatic this equals the model's reported cases (rho*symptomatic).
#              NB severe outcomes (hosp/deaths) are usually better ascertained than mild
#              cases, so their REPORTED values here are conservative lower bounds.
#
# Run order: source("MAYV_ca_lhs.R"); source("MAYV_ca_engine.R"); source(this)
# ============================================================
library(writexl)
if (!exists("fmtq")) source("ca_common.R")   # fmtq(v, d) -> "median (lo - hi)"

if (!file.exists("MAYV_ca_engine_results.rds"))
  stop("MAYV_ca_engine_results.rds not found -- run MAYV_ca_engine.R first.")
G <- readRDS("MAYV_ca_engine_results.rds")

rho_draw <- G$rho_draw
ob       <- G$outbreak                                   # sampled-R0 outbreak draws
base_pd  <- G$per_draw[["No vaccine (baseline)"]]; vac_pd <- G$per_draw[[G$vac_name]]
fb       <- G$fixed_base_pd; fv <- G$fixed_vac_pd        # fixed-R0 baseline / vaccine
all_draws <- seq_len(G$N_DRAWS)

dgt <- c(infections=0, reported=0, symptomatic=0, hospitalisations=1, deaths=2,
         n_nonhosp=0, n_chronic=0, yld_acute=0, yld_chronic=0, yld=0, yll=0, daly=0, doses=0)
d_of <- function(o) if (o %in% names(dgt)) dgt[[o]] else 0

burden_cols <- c("infections","symptomatic","hospitalisations","deaths",
                 "n_nonhosp","n_chronic","yld_acute","yld_chronic","yld","yll","daly")
vax_cols    <- c("symptomatic","hospitalisations","deaths","n_chronic","yld","yll","daly")

# --- no-vaccine burden: TRUE + REPORTED (= true x rho) for each outcome ----------
tr_burden <- function(mat, draws) do.call(rbind, lapply(burden_cols, function(o) {
  tv <- mat[draws, o]; rv <- tv * rho_draw[draws]
  data.frame(outcome = o,
             `true (median 95% UI)`     = fmtq(tv, d_of(o)),
             `reported (median 95% UI)` = fmtq(rv, d_of(o)),
             check.names = FALSE, row.names = NULL)
}))

# --- vaccine: averted (true + reported) and NNV (doses / averted) ----------------
tr_vaccine <- function(base, vac, draws) do.call(rbind, lapply(vax_cols, function(o) {
  av_t <- base[draws, o] - vac[draws, o]; av_r <- av_t * rho_draw[draws]
  doses <- vac[draws, "doses"]
  nnv_t <- doses / av_t; nnv_t[!is.finite(nnv_t) | nnv_t < 0] <- NA
  nnv_r <- doses / av_r; nnv_r[!is.finite(nnv_r) | nnv_r < 0] <- NA
  data.frame(outcome = o,
             `averted true`         = fmtq(av_t, d_of(o)),
             `averted reported`     = fmtq(av_r, d_of(o)),
             `NNV per true averted`     = fmtq(nnv_t, 0),
             `NNV per reported averted` = fmtq(nnv_r, 0),
             check.names = FALSE, row.names = NULL)
}))

# --- per-draw dumps --------------------------------------------------------------
draw_dump <- function(mat, attack, flagcol) data.frame(
  draw = all_draws, attack_pct = round(attack, 3), outbreak = flagcol,
  infections_true = round(mat[, "infections"]), reported_cases = round(mat[, "reported"]),
  symptomatic_true = round(mat[, "symptomatic"]), hospitalisations = round(mat[, "hospitalisations"], 1),
  deaths = round(mat[, "deaths"], 2), daly = round(mat[, "daly"], 1))

summary_df <- data.frame(
  item = c("R0 interpretation", "R0 sampled scenario", "R0 sampled range (Dodero limits)",
           "R0 fixed (representative)", "Monte Carlo draws", "P(outbreak takes off), sampled R0",
           "Outbreak threshold (attack > %)", "Sampled-R0 outbreak draws (n)",
           "Fixed-R0 draws (all take off)", "Evaluation window", "Eligible population (18-59)",
           "Coverage (median)", "VE disease-blocking (median)", "Vaccine mechanism",
           "Severity/DALY source", "TRUE means", "REPORTED means"),
  value = c("wet-season PEAK (truncated Lognormal)", G$R0_scenario, "[1.18, 3.51]",
            sprintf("%.1f", G$R0_FIX), G$N_DRAWS, sprintf("%.1f%%", 100*G$p_outbreak),
            G$OUTBREAK_ATTACK_THRESH, length(ob), G$N_DRAWS,
            sprintf("weeks %d-%d (2025-W40 .. 2026-W38)", min(G$EVAL_WIN), max(G$EVAL_WIN)),
            round(G$target_pop_elig), sprintf("%.0f%%", 100*median(G$cov_d)),
            sprintf("%.0f%%", 100*median(G$veb_d)), "Disease-blocking only, pre-outbreak",
            "Borrowed CHIKV (upper bound)", "full modelled epidemic quantity",
            "true x rho (surveillance-observed; severe outcomes = lower bound)"),
  check.names = FALSE)

sheets <- list(
  summary                = summary_df,
  sampledR0_burden       = tr_burden(base_pd, ob),          # no vaccine | outbreak
  sampledR0_vaccine      = tr_vaccine(base_pd, vac_pd, ob),
  fixedR0_burden         = tr_burden(fb, all_draws),        # no vaccine | fixed R0
  fixedR0_vaccine        = tr_vaccine(fb, fv, all_draws),
  per_draw_sampledR0     = draw_dump(base_pd, G$attack_base, all_draws %in% ob),
  per_draw_fixedR0       = draw_dump(fb, G$attack_fixed, rep(TRUE, G$N_DRAWS)))
write_xlsx(sheets, "caldas_mayv_engine_outputs.xlsx")
cat("Wrote caldas_mayv_engine_outputs.xlsx (sheets:", paste(names(sheets), collapse = ", "), ")\n")
cat(sprintf("  sampled-R0 (P=%.0f%%): true infections %s | reported cases %s\n", 100*G$p_outbreak,
            fmtq(base_pd[ob,"infections"],0), fmtq(base_pd[ob,"reported"],0)))
cat(sprintf("  fixed-R0 (%.1f):       true infections %s | reported cases %s\n", G$R0_FIX,
            fmtq(fb[,"infections"],0), fmtq(fb[,"reported"],0)))
