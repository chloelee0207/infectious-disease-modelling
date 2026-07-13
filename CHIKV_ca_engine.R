# ============================================================
# CHIKV_ca_engine.R -- Caldas Novas CHIKV UNIFIED Monte Carlo engine.
# ------------------------------------------------------------
# ONE Latin-hypercube uncertainty propagation from start to end. A single set of
# N_DRAWS parameter draws (transmission resampled from the LHS ensemble + vaccine
# efficacy, delivery, delay, SAMPLED coverage ~ Beta(30%,20-40), severity/CFR and
# every DALY input) is drawn ONCE; each draw is then run through the SEIRV for the
# baseline and every vaccination scenario. EVERY outcome is derived from those same
# runs, so they are mutually consistent draw-for-draw:
#     infections, symptomatic, hospitalisations, deaths        (burden)
#     hospitalised / non-hospitalised / chronic case counts    (severity phases)
#     YLD (acute, chronic), YLL, DALY                          (health loss)
#     doses delivered, and NNV = doses / burden averted        (per-draw ratio)
# NNV is NOT a separate experiment -- it is doses_i / (baseline_i - scenario_i) at
# that draw's sampled coverage, so its uncertainty is the SAME propagation as the
# burden/DALY numbers. Coverage is sampled (not a fixed sweep) per user request.
#
# Saves BOTH per-draw matrices and aggregated median/95% UI to
# CHIKV_ca_engine_results.rds. Severity-phase COUNTS are saved so a cost layer
# (cost to avert a hospitalisation, phase-specific costs, ICERs) can be added later
# as a cheap post-hoc multiply -- no SEIR re-run.
#
# The SEIRV simulator and DALY loader now live in ca_common.R (single source of
# truth); this engine is their only caller.
#
# Run order:  source("ca_common.R"); source("CHIKV_ca_lhs.R")  # ensemble (slow, once)
#             source("CHIKV_ca_engine.R")                        # this file
# ============================================================
library(dplyr); library(tidyr)
source("ca_common.R")   # fmtq, qs, burden, load_burden_params, load_caldas_age_cases,
                        # compute_age_weight, seirv_vaccinated, load_daly_params

# ---- knobs --------------------------------------------------
N_DRAWS    <- 1000
PHASE_MODE <- "hosp_severity"   # YLD severity axis (see daly_counts); matches CHIKV_ca_daly.R
SEV_SPLIT  <- c(mild_mod = 0.47, severe = 0.53)   # used only by PHASE_MODE == "acute_split"
set.seed(2030)

# ------------------------------------------------------------
# 0. Load the LHS ensemble
# ------------------------------------------------------------
if (!file.exists("CHIKV_ca_lhs_ensemble.rds"))
  stop("CHIKV_ca_lhs_ensemble.rds not found -- run CHIKV_ca_lhs.R first.")
E <- readRDS("CHIKV_ca_lhs_ensemble.rds")
N <- E$N; A <- E$A; age_df <- E$age_df; age_mid <- E$age_mid
observed_cases <- E$observed_cases; caldas_obs <- E$caldas_obs
week_1_cases <- E$week_1_cases; T_data <- E$T_weeks
n_ens <- nrow(E$beta); E0 <- rep(0, A)
cat(sprintf("Loaded ensemble: %d feasible transmission draws; running %d-draw unified MC.\n",
            n_ens, N_DRAWS))

# ------------------------------------------------------------
# 1. Horizon + per-draw immunity/seed
# ------------------------------------------------------------
EXTEND <- 26; T_sim <- T_data + EXTEND
extend_beta <- function(b) c(b, rep(b[length(b)], EXTEND))
draw_immunity <- function(foi) 1 - exp(-foi * age_mid)
draw_I0 <- function(foi, rho, gamma, ps) {
  Rimm <- draw_immunity(foi); sf <- N*(1-Rimm)/sum(N*(1-Rimm))
  round(((week_1_cases/rho/ps)/gamma) * sf)
}

# ------------------------------------------------------------
# 2. Vaccine programme + scenarios (52-week index) + 52-wk evaluation window
# ------------------------------------------------------------
target_age <- rep(0, A); target_age[c(4,5,6,7,8)] <- 1     # eligible = 18-59
target_pop_elig <- sum(N[target_age == 1])
immun_delay <- 2
idx_of <- function(yr, wk) caldas_obs$week_index[caldas_obs$Year == yr & caldas_obs$week == wk]
start_s1 <- idx_of(2026, 16)     # IXCHIQ real rollout
start_s2 <- idx_of(2026, 1)      # start of 2026
start_s3 <- idx_of(2025, 40)     # pre-outbreak
timings <- list("actual rollout" = start_s1, "start of 2026" = start_s2, "pre-outbreak" = start_s3)
arm_names <- c("Disease-blocking", "Disease + infection blocking")

# 52-epi-week evaluation window anchored at pre-outbreak implementation
# (2025-W40 -> 2026-W38); all burden/DALY/NNV accrued only within it.
EVAL_WIN <- start_s3:(start_s3 + 51)
stopifnot(length(EVAL_WIN) == 52, max(EVAL_WIN) <= T_sim)
cat(sprintf("52-week evaluation window: index %d-%d (2025-W40 -> 2026-W38).\n",
            min(EVAL_WIN), max(EVAL_WIN)))

scen <- list(list(name="No vaccine (baseline)", timing="No vaccine", arm="No vaccine",
                  start=NA_integer_, type="base"))
for (tn in names(timings)) for (an in arm_names)
  scen[[length(scen)+1]] <- list(name=paste0(tn," | ",an), timing=tn, arm=an,
                                 start=timings[[tn]], type=if(an=="Disease-blocking")"disb" else "both")
scen_names <- vapply(scen, function(s) s$name, character(1))
vac_names  <- setdiff(scen_names, "No vaccine (baseline)")

# ------------------------------------------------------------
# 3. Severity params + observed-age death correction + DALY params
# ------------------------------------------------------------
invisible(list2env(load_burden_params(A), globalenv()))     # ps_*, hosp_*, cfr_*, age_to_band, cfr_vec
obs_band_prop <- load_caldas_age_cases()$obs_band_prop
Rimm_base <- draw_immunity(E$base_foi)
I0_base   <- draw_I0(E$base_foi, E$base_rho, E$base_gamma, E$base_prop_symp)
out_base_pt <- seirv_vaccinated(T_sim, A, N, Rimm_base, I0_base, E0, extend_beta(E$base_beta),
                 E$base_sigma, E$base_gamma, E$base_rho, target_age, 0, 0.10, start_s3,
                 0, 0, immun_delay, prop_symp = E$base_prop_symp)
age_weight <- compute_age_weight(rowSums(out_base_pt$new_infections), obs_band_prop, age_to_band)
young_idx <- which(age_to_band <= 4); old_idx <- which(age_to_band >= 5)
dp <- load_daly_params()

# ------------------------------------------------------------
# 4. Outcome extractor: one SEIR run -> all outcomes within the 52-wk window.
#    Severity phases (hospitalised / non-hosp / chronic) are saved as COUNTS so a
#    cost layer can multiply unit costs onto them later without re-running.
# ------------------------------------------------------------
outcome_one <- function(out, covv, hosp_j, cfr_j, le_band,
                        dwmm, dwsv, dwch, dumm, dusv, duch, p14y, p14o, p90y, p90o) {
  infections <- sum(out$new_infections[, EVAL_WIN, drop = FALSE])
  symp_age   <- rowSums(out$new_symptomatic[, EVAL_WIN, drop = FALSE])
  symp_w  <- symp_age * age_weight
  symp_dw <- if (sum(symp_w) > 0) symp_w * (sum(symp_age)/sum(symp_w)) else symp_age
  st <- sum(symp_dw); sy <- sum(symp_dw[young_idx]); so <- sum(symp_dw[old_idx])

  n_hosp    <- st * hosp_j                     # hospitalised (= severe acute in hosp_severity)
  n_nonhosp <- st * (1 - hosp_j)               # non-hospitalised (mild/mod acute)
  n_chronic <- sy*(1-p14y)*(1-p90y) + so*(1-p14o)*(1-p90o)   # not recovered by 90d

  if (PHASE_MODE == "hosp_severity") {
    yld_ac  <- n_nonhosp*dwmm*dumm + n_hosp*dwsv*dusv
    yld_chr <- n_chronic * dwch * duch
  } else if (PHASE_MODE == "three_phase") {
    yld_ac  <- st*dwmm*dumm + (sy*(1-p14y) + so*(1-p14o))*dwsv*dusv
    yld_chr <- n_chronic * dwch * duch
  } else {                                     # acute_split
    yld_ac  <- st * (SEV_SPLIT["mild_mod"]*dwmm*dumm + SEV_SPLIT["severe"]*dwsv*dusv)
    yld_chr <- (sy*(1-p14y) + so*(1-p14o)) * dwch * duch
  }
  deaths <- sum(symp_dw * cfr_j)
  yll    <- sum(symp_dw * cfr_j * le_band[age_to_band])
  c(infections = infections, symptomatic = st,
    hospitalisations = n_hosp, deaths = deaths,
    n_nonhosp = n_nonhosp, n_chronic = n_chronic,
    yld_acute = unname(yld_ac), yld_chronic = unname(yld_chr),
    yld = unname(yld_ac + yld_chr), yll = yll, daly = unname(yld_ac + yld_chr) + yll,
    doses = target_pop_elig * covv)
}
OUTCOMES <- c("infections","symptomatic","hospitalisations","deaths",
              "n_nonhosp","n_chronic","yld_acute","yld_chronic","yld","yll","daly","doses")
NNV_OUT  <- c("symptomatic","hospitalisations","deaths","daly")   # NNV reported for these

# ------------------------------------------------------------
# 5. ONE Latin-hypercube design (42 stratified uniforms) + transmission resample.
#    cols: cov, ve, deliv, delay, hosp, cfrH(9), cfrN(9), DW(3), dur(3), LE(9), rec(4)
# ------------------------------------------------------------
beta_from_ci <- function(m, lo, hi) { v <- ((hi-lo)/(2*1.96))^2; k <- m*(1-m)/v-1; c(a=m*k, b=(1-m)*k) }
ve_ab    <- beta_from_ci(0.989, 0.967, 0.998)
cov_ab   <- beta_from_ci(0.30,  0.20,  0.40)   # SAMPLED coverage of eligible 18-59
deliv_ab <- beta_from_ci(0.10,  0.09,  0.11)

lhs_col <- function(n) (sample.int(n) - runif(n)) / n
K <- 42
U <- sapply(1:K, function(j) lhs_col(N_DRAWS)); col <- 0
nextU <- function(w = 1) { idx <- (col + 1):(col + w); col <<- col + w; U[, idx, drop = FALSE] }
cov_d   <- qbeta(nextU(), cov_ab["a"],   cov_ab["b"])
ve_d    <- qbeta(nextU(), ve_ab["a"],    ve_ab["b"])
deliv_d <- qbeta(nextU(), deliv_ab["a"], deliv_ab["b"])
delay_d <- 1 + round(2 * nextU())
hosp_d  <- qbeta(nextU(), hosp_a, hosp_b)
cfrH_d  <- qbeta(nextU(9), matrix(cfr_hosp_a, N_DRAWS, 9, byrow=TRUE), matrix(cfr_hosp_b, N_DRAWS, 9, byrow=TRUE))
cfrN_d  <- qbeta(nextU(9), matrix(cfr_nonh_a, N_DRAWS, 9, byrow=TRUE), matrix(cfr_nonh_b, N_DRAWS, 9, byrow=TRUE))
dwMM_d  <- qbeta(nextU(), dp$dw_mm$a,  dp$dw_mm$b)
dwSV_d  <- qbeta(nextU(), dp$dw_sev$a, dp$dw_sev$b)
dwCH_d  <- qbeta(nextU(), dp$dw_chr$a, dp$dw_chr$b)
duMM_d  <- qlnorm(nextU(), dp$du_mm$m,  dp$du_mm$s)
duSV_d  <- qlnorm(nextU(), dp$du_sev$m, dp$du_sev$s)
duCH_d  <- qlnorm(nextU(), dp$du_chr$m, dp$du_chr$s)
le_d    <- qlnorm(nextU(9), matrix(dp$le$m, N_DRAWS, 9, byrow=TRUE), matrix(dp$le$s, N_DRAWS, 9, byrow=TRUE))
p14y_d  <- qbeta(nextU(), dp$p14_y$a, dp$p14_y$b); p14o_d <- qbeta(nextU(), dp$p14_o$a, dp$p14_o$b)
p90y_d  <- qbeta(nextU(), dp$p90_y$a, dp$p90_y$b); p90o_d <- qbeta(nextU(), dp$p90_o$a, dp$p90_o$b)
stopifnot(col == K)
t_idx <- sample.int(n_ens, N_DRAWS, replace = TRUE)   # transmission draw per LHS row
rho_i <- E$rho[t_idx]                                 # per-draw reporting rate (for reported cases)

# ------------------------------------------------------------
# 6. Monte Carlo: per draw, baseline + all scenarios reuse the SAME parameters.
# ------------------------------------------------------------
per_draw <- setNames(lapply(scen_names, function(x)
  matrix(NA_real_, N_DRAWS, length(OUTCOMES), dimnames = list(NULL, OUTCOMES))), scen_names)
wk_symp  <- setNames(lapply(scen_names, function(x) matrix(NA_real_, N_DRAWS, T_sim)), scen_names)
wk_inf   <- setNames(lapply(scen_names, function(x) matrix(NA_real_, N_DRAWS, T_sim)), scen_names)

cat(sprintf("Running %d draws x %d scenarios (%s YLD)...\n", N_DRAWS, length(scen), PHASE_MODE))
for (i in 1:N_DRAWS) {
  ti <- t_idx[i]
  beta_ext <- extend_beta(E$beta[ti, ])
  Rimm_i <- draw_immunity(E$foi[ti]); I0_i <- draw_I0(E$foi[ti], E$rho[ti], E$gamma[ti], E$prop_symp[ti])
  cfr_j  <- (hosp_d[i]*cfrH_d[i, ] + (1-hosp_d[i])*cfrN_d[i, ])[age_to_band]
  for (s in scen) {
    if (s$type == "base") { covv <- 0; vi <- 0; vb <- 0; st <- start_s3 }
    else { covv <- cov_d[i]; vb <- ve_d[i]; vi <- if (s$type == "both") ve_d[i] else 0
           st <- min(s$start + delay_d[i], T_sim) }
    out <- seirv_vaccinated(T_sim, A, N, Rimm_i, I0_i, E0, beta_ext, E$sigma[ti], E$gamma[ti],
             E$rho[ti], target_age, covv, deliv_d[i], st, vi, vb, immun_delay, prop_symp = E$prop_symp[ti])
    per_draw[[s$name]][i, ] <- outcome_one(out, covv, hosp_d[i], cfr_j, le_d[i, ],
                                 dwMM_d[i], dwSV_d[i], dwCH_d[i], duMM_d[i], duSV_d[i], duCH_d[i],
                                 p14y_d[i], p14o_d[i], p90y_d[i], p90o_d[i])
    wk_symp[[s$name]][i, ] <- colSums(out$new_symptomatic)
    wk_inf [[s$name]][i, ] <- colSums(out$new_infections)
  }
  if (i %% 100 == 0) cat("  ", i, "/", N_DRAWS, "\n")
}

# ------------------------------------------------------------
# 7. Derived per-draw quantities: averted (vs baseline) and NNV = doses / averted.
#    Paired per draw, so NNV shares the burden/DALY uncertainty exactly.
# ------------------------------------------------------------
base_pd <- per_draw[["No vaccine (baseline)"]]
averted <- setNames(lapply(vac_names, function(nm)
  base_pd[, NNV_OUT, drop=FALSE] - per_draw[[nm]][, NNV_OUT, drop=FALSE]), vac_names)
nnv <- setNames(lapply(vac_names, function(nm) {
  av <- averted[[nm]]; d <- per_draw[[nm]][, "doses"]
  m <- d / av                                  # recycle doses across the 4 outcome cols
  m[!is.finite(m) | m < 0] <- NA
  m
}), vac_names)

# ------------------------------------------------------------
# 8. Aggregated median + 95% UI tables (per-draw matrices retained separately)
# ------------------------------------------------------------
agg <- function(mat, cols = colnames(mat)) do.call(rbind, lapply(cols, function(o) {
  q <- quantile(mat[, o], c(.5, .025, .975), na.rm = TRUE)
  data.frame(outcome = o, median = q[1], lo = q[2], hi = q[3], row.names = NULL)
}))
agg_burden <- do.call(rbind, lapply(scen_names, function(nm)
  cbind(scenario = nm, agg(per_draw[[nm]]))))
agg_averted <- do.call(rbind, lapply(vac_names, function(nm)
  cbind(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm), agg(averted[[nm]]))))
agg_nnv <- do.call(rbind, lapply(vac_names, function(nm)
  cbind(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm), agg(nnv[[nm]]))))

# ------------------------------------------------------------
# 9. Console summary (headline burden/DALY + pre-outbreak NNV)
# ------------------------------------------------------------
cat("\n=== Baseline burden / DALY (median, 95% UI) ===\n")
for (o in c("infections","symptomatic","hospitalisations","deaths","daly"))
  cat(sprintf("  %-16s %s\n", o, fmtq(base_pd[, o], if (o %in% c("deaths")) 1 else 0)))

cat("\n=== NNV, pre-outbreak (median, 95% UI) ===\n")
for (nm in grep("pre-outbreak", vac_names, value = TRUE)) {
  cat(" ", nm, "\n")
  for (o in NNV_OUT) cat(sprintf("    %-16s %s\n", o, fmtq(nnv[[nm]][, o], 1)))
}

# ------------------------------------------------------------
# 10. Save BOTH per-draw and aggregated (per user request)
# ------------------------------------------------------------
saveRDS(list(
  per_draw = per_draw, averted = averted, nnv = nnv,
  wk_symp = wk_symp, wk_inf = wk_inf, rho_i = rho_i,
  agg_burden = agg_burden, agg_averted = agg_averted, agg_nnv = agg_nnv,
  scen = scen, scen_names = scen_names, vac_names = vac_names,
  timings = timings, arm_names = arm_names,
  OUTCOMES = OUTCOMES, NNV_OUT = NNV_OUT, EVAL_WIN = EVAL_WIN,
  T_sim = T_sim, T_data = T_data, caldas_obs = caldas_obs, observed_cases = observed_cases,
  N_DRAWS = N_DRAWS, PHASE_MODE = PHASE_MODE, target_pop_elig = target_pop_elig,
  cov_d = cov_d, ve_d = ve_d),
  "CHIKV_ca_engine_results.rds")
cat("\nSaved CHIKV_ca_engine_results.rds (per-draw + aggregated; severity-phase counts included).\n")
