# ============================================================
# CHIKV_ca_engine.R -- Caldas Novas CHIKV Monte Carlo engine.
#
# One Latin hypercube propagation from transmission through to health loss. Each of
# N_DRAWS draws combines a re-fitted transmission trajectory resampled from
# CHIKV_ca_lhs_ensemble.rds with vaccine efficacy, delivery, delay, sampled coverage,
# severity, case fatality and every DALY input. That single draw is then run through
# the SEIRV for the no-vaccine baseline and every vaccination scenario, so all
# outcomes are mutually consistent draw for draw:
#     infections, symptomatic, hospitalisations, deaths
#     severity-phase counts (hospitalised / non-hospitalised / sub-acute / chronic)
#     YLD by phase, YLL, DALYs
#     doses delivered, and NNV = doses / burden averted at that draw's coverage
#
# The horizon is the 52-week observed window, 2025-W24 -> 2026-W22. Nothing is
# projected past the surveillance data, so no transmission rate is ever assumed.
# Every scenario, baseline included, is scored over this same fixed window; anything
# a scenario pushes beyond week 52 is out of scope for all scenarios alike.
#
# Saves per-draw matrices and aggregated median / 95% UI to
# CHIKV_ca_engine_results.rds. Severity-phase COUNTS are saved so the cost layer can
# multiply unit costs on afterwards without re-running the SEIR.
#
# Run order:  CHIKV_ca_lhs.R  (slow, once)  ->  this file  ->  outputs / costs
# ============================================================
library(dplyr); library(tidyr)
source("ca_common.R")   # fmtq, qs, burden, load_burden_params, load_caldas_age_cases,
                        # compute_age_weight, seirv_vaccinated, load_daly_params

# ---- knobs --------------------------------------------------
N_DRAWS    <- 1000
# YLD accrues by RECOVERY FUNNEL: each patient is counted exactly once, in the phase in
# which their illness resolved, using that group's total time-to-recovery as the duration.
# Phases are NOT nested -- dur_chronic (0.53 yr) is O'Driscoll's median time to arthralgia
# resolution measured FROM INFECTION (6.39 mo, 95% CI 5.48-7.66), so it already spans the
# acute and sub-acute period; adding those on top would double count.
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
T_sim <- T_data                 # 52 weeks; no projection beyond the observed data
# Immunity rebuilt exactly as the fit built it: truncated catalytic model, 12-yr cap.
exposure_age  <- E$exposure_age
draw_immunity <- function(foi) 1 - exp(-foi * exposure_age)
draw_I0 <- function(foi, rho, gamma, ps) {
  Rimm <- draw_immunity(foi); sf <- N*(1-Rimm)/sum(N*(1-Rimm))
  round(((week_1_cases/rho/ps)/gamma) * sf)
}

# ------------------------------------------------------------
# 2. Vaccine programme, scenarios and evaluation window
# ------------------------------------------------------------
target_age <- rep(0, A); target_age[c(4,5,6,7,8)] <- 1     # eligible = 18-59
target_pop_elig <- sum(N[target_age == 1])
immun_delay <- 2
idx_of <- function(yr, wk) caldas_obs$week_index[caldas_obs$Year == yr & caldas_obs$week == wk]
start_s1 <- idx_of(2026, 16)     # IXCHIQ real rollout
start_s2 <- idx_of(2026, 1)      # start of 2026
start_s3 <- idx_of(2025, 40)     # pre-outbreak
start_s0 <- 1                    # earliest: rollout completes well before the outbreak
timings <- list("earliest (2025-W24)" = start_s0, "actual rollout" = start_s1,
                "start of 2026" = start_s2, "pre-outbreak" = start_s3)
arm_names <- c("Disease-blocking", "Disease + infection blocking")

# Burden accrues over the full observed window, 2025-W24 -> 2026-W22.
EVAL_WIN <- seq_len(T_sim)
stopifnot(length(EVAL_WIN) == 52)
cat("Evaluation window: weeks 1-52 (2025-W24 -> 2026-W22), matching the data.\n")

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
out_base_pt <- seirv_vaccinated(T_sim, A, N, Rimm_base, I0_base, E0, E$base_beta,
                 E$base_sigma, E$base_gamma, E$base_rho, target_age, 0, 0.10, start_s3,
                 0, 0, immun_delay, prop_symp = E$base_prop_symp)
age_weight <- compute_age_weight(rowSums(out_base_pt$new_infections), obs_band_prop, age_to_band)
young_idx <- which(age_to_band <= 4); old_idx <- which(age_to_band >= 5)

# ------------------------------------------------------------
# 3b. Burden calculation AUDIT (point estimate, no-vaccine baseline, within EVAL_WIN).
#     A hand-checkable trace of how each number leads to the next. POINT estimate
#     (mean params) so every step ties out EXACTLY; the MC medians in the other tabs
#     are close but NOT identical (median of a product != product of medians).
#     Chain: infections -> symptomatic (x prop_symp) -> hospitalisations (x hosp_rate,
#     SINGLE rate) ; deaths are age-specific: sum over age of symptomatic_age x CFR_age.
# ------------------------------------------------------------
aud_inf_age  <- rowSums(out_base_pt$new_infections[,  EVAL_WIN, drop = FALSE])   # true infections by age
aud_symp_age <- rowSums(out_base_pt$new_symptomatic[, EVAL_WIN, drop = FALSE])   # = prop_symp * inf_age
aud_symp_w   <- aud_symp_age * age_weight
aud_symp_dw  <- if (sum(aud_symp_w) > 0) aud_symp_w * (sum(aud_symp_age)/sum(aud_symp_w)) else aud_symp_age
aud_inf  <- sum(aud_inf_age); aud_symp <- sum(aud_symp_age)
aud_hosp <- aud_symp * hosp_rate                    # SINGLE rate, all ages
aud_dth_age <- aud_symp_dw * cfr_vec; aud_dth <- sum(aud_dth_age)
aud_qb <- function(a, b) c(a/(a+b), qbeta(.025, a, b), qbeta(.975, a, b))
aud_ps <- aud_qb(ps_a, ps_b); aud_hr <- aud_qb(hosp_a, hosp_b)
burden_audit <- data.frame(
  step     = 1:5,
  quantity = c("True infections","True symptomatic","Reported cases","Hospitalisations","Deaths"),
  value    = c(round(aud_inf), round(aud_symp), round(E$base_rho*aud_symp), round(aud_hosp,1), round(aud_dth,2)),
  formula  = c("SEIR: sum(new_infections) in the 52-wk window",
               "true infections x prop_symp",
               "true symptomatic x rho   (= rho x prop_symp x infections)",
               "true symptomatic x hosp_rate   [SINGLE rate, ALL ages]",
               "sum over age of (symptomatic_age x CFR_age)   [age-specific]"),
  parameter = c("-",
                sprintf("prop_symp = %.3f (95%% UI %.3f-%.3f)", aud_ps[1], aud_ps[2], aud_ps[3]),
                sprintf("rho = %.3f (base; per-draw Beta(20,60) in the MC)", E$base_rho),
                sprintf("hosp_rate = %.4f (95%% UI %.4f-%.4f)", aud_hr[1], aud_hr[2], aud_hr[3]),
                "CFR by age -> see burden_audit_by_age tab"),
  stringsAsFactors = FALSE)
aud_band <- c("[0,10)","[10,20)","[20,30)","[30,40)","[40,50)","[50,60)","[60,70)","[70,80)","[80,90)")
burden_audit_by_age <- rbind(
  data.frame(age_group = as.character(age_df$age_group), cfr_band = aud_band[age_to_band],
             symptomatic_reweighted = round(aud_symp_dw, 1), CFR_death_per_symp = signif(cfr_vec, 3),
             deaths = round(aud_dth_age, 3), stringsAsFactors = FALSE),
  data.frame(age_group = "TOTAL", cfr_band = "-", symptomatic_reweighted = round(aud_symp, 1),
             CFR_death_per_symp = NA_real_, deaths = round(aud_dth, 3)))
cat(sprintf("Burden audit (point est.): infections %.0f -> symptomatic %.0f -> hosp %.1f -> deaths %.2f\n",
            aud_inf, aud_symp, aud_hosp, aud_dth))

dp <- load_daly_params()

# ------------------------------------------------------------
# 4. Outcome extractor: one SEIR run -> all outcomes within the 52-wk window.
#    Severity phases (hospitalised / non-hosp / chronic) are saved as COUNTS so a
#    cost layer can multiply unit costs onto them later without re-running.
# ------------------------------------------------------------
outcome_one <- function(out, covv, hosp_j, cfr_j, le_band,
                        dwmm, dwsv, dwch, dumm, dusv, dusb, duch, acy, aco, sby, sbo, chy, cho) {
  infections <- sum(out$new_infections[, EVAL_WIN, drop = FALSE])
  symp_age   <- rowSums(out$new_symptomatic[, EVAL_WIN, drop = FALSE])
  symp_w  <- symp_age * age_weight
  symp_dw <- if (sum(symp_w) > 0) symp_w * (sum(symp_age)/sum(symp_w)) else symp_age
  st <- sum(symp_dw); sy <- sum(symp_dw[young_idx]); so <- sum(symp_dw[old_idx])

  n_hosp    <- st * hosp_j                     # hospitalised (= severe acute in hosp_severity)
  n_nonhosp <- st * (1 - hosp_j)               # non-hospitalised (mild/mod acute)
  # Recovery funnel: resolved <=14d (acute) | 14d-3m (sub-acute) | >3m (chronic).
  # The three shares are MARGINAL and are renormalised to sum to 1 within each age class,
  # so every symptomatic case is counted in exactly one phase.
  sy_t <- acy + sby + chy; so_t <- aco + sbo + cho
  n_acute    <- sy*acy/sy_t + so*aco/so_t
  n_subacute <- sy*sby/sy_t + so*sbo/so_t
  n_chronic  <- sy*chy/sy_t + so*cho/so_t
  f_acute    <- n_acute / st

  yld_ac  <- (n_nonhosp*dwmm*dumm + n_hosp*dwsv*dusv) * f_acute  # resolved <=14d
  yld_sub <- n_subacute * dwch * dusb   # resolved 14d-3m (CHRONIC disability weight)
  yld_chr <- n_chronic  * dwch * duch   # resolved >3m

  deaths <- sum(symp_dw * cfr_j)
  yll    <- sum(symp_dw * cfr_j * le_band[age_to_band])
  yld_tot <- unname(yld_ac + yld_sub + yld_chr)
  c(infections = infections, symptomatic = st,
    hospitalisations = n_hosp, deaths = deaths,
    n_nonhosp = n_nonhosp, n_subacute = n_subacute, n_chronic = n_chronic,
    yld_acute = unname(yld_ac), yld_subacute = unname(yld_sub),
    yld_chronic = unname(yld_chr),
    yld = yld_tot, yll = yll, daly = yld_tot + yll,
    doses = target_pop_elig * covv)
}
OUTCOMES <- c("infections","symptomatic","hospitalisations","deaths",
              "n_nonhosp","n_subacute","n_chronic",
              "yld_acute","yld_subacute","yld_chronic","yld","yll","daly","doses")
# NNV reported for these. "infections" is meaningful only for infection-blocking arms;
# disease-blocking leaves infections unchanged, so its infection NNV is undefined (NA).
NNV_OUT  <- c("infections","symptomatic","hospitalisations","deaths","daly")

# ------------------------------------------------------------
# 5. ONE Latin-hypercube design (42 stratified uniforms) + transmission resample.
#    cols: cov, ve, deliv, delay, hosp, cfrH(9), cfrN(9), DW(3), dur(3), LE(9), rec(4)
# ------------------------------------------------------------
beta_from_ci <- function(m, lo, hi) { v <- ((hi-lo)/(2*1.96))^2; k <- m*(1-m)/v-1; c(a=m*k, b=(1-m)*k) }
# VE proxy = day-29 seroprotection from the VLA1553/Ixchiq trial, 263/266 seroprotected
# (Total column). Exact Beta posterior from the raw counts (uniform prior), which
# reproduces the published 96.7-99.8% Clopper-Pearson CI and carries the true trial n.
ve_ab    <- c(a = 263 + 1, b = 266 - 263 + 1)   # Beta(264, 4)
cov_ab   <- beta_from_ci(0.30,  0.20,  0.40)   # SAMPLED coverage of eligible 18-59
deliv_ab <- beta_from_ci(0.10,  0.09,  0.11)

lhs_col <- function(n) (sample.int(n) - runif(n)) / n
K <- 49
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
duSB_d  <- qlnorm(nextU(), dp$du_sub$m, dp$du_sub$s)
duCH_d  <- qlnorm(nextU(), dp$du_chr$m, dp$du_chr$s)
le_d    <- qlnorm(nextU(9), matrix(dp$le$m, N_DRAWS, 9, byrow=TRUE), matrix(dp$le$s, N_DRAWS, 9, byrow=TRUE))
# Marginal recovery shares: acute (<=14d), sub-acute (14d-3m), chronic (>3m = 6m+12m+30m)
acy_d <- qbeta(nextU(), dp$p14_y$a, dp$p14_y$b); aco_d <- qbeta(nextU(), dp$p14_o$a, dp$p14_o$b)
sby_d <- qbeta(nextU(), dp$p90_y$a, dp$p90_y$b); sbo_d <- qbeta(nextU(), dp$p90_o$a, dp$p90_o$b)
chy_d <- qbeta(nextU(), dp$p6_y$a,  dp$p6_y$b)  + qbeta(nextU(), dp$p12_y$a, dp$p12_y$b) +
         qbeta(nextU(), dp$p30_y$a, dp$p30_y$b)
cho_d <- qbeta(nextU(), dp$p6_o$a,  dp$p6_o$b)  + qbeta(nextU(), dp$p12_o$a, dp$p12_o$b) +
         qbeta(nextU(), dp$p30_o$a, dp$p30_o$b)
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
# Dose accounting per scenario: doses ADMINISTERED to eligible 18-59, and the subset
# reaching SUSCEPTIBLES. Wastage = 1 - on-target/administered (doses to already-immune
# or already-infected eligible people; higher for later-timing rollouts).
doses_deliv    <- setNames(lapply(scen_names, function(x) numeric(N_DRAWS)), scen_names)
doses_ontarget <- setNames(lapply(scen_names, function(x) numeric(N_DRAWS)), scen_names)

cat(sprintf("Running %d draws x %d scenarios...\n", N_DRAWS, length(scen)))
for (i in 1:N_DRAWS) {
  ti <- t_idx[i]
  Rimm_i <- draw_immunity(E$foi[ti]); I0_i <- draw_I0(E$foi[ti], E$rho[ti], E$gamma[ti], E$prop_symp[ti])
  cfr_j  <- (hosp_d[i]*cfrH_d[i, ] + (1-hosp_d[i])*cfrN_d[i, ])[age_to_band]
  for (s in scen) {
    if (s$type == "base") { covv <- 0; vi <- 0; vb <- 0; st <- start_s3 }
    else { covv <- cov_d[i]; vb <- ve_d[i]; vi <- if (s$type == "both") ve_d[i] else 0
           st <- min(s$start + delay_d[i], T_sim) }
    out <- seirv_vaccinated(T_sim, A, N, Rimm_i, I0_i, E0, E$beta[ti, ], E$sigma[ti], E$gamma[ti],
             E$rho[ti], target_age, covv, deliv_d[i], st, vi, vb, immun_delay, prop_symp = E$prop_symp[ti])
    per_draw[[s$name]][i, ] <- outcome_one(out, covv, hosp_d[i], cfr_j, le_d[i, ],
                                 dwMM_d[i], dwSV_d[i], dwCH_d[i], duMM_d[i], duSV_d[i],
                                 duSB_d[i], duCH_d[i],
                                 acy_d[i], aco_d[i], sby_d[i], sbo_d[i], chy_d[i], cho_d[i])
    wk_symp[[s$name]][i, ] <- colSums(out$new_symptomatic)
    wk_inf [[s$name]][i, ] <- colSums(out$new_infections)
    doses_deliv   [[s$name]][i] <- sum(out$total_used_age)          # administered
    doses_ontarget[[s$name]][i] <- sum(out$V_covered[, T_sim])      # reached susceptibles
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
  N_DRAWS = N_DRAWS, target_pop_elig = target_pop_elig,
  doses_deliv = doses_deliv, doses_ontarget = doses_ontarget,
  cov_d = cov_d, ve_d = ve_d,
  burden_audit = burden_audit, burden_audit_by_age = burden_audit_by_age),
  "CHIKV_ca_engine_results.rds")
cat("\nSaved CHIKV_ca_engine_results.rds (per-draw + aggregated; severity-phase counts included).\n")
