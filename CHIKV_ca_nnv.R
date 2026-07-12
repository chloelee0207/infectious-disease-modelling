# ============================================================
# CHIKV_ca_nnv.R -- Caldas Novas CHIKV Number Needed to Vaccinate (NNV),
#                   STANDALONE on the LHS ensemble.
# ------------------------------------------------------------
# Estimates the NNV to avert one (a) symptomatic case, (b) hospitalisation,
# (c) death, and (d) DALY for the PRE-OUTBREAK rollout (2025-W40), at a CENTRAL
# coverage of 30% with 10% / 90% extremes. All burden is accrued over the 52-epi-
# week window 2025-W40 -> 2026-W38 (matches CHIKV_ca_vacc_outputs_52.xlsx).
#
# NNV definition (population-level, captures herd effects -- the construction used
# by Hyolim Kang's CHIK_VIM reference this mirrors):
#     NNV(outcome) = (people vaccinated) / (outcome averted vs no-vaccine baseline)
# where "people vaccinated" = target population x coverage (eligible 18-59 pop x
# coverage; = Hyolim's n*coverage) and "outcome averted" = 52-wk baseline burden -
# 52-wk scenario burden. This generalises the individual-level NNV = 1/ARR to a
# transmission model by dividing the whole-population avertible burden into the
# doses that produced it, so NNV < 1 is possible (one dose averts >1 case via herd).
#
# The SEIR + scenario setup, burden extractor and DALY structure MIRROR
# CHIKV_ca_vacc.R and CHIKV_ca_daly.R (all standalone by design). Keep the SEIR /
# DALY blocks in sync with those engines if they change.
#
# Uncertainty: with RUN_MC = TRUE a 1000-draw Latin-hypercube MC gives 95% UI error
# bars via Hyolim's flip (median NNV = tot_vacc / averted_median; UI = tot_vacc /
# averted_{97.5%, 2.5%}). RUN_MC = FALSE falls back to a point estimate (median
# inputs); the figure/table adapt automatically.
#
# Run order:  source("ca_common.R"); source("CHIKV_ca_lhs.R")  # ensemble (slow, once)
#             source("CHIKV_ca_nnv.R")                          # this file
# ============================================================
library(dplyr); library(tidyr); library(ggplot2); library(writexl)
source("ca_common.R")   # fmtq, qs, burden, load_burden_params, load_caldas_age_cases, compute_age_weight

# ---- knobs --------------------------------------------------
COV_LEVELS  <- c(0.10, 0.30, 0.90)   # coverage sweep; 0.30 = central, 10%/90% extremes
COV_CENTRAL <- 0.30
RUN_MC      <- TRUE                   # TRUE -> 1000-draw LHS for 95% UI error bars (slow)
N_DRAWS     <- 1000
PHASE_MODE  <- "hosp_severity"        # DALY YLD axis (matches CHIKV_ca_daly.R default)
SEV_SPLIT   <- c(mild_mod = 0.47, severe = 0.53)   # used only by PHASE_MODE == "acute_split"
set.seed(2029)

# ------------------------------------------------------------
# 0. Load the LHS ensemble (same object the vacc / daly engines consume)
# ------------------------------------------------------------
if (!file.exists("CHIKV_ca_lhs_ensemble.rds"))
  stop("CHIKV_ca_lhs_ensemble.rds not found -- run CHIKV_ca_lhs.R first.")
E <- readRDS("CHIKV_ca_lhs_ensemble.rds")
N <- E$N; A <- E$A; age_df <- E$age_df; age_mid <- E$age_mid
observed_cases <- E$observed_cases; caldas_obs <- E$caldas_obs
week_1_cases <- E$week_1_cases; T_data <- E$T_weeks
n_ens <- nrow(E$beta)
E0 <- rep(0, A)
cat(sprintf("Loaded ensemble: %d feasible transmission draws.\n", n_ens))

# ------------------------------------------------------------
# 1. SEIRV simulator (age-structured, weekly) -- mirrors CHIKV_ca_vacc.R / _daly.R.
#    Returns total_used_age (doses delivered per age) so we can count who was
#    vaccinated -- the NNV denominator's numerator.
# ------------------------------------------------------------
seirv_vaccinated <- function(
    T_weeks, A, N, R_init_prop, I0, E0, base_beta, sigma, gamma, rho,
    target_age, total_coverage, weekly_delivery_speed, delay,
    VE_inf = 0.989, VE_block = 0, immun_delay = 2, prop_symp = 0.5242478, sub_steps = 7) {
  pmax0 <- function(x) pmax(0, x); N_total <- sum(N); dt <- 1/sub_steps
  S <- E <- I <- R <- V <- matrix(0, A, T_weeks)
  V_covered <- vacc_delayed <- coverage_frac <- matrix(0, A, T_weeks)
  new_infections <- new_symptomatic <- matrix(0, A, T_weeks)
  target_idx <- which(target_age == 1); target_pop <- sum(N[target_idx])
  total_supply <- target_pop * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed
  total_avail_age <- rep(0, A); total_avail_age[target_idx] <- total_supply * (N[target_idx]/target_pop)
  total_used_age <- rep(0, A); unvaccinated <- N
  S_now <- pmax0(N - I0 - E0 - R_init_prop*N); E_now <- E0; I_now <- I0
  R_now <- R_init_prop*N; V_now <- rep(0, A)
  for (t in 1:T_weeks) {
    prev_V_covered <- if (t > 1) V_covered[, t-1] else rep(0, A)
    if (t - immun_delay >= 1) {
      effective_dose <- vacc_delayed[, t-immun_delay]
      immunized <- round(VE_inf * effective_dose)
      V_covered[, t] <- prev_V_covered + effective_dose
    } else { immunized <- rep(0, A); V_covered[, t] <- prev_V_covered }
    S_now <- pmax0(S_now - immunized); V_now <- V_now + immunized
    coverage_frac[, t] <- V_covered[, t] / N
    if (t >= delay && target_pop > 0) {
      rem <- weekly_dose_total
      for (a in target_idx) {
        alloc <- min(ceiling(weekly_dose_total*(N[a]/target_pop)), rem,
                     unvaccinated[a], total_avail_age[a]-total_used_age[a])
        if (alloc > 0) {
          prop_S <- if (N[a] > 0) S_now[a]/N[a] else 0
          vacc_to_S <- round(alloc*prop_S)
          vacc_delayed[a, t] <- vacc_to_S
          total_used_age[a] <- total_used_age[a] + alloc
          unvaccinated[a] <- unvaccinated[a] - alloc; rem <- rem - alloc
        }
      }
    }
    new_I_week <- rep(0, A); beta_t <- base_beta[t]
    for (k in 1:sub_steps) {
      foi <- beta_t * sum(I_now)/N_total
      new_E <- foi*S_now*dt; new_I <- sigma*E_now*dt; new_R <- gamma*I_now*dt
      S_now <- pmax0(S_now-new_E); E_now <- pmax0(E_now+new_E-new_I)
      I_now <- pmax0(I_now+new_I-new_R); R_now <- pmax0(R_now+new_R)
      new_I_week <- new_I_week + new_I
    }
    S[,t]<-S_now; E[,t]<-E_now; I[,t]<-I_now; R[,t]<-R_now; V[,t]<-V_now
    new_infections[,t] <- new_I_week
    new_symptomatic[,t] <- prop_symp*new_I_week*(1 - VE_block*coverage_frac[,t])
  }
  list(new_infections=new_infections, new_symptomatic=new_symptomatic,
       new_reported=rho*new_symptomatic, V_covered=V_covered, total_used_age=total_used_age)
}

# ------------------------------------------------------------
# 2. Horizon + per-draw immunity/seed (mirrors CHIKV_ca_vacc.R / _daly.R)
# ------------------------------------------------------------
EXTEND <- 26; T_sim <- T_data + EXTEND
extend_beta <- function(b) c(b, rep(b[length(b)], EXTEND))
draw_immunity <- function(foi) 1 - exp(-foi * age_mid)
draw_I0 <- function(foi, rho, gamma, ps) {
  Rimm <- draw_immunity(foi); sf <- N*(1-Rimm)/sum(N*(1-Rimm))
  round(((week_1_cases/rho/ps)/gamma) * sf)
}

# ------------------------------------------------------------
# 3. Vaccine programme + scenarios (52-week index)
# ------------------------------------------------------------
target_age <- rep(0, A); target_age[c(4,5,6,7,8)] <- 1     # eligible = 18-59
target_pop_elig <- sum(N[target_age == 1])                 # eligible population (NNV numerator base)
immun_delay <- 2
idx_of <- function(yr, wk) caldas_obs$week_index[caldas_obs$Year == yr & caldas_obs$week == wk]
start_s1 <- idx_of(2026, 16)     # IXCHIQ real rollout
start_s2 <- idx_of(2026, 1)      # start of 2026
start_s3 <- idx_of(2025, 40)     # pre-outbreak
timings <- list("pre-outbreak" = start_s3)     # only the pre-outbreak rollout is of interest

# 52-epidemiological-week evaluation window: 2025-W40 (pre-outbreak vaccine
# implementation) -> 2026-W38. Matches CHIKV_ca_vacc_outputs_52.xlsx. All burden
# (and hence NNV) is accrued ONLY within this window; the pre-implementation weeks
# (1..16) are vaccine-independent and cancel in "averted", while cases the vaccine
# pushes past week 52 are correctly credited as averted-within-the-window.
EVAL_WIN <- start_s3:(start_s3 + 51)
stopifnot(length(EVAL_WIN) == 52, max(EVAL_WIN) <= T_sim)
cat(sprintf("52-week evaluation window: index %d-%d (2025-W40 -> 2026-W38).\n",
            min(EVAL_WIN), max(EVAL_WIN)))
arm_names <- c("Disease-blocking", "Disease + infection blocking")
scen <- list(list(name="No vaccine (baseline)", timing="No vaccine", arm="No vaccine",
                  start=NA_integer_, type="base"))
for (tn in names(timings)) for (an in arm_names)
  scen[[length(scen)+1]] <- list(name=paste0(tn," | ",an), timing=tn, arm=an,
                                 start=timings[[tn]], type=if(an=="Disease-blocking")"disb" else "both")
scen_names <- vapply(scen, function(s) s$name, character(1))
vac_names  <- setdiff(scen_names, "No vaccine (baseline)")

# ------------------------------------------------------------
# 4. Severity params (CFR/hosp) + observed-age (w_a) death correction
# ------------------------------------------------------------
invisible(list2env(load_burden_params(A), globalenv()))     # ps_*, hosp_*, cfr_*, age_to_band, cfr_vec
obs_band_prop <- load_caldas_age_cases()$obs_band_prop
I0_base   <- draw_I0(E$base_foi, E$base_rho, E$base_gamma, E$base_prop_symp)
Rimm_base <- draw_immunity(E$base_foi)
out_base_pt <- seirv_vaccinated(T_sim, A, N, Rimm_base, I0_base, E0, extend_beta(E$base_beta),
                 E$base_sigma, E$base_gamma, E$base_rho, target_age, 0, 0.10, start_s3,
                 0, 0, immun_delay, prop_symp = E$base_prop_symp)
age_weight <- compute_age_weight(rowSums(out_base_pt$new_infections), obs_band_prop, age_to_band)

young_idx <- which(age_to_band <= 4)     # groups mapping to [0,40)
old_idx   <- which(age_to_band >= 5)     # groups mapping to [40,90)

# ------------------------------------------------------------
# 5. DALY parameters from disease_progression.xlsx (same loader as CHIKV_ca_daly.R)
# ------------------------------------------------------------
load_daly_params <- function(dp_path = "disease_progression.xlsx",
                             dp_sheet = "disease_progression") {
  dp <- read_excel(dp_path, sheet = dp_sheet)
  names(dp)[1:10] <- c("parameter","group","median","ui_lo","ui_hi","dist","p1","v1","p2","v2")
  row <- function(param, grp = NULL) {
    d <- dp[dp$parameter == param, ]
    if (!is.null(grp)) d <- d[grepl(grp, d$group), ]
    d
  }
  beta_ab  <- function(param, grp = NULL) { d <- row(param, grp); list(a = d$v1, b = d$v2) }
  lnorm_ms <- function(param, grp = NULL) { d <- row(param, grp); list(m = d$v1, s = d$v2, med = d$median) }
  dw_mm  <- beta_ab("Disability weight for mild and moderate chikungunya")
  dw_sev <- beta_ab("Disability weight for severe chikungunya")
  dw_chr <- beta_ab("Disability weight for chronic chikungunya")
  du_mm  <- lnorm_ms("Duration of illness for mild and moderate chikungunya (years)")
  du_sev <- lnorm_ms("Duration of illness for severe chikungunya (years)")
  du_chr <- lnorm_ms("Duration of illness for chronic chikungunya (years)")
  le <- row("Remaining life-years")
  lo <- as.numeric(sub(".*\\[\\s*([0-9]+).*", "\\1", le$group)); le <- le[order(lo), ]
  le_ms <- list(m = le$v1, s = le$v2, med = le$median)
  stopifnot(length(le_ms$m) == 9)
  p14_y <- beta_ab("Probability of recovery within 14 days after onset of symptoms", "< 40")
  p14_o <- beta_ab("Probability of recovery within 14 days after onset of symptoms", "> 40")
  p90_y <- beta_ab("Probability of recovery within 90 days after acute period", "< 40")
  p90_o <- beta_ab("Probability of recovery within 90 days after acute period", "> 40")
  list(dw_mm=dw_mm, dw_sev=dw_sev, dw_chr=dw_chr,
       du_mm=du_mm, du_sev=du_sev, du_chr=du_chr, le=le_ms,
       p14_y=p14_y, p14_o=p14_o, p90_y=p90_y, p90_o=p90_o)
}
dp <- load_daly_params()

# YLD/YLL/DALY for one scenario run, one draw (identical structure to CHIKV_ca_daly.R)
daly_one <- function(symp_age, cfr_vec_j, hosp_j, le_band, dwmm, dwsv, dwch, dumm, dusv, duch,
                     p14y, p14o, p90y, p90o) {
  symp_w  <- symp_age * age_weight
  symp_dw <- if (sum(symp_w) > 0) symp_w * (sum(symp_age)/sum(symp_w)) else symp_age
  sy <- sum(symp_dw[young_idx]); so <- sum(symp_dw[old_idx]); st <- sy + so
  if (PHASE_MODE == "hosp_severity") {
    yld_ac  <- st * ((1-hosp_j)*dwmm*dumm + hosp_j*dwsv*dusv)
    yld_chr <- (sy*(1-p14y)*(1-p90y) + so*(1-p14o)*(1-p90o)) * dwch * duch
  } else if (PHASE_MODE == "three_phase") {
    yld_ac  <- st * dwmm * dumm +
               (sy*(1-p14y) + so*(1-p14o)) * dwsv * dusv
    yld_chr <- (sy*(1-p14y)*(1-p90y) + so*(1-p14o)*(1-p90o)) * dwch * duch
  } else {                       # acute_split
    yld_ac  <- st * (SEV_SPLIT["mild_mod"]*dwmm*dumm + SEV_SPLIT["severe"]*dwsv*dusv)
    yld_chr <- (sy*(1-p14y) + so*(1-p14o)) * dwch * duch
  }
  yld <- unname(yld_ac + yld_chr)
  deaths_age <- symp_dw * cfr_vec_j
  yll <- sum(deaths_age * le_band[age_to_band])
  c(daly = yld + yll, deaths = sum(deaths_age))
}

# ------------------------------------------------------------
# 6. POINT-ESTIMATE inputs (median of every distribution) -- default figure/table.
# ------------------------------------------------------------
ve_pt    <- 0.989                                    # vaccine efficacy (median, Hyolim Table 1)
deliv_pt <- 0.10; delay_pt <- 2                      # weekly delivery, deployment delay (medians)
# median DALY inputs
dwMM_pt <- qbeta(.5, dp$dw_mm$a,  dp$dw_mm$b)
dwSV_pt <- qbeta(.5, dp$dw_sev$a, dp$dw_sev$b)
dwCH_pt <- qbeta(.5, dp$dw_chr$a, dp$dw_chr$b)
duMM_pt <- dp$du_mm$med; duSV_pt <- dp$du_sev$med; duCH_pt <- dp$du_chr$med
le_pt   <- dp$le$med                                 # length 9
p14y_pt <- qbeta(.5, dp$p14_y$a, dp$p14_y$b); p14o_pt <- qbeta(.5, dp$p14_o$a, dp$p14_o$b)
p90y_pt <- qbeta(.5, dp$p90_y$a, dp$p90_y$b); p90o_pt <- qbeta(.5, dp$p90_o$a, dp$p90_o$b)

outcomes <- c("symptomatic","hospitalisations","deaths","daly")

# run one scenario at a given coverage; return the 4 burden outcomes + people vaccinated
run_scen <- function(s, cov, beta_ext, sigma, gamma, rho, ps, Rimm, I0i,
                     ve, deliv, delay, hosp_j, cfr_j,
                     dwmm, dwsv, dwch, dumm, dusv, duch, le_band, p14y, p14o, p90y, p90o) {
  if (s$type == "base") { covv<-0; vi<-0; vb<-0; st<-start_s3 }
  else { covv<-cov; vb<-ve; vi<-if (s$type=="both") ve else 0
         st<-min(s$start + delay, T_sim) }
  out <- seirv_vaccinated(T_sim, A, N, Rimm, I0i, E0, beta_ext, sigma, gamma, rho,
           target_age, covv, deliv, st, vi, vb, immun_delay, prop_symp = ps)
  symp_age <- rowSums(out$new_symptomatic[, EVAL_WIN, drop = FALSE])   # 52-wk window only
  dd <- daly_one(symp_age, cfr_j, hosp_j, le_band, dwmm, dwsv, dwch, dumm, dusv, duch,
                 p14y, p14o, p90y, p90o)
  # people vaccinated = target population x coverage (Hyolim's n*coverage; ~equal to
  # the engine's delivered-dose count sum(out$total_used_age) but exact & deterministic)
  c(symptomatic = sum(symp_age), hospitalisations = sum(symp_age) * hosp_j,
    deaths = unname(dd["deaths"]), daly = unname(dd["daly"]),
    vaccinated = target_pop_elig * covv)
}

# ------------------------------------------------------------
# 7. Point-estimate NNV: baseline once, then each coverage x scenario.
#    NNV(outcome) = vaccinated / (baseline_outcome - scenario_outcome).
# ------------------------------------------------------------
beta_base_ext <- extend_beta(E$base_beta)
pt_args <- list(beta_base_ext, E$base_sigma, E$base_gamma, E$base_rho, E$base_prop_symp,
                Rimm_base, I0_base, ve_pt, deliv_pt, delay_pt, hosp_rate, cfr_vec,
                dwMM_pt, dwSV_pt, dwCH_pt, duMM_pt, duSV_pt, duCH_pt, le_pt,
                p14y_pt, p14o_pt, p90y_pt, p90o_pt)
base_pt <- do.call(run_scen, c(list(scen[[1]], 0), pt_args))

nnv_pt <- do.call(rbind, lapply(COV_LEVELS, function(cov) {
  do.call(rbind, lapply(vac_names, function(nm) {
    s <- scen[[which(scen_names == nm)]]
    b <- do.call(run_scen, c(list(s, cov), pt_args))
    av <- base_pt[outcomes] - b[outcomes]          # averted (baseline - scenario)
    data.frame(coverage = cov, timing = s$timing, arm = s$arm, outcome = outcomes,
               averted = as.numeric(av),
               vaccinated = unname(b["vaccinated"]),
               nnv = unname(b["vaccinated"]) / as.numeric(av),
               row.names = NULL)
  }))
}))
nnv_pt$nnv[!is.finite(nnv_pt$nnv) | nnv_pt$nnv < 0] <- NA   # guard: no/negative aversion

# ------------------------------------------------------------
# 8. (Optional) 1000-draw LHS Monte Carlo for 95% UI error bars.
#    Mirrors CHIKV_ca_daly.R's design; adds people-vaccinated per draw and forms
#    the NNV ratio per draw so the UI captures the full joint uncertainty.
# ------------------------------------------------------------
nnv_mc <- NULL
if (RUN_MC) {
  cat(sprintf("Running %d-draw LHS for NNV 95%% UI (3 coverages x 6 scenarios)...\n", N_DRAWS))
  beta_from_ci <- function(m, lo, hi) { v <- ((hi-lo)/(2*1.96))^2; k <- m*(1-m)/v-1; c(a=m*k, b=(1-m)*k) }
  ve_ab    <- beta_from_ci(0.989, 0.967, 0.998)
  deliv_ab <- beta_from_ci(0.10,  0.09,  0.11)

  lhs_col <- function(n) (sample.int(n) - runif(n)) / n
  K <- 41   # 4 vaccine (ve,deliv,delay,hosp) + 9 cfrH + 9 cfrN + 3 DW + 3 dur + 9 LE + 4 recovery
  U <- sapply(1:K, function(j) lhs_col(N_DRAWS)); col <- 0
  nextU <- function(w = 1) { idx <- (col + 1):(col + w); col <<- col + w; U[, idx, drop = FALSE] }
  ve_d    <- qbeta(nextU(), ve_ab["a"], ve_ab["b"])
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
  t_idx   <- sample.int(n_ens, N_DRAWS, replace = TRUE)

  # Per-draw stores: averted burden [coverage, scenario, outcome] and people
  # vaccinated [coverage, scenario]. NNV is formed AFTER the MC from these, using
  # Hyolim Kang's construction (nnv = tot_vacc / averted; UI = tot_vacc / averted-CI
  # flipped), not a per-draw ratio -- see summary block below.
  key  <- function(cov, nm, o) paste(cov, nm, o, sep = "|")
  keyv <- function(cov, nm)    paste(cov, nm, sep = "|")
  av_store  <- new.env(parent = emptyenv())
  vac_store <- new.env(parent = emptyenv())
  for (cov in COV_LEVELS) for (nm in vac_names) {
    assign(keyv(cov, nm), rep(NA_real_, N_DRAWS), envir = vac_store)
    for (o in outcomes) assign(key(cov, nm, o), rep(NA_real_, N_DRAWS), envir = av_store)
  }

  for (i in 1:N_DRAWS) {
    ti <- t_idx[i]
    beta_ext <- extend_beta(E$beta[ti, ])
    Rimm_i <- draw_immunity(E$foi[ti]); I0_i <- draw_I0(E$foi[ti], E$rho[ti], E$gamma[ti], E$prop_symp[ti])
    cfr_j  <- (hosp_d[i]*cfrH_d[i, ] + (1-hosp_d[i])*cfrN_d[i, ])[age_to_band]
    args_i <- list(beta_ext, E$sigma[ti], E$gamma[ti], E$rho[ti], E$prop_symp[ti],
                   Rimm_i, I0_i, ve_d[i], deliv_d[i], delay_d[i], hosp_d[i], cfr_j,
                   dwMM_d[i], dwSV_d[i], dwCH_d[i], duMM_d[i], duSV_d[i], duCH_d[i], le_d[i, ],
                   p14y_d[i], p14o_d[i], p90y_d[i], p90o_d[i])
    base_i <- do.call(run_scen, c(list(scen[[1]], 0), args_i))
    for (cov in COV_LEVELS) for (nm in vac_names) {
      s <- scen[[which(scen_names == nm)]]
      b <- do.call(run_scen, c(list(s, cov), args_i))
      av <- base_i[outcomes] - b[outcomes]
      vv <- get(keyv(cov, nm), envir = vac_store); vv[i] <- unname(b["vaccinated"])
      assign(keyv(cov, nm), vv, envir = vac_store)
      for (j in seq_along(outcomes)) {
        v <- get(key(cov, nm, outcomes[j]), envir = av_store); v[i] <- as.numeric(av)[j]
        assign(key(cov, nm, outcomes[j]), v, envir = av_store)
      }
    }
    if (i %% 100 == 0) cat("  ", i, "/", N_DRAWS, "\n")
  }

  # NNV + 95% UI via Hyolim's construction: tot_vacc (median doses, ~fixed per
  # coverage) divided by the averted-burden median and its 95% bounds, flipped
  # (NNV decreases as averted increases). If the averted 2.5% bound <= 0 (vaccine
  # not reliably beneficial in the tail), the UPPER NNV bound is unbounded -> NA.
  nnv_mc <- do.call(rbind, lapply(COV_LEVELS, function(cov) {
    do.call(rbind, lapply(vac_names, function(nm) {
      s <- scen[[which(scen_names == nm)]]
      tot_vacc <- median(get(keyv(cov, nm), envir = vac_store), na.rm = TRUE)
      do.call(rbind, lapply(outcomes, function(o) {
        av <- get(key(cov, nm, o), envir = av_store)
        aq <- quantile(av, c(.5, .025, .975), na.rm = TRUE)   # averted median, lo, hi
        nnv_med <- if (aq[1] > 0) tot_vacc / aq[1] else NA_real_
        nnv_lo  <- if (aq[3] > 0) tot_vacc / aq[3] else NA_real_   # <- averted HI
        nnv_hi  <- if (aq[2] > 0) tot_vacc / aq[2] else NA_real_   # <- averted LO
        data.frame(coverage = cov, timing = s$timing, arm = s$arm, outcome = o,
                   tot_vacc = tot_vacc,
                   averted_med = aq[1], averted_lo = aq[2], averted_hi = aq[3],
                   nnv_med = nnv_med, nnv_lo = nnv_lo, nnv_hi = nnv_hi,
                   nnv_str = if (is.na(nnv_med)) "n/e"
                             else sprintf("%s (%s - %s)",
                                   formatC(round(nnv_med,1), big.mark=",", format="f", digits=1),
                                   formatC(round(nnv_lo,1),  big.mark=",", format="f", digits=1),
                                   if (is.na(nnv_hi)) "Inf" else formatC(round(nnv_hi,1), big.mark=",", format="f", digits=1)),
                   row.names = NULL)
      }))
    }))
  }))
}

# ------------------------------------------------------------
# 9. Assemble the plotting/table frame (MC medians if available, else point est.)
# ------------------------------------------------------------
out_labs <- c(symptomatic = "Symptomatic case", hospitalisations = "Hospitalisation",
              deaths = "Death", daly = "DALY")
cov_labs <- setNames(sprintf("%d%%", round(100*COV_LEVELS)), COV_LEVELS)

if (!is.null(nnv_mc)) {
  plt <- transform(nnv_mc, med = nnv_med, lo = nnv_lo, hi = nnv_hi)
} else {
  plt <- transform(nnv_pt, med = nnv, lo = NA_real_, hi = NA_real_)
}
plt$outcome  <- factor(out_labs[plt$outcome], levels = unname(out_labs))
plt$timing   <- factor(plt$timing, levels = names(timings))
plt$arm      <- factor(plt$arm, levels = arm_names)
plt$cov_lab  <- factor(cov_labs[as.character(plt$coverage)], levels = unname(cov_labs))

# ------------------------------------------------------------
# 10. Figure: NNV to avert each outcome (mirrors the reference layout).
#     rows = outcome, cols = vaccine-protection arm; x = timing; fill = coverage;
#     log10 y (NNV spans orders of magnitude across outcomes). Central 30% highlighted.
# ------------------------------------------------------------
cov_cols <- setNames(c("#3182bd", "#31a354", "#de2d26"), unname(cov_labs))   # 20 / 30 / 40 %
p_nnv <- ggplot(plt, aes(timing, med, fill = cov_lab)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45", linewidth = .3) +   # break-even: 1 dose = 1 averted
  geom_col(position = position_dodge(.8), width = .72, colour = "grey30", linewidth = .2) +
  { if (!is.null(nnv_mc))
      geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.8),
                    width = .25, linewidth = .35) } +
  facet_grid(outcome ~ arm, scales = "free_y", switch = "y") +
  scale_fill_manual(values = cov_cols, name = "Vaccine coverage") +
  scale_y_log10(breaks = scales::breaks_log(n = 6),
                labels = scales::label_number(drop0trailing = TRUE)) +
  labs(x = NULL, y = "NNV to avert one outcome",
       title = "Caldas Novas CHIKV Number Needed to Vaccinate") +
       # subtitle = sprintf("Pre-outbreak rollout, 52-wk window; central coverage %d%%, extremes %s%s",
       #                    round(100*COV_CENTRAL),
       #                    paste(sprintf("%d%%", round(100*setdiff(COV_LEVELS, COV_CENTRAL))), collapse=" / "),
       #                    if (is.null(nnv_mc)) " (point estimate)" else " (median + 95% UI)"),
       # caption = "Dashed line = NNV 1 (1 dose averts 1 case); NNV < 1 means each dose averts >1 case via herd protection") +
  theme_bw(11) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9),
        strip.placement = "outside", strip.background = element_rect(fill = "grey92", colour = NA),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),   # single timing -> no x label
        legend.position = "bottom", panel.grid.minor = element_blank())
print(p_nnv)
ggsave("CHIKV_ca_nnv.png", p_nnv, width = 7, height = 8.5, dpi = 120)

# ------------------------------------------------------------
# 11. Console + Excel output
# ------------------------------------------------------------
nnv_tbl <- if (!is.null(nnv_mc)) {
  transform(nnv_mc[, c("coverage","timing","arm","outcome","nnv_str")],
            outcome = out_labs[outcome])
} else {
  transform(nnv_pt, NNV = round(nnv, 1), averted = round(averted, 1),
            vaccinated = round(vaccinated),
            outcome = out_labs[outcome])[, c("coverage","timing","arm","outcome",
                                             "vaccinated","averted","NNV")]
}
cat("\n=== NNV to avert one outcome (", if (is.null(nnv_mc)) "point estimate" else "median, 95% UI", ") ===\n", sep="")
ow <- getOption("width"); options(width = 250); print(nnv_tbl, row.names = FALSE); options(width = ow)

notes <- data.frame(
  parameter = c("Location / pathogen", "Scenario", "Evaluation window", "NNV definition",
                "Numerator (people vaccinated)", "Denominator (outcomes averted)",
                "UI construction", "Outcomes", "Coverage sweep",
                "Vaccine efficacy", "Target group", "Uncertainty", "DALY structure"),
  value = c("Caldas Novas (Goias), CHIKV",
            "Pre-outbreak vaccine rollout (2025-W40)",
            sprintf("52 epi-weeks 2025-W40 -> 2026-W38 (index %d-%d); burden accrued only within window",
                    min(EVAL_WIN), max(EVAL_WIN)),
            "NNV = people vaccinated / outcomes averted vs no-vaccine baseline (population-level; incl. herd effects)",
            "target population x coverage (eligible 18-59 pop x coverage; = Hyolim n*coverage)",
            "52-wk baseline burden - 52-wk scenario burden, per outcome",
            "median NNV = tot_vacc / averted_median; 95% UI = tot_vacc / averted_{97.5%,2.5%} (flipped, per Hyolim)",
            "symptomatic case, hospitalisation, death, DALY",
            paste(sprintf("%d%%", round(100*COV_LEVELS)), collapse=", "),
            "Beta(98.9%, 96.7-99.8), median 0.989",
            "18-59 y (age_df groups 4-8)",
            if (is.null(nnv_mc)) "point estimate (median inputs)"
            else sprintf("%d-draw Latin-hypercube MC (vaccine + severity + DALY inputs, resampled transmission draws)", N_DRAWS),
            sprintf("%s YLD + undiscounted YLL (Hyolim Kang CHIK_VIM)", PHASE_MODE)),
  stringsAsFactors = FALSE)

sheets <- list(notes = notes, nnv = nnv_tbl)
if (!is.null(nnv_mc)) sheets$nnv_point_estimate <- transform(nnv_pt, NNV = round(nnv, 1),
  outcome = out_labs[outcome])[, c("coverage","timing","arm","outcome","vaccinated","averted","NNV")]
write_xlsx(sheets, "CHIKV_ca_nnv_outputs.xlsx")
cat("\nWrote CHIKV_ca_nnv_outputs.xlsx and CHIKV_ca_nnv.png\n")

saveRDS(list(nnv_pt = nnv_pt, nnv_mc = nnv_mc, plt = plt, COV_LEVELS = COV_LEVELS,
             outcomes = outcomes, out_labs = out_labs, RUN_MC = RUN_MC,
             N_DRAWS = if (RUN_MC) N_DRAWS else NA, PHASE_MODE = PHASE_MODE),
        "CHIKV_ca_nnv_results.rds")
cat("Saved CHIKV_ca_nnv_results.rds\n")
