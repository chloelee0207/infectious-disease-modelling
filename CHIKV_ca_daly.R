# ============================================================
# CHIKV_ca_daly.R -- Caldas Novas CHIKV DALYs, STANDALONE on the LHS ensemble.
# ------------------------------------------------------------
# Estimates YLD + YLL + DALYs (median, 95% UI) for the no-vaccine baseline and each
# vaccination scenario, plus DALYs AVERTED per timing x arm. Follows Hyolim Kang's
# CHIK_VIM DALY structure (github.com/hyolimkang/CHIK_VIM):
#     DALY = YLD + YLL
#     YLD  = sum over disease phases of  (cases in phase) x DW x duration
#     YLL  = deaths x remaining life-expectancy at age of death   (no discounting)
#
# Disease progression + all DALY inputs (with medians & 95% UIs) come from
# disease_progression.xlsx; the decision tree in chikv_decision_tree_A4_landscape.svg
# is the branch-probability map. Uncertainty is a single 1000-draw Latin-hypercube MC:
# each draw resamples one FEASIBLE transmission draw from the ensemble (beta_t, gamma,
# sigma, rho, prop_symp, FOI -> immunity) and layers an LHS row over the vaccine
# parameters AND every disease-progression input (DWs, durations, life-years, CFR,
# recovery probabilities), so calibration + programme + severity uncertainty propagate
# jointly. Point estimate is not needed here (the MC median is reported).
#
# This script MIRRORS the SEIR + scenario setup of CHIKV_ca_vacc.R (standalone by
# design, like that engine); it re-runs the scenarios to recover the AGE-STRATIFIED
# symptomatic that the engine's results RDS collapses to totals. Keep the SEIR block
# in sync with CHIKV_ca_vacc.R if that engine changes.
#
# Run order:  source("ca_common.R"); source("CHIKV_ca_lhs.R")  # ensemble (slow, once)
#             source("CHIKV_ca_daly.R")                          # this file
# ============================================================
library(dplyr); library(tidyr); library(ggplot2); library(writexl)
source("ca_common.R")   # fmtq, qs, burden, load_burden_params, load_caldas_age_cases, compute_age_weight

# ---- knobs --------------------------------------------------
N_DRAWS   <- 1000                 # LHS draws (per request)
PHASE_MODE <- "hosp_severity"     # Hyolim severity axis = hospitalisation: hospitalised -> severe DW,
                                  #   non-hosp -> mild/mod DW for the pre-chronic (acute/sub-acute)
                                  #   episode; chronic (not recovered by 90d) -> chronic DW.
                                  # "three_phase"  : acute=mild/mod (all), sub-acute=severe, chronic.
                                  # "acute_split"  : acute = SEV_SPLIT mild/mod + severe; + chronic (14d gate).
SEV_SPLIT <- c(mild_mod = 0.47, severe = 0.53)   # fixed acute severity split (SVG sketch); used only by "acute_split"
set.seed(2028)

# ------------------------------------------------------------
# 0. Load the LHS ensemble (same object CHIKV_ca_vacc.R consumes)
# ------------------------------------------------------------
if (!file.exists("CHIKV_ca_lhs_ensemble.rds"))
  stop("CHIKV_ca_lhs_ensemble.rds not found -- run CHIKV_ca_lhs.R first.")
E <- readRDS("CHIKV_ca_lhs_ensemble.rds")
N <- E$N; A <- E$A; age_df <- E$age_df; age_mid <- E$age_mid
observed_cases <- E$observed_cases; caldas_obs <- E$caldas_obs
week_1_cases <- E$week_1_cases; T_data <- E$T_weeks
n_ens <- nrow(E$beta)
E0 <- rep(0, A)
cat(sprintf("Loaded ensemble: %d feasible transmission draws; running %d-draw DALY LHS.\n",
            n_ens, N_DRAWS))

# ------------------------------------------------------------
# 1. SEIRV simulator (age-structured, weekly) -- mirrors CHIKV_ca_vacc.R
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
# 2. Horizon + per-draw immunity/seed (mirrors CHIKV_ca_vacc.R)
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
immun_delay <- 2
idx_of <- function(yr, wk) caldas_obs$week_index[caldas_obs$Year == yr & caldas_obs$week == wk]
start_s1 <- idx_of(2026, 16)     # IXCHIQ real rollout
start_s2 <- idx_of(2026, 1)      # start of 2026
start_s3 <- idx_of(2025, 40)     # pre-outbreak
timings <- list("actual rollout" = start_s1, "start of 2026" = start_s2, "pre-outbreak" = start_s3)
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

# age-class index sets from the decadal band map: <40 = bands 1-4, >=40 = bands 5-9
young_idx <- which(age_to_band <= 4)     # groups mapping to [0,40)
old_idx   <- which(age_to_band >= 5)     # groups mapping to [40,90)

# ------------------------------------------------------------
# 5. DALY parameters from disease_progression.xlsx (Lognormal/Beta, with 95% UIs)
#    Returns the distribution hyper-parameters used by the LHS draws below.
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

  # remaining life-years by decadal band [0,10)..[80,90) (length 9), ordered by lower bound
  le <- row("Remaining life-years")
  lo <- as.numeric(sub(".*\\[\\s*([0-9]+).*", "\\1", le$group)); le <- le[order(lo), ]
  le_ms <- list(m = le$v1, s = le$v2, med = le$median)
  stopifnot(length(le_ms$m) == 9)

  # acute (14d) and sub-acute (90d) recovery probabilities, by age class
  p14_y <- beta_ab("Probability of recovery within 14 days after onset of symptoms", "< 40")
  p14_o <- beta_ab("Probability of recovery within 14 days after onset of symptoms", "> 40")
  p90_y <- beta_ab("Probability of recovery within 90 days after acute period", "< 40")
  p90_o <- beta_ab("Probability of recovery within 90 days after acute period", "> 40")

  list(dw_mm=dw_mm, dw_sev=dw_sev, dw_chr=dw_chr,
       du_mm=du_mm, du_sev=du_sev, du_chr=du_chr, le=le_ms,
       p14_y=p14_y, p14_o=p14_o, p90_y=p90_y, p90_o=p90_o)
}
dp <- load_daly_params()

# ------------------------------------------------------------
# 6. Vaccine-parameter distributions (Hyolim Table 1; coverage median 30%)
# ------------------------------------------------------------
beta_from_ci <- function(m, lo, hi) { v <- ((hi-lo)/(2*1.96))^2; k <- m*(1-m)/v-1; c(a=m*k, b=(1-m)*k) }
ve_ab    <- beta_from_ci(0.989, 0.967, 0.998)
cov_ab   <- beta_from_ci(0.30,  0.20,  0.40)
deliv_ab <- beta_from_ci(0.10,  0.09,  0.11)

# ------------------------------------------------------------
# 7. Build the 1000-draw Latin-hypercube design
#    Columns: 4 vaccine + 1 hosp + 9 cfr_hosp + 9 cfr_nonhosp + 3 DW + 3 dur + 9 LE
#             + 4 recovery probs = 42 stratified uniforms, mapped through inverse CDFs.
#    Transmission draws are bootstrap-resampled from the feasible ensemble (n<1000).
# ------------------------------------------------------------
lhs_col <- function(n) (sample.int(n) - runif(n)) / n
K <- 42
U <- sapply(1:K, function(j) lhs_col(N_DRAWS))
col <- 0; nextU <- function(w = 1) { idx <- (col + 1):(col + w); col <<- col + w; U[, idx, drop = FALSE] }

ve_d    <- qbeta(nextU(), ve_ab["a"],    ve_ab["b"])
cov_d   <- qbeta(nextU(), cov_ab["a"],   cov_ab["b"])
deliv_d <- qbeta(nextU(), deliv_ab["a"], deliv_ab["b"])
delay_d <- 1 + round(2 * nextU())                                   # {1,2,3}, median 2
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
p14y_d  <- qbeta(nextU(), dp$p14_y$a, dp$p14_y$b)
p14o_d  <- qbeta(nextU(), dp$p14_o$a, dp$p14_o$b)
p90y_d  <- qbeta(nextU(), dp$p90_y$a, dp$p90_y$b)
p90o_d  <- qbeta(nextU(), dp$p90_o$a, dp$p90_o$b)
stopifnot(col == K)
t_idx   <- sample.int(n_ens, N_DRAWS, replace = TRUE)               # resampled transmission draw per LHS row

# ------------------------------------------------------------
# 8. YLD / YLL / DALY for one scenario run, one draw.
#    symp_age = TRUE symptomatic by model age group (length A). Reweighted to the
#    observed case-age split (age_weight, total-preserving) so the <40/>=40 chronic
#    proportions and the deaths-by-band both match CHIKV_ca_vacc.R's death calc.
# ------------------------------------------------------------
daly_one <- function(symp_age, cfr_vec_j, hosp_j, le_band, dwmm, dwsv, dwch, dumm, dusv, duch,
                     p14y, p14o, p90y, p90o) {
  symp_w  <- symp_age * age_weight
  symp_dw <- if (sum(symp_w) > 0) symp_w * (sum(symp_age)/sum(symp_w)) else symp_age
  sy <- sum(symp_dw[young_idx]); so <- sum(symp_dw[old_idx]); st <- sy + so

  if (PHASE_MODE == "hosp_severity") {
    # Hyolim severity axis = hospitalisation. Pre-chronic (acute/sub-acute) episode:
    # hospitalised (fraction hosp_j) accrue severe DW/duration, non-hosp accrue mild/mod.
    # Chronic = not recovered by 90d (after the sub-acute period) x chronic DW.
    yld_ac  <- st * ((1-hosp_j)*dwmm*dumm + hosp_j*dwsv*dusv)
    yld_sub <- 0
    yld_chr <- (sy*(1-p14y)*(1-p90y) + so*(1-p14o)*(1-p90o)) * dwch * duch
  } else if (PHASE_MODE == "three_phase") {
    # acute = mild/mod (ALL symptomatic); sub-acute = severe (not recovered by 14d);
    # chronic (not recovered by 90d). Recovery probs are age-class specific.
    yld_ac  <- st * dwmm * dumm
    yld_sub <- (sy*(1-p14y) + so*(1-p14o)) * dwsv * dusv
    yld_chr <- (sy*(1-p14y)*(1-p90y) + so*(1-p14o)*(1-p90o)) * dwch * duch
  } else {                       # "acute_split": acute severity split 47/53; + chronic (14d gate)
    yld_ac  <- st * (SEV_SPLIT["mild_mod"]*dwmm*dumm + SEV_SPLIT["severe"]*dwsv*dusv)
    yld_sub <- 0
    yld_chr <- (sy*(1-p14y) + so*(1-p14o)) * dwch * duch
  }
  yld <- yld_ac + yld_sub + yld_chr

  deaths_age <- symp_dw * cfr_vec_j                 # deaths by model age group
  yll        <- sum(deaths_age * le_band[age_to_band])
  c(yld_acute=unname(yld_ac), yld_subacute=unname(yld_sub), yld_chronic=unname(yld_chr),
    yld=unname(yld), yll=yll, daly=unname(yld)+yll, deaths=sum(deaths_age))
}

# ------------------------------------------------------------
# 9. Monte Carlo: 1000 draws x 7 scenarios
# ------------------------------------------------------------
metrics <- c("yld_acute","yld_subacute","yld_chronic","yld","yll","daly","deaths")
D <- setNames(lapply(scen_names, function(x) matrix(NA_real_, N_DRAWS, length(metrics),
              dimnames=list(NULL, metrics))), scen_names)

cat(sprintf("Running %d draws (7 scenarios each, %s YLD)...\n", N_DRAWS, PHASE_MODE))
for (i in 1:N_DRAWS) {
  ti <- t_idx[i]
  beta_ext <- extend_beta(E$beta[ti, ])
  Rimm_i <- draw_immunity(E$foi[ti]); I0_i <- draw_I0(E$foi[ti], E$rho[ti], E$gamma[ti], E$prop_symp[ti])
  cfr_j  <- (hosp_d[i]*cfrH_d[i, ] + (1-hosp_d[i])*cfrN_d[i, ])[age_to_band]
  for (s in scen) {
    if (s$type=="base") { covv<-0; vi<-0; vb<-0; st<-start_s3 }
    else { covv<-cov_d[i]; vb<-ve_d[i]; vi<-if(s$type=="both")ve_d[i] else 0
           st<-min(s$start+delay_d[i], T_sim) }
    out <- seirv_vaccinated(T_sim, A, N, Rimm_i, I0_i, E0, beta_ext, E$sigma[ti], E$gamma[ti],
             E$rho[ti], target_age, covv, deliv_d[i], st, vi, vb, immun_delay, prop_symp=E$prop_symp[ti])
    D[[s$name]][i, ] <- daly_one(rowSums(out$new_symptomatic), cfr_j, hosp_d[i], le_d[i, ],
                                 dwMM_d[i], dwSV_d[i], dwCH_d[i], duMM_d[i], duSV_d[i], duCH_d[i],
                                 p14y_d[i], p14o_d[i], p90y_d[i], p90o_d[i])
  }
  if (i %% 100 == 0) cat("  ", i, "/", N_DRAWS, "\n")
}

# ------------------------------------------------------------
# 10. Summaries: per-scenario burden + DALYs averted vs baseline
# ------------------------------------------------------------
base <- D[["No vaccine (baseline)"]]
cat("\n=== Baseline DALYs, no vaccine (median, 95% UI) ===\n")
for (m in c("yld","yll","daly","deaths"))
  cat(sprintf("  %-8s %s\n", m, fmtq(base[, m], if (m=="deaths") 1 else 0)))

daly_by_scenario <- do.call(rbind, lapply(scen_names, function(nm) {
  m <- D[[nm]]; b <- nm == "No vaccine (baseline)"
  data.frame(timing = if (b) "No vaccine" else sub(" \\|.*","",nm),
             arm    = if (b) "No vaccine" else sub(".*\\| ","",nm),
             YLD  = fmtq(m[,"yld"]), YLL = fmtq(m[,"yll"]), DALY = fmtq(m[,"daly"]),
             deaths = fmtq(m[,"deaths"], 1),
             YLD_acute = fmtq(m[,"yld_acute"]), YLD_subacute = fmtq(m[,"yld_subacute"]),
             YLD_chronic = fmtq(m[,"yld_chronic"]), row.names = NULL, check.names = FALSE)
}))

# averted (baseline - scenario), per draw -> preserves the paired correlation
av <- setNames(lapply(vac_names, function(nm) base - D[[nm]]), vac_names)
daly_averted <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- av[[nm]]
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
             YLD_averted = fmtq(m[,"yld"]), YLL_averted = fmtq(m[,"yll"]),
             DALY_averted = fmtq(m[,"daly"]),
             pct_DALY = sprintf("%.1f%%", 100*median(m[,"daly"]/base[,"daly"])),
             row.names = NULL, check.names = FALSE)
}))
cat("\n=== DALYs averted vs no-vaccine baseline (median, 95% UI) ===\n")
ow <- getOption("width"); options(width=250); print(daly_averted, row.names=FALSE); options(width=ow)

# ------------------------------------------------------------
# 11. Excel workbook
# ------------------------------------------------------------
notes <- data.frame(
  parameter = c("Location / pathogen", "DALY framework", "Discounting / age-weighting",
                "YLD phase structure", "YLL", "Uncertainty", "Transmission draws",
                "Severity / CFR", "Recovery probabilities", "Life expectancy", "Source sheet"),
  value = c("Caldas Novas (Goias), CHIKV",
            "DALY = YLD + YLL (Hyolim Kang CHIK_VIM structure)",
            "none (0% discount, no age-weighting; GBD 2010+)",
            switch(PHASE_MODE,
              hosp_severity = "severity = hospitalisation: hospitalised -> severe DW, non-hosp -> mild/mod DW (pre-chronic acute/sub-acute episode); chronic (not recovered by 90d) x DW 0.317",
              three_phase   = "acute = mild/mod (all symptomatic); sub-acute = severe (not recovered by 14d); chronic (not recovered by 90d) x DW 0.317",
              acute_split   = sprintf("acute = %.0f%% mild/mod + %.0f%% severe; chronic (not recovered by 14d) x DW 0.317", 100*SEV_SPLIT[1], 100*SEV_SPLIT[2])),
            "deaths (by decadal band) x remaining life-years at age of death",
            sprintf("%d-draw Latin-hypercube MC over vaccine + disease-progression inputs, paired with resampled transmission draws", N_DRAWS),
            sprintf("%d feasible ensemble draws (beta_t, gamma, sigma, rho, prop_symp, FOI), resampled with replacement", n_ens),
            "hosp among symptomatic + CFR hosp/non-hosp by decadal band (Beta); deaths reweighted to observed case-age split",
            "acute (14d) & sub-acute (90d) recovery, age <40 / >=40 (Beta)",
            "remaining life-years by decadal band (Lognormal)",
            "disease_progression.xlsx"),
  stringsAsFactors = FALSE)

sheets <- list(notes = notes, daly_by_scenario = daly_by_scenario, daly_averted = daly_averted)
write_xlsx(sheets, "CHIKV_ca_daly_outputs.xlsx")
cat("\nWrote CHIKV_ca_daly_outputs.xlsx (sheets:", paste(names(sheets), collapse=", "), ")\n")

# ------------------------------------------------------------
# 12. Figures: (a) baseline YLD/YLL composition; (b) DALYs averted by timing x arm
# ------------------------------------------------------------
comp <- data.frame(component = c("YLD (acute)","YLD (sub-acute)","YLD (chronic)","YLL"),
                   value = c(median(base[,"yld_acute"]), median(base[,"yld_subacute"]),
                             median(base[,"yld_chronic"]), median(base[,"yll"])))
comp$component <- factor(comp$component, levels = comp$component)
p_comp <- ggplot(comp, aes(x = "Baseline", y = value, fill = component)) +
  geom_col(width = .55) +
  scale_fill_manual(values = c("#c6dbef","#6baed6","#2171b5","#d6604d"), name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "DALYs (median)",
       title = "Caldas Novas CHIKV: baseline DALY composition (no vaccine)",
       caption = sprintf("%s YLD structure; YLL undiscounted; %d-draw LHS", PHASE_MODE, N_DRAWS)) +
  theme_bw(11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
print(p_comp); ggsave("CHIKV_ca_daly_composition.png", p_comp, width = 5.5, height = 5, dpi = 120)

av_long <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- av[[nm]]
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
             med = median(m[,"daly"]), lo = quantile(m[,"daly"],.025), hi = quantile(m[,"daly"],.975),
             row.names = NULL)
}))
av_long$timing <- factor(av_long$timing, levels = names(timings))
av_long$arm    <- factor(av_long$arm, levels = arm_names)
p_av <- ggplot(av_long, aes(timing, med, fill = arm)) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.8), width = .25, linewidth = .4) +
  scale_fill_manual(values = c("Disease-blocking"="#9ecae1","Disease + infection blocking"="#08519c"), name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "DALYs averted (vs no vaccine), median + 95% UI",
       title = "Caldas Novas CHIKV: DALYs averted by vaccination timing") +
  theme_bw(11) + theme(plot.title = element_text(face="bold", hjust=.5),
                       axis.text.x = element_text(angle = 20, hjust = 1),
                       legend.position = "bottom", panel.grid.minor = element_blank())
print(p_av); ggsave("CHIKV_ca_daly_averted.png", p_av, width = 8, height = 4.5, dpi = 120)

saveRDS(list(D = D, av = av, scen_names = scen_names, vac_names = vac_names,
             metrics = metrics, PHASE_MODE = PHASE_MODE, N_DRAWS = N_DRAWS),
        "CHIKV_ca_daly_results.rds")
cat("Saved figures CHIKV_ca_daly_composition.png, CHIKV_ca_daly_averted.png and CHIKV_ca_daly_results.rds\n")
