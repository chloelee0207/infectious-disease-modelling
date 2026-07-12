# ============================================================
# Caldas Novas CHIKV vaccine impact -- STANDALONE on the LHS ensemble.
# ------------------------------------------------------------
# Does NOT depend on CHIKV_ca_pre_vacc_optim.R. It loads CHIKV_ca_lhs_ensemble.rds
# (produced by CHIKV_ca_lhs.R), which carries, per feasible LHS draw, a re-fitted
# beta_t plus that draw's gamma, sigma, rho ~ Beta(20,60), prop_symp ~ Beta(35.84,32.56)
# and FOI (-> catalytic prior immunity). The vaccine Monte Carlo REUSES those draws
# (no re-fitting) and layers an LHS over the vaccine parameters + severity draws, so
# all uncertainty is propagated jointly. Point estimate = median inputs (rho 0.25).
#
# Two efficacy arms (Hyolim): Disease-blocking only (VE_inf 0, VE_block = ve) and
# Disease + infection blocking (VE_inf = ve = VE_block). Single shared efficacy per draw.
#
# Run order:  source("ca_common.R"); source("CHIKV_ca_lhs.R")  # makes the RDS (slow, once)
#             source("CHIKV_ca_vacc.R")                          # this file (fast)
# ============================================================
library(dplyr); library(tidyr); library(ggplot2)
source("ca_common.R")   # fmtq, burden, load_burden_params, load_caldas_age_cases, compute_age_weight

# ------------------------------------------------------------
# 0. Load the LHS ensemble
# ------------------------------------------------------------
if (!file.exists("CHIKV_ca_lhs_ensemble.rds"))
  stop("CHIKV_ca_lhs_ensemble.rds not found -- run CHIKV_ca_lhs.R first.")
E <- readRDS("CHIKV_ca_lhs_ensemble.rds")
N <- E$N; A <- E$A; age_df <- E$age_df; age_mid <- E$age_mid
observed_cases <- E$observed_cases; caldas_obs <- E$caldas_obs
week_1_cases <- E$week_1_cases; T_data <- E$T_weeks
n_draws <- nrow(E$beta)                 # feasible LHS draws
E0 <- rep(0, A)
cat(sprintf("Loaded ensemble: %d feasible draws, %d-week data window (obs total %d).\n",
            n_draws, T_data, sum(observed_cases)))

# ------------------------------------------------------------
# 1. SEIRV simulator with vaccination (age-structured, weekly)
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
# 2. Extended simulation horizon (so the flatter infection-blocking curve completes)
# ------------------------------------------------------------
EXTEND <- 26                            # weeks past the 52-week data window (adjustable)
T_sim  <- T_data + EXTEND
extend_beta <- function(b) c(b, rep(b[length(b)], EXTEND))   # hold beta flat past the data
sim_weeks  <- 1:T_sim

# per-draw immunity (catalytic FOI) and seed
draw_immunity <- function(foi) 1 - exp(-foi * age_mid)
draw_I0 <- function(foi, rho, gamma, ps) {
  Rimm <- draw_immunity(foi); sf <- N*(1-Rimm)/sum(N*(1-Rimm))
  round(((week_1_cases/rho/ps)/gamma) * sf)
}

# ------------------------------------------------------------
# 3. Vaccine programme + scenarios (52-week index; starts from caldas_obs)
# ------------------------------------------------------------
target_age <- rep(0, A); target_age[c(4,5,6,7,8)] <- 1     # eligible = 18-59 (age_df groups 4-8)
target_pop <- sum(N[target_age == 1])
immun_delay <- 2                                            # dose -> immunity (fixed, no UI)
idx_of <- function(yr, wk) caldas_obs$week_index[caldas_obs$Year == yr & caldas_obs$week == wk]
start_s1 <- idx_of(2026, 16)     # IXCHIQ real rollout (2026-W16)
start_s2 <- idx_of(2026, 1)      # start of 2026
start_s3 <- idx_of(2025, 40)     # pre-outbreak (2025-W40, per request)
timings <- list("actual rollout" = start_s1, "start of 2026" = start_s2, "pre-outbreak" = start_s3)
cat(sprintf("Scenario starts (week_index): actual rollout=%d, start of 2026=%d, pre-outbreak=%d\n",
            start_s1, start_s2, start_s3))
arm_names <- c("Disease-blocking", "Disease + infection blocking")

# scenario table: baseline + timing x arm
scen <- list(list(name="No vaccine (baseline)", timing="No vaccine", arm="No vaccine",
                  start=NA_integer_, type="base"))
for (tn in names(timings)) for (an in arm_names)
  scen[[length(scen)+1]] <- list(name=paste0(tn," | ",an), timing=tn, arm=an,
                                 start=timings[[tn]], type=if(an=="Disease-blocking")"disb" else "both")
scen_names <- vapply(scen, function(s) s$name, character(1))
vac_names  <- setdiff(scen_names, "No vaccine (baseline)")

# ------------------------------------------------------------
# 4. Burden severity params + observed-age (w_a) death correction
# ------------------------------------------------------------
invisible(list2env(load_burden_params(A), globalenv()))     # ps_*, hosp_*, cfr_*, age_to_band, cfr_vec
obs_band_prop <- load_caldas_age_cases()$obs_band_prop
# w_a from the POINT-estimate no-vaccine run (median inputs)
I0_base   <- draw_I0(E$base_foi, E$base_rho, E$base_gamma, E$base_prop_symp)
Rimm_base <- draw_immunity(E$base_foi)
beta_base_ext <- extend_beta(E$base_beta)
out_base_pt <- seirv_vaccinated(T_sim, A, N, Rimm_base, I0_base, E0, beta_base_ext,
                 E$base_sigma, E$base_gamma, E$base_rho, target_age, 0, 0.10, start_s3,
                 0, 0, immun_delay, prop_symp = E$base_prop_symp)
age_weight <- compute_age_weight(rowSums(out_base_pt$new_infections), obs_band_prop, age_to_band)
cat(sprintf("Observed-age death weight w_a range [%.2f, %.2f].\n", min(age_weight), max(age_weight)))

# ------------------------------------------------------------
# 5. Vaccine-parameter distributions (Hyolim Table 1; coverage median 30% per request)
# ------------------------------------------------------------
beta_from_ci <- function(m, lo, hi) { v <- ((hi-lo)/(2*1.96))^2; k <- m*(1-m)/v-1; c(a=m*k, b=(1-m)*k) }
ve_ab    <- beta_from_ci(0.989, 0.967, 0.998)   # efficacy (shared VE_inf & VE_block)
cov_ab   <- beta_from_ci(0.30,  0.20,  0.40)    # coverage of eligible (median 30%, 95% 20-40)
deliv_ab <- beta_from_ci(0.10,  0.09,  0.11)    # weekly delivery speed
ve_med <- 0.989; cov_med <- 0.30; deliv_med <- 0.10; delay_med <- 2   # deployment delay median (1-3)

# run one scenario given calibration + vaccine params (returns burden vector)
run_scen <- function(s, beta_ext, sigma, gamma, rho, ps, Rimm, I0i,
                     ve, cov, deliv, delay, hr, cv) {
  if (s$type == "base") { covv<-0; vi<-0; vb<-0; st<-start_s3 }
  else { covv<-cov; vb<-ve; vi<-if (s$type=="both") ve else 0
         st<-min(s$start + delay, T_sim) }
  out <- seirv_vaccinated(T_sim, A, N, Rimm, I0i, E0, beta_ext, sigma, gamma, rho,
           target_age, covv, deliv, st, vi, vb, immun_delay, prop_symp = ps)
  c(burden(out, hr, cv, age_weight), reported = rho * sum(out$new_symptomatic))
}

# ------------------------------------------------------------
# 6. Point estimate (median inputs) -- weekly symptomatic curves for the plots too
# ------------------------------------------------------------
outcomes <- c("infections","symptomatic","hospitalisations","deaths")
pt_burden <- t(sapply(scen, function(s)
  run_scen(s, beta_base_ext, E$base_sigma, E$base_gamma, E$base_rho, E$base_prop_symp,
           Rimm_base, I0_base, ve_med, cov_med, deliv_med, delay_med, hosp_rate, cfr_vec)))
rownames(pt_burden) <- scen_names

# ------------------------------------------------------------
# 7. Monte Carlo: reuse the LHS draws, layer an LHS over the vaccine params + severity
# ------------------------------------------------------------
set.seed(2027)
lhs_col <- function(n) (sample.int(n) - runif(n))/n
Uv <- sapply(1:4, function(j) lhs_col(n_draws))
ve_d    <- qbeta(Uv[,1], ve_ab["a"],   ve_ab["b"])
cov_d   <- qbeta(Uv[,2], cov_ab["a"],  cov_ab["b"])
deliv_d <- qbeta(Uv[,3], deliv_ab["a"],deliv_ab["b"])
delay_d <- 1 + round(2 * Uv[,4])                 # {1,2,3}, median 2

# per-draw weekly symptomatic (summed over age) for ribbon plots -- store the 3 both-arm
# timings + baseline (disease-blocking overlaps baseline for infections)
bmat <- setNames(lapply(scen_names, function(x) matrix(NA_real_, n_draws, 5,
                  dimnames=list(NULL, c(outcomes,"reported")))), scen_names)
wk_symp <- setNames(lapply(scen_names, function(x) matrix(NA_real_, n_draws, T_sim)), scen_names)

cat("Running", n_draws, "draws (7 scenarios each, extended horizon)...\n")
for (i in 1:n_draws) {
  beta_ext <- extend_beta(E$beta[i,])
  Rimm_i <- draw_immunity(E$foi[i])
  I0_i   <- draw_I0(E$foi[i], E$rho[i], E$gamma[i], E$prop_symp[i])
  hosp_i <- rbeta(1, hosp_a, hosp_b)
  cfr_i  <- (hosp_i*rbeta(9,cfr_hosp_a,cfr_hosp_b) + (1-hosp_i)*rbeta(9,cfr_nonh_a,cfr_nonh_b))[age_to_band]
  for (s in scen) {
    if (s$type=="base") { covv<-0; vi<-0; vb<-0; st<-start_s3 }
    else { covv<-cov_d[i]; vb<-ve_d[i]; vi<-if(s$type=="both")ve_d[i] else 0
           st<-min(s$start+delay_d[i], T_sim) }
    out <- seirv_vaccinated(T_sim, A, N, Rimm_i, I0_i, E0, beta_ext, E$sigma[i], E$gamma[i],
             E$rho[i], target_age, covv, deliv_d[i], st, vi, vb, immun_delay, prop_symp=E$prop_symp[i])
    bmat[[s$name]][i,] <- c(burden(out, hosp_i, cfr_i, age_weight), rho=E$rho[i]*sum(out$new_symptomatic))
    wk_symp[[s$name]][i,] <- colSums(out$new_symptomatic)
  }
  if (i %% 100 == 0) cat("  ", i, "/", n_draws, "\n")
}

# reported per draw = rho[i] * true (anchored); rho already per-draw in the ensemble
rho_i <- E$rho

# ------------------------------------------------------------
# 8. Summaries: baseline true/reported; per-scenario totals; averted vs baseline
# ------------------------------------------------------------
base_true <- bmat[["No vaccine (baseline)"]][, outcomes, drop=FALSE]
cat("\n=== Baseline burden, no vaccine: TRUE vs REPORTED (median, 95% UI) ===\n")
base_tbl <- data.frame(outcome=character(), true=character(), reported=character())
for (o in outcomes) {
  d <- if (o=="deaths") 1 else 0
  cat(sprintf("  %-16s  %-26s  %-26s\n", o, fmtq(base_true[,o], d), fmtq(rho_i*base_true[,o], d)))
  base_tbl <- rbind(base_tbl, data.frame(outcome=o, true=fmtq(base_true[,o],d), reported=fmtq(rho_i*base_true[,o],d)))
}

# averted (baseline - scenario) per draw
av <- setNames(lapply(vac_names, function(nm) base_true - bmat[[nm]][, outcomes, drop=FALSE]), vac_names)
mc_tbl <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- av[[nm]]
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
             Infections=fmtq(m[,"infections"]), Symptomatic=fmtq(m[,"symptomatic"]),
             Hospitalisations=fmtq(m[,"hospitalisations"],1), Deaths=fmtq(m[,"deaths"],2),
             pct_symp=sprintf("%.1f%%",100*median(m[,"symptomatic"]/base_true[,"symptomatic"])),
             row.names=NULL)
}))
cat("\n=== Averted vs no-vaccine baseline (median, 95% UI) ===\n")
old_w<-getOption("width"); options(width=250); print(mc_tbl, row.names=FALSE); options(width=old_w)

saveRDS(list(bmat=bmat, av=av, wk_symp=wk_symp, base_true=base_true, rho_i=rho_i,
             scen=scen, scen_names=scen_names, vac_names=vac_names, timings=timings,
             outcomes=outcomes, T_sim=T_sim, T_data=T_data, caldas_obs=caldas_obs,
             observed_cases=observed_cases, base_tbl=base_tbl, mc_tbl=mc_tbl,
             pt_burden=pt_burden, EXTEND=EXTEND),
        "CHIKV_ca_vacc_results.rds")
cat("\nEngine done. Source CHIKV_ca_outputs.R for figures + Excel (loads CHIKV_ca_vacc_results.rds).\n")
