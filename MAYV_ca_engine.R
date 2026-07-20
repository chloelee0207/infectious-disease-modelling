# ============================================================
# MAYV_ca_engine.R -- Caldas Novas MAYV UNIFIED Monte Carlo engine.
# ------------------------------------------------------------
# MAYV's own engine (the analogue of CHIKV_ca_engine.R; the CHIKV engine/LHS are
# left untouched). ONE uncertainty propagation -> burden + DALY + NNV, consistent
# draw-for-draw. It CONSUMES MAYV_ca_lhs_ensemble.rds (the rainfall-envelope outbreak
# at sampled wet-season-PEAK R0 ~ Lognormal[high], flat Lima-2021 immunity), and
# layers vaccine + severity + DALY draws on top.
#
# KEY MAYV FRAMING -- CONDITION ON AN OUTBREAK. Under seasonal-peak scaling a single
# introduction only self-sustains when peak R0 is high (~2.5-3.5); most high-scenario
# draws fizzle. So we report BOTH: P(outbreak) = share of draws that take off
# (baseline attack > OUTBREAK_ATTACK_THRESH), and the burden/DALY/NNV distributions
# CONDITIONAL on taking off (over the outbreak draws only). This is the honest read of
# a steep-threshold, single-introduction system (seeding kept deliberately conservative).
#
# BORROWED severity/DALY (no MAYV-specific data): CHIKV disease-progression params via
# load_burden_params()/load_daly_params() in ca_common.R -- a CHIKV-equivalent UPPER
# BOUND on MAYV severity (the established MAYV convention, see MAYV_ca_vacc.R).
#
# VACCINE: DISEASE-BLOCKING ONLY (VE_inf = 0 -> infections identical across arms; only
# symptomatic/hosp/deaths/DALY move), pre-outbreak campaign, coverage/VE_block/delivery
# sampled. age_weight = 1 (uniform): no observed MAYV age distribution to correct toward.
#
# Run order: source("MAYV_ca_lhs.R")  (ensemble, once);  source("MAYV_ca_engine.R")
# ============================================================
library(dplyr); library(ggplot2)
source("ca_common.R")   # fmtq, load_burden_params, load_daly_params

# ---- knobs --------------------------------------------------
PHASE_MODE             <- "hosp_severity"   # YLD severity axis (matches CHIKV engine)
OUTBREAK_ATTACK_THRESH <- 1.0               # % of susceptibles infected -> "took off"
R0_FIX                 <- 3.0               # fixed peak R0 for the representative-outbreak plot (6c); NA to skip
set.seed(2031)

# ------------------------------------------------------------
# SEIRV simulator (disease-blocking capable, seed_week-aware, FLAT immunity).
# Lifted from MAYV_ca_vacc.R so the engine is self-contained. coverage = 0 reproduces
# the pre-vacc baseline exactly. With VE_inf = 0 (our only mode) infections are
# vaccine-invariant; the vaccine only scales new_symptomatic by (1 - VE_block*coverage).
# ------------------------------------------------------------
seirv_vaccinated_MAYV <- function(
    T_weeks, A, N, R_init_prop, I0, base_beta, sigma, gamma, rho,
    target_age, total_coverage, weekly_delivery_speed, delay,
    VE_inf = 0, VE_block = 0, immun_delay = 2, prop_symp = 0.5242478,
    sub_steps = 7, E0 = rep(0, A), seed_week = 1) {
  pmax0 <- function(x) pmax(0, x); N_total <- sum(N); dt <- 1/sub_steps
  S <- E <- I <- R <- V <- matrix(0, A, T_weeks)
  V_covered <- vacc_delayed <- coverage_frac <- matrix(0, A, T_weeks)
  new_infections <- new_symptomatic <- matrix(0, A, T_weeks)
  target_idx <- which(target_age == 1); target_pop <- sum(N[target_idx])
  total_supply <- target_pop * total_coverage; weekly_dose_total <- total_supply * weekly_delivery_speed
  total_avail_age <- rep(0, A); total_avail_age[target_idx] <- total_supply * (N[target_idx]/target_pop)
  total_used_age <- rep(0, A); unvaccinated <- N
  S_now <- pmax0(N - E0 - R_init_prop*N); E_now <- E0; I_now <- rep(0, A); R_now <- R_init_prop*N; V_now <- rep(0, A)
  for (t in 1:T_weeks) {
    prev_Vc <- if (t > 1) V_covered[, t-1] else rep(0, A)
    if (t - immun_delay >= 1) {
      eff_dose <- vacc_delayed[, t-immun_delay]; immunized <- round(VE_inf*eff_dose)
      V_covered[, t] <- prev_Vc + eff_dose
    } else { immunized <- rep(0, A); V_covered[, t] <- prev_Vc }
    S_now <- pmax0(S_now - immunized); V_now <- V_now + immunized
    coverage_frac[, t] <- V_covered[, t]/N
    if (t >= delay && target_pop > 0) {
      rem <- weekly_dose_total
      for (a in target_idx) {
        alloc <- min(ceiling(weekly_dose_total*(N[a]/target_pop)), rem, unvaccinated[a],
                     total_avail_age[a]-total_used_age[a])
        if (alloc > 0) {
          prop_S <- if (N[a] > 0) S_now[a]/N[a] else 0; vacc_to_S <- round(alloc*prop_S)
          vacc_delayed[a, t] <- vacc_to_S; total_used_age[a] <- total_used_age[a]+alloc
          unvaccinated[a] <- unvaccinated[a]-alloc; rem <- rem-alloc
        }
      }
    }
    if (t == seed_week) { I_now <- I_now + I0; S_now <- pmax0(S_now - I0) }
    new_I_week <- rep(0, A); beta_t <- base_beta[t]
    for (k in 1:sub_steps) {
      foi <- beta_t*sum(I_now)/N_total
      new_E <- foi*S_now*dt; new_I <- sigma*E_now*dt; new_R <- gamma*I_now*dt
      S_now <- pmax0(S_now-new_E); E_now <- pmax0(E_now+new_E-new_I)
      I_now <- pmax0(I_now+new_I-new_R); R_now <- pmax0(R_now+new_R); new_I_week <- new_I_week+new_I
    }
    S[,t]<-S_now; E[,t]<-E_now; I[,t]<-I_now; R[,t]<-R_now; V[,t]<-V_now
    new_infections[,t] <- new_I_week
    new_symptomatic[,t] <- prop_symp*new_I_week*(1 - VE_block*coverage_frac[,t])
  }
  list(new_infections = new_infections, coverage_frac = coverage_frac)
}

# ------------------------------------------------------------
# 0. Load the MAYV LHS ensemble (transmission draws)
# ------------------------------------------------------------
if (!file.exists("MAYV_ca_lhs_ensemble.rds"))
  stop("MAYV_ca_lhs_ensemble.rds not found -- run MAYV_ca_lhs.R first.")
E <- readRDS("MAYV_ca_lhs_ensemble.rds")
N <- E$N; A <- E$A; age_df <- E$age_df; season <- E$season
seed_week <- E$seed_week; E0 <- E$E0; I0_total <- E$I0_total
T_weeks <- E$T_weeks; weeks <- E$weeks; x_ticks <- E$x_ticks; year_break <- E$year_break
n_ens <- length(E$R0); N_DRAWS <- n_ens
cat(sprintf("Loaded MAYV ensemble: %d draws | R0 scenario '%s' (peak=%s) | seed wk %d\n",
            n_ens, E$R0_scenario, E$r0_is_peak, seed_week))

EVAL_WIN <- 1:T_weeks   # outbreak resolves within the 52-week rainfall window

# ------------------------------------------------------------
# 1. Severity + DALY params (borrowed CHIKV) + eligibility + uniform age weight
# ------------------------------------------------------------
invisible(list2env(load_burden_params(A), globalenv()))   # ps_*, hosp_*, cfr_*, age_to_band, cfr_vec
dp <- load_daly_params()
age_weight <- rep(1, A)                                    # uniform: no observed MAYV age split
young_idx  <- which(age_to_band <= 4); old_idx <- which(age_to_band >= 5)

target_age <- rep(0, A); target_age[c(4,5,6,7,8)] <- 1     # eligible adults 18-59 (as CHIKV/MAYV vacc)
target_pop_elig <- sum(N[target_age == 1])
immun_delay <- 2
start_pre   <- 1        # pre-outbreak campaign from the window open (2025-W40), before the seed (wk 5)
stopifnot(start_pre < seed_week)

# ------------------------------------------------------------
# 2. Outcome extractor: symptomatic-by-age (already vaccine-adjusted) -> all outcomes.
#    Severity phases saved as COUNTS so a cost layer can multiply later (no re-run).
# ------------------------------------------------------------
OUTCOMES <- c("infections","reported","symptomatic","hospitalisations","deaths",
              "n_nonhosp","n_chronic","yld_acute","yld_chronic","yld","yll","daly","doses")
NNV_OUT  <- c("reported","symptomatic","hospitalisations","deaths","daly")

outcome_one <- function(symp_age, infections, doses, rho, hosp_j, cfr_j, le_band,
                        dwmm, dwsv, dwch, dumm, dusv, duch, p14y, p14o, p90y, p90o) {
  symp_age <- symp_age * age_weight                        # uniform (=1) here
  st <- sum(symp_age); sy <- sum(symp_age[young_idx]); so <- sum(symp_age[old_idx])
  reported <- rho * st                                     # rho * symptomatic (surveillance-visible)
  n_hosp    <- st * hosp_j
  n_nonhosp <- st * (1 - hosp_j)
  n_chronic <- sy*(1-p14y)*(1-p90y) + so*(1-p14o)*(1-p90o)
  if (PHASE_MODE == "hosp_severity") {
    yld_ac <- n_nonhosp*dwmm*dumm + n_hosp*dwsv*dusv
  } else {                                                  # three_phase
    yld_ac <- st*dwmm*dumm + (sy*(1-p14y) + so*(1-p14o))*dwsv*dusv
  }
  yld_chr <- n_chronic * dwch * duch
  deaths  <- sum(symp_age * cfr_j)
  yll     <- sum(symp_age * cfr_j * le_band[age_to_band])
  c(infections = infections, reported = reported, symptomatic = st,
    hospitalisations = n_hosp, deaths = deaths,
    n_nonhosp = n_nonhosp, n_chronic = n_chronic,
    yld_acute = unname(yld_ac), yld_chronic = unname(yld_chr), yld = unname(yld_ac+yld_chr),
    yll = yll, daly = unname(yld_ac+yld_chr) + yll, doses = doses)
}

# ------------------------------------------------------------
# 3. ONE Latin-hypercube design for the layered (vaccine + severity + DALY) draws.
#    Transmission is the paired ensemble draw i (R0/gamma/sigma/rho/prop_symp/immune).
# ------------------------------------------------------------
beta_from_ci <- function(m, lo, hi) { v <- ((hi-lo)/(2*1.96))^2; k <- m*(1-m)/v-1; c(a=m*k, b=(1-m)*k) }
cov_ab <- beta_from_ci(0.30, 0.20, 0.40)     # coverage of eligible 18-59
veb_ab <- beta_from_ci(0.50, 0.25, 0.75)     # disease-blocking efficacy (hypothetical; 25-75%)
del_ab <- beta_from_ci(0.10, 0.09, 0.11)     # weekly delivery speed

lhs_col <- function(n) (sample.int(n) - runif(n)) / n
K <- 41; U <- sapply(1:K, function(j) lhs_col(N_DRAWS)); col <- 0
nextU <- function(w = 1) { idx <- (col+1):(col+w); col <<- col + w; U[, idx, drop = FALSE] }
cov_d  <- qbeta(nextU(), cov_ab["a"], cov_ab["b"])
veb_d  <- qbeta(nextU(), veb_ab["a"], veb_ab["b"])
del_d  <- qbeta(nextU(), del_ab["a"], del_ab["b"])
hosp_d <- qbeta(nextU(), hosp_a, hosp_b)
cfrH_d <- qbeta(nextU(9), matrix(cfr_hosp_a, N_DRAWS, 9, byrow=TRUE), matrix(cfr_hosp_b, N_DRAWS, 9, byrow=TRUE))
cfrN_d <- qbeta(nextU(9), matrix(cfr_nonh_a, N_DRAWS, 9, byrow=TRUE), matrix(cfr_nonh_b, N_DRAWS, 9, byrow=TRUE))
dwMM_d <- qbeta(nextU(), dp$dw_mm$a, dp$dw_mm$b); dwSV_d <- qbeta(nextU(), dp$dw_sev$a, dp$dw_sev$b)
dwCH_d <- qbeta(nextU(), dp$dw_chr$a, dp$dw_chr$b)
duMM_d <- qlnorm(nextU(), dp$du_mm$m, dp$du_mm$s); duSV_d <- qlnorm(nextU(), dp$du_sev$m, dp$du_sev$s)
duCH_d <- qlnorm(nextU(), dp$du_chr$m, dp$du_chr$s)
le_d   <- qlnorm(nextU(9), matrix(dp$le$m, N_DRAWS, 9, byrow=TRUE), matrix(dp$le$s, N_DRAWS, 9, byrow=TRUE))
p14y_d <- qbeta(nextU(), dp$p14_y$a, dp$p14_y$b); p14o_d <- qbeta(nextU(), dp$p14_o$a, dp$p14_o$b)
p90y_d <- qbeta(nextU(), dp$p90_y$a, dp$p90_y$b); p90o_d <- qbeta(nextU(), dp$p90_o$a, dp$p90_o$b)
stopifnot(col == K)

# ------------------------------------------------------------
# 4. Monte Carlo: one SEIRV run per draw (VE_inf = 0 -> infections vaccine-invariant),
#    so baseline & vaccine symptomatic both come from the SAME run's new_infections.
# ------------------------------------------------------------
scen_names <- c("No vaccine (baseline)", "Pre-outbreak | Disease-blocking")
vac_name   <- "Pre-outbreak | Disease-blocking"
per_draw <- setNames(lapply(scen_names, function(x)
  matrix(NA_real_, N_DRAWS, length(OUTCOMES), dimnames = list(NULL, OUTCOMES))), scen_names)
attack_base <- numeric(N_DRAWS)
wk_base <- wk_vacc <- matrix(NA_real_, N_DRAWS, T_weeks)   # weekly symptomatic (for the epicurve)

cat(sprintf("Running %d draws (baseline + pre-outbreak disease-blocking), conditioning threshold attack > %.1f%%...\n",
            N_DRAWS, OUTBREAK_ATTACK_THRESH))
for (i in 1:N_DRAWS) {
  R0i<-E$R0[i]; gi<-E$gamma[i]; si<-E$sigma[i]; ri<-E$rho[i]; psi<-E$prop_symp[i]; immi<-E$immune_frac[i]
  Rimm <- rep(immi, A); sus <- N*(1-Rimm); I0i <- I0_total * sus/sum(sus)
  base_beta <- R0i * gi * season
  run <- seirv_vaccinated_MAYV(T_weeks, A, N, Rimm, I0i, base_beta, si, gi, ri,
           target_age, cov_d[i], del_d[i], start_pre, VE_inf = 0, VE_block = veb_d[i],
           immun_delay = immun_delay, prop_symp = psi, E0 = E0, seed_week = seed_week)
  ninf <- run$new_infections; covf <- run$coverage_frac
  inf_tot <- sum(ninf[, EVAL_WIN, drop = FALSE])
  attack_base[i] <- 100 * inf_tot / sum(sus)

  symp_base <- rowSums((psi * ninf)[, EVAL_WIN, drop = FALSE])                      # no vaccine
  symp_vacc <- rowSums((psi * ninf * (1 - veb_d[i]*covf))[, EVAL_WIN, drop = FALSE]) # disease-blocked
  wk_base[i, ] <- colSums(psi * ninf)                                              # weekly (all ages)
  wk_vacc[i, ] <- colSums(psi * ninf * (1 - veb_d[i]*covf))
  cfr_j <- (hosp_d[i]*cfrH_d[i,] + (1-hosp_d[i])*cfrN_d[i,])[age_to_band]
  args_daly <- list(le_d[i,], dwMM_d[i], dwSV_d[i], dwCH_d[i], duMM_d[i], duSV_d[i], duCH_d[i],
                    p14y_d[i], p14o_d[i], p90y_d[i], p90o_d[i])
  per_draw[["No vaccine (baseline)"]][i, ] <-
    do.call(outcome_one, c(list(symp_base, inf_tot, 0, ri, hosp_d[i], cfr_j), args_daly))
  per_draw[[vac_name]][i, ] <-
    do.call(outcome_one, c(list(symp_vacc, inf_tot, target_pop_elig*cov_d[i], ri, hosp_d[i], cfr_j), args_daly))
  if (i %% 100 == 0) cat("  ", i, "/", N_DRAWS, "\n")
}

# ------------------------------------------------------------
# 5. Condition on an outbreak; averted + NNV over outbreak draws (paired per draw)
# ------------------------------------------------------------
outbreak <- which(attack_base > OUTBREAK_ATTACK_THRESH)
p_outbreak <- length(outbreak) / N_DRAWS
cat(sprintf("\nP(outbreak takes off | high-R0 prior) = %.1f%%  (%d / %d draws, attack > %.1f%%)\n",
            100*p_outbreak, length(outbreak), N_DRAWS, OUTBREAK_ATTACK_THRESH))
if (length(outbreak) < 20)
  cat("  NOTE: few outbreak draws -- conditional estimates are noisy; consider a wider high prior or lower threshold.\n")

base_pd <- per_draw[["No vaccine (baseline)"]]
vac_pd  <- per_draw[[vac_name]]
averted <- base_pd[, NNV_OUT, drop=FALSE] - vac_pd[, NNV_OUT, drop=FALSE]
nnv     <- vac_pd[, "doses"] / averted                     # doses recycled across the 4 cols
nnv[!is.finite(nnv) | nnv < 0] <- NA

# Conditional (outbreak-only) aggregation
q3  <- function(x) quantile(x, c(.5,.025,.975), na.rm = TRUE)
aggc <- function(mat) do.call(rbind, lapply(colnames(mat), function(o) {
  q <- q3(mat[outbreak, o]); data.frame(outcome=o, median=q[1], lo=q[2], hi=q[3], row.names=NULL) }))
agg_burden_cond  <- cbind(scenario = "baseline",           aggc(base_pd))
agg_averted_cond <- cbind(scenario = vac_name,             aggc(averted))
agg_nnv_cond     <- cbind(scenario = vac_name,             aggc(nnv))

# ------------------------------------------------------------
# 6. Console summary (conditional on taking off)
# ------------------------------------------------------------
fmtq_ <- function(v, d=0) fmtq(v[outbreak], d)
cat("\n=== Baseline burden | outbreak (median, 95% UI over outbreak draws) ===\n")
for (o in c("infections","symptomatic","hospitalisations","deaths","daly"))
  cat(sprintf("  %-16s %s\n", o, fmtq_(base_pd[, o], if (o=="deaths") 1 else 0)))
cat("\n=== Pre-outbreak disease-blocking vaccine, AVERTED | outbreak ===\n")
for (o in NNV_OUT) cat(sprintf("  %-16s %s\n", o, fmtq_(averted[, o], if (o=="deaths") 1 else 0)))
cat("\n=== NNV (doses per burden averted) | outbreak ===\n")
for (o in NNV_OUT) cat(sprintf("  %-16s %s\n", o, fmtq(nnv[outbreak, o], 0)))

# ------------------------------------------------------------
# 6b. Epidemic curve (CHIKV-style): symptomatic cases, baseline vs pre-outbreak
#     disease-blocking, plotted CONDITIONAL on the outbreak taking off. If too few
#     outbreak draws (e.g. the low-R0 scenario), fall back to ALL draws -- which then
#     correctly shows a flat ~no-outbreak curve rather than a spurious hump.
# ------------------------------------------------------------
use_cond <- length(outbreak) >= 20
draw_set <- if (use_cond) outbreak else seq_len(N_DRAWS)
bandq    <- function(M) apply(M[draw_set, , drop = FALSE], 2, quantile, c(.025,.5,.975), na.rm = TRUE)
bb <- bandq(wk_base); bv <- bandq(wk_vacc)
wk_num   <- function(idx) ifelse(idx <= 14, idx + 39, idx - 14)     # 2025-W40..W53 | 2026-W01..W38
tick_idx <- c(5, 10, 15, 22, 34, 46); xt <- data.frame(i = tick_idx, w = wk_num(tick_idx))
red <- 100 * (base_pd[draw_set,"symptomatic"] - vac_pd[draw_set,"symptomatic"]) / base_pd[draw_set,"symptomatic"]
rq  <- quantile(red, c(.5,.025,.975), na.rm = TRUE)
lab <- sprintf("%% symptomatic reduction\nDisease-blocking: %.1f%% (%.1f-%.1f%%)", rq[1], rq[2], rq[3])
roll_end <- start_pre + max(1, round(1 / mean(del_d))) - 0.5       # ~vaccine rollout window
pdf_df <- data.frame(week = 1:T_weeks, b_lo=bb[1,], b_md=bb[2,], b_hi=bb[3,],
                     v_lo=bv[1,], v_md=bv[2,], v_hi=bv[3,])
ttl <- sprintf("MAYV symptomatic cases (2025-W40 - 2026-W38) | R0 '%s'%s", E$R0_scenario,
               if (use_cond) sprintf(", conditional on outbreak (%d/%d)", length(outbreak), N_DRAWS)
               else ", all draws (~no outbreak)")
p_epi <- ggplot(pdf_df, aes(week)) +
  annotate("rect", xmin = start_pre-0.5, xmax = roll_end, ymin = -Inf, ymax = Inf, fill = "#cdebc5", alpha = 0.5) +
  geom_vline(xintercept = E$year_break, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = b_lo, ymax = b_hi), fill = "grey55",  alpha = 0.25) +
  geom_ribbon(aes(ymin = v_lo, ymax = v_hi), fill = "#4292c6", alpha = 0.22) +
  geom_line(aes(y = b_md), colour = "grey30",  linewidth = 1) +
  geom_line(aes(y = v_md), colour = "#2171b5", linewidth = 1) +
  annotate("text", x = T_weeks, y = Inf, label = lab, hjust = 1, vjust = 1.2, size = 3.1) +
  scale_x_continuous(breaks = xt$i, labels = xt$w) +
  labs(x = "Week", y = "Predicted symptomatic cases", title = ttl) +
  theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
                       panel.grid.minor = element_blank())
epi_fn <- sprintf("MAYV_ca_symptomatic_%s.png", E$R0_scenario)
ggsave(epi_fn, p_epi, width = 8, height = 4.5, dpi = 120)
cat(sprintf("Saved epidemic-curve plot: %s (grey = no vaccine, blue = disease-blocking)\n", epi_fn))

# reusable BEFORE-VACCINE (baseline-only) curve: median + 95% band, no vaccine arm
base_curve_plot <- function(M, draws, ytitle, ttl, fn) {
  bq <- apply(M[draws, , drop = FALSE], 2, quantile, c(.025,.5,.975), na.rm = TRUE)
  dfp <- data.frame(week = 1:T_weeks, lo = bq[1,], md = bq[2,], hi = bq[3,])
  p <- ggplot(dfp, aes(week)) +
    geom_vline(xintercept = E$year_break, linetype = "dashed", colour = "grey55") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey55", alpha = 0.28) +
    geom_line(aes(y = md), colour = "grey20", linewidth = 1) +
    scale_x_continuous(breaks = xt$i, labels = xt$w) +
    labs(x = "Week", y = ytitle, title = ttl) +
    theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
                         panel.grid.minor = element_blank())
  ggsave(fn, p, width = 8, height = 4.5, dpi = 120); cat("Saved baseline (no-vaccine) plot:", fn, "\n")
}
base_curve_plot(wk_base, draw_set, "Predicted symptomatic cases (no vaccine)",
  sprintf("MAYV symptomatic, NO vaccine (2025-W40 - 2026-W38) | R0 '%s'%s", E$R0_scenario,
          if (use_cond) sprintf(", conditional on outbreak (%d/%d)", length(outbreak), N_DRAWS) else ", all draws"),
  sprintf("MAYV_ca_baseline_%s.png", E$R0_scenario))

# ------------------------------------------------------------
# 6c. FIXED-R0 representative outbreak. Pin peak R0 = R0_FIX and REUSE every other draw
#     (immunity, gamma, sigma, rho, prop_symp, coverage, VE) -- no new LHS. The band is
#     then PARAMETER uncertainty only (R0 removed), so it reads like the CHIKV plot:
#     one clean "substantial urban outbreak" + the vaccine's effect on it, rather than
#     the R0-dominated fan of the propagated/conditional plot above.
# ------------------------------------------------------------
fixed_base_pd <- fixed_vac_pd <- NULL; attack_fixed <- NULL
if (!is.na(R0_FIX)) {
  wkb <- wkv <- wk_inf_f <- wk_rep_f <- matrix(NA_real_, N_DRAWS, T_weeks)
  fixed_base_pd <- fixed_vac_pd <- matrix(NA_real_, N_DRAWS, length(OUTCOMES), dimnames = list(NULL, OUTCOMES))
  attack_fixed  <- numeric(N_DRAWS)
  for (i in 1:N_DRAWS) {
    gi<-E$gamma[i]; si<-E$sigma[i]; ri<-E$rho[i]; psi<-E$prop_symp[i]; immi<-E$immune_frac[i]
    Rimm <- rep(immi, A); sus <- N*(1-Rimm); I0i <- I0_total * sus/sum(sus)
    run <- seirv_vaccinated_MAYV(T_weeks, A, N, Rimm, I0i, R0_FIX * gi * season, si, gi, ri,
             target_age, cov_d[i], del_d[i], start_pre, VE_inf = 0, VE_block = veb_d[i],
             immun_delay = immun_delay, prop_symp = psi, E0 = E0, seed_week = seed_week)
    ninf <- run$new_infections; covf <- run$coverage_frac
    wkb[i, ] <- colSums(psi * ninf); wkv[i, ] <- colSums(psi * ninf * (1 - veb_d[i]*covf))
    wk_inf_f[i, ] <- colSums(ninf)          # weekly TRUE infections (no vaccine)
    wk_rep_f[i, ] <- ri * wkb[i, ]          # weekly REPORTED = rho * symptomatic (no vaccine)
    # full burden outcomes at fixed R0 (every fixed-R0 draw is an outbreak, so no conditioning)
    symp_b <- rowSums((psi*ninf)[, EVAL_WIN, drop=FALSE])
    symp_v <- rowSums((psi*ninf*(1 - veb_d[i]*covf))[, EVAL_WIN, drop=FALSE])
    inf_t  <- sum(ninf[, EVAL_WIN, drop=FALSE]); attack_fixed[i] <- 100*inf_t/sum(sus)
    cfr_j  <- (hosp_d[i]*cfrH_d[i,] + (1-hosp_d[i])*cfrN_d[i,])[age_to_band]
    ad <- list(le_d[i,], dwMM_d[i], dwSV_d[i], dwCH_d[i], duMM_d[i], duSV_d[i], duCH_d[i],
               p14y_d[i], p14o_d[i], p90y_d[i], p90o_d[i])
    fixed_base_pd[i,] <- do.call(outcome_one, c(list(symp_b, inf_t, 0, ri, hosp_d[i], cfr_j), ad))
    fixed_vac_pd[i,]  <- do.call(outcome_one, c(list(symp_v, inf_t, target_pop_elig*cov_d[i], ri, hosp_d[i], cfr_j), ad))
  }
  # BEFORE-VACCINE plots at fixed R0: (i) baseline symptomatic, (ii) true infections vs reported
  base_curve_plot(wkb, 1:N_DRAWS, "Predicted symptomatic cases (no vaccine)",
    sprintf("MAYV symptomatic, NO vaccine | representative outbreak, fixed R0 = %.1f", R0_FIX),
    sprintf("MAYV_ca_baseline_fixedR0_%.1f.png", R0_FIX))
  ibq <- apply(wk_inf_f, 2, quantile, c(.025,.5,.975)); rbq <- apply(wk_rep_f, 2, quantile, c(.025,.5,.975))
  dfir <- data.frame(week = 1:T_weeks, i_lo=ibq[1,], i_md=ibq[2,], i_hi=ibq[3,],
                     r_lo=rbq[1,], r_md=rbq[2,], r_hi=rbq[3,])
  p_ir <- ggplot(dfir, aes(week)) +
    geom_vline(xintercept = E$year_break, linetype = "dashed", colour = "grey55") +
    geom_ribbon(aes(ymin=i_lo, ymax=i_hi, fill="True infections"), alpha=0.25) +
    geom_ribbon(aes(ymin=r_lo, ymax=r_hi, fill="Reported"), alpha=0.30) +
    geom_line(aes(y=i_md, colour="True infections"), linewidth=1) +
    geom_line(aes(y=r_md, colour="Reported"), linewidth=1) +
    scale_colour_manual(name=NULL, values=c("True infections"="#3182bd", "Reported"="#d6604d")) +
    scale_fill_manual(name=NULL, values=c("True infections"="#a8d1e7", "Reported"="#f4a582")) +
    scale_x_continuous(breaks=xt$i, labels=xt$w) +
    labs(x="Week", y="Weekly cases (no vaccine)",
         title=sprintf("MAYV true infections vs reported, NO vaccine | fixed R0 = %.1f", R0_FIX)) +
    theme_bw(12) + theme(plot.title=element_text(face="bold", hjust=0.5, size=10.5),
                         legend.position="inside", legend.position.inside=c(0.98,0.98),
                         legend.justification=c(1,1), panel.grid.minor=element_blank())
  ggsave(sprintf("MAYV_ca_baseline_true_vs_reported_fixedR0_%.1f.png", R0_FIX), p_ir, width=8, height=4.5, dpi=120)
  cat(sprintf("Saved baseline true-vs-reported plot: MAYV_ca_baseline_true_vs_reported_fixedR0_%.1f.png\n", R0_FIX))
  bqf <- function(M) apply(M, 2, quantile, c(.025,.5,.975), na.rm = TRUE)
  bF <- bqf(wkb); vF <- bqf(wkv)
  redF <- 100 * (rowSums(wkb) - rowSums(wkv)) / rowSums(wkb); rqF <- quantile(redF, c(.5,.025,.975))
  labF <- sprintf("%% symptomatic reduction\nDisease-blocking: %.1f%% (%.1f-%.1f%%)", rqF[1], rqF[2], rqF[3])
  dfF <- data.frame(week = 1:T_weeks, b_lo=bF[1,], b_md=bF[2,], b_hi=bF[3,], v_lo=vF[1,], v_md=vF[2,], v_hi=vF[3,])
  p_fix <- ggplot(dfF, aes(week)) +
    annotate("rect", xmin = start_pre-0.5, xmax = roll_end, ymin = -Inf, ymax = Inf, fill = "#cdebc5", alpha = 0.5) +
    geom_vline(xintercept = E$year_break, linetype = "dashed", colour = "grey55") +
    geom_ribbon(aes(ymin = b_lo, ymax = b_hi), fill = "grey55",  alpha = 0.25) +
    geom_ribbon(aes(ymin = v_lo, ymax = v_hi), fill = "#4292c6", alpha = 0.22) +
    geom_line(aes(y = b_md), colour = "grey30",  linewidth = 1) +
    geom_line(aes(y = v_md), colour = "#2171b5", linewidth = 1) +
    annotate("text", x = T_weeks, y = Inf, label = labF, hjust = 1, vjust = 1.2, size = 3.1) +
    scale_x_continuous(breaks = xt$i, labels = xt$w) +
    labs(x = "Week", y = "Predicted symptomatic cases",
         title = sprintf("MAYV symptomatic cases (2025-W40 - 2026-W38) | representative outbreak, fixed R0 = %.1f", R0_FIX)) +
    theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
                         panel.grid.minor = element_blank())
  fix_fn <- sprintf("MAYV_ca_symptomatic_fixedR0_%.1f.png", R0_FIX)
  ggsave(fix_fn, p_fix, width = 8, height = 4.5, dpi = 120)
  cat(sprintf("Saved fixed-R0 plot: %s (median baseline symptomatic total = %.0f; band = parameter uncertainty only)\n",
              fix_fn, median(rowSums(wkb))))
}

# ------------------------------------------------------------
# 7. Save per-draw + aggregated + the outbreak index (severity-phase counts included)
# ------------------------------------------------------------
saveRDS(list(
  per_draw = per_draw, averted = averted, nnv = nnv,
  attack_base = attack_base, outbreak = outbreak, p_outbreak = p_outbreak,
  OUTBREAK_ATTACK_THRESH = OUTBREAK_ATTACK_THRESH,
  fixed_base_pd = fixed_base_pd, fixed_vac_pd = fixed_vac_pd, attack_fixed = attack_fixed,
  R0_FIX = R0_FIX, rho_draw = E$rho,
  agg_burden_cond = agg_burden_cond, agg_averted_cond = agg_averted_cond, agg_nnv_cond = agg_nnv_cond,
  scen_names = scen_names, vac_name = vac_name, OUTCOMES = OUTCOMES, NNV_OUT = NNV_OUT,
  N_DRAWS = N_DRAWS, PHASE_MODE = PHASE_MODE, R0_scenario = E$R0_scenario,
  target_pop_elig = target_pop_elig, cov_d = cov_d, veb_d = veb_d, EVAL_WIN = EVAL_WIN),
  "MAYV_ca_engine_results.rds")
cat("\nSaved MAYV_ca_engine_results.rds (per-draw + conditional aggregates; severity-phase counts included).\n")
