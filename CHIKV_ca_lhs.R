# ============================================================
# Caldas Novas CHIKV — full uncertainty propagation (self-contained).
#
# Propagates FIVE uncertain inputs through the age-structured SEIR fit:
#     gamma (recovery rate), sigma (latent rate), rho (reporting rate, Beta(20,60)),
#     prop_symp (symptomatic fraction), FOI (long-term force of infection ->
#     age-specific prior immunity). Note the fit only sees rho*prop_symp, so prop_symp
#     mainly widens the TRUE-infection iceberg (symptomatic/hosp/deaths anchor at obs/rho).
# For each Latin-Hypercube draw we RE-FIT the beta-spline + theta, so the final
# R0 / attack-rate / infection / beta bands carry all four sources of uncertainty
# instead of the fixed point values used in the point-estimate fit.
#
# This script copies the model machinery so it does not depend on the state of
# CHIKV_ca_pre_vacc_optim.R. Modelling choices are set explicitly below.
# ============================================================
setwd("/Users/chloelee/Documents/R/summer_project")
suppressMessages({library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(splines)})

# ------------------------------------------------------------
# 1. SEIR simulator (weekly, age-structured, 7 daily sub-steps)
# ------------------------------------------------------------
seir_baseline <- function(T_weeks, A, N, R_init_prop, I0, base_beta, sigma, gamma, rho,
                          prop_symp = 0.5242478, sub_steps = 7, E0 = rep(0, A)) {
  pmax0 <- function(x) pmax(0, x); N_total <- sum(N); dt <- 1/sub_steps
  S <- E <- I <- R <- matrix(0, A, T_weeks)
  new_infections <- new_symptomatic <- matrix(0, A, T_weeks)
  S_now <- pmax0(N - I0 - E0 - R_init_prop*N); E_now <- E0; I_now <- I0; R_now <- R_init_prop*N
  for (t in 1:T_weeks) {
    new_I_week <- rep(0, A); beta_t <- base_beta[t]
    for (k in 1:sub_steps) {
      foi <- beta_t * sum(I_now)/N_total
      new_E <- foi*S_now*dt; new_I <- sigma*E_now*dt; new_R <- gamma*I_now*dt
      S_now <- pmax0(S_now - new_E); E_now <- pmax0(E_now + new_E - new_I)
      I_now <- pmax0(I_now + new_I - new_R); R_now <- pmax0(R_now + new_R)
      new_I_week <- new_I_week + new_I
    }
    S[,t]<-S_now; E[,t]<-E_now; I[,t]<-I_now; R[,t]<-R_now
    new_infections[,t] <- new_I_week; new_symptomatic[,t] <- prop_symp*new_I_week
  }
  list(new_infections=new_infections, new_symptomatic=new_symptomatic, new_reported=rho*new_symptomatic)
}

# beta_t = exp(spline) over the active window, held flat afterwards (uses globals)
make_beta_t <- function(coefs) {
  beta_active <- as.numeric(exp(basis_full %*% coefs))
  beta_full   <- rep(beta_active[active_weeks], T_weeks)
  beta_full[seq_len(active_weeks)] <- beta_active
  beta_full
}

# negative log-posterior (negbin likelihood, peak-weighted, + lognormal beta prior)
neg_log_lik <- function(params, observed, T_weeks, df, A, N, R_init_prop, I0,
                        sigma, gamma, prop_symp, rho, E0) {
  coefs <- params[1:df]; theta <- 1 + exp(params[df+1])
  beta_t <- make_beta_t(coefs)
  if (any(!is.finite(beta_t)) || any(beta_t>20) || any(beta_t<1e-6)) return(1e10)
  out <- tryCatch(seir_baseline(T_weeks,A,N,R_init_prop,I0,beta_t,sigma,gamma,rho,prop_symp,E0=E0),
                  error=function(e) NULL)
  if (is.null(out)) return(1e10)
  predicted <- pmax(colSums(out$new_reported), 1e-6)
  wts <- 1 + peak_emphasis * (observed/max(observed))
  ll  <- sum(wts * dnbinom(observed, mu=predicted, size=theta, log=TRUE))
  if (!is.finite(ll)) return(1e10)
  lp <- if (is.finite(prior_logsd)) sum(dnorm(log(beta_t[1:active_weeks]), prior_logmean, prior_logsd, log=TRUE)) else 0
  -(ll + lp)
}

# ------------------------------------------------------------
# 2. Data: age structure (with age midpoints for the FOI model) + weekly cases
# ------------------------------------------------------------
normalize <- function(x) tolower(gsub("[áàâã]","a",gsub("[éèê]","e",gsub("[íì]","i",
              gsub("[óòôõ]","o",gsub("[úù]","u",gsub("ç","c",x)))))))
age_df <- read_excel("population.xlsx", sheet="prop_immune")
age_df <- as.data.frame(age_df[normalize(age_df$municipality)==normalize("Caldas Novas"),])
age_mid <- age_df$age_midpoint
pop_2022_total <- 98622; pop_2025_total <- 106820
N <- age_df$pop_num * (pop_2025_total/pop_2022_total)      # grow 2022 -> 2025
A <- nrow(age_df)

# Observed weekly cases: canonical 52-week loader from ca_common.R (single source of
# truth, shared with the fit and the age-stratified script). Window 2025-W24 -> 2026-W22.
if (!exists("load_caldas_age_cases")) source("ca_common.R")
ca_cases       <- load_caldas_age_cases()
caldas_obs     <- ca_cases$caldas_obs
observed_cases <- ca_cases$observed_cases
T_weeks <- length(observed_cases)
stopifnot(T_weeks == 52, sum(observed_cases) == 8204)
weeks <- 1:T_weeks; obs_total <- sum(observed_cases); week_1_cases <- observed_cases[1]

# ------------------------------------------------------------
# 3. Fixed modelling choices (explicit, not inherited)
# ------------------------------------------------------------
prop_symp     <- 0.5242478
df_spline     <- 5
active_weeks  <- 49        # hold beta flat from 2026-W19 (index 49 in the 52-week window)
prior_logmean <- log(0.54); prior_logsd <- 0.70
peak_emphasis <- 10
E0            <- rep(0, A)
basis_full    <- ns(seq_len(active_weeks), df=df_spline, intercept=TRUE)
peak_idx      <- which.max(observed_cases)
x_ticks   <- caldas_obs |> filter((Year==2025 & week %in% c(30,40,50))|(Year==2026 & week %in% c(10,20)))
year_break<- mean(c(max(caldas_obs$week_index[caldas_obs$Year==2025]), min(caldas_obs$week_index[caldas_obs$Year==2026])))

# ------------------------------------------------------------
# 4. One re-fit at given (FOI, gamma, sigma, rho): recompute immunity -> refit beta+theta
# ------------------------------------------------------------
cfl <- function(v) coef(lm(v ~ basis_full - 1))
gen_start <- c(cfl(log(c(seq(1.0,2.2,length.out=peak_idx), seq(2.2,0.5,length.out=active_weeks-peak_idx)))), log(50))
refit <- function(foi, g, s, r, ps, start) {
  Rimm <- 1 - exp(-foi * age_mid)                          # age-specific prior immunity (catalytic model)
  pool <- sum(N * (1 - Rimm)); sfrac <- (N*(1-Rimm))/pool
  I0i  <- round(((week_1_cases/r/ps)/g) * sfrac)            # seed depends on rho, gamma, prop_symp, immunity
  nll  <- function(par) neg_log_lik(par, observed_cases, T_weeks, df_spline, A, N, Rimm,
                                    I0i, s, g, ps, r, E0)
  f <- tryCatch(optim(start, nll, method="BFGS", control=list(maxit=900)), error=function(e) NULL)
  if (is.null(f)) return(NULL)
  bt <- make_beta_t(f$par[1:df_spline])
  o  <- seir_baseline(T_weeks,A,N,Rimm,I0i,bt,s,g,r,ps,E0=E0)
  list(par=f$par, beta=bt, pred=colSums(o$new_reported), inf=colSums(o$new_infections),
       R0=bt/g, attack=100*sum(o$new_infections)/pool, total=sum(colSums(o$new_reported)),
       pool=pool, immune=100*sum(Rimm*N)/sum(N))
}

# ------------------------------------------------------------
# 5. Samplers: gamma/sigma/rho from model_calibration.xlsx, FOI from its 95% UI
# ------------------------------------------------------------
cal <- as.data.frame(read_excel("model_calibration.xlsx", sheet=1))
cal <- cal[cal$Group=="CHIKV" & !is.na(cal$Median),]
row_for <- function(k) cal[grepl(k, cal$Parameter, ignore.case=TRUE),][1,]
sd_of   <- function(r) (r[["95% UI upper"]]-r[["95% UI lower"]])/(2*1.96)
gr<-row_for("gamma"); sr<-row_for("sigma"); rr<-row_for("reporting")
g_m<-gr$Median; g_sd<-sd_of(gr)                 # gamma  (rate)
p_m<-sr$Median; p_sd<-sd_of(sr)                 # latent PERIOD -> sigma = 1/period
r_m<-rr$Median; r_sd<-sd_of(rr)
# rho ~ Beta(20, 60): Hyolim's stated generative prior (mean 0.25, 95% ~0.162-0.350),
# wider than the model_calibration.xlsx posterior. Point estimate uses the median 0.25.
rab <- c(a = 20, b = 60)
foi_med<-0.008; foi_lo<-0.003; foi_hi<-0.020                        # long-term average FOI (from serology)
foi_mlog<-log(foi_med); foi_slog<-(log(foi_hi)-log(foi_lo))/(2*1.96)# FOI ~ Lognormal
# prop_symp ~ Beta(35.84, 32.56): symptomatic fraction among infections (Hyolim Table S4 /
# disease_progression.xlsx, median 0.524, 95% 0.406-0.641). Point estimate uses the median.
ps_a <- 35.84; ps_b <- 32.56
cat(sprintf("Samplers: gamma~N(%.3f,%.3f) latent~N(%.3f,%.3f) rho~Beta(%.1f,%.1f) FOI~logN(med %.3f, 95%%[%.3f,%.3f]) prop_symp~Beta(%.2f,%.2f)\n",
            g_m,g_sd,p_m,p_sd,rab["a"],rab["b"],foi_med,foi_lo,foi_hi,ps_a,ps_b))

# ------------------------------------------------------------
# 6. Baseline at the median inputs (for the dashed reference)
# ------------------------------------------------------------
base <- refit(foi_med, 0.54, 1/0.60, 0.25, prop_symp, gen_start); warm <- base$par
pk   <- which.max(base$beta)
cat(sprintf("Baseline (median inputs): immune %.1f%% | R0 at peak %.2f | attack %.1f%% | total %.0f\n",
            base$immune, (base$beta/0.54)[pk], base$attack, base$total))

# ---- Point-estimate fit diagnostics (median inputs: rho 0.25, gamma 0.54, latent 0.60 wk) ----
base_theta <- 1 + exp(base$par[df_spline + 1])
base_ll    <- sum(dnbinom(observed_cases, mu = pmax(base$pred, 1e-6), size = base_theta, log = TRUE))
base_k     <- df_spline + 1                                   # spline coefs + log_theta (rho fixed)
cat("\n--- Baseline point-estimate fit diagnostics ---\n")
cat(sprintf("Best rho:   0.25 (fixed) | Best theta: %.2f\n", base_theta))
cat(sprintf("beta_t range: [%.3f, %.3f] | R0 = beta/gamma range: [%.2f, %.2f]\n",
            min(base$beta), max(base$beta), min(base$R0), max(base$R0)))
cat(sprintf("Predicted total reported: %.0f | Observed total reported: %d  (%+.1f%%)\n",
            base$total, obs_total, 100*(base$total - obs_total)/obs_total))
cat(sprintf("Predicted peak week: %d (%.0f cases) | Observed peak week: %d (%d cases)\n",
            which.max(base$pred), max(base$pred), which.max(observed_cases), max(observed_cases)))
cat(sprintf("Log-likelihood: %.2f | k: %d | n: %d | AIC: %.2f | BIC: %.2f\n\n",
            base_ll, base_k, T_weeks, -2*base_ll + 2*base_k, -2*base_ll + base_k*log(T_weeks)))

# ------------------------------------------------------------
# 7. Latin Hypercube (5 inputs), re-fit each; feasibility handled by filtering
# ------------------------------------------------------------
set.seed(2024); n <- 1000        # match the 1000 LHS draws used in the reference study
lhs_col <- function(n) (sample.int(n)-runif(n))/n
U   <- sapply(1:5, function(j) lhs_col(n))
foi <- qlnorm(U[,1], foi_mlog, foi_slog)
gam <- qnorm (U[,2], g_m, g_sd)
sig <- 1/qnorm(U[,3], p_m, p_sd)
rho <- qbeta (U[,4], rab["a"], rab["b"])
psy <- qbeta (U[,5], ps_a, ps_b)                 # prop_symp per draw

beta_mat <- pred_mat <- inf_mat <- r0_mat <- matrix(NA_real_, n, T_weeks)
R0peak <- attack <- totrep <- immune <- loglik <- rep(NA_real_, n)
cat("Re-fitting", n, "LHS draws (5 inputs)...\n")
for (i in 1:n) {
  d <- tryCatch(refit(foi[i], gam[i], sig[i], rho[i], psy[i], warm), error=function(e) NULL)
  if (is.null(d)) next
  beta_mat[i,]<-d$beta; pred_mat[i,]<-d$pred; inf_mat[i,]<-d$inf; r0_mat[i,]<-d$R0
  R0peak[i]<-d$R0[pk]; attack[i]<-d$attack; totrep[i]<-d$total; immune[i]<-d$immune
  loglik[i]<-sum(dnbinom(observed_cases, mu=pmax(d$pred,1e-6), size=1+exp(d$par[df_spline+1]), log=TRUE))
  if (i %% 50 == 0) cat("  ", i, "/", n, "\n")
}
# Keep feasible + converged draws (attack < 95% and total within 10% of observed).
ok <- which(!is.na(totrep) & attack < 95 & abs(totrep-obs_total)/obs_total < 0.10)
cat(sprintf("Kept %d / %d draws (dropped %d infeasible/non-converged, mostly the high-FOI+low-rho corner the data rule out).\n",
            length(ok), n, n-length(ok)))

# ------------------------------------------------------------
# 8. Summaries (baseline point vs propagated median [95% UI]) + plots
# ------------------------------------------------------------
q3   <- function(x) quantile(x, c(.5,.025,.975), na.rm=TRUE)
band <- function(M) apply(M[ok,,drop=FALSE], 2, quantile, c(.025,.5,.975), na.rm=TRUE)
# Both totals come FROM THE MODEL (per draw), so they are internally consistent:
#   reported = sum(new_reported) = rho * prop_symp * sum(new_infections) = true infections scaled.
rep_tot  <- totrep[ok]                 # model reported cases  (sum new_reported per draw)
true_inf <- rowSums(inf_mat)[ok]       # model true infections (sum new_infections per draw)
cat("\n=========== PROPAGATED RESULTS (", length(ok), " feasible draws, 4 inputs) ===========\n", sep="")
cat(sprintf("Prior immunity:       baseline %.1f%%  -> propagated %.1f%% [%.1f%%, %.1f%%]\n",
            base$immune, q3(immune[ok])[1], q3(immune[ok])[2], q3(immune[ok])[3]))
cat(sprintf("R0 at peak (%s): baseline %.2f -> propagated %.2f [%.2f, %.2f]\n",
            caldas_obs$week_label[pk], (base$beta/0.54)[pk], q3(R0peak[ok])[1], q3(R0peak[ok])[2], q3(R0peak[ok])[3]))
cat(sprintf("Attack rate:          baseline %.1f%%  -> propagated %.1f%% [%.1f%%, %.1f%%]\n",
            base$attack, q3(attack[ok])[1], q3(attack[ok])[2], q3(attack[ok])[3]))
cat(sprintf("Reported cases total (model):  median %s [%s, %s]   (observed %s)\n",
            format(round(q3(rep_tot)[1]),big.mark=","), format(round(q3(rep_tot)[2]),big.mark=","),
            format(round(q3(rep_tot)[3]),big.mark=","), format(obs_total,big.mark=",")))
cat(sprintf("True infections total (model): median %s [%s, %s]\n",
            format(round(q3(true_inf)[1]),big.mark=","), format(round(q3(true_inf)[2]),big.mark=","),
            format(round(q3(true_inf)[3]),big.mark=",")))
aic_d <- -2*loglik + 2*(df_spline+1); bic_d <- -2*loglik + (df_spline+1)*log(T_weeks)
cat(sprintf("Goodness of fit across draws: logLik %.1f [%.1f, %.1f] | AIC %.1f [%.1f, %.1f] | BIC %.1f [%.1f, %.1f]\n",
            q3(loglik[ok])[1],q3(loglik[ok])[2],q3(loglik[ok])[3],
            q3(aic_d[ok])[1],q3(aic_d[ok])[2],q3(aic_d[ok])[3],
            q3(bic_d[ok])[1],q3(bic_d[ok])[2],q3(bic_d[ok])[3]))

bb <- band(beta_mat); r0b <- band(r0_mat); ib <- band(inf_mat); pb <- band(pred_mat)
save_band <- function(file, dfp, ytitle, ttl, add_dashed=NULL){
  g <- ggplot(dfp) + geom_vline(xintercept=year_break, linetype="dashed", colour="grey60") +
    geom_ribbon(aes(week, ymin=lo, ymax=hi), fill="#a8d1e7", alpha=.5) +
    geom_line(aes(week, med), colour="#3182bd", linewidth=1) +
    scale_x_continuous(breaks=x_ticks$week_index, labels=x_ticks$week) +
    labs(x="Week", y=ytitle, title=ttl) + theme_bw(12)
  if(!is.null(add_dashed)) g <- g + geom_line(data=add_dashed, aes(week,y), colour="#d6604d", linewidth=.8, linetype="dashed")
  ggsave(file, g, width=7.5, height=4.4, dpi=110)
}
save_band("CHIKV_ca_prop_beta.png", data.frame(week=weeks, lo=bb[1,], med=bb[2,], hi=bb[3,]),
          expression(beta[t]), "Caldas Novas beta(t): propagated 95% band (all 4 inputs)",
          data.frame(week=weeks, y=base$beta))
save_band("CHIKV_ca_prop_R0.png", data.frame(week=weeks, lo=r0b[1,], med=r0b[2,], hi=r0b[3,]),
          "R0(t) = beta/gamma", "Caldas Novas R0(t): propagated 95% band vs baseline (dashed)",
          data.frame(week=weeks, y=base$beta/0.54))
ggsave("CHIKV_ca_prop_infections.png",
  ggplot(data.frame(week=weeks, lo=ib[1,], med=ib[2,], hi=ib[3,], rep=pb[2,], obs=observed_cases)) +
    geom_ribbon(aes(week, ymin=lo, ymax=hi), fill="#a8d1e7", alpha=.5) +
    geom_line(aes(week, med), colour="#3182bd", linewidth=1) +
    geom_line(aes(week, rep), colour="#d6604d", linewidth=.9) + geom_point(aes(week, obs), size=1) +
    scale_x_continuous(breaks=x_ticks$week_index, labels=x_ticks$week) +
    labs(x="Week", y="Weekly cases", title="Caldas Novas: true infections (band) vs reported (coral/dots)") +
    theme_bw(12), width=7.5, height=4.4, dpi=110)

write.csv(data.frame(draw=1:n, FOI=foi, gamma=gam, sigma=sig, rho=rho, prop_symp=psy,
                     immune_pct=immune, R0_peak=R0peak, attack_pct=attack, total_reported=totrep,
                     feasible=(seq_len(n) %in% ok)), "CHIKV_ca_lhs_draws.csv", row.names=FALSE)
cat("\nSaved CHIKV_ca_prop_beta.png, CHIKV_ca_prop_R0.png, CHIKV_ca_prop_infections.png, CHIKV_ca_lhs_draws.csv\n")

# ------------------------------------------------------------
# 9. Export the FEASIBLE-draw ensemble for the standalone vaccine model.
#    CHIKV_ca_vacc.R loads this RDS instead of depending on CHIKV_ca_pre_vacc_optim.R
#    or re-running the 1000 refits. It carries all four calibration uncertainties
#    (FOI -> immunity, gamma, sigma, rho ~ Beta(20,60)); the vaccine MC iterates over
#    these draws and layers vaccine-parameter draws on top. The point estimate uses
#    the median inputs (rho 0.25). Immunity Rimm and the seed I0 are recomputed per
#    draw in the vaccine script from foi/gamma/rho (Rimm = 1 - exp(-foi * age_mid)).
ca_lhs_ensemble <- list(
  beta   = beta_mat[ok, , drop = FALSE],       # n_ok x T_weeks fitted beta_t per draw
  gamma  = gam[ok], sigma = sig[ok], rho = rho[ok], foi = foi[ok], prop_symp = psy[ok],
  base_beta = base$beta, base_rho = 0.25, base_foi = foi_med,
  base_gamma = 0.54, base_sigma = 1/0.60, base_prop_symp = prop_symp,
  N = N, A = A, age_mid = age_mid, age_df = age_df,
  week_1_cases = week_1_cases,
  T_weeks = T_weeks, observed_cases = observed_cases,
  caldas_obs = caldas_obs, weeks = weeks, x_ticks = x_ticks, year_break = year_break
)
saveRDS(ca_lhs_ensemble, "CHIKV_ca_lhs_ensemble.rds")
cat(sprintf("Saved CHIKV_ca_lhs_ensemble.rds (%d feasible draws) for the standalone vaccine model.\n",
            length(ok)))
