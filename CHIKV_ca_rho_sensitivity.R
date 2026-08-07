# ============================================================
# CHIKV_ca_rho_sensitivity.R
# ------------------------------------------------------------
# A VARIANT of CHIKV_ca_lhs.R, run as a sensitivity analysis. Two departures from the
# base case, both deliberate:
#
#   1. NO FLAT-HOLD. The base case estimates beta_t over weeks 1-49 and holds it flat to
#      week 52 (active_weeks = 49). Here the spline is estimated over all 52 weeks, so
#      the tail is free to move and any rise there is the fit's own doing.
#
#   2. THREE REPORTING-RATE SETTINGS, to see how much of the tail behaviour is really
#      about how many infections sit behind the 8,204 observed cases:
#        a. Beta(20,60)  -- the base case: median 25%, 95% UI 20.1-32.5%, SAMPLED
#        b. Lognormal matched to 7.02-38.87% (stated median 15.40%)  -- SAMPLED
#        d. Beta matched to the SAME two endpoints, for comparability with (a) -- SAMPLED
#        e. Beta matched to the MEDIAN and the upper bound instead              -- SAMPLED
#        c. 38.87% as a POINT ESTIMATE -- not sampled; everything else still is
#
# Everything else (FOI, gamma, sigma, prop_symp, the penalty, the peak weighting, the
# feasibility filter) matches CHIKV_ca_lhs.R. The base script is sourced with DEFS_ONLY
# and is NOT modified; the spline settings are overridden here.
#
# BIC uses the UNWEIGHTED negative binomial log-likelihood at the fitted parameters,
# as CHIKV_ca_lhs.R does -- the peak weighting is a fitting device, not a likelihood.
#
# Output: CHIKV_ca_rho_sensitivity.png and .xlsx
# ============================================================
DEFS_ONLY <- TRUE
suppressMessages({library(ggplot2); library(dplyr); library(writexl); library(splines)
                  source("CHIKV_ca_lhs.R")})

N_DRAWS <- if (exists("RHO_SENS_DRAWS")) RHO_SENS_DRAWS else 500

# ---- spline window --------------------------------------------------------------
# Default 52 = the whole window, removing the base case's flat-hold. Set
# RHO_SENS_ACTIVE_WEEKS <- 49 to keep the BASE-CASE spline, which isolates the reporting
# rate as the only thing that differs from the main analysis.
active_weeks <<- if (exists("RHO_SENS_ACTIVE_WEEKS")) RHO_SENS_ACTIVE_WEEKS else T_weeks
basis_full   <<- ns(seq_len(active_weeks), df = df_spline, intercept = TRUE)
gen_start    <<- c(cfl(log(c(seq(1.0, 2.2, length.out = peak_idx),
                             seq(2.2, 0.5, length.out = active_weeks - peak_idx)))), log(50))
cat(sprintf("Spline now active over all %d weeks (base case: 49).\n", active_weeks))

# ---- reporting-rate settings ------------------------------------------------------
# The municipal estimate 15.40% (7.02-38.87%) is right-skewed and NOT Beta-shaped: the
# best Beta misses the lower bound badly (3.9% against 7.0%), and the lower bound is what
# decides how many draws survive the attack-rate filter. A Lognormal matched to the two
# 95% endpoints reproduces them exactly. Its implied median is 16.5%, slightly above the
# stated 15.40%, because the quoted triple is not itself lognormal -- the interval is
# matched in preference to the median, since the tails drive feasibility here.
# This mirrors how MAYV_ca_lhs.R handles the skewed Lima seroprevalence prior.
low_ml <- (log(0.0702) + log(0.3887))/2
low_sl <- (log(0.3887) - log(0.0702))/(2*1.96)
cat(sprintf("Lognormal for 15.40%% (7.02-38.87%%): meanlog %.4f, sdlog %.4f -> quantiles %.4f / %.4f / %.4f\n",
            low_ml, low_sl, qlnorm(.025, low_ml, low_sl), qlnorm(.5, low_ml, low_sl),
            qlnorm(.975, low_ml, low_sl)))

SCEN <- list(
  list(key = "base",  lab = "Beta(20,60): 25% (20.1-32.5)",      sampled = TRUE,
       a = 20, b = 60),
  list(key = "low",   lab = "Lognormal: 7.02-38.87 (median 16.5)", sampled = TRUE,
       ml = low_ml, sl = low_sl),
  list(key = "point", lab = "38.87% point estimate (not sampled)", sampled = FALSE,
       fixed = 0.3887),
  # A Beta has two parameters, so it can match two of the three target quantiles but not
  # all three. Matching BOTH 95% endpoints (as the lognormal does) is the like-for-like
  # comparison; it costs a median of 19.7% against the stated 15.4%, worse than the
  # lognormal's 16.5%, because a Beta on (0,1) cannot be skewed enough for an interval
  # with hi/lo = 5.5. Included so the prior family can be chosen on evidence.
  list(key = "beta_low", lab = "Beta(4.69,18.17): 7.02-38.87 (median 19.7)", sampled = TRUE,
       a = 4.6898, b = 18.1660),
  # The other Beta the target admits: median 15.40% and upper 38.87% matched exactly, at
  # the cost of the lower bound (2.5th percentile 3.17% against the stated 7.02%). It
  # therefore places 13.6% of its mass below the lower bound the estimate reports, versus
  # 2.5% for the endpoint-matched version -- included so that cost can be seen in the
  # results rather than argued from the prior alone.
  list(key = "beta_low2", lab = "Beta(2.52,12.39): median 15.4 + upper 38.87", sampled = TRUE,
       a = 2.5181, b = 12.3863))

# Optionally run a subset, e.g. RHO_SENS_KEYS <- "beta_low" to add one scenario.
if (exists("RHO_SENS_KEYS"))
  SCEN <- Filter(function(x) x$key %in% RHO_SENS_KEYS, SCEN)

# ---- LHS over the inputs that stay sampled in every scenario ----------------------
# Same priors as CHIKV_ca_lhs.R. rho gets its own column, used only where sampled.
set.seed(2027)
lhs_col <- function(n) (sample.int(n) - runif(n)) / n
U <- sapply(1:5, function(j) lhs_col(N_DRAWS))
foi_lo <- 0.003; foi_hi <- 0.020
foi_ml <- (log(foi_lo) + log(foi_hi))/2; foi_sl <- (log(foi_hi) - log(foi_lo))/(2*1.96)
foi_d  <- qlnorm(U[,1], foi_ml, foi_sl)
gam_d  <- qnorm(U[,2], 0.54, 0.0714);  gam_d <- pmax(gam_d, 0.2)
lat_d  <- qnorm(U[,3], 0.60, 0.05);    sig_d <- 1/pmax(lat_d, 0.2)
psy_d  <- qbeta(U[,4], 35.84, 32.56)

peak_rise <- function(b) { pk <- which.max(b); tl <- b[pk:length(b)]
  tr <- which.min(tl); 100 * (max(tl[tr:length(tl)]) / min(tl) - 1) }

run_scen <- function(sc) {
  rho_d <- if (!sc$sampled) rep(sc$fixed, N_DRAWS)
           else if (!is.null(sc$a)) qbeta(U[,5], sc$a, sc$b)
           else qlnorm(U[,5], sc$ml, sc$sl)
  cat(sprintf("\n-- %s: %d draws --\n", sc$lab, N_DRAWS))
  out <- vector("list", N_DRAWS); bmat <- matrix(NA_real_, N_DRAWS, T_weeks)
  for (i in seq_len(N_DRAWS)) {
    f <- tryCatch(refit(foi_d[i], gam_d[i], sig_d[i], rho_d[i], psy_d[i], gen_start),
                  error = function(e) NULL)
    if (is.null(f)) next
    th <- 1 + exp(f$par[df_spline + 1])
    ll <- sum(dnbinom(observed_cases, mu = pmax(f$pred, 1e-6), size = th, log = TRUE))
    Rimm <- 1 - exp(-foi_d[i] * exposure_age)
    out[[i]] <- data.frame(
      draw = i, rho = rho_d[i],
      immune = sum(Rimm * N), susceptible = f$pool,
      infections = sum(f$inf), attack = f$attack, totrep = f$total,
      peak_beta = max(as.numeric(f$beta)), rise = peak_rise(as.numeric(f$beta)),
      bic = -2*ll + (df_spline + 1)*log(T_weeks))
    bmat[i, ] <- as.numeric(f$beta)
    if (i %% 100 == 0) cat(sprintf("   %d / %d\n", i, N_DRAWS))
  }
  d <- do.call(rbind, out)
  # same feasibility filter as the base script
  keep <- with(d, attack < 100 & abs(totrep - obs_total)/obs_total < 0.10)
  list(draws = d, keep_idx = d$draw[keep], n_keep = sum(keep), beta = bmat, lab = sc$lab)
}
RES <- lapply(SCEN, run_scen)
names(RES) <- vapply(SCEN, `[[`, character(1), "key")

# ---- summary table ----------------------------------------------------------------
q3  <- function(x) quantile(x, c(.5, .025, .975), na.rm = TRUE)
fmt <- function(x, d = 0) { q <- q3(x); sprintf("%s (%s-%s)",
        formatC(round(q[1],d), big.mark=",", format="f", digits=d),
        formatC(round(q[2],d), big.mark=",", format="f", digits=d),
        formatC(round(q[3],d), big.mark=",", format="f", digits=d)) }
tbl <- do.call(rbind, lapply(RES, function(r) {
  a <- r$draws                                   # all draws, before the feasibility filter
  d <- r$draws[r$draws$draw %in% r$keep_idx, ]   # retained
  data.frame(reporting_rate = r$lab,
             draws_retained = sprintf("%d / %d", r$n_keep, N_DRAWS),
             total_population = formatC(round(sum(N)), big.mark = ","),
             # Immunity is 1 - exp(-FOI * min(age,12)), a function of FOI ALONE -- rho
             # does not enter it. Over ALL draws it is therefore identical in every
             # scenario (the FOI column is shared). It differs across the RETAINED
             # subsets only because the feasibility filter selects on FOI: high-FOI draws
             # have a smaller susceptible pool, so they hit the attack-rate ceiling and
             # are dropped. Both are reported so that selection is visible, not hidden.
             immune_all_draws = fmt(a$immune),
             immune_retained = fmt(d$immune),
             susceptible_all_draws = fmt(a$susceptible),
             susceptible_retained = fmt(d$susceptible),
             true_infections = fmt(d$infections),
             attack_rate_pct = fmt(d$attack, 1),
             peak_beta = fmt(d$peak_beta, 3),
             rise_after_trough_pct = fmt(d$rise, 1),
             pct_draws_rise_over_5 = sprintf("%.1f%%", 100*mean(d$rise > 5)),
             BIC = fmt(d$bic, 1),
             stringsAsFactors = FALSE) }))
print(tbl, row.names = FALSE)

# ---- beta figure ------------------------------------------------------------------
band <- do.call(rbind, lapply(RES, function(r) {
  b <- r$beta[r$keep_idx, , drop = FALSE]
  q <- apply(b, 2, quantile, c(.025, .5, .975), na.rm = TRUE)
  data.frame(week = seq_len(T_weeks), lo = q[1,], med = q[2,], hi = q[3,],
             rho = r$lab, stringsAsFactors = FALSE) }))
notes <- data.frame(item = c("Departure 1", "Departure 2", "Held as in the base case",
                             "Feasibility filter", "BIC", "Beta fit", "Draws",
                             "Choice of prior family", "Immune and susceptible columns",
                             "Reading the rise column"),
  detail = c(
  sprintf("beta_t is spline-estimated over %d of the %d weeks (base case: 49, holding beta flat to week 52). At 49 the reporting rate is the ONLY departure from the main analysis; at 52 the flat-hold is removed as well.", active_weeks, T_weeks),
  "Three reporting-rate settings: Beta(20,60) sampled (base case); a Beta fitted to 15.40% (7.02-38.87%) sampled; and 38.87% as a fixed point estimate with everything else still sampled.",
  "FOI ~ Lognormal(95% 0.003-0.020), gamma ~ N(0.54, 0.0714), latent ~ N(0.60, 0.05), prop_symp ~ Beta(35.84, 32.56), ridge penalty on log beta, peak weighting 1 + 10*(y/max y).",
  "attack < 100% and predicted reported total within 10% of the observed 8,204, as in CHIKV_ca_lhs.R.",
  "-2*logL + k*log(52) with k = df_spline + 1 = 6, using the UNWEIGHTED negative binomial log-likelihood: the peak weighting is a fitting device, not a likelihood.",
  sprintf("Lognormal(meanlog %.4f, sdlog %.4f) matched to the 95%% endpoints 7.02-38.87%%; implied median 16.5%% vs the stated 15.40%%. A Beta was rejected: the best fit put the 2.5th percentile at 3.9%% against 7.0%%, and the lower tail is what decides feasibility.", low_ml, low_sl),
  sprintf("%d per scenario (the base pipeline uses 1000; reduced here because each draw is a full re-fit).", N_DRAWS),
  "Two Beta parameterisations of the municipal estimate are included because a Beta has two parameters and the target has three constraints: Beta(4.69,18.17) matches both 95% endpoints (median falls to 19.7%), Beta(2.52,12.39) matches the median and upper bound (2.5th percentile falls to 3.17%, so 13.6% of its mass sits below the stated 7.02% lower bound). The Lognormal matches both endpoints with a median of 16.5%.",
  "Prior immunity is 1 - exp(-FOI * min(age, 12)) and the susceptible pool is its complement, so BOTH depend on FOI alone -- the reporting rate does not enter either. Over all draws they are identical in every scenario, at 8,709 (3,478-20,841) immune. The retained columns differ only because the feasibility filter selects on FOI: draws with high FOI have a smaller susceptible pool, hit the attack-rate ceiling, and are dropped, so the retained subset is skewed towards lower immunity. Read the _retained columns as a diagnostic of the filter, not as an effect of rho.",
  "The MEDIAN rise after the post-peak trough is 0.0% in all three scenarios -- most draws give a monotonically falling tail even with the spline free over all 52 weeks. The rise is a TAIL phenomenon confined to low-rho draws, so pct_draws_rise_over_5 (the share of retained draws rising by more than 5%) is the informative statistic, not the median."),
  stringsAsFactors = FALSE)
# When only a subset of scenarios is run (RHO_SENS_KEYS), MERGE into the existing
# workbook rather than overwriting it, so a scenario can be added without re-fitting the
# others. The summary and beta_bands sheets are rebuilt from the union.
XL <- if (exists("RHO_SENS_XLSX")) RHO_SENS_XLSX else "CHIKV_ca_rho_sensitivity.xlsx"
new_sheets <- setNames(lapply(RES, `[[`, "draws"), paste0("draws_", names(RES)))
if (file.exists(XL) && exists("RHO_SENS_KEYS")) {
  old <- setNames(lapply(readxl::excel_sheets(XL),
                         function(z) as.data.frame(readxl::read_xlsx(XL, z))),
                  readxl::excel_sheets(XL))
  keep_old <- old[setdiff(grep("^draws_", names(old), value = TRUE), names(new_sheets))]
  sheets <- c(keep_old, new_sheets)
  tbl  <- rbind(old$summary[!(old$summary$reporting_rate %in% tbl$reporting_rate), , drop = FALSE], tbl)
  band <- rbind(old$beta_bands[!(old$beta_bands$rho %in% band$rho), , drop = FALSE], band)
} else sheets <- new_sheets
write_xlsx(c(list(summary = tbl, notes = notes, beta_bands = band), sheets), XL)

# Figure is drawn from the MERGED band, so a subset run still redraws every scenario.
tk  <- c(1, 10, 20, 30, 40, 52)
lab <- sprintf("%d-W%02d", caldas_obs$Year[tk], caldas_obs$week[tk])
# The workbook keeps all five scenarios; the FIGURE shows the three that are
# substantively distinct -- the base prior, one Beta parameterisation of the municipal
# estimate, and the point estimate. The other two low-rho priors sit on top of
# Beta(2.52,12.39) and only crowd the legend.
PLOT_KEYS <- c("Beta(20,60): 25% (20.1-32.5)",
               "Beta(2.52,12.39): median 15.4 + upper 38.87",
               "38.87% point estimate (not sampled)")
band <- band[band$rho %in% PLOT_KEYS, ]
band$rho <- factor(band$rho, levels = PLOT_KEYS)
p <- ggplot(band, aes(week, med, colour = rho, fill = rho)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .14, colour = NA) +
  geom_line(linewidth = .9) +
  scale_x_continuous(breaks = tk, labels = lab) +
  coord_cartesian(ylim = c(0, 2.6)) +
  labs(x = "Week", y = expression(paste("Weekly transmission rate ", beta[t])),
       colour = NULL, fill = NULL, title = "Fitted beta_t by reporting rate") +
  theme_bw(12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12))
ggsave("CHIKV_ca_rho_sensitivity.png", p, width = 9.5, height = 5.6, dpi = 150)


cat(sprintf("\nWrote %s and %s\n", XL, sub("\\.xlsx$", ".png", XL)))
