# ============================================================
# ca_common.R -- shared helpers for the Caldas Novas vaccination scripts
# ------------------------------------------------------------
# Single source of truth for code that was previously duplicated verbatim in
# CHIKV_ca_vacc.R and MAYV_ca_vacc.R. Source it once near the top of each:
#     source("ca_common.R")
# Provides:
#   week_to_index()       calendar (year, epi-week) -> within-window index
#   fmtq()                "median (2.5% - 97.5%)" formatter for MC draws
#   burden()              summarise an SEIR run -> infections/symp/hosp/deaths
#   load_burden_params()  read disease_progression.xlsx (Hyolim Table S4) and
#                         return the Beta(alpha,beta) severity parameters
# ============================================================
library(readxl)

# Calendar -> within-window index for the 2025-W23 -> 2026-W22 fit window:
#   2025-W23 = index 1  (2025 weeks: index = week - 22)
#   2026-W01 = index 31 (2026 weeks: index = 30 + week)
week_to_index <- function(year, week) ifelse(year == 2025, week - 22L, 30L + week)

# Summarise a Monte Carlo draw vector as "median (2.5% - 97.5%)" with d decimals.
fmtq <- function(v, d = 0) {
  q <- quantile(v, c(.5, .025, .975), na.rm = TRUE)
  f <- function(x) formatC(round(x, d), big.mark = ",", format = "f", digits = d)
  sprintf("%s (%s - %s)", f(q[1]), f(q[2]), f(q[3]))
}

# Burden extractor: infections, symptomatic, hospitalisations, deaths.
# hr (hosp rate) and cv (CFR-by-age) default to the means but accept per-draw values.
burden <- function(out, hr = hosp_rate, cv = cfr_vec) {
  symp_age <- rowSums(out$new_symptomatic)
  symp     <- sum(symp_age)
  c(infections       = sum(out$new_infections),
    symptomatic      = symp,
    hospitalisations = symp * hr,
    deaths           = sum(symp_age * cv))
}

# Read the disease-progression Beta(alpha, beta) hyperparameters (Hyolim Table S4,
# disease_progression.xlsx) and return the severity parameters used by burden() and
# the Monte Carlo draws. `A` is the number of model age groups (for the age->band
# length check). Returns a named list; the caller typically unpacks it into globals:
#     list2env(load_burden_params(A), globalenv())
# so downstream code can reference ps_a, hosp_a, cfr_vec, age_to_band, ... directly.
load_burden_params <- function(A,
                               dp_path  = "disease_progression.xlsx",
                               dp_sheet = "disease_progression") {
  dp <- read_excel(dp_path, sheet = dp_sheet)
  names(dp)[1:10] <- c("parameter", "group", "median", "ui_lo", "ui_hi",
                       "dist", "p1", "alpha", "p2", "beta")
  stopifnot(all(dp$dist == "Beta"))

  # Pull alpha/beta for rows matching an exact Parameter (+ optional Group regex),
  # ordered by the lower age bound when the Group is an "Age [lo, hi)" band.
  get_beta_ab <- function(param, group_regex = NULL) {
    d <- dp[dp$parameter == param, ]
    if (!is.null(group_regex)) d <- d[grepl(group_regex, d$group), ]
    lo <- suppressWarnings(as.numeric(sub(".*\\[\\s*([0-9]+).*", "\\1", d$group)))
    if (!all(is.na(lo))) d <- d[order(lo), ]
    list(a = d$alpha, b = d$beta)
  }

  # prop_symp (symptomatic among infections), "Overall"
  ps <- get_beta_ab("Probability of symptomatic cases among infections", "^Overall$")
  ps_a <- ps$a; ps_b <- ps$b

  # hospitalisation among symptomatic (single, not age-specific)
  hp <- get_beta_ab("Probability of hospitalisation among symptomatic cases")
  hosp_a <- hp$a; hosp_b <- hp$b
  hosp_rate <- hosp_a / (hosp_a + hosp_b)

  # case fatality by decadal band [0,10)..[80,90) (length 9), hospitalised vs not
  ch <- get_beta_ab("Probability of death among hospitalised cases")
  cn <- get_beta_ab("Probability of death among non-hospitalised cases")
  cfr_hosp_a <- ch$a; cfr_hosp_b <- ch$b
  cfr_nonh_a <- cn$a; cfr_nonh_b <- cn$b
  stopifnot(length(ps_a) == 1, length(hosp_a) == 1,
            length(cfr_hosp_a) == 9, length(cfr_nonh_a) == 9)

  cfr_h_mean <- cfr_hosp_a / (cfr_hosp_a + cfr_hosp_b)
  cfr_n_mean <- cfr_nonh_a / (cfr_nonh_a + cfr_nonh_b)
  cfr_band   <- hosp_rate * cfr_h_mean + (1 - hosp_rate) * cfr_n_mean  # death per symptomatic, by band

  age_to_band <- c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9)   # 12 model groups -> 9 decadal bands
  stopifnot(length(age_to_band) == A)
  cfr_vec <- cfr_band[age_to_band]                       # death per symptomatic case, length A

  list(ps_a = ps_a, ps_b = ps_b,
       hosp_a = hosp_a, hosp_b = hosp_b, hosp_rate = hosp_rate,
       cfr_hosp_a = cfr_hosp_a, cfr_hosp_b = cfr_hosp_b,
       cfr_nonh_a = cfr_nonh_a, cfr_nonh_b = cfr_nonh_b,
       cfr_h_mean = cfr_h_mean, cfr_n_mean = cfr_n_mean, cfr_band = cfr_band,
       age_to_band = age_to_band, cfr_vec = cfr_vec)
}
