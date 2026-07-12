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
suppressMessages({library(dplyr); library(tidyr)})

# Calendar -> within-window index for the 2025-W24 -> 2026-W22 fit window.
# The 2025 epi calendar has a Semana 53, so 2025 contributes 30 weeks (W24..W53):
#   2025-W24 = index 1  ; 2025-W53 = index 30   (2025 weeks: index = week - 23)
#   2026-W01 = index 31                          (2026 weeks: index = 30 + week)
week_to_index <- function(year, week) ifelse(year == 2025, week - 23L, 30L + week)

# Summarise a Monte Carlo draw vector as "median (2.5% - 97.5%)" with d decimals.
fmtq <- function(v, d = 0) {
  q <- quantile(v, c(.5, .025, .975), na.rm = TRUE)
  f <- function(x) formatC(round(x, d), big.mark = ",", format = "f", digits = d)
  sprintf("%s (%s - %s)", f(q[1]), f(q[2]), f(q[3]))
}

# Numeric median + 95% quantiles as c(median, 2.5%, 97.5%) for a MC draw vector.
# (fmtq() is the string formatter; qs() returns the raw numbers for further use.)
qs <- function(v) as.numeric(quantile(v, c(.5, .025, .975), na.rm = TRUE))

# Optional age re-weighting for the (age-sensitive) death calc. Defaults to 1 (no
# correction). The pipeline sets age_weight to observed_prop / model_prop by age so
# deaths use the DATA-based age distribution instead of the population-structure one.
# Guarded with exists() so RE-sourcing ca_common.R (e.g. from CHIKV_outputs.R) does
# not clobber an age_weight the engine already set.
if (!exists("age_weight")) age_weight <- 1

# Burden extractor: infections, symptomatic, hospitalisations, deaths.
# hr (hosp rate) and cv (CFR-by-age) default to the means but accept per-draw values.
# w re-weights the age distribution for DEATHS only, redistributing symptomatic cases
# across ages while PRESERVING the total (so infections/symptomatic totals and the
# single-rate hospitalisations are unchanged; only deaths move).
burden <- function(out, hr = hosp_rate, cv = cfr_vec, w = age_weight) {
  symp_age <- rowSums(out$new_symptomatic)
  symp     <- sum(symp_age)
  symp_w   <- symp_age * w
  sw       <- sum(symp_w)
  symp_dw  <- if (sw > 0) symp_w * (symp / sw) else symp_age   # renormalise: keep total
  c(infections       = sum(out$new_infections),
    symptomatic      = symp,
    hospitalisations = symp * hr,
    deaths           = sum(symp_dw * cv))
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

# ------------------------------------------------------------
# Canonical Caldas Novas age-stratified case loader (the "ca_combined" SINAN sheet).
# Single source of truth for the outbreak-window cases used by BOTH the fit
# (CHIKV_ca_pre_vacc_optim.R) and the age-stratified script (weekly_age_stratified.R),
# so neither has to depend on the older plain weekly_all series (weekly_case.R).
# Window 2025-W24 -> 2026-W22 = 52 weeks (2025 has an epi Semana 53). Missing zero-case
# weeks (2025-W33 & W40) are zero-filled on a canonical contiguous grid.
# Returns:
#   observed_cases  weekly totals (length 52), ordered by week_index
#   caldas_obs      week grid + totals (cols: Year, week, week_index, week_label, cases)
#   ca_age          long age x week table (cols: week_index, week_label, Year, week,
#                   age_group, cases)
#   age_totals      cases summed by ca_combined age group
#   obs_band_prop   observed case proportion across the 9 decadal CFR bands (length 9)
#   age_levels      the ca_combined age-group labels, in order
# ------------------------------------------------------------
load_caldas_age_cases <- function(path = "weekly_case.xlsx", sheet = "ca_combined") {
  age_levels <- c("<1 Ano", "1-4", "5-9", "10-14", "15-19",
                  "20-39", "40-59", "60-64", "65-69", "70-79", "80 e +")
  ca_long <- read_excel(path, sheet = sheet) |>
    rename(Year = Ano, semana = Semana) |>
    mutate(Year = as.integer(Year), week = as.integer(sub("Semana ", "", semana))) |>
    filter((Year == 2025 & week >= 24) | (Year == 2026 & week <= 22)) |>
    pivot_longer(all_of(age_levels), names_to = "age_group", values_to = "cases") |>
    mutate(cases = ifelse(is.na(cases), 0, cases))

  week_grid <- bind_rows(tibble(Year = 2025L, week = 24:53),
                         tibble(Year = 2026L, week = 1:22)) |>
    arrange(Year, week) |>
    mutate(week_index = row_number(),
           week_label = paste0(Year, "-W", sprintf("%02d", week)))

  caldas_obs <- week_grid |>
    left_join(ca_long |> group_by(Year, week) |>
                summarise(cases = sum(cases), .groups = "drop"),
              by = c("Year", "week")) |>
    mutate(cases = ifelse(is.na(cases), 0, cases))

  ca_age <- tidyr::expand_grid(week_grid, age_group = factor(age_levels, age_levels)) |>
    left_join(ca_long |> mutate(age_group = factor(age_group, age_levels)),
              by = c("Year", "week", "age_group")) |>
    mutate(cases = ifelse(is.na(cases), 0, cases)) |>
    dplyr::select(week_index, week_label, Year, week, age_group, cases)

  age_totals <- ca_age |> group_by(age_group) |>
    summarise(cases = sum(cases), .groups = "drop")

  # observed case proportion across the 9 decadal CFR bands; two-band groups split evenly
  group_to_bands <- list("<1 Ano" = 1, "1-4" = 1, "5-9" = 1, "10-14" = 2, "15-19" = 2,
                         "20-39" = c(3, 4), "40-59" = c(5, 6),
                         "60-64" = 7, "65-69" = 7, "70-79" = 8, "80 e +" = 9)
  obs_band <- numeric(9)
  for (g in age_levels) {
    b <- group_to_bands[[g]]
    n <- age_totals$cases[age_totals$age_group == g]
    obs_band[b] <- obs_band[b] + n / length(b)
  }

  list(observed_cases = caldas_obs$cases, caldas_obs = caldas_obs,
       ca_age = ca_age, age_totals = age_totals,
       obs_band_prop = obs_band / sum(obs_band), age_levels = age_levels)
}

# Age re-weighting vector w_a (length A) that maps the model's baseline infection age
# split to the OBSERVED case age split, per decadal band. Feed the model's baseline
# infections-by-age (length A) plus obs_band_prop + age_to_band from load_burden_params.
# Deaths computed with this w reproduce the observed age distribution (see burden()).
compute_age_weight <- function(inf_age_model, obs_band_prop, age_to_band) {
  mod_band <- as.numeric(tapply(inf_age_model, age_to_band, sum)[as.character(seq_along(obs_band_prop))])
  mod_band[is.na(mod_band)] <- 0
  mod_prop <- mod_band / sum(mod_band)
  w_band   <- ifelse(mod_prop > 0, obs_band_prop / mod_prop, 1)
  w_band[age_to_band]
}
