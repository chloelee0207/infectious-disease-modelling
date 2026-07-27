# ============================================================
# CHIKV_ca_owsa.R -- Caldas Novas CHIKV one-way sensitivity analysis.
#
# DETERMINISTIC. Every input is held at its central value and ONE parameter at a time
# is moved to its lower and upper bound; the change in vaccine impact is the tornado.
# This is the deterministic complement to the probabilistic (LHS) analysis in
# CHIKV_ca_lhs.R / CHIKV_ca_engine.R, which vary everything at once.
#
# Two classes of parameter:
#   FIT parameters (FOI, rho) change prior immunity or the case scaling, so beta must
#     be RE-FITTED to the same 8,204 observed cases at each bound.
#   CAMPAIGN parameters (VE, coverage, delivery, delay, time-to-immunity) leave the
#     epidemic fit untouched; only the SEIRV re-runs.
#
# Outcomes, all accrued over the 52-week observed window (2025-W24 -> 2026-W22):
#   symptomatic cases averted, deaths averted, DALYs averted, direct medical cost
#   averted (2019 BRL), and % symptomatic reduction.
#
# Run order:  CHIKV_ca_lhs.R (once, for the ensemble) -> this file
# ============================================================
suppressMessages({library(dplyr); library(ggplot2); library(writexl)})

DEFS_ONLY <- TRUE                 # load the fit machinery without the 1000-draw run
source("CHIKV_ca_lhs.R")          # refit(), seir_baseline(), data, fixed choices
source("ca_common.R")             # seirv_vaccinated, load_burden_params, load_daly_params

invisible(list2env(load_burden_params(A), globalenv()))
dp <- load_daly_params()
obs_band_prop <- load_caldas_age_cases()$obs_band_prop
young_idx <- which(age_to_band <= 4); old_idx <- which(age_to_band >= 5)

# ------------------------------------------------------------
# 1. Central values and one-way bounds
# ------------------------------------------------------------
# Central = the same point estimate the LHS uses for its reference fit. Bounds are the
# 95% interval of each parameter's propagated distribution, except delay and
# time-to-immunity, which are varied by +/- 1 week (their sampled range).
BASE <- list(foi = 0.008, rho = 0.25, ve = 263/266, cov = 0.30,
             deliv = 0.10, delay = 2, immun = 2)

BOUNDS <- list(
  foi   = c(0.003, 0.020),                        # serocatalytic 95% UI
  rho   = unname(qbeta(c(.025, .975), 20, 60)),   # Beta(20,60)
  ve    = unname(qbeta(c(.025, .975), 264, 4)),   # Beta(264,4) from 263/266 seroprotected
  cov   = c(0.20, 0.40),
  deliv = c(0.09, 0.11),
  delay = c(1, 3),
  immun = c(1, 3))

PAR_LAB <- c(foi   = "FOI",
             rho   = "Reporting rate",
             ve    = "Vaccine efficacy",
             cov   = "Vaccine coverage",
             deliv = "Weekly delivery speed",
             delay = "Delay in deployment",
             immun = "Time to immunity")
REFIT_PARS <- c("foi", "rho")     # these change the epidemic fit

GAMMA <- 0.54; SIGMA <- 1/0.60; PSYMP <- prop_symp   # held fixed (see header of the workbook)
start_pre <- caldas_obs$week_index[caldas_obs$Year == 2025 & caldas_obs$week == 40]
target_age <- rep(0, A); target_age[c(4,5,6,7,8)] <- 1          # eligible 18-59
target_pop_elig <- sum(N[target_age == 1])
EVAL_WIN <- seq_len(T_weeks)

# ------------------------------------------------------------
# 2. Deterministic outcome extractor (mirrors CHIKV_ca_engine.R outcome_one at
#    central severity/DALY values; recovery funnel, each case counted once).
# ------------------------------------------------------------
mean_ab <- function(x) x$a / (x$a + x$b)
DW  <- c(mm = mean_ab(dp$dw_mm), sev = mean_ab(dp$dw_sev), chr = mean_ab(dp$dw_chr))
DUR <- c(mm = dp$du_mm$med, sev = dp$du_sev$med, sub = dp$du_sub$med, chr = dp$du_chr$med)
LE  <- dp$le$med
REC <- c(acy = mean_ab(dp$p14_y), aco = mean_ab(dp$p14_o),
         sby = mean_ab(dp$p90_y), sbo = mean_ab(dp$p90_o),
         chy = mean_ab(dp$p6_y) + mean_ab(dp$p12_y) + mean_ab(dp$p30_y),
         cho = mean_ab(dp$p6_o) + mean_ab(dp$p12_o) + mean_ab(dp$p30_o))
CFR <- cfr_vec                                            # death per symptomatic, by age

# Direct medical cost per case at median unit costs (Goncalves 2024; see CHIKV_ca_costs.R)
CO <- as.data.frame(readxl::read_excel("costs.xlsx", sheet = "costs"))
names(CO)[1:5] <- c("parameter", "median", "lo", "hi", "dist")
um <- function(nm) CO$median[match(nm, CO$parameter)]
appt <- um("Medical appointment cost (2019 R$)")
PC <- c(
  acute = 2*appt + um("Minor/moderate pain as share of acute phase cases") *
            (um("Dipyrone - cost per case (2019 R$)") + um("Acetaminophen - cost per case (2019 R$)")) +
          (1 - um("Minor/moderate pain as share of acute phase cases")) *
            ((um("Dipyrone - cost per case (2019 R$)") + um("Acetaminophen - cost per case (2019 R$)"))/2 +
             (um("Tramadol - cost per case (2019 R$)") + um("Codeine - cost per case (2019 R$)") +
              um("Oxycodone - cost per case (2019 R$)"))/3),
  sub   = 3*appt + um("Arthritis as share of sub-acute phase cases")*um("Prednisone - cost per case (2019 R$)") +
          (1 - um("Arthritis as share of sub-acute phase cases")) *
            (((um("Amitriptyline - cost per case (2019 R$)") + um("Gabapentin - cost per case (2019 R$)"))/2 +
               um("Ibuprofen - cost per case (2019 R$)"))/2),
  chr   = 3*appt + um("Mild illness as share of chronic phase cases")*um("Hydroxychloroquine - cost per case (2019 R$)") +
          (1 - um("Mild illness as share of chronic phase cases")) *
            (um("Methotrexate - cost per case (2019 R$)") + um("Folic acid - cost per case (2019 R$)")),
  hosp  = um("Public share of inpatient admissions")*um("Average inpatient stay cost (2019 R$) - public sector") +
          (1 - um("Public share of inpatient admissions"))*um("Average inpatient stay cost (2019 R$) - private sector"))

outcomes_det <- function(out, age_weight) {
  symp_age <- rowSums(out$new_symptomatic[, EVAL_WIN, drop = FALSE])
  symp_w   <- symp_age * age_weight
  symp_dw  <- if (sum(symp_w) > 0) symp_w * (sum(symp_age)/sum(symp_w)) else symp_age
  st <- sum(symp_dw); sy <- sum(symp_dw[young_idx]); so <- sum(symp_dw[old_idx])

  n_hosp <- st * hosp_rate; n_nonhosp <- st - n_hosp
  syt <- REC["acy"] + REC["sby"] + REC["chy"]; sot <- REC["aco"] + REC["sbo"] + REC["cho"]
  n_acute <- sy*REC["acy"]/syt + so*REC["aco"]/sot     # resolved <= 14d
  n_sub   <- sy*REC["sby"]/syt + so*REC["sbo"]/sot     # resolved 14d - 3m
  n_chr   <- sy*REC["chy"]/syt + so*REC["cho"]/sot     # still ill > 3m

  yld <- unname((n_nonhosp*DW["mm"]*DUR["mm"] + n_hosp*DW["sev"]*DUR["sev"]) * (n_acute/st) +
                n_sub*DW["chr"]*DUR["sub"] + n_chr*DW["chr"]*DUR["chr"])
  deaths <- sum(symp_dw * CFR)
  yll    <- sum(symp_dw * CFR * LE[age_to_band])
  # cost: outpatient on NON-hospitalised (nested phases) + inpatient per admission
  cost <- n_hosp*PC["hosp"] + n_nonhosp*PC["acute"] +
          n_nonhosp*((n_sub + n_chr)/st)*PC["sub"] + n_nonhosp*(n_chr/st)*PC["chr"]
  c(infections = sum(out$new_infections[, EVAL_WIN]), symptomatic = st,
    hospitalisations = n_hosp, deaths = deaths, daly = yld + yll, cost = unname(cost))
}

# ------------------------------------------------------------
# 3. One full scenario: re-fit if needed, then baseline + both vaccine arms
# ------------------------------------------------------------
fit_cache <- list()
get_fit <- function(foi, rho) {
  k <- sprintf("%.6f_%.6f", foi, rho)
  if (is.null(fit_cache[[k]])) {
    f <- refit(foi, GAMMA, SIGMA, rho, PSYMP, gen_start)
    if (is.null(f)) stop("re-fit failed at FOI=", foi, " rho=", rho)
    fit_cache[[k]] <<- f
  }
  fit_cache[[k]]
}

run_scenario <- function(p) {
  f    <- get_fit(p$foi, p$rho)
  Rimm <- 1 - exp(-p$foi * exposure_age)
  sfrac <- (N*(1-Rimm))/sum(N*(1-Rimm))
  I0i  <- round(((week_1_cases/p$rho/PSYMP)/GAMMA) * sfrac)
  st   <- min(start_pre + p$delay, T_weeks)
  sim  <- function(cov, vi, vb)
    seirv_vaccinated(T_weeks, A, N, Rimm, I0i, E0, f$beta, SIGMA, GAMMA, p$rho,
                     target_age, cov, p$deliv, st, vi, vb, p$immun, prop_symp = PSYMP)
  base <- sim(0, 0, 0)
  aw   <- compute_age_weight(rowSums(base$new_infections), obs_band_prop, age_to_band)
  ob   <- outcomes_det(base, aw)
  list(base = ob,
       disb = ob - outcomes_det(sim(p$cov, 0,    p$ve), aw),   # averted, disease-blocking
       both = ob - outcomes_det(sim(p$cov, p$ve, p$ve), aw))   # averted, + infection-blocking
}

# ------------------------------------------------------------
# 4. Run base case, then each parameter at both bounds
# ------------------------------------------------------------
cat("OWSA: base case + ", 2*length(BOUNDS), " one-way runs...\n", sep = "")
base_run <- run_scenario(BASE)
rows <- list()
for (nm in names(BOUNDS)) {
  for (side in c("lower", "upper")) {
    p <- BASE; p[[nm]] <- BOUNDS[[nm]][if (side == "lower") 1 else 2]
    r <- run_scenario(p)
    for (arm in c("disb", "both"))
      rows[[length(rows)+1]] <- data.frame(
        parameter = unname(PAR_LAB[nm]), par_key = nm, bound = side,
        value = p[[nm]],
        arm = if (arm == "disb") "Disease-blocking" else "Disease + infection blocking",
        symptomatic = r[[arm]]["symptomatic"], deaths = r[[arm]]["deaths"],
        daly = r[[arm]]["daly"], cost = r[[arm]]["cost"],
        pct_symp = 100 * r[[arm]]["symptomatic"] / r$base["symptomatic"],
        row.names = NULL)
  }
  cat("  ", nm, " done\n", sep = "")
}
owsa <- do.call(rbind, rows)

base_tbl <- do.call(rbind, lapply(c("disb","both"), function(arm) data.frame(
  arm = if (arm=="disb") "Disease-blocking" else "Disease + infection blocking",
  symptomatic = base_run[[arm]]["symptomatic"], deaths = base_run[[arm]]["deaths"],
  daly = base_run[[arm]]["daly"], cost = base_run[[arm]]["cost"],
  pct_symp = 100*base_run[[arm]]["symptomatic"]/base_run$base["symptomatic"], row.names = NULL)))

# ------------------------------------------------------------
# 5. Tornado: symptomatic cases averted, ordered by swing
# ------------------------------------------------------------
tornado <- function(outcome, xlab, file) {
  arm_lv <- c("Disease-blocking", "Disease + infection blocking")
  d <- owsa; d$val <- d[[outcome]]
  b <- setNames(base_tbl[[outcome]], base_tbl$arm)
  d$base <- b[d$arm]; d$arm <- factor(d$arm, levels = arm_lv)
  sw <- d |> group_by(arm, parameter) |> summarise(swing = diff(range(val)), .groups = "drop") |>
        group_by(parameter) |> summarise(swing = max(swing), .groups = "drop")
  d$parameter <- factor(d$parameter, levels = sw$parameter[order(sw$swing)])
  bl <- data.frame(arm = factor(names(b), levels = arm_lv), base = as.numeric(b))
  p <- ggplot(d, aes(y = parameter)) +
    geom_vline(data = bl, aes(xintercept = base), linetype = "dashed", colour = "grey45") +
    geom_segment(aes(x = base, xend = val, yend = parameter, colour = bound),
                 linewidth = 6, alpha = .85) +
    facet_wrap(~ arm, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = c(lower = "#d6604d", upper = "#4393c3"),
                        labels = c(lower = "Lower", upper = "Upper"), name = "Bound") +
    scale_x_continuous(labels = scales::comma) +
    labs(x = xlab, y = NULL) +
    theme_bw(11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                         strip.text = element_text(face = "bold"))
  ggsave(file, p, width = 11, height = 4.6, dpi = 130); p
}
p_symp <- tornado("symptomatic", "Symptomatic cases averted", "CHIKV_ca_owsa_symptomatic.png")
p_daly <- tornado("daly",        "DALYs averted",                              "CHIKV_ca_owsa_daly.png")
p_cost <- tornado("cost",        "Direct medical cost averted (2019 BRL)",     "CHIKV_ca_owsa_cost.png")
print(p_symp)

# ------------------------------------------------------------
# 6. Export
# ------------------------------------------------------------
notes <- data.frame(item = c(
  "Analysis", "Window", "Vaccination", "Central values", "Bounds",
  "Fixed (not varied)", "Re-fitted parameters", "Outcomes"),
  detail = c(
  "Deterministic one-way sensitivity: all inputs at central values, one varied at a time.",
  sprintf("52 weeks, 2025-W24 -> 2026-W22 (indices %d-%d).", min(EVAL_WIN), max(EVAL_WIN)),
  sprintf("Pre-outbreak campaign, week index %d, eligible 18-59 (n = %s).",
          start_pre, format(round(target_pop_elig), big.mark = ",")),
  paste(sprintf("%s = %.4g", names(BASE), unlist(BASE)), collapse = "; "),
  paste(sprintf("%s [%.4g, %.4g]", names(BOUNDS), sapply(BOUNDS, `[`, 1), sapply(BOUNDS, `[`, 2)), collapse = "; "),
  paste("gamma, sigma, prop_symp. gamma is absorbed by the beta re-fit (R0 = beta/gamma",
        "against the same cases); sigma shifts peak timing not size; prop_symp cancels",
        "because the fit anchors on rho x prop_symp x infections = 8,204."),
  "FOI and rho change prior immunity / case scaling, so beta is re-fitted at each bound.",
  "Averted = baseline - scenario, both inside the window."),
  stringsAsFactors = FALSE)
write_xlsx(list(notes = notes, base_case = base_tbl, owsa = owsa), "CHIKV_ca_owsa.xlsx")
saveRDS(list(owsa = owsa, base = base_tbl, BASE = BASE, BOUNDS = BOUNDS,
             base_run = base_run), "CHIKV_ca_owsa.rds")

cat("\n=== Base case, pre-outbreak rollout ===\n"); print(base_tbl, row.names = FALSE)
cat("\n=== Swing in symptomatic cases averted (both-blocking) ===\n")
print(owsa |> filter(arm == "Disease + infection blocking") |>
        group_by(parameter) |> summarise(low = min(symptomatic), high = max(symptomatic),
                                         swing = high - low, .groups = "drop") |>
        arrange(desc(swing)) |> as.data.frame(), row.names = FALSE, digits = 5)
cat("\nWrote CHIKV_ca_owsa.xlsx, CHIKV_ca_owsa.rds and 3 tornado figures.\n")
