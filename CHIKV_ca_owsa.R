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
# 6. Scenario surface: campaign START WEEK x vaccine COVERAGE.
# The one-way tornado varies deployment DELAY by a week or two; this asks the bigger
# question the delay parameter cannot -- what does the campaign DATE cost? Ixchiq was
# intended as a pre-outbreak prophylactic (2025-W40); in Caldas Novas it was announced on
# 18 April 2026 = 2026-W16, after the modelled peak.
#
# Deterministic, like the rest of this script: 52 start weeks x 10 coverage levels x 2
# arms is a scenario surface, and re-propagating 1,040 cells would obscure the pattern
# without changing its shape. FOI and rho stay at their central values, so the fit is
# reused and the baseline epidemic is identical in every cell -- only the campaign moves.
# ------------------------------------------------------------
idx_of <- function(y, w) caldas_obs$week_index[caldas_obs$Year == y & caldas_obs$week == w]
WK_PRE <- start_pre                       # intended: pre-outbreak, 2025-W40
WK_ACT <- idx_of(2026, 16)                # actual: announced 18 April 2026

sw_fit <- get_fit(BASE$foi, BASE$rho)     # cached; the tornado has already built it
sw_Rimm  <- 1 - exp(-BASE$foi * exposure_age)
sw_sfrac <- (N*(1-sw_Rimm))/sum(N*(1-sw_Rimm))
sw_I0i   <- round(((week_1_cases/BASE$rho/PSYMP)/GAMMA) * sw_sfrac)
# The engine's age reweighting rescales to preserve the total, so for a TOTAL the raw
# sum is the same number and no weighting is needed here.
# The y-axis is the week the campaign is DECIDED. Dosing begins BASE$delay weeks later,
# exactly as run_scenario does (st <- start_pre + delay), so a cell here means the same
# thing as the corresponding scenario in the tornado and in the engine.
sw_sim <- function(cov, start, vi, vb)
  sum(seirv_vaccinated(T_weeks, A, N, sw_Rimm, sw_I0i, E0, sw_fit$beta, SIGMA, GAMMA,
                       BASE$rho, target_age, cov, BASE$deliv,
                       min(start + BASE$delay, T_weeks), vi, vb, BASE$immun,
                       prop_symp = PSYMP)$new_symptomatic[, EVAL_WIN])
SW_BASE <- sw_sim(0, WK_PRE, 0, 0)        # start week is irrelevant at coverage 0

SW_WEEKS <- seq_len(T_weeks); SW_COVS <- seq(0.10, 1.00, by = 0.10)
SW_ARMS  <- c("Disease-blocking", "Disease + infection blocking")
cat(sprintf("Scenario surface: %d start weeks x %d coverage levels x 2 arms...\n",
            length(SW_WEEKS), length(SW_COVS)))
surface <- do.call(rbind, lapply(SW_WEEKS, function(w) do.call(rbind, lapply(SW_COVS, function(cv)
  data.frame(week = w, coverage = cv, arm = factor(SW_ARMS, levels = SW_ARMS),
             pct_averted = 100 * (SW_BASE - c(sw_sim(cv, w, 0, BASE$ve),
                                              sw_sim(cv, w, BASE$ve, BASE$ve))) / SW_BASE,
             stringsAsFactors = FALSE)))))
surface$epi_week <- sprintf("%d-W%02d", caldas_obs$Year[surface$week],
                            caldas_obs$week[surface$week])

# fold loss between the intended and the actual campaign date, at each coverage level
intended_vs_actual <- surface |>
  filter(week %in% c(WK_PRE, WK_ACT)) |>
  mutate(timing = ifelse(week == WK_PRE, "intended_pre_outbreak", "actual_2026W16")) |>
  select(arm, coverage, timing, pct_averted) |>
  tidyr::pivot_wider(names_from = timing, values_from = pct_averted) |>
  mutate(fold_loss = intended_pre_outbreak / actual_2026W16) |> as.data.frame()


# The surface itself is deterministic, so its cells carry no interval. The two marked
# dates DO have one: both are engine scenarios, so their uncertainty is already
# propagated over the 1000-draw ensemble. The fold loss is formed per draw, keeping the
# two timings paired.
timing_propagated <- NULL
if (file.exists("CHIKV_ca_engine_results.rds")) {
  Ge <- readRDS("CHIKV_ca_engine_results.rds")
  be <- Ge$per_draw[["No vaccine (baseline)"]][, "symptomatic"]
  qf <- function(x, d) sprintf("%.*f (%.*f-%.*f)", d, median(x), d,
                               quantile(x, .025), d, quantile(x, .975))
  timing_propagated <- do.call(rbind, lapply(SW_ARMS, function(a) {
    pre <- 100*(be - Ge$per_draw[[paste0("pre-outbreak | ", a)]][, "symptomatic"])/be
    act <- 100*(be - Ge$per_draw[[paste0("actual rollout | ", a)]][, "symptomatic"])/be
    data.frame(arm = a, coverage = "sampled, Beta(30%, 20-40%)",
               intended_pre_outbreak_pct = qf(pre, 2),
               actual_2026W16_pct        = qf(act, 3),
               fold_loss                 = qf(pre/act, 0),
               stringsAsFactors = FALSE) }))
  cat("\n=== Propagated (1000 draws), the two marked campaign dates ===\n")
  print(timing_propagated, row.names = FALSE)
} else cat("NOTE: CHIKV_ca_engine_results.rds not found -- propagated timing sheet omitted.\n")

sw_ticks <- c(1, 10, 20, 30, 40, 52)
p_surface <- ggplot(surface, aes(100*coverage, week, fill = pct_averted)) +
  geom_raster(interpolate = TRUE) +
  geom_hline(yintercept = WK_PRE, colour = "grey15", linetype = "22", linewidth = .5) +
  geom_hline(yintercept = WK_ACT, colour = "grey15", linewidth = .5) +
  annotate("text", x = 12, y = WK_PRE + 2.4, hjust = 0, size = 4, colour = "grey15",
           label = "Modelled: pre-outbreak, 2025-W40") +
  annotate("text", x = 12, y = WK_ACT + 2.4, hjust = 0, size = 4, colour = "grey15",
           label = "Actual: announced 18 Apr 2026, 2026-W16") +
  facet_wrap(~ arm, nrow = 1) +
  scale_fill_gradientn(colours = c("#3b7fb6","#7fcdbb","#d9f0a3","#fee391","#fc8d59","#d73027"),
                       name = "% of symptomatic\ncases averted") +
  scale_x_continuous(breaks = seq(10, 100, 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = sw_ticks, expand = c(0, 0),
                     labels = sprintf("%d\n(%d-W%02d)", sw_ticks,
                                      caldas_obs$Year[sw_ticks], caldas_obs$week[sw_ticks])) +
  labs(x = "Vaccine coverage (% of eligible population aged 18-59)",
       y = "Start of vaccination campaign (week index)") +
  theme_bw(11) +
  theme(strip.text = element_text(face = "bold", size = 12),
        panel.grid = element_blank(), legend.position = "right")
ggsave("CHIKV_ca_timing_coverage.png", p_surface, width = 11, height = 5.6, dpi = 150)
ggsave("CHIKV_ca_timing_coverage.pdf", p_surface, width = 11, height = 5.6)
# 600-dpi companion for the manuscript (same size/layout as the 150-dpi version).
ggsave("CHIKV_ca_timing_coverage_600dpi.png", p_surface, width = 11, height = 5.6, dpi = 600)

# ------------------------------------------------------------
# 7. Export
# ------------------------------------------------------------
notes <- data.frame(item = c(
  "Analysis", "Window", "Vaccination", "Central values", "Bounds",
  "Fixed (not varied)", "Re-fitted parameters", "Outcomes",
  "Scenario surface", "Surface coverage axis", "Surface timing axis", "Intervals",
  "Interpreting the fold loss"),
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
  "Averted = baseline - scenario, both inside the window.",
  sprintf("Campaign start week x coverage, %d x %d x 2 arms, deterministic at the central set. Intended pre-outbreak date 2025-W40 (week %d) vs actual announcement 18 April 2026 = 2026-W16 (week %d). See the surface and intended_vs_actual sheets.", T_weeks, length(SW_COVS), WK_PRE, WK_ACT),
  "Coverage on the surface is of the ELIGIBLE 18-59 group (61.5% of the population), not of the whole population.",
  sprintf("The surface y-axis is the week the campaign is DECIDED; dosing begins %d weeks later (BASE$delay), as in run_scenario and the engine.", BASE$delay),
  "Surface cells are deterministic and carry NO interval. The two marked dates are engine scenarios, so their 95%% UIs are propagated over the 1000-draw ensemble -- see intended_vs_actual_propagated. Those medians differ slightly from the deterministic cells because the engine also samples coverage, efficacy and deployment delay.",
  "The fold loss is UNSTABLE and should not be quoted precisely: its denominator is ~0.03% of cases, so small changes in assumptions move it by a multiple. Adding the 2-week deployment delay alone takes it from ~193x to ~590x, because two weeks costs nothing pre-outbreak (99.5% of cases still ahead) but 68% of the remaining benefit at 2026-W16 (cases left after immunity fall 1,721 -> 806, and the 10-week rollout is truncated by the window end). Prefer reporting the two percentages, which are stable, and describing the loss as two to three orders of magnitude."),
  stringsAsFactors = FALSE)
write_xlsx(list(notes = notes, base_case = base_tbl, owsa = owsa,
                timing_coverage_surface = surface,
                intended_vs_actual = intended_vs_actual,
                intended_vs_actual_propagated = timing_propagated), "CHIKV_ca_owsa.xlsx")
saveRDS(list(owsa = owsa, base = base_tbl, BASE = BASE, BOUNDS = BOUNDS,
             base_run = base_run, surface = surface,
             intended_vs_actual = intended_vs_actual), "CHIKV_ca_owsa.rds")

cat("\n=== Base case, pre-outbreak rollout ===\n"); print(base_tbl, row.names = FALSE)
cat("\n=== Swing in symptomatic cases averted (both-blocking) ===\n")
print(owsa |> filter(arm == "Disease + infection blocking") |>
        group_by(parameter) |> summarise(low = min(symptomatic), high = max(symptomatic),
                                         swing = high - low, .groups = "drop") |>
        arrange(desc(swing)) |> as.data.frame(), row.names = FALSE, digits = 5)
cat("\n=== Campaign date: intended 2025-W40 vs actual 2026-W16 ===\n")
print(intended_vs_actual |> mutate(coverage = 100*coverage,
        across(where(is.numeric), \(x) round(x, 2))) |> as.data.frame(), row.names = FALSE)
cat("\nWrote CHIKV_ca_owsa.xlsx, CHIKV_ca_owsa.rds, 3 tornado figures",
    "and CHIKV_ca_timing_coverage.png/.pdf\n")
