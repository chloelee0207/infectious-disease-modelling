# ============================================================
# CHIKV_outputs.R -- Caldas Novas presentation & export layer
# ------------------------------------------------------------
# The MODEL/SCENARIO/Monte-Carlo work lives in CHIKV_ca_vacc.R (the "engine").
# This script is PURELY presentation: it takes the engine's in-memory result
# objects and (1) draws the figures and (2) writes the Excel workbook.
#
# Run order in one R session:
#     source("ca_common.R")              # shared helpers
#     source("CHIKV_ca_pre_vacc_optim.R")# fit -> best_beta_t, N, ... (slow)
#     source("CHIKV_ca_vacc.R")          # scenarios + Monte Carlo (the engine)
#     source("CHIKV_outputs.R")          # <- this file: figures + Excel
#
# (The old coverage-based Sete-Lagoas burden Monte Carlo that used a global
#  `delay` has been retired -- that design does not match the Caldas timing
#  analysis, which passes delay = start_week per scenario. Its code remains in
#  git history if ever needed.)
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
if (file.exists("ca_common.R")) source("ca_common.R")   # fmtq(), etc.

# ------------------------------------------------------------
# 0. Guard: the engine must have been run first
# ------------------------------------------------------------
.needed <- c("scenarios", "scen_meta", "summary_tbl", "av_draws", "vac_names",
             "outcomes", "timings", "arms", "T_weeks", "caldas_obs",
             "base_tbl", "mc_tbl", "burden_mat",
             "wk_symp", "base_draws", "weekly_delivery_speed",
             "observed_cases", "best_rho", "ixchiq_efficacy",
             "total_coverage", "immun_delay", "rho_draws")
.missing <- .needed[!vapply(.needed, exists, logical(1))]
if (length(.missing)) {
  stop("CHIKV_outputs.R needs objects from the engine. Run CHIKV_ca_vacc.R first.\n",
       "  Missing: ", paste(.missing, collapse = ", "))
}

# ------------------------------------------------------------
# 1. Presentation styling (colours + shared calendar x-axis)
# ------------------------------------------------------------
timing_cols <- c("actual rollout" = "#f4a582",
                 "start of 2026"  = "#2166ac",
                 "pre-outbreak"   = "#b2182b")

x_ticks <- caldas_obs |>
  filter((Year == 2025 & week %in% c(30, 40, 50)) |
           (Year == 2026 & week %in% c(10, 20)))
year_break <- mean(c(max(caldas_obs$week_index[caldas_obs$Year == 2025]),
                     min(caldas_obs$week_index[caldas_obs$Year == 2026])))

# ============================================================
# 2. Figure (a): weekly TRUE infections, baseline + 3 timings
#    (infection-blocking arm; disease-blocking does not change infections)
# ============================================================
inf_keys <- c("No vaccine (baseline)",
              paste0(names(timings), " | Disease + infection blocking"))
inf_df <- do.call(rbind, lapply(inf_keys, function(nm) {
  lbl <- if (nm == "No vaccine (baseline)") nm else scen_meta$timing[scen_meta$name == nm]
  data.frame(week = 1:T_weeks, infections = colSums(scenarios[[nm]]$new_infections), scenario = lbl)
}))
inf_df$scenario <- factor(inf_df$scenario,
                          levels = c("No vaccine (baseline)", names(timings)))

p_inf <- ggplot(inf_df, aes(week, infections, colour = scenario)) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey60") +
  # geom_vline(xintercept = unlist(timings), linetype = "dotted",
  #            colour = timing_cols[names(timings)], linewidth = 0.5) +
  annotate("text", x = year_break, y = 0, label = "2026", angle = 90,
           vjust = -0.5, hjust = -11.5, fontface = "bold", size = 3.5) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("No vaccine (baseline)" = "grey40", timing_cols), name = NULL) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week", y = "Weekly infections (true)",
       title = "Caldas Novas CHIKV cases - disease + infection blocking") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())
print(p_inf)
ggsave("ca_vacc_infections.png", p_inf, width = 8, height = 4.5, dpi = 120)

# ============================================================
# 3. Figure (b): averted burden by timing x arm (point estimate)
# ============================================================
bar_df <- summary_tbl |>
  dplyr::select(timing, arm, Infections = inf_averted, Symptomatic = symp_averted,
                Hospitalisations = hosp_averted, Deaths = deaths_averted) |>
  pivot_longer(c(Infections, Symptomatic, Hospitalisations, Deaths),
               names_to = "outcome", values_to = "averted")
bar_df$timing  <- factor(bar_df$timing, levels = names(timings))
bar_df$arm     <- factor(bar_df$arm, levels = names(arms))   # disease-blocking first
bar_df$outcome <- factor(bar_df$outcome,
                         levels = c("Infections", "Symptomatic", "Hospitalisations", "Deaths"))

p_bar <- ggplot(bar_df, aes(timing, averted, fill = arm)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ outcome, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Disease-blocking" = "#9ecae1",
                               "Disease + infection blocking" = "#08519c"),
                    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Averted (vs no vaccine)",
       title = "Caldas Novas CHIKV cases - burden averted") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom", panel.grid.minor = element_blank())
# print(p_bar)
ggsave("ca_vacc_averted.png", p_bar, width = 10, height = 4.2, dpi = 120)

# ============================================================
# 4. Figure (c): averted burden with 95% UI (Monte Carlo)
# ============================================================
# Averted is TRUE-scale (real cases prevented), so carry the rho ~ Beta(20,60)
# uncertainty by scaling each draw by best_rho / rho_draws -- matching mc_tbl and the
# vaccinated_true_reported totals (so baseline - scenario = averted reconciles per draw).
f_rho <- best_rho / rho_draws
mc_long <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- av_draws[[nm]] * f_rho
  data.frame(timing  = scen_meta$timing[scen_meta$name == nm],
             arm     = scen_meta$arm[scen_meta$name == nm],
             outcome = outcomes,
             med = apply(m, 2, median, na.rm = TRUE),
             lo  = apply(m, 2, quantile, .025, na.rm = TRUE),
             hi  = apply(m, 2, quantile, .975, na.rm = TRUE),
             row.names = NULL)
}))
mc_long$timing  <- factor(mc_long$timing, levels = names(timings))
mc_long$arm     <- factor(mc_long$arm, levels = names(arms))   # disease-blocking first
mc_long$outcome <- factor(mc_long$outcome,
                          levels = c("infections", "symptomatic", "hospitalisations", "deaths"),
                          labels = c("Infections", "Symptomatic", "Hospitalisations", "Deaths"))

p_mc <- ggplot(mc_long, aes(timing, med, fill = arm)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8),
                width = 0.25, linewidth = 0.4) +
  facet_wrap(~ outcome, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Disease-blocking" = "#9ecae1",
                               "Disease + infection blocking" = "#08519c"), name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Burden averted, median + 95% UI",
       title = "Caldas Novas CHIKV cases - burden averted") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom", panel.grid.minor = element_blank())
print(p_mc)
ggsave("ca_vacc_averted_mc.png", p_mc, width = 10, height = 4.2, dpi = 120)

# ============================================================
# 4b. Figure (d): epidemic-curve with 95% UI ribbons, at 30% coverage
#     One plot per vaccination timing, mirroring the old per-coverage design.
#     y = predicted symptomatic cases; grey = no vaccine, blue = disease-blocking,
#     red = disease + infection blocking. Shaded band = vaccine rollout window.
# ============================================================
scen_cols <- c("No vaccination"               = "grey60",
               "Disease-blocking"             = "#4393c3",
               "Disease + infection blocking" = "#d6604d")

# paired % reduction in symptomatic cases (median, 95% UI) for one scenario
pct_reduction <- function(av_symp, base_symp) {
  q <- quantile(100 * av_symp / base_symp, c(.5, .025, .975), na.rm = TRUE)
  sprintf("%.1f%% (%.1f-%.1f%%)", q[1], q[2], q[3])
}

# weekly 2.5 / 50 / 97.5% across draws for one n_draws x T_weeks matrix
ribbon_df <- function(wk, scenario) {
  q <- apply(wk, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  data.frame(week = 1:T_weeks, lo = q[1, ], med = q[2, ], hi = q[3, ], scenario = scenario)
}

make_epicurve <- function(tn) {
  disb_nm <- paste0(tn, " | Disease-blocking")
  both_nm <- paste0(tn, " | Disease + infection blocking")
  pdf <- rbind(
    ribbon_df(wk_symp[["No vaccine (baseline)"]], "No vaccination"),
    ribbon_df(wk_symp[[disb_nm]],                 "Disease-blocking"),
    ribbon_df(wk_symp[[both_nm]],                 "Disease + infection blocking")
  )
  pdf$scenario <- factor(pdf$scenario, levels = names(scen_cols))

  lab <- paste0("% reduction (symptomatic):\n",
                "Disease-blocking: ",
                pct_reduction(av_draws[[disb_nm]][, "symptomatic"], base_draws[, "symptomatic"]), "\n",
                "Disease + infection blocking: ",
                pct_reduction(av_draws[[both_nm]][, "symptomatic"], base_draws[, "symptomatic"]))

  roll_start <- timings[[tn]]
  roll_end   <- roll_start + 1 / weekly_delivery_speed

  ggplot(pdf, aes(week, med, colour = scenario, fill = scenario)) +
    annotate("rect", xmin = roll_start, xmax = roll_end, ymin = -Inf, ymax = Inf,
             fill = "#3a7d3a", alpha = 0.10) +
    geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey70") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.02, vjust = 1.2,
             size = 2.8, lineheight = 0.95) +
    scale_colour_manual(values = scen_cols, aesthetics = c("colour", "fill")) +
    scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Week", y = "Predicted symptomatic CHIKV cases", colour = NULL, fill = NULL,
         title = paste0("Caldas Novas 2025-2026: symptomatic cases at 30% coverage"),
         caption = sprintf("Shaded band = vaccine deployment window (weeks %g-%g)",
                           roll_start, roll_end)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

epi_plots <- setNames(lapply(names(timings), make_epicurve), names(timings))
for (tn in names(timings)) {
  print(epi_plots[[tn]])
  fn <- sprintf("ca_vacc_epicurve_%s.png", gsub("[^a-z0-9]+", "_", tolower(tn)))
  ggsave(fn, epi_plots[[tn]], width = 8, height = 5, dpi = 120)
}

# ============================================================
# 5. Excel workbook: all outcome tables in one .xlsx
# ============================================================
# NOTE on deaths: all death figures below are OBSERVED-AGE corrected -- the engine
# (CHIKV_ca_vacc.R) sets age_weight = observed/model infection proportion by age, so
# deaths use the real ca_combined age distribution rather than the population
# structure. Infections/symptomatic/hospitalisations are unaffected. The `notes`
# sheet records the live run parameters, so the workbook always documents the run
# that produced it (any change in CHIKV_ca_vacc.R is reflected automatically).

# burden_mat (7 scenarios x 4 outcomes) -> tidy data frame with a scenario column
scenario_totals <- data.frame(scenario = rownames(burden_mat),
                              round(as.data.frame(burden_mat), 1),
                              row.names = NULL, check.names = FALSE)

# Per-scenario TOTAL burden (true & reported), median + 95% UI -- the counterpart to
# the averted table so you can cross-check totals vs averted. The per-draw scenario
# total is base_draws - av_draws (paired), which exactly recovers each draw's total;
# you cannot subtract the averted medians from the baseline medians because each UI is
# a quantile over the MC draws.
# Reporting-rate convention MATCHES caldas_baseline_burden (base_tbl in CHIKV_ca_vacc.R):
#   REPORTED = best_rho * total          -- anchored, rho-independent (narrow UI)
#   TRUE     = REPORTED / rho_draws      -- iceberg, carries rho ~ Beta(20,60) uncertainty
# so the "No vaccine" rows here reproduce caldas_baseline_burden exactly. Deaths w_a-corrected.
vac_scen <- c("No vaccine (baseline)", vac_names)
vaccinated_true_reported <- do.call(rbind, lapply(vac_scen, function(nm) {
  tot <- if (nm == "No vaccine (baseline)") base_draws else base_draws - av_draws[[nm]]
  data.frame(
    timing   = scen_meta$timing[scen_meta$name == nm],
    arm      = scen_meta$arm[scen_meta$name == nm],
    outcome  = outcomes,
    true     = vapply(outcomes, function(o) fmtq(best_rho * tot[, o] / rho_draws, if (o == "deaths") 1 else 0), character(1)),
    reported = vapply(outcomes, function(o) fmtq(best_rho * tot[, o],             if (o == "deaths") 1 else 0), character(1)),
    row.names = NULL, stringsAsFactors = FALSE)
}))

# Run parameters / assumptions, pulled live from the engine objects.
wa_on <- exists("age_weight") && !isTRUE(all.equal(age_weight, 1))
notes <- data.frame(
  parameter = c("Outbreak window", "Total reported cases (observed)",
                "Reporting rate (rho)", "Vaccine efficacy (IXCHIQ)",
                "Target coverage of 18-59", "Weekly delivery speed",
                "Dose-to-immunity delay (weeks)",
                "Deaths: observed-age (w_a) correction",
                "Start week_index: actual rollout",
                "Start week_index: start of 2026",
                "Start week_index: pre-outbreak"),
  value = c(sprintf("2025-W23 -> 2026-W22 (%d weeks)", T_weeks),
            format(sum(observed_cases), big.mark = ","),
            sprintf("%.2f", best_rho),
            sprintf("%.3f", ixchiq_efficacy),
            sprintf("%.0f%%", 100 * total_coverage),
            sprintf("%.2f (%g-week rollout)", weekly_delivery_speed, round(1 / weekly_delivery_speed)),
            as.character(immun_delay),
            if (wa_on) sprintf("ON (w range %.2f-%.2f)", min(age_weight), max(age_weight))
                       else "OFF (population structure)",
            as.character(timings[["actual rollout"]]),
            as.character(timings[["start of 2026"]]),
            as.character(timings[["pre-outbreak"]])),
  stringsAsFactors = FALSE)

sheets <- list(
  "notes"                    = notes,         # live run parameters / assumptions
  "baseline_true_reported"   = base_tbl,      # outcome, true, reported (median, 95% UI)
  "vaccinated_true_reported" = vaccinated_true_reported,  # per-scenario TOTALS, true & reported
  "averted_point"            = summary_tbl,   # timing, arm, *_averted, pct_* (numeric)
  "averted_MC_95UI"          = mc_tbl,        # timing, arm, outcome cols (median, 95% UI)
  "scenario_totals"          = scenario_totals# every scenario's total burden (point est.)
)
write_xlsx(sheets, "caldas_vacc_outputs.xlsx")
cat("Wrote caldas_vacc_outputs.xlsx (sheets:",
    paste(names(sheets), collapse = ", "), ")\n")
cat("Saved figures: ca_vacc_infections.png, ca_vacc_averted.png, ca_vacc_averted_mc.png,\n")
cat("               ca_vacc_epicurve_{actual_rollout,start_of_2026,pre_outbreak}.png\n")
