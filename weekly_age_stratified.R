library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
source("ca_common.R")   # load_caldas_age_cases(), load_burden_params(), fmtq(), ...

# ----------------------------------------------------------------------------
# Caldas Novas CHIKV weekly cases STRATIFIED BY AGE GROUP.
# Source: "ca_combined" sheet (SINAN download, age-group columns + Total).
# Outbreak window: 2025-W24 .. 2026-W22 = 52 weeks (30 in 2025, 22 in 2026).
# The 2025 epidemiological calendar has a Semana 53, so 2025 contributes weeks
# W24..W53 (30 weeks); 2025-W23 is excluded to give a clean 52-week window.
# ----------------------------------------------------------------------------

# Age-stratified cases from the canonical shared loader in ca_common.R. This is the
# SAME "ca_combined" data the fit (CHIKV_ca_pre_vacc_optim.R) uses -- single source of
# truth, so this script and the model can never drift apart.
age_levels <- c("<1 Ano", "1-4", "5-9", "10-14", "15-19",
                "20-39", "40-59", "60-64", "65-69", "70-79", "80 e +")
cc       <- load_caldas_age_cases()
ca_age   <- cc$ca_age     |> rename(epi_week = week)                     # age x week (zero-filled)
ca_total <- cc$caldas_obs |> rename(epi_week = week, tot_cases = cases)  # weekly totals
observed_cases <- cc$observed_cases
T_weeks <- length(observed_cases)

stopifnot(T_weeks == 52, sum(observed_cases) == 8204)
cat("Caldas Novas age-stratified CHIKV (2025-W24 to 2026-W22)\n")
cat("  weeks:", T_weeks, " total cases:", sum(observed_cases), "\n")
cat("  peak week:", ca_total$week_label[which.max(observed_cases)],
    "with", max(observed_cases), "cases\n")

# Cases by age group over the whole window
age_totals <- ca_age |>
  group_by(age_group) |>
  summarise(cases = sum(cases), .groups = "drop") |>
  mutate(pct = cases / sum(cases) * 100)
cat("\nCases by age group:\n")
print(as.data.frame(age_totals), row.names = FALSE)

# ----------------------------------------------------------------------------
# Severity per age group -- hospitalisations & deaths with 95% uncertainty.
# Parameters from Hyolim Table S4 (disease_progression.xlsx, via ca_common.R):
#   - hospitalisation among symptomatic: single Beta rate (~4.0%, not age-specific)
#   - death among hospitalised cases   : Beta by 9 decadal bands [0,10)..[80,90)
#   - death among non-hospitalised     : Beta by the same 9 bands
# The ca_combined age groups don't align 1:1 with the decadal bands, so each is
# mapped to the band(s) it covers; groups spanning two bands (20-39, 40-59) take
# the mean of those bands. Uncertainty is propagated by Monte Carlo draws from
# the Beta(alpha, beta) hyperparameters.
# ----------------------------------------------------------------------------
bp <- load_burden_params(12)          # exposes hosp_a/b and per-band cfr Beta params (ca_common.R)

set.seed(42)
N_mc <- 50000
hosp_draw  <- rbeta(N_mc, bp$hosp_a, bp$hosp_b)                               # single hosp rate
cfr_h_draw <- sapply(1:9, function(j) rbeta(N_mc, bp$cfr_hosp_a[j], bp$cfr_hosp_b[j]))  # N x 9
cfr_n_draw <- sapply(1:9, function(j) rbeta(N_mc, bp$cfr_nonh_a[j], bp$cfr_nonh_b[j]))  # N x 9

# ca_combined age group -> decadal band index/indices (1=[0,10) .. 9=[80,90))
group_to_bands <- list(
  "<1 Ano" = 1, "1-4" = 1, "5-9" = 1,
  "10-14" = 2, "15-19" = 2,
  "20-39" = c(3, 4), "40-59" = c(5, 6),
  "60-64" = 7, "65-69" = 7, "70-79" = 8, "80 e +" = 9
)
qs <- function(x) as.numeric(quantile(x, c(.5, .025, .975)))   # median, lo, hi

# Per-group Monte Carlo, accumulating totals across groups on the same draws
tot_hosp <- tot_dh <- tot_dtot <- numeric(N_mc)
sev_list <- lapply(names(group_to_bands), function(g) {
  b   <- group_to_bands[[g]]
  n_g <- age_totals$cases[age_totals$age_group == g]
  cfr_h_g <- if (length(b) == 1) cfr_h_draw[, b] else rowMeans(cfr_h_draw[, b])
  cfr_n_g <- if (length(b) == 1) cfr_n_draw[, b] else rowMeans(cfr_n_draw[, b])
  hosp_g  <- n_g * hosp_draw                                   # hospitalisations
  dh_g    <- hosp_g * cfr_h_g                                  # deaths among hospitalised
  dtot_g  <- n_g * (hosp_draw * cfr_h_g + (1 - hosp_draw) * cfr_n_g)  # total deaths
  tot_hosp <<- tot_hosp + hosp_g; tot_dh <<- tot_dh + dh_g; tot_dtot <<- tot_dtot + dtot_g
  qh <- qs(hosp_g); qd <- qs(dh_g); qt <- qs(dtot_g); qr <- qs(cfr_h_g)
  data.frame(age_group = g, symptomatic = n_g,
             death_per_hosp_med = qr[1], death_per_hosp_lo = qr[2], death_per_hosp_hi = qr[3],
             hosp_med = qh[1], hosp_lo = qh[2], hosp_hi = qh[3],
             deaths_hosp_med = qd[1], deaths_hosp_lo = qd[2], deaths_hosp_hi = qd[3],
             deaths_total_med = qt[1], deaths_total_lo = qt[2], deaths_total_hi = qt[3])
})
severity_tab <- do.call(rbind, sev_list)

# TOTAL row (from the summed draws, so the interval reflects correlated hosp rate)
qh <- qs(tot_hosp); qd <- qs(tot_dh); qt <- qs(tot_dtot)
severity_tab <- rbind(severity_tab, data.frame(
  age_group = "TOTAL", symptomatic = sum(severity_tab$symptomatic),
  death_per_hosp_med = NA, death_per_hosp_lo = NA, death_per_hosp_hi = NA,
  hosp_med = qh[1], hosp_lo = qh[2], hosp_hi = qh[3],
  deaths_hosp_med = qd[1], deaths_hosp_lo = qd[2], deaths_hosp_hi = qd[3],
  deaths_total_med = qt[1], deaths_total_lo = qt[2], deaths_total_hi = qt[3]))

# Add the (non-age-specific) hospitalisation rate columns for the export
severity_tab$hosp_rate_med <- qs(hosp_draw)[1]
severity_tab$hosp_rate_lo  <- qs(hosp_draw)[2]
severity_tab$hosp_rate_hi  <- qs(hosp_draw)[3]

# Pretty console print: "median (lo - hi)"
mlh  <- function(m, l, h, d = 1) ifelse(is.na(m), "-",
          sprintf(paste0("%.", d, "f (%.", d, "f - %.", d, "f)"), m, l, h))
pct3 <- function(m, l, h) ifelse(is.na(m), "-",
          sprintf("%.2f%% (%.2f - %.2f)", 100*m, 100*l, 100*h))
cat("\nSeverity per age group (hospitalisation rate = 4.0%, not age-specific):\n")
print(data.frame(
  age_group        = severity_tab$age_group,
  symptomatic      = severity_tab$symptomatic,
  `death|hosp`     = pct3(severity_tab$death_per_hosp_med, severity_tab$death_per_hosp_lo, severity_tab$death_per_hosp_hi),
  hospitalisations = mlh(severity_tab$hosp_med, severity_tab$hosp_lo, severity_tab$hosp_hi),
  `deaths (hosp)`  = mlh(severity_tab$deaths_hosp_med, severity_tab$deaths_hosp_lo, severity_tab$deaths_hosp_hi, 2),
  `deaths (total)` = mlh(severity_tab$deaths_total_med, severity_tab$deaths_total_lo, severity_tab$deaths_total_hi, 2),
  check.names = FALSE), row.names = FALSE)

# ---- Export: hospitalisations & deaths per age group (median + 95% UI) ----
export <- severity_tab |>
  dplyr::select(age_group, symptomatic,
                hosp_rate_med, hosp_rate_lo, hosp_rate_hi,
                death_per_hosp_med, death_per_hosp_lo, death_per_hosp_hi,
                hosp_med, hosp_lo, hosp_hi,
                deaths_hosp_med, deaths_hosp_lo, deaths_hosp_hi,
                deaths_total_med, deaths_total_lo, deaths_total_hi)
write.csv(export, "CHIKV_ca_severity_by_age.csv", row.names = FALSE)
if (requireNamespace("writexl", quietly = TRUE))
  writexl::write_xlsx(list(severity_by_age = export), "CHIKV_ca_severity_by_age.xlsx")
cat("\nWrote CHIKV_ca_severity_by_age.csv",
    if (requireNamespace("writexl", quietly = TRUE)) "and CHIKV_ca_severity_by_age.xlsx" else "", "\n")

# ----------------------------------------------------------------------------
# True number of infections (under-ascertainment correction)
#   true infections = reported cases / P(symptomatic) / P(reported)
# P(symptomatic): prop_symp ~ Beta(ps_a, ps_b) (Hyolim; mean 52.4%, 95% 40.6-64.0);
#                 already loaded in bp above.
# P(reported)   : reporting rate rho ~ Beta(20, 60) (Hyolim supplement; mean 25%,
#                 95% 16.2-35.0). NB this is her stated generative Beta and is a
#                 touch wider than the paper's headline "25% (20.1-32.5)"; using it
#                 keeps this consistent with CHIKV_ca_vacc.R. To honour the headline
#                 UI instead, swap in rbeta(N_mc, 49.76, 141.17).
# Uncertainty in both is propagated by Monte Carlo through the two divisions.
# ----------------------------------------------------------------------------
ps_draw  <- rbeta(N_mc, bp$ps_a, bp$ps_b)
rho_draw <- rbeta(N_mc, 20, 60)
true_inf_draw <- sum(observed_cases) / ps_draw / rho_draw

ti_point <- sum(observed_cases) / (bp$ps_a / (bp$ps_a + bp$ps_b)) / 0.25
ti_q     <- qs(true_inf_draw)                 # median, 2.5%, 97.5%
cat(sprintf("\nTrue infections = %s reported cases / prop_symp / reporting_rate\n",
            format(sum(observed_cases), big.mark = ",")))
cat(sprintf("  point estimate (0.524, 0.25) : %s\n", format(round(ti_point), big.mark = ",")))
cat(sprintf("  Monte Carlo median (95%% UI)  : %s (%s - %s)\n",
            format(round(ti_q[1]), big.mark = ","),
            format(round(ti_q[2]), big.mark = ","),
            format(round(ti_q[3]), big.mark = ",")))

# ---- Axis helpers (mirror weekly_case.R) ----
x_ticks <- ca_total |>
  filter((Year == 2025 & epi_week %in% c(30, 40, 50)) |
         (Year == 2026 & epi_week %in% c(10, 20))) |>
  mutate(label = as.character(epi_week))

year_break <- mean(c(
  max(ca_total$week_index[ca_total$Year == 2025]),
  min(ca_total$week_index[ca_total$Year == 2026])
))

# ---- Plot 1: weekly total ----
caldas_weekly <- ggplot(ca_total, aes(x = week_index, y = tot_cases)) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  annotate("text", x = year_break, y = 0,
           label = "2026", angle = 90, vjust = -0.4, hjust = -9.3,
           fontface = "bold", size = 3.6, colour = "grey40") +
  geom_line(linewidth = 0.6) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$label) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week", y = "Reported CHIKV cases",
       title = "Observed Weekly CHIKV cases in Caldas Novas") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_line(linetype = "dotted", colour = "grey80"),
        panel.grid.minor = element_blank())

print(caldas_weekly)
ggsave("CHIKV_ca_weekly.png", caldas_weekly, width = 9, height = 5, dpi = 120)

# ---- Plot 2: weekly cases stacked by age group ----
caldas_age <- ggplot(ca_age, aes(x = week_index, y = cases, fill = age_group)) +
  geom_col(width = 1) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey30") +
  annotate("text", x = year_break, y = 0,
           label = "2026", angle = 90, vjust = -0.4, hjust = -9.3,
           fontface = "bold", size = 3.6, colour = "grey40") +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$label) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_viridis_d(option = "turbo", name = "Age group") +
  labs(x = "Week", y = "Reported CHIKV cases",
       title = "Observed Weekly CHIKV cases in Caldas Novas by age group") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_line(linetype = "dotted", colour = "grey80"),
        panel.grid.minor = element_blank())

print(caldas_age)
ggsave("CHIKV_ca_weekly_age.png", caldas_age, width = 9, height = 4.5, dpi = 120)

# ============================================================================
# SENSE CHECK: expected CHIKV deaths in Caldas Novas via the Goias state CFR
# ----------------------------------------------------------------------------
# We cannot observe how many people died of chikungunya in Caldas Novas in
# 2025-2026, so we anchor to Goias state and scale down:
#   1. Goias 2024 chikungunya cases      -- weekly_all, State == Goias, Year 2024
#   2. Goias 2024 arbovirus/VHF deaths   -- mortality_2024.xlsx, CID-BR-10 row
#      "020 Outras febres por arbovirus e febres hemorragicas virais"
#   3. CFR = deaths / cases  (assumes ALL those arbovirus deaths are chikungunya)
#   4. Goias cases over the outbreak window (2025-W24 .. 2026-W22)
#   5. Caldas Novas population as a share of Goias  -- population.xlsx, 2025 sheet
#   6. Estimated Caldas deaths = CFR x Goias-window cases x population share
# ============================================================================
wa <- suppressMessages(read_excel("weekly_case.xlsx", sheet = "weekly_all")) |>
  filter(State == "Goiás") |>
  pivot_longer(starts_with("Week"), names_to = "w", values_to = "c") |>
  mutate(week = as.integer(sub("Week ", "", w)),
         c = ifelse(is.na(c), 0, c))

goias_cases_2024 <- sum(wa$c[wa$Year == 2024])
# Outbreak window 2025-W24 .. 2026-W22 (INCLUDING 2025-W53), matching the Caldas
# analysis window exactly (2025 has an epidemiological Semana 53).
goias_win <- sum(wa$c[(wa$Year == 2025 & wa$week >= 24) |
                      (wa$Year == 2026 & wa$week <= 22)])
# For reference only: the same window excluding 2025-W53.
goias_win_excl53 <- sum(wa$c[(wa$Year == 2025 & wa$week >= 24 & wa$week <= 52) |
                             (wa$Year == 2026 & wa$week <= 22)])

# --- Goias arbovirus/VHF deaths 2024 (CID-BR-10 code 020) ---
mort <- suppressMessages(read_excel("mortality_2024.xlsx", sheet = 1, col_names = FALSE))[[1]]
hdr  <- mort[grepl('RO";"AC', mort)][1]                 # header line with state codes
states <- gsub('"', '', strsplit(hdr, ";")[[1]])        # [1]=label, then RO..DF, Total
row020 <- mort[grepl("020 Out febres", mort)][1]
vals   <- strsplit(row020, ";")[[1]]                    # [1]=label, aligned with states
goias_arbo_deaths_2024 <- as.numeric(gsub("-", "0", vals[which(states == "GO")]))

goias_cfr <- goias_arbo_deaths_2024 / goias_cases_2024

# --- Caldas Novas population share of Goias (IBGE 2025 estimate) ---
popn <- suppressWarnings(suppressMessages(
          read_excel("population.xlsx", sheet = "2025", skip = 1)))
names(popn)[1:5] <- c("uf", "cod_uf", "cod_munic", "municipio", "pop")
popn$pop <- as.numeric(popn$pop)
goias_pop  <- sum(popn$pop[popn$uf == "GO"], na.rm = TRUE)
caldas_pop <- popn$pop[popn$uf == "GO" & grepl("Caldas Novas", popn$municipio)]
pop_share  <- caldas_pop / goias_pop

# --- Estimate (user's approach: scale state deaths by population share) ---
goias_deaths_win <- goias_cfr * goias_win
est_caldas_deaths <- goias_deaths_win * pop_share

cat("\n================ SENSE CHECK: deaths in Caldas Novas ================\n")
cat(sprintf("Goias chikungunya cases 2024              : %s\n", format(goias_cases_2024, big.mark = ",")))
cat(sprintf("Goias arbovirus/VHF deaths 2024 (CID 020) : %d\n", goias_arbo_deaths_2024))
cat(sprintf("Implied Goias CFR (2024)                  : %.3f%%\n", 100 * goias_cfr))
cat(sprintf("Goias cases, window 2025-W24..2026-W22    : %s\n",
            format(goias_win, big.mark = ",")))
cat(sprintf("Caldas Novas / Goias population share     : %s / %s = %.3f%%\n",
            format(caldas_pop, big.mark = ","), format(goias_pop, big.mark = ","), 100 * pop_share))
cat(sprintf("Implied Goias deaths in window            : %.1f\n", goias_deaths_win))
cat(sprintf(">> Estimated Caldas Novas CHIKV deaths     : %.2f  (population-share scaling)\n", est_caldas_deaths))

# --- Alternative: apply the CFR to Caldas's OWN observed cases -------------
# Caldas Novas holds a large share of Goias' window cases, so population-share
# scaling understates its deaths. Applying the state CFR directly to Caldas'
# own case count is an upper-bound style comparison.
caldas_case_share <- sum(observed_cases) / goias_win
est_caldas_deaths_cases <- goias_cfr * sum(observed_cases)
cat(sprintf("\nContext: Caldas holds %.0f%% of Goias' window cases (%s of %s).\n",
            100 * caldas_case_share, format(sum(observed_cases), big.mark = ","),
            format(goias_win, big.mark = ",")))
cat(sprintf("If the Goias CFR is applied to Caldas' OWN %s cases -> %.1f deaths.\n",
            format(sum(observed_cases), big.mark = ","), est_caldas_deaths_cases))
cat(sprintf("Model-based estimate (this script)         : %.1f total deaths (95%% UI %.1f-%.1f)\n",
            qs(tot_dtot)[1], qs(tot_dtot)[2], qs(tot_dtot)[3]))
cat("=====================================================================\n")
