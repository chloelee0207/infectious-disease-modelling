library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(splines)

# NOTE: this script assumes you have already run seir_mle_fit.R and have in
# memory:
#   - best_beta_t        (fitted Fortaleza weekly beta_t, length 52)
#   - best_rho           (fitted reporting rate)
#   - N, R_init_prop, A, age_df  (Fortaleza population and immunity)
#   - I0, sigma, gamma, prop_symp, T_weeks, observed_cases
#
# If you're running this script standalone, source seir_mle_fit.R first or
# re-load those objects.

# ------------------------------------------------------------
# 1. SEIRV function with vaccination
# ------------------------------------------------------------
seirv_vaccinated <- function(
    T_weeks, 
    A, 
    N, 
    R_init_prop, 
    I0,
    base_beta, 
    sigma, 
    gamma, 
    rho,
    # Vaccine parameters:
    target_age,            # length-A binary vector; 1 = eligible
    total_coverage,        # final coverage of target pop (e.g., 0.30)
    weekly_delivery_speed, # fraction of supply delivered per week (e.g., 0.10)
    delay,                 # week at which vaccination begins (e.g., 2)
    VE_inf  = 0.989,       # infection-blocking efficacy
    VE_block = 0,          # disease-blocking efficacy
    immun_delay = 2,       # weeks from dose to immunity (Hyolim uses 2)
    prop_symp = 0.5242478,
    sub_steps = 7
) {
  
  pmax0 <- function(x) pmax(0, x)
  N_total <- sum(N)
  dt      <- 1 / sub_steps
  
  # Storage at WEEKLY resolution
  S <- E <- I <- R <- V <- matrix(0, nrow = A, ncol = T_weeks)
  V_covered <- matrix(0, nrow = A, ncol = T_weeks)
  vacc_delayed <- matrix(0, nrow = A, ncol = T_weeks)
  coverage_frac <- matrix(0, nrow = A, ncol = T_weeks)
  new_infections <- new_symptomatic <- matrix(0, nrow = A, ncol = T_weeks)
  
  # Vaccine supply machinery
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  total_supply      <- target_pop * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  total_avail_age <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  unvaccinated   <- N
  
  # Initial conditions at t = 1
  S_now <- N - I0 - R_init_prop * N
  E_now <- rep(0, A)
  I_now <- I0
  R_now <- R_init_prop * N
  V_now <- rep(0, A)
  
  S[, 1] <- S_now;  E[, 1] <- E_now;  I[, 1] <- I_now
  R[, 1] <- R_now;  V[, 1] <- V_now
  
  for (t in 2:T_weeks) {
    
    # ---- (a) People vaccinated `immun_delay` weeks ago become immune now
    if (t - immun_delay >= 1) {
      effective_dose <- vacc_delayed[, t - immun_delay]
      immunized      <- round(VE_inf * effective_dose)
      V_covered[, t] <- V_covered[, t - 1] + effective_dose
    } else {
      immunized      <- rep(0, A)
      V_covered[, t] <- V_covered[, t - 1]
    }
    
    # Move immunized people S -> V
    S_now <- pmax0(S_now - immunized)
    V_now <- V_now + immunized
    
    coverage_frac[, t] <- V_covered[, t] / N
    
    # ---- (b) Allocate this week's doses (if past delay)
    if (t >= delay && target_pop > 0) {
      rem <- weekly_dose_total
      for (a in target_idx) {
        alloc <- min(
          ceiling(weekly_dose_total * (N[a] / target_pop)),
          rem,
          unvaccinated[a],
          total_avail_age[a] - total_used_age[a]
        )
        if (alloc > 0) {
          prop_S    <- if (N[a] > 0) S_now[a] / N[a] else 0
          vacc_to_S <- round(alloc * prop_S)  # only S can be successfully vacc
          vacc_delayed[a, t] <- vacc_to_S
          total_used_age[a]  <- total_used_age[a] + alloc
          unvaccinated[a]    <- unvaccinated[a] - alloc
          rem                <- rem - alloc
        }
      }
    }
    
    # ---- (c) SEIR dynamics with sub-steps
    new_I_week <- rep(0, A)
    beta_t <- base_beta[t - 1]
    
    for (k in 1:sub_steps) {
      foi <- beta_t * sum(I_now) / N_total
      
      new_E <- foi   * S_now * dt
      new_I <- sigma * E_now * dt
      new_R <- gamma * I_now * dt
      
      S_now <- pmax0(S_now - new_E)
      E_now <- pmax0(E_now + new_E - new_I)
      I_now <- pmax0(I_now + new_I - new_R)
      R_now <- pmax0(R_now + new_R)
      
      new_I_week <- new_I_week + new_I
    }
    
    S[, t] <- S_now;  E[, t] <- E_now;  I[, t] <- I_now
    R[, t] <- R_now;  V[, t] <- V_now
    
    new_infections[, t]  <- new_I_week
    # Symptomatic cases reduced by disease-blocking efficacy in vaccinated
    new_symptomatic[, t] <- prop_symp * new_I_week *
      (1 - VE_block * coverage_frac[, t])
  }
  
  list(
    S = S, 
    E = E, 
    I = I, 
    R = R, 
    V = V,
    V_covered = V_covered,
    coverage_frac = coverage_frac,
    vacc_delayed  = vacc_delayed,
    total_supply  = total_supply,
    total_used_age = total_used_age,
    new_infections  = new_infections,
    new_symptomatic = new_symptomatic,
    new_reported    = rho * new_symptomatic
  )
}

# ============================================================
# 2a. Compute Fortaleza's share of allocated doses
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)

# Load weekly case data
raw <- read_excel("weekly_case.xlsx", sheet = "weekly_case")

# Compute total annual cases per municipality
totals <- raw |>
  pivot_longer(cols = starts_with("Week"),
               names_to  = "week",
               values_to = "cases") |>
  mutate(cases = ifelse(is.na(cases), 0, cases)) |>
  group_by(Code, Municipality, State) |>
  summarise(total_cases = sum(cases), .groups = "drop") |>
  arrange(desc(total_cases))

# Take top 10 by total cases, then exclude the highest as outlier
top10 <- totals |> slice_head(n = 10)

cat("Top 10 municipalities by total cases:\n")
print(top10)

# Compute Fortaleza's share of cases among the 9 remaining municipalities
fortaleza_total <- top10 |>
  filter(Code == 230440) |>
  pull(total_cases)

pool_total <- sum(top10$total_cases)

fortaleza_share <- fortaleza_total / pool_total

cat("\nFortaleza total cases:        ", fortaleza_total, "\n")
cat("Pool total:        ", pool_total, "\n")
cat("Fortaleza share:              ", round(100 * fortaleza_share, 1), "%\n")

# ------------------------------------------------------------
# 2b. Set vaccine parameters for the Brazil PVS scenario
# ------------------------------------------------------------
# Brazil PVS programme: 500,000 doses across 10 participating municipalities.
# Fortaleza represents ~41% of cases among the participating municipalities
# (Brazil MoH SINAN 2022; case-burden-weighted allocation assumed).
# Eligible population: adults 18-59
age_group <- as.character(age_df$age_group)
target_age <- rep(0, A)
target_age[c(4, 5, 6, 7, 8)] <- 1   
cat("Target age groups for vaccination:\n")
print(data.frame(age_group = age_group, eligible = target_age))
cat("Target population size:", sum(N[target_age == 1]), "\n")

national_supply  <- 500000
fortaleza_doses  <- national_supply * fortaleza_share

target_pop <- sum(N[target_age == 1])
total_coverage <- fortaleza_doses / target_pop

# Delivery speed: Hyolim used 10% per week
weekly_delivery_speed <- 0.10

# Delay: Hyolim used 2 weeks after outbreak onset
delay <- 2

# IXCHIQ efficacy: 98.9% per Hyolim's paper
ixchiq_efficacy <- 0.989

# Immune delay: 2 weeks from dose to immunity (Hyolim's default; IXCHIQ
# label suggests 14 days)
immun_delay <- 2

cat("Total doses to Fortaleza:    ", round(fortaleza_doses), "\n")
cat("Eligible 18-59 population:   ", target_pop, "\n")
cat("Implied coverage:            ", round(100 * total_coverage, 1), "%\n")
cat("  Weekly delivery:     ", round(sum(N[target_age == 1]) * total_coverage *
                                       weekly_delivery_speed), "doses/week\n")
cat("  Weeks to deplete:    ", round(1 / weekly_delivery_speed), "\n")
cat("  Rollout begins week: ", delay, "\n")

# ------------------------------------------------------------
# 3. Run the scenarios
# ------------------------------------------------------------
# Baseline: no vaccination
out_baseline <- seirv_vaccinated(
  T_weeks = T_weeks, 
  A = A, 
  N = N, 
  R_init_prop = R_init_prop,
  I0 = I0, 
  base_beta = best_beta_t,
  sigma = sigma, 
  gamma = gamma, 
  rho = best_rho,
  target_age = target_age,
  total_coverage = 0,    # zero coverage = no vaccination
  weekly_delivery_speed = weekly_delivery_speed,
  delay = delay,
  VE_inf = 0,
  VE_block = 0,
  immun_delay = immun_delay
)

# Scenario 1: disease-blocking ONLY (conservative floor)
out_disblock <- seirv_vaccinated(
  T_weeks = T_weeks, 
  A = A, 
  N = N, 
  R_init_prop = R_init_prop,
  I0 = I0, 
  base_beta = best_beta_t,
  sigma = sigma, 
  gamma = gamma, 
  rho = best_rho,
  target_age = target_age,
  total_coverage = total_coverage,
  weekly_delivery_speed = weekly_delivery_speed,
  delay = delay,
  VE_inf = 0,
  VE_block = ixchiq_efficacy,
  immun_delay = immun_delay
)

# Scenario 2: disease + infection blocking (optimistic ceiling)
out_both <- seirv_vaccinated(
  T_weeks = T_weeks, 
  A = A, 
  N = N, 
  R_init_prop = R_init_prop,
  I0 = I0, 
  base_beta = best_beta_t,
  sigma = sigma, 
  gamma = gamma, 
  rho = best_rho,
  target_age = target_age,
  total_coverage = total_coverage,
  weekly_delivery_speed = weekly_delivery_speed,
  delay = delay,
  VE_inf = ixchiq_efficacy,
  VE_block = ixchiq_efficacy,
  immun_delay = immun_delay
)

# ------------------------------------------------------------
# 4. Summary
# ------------------------------------------------------------
baseline_total <- sum(out_baseline$new_reported)
disblock_total <- sum(out_disblock$new_reported)
both_total     <- sum(out_both$new_reported)

averted_disblock     <- baseline_total - disblock_total
averted_both         <- baseline_total - both_total
incremental_infblock <- averted_both - averted_disblock   # extra benefit from infection-blocking

cat("\n--- Results ---\n")
cat("Baseline (no vaccine):              ", round(baseline_total), "reported cases\n")
cat("Disease-blocking only:              ", round(disblock_total), "reported cases\n")
cat("Disease + infection blocking:       ", round(both_total),     "reported cases\n")

cat("\nCases averted vs baseline:\n")
cat("  Disease-blocking only:    ", round(averted_disblock),
    " (", round(100 * averted_disblock / baseline_total, 1), "%)\n", sep = "")
cat("  Disease + infection:      ", round(averted_both),
    " (", round(100 * averted_both / baseline_total, 1), "%)\n", sep = "")

cat("\nIncremental benefit of infection-blocking",
    "(beyond disease-blocking floor):\n")
cat("  Additional cases averted: ", round(incremental_infblock),
    " (", round(100 * incremental_infblock / baseline_total, 1),
    "% of baseline)\n", sep = "")

# Coverage check
final_coverage <- max(colSums(out_both$V_covered) /
                        sum(N[target_age == 1]))
cat("\nFinal coverage achieved:", round(100 * final_coverage, 1), "%\n")

total_doses_allocated <- sum(out_both$total_used_age)
total_successful_vacc <- sum(out_both$V_covered[, T_weeks])
waste_fraction <- 1 - total_successful_vacc / total_doses_allocated

cat("Total doses allocated: ", round(total_doses_allocated), "\n")
cat("Successfully vaccinated:", round(total_successful_vacc), "\n")
cat("Dose wastage fraction: ", round(100 * waste_fraction, 1), "%\n")

# ------------------------------------------------------------
# 5. Plot the three scenarios
# ------------------------------------------------------------
plot_df <- data.frame(
  week = rep(1:T_weeks, 3),
  cases = c(colSums(out_baseline$new_reported),
            colSums(out_disblock$new_reported),
            colSums(out_both$new_reported)),
  scenario = rep(c("Baseline (no vaccine)",
                   "Disease-blocking only",
                   "Disease + infection blocking"), each = T_weeks)
)

plot_df$scenario <- factor(plot_df$scenario,
                           levels = c("Baseline (no vaccine)",
                                      "Disease-blocking only",
                                      "Disease + infection blocking"))

ggplot(plot_df, aes(x = week, y = cases, colour = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(data = data.frame(week = 1:T_weeks,
                               cases = observed_cases),
             aes(x = week, y = cases),
             colour = "black", size = 1.2, inherit.aes = FALSE) +
  scale_colour_manual(values = c("Baseline (no vaccine)"        = "grey",
                                 "Disease-blocking only"        = "#a8d1e7",
                                 "Disease + infection blocking" = "#f4a582")) +
  scale_x_continuous(breaks = seq(0, 52, 10)) +
  labs(x = "Week",
       y = "Reported CHIKV cases",
       title = "Fortaleza 2022: vaccinated vs baseline scenarios",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 11),
    legend.position  = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )