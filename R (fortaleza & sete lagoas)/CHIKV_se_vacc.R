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
    E0,
    base_beta,
    sigma, 
    gamma, 
    rho,
    # Vaccine parameters:
    target_age,            # length-A binary vector; 1 = eligible
    total_coverage,        # final coverage of target pop
    weekly_delivery_speed, # fraction of supply delivered per week
    delay,                 # week at which vaccination begins
    VE_inf  = 0.989,       # infection-blocking efficacy
    VE_block = 0,          # disease-blocking efficacy
    immun_delay = 2,       # weeks from dose to immunity
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
  
  # Working state at the start of week 1. Seed E0 (exposed) as well as I0 so the
  # E->I flow (sigma * E0) reproduces week-1 incidence, exactly as the fitting
  # model does. Without seeding E0 the epidemic spins up from zero and week 1 = 0.
  S_now <- pmax0(N - I0 - E0 - R_init_prop * N)
  E_now <- E0
  I_now <- I0
  R_now <- R_init_prop * N
  V_now <- rep(0, A)

  # Loop from week 1 so incidence is recorded for every week (including week 1).
  # At t = 1 the vaccination steps below are inert (delay >= 2, immun_delay >= 2),
  # so the baseline scenario reproduces the fitted predicted curve exactly.
  for (t in 1:T_weeks) {

    # ---- (a) People vaccinated `immun_delay` weeks ago become immune now
    prev_V_covered <- if (t > 1) V_covered[, t - 1] else rep(0, A)
    if (t - immun_delay >= 1) {
      effective_dose <- vacc_delayed[, t - immun_delay]
      immunized      <- round(VE_inf * effective_dose)
      V_covered[, t] <- prev_V_covered + effective_dose
    } else {
      immunized      <- rep(0, A)
      V_covered[, t] <- prev_V_covered
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
    beta_t <- base_beta[t]
    
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
# 2. Set vaccine parameters for the target-coverage scenario
# ============================================================
# Brazil's IXCHIQ rollout has been extended beyond the original 10
# municipalities, and the number of doses allocated to Sete Lagoas has not
# been publicly announced. Rather than back-calculating doses from a national
# supply (the previous Fortaleza approach), we model explicit *target
# coverage* of the eligible population (adults 18-59): a low (20%) and a high
# (40%) bound, bracketing the plausible programme reach.

# Eligible population: adults 18-59
age_group <- as.character(age_df$age_group)
target_age <- rep(0, A)
target_age[c(4, 5, 6, 7, 8)] <- 1
cat("Target age groups for vaccination:\n")
print(data.frame(age_group = age_group, eligible = target_age))

target_pop <- sum(N[target_age == 1])
cat("Eligible 18-59 population:", target_pop, "\n")

# Target coverage of the eligible population (low / high bounds)
coverage_low  <- 0.20
coverage_high <- 0.40

# Delivery speed: 10% of supply per week
weekly_delivery_speed <- 0.10

# Delay: vaccination begins 2 weeks after outbreak onset
delay <- 2

# IXCHIQ efficacy: 98.9% per Hyolim's paper
ixchiq_efficacy <- 0.989

# Immune delay: 2 weeks from dose to immunity (Hyolim's default; IXCHIQ
# label suggests 14 days)
immun_delay <- 2

for (cov in c(coverage_low, coverage_high)) {
  cat(sprintf("\nCoverage target: %.0f%% of eligible 18-59\n", 100 * cov))
  cat("  Doses (max, to reach target):", round(target_pop * cov), "\n")
  cat("  Weekly delivery:             ",
      round(target_pop * cov * weekly_delivery_speed), "doses/week\n")
}
cat("\n  Weeks to deplete:    ", round(1 / weekly_delivery_speed), "\n")
cat("  Rollout begins week: ", delay, "\n")

# ------------------------------------------------------------
# 3. Run the scenarios
# ------------------------------------------------------------
# Helper: run one scenario at a given coverage + efficacy combination.
run_scenario <- function(coverage, VE_inf, VE_block) {
  seirv_vaccinated(
    T_weeks = T_weeks,
    A = A,
    N = N,
    R_init_prop = R_init_prop,
    I0 = I0,
    E0 = E0,
    base_beta = best_beta_t,
    sigma = sigma,
    gamma = gamma,
    rho = best_rho,
    target_age = target_age,
    total_coverage = coverage,
    weekly_delivery_speed = weekly_delivery_speed,
    delay = delay,
    VE_inf = VE_inf,
    VE_block = VE_block,
    immun_delay = immun_delay
  )
}

# Baseline: no vaccination
out_baseline <- run_scenario(0, VE_inf = 0, VE_block = 0)

# Disease-blocking ONLY (conservative floor) at 20% and 40% coverage
out_disblock_low  <- run_scenario(coverage_low,  VE_inf = 0, VE_block = ixchiq_efficacy)
out_disblock_high <- run_scenario(coverage_high, VE_inf = 0, VE_block = ixchiq_efficacy)

# Disease + infection blocking (optimistic ceiling) at 20% and 40% coverage
out_both_low  <- run_scenario(coverage_low,  VE_inf = ixchiq_efficacy, VE_block = ixchiq_efficacy)
out_both_high <- run_scenario(coverage_high, VE_inf = ixchiq_efficacy, VE_block = ixchiq_efficacy)

# ------------------------------------------------------------
# 4. Summary
# ------------------------------------------------------------
scenarios <- list(
  "Baseline (no vaccine)"             = out_baseline,
  "Disease-blocking only, 20%"        = out_disblock_low,
  "Disease-blocking only, 40%"        = out_disblock_high,
  "Disease + infection blocking, 20%" = out_both_low,
  "Disease + infection blocking, 40%" = out_both_high
)

baseline_total <- sum(out_baseline$new_reported)

cat("\n--- Results (reported cases over the season) ---\n")
for (nm in names(scenarios)) {
  tot     <- sum(scenarios[[nm]]$new_reported)
  averted <- baseline_total - tot
  if (nm == "Baseline (no vaccine)") {
    cat(sprintf("  %-34s %8.0f\n", nm, tot))
  } else {
    cat(sprintf("  %-34s %8.0f   averted %7.0f (%4.1f%%)\n",
                nm, tot, averted, 100 * averted / baseline_total))
  }
}

# Incremental benefit of infection-blocking beyond the disease-blocking floor
for (cov in c("20%", "40%")) {
  db <- scenarios[[paste("Disease-blocking only,", cov)]]
  bo <- scenarios[[paste("Disease + infection blocking,", cov)]]
  inc <- sum(db$new_reported) - sum(bo$new_reported)
  cat(sprintf("\nIncremental benefit of infection-blocking at %s coverage: %0.0f (%0.1f%% of baseline)\n",
              cov, inc, 100 * inc / baseline_total))
}

# Coverage + dose diagnostics for each coverage level
for (nm in c("Disease + infection blocking, 20%",
             "Disease + infection blocking, 40%")) {
  out <- scenarios[[nm]]
  final_coverage <- max(colSums(out$V_covered) / target_pop)
  total_doses_allocated <- sum(out$total_used_age)
  total_successful_vacc <- sum(out$V_covered[, T_weeks])
  waste_fraction <- 1 - total_successful_vacc / total_doses_allocated
  cat(sprintf("\n[%s]\n", nm))
  cat("  Final coverage achieved:", round(100 * final_coverage, 1), "%\n")
  cat("  Total doses allocated:  ", round(total_doses_allocated), "\n")
  cat("  Successfully vaccinated: ", round(total_successful_vacc), "\n")
  cat("  Dose wastage fraction:  ", round(100 * waste_fraction, 1), "%\n")
}

# ------------------------------------------------------------
# 5. Plot the three scenarios
# ------------------------------------------------------------
plot_df <- do.call(rbind, lapply(names(scenarios), function(nm) {
  data.frame(week     = 1:T_weeks,
             cases    = colSums(scenarios[[nm]]$new_reported),
             scenario = nm)
}))

plot_df$scenario <- factor(plot_df$scenario, levels = names(scenarios))

ggplot(plot_df, aes(x = week, y = cases, colour = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(data = data.frame(week = 1:T_weeks,
                               cases = observed_cases),
             aes(x = week, y = cases),
             colour = "black", size = 1.2, inherit.aes = FALSE) +
  scale_colour_manual(values = c(
    "Baseline (no vaccine)"             = "grey40",
    "Disease-blocking only, 20%"        = "#a8d1e7",
    "Disease-blocking only, 40%"        = "#2166ac",
    "Disease + infection blocking, 20%" = "#f4a582",
    "Disease + infection blocking, 40%" = "#b2182b")) +
  scale_x_continuous(breaks = seq(0, 52, 10)) +
  labs(x = "Week",
       y = "Reported CHIKV cases",
       title = "Sete Lagoas 2024: vaccinated vs baseline scenarios",
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