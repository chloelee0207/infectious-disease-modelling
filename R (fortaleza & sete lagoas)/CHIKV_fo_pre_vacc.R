# install.packages("readxl")  # run once if not already installed
library(readxl)

seir_baseline <- function(
    T_weeks,        # number of weeks to simulate
    A,              # number of age groups
    N,              # vector length A: population by age group
    R_init_prop,    # vector length A: proportion already immune by age
    I0,             # vector length A: initial infections by age
    base_beta,      # vector length T_weeks: weekly transmission rate
    sigma,          # 1 / latent period (per week)
    gamma,          # recovery rate (per week)
    rho,            # reporting rate
    prop_symp = 0.5242478,  # proportion symptomatic
    sub_steps = 7       # sub-steps per week (7 = daily updates)
) {
  
  pmax0 <- function(x) pmax(0, x)
  N_total <- sum(N)
  dt      <- 1 / sub_steps   # length of each sub-step in "weeks"
  
  # Storage at WEEKLY resolution
  S <- E <- I <- R <- matrix(0, nrow = A, ncol = T_weeks)
  new_infections <- new_symptomatic <- matrix(0, nrow = A, ncol = T_weeks)
  
  # Working state at fine resolution
  S_now <- N - I0 - R_init_prop * N
  E_now <- rep(0, A)
  I_now <- I0
  R_now <- R_init_prop * N
  
  # Record t = 1
  S[, 1] <- S_now;  E[, 1] <- E_now;  I[, 1] <- I_now;  R[, 1] <- R_now
  
  for (t in 2:T_weeks) {
    
    # Track new infections accumulated within this week
    new_I_week <- rep(0, A)
    
    # beta is constant across the sub-steps within a week
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
    
    S[, t] <- S_now;  E[, t] <- E_now;  I[, t] <- I_now;  R[, t] <- R_now
    new_infections[, t]  <- new_I_week
    new_symptomatic[, t] <- prop_symp * new_I_week
  }
  
  list(
    S = S, E = E, I = I, R = R,
    new_infections  = new_infections,
    new_symptomatic = new_symptomatic,
    new_reported    = rho * new_symptomatic   # reported = rho * symptomatic
  )
}

municipality_name <- "fortaleza"      # match value in population.xlsx
state_name        <- "Ceará"          # match value in beta_hyolim.xls

# Load Hyolim's beta_t posteriors and filter for Ceará
beta_all <- read_excel("beta_hyolim.xlsx", sheet = "beta")
normalize <- function(x) tolower(gsub("[áàâã]", "a",
                                      gsub("[éèê]",   "e",
                                           gsub("[íì]",    "i",
                                                gsub("[óòôõ]", "o",
                                                     gsub("[úù]",   "u",
                                                          gsub("ç",      "c", x)))))))

beta_state <- beta_all[normalize(beta_all$region) == normalize(state_name), ]
stopifnot(nrow(beta_state) > 0)

beta_state$week <- as.integer(sub("beta_week", "", beta_state$parameter))
beta_state      <- beta_state[order(beta_state$week), ]

base_beta <- beta_state$median
T_weeks   <- length(base_beta)
stopifnot(T_weeks == 52)

# Load municipality age structure and baseline immunity
age_df <- read_excel("population.xlsx", sheet = "prop_immune")
age_df <- age_df[normalize(age_df$municipality) == normalize(municipality_name), ]
stopifnot(nrow(age_df) > 0)
age_df$age_group <- factor(age_df$age_group, levels = age_df$age_group)

age_group   <- as.character(age_df$age_group)
N           <- age_df$pop_num
R_init_prop <- age_df$prop_ever_infected
age_midpt   <- age_df$age_midpoint

A <- nrow(age_df)
stopifnot(length(N) == A, length(R_init_prop) == A)

cat("Running model for: ", municipality_name, " (", state_name, ")\n", sep = "")
cat("Loaded", A, "age groups\n")
cat("Total population: ", sum(N), "\n")
cat("Population-weighted baseline immunity: ",
    round(sum(R_init_prop * N) / sum(N) * 100, 1), "%\n")

# Disease and reporting parameters
sigma <- 1 / 0.60    # 1 / latent period (weeks)
gamma <- 0.54        # recovery rate (per week)
rho   <- 0.18        # reporting rate (national prior; replace with Ceará posterior)

# Initial infections
week_1_cases <- 10
rho_for_seed <- 0.10  # only used to back-calculate the seed
prop_symp    <- 0.5242478
I0_total <- week_1_cases / rho_for_seed / prop_symp
susceptible_pop <- N * (1 - R_init_prop)
I0 <- round(I0_total * susceptible_pop / sum(susceptible_pop))

cat("Back-calculated total initial infections: ", round(I0_total, 1), "\n")
cat("Seeded I0 by age group:\n")
print(setNames(I0, age_group))

# ---- Run ----------------------------------------------------------
out <- seir_baseline(
  T_weeks     = T_weeks,
  A           = A,
  N           = N,
  R_init_prop = R_init_prop,
  I0          = I0,
  base_beta   = base_beta,
  sigma       = sigma,
  gamma       = gamma,
  rho         = rho
)

# # ---- Summary ------------------------------------------------------
cat("\nTotal infections:                  ", round(sum(out$new_infections)),  "\n")
cat("Total symptomatic cases (true):    ", round(sum(out$new_symptomatic)), "\n")
cat("Total reported cases (predicted):  ", round(sum(out$new_reported)),    "\n")
# 
# # ---- Plots: all line charts --------------------------------------
# weeks   <- 1:T_weeks
# palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
#              "#000000")
# 
# par(mfrow = c(4, 1), mar = c(4, 4.5, 2, 1))
# 
# # (1) Input beta_t
# plot(weeks, base_beta, type = "l", lwd = 2, col = "grey30",
#      xlab = "Week", ylab = expression(beta[t]),
#      main = "Ceará weekly transmission rate (Hyolim's posterior median)")
# 
# # (2) Total predicted reported cases (compare with SINAN later)
# plot(weeks, colSums(out$new_reported), type = "l", lwd = 2, col = "darkorange",
#      xlab = "Week", ylab = "Reported cases",
#      main = "Predicted weekly reported CHIKV cases, Fortaleza 2022")
# 
# # (3) True symptomatic cases by age (NOT underreporting-adjusted)
# matplot(weeks, t(out$new_symptomatic), type = "l", lty = 1, lwd = 1.8,
#         col = palette[1:A],
#         xlab = "Week", ylab = "True symptomatic cases",
#         main = "Weekly TRUE symptomatic cases by age group (no underreporting)")
# legend("topright", legend = age_group, col = palette[1:A], lty = 1, lwd = 1.8,
#        cex = 0.65, ncol = 2, bty = "n")
# 
# # (4) Reported cases by age (underreporting-adjusted; compare to SINAN by age)
# matplot(weeks, t(out$new_reported), type = "l", lty = 1, lwd = 1.8,
#         col = palette[1:A],
#         xlab = "Week", ylab = "Reported cases",
#         main = "Weekly REPORTED cases by age group (rho = 0.25 applied)")
# legend("topright", legend = age_group, col = palette[1:A], lty = 1, lwd = 1.8,
#        cex = 0.65, ncol = 2, bty = "n")
# 
# par(mfrow = c(1, 1))

# ============================================================
# Overlay predicted vs. observed reported cases (Fortaleza)
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Observed: load and reshape Fortaleza SINAN data --------------
raw <- read_excel("weekly_case.xlsx", sheet = "weekly_case")

fortaleza_obs <- raw |>
  filter(Code == 230440) |>
  pivot_longer(
    cols      = starts_with("Week"),
    names_to  = "week",
    values_to = "cases"
  ) |>
  mutate(
    week  = as.integer(sub("Week ", "", week)),
    cases = ifelse(is.na(cases), 0, cases),
    source = "Observed (SINAN)"
  ) |>
  select(week, cases, source)

stopifnot(sum(fortaleza_obs$cases) == 29660)

# ---- Predicted: weekly reported cases from the SEIR ---------------
predicted <- data.frame(
  week   = 1:T_weeks,
  cases  = colSums(out$new_reported),
  source = "Predicted (SEIR)"
)

# ---- Combine and plot --------------------------------------------
combined <- bind_rows(fortaleza_obs, predicted)

ggplot(combined, aes(x = week, y = cases, colour = source)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = fortaleza_obs, size = 1.2, colour = "black") +
  scale_colour_manual(values = c("Observed (SINAN)" = "black",
                                 "Predicted (SEIR)" = "#f4a582")) +
  scale_x_continuous(breaks = seq(0, 52, 10)) +
  labs(
    x = "Epidemiological week (2022)",
    y = "Reported CHIKV cases",
    colour = NULL,
    title = "Fortaleza 2022: predicted vs. observed weekly reported cases"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.major = element_line(linetype = "dotted", colour = "grey80"),
    panel.grid.minor = element_blank()
  )

# # ---- Quick fit diagnostics ---------------------------------------
# fit_summary <- combined |>
#   pivot_wider(names_from = source, values_from = cases) |>
#   rename(observed = `Observed (SINAN)`, predicted = `Predicted (SEIR)`)
# 
# cat("\n--- Fit summary ---\n")
# cat("Observed total: ",  sum(fit_summary$observed),  "\n")
# cat("Predicted total:",  round(sum(fit_summary$predicted)), "\n")
# cat("Ratio (pred/obs):", round(sum(fit_summary$predicted) / sum(fit_summary$observed), 2), "\n")
# cat("Observed peak week: ", fit_summary$week[which.max(fit_summary$observed)],
#     " (", max(fit_summary$observed), " cases)\n", sep = "")
# cat("Predicted peak week:", fit_summary$week[which.max(fit_summary$predicted)],
#     " (", round(max(fit_summary$predicted)), " cases)\n", sep = "")