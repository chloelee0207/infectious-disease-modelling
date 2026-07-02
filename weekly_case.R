setwd("/Users/chloelee/Documents/R/summer_project")

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

raw <- read_excel("weekly_case.xlsx", sheet = "weekly_all")

caldas <- raw |>
  filter(Code == 520450 & Year %in% c(2025, 2026)) |>
  pivot_longer(
    cols      = starts_with("Week"),
    names_to  = "epi_week",
    values_to = "tot_cases"
  ) |>
  mutate(
    epi_week  = as.integer(sub("Week ", "", epi_week)),
    tot_cases = ifelse(is.na(tot_cases), 0, tot_cases)
  ) |>
  filter((Year == 2025 & epi_week >= 23) | (Year == 2026 & epi_week <= 22)) |>
  arrange(Year, epi_week) |>
  mutate(
    week_index = row_number(),                                 
    week_label = paste0(Year, "-W", sprintf("%02d", epi_week))
  ) |>
  dplyr::select(week_index, week_label, Year, epi_week, tot_cases)

observed_cases <- caldas$tot_cases
T_weeks <- length(observed_cases)
stopifnot(T_weeks == 52, sum(observed_cases) == 8085)
cat("Caldas Novas outbreak (2025-W23 to 2026-W22) reported CHIKV cases:",
    sum(observed_cases), "\n")
cat("Peak week:", caldas$week_label[which.max(observed_cases)],
    "with", max(observed_cases), "cases\n")

x_ticks <- caldas |>
  filter((Year == 2025 & epi_week %in% c(30, 40, 50)) |
         (Year == 2026 & epi_week %in% c(10, 20))) |>
  mutate(label = as.character(epi_week))   # plain epi-week numbers (reset at the year boundary)

year_break <- mean(c(
  max(caldas$week_index[caldas$Year == 2025]),
  min(caldas$week_index[caldas$Year == 2026])
))

# plot
caldas_weekly <- ggplot(caldas, aes(x = week_index, y = tot_cases)) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  annotate("text", x = year_break, y = 0,
           label = "2026", angle = 90, vjust = -0.4, hjust = -11,
           fontface = "bold", size = 3.6, colour = "grey40") +
  geom_line(linewidth = 0.6) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$label) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week",
       y = "Reported CHIKV cases",
       title = "Weekly CHIKV cases in Caldas Novas (2025-2026)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_line(linetype = "dotted", colour = "grey80"),
    panel.grid.minor = element_blank()
  )

print(caldas_weekly)
# ggsave("caldas_weekly.png", caldas_weekly, width = 8, height = 4.5, dpi = 120)
