# ============================================================
# Build the MAYV transmission seasonality envelope from CHIRPS rainfall.
#
# SOURCE: rainfall_caldas.csv (Caldas Novas rows extracted from the 129 MB HDX
# CHIRPS workbook). We use rfh_avg = the FIXED long-term climatological normal
# rainfall per dekad (identical across years; NOT a trailing mean), so the envelope
# represents a CLIMATOLOGICALLY TYPICAL year and has NO gap for the un-observed
# future weeks (mid-2026). This suits a HYPOTHETICAL MAYV outbreak: we pair it with
# the typical seasonal rhythm, not one year's realised weather.
#
# WINDOW: 2025-W40 -> 2026-W38 (52 epi weeks; 2025 carries an epi-week 53, so
# 2025-W40..W53 = 14 wks and 2026-W01..W38 = 38 wks).
#
# LAG: envelope = normalised rainfall shifted by LAG_WEEKS. Set to 0 by DEFAULT.
# Justification (see the CHIKV timing analysis): the envelope IS beta, and comparing
# CLIMATOLOGICAL rainfall to the fitted CHIKV beta_t (the empirical transmission
# signal, since we assume CHIKV & MAYV share the vector season) gives cross-corr 0.90
# at lag 0, falling monotonically -> big lags are ruled out, ~0 fits best. The ~10-week
# gap you SEE between rain onset and case take-off is the SEIR's own seed->outbreak
# GROWTH time, not a transmission lag; the model reproduces it, so we must NOT bake it
# into the envelope (that would double-count). A small lag (<=3 wk) is within noise;
# change LAG_WEEKS to explore it.
#
# FORM: linear in rainfall, normalised to mean 1 (drop-in for base_beta = R0*gamma*season).
# A saturating/optimum transform is a later refinement.
# ============================================================
setwd("/Users/chloelee/Documents/R/summer_project")
suppressMessages({library(dplyr); library(ggplot2)})

LAG_WEEKS <- 0                                   # rainfall -> transmission lag (weeks)

# ------------------------------------------------------------
# 1. Window: 52 epi weeks 2025-W40 -> 2026-W38, with MMWR-Sunday midpoints
# ------------------------------------------------------------
win <- rbind(data.frame(Year = 2025, week = 40:53),
             data.frame(Year = 2026, week = 1:38))
win$week_index <- seq_len(nrow(win))
win$week_label <- sprintf("%d-W%02d", win$Year, win$week)
stopifnot(nrow(win) == 52)
# MMWR W01 Sundays: 2025 -> 2024-12-29, 2026 -> 2025-12-28. Use the week midpoint (+3 d).
w01 <- ifelse(win$Year == 2025, as.Date("2024-12-29"), as.Date("2025-12-28"))
win$date <- as.Date(w01, origin = "1970-01-01") + (win$week - 1) * 7 + 3

# ------------------------------------------------------------
# 2. Climatological rainfall grid (rfh_avg by dekad-of-year), spanning the window
# ------------------------------------------------------------
r <- read.csv("rainfall_caldas.csv"); r$date <- as.Date(r$date)
r$md <- format(r$date, "%m-%d")
clim <- r %>% group_by(md) %>% summarise(rfh_avg = first(rfh_avg), .groups = "drop")
# Place each dekad-of-year at its ~mid-date in 2025 and 2026 (margin covers the window
# edges and any small negative lag shift).
grid <- do.call(rbind, lapply(c(2025, 2026), function(y)
  data.frame(date = as.Date(paste0(y, "-", clim$md)) + 4, rain = clim$rfh_avg)))
grid <- grid[order(grid$date), ]

# ------------------------------------------------------------
# 3. Interpolate to weekly, apply lag, normalise to mean 1
# ------------------------------------------------------------
win$rain   <- approx(grid$date, grid$rain, xout = win$date - LAG_WEEKS * 7, rule = 2)$y
season     <- win$rain / mean(win$rain)          # mean-1 envelope
win$season <- season
stopifnot(abs(mean(season) - 1) < 1e-9, length(season) == 52)

cat(sprintf("Rainfall envelope: lag %d wk | peak %s | range [%.2f, %.2f] | mean %.3f\n",
            LAG_WEEKS, win$week_label[which.max(season)], min(season), max(season), mean(season)))

# ------------------------------------------------------------
# 4. Save (drop-in for the old caldas_beta_season.rds) + diagnostic plot
# ------------------------------------------------------------
saveRDS(season, "caldas_rain_season.rds")
write.csv(win[, c("week_index", "week_label", "rain", "season")],
          "caldas_rain_season.csv", row.names = FALSE)

year_break <- 14.5   # between 2025-W53 (idx 14) and 2026-W01 (idx 15)
tick_idx   <- c(5, 10, 15, 22, 34, 46)
x_ticks    <- data.frame(week_index = tick_idx, week_label = win$week_label[tick_idx])
p <- ggplot(win, aes(week_index, season)) +
  annotate("rect", xmin = min(which(season >= 1)) - 0.5, xmax = max(which(season >= 1)) + 0.5,
           ymin = -Inf, ymax = Inf, fill = "#cfe6f2", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey50") +
  geom_line(colour = "#2c7fb8", linewidth = 1) +
  scale_x_continuous(breaks = x_ticks$week_index, labels = x_ticks$week_label) +
  labs(x = "Epi week", y = "Transmission envelope (mean = 1)",
       title = "Caldas Novas MAYV seasonal envelope from CHIRPS climatological rainfall",
       subtitle = sprintf("rfh_avg, dekad->week, lag %d wk, normalised; shaded = above-mean (wet)", LAG_WEEKS)) +
  theme_bw(12) + theme(plot.title = element_text(face = "bold", hjust = 0.5),
                       plot.subtitle = element_text(hjust = 0.5, size = 9),
                       panel.grid.minor = element_blank())
ggsave("caldas_rain_season.png", p, width = 8, height = 4.2, dpi = 120)
cat("Saved caldas_rain_season.rds, caldas_rain_season.csv, caldas_rain_season.png\n")
print(win[, c("week_index", "week_label", "rain", "season")], row.names = FALSE)
