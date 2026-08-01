# ============================================================
# combined_outputs.R -- every figure that spans BOTH pathogens.
#
# Each of the four pathogen scripts owns its own single-pathogen figures; anything that
# puts Chikungunya next to Mayaro lives here, so there is one place to change shared
# styling and one place that knows how the two models line up.
#
# Produces:
#   1. CHIKV_MAYV_epicurves.png       epidemic curves, side by side
#   2. combined_residual_burden.png   burden as % of no vaccination
#   3. CHIKV_MAYV_owsa.png            tornado, CHIKV over MAYV
#   4. combined_master.png            A epicurves / B DALYs / C deaths / D healthcare cost
#
# Pure presentation: reads the .rds files the model scripts write and does no SEIR or
# Monte Carlo of its own, so every value here shares the engines' per-draw propagation.
#
# Panel geometry. The two blocks are sized for what they have to show, NOT to line up
# with each other: the epidemic curves get equal halves of the width, while the burden
# blocks are only two bars wide and are inset so they do not sprawl. A therefore does not
# share a column grid with B-D, so each block is labelled Chikungunya / Mayaro in its own
# right.
#
# Mayaro has no deaths row: MAYV_ZERO_DEATHS = TRUE fixes the MAYV CFR at 0, so the
# ratio is undefined rather than zero. The facet level is kept (drop = FALSE) to hold the
# row grid aligned, then the panel grob is deleted so no box or gridlines are drawn.
#
# Run order:  CHIKV_ca_engine.R -> CHIKV_ca_costs.R -> CHIKV_ca_outputs.R
#             MAYV_ca_engine.R  -> MAYV_ca_costs.R  -> MAYV_ca_outputs.R
#             (CHIKV_ca_owsa.R, MAYV_ca_owsa.R for figure 3)  ->  this file
# ============================================================
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})

# Which MAYV scenario to show. Reads the scenario-TAGGED results file so the figures
# never depend on which scenario happened to run last.
if (!exists("MAYV_EPI_SCENARIO")) MAYV_EPI_SCENARIO <- "high"   # "high" R0 2.50 | "low" R0 1.20

ARMS      <- c("Disease-blocking", "Disease + infection blocking")
OUT_LV    <- c("Cumulative DALYs", "Cumulative deaths", "Healthcare cost")
FILL      <- c("No vaccination" = "grey60", "Vaccination" = "#4e79a7")
scen_cols <- c("No vaccination" = "grey55", "Disease-blocking" = "#4393c3",
               "Disease + infection blocking" = "#d6604d")

need1 <- function(f, how) if (!file.exists(f)) stop("missing ", f, " -- ", how) else f
G <- readRDS(need1("CHIKV_ca_engine_results.rds", "run CHIKV_ca_engine.R"))
M <- readRDS(need1(sprintf("MAYV_ca_engine_results_%s.rds", MAYV_EPI_SCENARIO),
                   sprintf("run MAYV_ca_lhs.R + MAYV_ca_engine.R with R0_SCENARIO <- '%s'",
                           MAYV_EPI_SCENARIO)))
chik <- readRDS(need1("CHIKV_ca_residual_burden.rds", "run CHIKV_ca_outputs.R"))
mayv <- readRDS(need1("MAYV_ca_residual_burden.rds",  "run MAYV_ca_outputs.R"))
stopifnot(!is.null(M$wk_base))                 # re-run MAYV_ca_engine.R to store curves

T_sim <- G$T_sim
mayv_lab <- sprintf("Mayaro (fixed R0 = %.2f)", M$R0_fixed)

# ------------------------------------------------------------
# 1. Epidemic curves.
# CHIKV carries both vaccine arms; MAYV is disease-blocking only, so it has two curves.
# MAYV uses a FIXED peak R0 per scenario, so every draw is the same transmission regime
# and ALL draws are plotted -- there is no take-off conditioning to apply.
# ------------------------------------------------------------
qband <- function(mat, rows, lab, measure) {
  q <- apply(mat[rows, , drop = FALSE], 2, quantile, c(.025, .5, .975), na.rm = TRUE)
  data.frame(week = seq_len(ncol(mat)), lo = q[1,], med = q[2,], hi = q[3,],
             scenario = lab, measure = measure)
}
pair <- function(mat, rho, rows, lab)
  rbind(qband(mat, rows, lab, "True symptomatic"),
        qband(mat * rho, rows, lab, "Reported"))

ch_rows <- seq_len(nrow(G$wk_symp[["No vaccine (baseline)"]]))
mv_rows <- seq_len(M$N_DRAWS)                  # ALL draws: R0 fixed, so no conditioning
epi_ch <- rbind(
  pair(G$wk_symp[["No vaccine (baseline)"]],                    G$rho_i, ch_rows, "No vaccination"),
  pair(G$wk_symp[["pre-outbreak | Disease-blocking"]],          G$rho_i, ch_rows, ARMS[1]),
  pair(G$wk_symp[["pre-outbreak | Disease + infection blocking"]], G$rho_i, ch_rows, ARMS[2]))
epi_mv <- rbind(pair(M$wk_base, M$rho_draw, mv_rows, "No vaccination"),
                pair(M$wk_vacc, M$rho_draw, mv_rows, ARMS[1]))
for (nm in c("epi_ch", "epi_mv")) {
  d <- get(nm); d$measure  <- factor(d$measure, levels = c("True symptomatic", "Reported"))
  d$scenario <- factor(d$scenario, levels = names(scen_cols)); assign(nm, d)
}

pct_str <- function(b, v) { p <- 100*(b - v)/b; q <- quantile(p, c(.5, .025, .975), na.rm = TRUE)
  sprintf("%.1f %% (95%% UI %.1f-%.1f %%)", q[1], q[2], q[3]) }
cb <- G$per_draw[["No vaccine (baseline)"]][, "symptomatic"]
mb <- M$per_draw[["No vaccine (baseline)"]][mv_rows, "symptomatic"]
ann_ch <- paste0("% Reduction in predicted symptomatic cases\nDisease blocking only: ",
                 pct_str(cb, G$per_draw[["pre-outbreak | Disease-blocking"]][, "symptomatic"]),
                 "\nDisease & infection blocking: ",
                 pct_str(cb, G$per_draw[["pre-outbreak | Disease + infection blocking"]][, "symptomatic"]))
ann_mv <- paste0("% Reduction in predicted symptomatic cases\nDisease blocking only: ",
                 pct_str(mb, M$per_draw[[M$vac_name]][mv_rows, "symptomatic"]))

# Rollout band = DOSING START (campaign week + median deployment delay) through the end
# of the ~1/delivery-rate week rollout. Ixchiq is deployed once, so both pathogens use
# the same definition and the two bands coincide.
ch_dose <- G$timings[["pre-outbreak"]] + 2                       # median delay 2 wk
mv_dose <- if (!is.null(M$dose_start))   M$dose_start            else ch_dose
mv_len  <- if (!is.null(M$deliv_median)) round(1/M$deliv_median) else 10

epicurve <- function(d, strip, ann, xmin, xmax, ylab, ann_size = 4) {
  d$panel <- factor(strip)
  ggplot(d, aes(week, med, colour = scenario, fill = scenario)) +
    annotate("rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
             fill = "#3a7d3a", alpha = .12) +
    geom_area(data = subset(d, measure == "True symptomatic"),
              position = "identity", alpha = .25, colour = NA) +
    geom_line(data = subset(d, measure == "Reported"),
              aes(linetype = measure), linewidth = .8) +
    annotate("text", x = -Inf, y = Inf, label = ann, hjust = -0.02, vjust = 1.15,
             size = ann_size, lineheight = 1.05) +
    expand_limits(y = 0) +
    facet_wrap(~ panel) +
    scale_colour_manual(values = scen_cols, aesthetics = c("colour", "fill"), drop = FALSE) +
    scale_linetype_manual(values = c("Reported" = "dotted"),
                          labels = c("Reported" = "Reported symptomatic cases"), name = NULL) +
    guides(fill   = guide_legend(order = 1, nrow = 2, byrow = TRUE,
                                 override.aes = list(linetype = 0, alpha = .55)),
           colour = guide_legend(order = 1, nrow = 2, byrow = TRUE,
                                 override.aes = list(linetype = 0)),
           linetype = guide_legend(order = 2, override.aes = list(colour = "grey25"))) +
    scale_x_continuous(breaks = c(1, seq(10, T_sim, by = 10))) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Week (index, 1 = 2025-W24)", y = ylab, colour = NULL, fill = NULL) +
    theme_bw(11) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 10),
          # the four keys overrun an 11in figure on one line, so the scenario guide
          # wraps to two rows rather than being shrunk until it is unreadable
          legend.text = element_text(size = 9), legend.key.size = unit(.9, "lines"),
          legend.justification = "center",
          panel.grid.minor = element_blank())
}

pE_c <- epicurve(epi_ch, "Chikungunya", ann_ch, ch_dose, ch_dose + 10, "Predicted symptomatic cases")
# the MAYV panel is a third of the width, so its annotation needs a smaller size
pE_m <- epicurve(epi_mv, mayv_lab, ann_mv, mv_dose, mv_dose + mv_len, NULL, ann_size = 3)

p_epi <- pE_c + pE_m + plot_layout(ncol = 2, widths = c(2, 1), guides = "collect") &
  theme(legend.position = "bottom")
ggsave("CHIKV_MAYV_epicurves.png", p_epi, width = 11, height = 4.8, dpi = 130)
cat(sprintf("Saved CHIKV_MAYV_epicurves.png (MAYV = '%s', fixed R0 = %.2f, all %d draws).\n",
            MAYV_EPI_SCENARIO, M$R0_fixed, M$N_DRAWS))

# ------------------------------------------------------------
# 2. Residual burden as a % of no vaccination.
# One block per outcome so each can carry its own panel letter; the CHIKV block holds
# both arms and the MAYV block one, at widths 2:1.
# ------------------------------------------------------------
chik$outcome <- factor(as.character(chik$outcome), levels = OUT_LV)
mayv$outcome <- factor(as.character(mayv$outcome), levels = OUT_LV)   # deaths kept as an
mayv$arm     <- factor(ARMS[1], levels = ARMS[1])                     # empty level
chik$arm     <- factor(as.character(chik$arm), levels = ARMS)

burden <- function(d, strip, ylab, x_axis, y_axis, tag) {
  ggplot(d, aes(scenario, med, fill = scenario)) +
    geom_col(width = .6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, linewidth = .35) +
    facet_grid(~ arm, drop = FALSE) +
    scale_fill_manual(values = FILL, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 105), breaks = seq(0, 100, 25)) +
    labs(tag = tag, x = NULL, y = ylab) +
    theme_bw(11) +
    theme(plot.tag = element_text(face = "bold", size = 13),
          plot.tag.position = c(0, 1),
          axis.title.y = element_text(size = 10),
          # top margin gives each panel letter its own band, so B, C and D stay legible
          # instead of being squeezed against the row above
          plot.margin = margin(t = 22, r = 5.5, b = 5.5, l = 5.5),
          axis.text.x  = if (x_axis) element_text(size = 8.5) else element_blank(),
          axis.ticks.x = if (x_axis) element_line() else element_blank(),
          axis.text.y  = if (y_axis) element_text() else element_blank(),
          axis.ticks.y = if (y_axis) element_line() else element_blank(),
          strip.text   = if (is.null(strip)) element_blank()
                         else element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank(),
          # the bars are already named on the x axis, so a fill legend adds nothing
          legend.position = "none")
}

# strip = NULL blanks the arm strip on the lower blocks: C and D repeat B's columns.
row_burden <- function(o, tag, ylab = NULL, strip = FALSE, x_axis = FALSE) {
  st <- if (strip) TRUE else NULL
  c_blk <- burden(subset(chik, outcome == o), st, ylab, x_axis, TRUE, tag)
  if (o == "Cumulative deaths") return(list(c = c_blk, m = NULL))     # MAYV: CFR fixed at 0
  list(c = c_blk, m = burden(subset(mayv, outcome == o), st, NULL, x_axis, FALSE, NULL))
}

# Disease names head the whole B-D block, above the panel letters. They sit on their own
# void row rather than as titles on row B so that they read before the B, and so C and D
# inherit them without the names being repeated three times.
head_lab <- function(txt) ggplot() + labs(title = txt) + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = .5,
                                  margin = margin(t = 0, b = 2)),
        plot.margin = margin(0, 0, 0, 0))

YLAB <- "% of no vaccination"
rB <- row_burden("Cumulative DALYs",  "B", "Cumulative DALYs", strip = TRUE)
rC <- row_burden("Cumulative deaths", "C", "Cumulative deaths")
rD <- row_burden("Healthcare cost",   "D", "Healthcare cost", x_axis = TRUE)

# standalone version of the burden figure, three outcomes stacked
burden_rows <- (head_lab("Chikungunya") + head_lab("Mayaro") +
                  plot_layout(widths = c(2, 1))) /
               (rB$c + rB$m + plot_layout(widths = c(2, 1))) /
               (rC$c + plot_spacer() + plot_layout(widths = c(2, 1))) /
               (rD$c + rD$m + plot_layout(widths = c(2, 1))) +
               plot_layout(heights = c(.01, 1, 1, 1))

# One rotated y title for the whole B-D block. Put on a single row it would be taller
# than that row and get clipped, so it is attached to the wrapped block instead.
with_ylab <- function(blk) wrap_elements(patchworkGrob(blk)) +
  labs(tag = YLAB) +
  theme(plot.tag = element_text(size = 11, angle = 90), plot.tag.position = "left",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

p_burden <- with_ylab(burden_rows)
ggsave("combined_residual_burden.png", p_burden, width = 9.5, height = 9, dpi = 150)
cat("Saved combined_residual_burden.png (B DALYs / C deaths / D healthcare cost).\n")

# ------------------------------------------------------------
# 3. One-way sensitivity: CHIKV (A) over MAYV (B).
# MAYV has only the disease-blocking arm, so its panel occupies the left-hand column
# with the right-hand cell left empty.
# ------------------------------------------------------------
if (!all(file.exists("CHIKV_ca_owsa.rds", "MAYV_ca_owsa.rds"))) {
  cat("Skipped CHIKV_MAYV_owsa.png (run CHIKV_ca_owsa.R and MAYV_ca_owsa.R).\n")
} else {
  CO <- readRDS("CHIKV_ca_owsa.rds"); MO <- readRDS("MAYV_ca_owsa.rds")

  ca <- CO$owsa; ca$val <- ca$symptomatic
  bA <- setNames(CO$base$symptomatic, CO$base$arm)
  ca$base <- bA[ca$arm]; ca$arm <- factor(ca$arm, levels = ARMS)
  swA <- ca |> group_by(parameter) |> summarise(s = max(abs(val - base)), .groups = "drop")
  ca$parameter <- factor(ca$parameter, levels = swA$parameter[order(swA$s)])
  blA <- data.frame(arm = factor(names(bA), levels = ARMS), base = as.numeric(bA))

  mb <- MO$owsa; mb$val <- mb$averted; mb$base <- MO$base$averted
  mb$arm <- factor(ARMS[1], levels = ARMS[1])
  swB <- mb |> group_by(parameter) |> summarise(s = max(abs(val - base)), .groups = "drop")
  mb$parameter <- factor(mb$parameter, levels = swB$parameter[order(swB$s)])
  blB <- data.frame(arm = factor(ARMS[1], levels = ARMS[1]), base = MO$base$averted)

  tor <- function(d, bl, ttl, legend) {
    ggplot(d, aes(y = parameter)) +
      geom_vline(data = bl, aes(xintercept = base), linetype = "dashed", colour = "grey45") +
      geom_segment(aes(x = base, xend = val, yend = parameter, colour = bound),
                   linewidth = 5.5, alpha = .85) +
      facet_wrap(~ arm, nrow = 1, scales = "free_x") +
      scale_colour_manual(values = c(lower = "#d6604d", upper = "#4393c3"),
                          labels = c(lower = "Lower", upper = "Upper"), name = "Bound") +
      scale_x_continuous(labels = scales::comma) +
      labs(title = ttl, x = "Symptomatic cases averted", y = NULL) +
      theme_bw(11) +
      theme(text = element_text(size = 14),
            legend.position = if (legend) "bottom" else "none",
            plot.title = element_text(face = "bold", size = 12),
            strip.text = element_text(face = "bold", size = 10),
            panel.grid.minor = element_blank())
  }
  g <- tor(ca, blA, "A   Chikungunya", FALSE) + tor(mb, blB, "B   Mayaro", TRUE) +
    plot_layout(design = "AA\nB#", heights = c(1, 1.05))   # A spans, B is left column only
  ggsave("CHIKV_MAYV_owsa.png", g, width = 11, height = 8.4, dpi = 130)
  cat("Saved CHIKV_MAYV_owsa.png (A = CHIKV, B = MAYV).\n")
}

# ------------------------------------------------------------
# 4. Master figure: A epidemic curves, B DALYs, C deaths, D healthcare cost.
# The disease labels sit once, on the strips over row A; B-D inherit those columns.
# ------------------------------------------------------------
pA_c <- pE_c + labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 13), plot.tag.position = c(0, 1))
# Both epicurve panels carry identical keys, so MAYV's copy is dropped at the SCALE
# level -- a theme(legend.position = "none") would be undone by the shared `&` theme.
pA_m <- pE_m + guides(fill = "none", colour = "none", linetype = "none")

# equal halves: the Mayaro curve needs as much room to be read as the Chikungunya one
row_A <- pA_c + pA_m + plot_layout(widths = c(1, 1)) &
  theme(legend.position = "bottom", plot.margin = margin(t = 5.5, r = 5.5, b = 2, l = 5.5))

# the burden panels are two bars wide, so the block is inset rather than stretched
row_BCD <- plot_spacer() + with_ylab(burden_rows) + plot_spacer() +
  plot_layout(widths = c(.08, 1, .08))

master <- row_A / row_BCD + plot_layout(heights = c(1.3, 2.7))

ggsave("combined_master.png", master, width = 11, height = 12.5, dpi = 150)
ggsave("combined_master.pdf", master, width = 11, height = 12.5)
cat("Saved combined_master.png / .pdf (A epicurves, B DALYs, C deaths, D healthcare cost).\n")
