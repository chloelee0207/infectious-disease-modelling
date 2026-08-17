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
#   4. combined_master.png            A epicurves / B DALYs / C deaths / D hospitalisation cost
#   5. combined_per_100k_doses.xlsx   benefit per 100,000 doses, CHIKV beside MAYV
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
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(writexl)
                  library(readxl)})

# Which MAYV scenario to show. Reads the scenario-TAGGED results file so the figures
# never depend on which scenario happened to run last.
if (!exists("MAYV_EPI_SCENARIO")) MAYV_EPI_SCENARIO <- "high"   # "high" R0 2.1-2.9 | "low" R0 1.1-1.3

ARMS      <- c("Disease-blocking", "Disease + infection blocking")
OUT_LV    <- c("Cumulative DALYs", "Cumulative deaths", "Hospitalisation cost")  # as stored in the .rds
# Strip text for the B/C/D rows -- EDIT HERE to rename them. Must stay in OUT_LV order;
# the .rds names above are the lookup keys and must not be changed.
OUT_LAB   <- c("DALYs", "Deaths", "Hospitalisation costs")
FILL      <- c("No vaccination" = "grey60", "Vaccination" = "#4e79a7")
scen_cols <- c("No vaccination" = "grey55", "Disease-blocking" = "#4393c3",
               "Disease + infection blocking" = "#d6604d")

# ------------------------------------------------------------
# FONT SIZES -- every text size in the combined figures, in one place. Tune here.
# All values are points EXCEPT epi_annot_*, which is the annotate()/geom_text() scale
# (roughly points / 2.845), because that text is drawn inside the panel as data.
# ------------------------------------------------------------
FS <- list(
  panel_letter   = 15,  # the bold A / B / C / D (was 14 for A, 13 for B-D; unified)
  disease_header = 13,  # "Chikungunya" / "Mayaro" headings above the B-D block
  block_ylab     = 13,  # rotated "Cumulative burden (% of no vaccination)"

  # ---- panel A, the epidemic curves
  epi_base       = 12,  # theme_bw() base; anything not named below inherits from it
  epi_axis_text  = 12,  # week numbers and case counts
  epi_axis_title = 13,  # "Week (index, 1 = 2025-W24)", "Predicted symptomatic cases"
  epi_strip      = 12,  # grey strips: "Chikungunya", "Mayaro (fixed R0 = ...)"
  epi_annot_chik = 4,   # in-panel "% Reduction in predicted symptomatic cases", CHIKV
  epi_annot_mayv = 4,   # the same text in the MAYV panel  (annotate scale, see above)
  epi_legend     = 11,  # legend keys under panel A

  # ---- panels B, C, D, the burden bars
  bur_base       = 14,  # theme_bw() base for the burden blocks
  bur_axis_x     = 12,  # "No vaccination" / "Vaccination"
  bur_axis_y     = 12,  # 0% - 100%
  bur_strip_x    = 12,  # top grey strips: the vaccine arm names
  bur_strip_y    = 12,  # right grey strips: DALYs / Deaths / Hospitalisation costs

  # ---- the separate tornado figure, CHIKV_MAYV_owsa.png
  owsa_base      = 14,
  owsa_title     = 12,
  owsa_strip     = 10
)

need1 <- function(f, how) if (!file.exists(f)) stop("missing ", f, " -- ", how) else f
G <- readRDS(need1("CHIKV_ca_engine_results.rds", "run CHIKV_ca_engine.R"))
M <- readRDS(need1(sprintf("MAYV_ca_engine_results_%s.rds", MAYV_EPI_SCENARIO),
                   sprintf("run MAYV_ca_lhs.R + MAYV_ca_engine.R with R0_SCENARIO <- '%s'",
                           MAYV_EPI_SCENARIO)))
chik <- readRDS(need1("CHIKV_ca_residual_burden.rds", "run CHIKV_ca_outputs.R"))
mayv <- readRDS(need1("MAYV_ca_residual_burden.rds",  "run MAYV_ca_outputs.R"))
stopifnot(!is.null(M$wk_base))                 # re-run MAYV_ca_engine.R to store curves

T_sim <- G$T_sim
mayv_lab <- if (isTRUE(M$R0_sampled))
  sprintf("Mayaro (R0 %.1f-%.1f)", M$R0_lo, M$R0_hi) else
  sprintf("Mayaro (fixed R0 = %.2f)", M$R0_fixed)

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

# ------------------------------------------------------------
# Figure saving. Every manuscript figure is written TWICE:
#   <name>.png         screen resolution, for drafting and quick review
#   <name>_600dpi.png  600 dpi, for submission
# Both are kept deliberately: a Word document with several 600-dpi images gets very large,
# so the screen-resolution copy is the fallback if the .docx becomes unwieldy. The figures
# are identical in physical size (inches) and layout -- only the pixel density differs, so
# they are interchangeable in the document.
# ------------------------------------------------------------
DPI_PUB <- 600
save_fig <- function(file, plot, width, height, dpi) {
  ggsave(file, plot, width = width, height = height, dpi = dpi)
  hi <- sub("[.]png$", "_600dpi.png", file)
  ggsave(hi, plot, width = width, height = height, dpi = DPI_PUB)
  cat(sprintf("  saved %s (%d dpi, %.1f MB) and %s (%d dpi, %.1f MB)\n",
              file, dpi, file.size(file)/1e6, hi, DPI_PUB, file.size(hi)/1e6))
  invisible(NULL)
}

epicurve <- function(d, strip, ann, xmin, xmax, ylab, ann_size = FS$epi_annot_chik) {
  d$panel <- factor(strip)
  # Annotation placement. Two traps here:
  #   * vjust > 1 offsets the block by a fraction of its OWN height, so the 2-line Mayaro
  #     label sat higher than the 3-line chikungunya one. Use vjust = 1 (tops aligned).
  #   * vjust = 1 at y = Inf then sits flush against the frame. Deriving a y from
  #     max(d$hi) does NOT fix it: d holds several measures/scenarios, so its maximum is
  #     not the plotted range and the axis blows up.
  # Simplest robust fix: keep y = Inf / vjust = 1 and prepend a blank line. Line height is
  # identical in both panels, so both blocks drop by exactly one line -> equal gap, and it
  # cannot interact with the data range.
  TOP_EXP <- 0.20                 # y-scale headroom above the data, for the annotation
  ann <- paste0("\n", ann)
  ggplot(d, aes(week, med, colour = scenario, fill = scenario)) +
    annotate("rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
             fill = "#3a7d3a", alpha = .12) +
    geom_area(data = subset(d, measure == "True symptomatic"),
              position = "identity", alpha = .25, colour = NA) +
    geom_line(data = subset(d, measure == "Reported"),
              aes(linetype = measure), linewidth = .8) +
    # Anchored at week 1 rather than x = -Inf: with -Inf the hjust offset scales with
    # the text width, so the two panels get different left margins. A data-unit anchor
    # gives the same small gap in both.
    # vjust = 1 (NOT 1.3): vjust > 1 offsets the block by a FRACTION OF ITS OWN HEIGHT,
    # so the 2-line Mayaro label was pushed down less than the 3-line chikungunya one and
    # sat visibly higher. vjust = 1 tops-aligns both regardless of line count; the
    # headroom comes from the y-scale expansion below instead.
    annotate("text", x = 1, y = Inf, label = ann, hjust = 0, vjust = 1,
             size = ann_size, lineheight = 1.15) +
    expand_limits(y = 0) +
    facet_wrap(~ panel) +
    scale_colour_manual(values = scen_cols, aesthetics = c("colour", "fill"), drop = FALSE) +
    scale_linetype_manual(values = c("Reported" = "dotted"),
                          labels = c("Reported" = "Reported symptomatic cases"), name = NULL) +
    guides(fill   = guide_legend(order = 1, nrow = 1,
                                 override.aes = list(linetype = 0, alpha = .55)),
           colour = guide_legend(order = 1, nrow = 1,
                                 override.aes = list(linetype = 0)),
           linetype = guide_legend(order = 2, nrow = 1,
                                   override.aes = list(colour = "grey25"))) +
    scale_x_continuous(breaks = c(1, seq(10, T_sim, by = 10))) +
    # extra headroom at the top for the (now top-aligned) annotation
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0.05, TOP_EXP))) +
    labs(x = "Week (index, 1 = 2025-W24)", y = ylab, colour = NULL, fill = NULL) +
    theme_bw(FS$epi_base) +
    theme(legend.position = "bottom",
          # the two guides STACK vertically: the three vaccination scenarios on one row,
          # then "Reported symptomatic cases" on its own row underneath
          legend.box = "vertical",
          legend.box.just = "center",
          legend.spacing.y = unit(1, "pt"),
          strip.text = element_text(face = "bold", size = FS$epi_strip),
          legend.text = element_text(size = FS$epi_legend),
          legend.key.size = unit(.9, "lines"),
          axis.text  = element_text(size = FS$epi_axis_text),
          axis.title = element_text(size = FS$epi_axis_title),
          legend.justification = "center",
          panel.grid.minor = element_blank())
}

# Wrap guard only. At annotation size 4 the longest line (57 characters) fits a
# half-width panel, so 60 leaves the natural three-line / two-line structure intact and
# only catches a string that grows unexpectedly. strwrap breaks on spaces, so a wrap can
# never land inside a number.
wrap_ann <- function(x, w) paste(vapply(strsplit(x, "\n")[[1]],
  function(l) paste(strwrap(l, width = w), collapse = "\n"), character(1)), collapse = "\n")

pE_c <- epicurve(epi_ch, "Chikungunya", wrap_ann(ann_ch, 60), ch_dose, ch_dose + 10,
                 "Predicted symptomatic cases")
pE_m <- epicurve(epi_mv, mayv_lab, wrap_ann(ann_mv, 60), mv_dose, mv_dose + mv_len, NULL)

# Equal panel widths, and ONE legend: the two panels carry identical keys, so MAYV's copy
# is dropped at the SCALE level -- a theme(legend.position = "none") would be undone by
# the shared `&` theme that positions the collected guide.
p_epi <- pE_c + (pE_m + guides(fill = "none", colour = "none", linetype = "none")) +
  # guides = "collect" lifts the legend out of the first panel and gives it its own
  # strip across the bottom of BOTH panels; without it the legend belongs to Chikungunya
  # and sits left-aligned under that panel alone.
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  # legend.box MUST be set here as well as inside epicurve(): with guides = "collect"
  # patchwork lays the collected guides out using the PATCHWORK-level theme, so a
  # per-subplot legend.box is ignored. "vertical" stacks the two guide boxes -> the three
  # vaccination scenarios on one row, "Reported symptomatic cases" on a second.
  theme(legend.position = "bottom", legend.justification = "center",
        legend.box = "vertical", legend.box.just = "center",
        legend.spacing.y = unit(1, "pt"))
save_fig("CHIKV_MAYV_epicurves.png", p_epi, width = 11, height = 4.8, dpi = 130)
cat(sprintf("Saved CHIKV_MAYV_epicurves.png (MAYV = '%s', fixed R0 = %.2f, all %d draws).\n",
            MAYV_EPI_SCENARIO, M$R0_fixed, M$N_DRAWS))

# ------------------------------------------------------------
# 2. Residual burden as a % of no vaccination.
# One block per outcome so each can carry its own panel letter; the CHIKV block holds
# both arms and the MAYV block one, at widths 2:1.
# ------------------------------------------------------------
chik$outcome <- factor(as.character(chik$outcome), levels = OUT_LV, labels = OUT_LAB)
mayv$outcome <- factor(as.character(mayv$outcome), levels = OUT_LV, labels = OUT_LAB)
mayv$arm     <- factor(ARMS[1], levels = ARMS[1])       # deaths kept as an empty level
chik$arm     <- factor(as.character(chik$arm), levels = ARMS)

burden <- function(d, strip, ylab, x_axis, y_axis, tag, row_strip = FALSE) {
  ggplot(d, aes(scenario, med, fill = scenario)) +
    geom_col(width = .6, na.rm = TRUE) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, linewidth = .35, na.rm = TRUE) +
    facet_grid(outcome ~ arm, drop = FALSE) +
    scale_fill_manual(values = FILL, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 105), breaks = seq(0, 100, 25)) +
    labs(tag = tag, x = NULL, y = ylab) +
    theme_bw(FS$bur_base) +
    theme(plot.tag = element_text(face = "bold", size = FS$panel_letter),
          plot.tag.position = c(0, 1),
          # top margin gives each panel letter its own band, so B, C and D stay legible
          # instead of being squeezed against the row above
          plot.margin = margin(t = 22, r = 5.5, b = 5.5, l = 5.5),
          axis.text.x  = if (x_axis) element_text(size = FS$bur_axis_x) else element_blank(),
          axis.ticks.x = if (x_axis) element_line() else element_blank(),
          axis.text.y  = if (y_axis) element_text(size = FS$bur_axis_y) else element_blank(),
          axis.ticks.y = if (y_axis) element_line() else element_blank(),
          strip.text.x = if (is.null(strip)) element_blank()
                         else element_text(face = "bold", size = FS$bur_strip_x),
          # the outcome strip rides on the RIGHT-hand (MAYV) block only, so it sits at
          # the far edge of the figure and cannot collide with the B/C/D letters
          strip.text.y = if (row_strip) element_text(face = "bold", size = FS$bur_strip_y,
                                                    angle = -90)
                         else element_blank(),
          panel.grid.minor = element_blank(),
          # the bars are already named on the x axis, so a fill legend adds nothing
          legend.position = "none")
}

# strip = NULL blanks the arm strip on the lower blocks: C and D repeat B's columns.
# drop = FALSE is needed for the ARM dimension (MAYV lacks the infection-blocking level),
# but it would also render all three outcome rows in every block, so the outcome factor is
# reduced to the single level this row shows.
one_row <- function(d, o) { d <- subset(d, outcome == o); d$outcome <- droplevels(d$outcome); d }

# The Mayaro deaths cell has no bars (CFR fixed at 0), but the row still needs its
# right-hand strip, so it gets a data-free block with every panel element blanked.
blank_row <- function(o) {
  d <- data.frame(outcome = factor(o), arm = factor(ARMS[1], levels = ARMS[1]),
                  scenario = factor("No vaccination", levels = names(FILL)),
                  med = NA_real_, lo = NA_real_, hi = NA_real_)
  burden(d, NULL, NULL, FALSE, FALSE, NULL, row_strip = TRUE) +
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank())
}

row_burden <- function(o, tag, strip = FALSE, x_axis = FALSE) {
  st <- if (strip) TRUE else NULL
  c_blk <- burden(one_row(chik, o), st, NULL, x_axis, TRUE, tag)
  m_blk <- if (o == OUT_LAB[2]) blank_row(o)                          # MAYV CFR is 0
           else burden(one_row(mayv, o), st, NULL, x_axis, FALSE, NULL, row_strip = TRUE)
  list(c = c_blk, m = m_blk)
}

# Disease names head the whole B-D block, above the panel letters. They sit on their own
# void row rather than as titles on row B so that they read before the B, and so C and D
# inherit them without the names being repeated three times.
head_lab <- function(txt) ggplot() + labs(title = txt) + theme_void() +
  theme(plot.title = element_text(face = "bold", size = FS$disease_header, hjust = .5,
                                  margin = margin(t = 0, b = 2)),
        plot.margin = margin(0, 0, 0, 0))

YLAB <- "Cumulative burden (% of no vaccination)"
rB <- row_burden(OUT_LAB[1], "B", strip = TRUE)
rC <- row_burden(OUT_LAB[2], "C")
rD <- row_burden(OUT_LAB[3], "D", x_axis = TRUE)

# standalone version of the burden figure, three outcomes stacked
burden_rows <- (head_lab("Chikungunya") + head_lab("Mayaro") +
                  plot_layout(widths = c(2, 1))) /
               (rB$c + rB$m + plot_layout(widths = c(2, 1))) /
               (rC$c + rC$m + plot_layout(widths = c(2, 1))) /
               (rD$c + rD$m + plot_layout(widths = c(2, 1))) +
               plot_layout(heights = c(.01, 1, 1, 1))

# One rotated y title for the whole B-D block. Put on a single row it would be taller
# than that row and get clipped, so it is attached to the wrapped block instead.
with_ylab <- function(blk) wrap_elements(patchworkGrob(blk)) +
  labs(tag = YLAB) +
  theme(plot.tag = element_text(size = FS$block_ylab, angle = 90), plot.tag.position = "left",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

p_burden <- with_ylab(burden_rows)
save_fig("combined_residual_burden.png", p_burden, width = 9.5, height = 9, dpi = 150)
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
      theme_bw(FS$owsa_base) +
      theme(text = element_text(size = FS$owsa_base),
            legend.position = if (legend) "bottom" else "none",
            plot.title = element_text(face = "bold", size = FS$owsa_title),
            strip.text = element_text(face = "bold", size = FS$owsa_strip),
            panel.grid.minor = element_blank())
  }
  g <- tor(ca, blA, "A   Chikungunya", FALSE) + tor(mb, blB, "B   Mayaro", TRUE) +
    plot_layout(design = "AA\nB#", heights = c(1, 1.05))   # A spans, B is left column only
  save_fig("CHIKV_MAYV_owsa.png", g, width = 11, height = 8.4, dpi = 130)
  cat("Saved CHIKV_MAYV_owsa.png (A = CHIKV, B = MAYV).\n")
}

# ------------------------------------------------------------
# 4. Master figure: A epidemic curves, B DALYs, C deaths, D healthcare cost.
# The disease labels sit once, on the strips over row A; B-D inherit those columns.
# ------------------------------------------------------------
pA_c <- pE_c + labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = FS$panel_letter),
        plot.tag.position = c(0, 1))
# Both epicurve panels carry identical keys, so MAYV's copy is dropped at the SCALE
# level -- a theme(legend.position = "none") would be undone by the shared `&` theme.
# Rebuilt rather than reused from pE_m so the master figure's strip is a plain "Mayaro":
# in a multi-panel figure the R0 belongs in the caption, not repeated in a panel header.
# The standalone CHIKV_MAYV_epicurves.png keeps "Mayaro (fixed R0 = x.xx)" because that
# figure is often viewed on its own, where the scenario would otherwise be ambiguous.
pA_m <- epicurve(epi_mv, "Mayaro", wrap_ann(ann_mv, 60), mv_dose, mv_dose + mv_len, NULL) +
  guides(fill = "none", colour = "none", linetype = "none")

# equal halves: the Mayaro curve needs as much room to be read as the Chikungunya one
# guides = "collect" lifts the legend out of the Chikungunya panel and gives it its own
# strip across the bottom of BOTH panels, so it is centred on the row rather than
# left-aligned under the first plot. Note that a collected guide is drawn from the
# PATCHWORK's theme, so its size must be set here with `&`, not inside epicurve().
row_A <- pA_c + pA_m + plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "center",
        # same two-row legend as the standalone epicurve figure (see section 1)
        legend.box = "vertical", legend.box.just = "center",
        legend.spacing.y = unit(1, "pt"),
        legend.text = element_text(size = FS$epi_legend), legend.key.size = unit(1.1, "lines"),
        plot.margin = margin(t = 5.5, r = 5.5, b = 2, l = 5.5))

# the burden panels are two bars wide, so the block is inset rather than stretched
row_BCD <- plot_spacer() + with_ylab(burden_rows) + plot_spacer() +
  plot_layout(widths = c(.08, 1, .08))

# Burden block gets a larger share and the figure is taller overall. The right-hand strip
# text is rotated, so its available length IS the row height -- "Hospitalisation costs" is
# the longest label and sets the requirement.
master <- row_A / row_BCD + plot_layout(heights = c(1.15, 3.1))

save_fig("combined_master.png", master, width = 13, height = 15, dpi = 150)
ggsave("combined_master.pdf", master, width = 13, height = 15)
cat("Saved combined_master.png / .pdf (A epicurves, B DALYs, C deaths, D healthcare cost).\n")


# ------------------------------------------------------------
# 5. Benefit per 100,000 doses, CHIKV beside MAYV.
# Ixchiq is deployed ONCE, so both models cost the same doses: 19,584 (CHIKV) vs 19,589
# (MAYV), a 0.02% difference from sampling alone. Normalising by doses therefore puts the
# two pathogens in a shared unit WITHOUT dividing one by the other -- the Mayaro figure is
# the additional return on doses already being bought for chikungunya, expressed as an
# absolute quantity rather than a ratio.
#
# The combined column adds the two per draw. CHIKV and MAYV draws are NOT paired (separate
# models, separate LHS designs), so one is randomly permuted before adding: that samples
# the sum correctly under the assumption the two epidemics are independent, which is
# stated rather than hidden. The combined total is conditional on BOTH outbreaks occurring
# as modelled -- the Mayaro one is hypothetical.
# ------------------------------------------------------------
per100k <- function(pd, base, o, doses) {
  r <- 1e5 * (base[, o] - pd[, o]) / doses; r[!is.finite(r) | r < 0] <- NA; r
}
OUT100 <- c(infections = "Infections", symptomatic = "Symptomatic",
            hospitalisations = "Hospitalisations", deaths = "Deaths", daly = "DALYs")
gb <- G$per_draw[["No vaccine (baseline)"]]; mbase <- M$per_draw[["No vaccine (baseline)"]]
ch_db <- G$per_draw[["pre-outbreak | Disease-blocking"]]
ch_di <- G$per_draw[["pre-outbreak | Disease + infection blocking"]]
mv_db <- M$per_draw[[M$vac_name]]
set.seed(20260802); perm <- sample.int(nrow(mv_db))     # independence pairing, reproducible

fmt100 <- function(v, d = 0) { q <- quantile(v, c(.5,.025,.975), na.rm = TRUE)
  sprintf("%s (%s-%s)", formatC(round(q[1],d), big.mark=",", format="f", digits=d),
          formatC(round(q[2],d), big.mark=",", format="f", digits=d),
          formatC(round(q[3],d), big.mark=",", format="f", digits=d)) }

rows <- lapply(names(OUT100), function(o) {
  a <- per100k(ch_db, gb, o, ch_db[, "doses"])
  b <- per100k(ch_di, gb, o, ch_di[, "doses"])
  m <- per100k(mv_db, mbase, o, mv_db[, "doses"])
  data.frame(outcome = OUT100[[o]],
             CHIKV_disease_blocking = fmt100(a),
             CHIKV_disease_and_infection = fmt100(b),
             MAYV_disease_blocking = fmt100(m),
             combined_CHIKV_db_plus_MAYV = fmt100(a + m[perm]),
             stringsAsFactors = FALSE) })

# healthcare cost, if the cost layer has been run
cc <- if (file.exists("CHIKV_ca_costs.rds")) readRDS("CHIKV_ca_costs.rds") else NULL
cm <- if (file.exists("MAYV_ca_costs.rds"))  readRDS("MAYV_ca_costs.rds")  else NULL
if (!is.null(cc) && !is.null(cm)) {
  K <- 1.47761 / 5.1257                                   # BRL2019 -> USD2026, as in the cost scripts
  # Costs are split into INPATIENT and OUTPATIENT rather than reported as one total.
  # MAYV outpatient care is deliberately not costed (no Mayaro treatment flowchart, drugs
  # or unit costs exist), so those cells are NA and read "not estimated" -- never 0, which
  # would assert that no outpatient cost is incurred.
  cost100 <- function(cp, nm, doses, comp) {
    base <- cp[["No vaccine (baseline)"]]; scen <- cp[[nm]]
    a <- if (identical(comp, "outpatient")) {
           (rowSums(base[, c("out_acute","out_subacute","out_chronic")]) -
            rowSums(scen[, c("out_acute","out_subacute","out_chronic")])) * K
         } else (base[, comp] - scen[, comp]) * K
    r <- 1e5 * a / doses; r[!is.finite(r)] <- NA; r }
  f100 <- function(v) if (all(is.na(v))) "not estimated" else fmt100(v)
  mv_nm <- names(cm$cost_pd)[2]
  for (cmp in c("hosp_inpatient", "outpatient")) {
    a <- cost100(cc$cost_pd, "pre-outbreak | Disease-blocking", ch_db[, "doses"], cmp)
    b <- cost100(cc$cost_pd, "pre-outbreak | Disease + infection blocking", ch_di[, "doses"], cmp)
    m <- cost100(cm$cost_pd, mv_nm, mv_db[, "doses"], cmp)
    comb <- if (all(is.na(m))) "CHIKV only (MAYV not estimated)" else f100(a + m[perm])
    rows[[length(rows)+1]] <- data.frame(
      outcome = if (cmp == "hosp_inpatient") "Hospitalisation cost (US$ 2026)"
                else "Outpatient cost (US$ 2026)",
      CHIKV_disease_blocking = f100(a), CHIKV_disease_and_infection = f100(b),
      MAYV_disease_blocking = f100(m), combined_CHIKV_db_plus_MAYV = comb,
      stringsAsFactors = FALSE)
  }
} else cat("NOTE: cost layers not found -- healthcare-cost row omitted from the per-100k table.\n")

per100k_tbl <- do.call(rbind, rows)
p100_notes <- data.frame(item = c("Unit", "Why doses", "Why not a ratio", "Combined column",
                                  "Deaths", "Doses", "Cost split"),
  detail = c(
  "Outcomes averted per 100,000 doses administered, median (95% UI), computed per draw.",
  "Ixchiq is deployed once, so both models consume the same doses. Doses are the shared input, which makes per-dose benefit comparable across pathogens without either model being divided by the other.",
  "Dividing Mayaro cases averted by chikungunya cases averted compares two different models on a scale with no natural zero. Per-dose benefit has a natural zero (no benefit) and is additive.",
  "CHIKV disease-blocking + MAYV, added per draw with MAYV randomly permuted, since the two models' draws are not paired. Assumes the two epidemics are independent, and is conditional on BOTH occurring as modelled -- the Mayaro outbreak is hypothetical.",
  "MAYV deaths are identically zero (MAYV_ZERO_DEATHS = TRUE fixes the CFR at 0), so the Mayaro column is 0 by construction.",
  "Costs are split into inpatient and outpatient. MAYV outpatient is NOT costed -- no Mayaro treatment flowchart, drug list or unit costs exist, and borrowing the chikungunya regimen would manufacture a figure. Those cells read 'not estimated' rather than 0, and the combined outpatient row is therefore CHIKV only.",
  sprintf("CHIKV %s, MAYV %s (median), a %.2f%% difference from sampling alone.",
          format(round(median(ch_db[, "doses"])), big.mark=","),
          format(round(median(mv_db[, "doses"])), big.mark=","),
          100*abs(median(mv_db[,"doses"])/median(ch_db[,"doses"]) - 1))),
  stringsAsFactors = FALSE)
writexl::write_xlsx(list(averted_per_100k_doses = per100k_tbl, notes = p100_notes),
                    "combined_per_100k_doses.xlsx")
cat("Saved combined_per_100k_doses.xlsx (CHIKV beside MAYV, per 100,000 doses).\n")
print(per100k_tbl, row.names = FALSE)


# ------------------------------------------------------------
# 6. Proportion still symptomatic at three months -- supplementary table.
#
# This is a PARAMETER-level table, not a Caldas Novas result: it describes the severity
# inputs the two models draw from, so it carries no dependence on outbreak size, age mix
# or vaccination. Sourced from disease_progression.xlsx so it can never drift from what
# the engines actually sample.
#
# The comparable quantity across the two pathogens is the proportion of symptomatic cases
# STILL SYMPTOMATIC AT ~3 MONTHS:
#   Chikungunya  1 - (recovered by 14 d) - (recovered 14 d - 90 d) = sum of the 6 m, 12 m
#                and 30 m rows of O'Driscoll et al. 2021, whose five rows are MARGINAL
#                shares of one cohort (they sum to 1.000 age-wise) and are renormalised
#                to sum to exactly 1, as in the engines.
#   Mayaro       11/16 with arthralgia at the 3-month visit (Halsey et al. 2015), the
#                third component of Dirichlet(1, 4, 11); its exact marginal is Beta(11, 5).
#
# Summed PER DRAW, never by adding medians: medians are not additive, and adding the three
# components' UI endpoints would assume they move in lockstep, giving an interval ~1.6x too
# wide (the three chikungunya Betas are drawn independently, so variances add, not SDs).
# Means happen to be additive, and these Betas are near-symmetric, so summing the medians
# would be within ~0.15% here -- but the interval would not be.
# ------------------------------------------------------------
source("ca_common.R")
dpar <- load_daly_params()
set.seed(20260815); N_MC <- 2e5
rb <- function(p) rbeta(N_MC, p$a, p$b)

chronic_share <- function(p14, p90, p6, p12, p30) {
  a <- rb(p14); s <- rb(p90); ch <- rb(p6) + rb(p12) + rb(p30)
  tot <- a + s + ch                       # renormalise exactly as the engines do
  list(acute = a/tot, subacute = s/tot, chronic = ch/tot, raw_total = tot)
}
cs_y <- chronic_share(dpar$p14_y, dpar$p90_y, dpar$p6_y, dpar$p12_y, dpar$p30_y)
cs_o <- chronic_share(dpar$p14_o, dpar$p90_o, dpar$p6_o, dpar$p12_o, dpar$p30_o)
mv_chronic <- rbeta(N_MC, 11, 5)          # exact marginal of Dirichlet(1, 4, 11)

pc <- function(v, d = 1) { q <- quantile(v, c(.5, .025, .975))
  sprintf("%.*f%% (%.*f-%.*f)", d, 100*q[1], d, 100*q[2], d, 100*q[3]) }

chronic_tbl <- data.frame(
  pathogen  = c("Chikungunya", "Chikungunya", "Mayaro"),
  age_group = c("< 40 years", ">= 40 years", "All ages (not stratified)"),
  still_symptomatic_3m = c(pc(cs_y$chronic), pc(cs_o$chronic), pc(mv_chronic)),
  distribution = c("Sum of Beta(799.3, 3307.0), Beta(77.6, 836.7) and Beta(7.6, 1075.0), renormalised",
                   "Sum of Beta(5998.9, 21654.9), Beta(184.7, 1216.1) and Beta(15.8, 516.3), renormalised",
                   "Beta(11, 5), the chronic margin of Dirichlet(1, 4, 11)"),
  source = c("O'Driscoll et al. 2021 IJID (pooled cohorts)",
             "O'Driscoll et al. 2021 IJID (pooled cohorts)",
             "Halsey et al. 2015 (Peru, n = 16)"),
  stringsAsFactors = FALSE)

# Full three-phase partition, for readers who want to see that the rows sum to 1.
partition_tbl <- data.frame(
  pathogen  = c(rep("Chikungunya", 2), "Mayaro"),
  age_group = c("< 40 years", ">= 40 years", "All ages (not stratified)"),
  resolved_by_14d      = c(pc(cs_y$acute),    pc(cs_o$acute),    pc(rbeta(N_MC, 1, 15))),
  resolved_14d_to_3m   = c(pc(cs_y$subacute), pc(cs_o$subacute), pc(rbeta(N_MC, 4, 12))),
  still_symptomatic_3m = c(pc(cs_y$chronic),  pc(cs_o$chronic),  pc(mv_chronic)),
  stringsAsFactors = FALSE)

chronic_notes <- data.frame(item = c(
  "Quantity", "Chikungunya derivation", "Mayaro derivation", "Marginal, not conditional",
  "Why medians are not summed", "Interval widths are not comparable",
  "Mayaro is prevalence, not survival", "Mayaro acute bin", "No pooled value"),
  detail = c(
  "Proportion of symptomatic cases still reporting arthralgia at approximately three months from symptom onset. Median (95% UI) over 200,000 Monte Carlo draws.",
  "One minus the cumulative probability of recovery by 90 days, i.e. the sum of the 6-month, 12-month and 30-month recovery rows, renormalised so the three phases sum to 1. The five O'Driscoll rows sum to 1.000 (<40) and 1.000 (>=40) using the published Beta hyperparameters, confirming they partition one cohort.",
  "Proportion with arthralgia at the 3-month visit, 11 of 16 patients. Drawn as the chronic component of a Dirichlet(1, 4, 11) so the three phases stay correlated and sum to 1 in every draw; the exact marginal is Beta(11, 5).",
  "Despite the 'within X after Y period' wording of the source tables, these are marginal shares of recovery time, not conditional recovery probabilities. Reading row 2 as conditional on not having recovered in row 1 gives a different number (26.7% rather than 23.9% for Mayaro).",
  "Components are summed per draw, not by adding medians. Medians are not additive in general; adding the three components' UI endpoints would assume perfect positive correlation and yields an interval about 1.6 times too wide, since the chikungunya Betas are drawn independently.",
  "The chikungunya intervals are narrow because O'Driscoll's Beta hyperparameters encode very large pooled sample sizes; they describe precision in the source meta-analysis, not transportability to Caldas Novas. The Mayaro interval reflects n = 16. The widths differ by roughly 17-fold and should not be read as relative confidence.",
  "Halsey reports arthralgia prevalence falling to 3/16 at day 20 and rising to 11/16 at 3 months, a biphasic pattern the authors treat as real. Recovery is therefore not absorbing and the 3-month figure is a POINT PREVALENCE, not the probability of never having recovered.",
  "Halsey's first follow-up visit was day 20, so the 14-day boundary was never observed. The acute share is 'no arthralgia at the acute visit' (1 of 16), which is why it is far below the chikungunya values.",
  "No age-pooled chikungunya value is given: the source stratifies at 40 years and a pooled figure would depend on the age distribution of symptomatic cases, making it a model output rather than a parameter."),
  stringsAsFactors = FALSE)

write_xlsx(list(still_symptomatic_3m = chronic_tbl,
                full_partition       = partition_tbl,
                notes                = chronic_notes),
           "combined_chronic_share.xlsx")
cat("\nSaved combined_chronic_share.xlsx (proportion still symptomatic at 3 months).\n")
print(chronic_tbl[, 1:3], row.names = FALSE)
cat("\nFull partition:\n"); print(partition_tbl, row.names = FALSE)
cat(sprintf("\nO'Driscoll rows sum to %.4f (<40) and %.4f (>=40) before renormalisation.\n",
            median(cs_y$raw_total), median(cs_o$raw_total)))
