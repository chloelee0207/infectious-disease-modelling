# ============================================================
# combine_residual_burden.R -- merge the CHIKV and MAYV residual-burden charts
# into one three-column figure:
#
#   A Chikungunya : Disease-blocking | Disease + infection blocking
#   B Mayaro      : Disease-blocking
#
# Rows are the three outcomes; row strips sit on the far right of the whole figure.
#
# Reads the summarised frames written by the two outputs scripts, so it does no
# arithmetic of its own -- every value here is already per-draw propagated.
#
# The Mayaro deaths cell is blank. MAYV_ZERO_DEATHS = TRUE fixes the MAYV CFR at 0, so
# that ratio is undefined rather than zero. The facet level is kept (drop = FALSE) only
# to hold the row grid aligned with Chikungunya; the panel grob is then deleted from the
# assembled gtable, so no box, border or gridlines are drawn. The right-hand
# "Cumulative deaths" strip stays, because it labels that row for the Chikungunya block.
#
# Run order: CHIKV_ca_outputs.R and MAYV_ca_outputs.R -> this file
# ============================================================
library(ggplot2); library(patchwork)

need <- c("CHIKV_ca_residual_burden.rds", "MAYV_ca_residual_burden.rds")
miss <- need[!file.exists(need)]
if (length(miss))
  stop("missing ", paste(miss, collapse = " and "),
       " -- run CHIKV_ca_outputs.R / MAYV_ca_outputs.R first.")
chik <- readRDS(need[1]); mayv <- readRDS(need[2])

OUT_LV  <- c("Cumulative DALYs", "Cumulative deaths", "Healthcare cost")
FILL    <- c("No vaccination" = "grey60", "Vaccination" = "#4e79a7")
chik$outcome <- factor(as.character(chik$outcome), levels = OUT_LV)
mayv$outcome <- factor(as.character(mayv$outcome), levels = OUT_LV)  # deaths kept as an
mayv$arm     <- factor("Disease-blocking")                           # empty level

# common panel body; only the strip/axis furniture differs between the two blocks
panel <- function(d, title) {
  ggplot(d, aes(scenario, med, fill = scenario)) +
    geom_col(width = .6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, linewidth = .35) +
    scale_fill_manual(values = FILL, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 105), breaks = seq(0, 100, 25)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(11) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
          axis.text.x = element_text(size = 8.5),
          strip.text.x = element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

p_chik <- panel(chik, "A  Chikungunya") +
  facet_grid(outcome ~ arm) +
  theme(strip.text.y = element_blank())          # row labels live on the MAYV block

p_mayv <- panel(mayv, "B  Mayaro") +
  facet_grid(outcome ~ arm, drop = FALSE) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(face = "bold", size = 9))

# widths 2:1 so every panel column is the same physical width
# suppressWarnings: the empty deaths facet has no data for ggplot to take a range over.
p <- suppressWarnings(patchworkGrob(
  p_chik + p_mayv + plot_layout(ncol = 2, widths = c(2, 1), guides = "collect") &
    theme(legend.position = "bottom")))

# Blank the Mayaro deaths panel. The Mayaro block is the only sub-gtable with a single
# facet column (3 panels vs Chikungunya's 6); within it facet_grid names the rows
# panel-1-<row> in OUT_LV order, so deaths is panel-1-2.
blk <- which(vapply(p$grobs, function(gr) inherits(gr, "gtable") &&
                      sum(grepl("^panel-", gr$layout$name)) == 3L, logical(1)))
stopifnot(length(blk) == 1)
mg <- p$grobs[[blk]]
j  <- which(mg$layout$name == "panel-1-2")
stopifnot(length(j) == 1)
mg$grobs[[j]] <- grid::nullGrob()
p$grobs[[blk]] <- mg

# shared y-axis title, rotated, outside both blocks
fig <- wrap_elements(p) +
  labs(tag = "Cumulative burden (% of no vaccination)") +
  theme(plot.tag = element_text(size = 11, angle = 90),
        plot.tag.position = "left")

ggsave("combined_residual_burden.png", fig, width = 10, height = 7.4, dpi = 200)
ggsave("combined_residual_burden.pdf", fig, width = 10, height = 7.4)
cat("Saved combined_residual_burden.png / .pdf",
    "(A Chikungunya: 2 arms; B Mayaro: disease-blocking, deaths cell blank).\n")
