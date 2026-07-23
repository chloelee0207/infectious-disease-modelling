# ============================================================
# CHIKV_ca_outputs.R -- Caldas Novas CHIKV presentation and export layer.
#
# Reads CHIKV_ca_engine_results.rds and writes the burden workbook and figures.
# Pure presentation: no SEIR or Monte Carlo of its own, so every number here shares
# the engine's per-draw propagation.
#
# Horizon is the 52-week observed window, 2025-W24 -> 2026-W22.
#
# Run order:  CHIKV_ca_engine.R  ->  this file
# ============================================================
library(dplyr); library(tidyr); library(ggplot2); library(writexl)
if (!exists("fmtq")) source("ca_common.R")

if (!file.exists("CHIKV_ca_engine_results.rds"))
  stop("CHIKV_ca_engine_results.rds not found -- run CHIKV_ca_engine.R first.")
G <- readRDS("CHIKV_ca_engine_results.rds")
bmat <- G$per_draw; wk_symp <- G$wk_symp; wk_inf <- G$wk_inf; rho_i <- G$rho_i
scen_names <- G$scen_names; vac_names <- G$vac_names
timings <- G$timings; T_sim <- G$T_sim; T_data <- G$T_data
arm_names <- G$arm_names; nnv <- G$nnv; NNV_OUT <- G$NNV_OUT; cov_d <- G$cov_d
caldas_obs <- G$caldas_obs; observed_cases <- G$observed_cases
outcomes <- c("infections","symptomatic","hospitalisations","deaths")   # burden subset of OUTCOMES
stopifnot(T_sim == T_data)                        # horizon is the observed window
base_true <- bmat[["No vaccine (baseline)"]][, outcomes, drop = FALSE]

# averted (baseline - scenario), per draw, for the burden outcomes (computed locally
# so infections is included -- G$averted carries only the NNV outcomes).
av <- setNames(lapply(vac_names, function(nm)
  base_true - bmat[[nm]][, outcomes, drop = FALSE]), vac_names)

# baseline true-vs-reported and averted-MC tables (built here; the engine stores raw draws)
base_tbl <- do.call(rbind, lapply(outcomes, function(o) {
  d <- if (o == "deaths") 1 else 0
  data.frame(outcome = o, true = fmtq(base_true[,o], d), reported = fmtq(rho_i*base_true[,o], d))
}))
mc_tbl <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- av[[nm]]
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
             Infections = fmtq(m[,"infections"]), Symptomatic = fmtq(m[,"symptomatic"]),
             Hospitalisations = fmtq(m[,"hospitalisations"],1), Deaths = fmtq(m[,"deaths"],2),
             pct_symp = sprintf("%.1f%%", 100*median(m[,"symptomatic"]/base_true[,"symptomatic"], na.rm=TRUE)),
             row.names = NULL)
}))

# Outcomes averted per 100,000 doses (= 1e5 x averted / doses, per draw). Scale-free,
# so it is comparable across settings of different population size. Infections are
# averted only by the infection-blocking arm, so disease-blocking shows NA there.
per100k_outcomes <- c(infections="Infections", symptomatic="Symptomatic",
                      hospitalisations="Hospitalisations", deaths="Deaths", daly="DALYs")
mc_per100k <- do.call(rbind, lapply(vac_names, function(nm) {
  doses <- bmat[[nm]][,"doses"]
  base_o <- bmat[["No vaccine (baseline)"]]; scen_o <- bmat[[nm]]
  cells <- lapply(names(per100k_outcomes), function(o) {
    r <- 1e5 * (base_o[,o] - scen_o[,o]) / doses; r[!is.finite(r) | r < 0] <- NA
    fmtq(r, 0)
  })
  setNames(data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
                      cells, row.names = NULL, check.names = FALSE),
           c("timing","arm", unname(per100k_outcomes)))
}))

# ------------------------------------------------------------
# 1. Calendar axis for the 52 observed weeks
# ------------------------------------------------------------
grid <- data.frame(week_index = 1:T_sim)
grid$Year <- caldas_obs$Year[match(grid$week_index, caldas_obs$week_index)]
grid$wk   <- caldas_obs$week[match(grid$week_index, caldas_obs$week_index)]
tick_specs <- list(c(2025,30), c(2025,40), c(2025,50), c(2026,10), c(2026,20))
x_breaks <- sapply(tick_specs, function(s) { j <- which(grid$Year==s[1] & grid$wk==s[2]); if (length(j)) j[1] else NA })
x_labs   <- sapply(tick_specs, function(s) s[2])
keep <- !is.na(x_breaks); x_breaks <- x_breaks[keep]; x_labs <- x_labs[keep]
year_break <- mean(c(max(which(grid$Year==2025)), min(which(grid$Year==2026))))
grid$week_label <- sprintf("%d-W%02d", grid$Year, grid$wk)

# ------------------------------------------------------------
# 1b. Baseline (no-vaccine) weekly REPORTED cases.
# wk_symp holds per-draw weekly TRUE symptomatic; reported = rho_i * true (rho_i
# recycles down the columns, so row i is scaled by that draw's reporting rate).
# ------------------------------------------------------------
rep_mat <- wk_symp[["No vaccine (baseline)"]] * rho_i
rq <- apply(rep_mat, 2, quantile, c(.025, .5, .975), na.rm = TRUE)
weekly_reported <- data.frame(
  week_index  = 1:T_sim,
  week_label  = grid$week_label,
  observed    = observed_cases,
  pred_median = round(rq[2, ], 1),
  pred_lo     = round(rq[1, ], 1),
  pred_hi     = round(rq[3, ], 1))
write.csv(weekly_reported, "CHIKV_ca_vacc_weekly_reported.csv", row.names = FALSE)

# ------------------------------------------------------------
# 2. Excel workbook
# ------------------------------------------------------------
# per-scenario TOTALS, true & reported (median, 95% UI). rho already per-draw in the
# ensemble, so reported = rho_i * true carries the reporting-rate uncertainty jointly.
vtr <- do.call(rbind, lapply(scen_names, function(nm) {
  m <- bmat[[nm]]; base <- nm == "No vaccine (baseline)"
  data.frame(timing = if (base) "No vaccine" else sub(" \\|.*","",nm),
             arm    = if (base) "No vaccine" else sub(".*\\| ","",nm),
             outcome = outcomes,
             true     = sapply(outcomes, function(o) fmtq(m[,o],          if (o=="deaths") 1 else 0)),
             reported = sapply(outcomes, function(o) fmtq(rho_i*m[,o],    if (o=="deaths") 1 else 0)),
             row.names = NULL)
}))
scenario_totals <- data.frame(scenario = scen_names,
                              do.call(rbind, lapply(scen_names, function(nm)
                                round(apply(bmat[[nm]][, outcomes, drop=FALSE], 2, median, na.rm=TRUE), 1))),
                              row.names = NULL, check.names = FALSE)
notes <- data.frame(
  parameter = c("Data window", "Simulation horizon", "Observed reported cases",
                "Reporting rate rho", "Symptomatic fraction (prop_symp)",
                "Prior immunity model", "Vaccine efficacy (shared VE_inf/VE_block)",
                "Coverage of 18-59", "Weekly delivery", "Deployment delay",
                "Uncertainty", "Feasible LHS draws"),
  value = c(sprintf("2025-W24 -> 2026-W22 (%d wk)", T_data),
            sprintf("%d weeks; no projection beyond the observed data", T_sim),
            format(sum(observed_cases), big.mark=","),
            "Beta(20,60), median 0.25 (point est.)",
            "Beta(35.84,32.56), median 0.524",
            "truncated catalytic: 1 - exp(-FOI*min(age,12)), FOI ~ Lognormal",
            "Beta(98.9%, 96.7-99.8)",
            "Beta(30%, 95% 20-40)",
            "Beta(10%, 95% 9-11)",
            "median 2 wk (sampled 1-3), shifts scenario start",
            "5-input LHS (FOI,gamma,sigma,rho,prop_symp) re-fit per draw + vaccine LHS",
            as.character(nrow(base_true))),
  stringsAsFactors = FALSE)
notes <- rbind(notes, data.frame(parameter = "Scope of the window", value = paste(
  "Burden for EVERY scenario, including the no-vaccine baseline, is counted only within",
  "the 52-week observed window (2025-W24 -> 2026-W22). The model is not projected past",
  "the data, since that needs an assumed transmission rate beyond the last observation.",
  "Averted = baseline - scenario, both inside this window. For infection-blocking arms",
  "the within-window reduction combines cases prevented outright with cases the vaccine",
  "delays past 2026-W22; neither is counted after the window, by design."),
  stringsAsFactors = FALSE))

# doses actually delivered to the eligible 18-59, and dose wastage, per vaccine
# scenario. wastage = doses administered to already-immune / already-infected eligible
# people (cannot benefit) = 1 - on-target (reached susceptibles) / administered; it
# rises for later-timing rollouts, when more of the eligible are already infected.
pct3 <- function(x) { q <- quantile(x, c(.5,.025,.975), na.rm = TRUE)
  sprintf("%.1f%% (%.1f - %.1f%%)", 100*q[1], 100*q[2], 100*q[3]) }
doses_wastage <- do.call(rbind, lapply(vac_names, function(nm) {
  del <- G$doses_deliv[[nm]]; ont <- G$doses_ontarget[[nm]]
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
             `doses delivered` = fmtq(del, 0),
             `doses on-target (to susceptibles)` = fmtq(ont, 0),
             `wastage %` = pct3(1 - ont/del), check.names = FALSE, row.names = NULL)
}))

sheets <- list(notes = notes, baseline_true_reported = base_tbl,
               vaccinated_true_reported = vtr, averted_MC_95UI = mc_tbl,
               averted_per_100k_doses = mc_per100k, doses_wastage = doses_wastage,
               scenario_totals = scenario_totals, weekly_reported = weekly_reported,
               burden_audit = G$burden_audit, burden_audit_by_age = G$burden_audit_by_age)
write_xlsx(sheets, "CHIKV_ca_vacc_outputs.xlsx")
cat("Wrote CHIKV_ca_vacc_outputs.xlsx  (sheets:", paste(names(sheets), collapse=", "), ")\n")

# ------------------------------------------------------------
# 3. Figure (a): epidemic-curve ribbons, per timing
# ------------------------------------------------------------
scen_cols <- c("No vaccination"="grey55", "Disease-blocking"="#4393c3",
               "Disease + infection blocking"="#d6604d")
# TRUE symptomatic (solid, with ribbon) and REPORTED = rho_i * symptomatic (dotted).
band <- function(nm, lab) {
  qt <- apply(wk_symp[[nm]],          2, quantile, c(.025,.5,.975), na.rm=TRUE)
  qr <- apply(wk_symp[[nm]] * rho_i,  2, quantile, c(.025,.5,.975), na.rm=TRUE)
  rbind(
    data.frame(week=1:T_sim, lo=qt[1,], med=qt[2,], hi=qt[3,], scenario=lab, measure="True symptomatic"),
    data.frame(week=1:T_sim, lo=qr[1,], med=qr[2,], hi=qr[3,], scenario=lab, measure="Reported"))
}
pct_of <- function(nm) sprintf("%.1f%% (%.1f-%.1f%%)",
  100*median(av[[nm]][,"symptomatic"]/base_true[,"symptomatic"]),
  100*quantile(av[[nm]][,"symptomatic"]/base_true[,"symptomatic"], .025, na.rm=TRUE),
  100*quantile(av[[nm]][,"symptomatic"]/base_true[,"symptomatic"], .975, na.rm=TRUE))

make_epicurve <- function(tn) {
  disb <- paste0(tn," | Disease-blocking"); both <- paste0(tn," | Disease + infection blocking")
  pdf <- rbind(band("No vaccine (baseline)","No vaccination"),
               band(disb,"Disease-blocking"), band(both,"Disease + infection blocking"))
  pdf$scenario <- factor(pdf$scenario, levels=names(scen_cols))
  pdf$measure  <- factor(pdf$measure, levels=c("True symptomatic","Reported"))
  roll0 <- timings[[tn]] + 2; roll1 <- roll0 + 10   # median deploy delay + 10-wk rollout
  lab <- paste0("% symptomatic reduction:\nDisease-blocking: ", pct_of(disb),
                "\nDisease + infection blocking: ", pct_of(both))
  ggplot(pdf, aes(week, med, colour=scenario, fill=scenario)) +
    annotate("rect", xmin=roll0, xmax=roll1, ymin=-Inf, ymax=Inf, fill="#3a7d3a", alpha=.10) +
    geom_vline(xintercept=year_break, linetype="dashed", colour="grey60") +
    # geom_vline(xintercept=T_data+0.5, linetype="dotted", colour="grey60") +
    geom_ribbon(data=subset(pdf, measure=="True symptomatic"), aes(ymin=lo, ymax=hi), alpha=.18, colour=NA) +
    geom_line(aes(linetype=measure), linewidth=.9) +
    annotate("text", x=Inf, y=Inf, label=lab, hjust=1.02, vjust=1.2, size=3.5, lineheight=.95) +
    scale_colour_manual(values=scen_cols, aesthetics=c("colour","fill")) +
    scale_linetype_manual(values=c("True symptomatic"="solid","Reported"="dotted"), name=NULL) +
    scale_x_continuous(breaks=x_breaks, labels=x_labs) +
    scale_y_continuous(labels=scales::comma) +
    labs(x="Week", y="Predicted symptomatic cases", colour=NULL, fill=NULL,
         title=paste0("CHIKV symptomatic cases at 30% coverage")) +
         # caption=sprintf("Green = vaccination rollout window; band = 95%% UI", T_data)) +
    theme_bw(11) + theme(legend.position="bottom", plot.title=element_text(face="bold"),
                         panel.grid.minor=element_blank())
}
for (tn in names(timings)) {
  p <- make_epicurve(tn); print(p)
  ggsave(sprintf("CHIKV_ca_vacc_epicurve_%s.png", gsub("[^a-z0-9]+","_",tolower(tn))), p, width=8, height=5, dpi=120)
}

# ------------------------------------------------------------
# 4. Figure (b): averted burden by timing x arm (median + 95% UI)
# ------------------------------------------------------------
mc_long <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- av[[nm]]
  data.frame(timing=sub(" \\|.*","",nm), arm=sub(".*\\| ","",nm), outcome=outcomes,
             med=apply(m,2,median,na.rm=TRUE), lo=apply(m,2,quantile,.025,na.rm=TRUE),
             hi=apply(m,2,quantile,.975,na.rm=TRUE), row.names=NULL)
}))
mc_long$timing  <- factor(mc_long$timing, levels=names(timings))
mc_long$arm     <- factor(mc_long$arm, levels=c("Disease-blocking","Disease + infection blocking"))
mc_long$outcome <- factor(mc_long$outcome, levels=outcomes,
                          labels=c("Infections","Symptomatic","Hospitalisations","Deaths"))
p_bar <- ggplot(mc_long, aes(timing, med, fill=arm)) +
  geom_col(position=position_dodge(.8), width=.7) +
  geom_errorbar(aes(ymin=lo, ymax=hi), position=position_dodge(.8), width=.25, linewidth=.4) +
  facet_wrap(~outcome, scales="free_y", nrow=1) +
  scale_fill_manual(values=c("Disease-blocking"="#9ecae1","Disease + infection blocking"="#08519c"), name=NULL) +
  scale_y_continuous(labels=scales::comma) +
  labs(x=NULL, y="Averted (vs no vaccine), median + 95% UI",
       title="Caldas Novas CHIKV: burden averted by vaccination timing") +
  theme_bw(11) + theme(plot.title=element_text(face="bold", hjust=.5),
                       axis.text.x=element_text(angle=30, hjust=1, size=8),
                       legend.position="bottom", panel.grid.minor=element_blank())
print(p_bar); ggsave("CHIKV_ca_vacc_averted_mc.png", p_bar, width=10, height=4.2, dpi=120)

# ------------------------------------------------------------
# 5. Figure (c): baseline observed vs predicted reported cases.
# ------------------------------------------------------------
pred <- data.frame(week = 1:T_sim, lo = rq[1, ], med = rq[2, ], hi = rq[3, ])
obs  <- data.frame(week = 1:T_data, cases = observed_cases)
p_fit <- ggplot(pred, aes(week, med)) +
  geom_vline(xintercept = year_break, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#4393c3", alpha = .20) +
  geom_line(colour = "#2166ac", linewidth = .9) +
  geom_point(data = obs, aes(week, cases), inherit.aes = FALSE, size = 1.3, colour = "grey20") +
  scale_x_continuous(breaks = x_breaks, labels = x_labs) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week", y = "Reported CHIKV cases", colour = NULL,
       title = "Caldas Novas CHIKV: observed vs predicted reported cases (no vaccine)",
       caption = "Dots = observed; line = median, band = 95% UI over LHS draws; dashed = year boundary") +
  theme_bw(11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
print(p_fit); ggsave("CHIKV_ca_vacc_fit_observed.png", p_fit, width = 9, height = 5, dpi = 120)

# ------------------------------------------------------------
# 6. DALYs: workbook + composition and averted figures.
# YLD (by phase), YLL and DALY are per-draw engine outcomes, consistent draw-for-draw
# with the burden totals. DALY = YLD + YLL, undiscounted.
# ------------------------------------------------------------
base_daly <- bmat[["No vaccine (baseline)"]]
daly_by_scenario <- do.call(rbind, lapply(scen_names, function(nm) {
  m <- bmat[[nm]]; b <- nm == "No vaccine (baseline)"
  data.frame(timing = if (b) "No vaccine" else sub(" \\|.*","",nm),
             arm    = if (b) "No vaccine" else sub(".*\\| ","",nm),
             YLD = fmtq(m[,"yld"]), YLL = fmtq(m[,"yll"]), DALY = fmtq(m[,"daly"]),
             YLD_acute = fmtq(m[,"yld_acute"]), YLD_subacute = fmtq(m[,"yld_subacute"]),
             YLD_chronic = fmtq(m[,"yld_chronic"]), row.names = NULL, check.names = FALSE)
}))
daly_averted <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- G$averted[[nm]]
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
             DALY_averted = fmtq(m[,"daly"]),
             pct_DALY = sprintf("%.1f%%", 100*median(m[,"daly"]/base_daly[,"daly"], na.rm=TRUE)),
             row.names = NULL, check.names = FALSE)
}))
write_xlsx(list(daly_by_scenario = daly_by_scenario, daly_averted = daly_averted),
           "CHIKV_ca_daly_outputs.xlsx")

# (a) baseline DALY composition
comp <- data.frame(component = factor(c("YLD (acute)","YLD (sub-acute)","YLD (chronic)","YLL"),
                     levels = c("YLD (acute)","YLD (sub-acute)","YLD (chronic)","YLL")),
                   value = c(median(base_daly[,"yld_acute"]), median(base_daly[,"yld_subacute"]),
                             median(base_daly[,"yld_chronic"]), median(base_daly[,"yll"])))
p_comp <- ggplot(comp, aes("Baseline", value, fill = component)) +
  geom_col(width = .55) +
  scale_fill_manual(values = c("#c6dbef","#6baed6","#2171b5","#d6604d"), name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "DALYs (median)", title = "Caldas Novas CHIKV: baseline DALY composition (no vaccine)") +
  theme_bw(11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
print(p_comp); ggsave("CHIKV_ca_daly_composition.png", p_comp, width = 5.5, height = 5, dpi = 120)

# (b) DALYs averted by timing x arm
dav_long <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- G$averted[[nm]]
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm),
             med = median(m[,"daly"], na.rm=TRUE),
             lo = quantile(m[,"daly"],.025,na.rm=TRUE), hi = quantile(m[,"daly"],.975,na.rm=TRUE),
             row.names = NULL)
}))
dav_long$timing <- factor(dav_long$timing, levels = names(timings))
dav_long$arm    <- factor(dav_long$arm, levels = arm_names)
p_dav <- ggplot(dav_long, aes(timing, med, fill = arm)) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.8), width = .25, linewidth = .4) +
  scale_fill_manual(values = c("Disease-blocking"="#9ecae1","Disease + infection blocking"="#08519c"), name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "DALYs averted (vs no vaccine), median + 95% UI",
       title = "Caldas Novas CHIKV: DALYs averted by vaccination timing") +
  theme_bw(11) + theme(plot.title = element_text(face="bold", hjust=.5),
                       axis.text.x = element_text(angle = 20, hjust = 1),
                       legend.position = "bottom", panel.grid.minor = element_blank())
print(p_dav); ggsave("CHIKV_ca_daly_averted.png", p_dav, width = 8, height = 4.5, dpi = 120)

# ------------------------------------------------------------
# 7. Number Needed to Vaccinate: workbook + figure (pre-outbreak rollout).
# NNV = doses / burden averted, per draw, so it shares the burden/DALY propagation.
# NNV < 1 means one dose averts more than one case via herd protection.
# ------------------------------------------------------------
out_labs <- c(infections = "Infection", symptomatic = "Symptomatic case",
              hospitalisations = "Hospitalisation", deaths = "Death", daly = "DALY")
nnv_tbl <- do.call(rbind, lapply(vac_names, function(nm) {
  m <- nnv[[nm]]
  vals <- setNames(as.list(vapply(NNV_OUT, function(o) fmtq(m[, o], 1), character(1))),
                   unname(out_labs[NNV_OUT]))
  data.frame(timing = sub(" \\|.*","",nm), arm = sub(".*\\| ","",nm), vals,
             row.names = NULL, check.names = FALSE)
}))
write_xlsx(list(nnv = nnv_tbl), "CHIKV_ca_nnv_outputs.xlsx")

focus_names <- grep("pre-outbreak", vac_names, value = TRUE)
nnv_plt <- do.call(rbind, lapply(focus_names, function(nm) {
  m <- nnv[[nm]]
  do.call(rbind, lapply(NNV_OUT, function(o) {
    q <- quantile(m[, o], c(.5, .025, .975), na.rm = TRUE)
    data.frame(arm = sub(".*\\| ","",nm), outcome = o, med = q[1], lo = q[2], hi = q[3], row.names = NULL)
  }))
}))
nnv_plt$outcome <- factor(out_labs[nnv_plt$outcome], levels = unname(out_labs))
nnv_plt$arm     <- factor(nnv_plt$arm, levels = arm_names)
p_nnv <- ggplot(nnv_plt, aes(arm, med, fill = arm)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45", linewidth = .3) +
  geom_col(width = .65, colour = "grey30", linewidth = .2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .2, linewidth = .35) +
  facet_wrap(~ outcome, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Disease-blocking"="#4393c3","Disease + infection blocking"="#d6604d"), name = NULL) +
  scale_y_log10(breaks = scales::breaks_log(n = 6), labels = scales::label_number(drop0trailing = TRUE)) +
  labs(x = NULL, y = "NNV to avert one outcome (pre-outbreak rollout)",
       title = "Number Needed to Vaccinate (NNV)") +
  theme_bw(11) + theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.grid.minor = element_blank())
print(p_nnv); ggsave("CHIKV_ca_nnv.png", p_nnv, width = 11, height = 4.6, dpi = 120)

cat(sprintf("Saved figures: epicurve_{%s}, averted_mc, fit_observed, daly_composition, daly_averted, nnv\n",
            paste(gsub("[^a-z0-9]+","_",tolower(names(timings))), collapse=", ")))
cat("Wrote CHIKV_ca_vacc_outputs.xlsx, CHIKV_ca_daly_outputs.xlsx, CHIKV_ca_nnv_outputs.xlsx,",
    "CHIKV_ca_vacc_weekly_reported.csv\n")
