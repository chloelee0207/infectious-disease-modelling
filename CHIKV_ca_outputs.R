# ============================================================
# CHIKV_ca_outputs.R -- Caldas Novas CHIKV presentation & export layer.
# ------------------------------------------------------------
# Loads CHIKV_ca_vacc_results.rds (produced by CHIKV_ca_vacc.R, the standalone LHS-based
# engine) and produces the figures + the Excel workbook. Pure presentation.
#
# Run order:  source("ca_common.R"); source("CHIKV_ca_lhs.R")   # ensemble (slow, once)
#             source("CHIKV_ca_vacc.R")                          # engine -> CHIKV_ca_vacc_results.rds
#             source("CHIKV_ca_outputs.R")                       # <- this file
# ============================================================
library(dplyr); library(tidyr); library(ggplot2); library(writexl)
if (!exists("fmtq")) source("ca_common.R")

if (!file.exists("CHIKV_ca_vacc_results.rds"))
  stop("CHIKV_ca_vacc_results.rds not found -- run CHIKV_ca_vacc.R first.")
R <- readRDS("CHIKV_ca_vacc_results.rds")
bmat <- R$bmat; av <- R$av; wk_symp <- R$wk_symp; base_true <- R$base_true
rho_i <- R$rho_i; scen_names <- R$scen_names; vac_names <- R$vac_names
timings <- R$timings; outcomes <- R$outcomes; T_sim <- R$T_sim; T_data <- R$T_data
caldas_obs <- R$caldas_obs; observed_cases <- R$observed_cases

# ------------------------------------------------------------
# 1. Full-horizon calendar axis (data weeks 1..52 + extension 53..T_sim)
# ------------------------------------------------------------
grid <- data.frame(week_index = 1:T_sim)
grid$Year <- ifelse(grid$week_index <= T_data,
                    caldas_obs$Year[match(grid$week_index, caldas_obs$week_index)], 2026L)
grid$wk   <- ifelse(grid$week_index <= T_data,
                    caldas_obs$week[match(grid$week_index, caldas_obs$week_index)],
                    22L + (grid$week_index - T_data))
tick_specs <- list(c(2025,30), c(2025,40), c(2025,50), c(2026,10),
                   c(2026,20), c(2026,30), c(2026,40))
x_breaks <- sapply(tick_specs, function(s) { j <- which(grid$Year==s[1] & grid$wk==s[2]); if (length(j)) j[1] else NA })
x_labs   <- sapply(tick_specs, function(s) s[2])
keep <- !is.na(x_breaks); x_breaks <- x_breaks[keep]; x_labs <- x_labs[keep]
year_break <- mean(c(max(which(grid$Year==2025)), min(which(grid$Year==2026))))

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
scenario_totals <- data.frame(scenario = rownames(R$pt_burden),
                              round(as.data.frame(R$pt_burden[, outcomes, drop=FALSE]), 1),
                              row.names = NULL, check.names = FALSE)
notes <- data.frame(
  parameter = c("Data window", "Simulation horizon", "Observed reported cases",
                "Reporting rate rho", "Symptomatic fraction (prop_symp)",
                "Prior immunity model", "Vaccine efficacy (shared VE_inf/VE_block)",
                "Coverage of 18-59", "Weekly delivery", "Deployment delay",
                "Uncertainty", "Feasible LHS draws"),
  value = c(sprintf("2025-W24 -> 2026-W22 (%d wk), fit", T_data),
            sprintf("%d weeks (data + %d extension)", T_sim, R$EXTEND),
            format(sum(observed_cases), big.mark=","),
            "Beta(20,60), median 0.25 (point est.)",
            "Beta(35.84,32.56), median 0.524",
            "catalytic FOI (1 - exp(-FOI*age)), FOI ~ Lognormal",
            "Beta(98.9%, 96.7-99.8)",
            "Beta(30%, 95% 20-40)",
            "Beta(10%, 95% 9-11)",
            "median 2 wk (sampled 1-3), shifts scenario start",
            "5-input LHS (FOI,gamma,sigma,rho,prop_symp) re-fit per draw + vaccine LHS",
            as.character(nrow(base_true))),
  stringsAsFactors = FALSE)

sheets <- list(notes = notes, baseline_true_reported = R$base_tbl,
               vaccinated_true_reported = vtr, averted_MC_95UI = R$mc_tbl,
               scenario_totals = scenario_totals)
write_xlsx(sheets, "CHIKV_ca_vacc_outputs.xlsx")
cat("Wrote CHIKV_ca_vacc_outputs.xlsx (sheets:", paste(names(sheets), collapse=", "), ")\n")

# ------------------------------------------------------------
# 3. Figure (a): epidemic-curve ribbons, per timing (extended horizon)
# ------------------------------------------------------------
scen_cols <- c("No vaccination"="grey55", "Disease-blocking"="#4393c3",
               "Disease + infection blocking"="#d6604d")
band <- function(nm, lab) {
  q <- apply(wk_symp[[nm]], 2, quantile, c(.025,.5,.975), na.rm=TRUE)
  data.frame(week=1:T_sim, lo=q[1,], med=q[2,], hi=q[3,], scenario=lab)
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
  roll0 <- timings[[tn]] + 2; roll1 <- roll0 + 10   # median deploy delay + 10-wk rollout
  lab <- paste0("% symptomatic reduction:\nDisease-blocking: ", pct_of(disb),
                "\nDisease + infection blocking: ", pct_of(both))
  ggplot(pdf, aes(week, med, colour=scenario, fill=scenario)) +
    annotate("rect", xmin=roll0, xmax=roll1, ymin=-Inf, ymax=Inf, fill="#3a7d3a", alpha=.10) +
    geom_vline(xintercept=year_break, linetype="dashed", colour="grey60") +
    geom_vline(xintercept=T_data+0.5, linetype="dotted", colour="grey60") +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=.18, colour=NA) + geom_line(linewidth=.9) +
    annotate("text", x=Inf, y=Inf, label=lab, hjust=1.02, vjust=1.2, size=2.7, lineheight=.95) +
    scale_colour_manual(values=scen_cols, aesthetics=c("colour","fill")) +
    scale_x_continuous(breaks=x_breaks, labels=x_labs) +
    scale_y_continuous(labels=scales::comma) +
    labs(x="Week", y="Predicted symptomatic CHIKV cases", colour=NULL, fill=NULL,
         title=paste0("Caldas Novas CHIKV: symptomatic cases at 30% coverage - ", tn),
         caption=sprintf("Green = rollout window; dotted = end of observed data (wk %d); band = 95%% UI over LHS draws", T_data)) +
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

cat("Saved figures: CHIKV_ca_vacc_epicurve_{actual_rollout,start_of_2026,pre_outbreak}.png, CHIKV_ca_vacc_averted_mc.png\n")
