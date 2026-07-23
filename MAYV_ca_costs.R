# ============================================================
# MAYV_ca_costs.R -- Caldas Novas MAYV direct MEDICAL cost layer.
# ------------------------------------------------------------
# The MAYV twin of CHIKV_ca_costs.R. Post-hoc multiply onto
# MAYV_ca_engine_results.rds. No SEIR re-run: the engine already stores per-draw
# case counts, so costs are applied draw-by-draw and the cost-parameter uncertainty
# is propagated jointly with the epidemic uncertainty.
#
# *** BORROWED COSTS -- READ THIS ***
# There is NO published Mayaro cost-of-illness study. As with the severity/DALY
# parameters, this layer BORROWS the CHIKV unit costs and treatment regimens:
#   Goncalves et al. 2024 (Rev Bras Epidemiol 27:e240026), Rio de Janeiro 2019, 2019 BRL,
#   regimens per the Brazilian MoH / SES-RJ chikungunya flowchart (10/01/19).
# MAYV and CHIKV are both alphaviruses producing an acute febrile arthralgia that
# can persist, so the CHIKV care pathway is the closest available proxy -- but these
# are "CHIKV-equivalent" costs, not measured MAYV costs. MAYV is generally milder,
# so treat them as an UPPER bound. The goncalves_check sheet confirms the per-case
# formulas reproduce the published CHIKV phase totals to <0.002%, i.e. it validates
# the borrowed FORMULAS (not MAYV itself).
#
# NOTE ON DENOMINATORS (same convention as the CHIKV layer):
#   Outpatient costs are applied to NON-HOSPITALISED modelled symptomatic cases only,
#   and inpatient costs to hospitalised cases -- mutually exclusive, no double count.
#   The denominator is MODELLED symptomatic (true burden); a reported-case equivalent
#   (x rho) is given in case_counts.
#
# Phase counts are NESTED (cumulative) for costs, unlike the DALY layer which is a
# disjoint exit funnel. A patient who becomes chronic also consumed acute and
# sub-acute care, so they pay at every phase they passed through.
#
# CONDITIONING: the MAYV epidemic only takes off in a minority of draws, so (exactly
# as the MAYV engine reports burden) every summary here is CONDITIONAL ON AN OUTBREAK
# (G$outbreak). A fixed-R0 representative outbreak is costed separately.
#
# DEATHS: MAYV has no confirmed attributable death, so the engine sets CFR = 0. This
# layer is DIRECT MEDICAL cost only and never costed deaths anyway, so nothing changes.
#
# Run after: MAYV_ca_engine.R
# ============================================================
library(dplyr); library(writexl); library(readxl)

set.seed(4042)
G  <- readRDS("MAYV_ca_engine_results.rds")
ND <- G$N_DRAWS
scen_names <- G$scen_names
ok <- G$outbreak                      # integer indices of draws that took off
stopifnot(length(ok) > 0)

# ------------------------------------------------------------
# 1. Unit-cost parameters: fit Gamma / Beta to median + 95% UI
# ------------------------------------------------------------
CO <- as.data.frame(read_excel("costs.xlsx", sheet = "costs"))
names(CO)[1:5] <- c("parameter", "median", "lo", "hi", "dist")
CO <- CO[!is.na(CO$parameter), ]
P <- function(nm) {
  r <- CO[CO$parameter == nm, ]
  if (nrow(r) != 1) stop("costs.xlsx: parameter not found / not unique: ", nm)
  r
}
gam_fit  <- function(m, lo, hi) { s <- (hi-lo)/(2*1.96); c(shape = (m/s)^2, rate = m/s^2) }
beta_fit <- function(m, lo, hi) { v <- ((hi-lo)/(2*1.96))^2; k <- m*(1-m)/v - 1; c(a = m*k, b = (1-m)*k) }

lhs_col <- function(n) (sample.int(n) - runif(n)) / n
PARS <- c("Medical appointment cost (2019 R$)",
          "Minor/moderate pain as share of acute phase cases",
          "Dipyrone - cost per case (2019 R$)", "Acetaminophen - cost per case (2019 R$)",
          "Tramadol - cost per case (2019 R$)", "Codeine - cost per case (2019 R$)",
          "Oxycodone - cost per case (2019 R$)",
          "Arthritis as share of sub-acute phase cases",
          "Prednisone - cost per case (2019 R$)", "Ibuprofen - cost per case (2019 R$)",
          "Amitriptyline - cost per case (2019 R$)", "Gabapentin - cost per case (2019 R$)",
          "Mild illness as share of chronic phase cases",
          "Hydroxychloroquine - cost per case (2019 R$)",
          "Methotrexate - cost per case (2019 R$)", "Folic acid - cost per case (2019 R$)",
          "Public share of inpatient admissions",
          "Average inpatient stay cost (2019 R$) - public sector",
          "Average inpatient stay cost (2019 R$) - private sector")
U <- sapply(seq_along(PARS), function(j) lhs_col(ND))
draws <- list(); fits <- list()
for (j in seq_along(PARS)) {
  r <- P(PARS[j])
  if (r$dist == "Beta") { ab <- beta_fit(r$median, r$lo, r$hi); fits[[PARS[j]]] <- ab
                          draws[[PARS[j]]] <- qbeta(U[, j], ab["a"], ab["b"]) }
  else                  { sr <- gam_fit(r$median, r$lo, r$hi);  fits[[PARS[j]]] <- sr
                          draws[[PARS[j]]] <- qgamma(U[, j], shape = sr["shape"], rate = sr["rate"]) }
}
d <- function(nm) draws[[nm]]
appt <- d(PARS[1]); mm_sh <- d(PARS[2]); dip <- d(PARS[3]); ace <- d(PARS[4])
tra  <- d(PARS[5]); cod   <- d(PARS[6]); oxy <- d(PARS[7]); art_sh <- d(PARS[8])
pred <- d(PARS[9]); ibu   <- d(PARS[10]); ami <- d(PARS[11]); gab <- d(PARS[12])
mild_sh <- d(PARS[13]); hcq <- d(PARS[14]); mtx <- d(PARS[15]); fol <- d(PARS[16])
pub_sh  <- d(PARS[17]); pub_c <- d(PARS[18]); priv_c <- d(PARS[19])

# ------------------------------------------------------------
# 2. Cost per case, per phase (vectors of length ND) -- CHIKV regimens (borrowed)
#      acute mild/moderate  : dipyrone AND paracetamol alternating -> SUM
#      acute severe         : dipyrone OR paracetamol (mean) + one opioid (mean of 3)
#      sub-acute arthritis  : prednisone
#      sub-acute no-arthritis: mean of [nerve-pain arm mean(ami,gab)] and [ibuprofen]
#      chronic mild         : hydroxychloroquine
#      chronic moderate/int : methotrexate + folic acid
# ------------------------------------------------------------
pc_acute <- 2*appt + mm_sh*(dip + ace) + (1 - mm_sh)*((dip + ace)/2 + (tra + cod + oxy)/3)
pc_sub   <- 3*appt + art_sh*pred + (1 - art_sh)*(((ami + gab)/2 + ibu)/2)
pc_chr   <- 3*appt + mild_sh*hcq + (1 - mild_sh)*(mtx + fol)
pc_hosp  <- pub_sh*pub_c + (1 - pub_sh)*priv_c

# ------------------------------------------------------------
# 3. Per-draw costs. Helper works on any per-draw matrix (sampled-R0 or fixed-R0).
#    NESTED phase entry shares (the engine stores DISJOINT exit groups).
# ------------------------------------------------------------
cost_of <- function(pd) {
  symp <- pd[, "symptomatic"]; hosp <- pd[, "hospitalisations"]
  nonh <- symp - hosp
  safe  <- function(x, y) ifelse(y > 0, x / y, 0)          # near-zero outbreaks
  p_sub <- safe(pd[, "n_subacute"] + pd[, "n_chronic"], symp)  # still ill past 14d
  p_chr <- safe(pd[, "n_chronic"], symp)                       # still ill past 3 months
  n_sub <- nonh * p_sub; n_chr <- nonh * p_chr
  ch <- hosp*pc_hosp; ca <- nonh*pc_acute; cs <- n_sub*pc_sub; cc <- n_chr*pc_chr
  list(cost = cbind(hosp_inpatient = ch, out_acute = ca, out_subacute = cs,
                    out_chronic = cc, total_direct_medical = ch + ca + cs + cc),
       count = cbind(symptomatic = symp, reported = pd[, "reported"], hospitalised = hosp,
                     non_hospitalised = nonh, entered_subacute = n_sub,
                     entered_chronic = n_chr))
}
COMP <- c("hosp_inpatient", "out_acute", "out_subacute", "out_chronic", "total_direct_medical")
cost_pd  <- setNames(lapply(scen_names, function(s) cost_of(G$per_draw[[s]])$cost),  scen_names)
count_pd <- setNames(lapply(scen_names, function(s) cost_of(G$per_draw[[s]])$count), scen_names)

# ------------------------------------------------------------
# 4. Summaries -- CONDITIONAL ON AN OUTBREAK (rows `ok`)
# ------------------------------------------------------------
q3  <- function(x) c(median(x), quantile(x, .025), quantile(x, .975))
fmt <- function(x, dp = 0) sprintf("%s (%s - %s)",
        formatC(round(x[1], dp), big.mark = ",", format = "f", digits = dp),
        formatC(round(x[2], dp), big.mark = ",", format = "f", digits = dp),
        formatC(round(x[3], dp), big.mark = ",", format = "f", digits = dp))

sh_costs <- do.call(rbind, lapply(scen_names, function(s) {
  m <- cost_pd[[s]][ok, , drop = FALSE]
  data.frame(scenario = s,
             hosp_inpatient       = fmt(q3(m[, "hosp_inpatient"])),
             out_acute            = fmt(q3(m[, "out_acute"])),
             out_subacute         = fmt(q3(m[, "out_subacute"])),
             out_chronic          = fmt(q3(m[, "out_chronic"])),
             TOTAL_direct_medical = fmt(q3(m[, "total_direct_medical"])),
             stringsAsFactors = FALSE) }))

sh_counts <- do.call(rbind, lapply(scen_names, function(s) {
  m <- count_pd[[s]][ok, , drop = FALSE]
  data.frame(scenario = s,
             symptomatic      = fmt(q3(m[, "symptomatic"])),
             reported         = fmt(q3(m[, "reported"])),
             hospitalised      = fmt(q3(m[, "hospitalised"]), 1),
             non_hospitalised = fmt(q3(m[, "non_hospitalised"])),
             entered_subacute = fmt(q3(m[, "entered_subacute"])),
             entered_chronic  = fmt(q3(m[, "entered_chronic"])),
             stringsAsFactors = FALSE) }))

# averted vs baseline, paired per draw (so the UI keeps the correlation)
base <- cost_pd[["No vaccine (baseline)"]]
vac  <- setdiff(scen_names, "No vaccine (baseline)")
sh_averted <- do.call(rbind, lapply(vac, function(s) {
  a  <- (base - cost_pd[[s]])[ok, , drop = FALSE]
  bs <- base[ok, , drop = FALSE]
  dz <- G$per_draw[[s]][ok, "doses"]
  data.frame(scenario = s,
             hosp_inpatient       = fmt(q3(a[, "hosp_inpatient"])),
             out_acute            = fmt(q3(a[, "out_acute"])),
             out_subacute         = fmt(q3(a[, "out_subacute"])),
             out_chronic          = fmt(q3(a[, "out_chronic"])),
             TOTAL_cost_averted   = fmt(q3(a[, "total_direct_medical"])),
             pct_of_baseline      = fmt(q3(100*a[, "total_direct_medical"] /
                                            bs[, "total_direct_medical"]), 1),
             cost_averted_per_dose= fmt(q3(a[, "total_direct_medical"] / dz), 2),
             stringsAsFactors = FALSE) }))

# ------------------------------------------------------------
# 5. Unit-cost audit: input vs fitted distribution vs realised draws
# ------------------------------------------------------------
sh_units <- do.call(rbind, lapply(PARS, function(nm) {
  r <- P(nm); f <- fits[[nm]]; x <- draws[[nm]]
  data.frame(parameter = nm, input_median = r$median, input_lo = r$lo, input_hi = r$hi,
             distribution = r$dist,
             fitted_par1 = round(unname(f[1]), 4), fitted_par2 = round(unname(f[2]), 4),
             drawn_median_95UI = fmt(q3(x), 3), stringsAsFactors = FALSE) }))

sh_percase <- data.frame(
  phase = c("Acute (per non-hospitalised case)", "Sub-acute (per case entering phase)",
            "Chronic (per case entering phase)", "Inpatient (per admission)"),
  formula = c("2*appt + mm*(dipyrone+acetaminophen) + (1-mm)*[(dipyrone+acetaminophen)/2 + mean(tramadol,codeine,oxycodone)]",
              "3*appt + arth*prednisone + (1-arth)*mean( mean(amitriptyline,gabapentin) , ibuprofen )",
              "3*appt + mild*hydroxychloroquine + (1-mild)*(methotrexate + folic acid)",
              "public_share*public_stay_cost + (1-public_share)*private_stay_cost"),
  cost_per_case_BRL2019 = c(fmt(q3(pc_acute), 2), fmt(q3(pc_sub), 2),
                            fmt(q3(pc_chr), 2), fmt(q3(pc_hosp), 2)),
  stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 6. Deterministic point-estimate chain (hand-checkable, baseline | outbreak)
# ------------------------------------------------------------
M <- function(nm) P(nm)$median
m_appt<-M(PARS[1]); m_mm<-M(PARS[2]); m_dip<-M(PARS[3]); m_ace<-M(PARS[4])
m_tra<-M(PARS[5]); m_cod<-M(PARS[6]); m_oxy<-M(PARS[7]); m_art<-M(PARS[8])
m_pred<-M(PARS[9]); m_ibu<-M(PARS[10]); m_ami<-M(PARS[11]); m_gab<-M(PARS[12])
m_mild<-M(PARS[13]); m_hcq<-M(PARS[14]); m_mtx<-M(PARS[15]); m_fol<-M(PARS[16])
m_pub<-M(PARS[17]); m_pubc<-M(PARS[18]); m_privc<-M(PARS[19])
mpc_a <- 2*m_appt + m_mm*(m_dip+m_ace) + (1-m_mm)*((m_dip+m_ace)/2 + (m_tra+m_cod+m_oxy)/3)
mpc_s <- 3*m_appt + m_art*m_pred + (1-m_art)*(((m_ami+m_gab)/2 + m_ibu)/2)
mpc_c <- 3*m_appt + m_mild*m_hcq + (1-m_mild)*(m_mtx+m_fol)
mpc_h <- m_pub*m_pubc + (1-m_pub)*m_privc
b <- G$per_draw[["No vaccine (baseline)"]][ok, , drop = FALSE]
B_symp <- median(b[, "symptomatic"]); B_hosp <- median(b[, "hospitalisations"])
B_nonh <- B_symp - B_hosp
B_psub <- median((b[, "n_subacute"] + b[, "n_chronic"]) / b[, "symptomatic"])
B_pchr <- median(b[, "n_chronic"] / b[, "symptomatic"])
sh_audit <- data.frame(
  step = 1:9,
  quantity = c("Symptomatic cases (median | outbreak)", "Hospitalised", "Non-hospitalised",
               "Share entering sub-acute", "Share entering chronic",
               "Inpatient cost", "Outpatient acute cost", "Outpatient sub-acute cost",
               "Outpatient chronic cost"),
  formula = c("from engine (conditional on outbreak)", "from engine", "symptomatic - hospitalised",
              "(n_subacute + n_chronic) / symptomatic", "n_chronic / symptomatic",
              sprintf("%.1f hospitalised x %.2f per admission", B_hosp, mpc_h),
              sprintf("%.0f non-hosp x %.2f per case", B_nonh, mpc_a),
              sprintf("%.0f x %.4f x %.2f per case", B_nonh, B_psub, mpc_s),
              sprintf("%.0f x %.4f x %.2f per case", B_nonh, B_pchr, mpc_c)),
  value = round(c(B_symp, B_hosp, B_nonh, B_psub, B_pchr, B_hosp*mpc_h, B_nonh*mpc_a,
                  B_nonh*B_psub*mpc_s, B_nonh*B_pchr*mpc_c), 4),
  stringsAsFactors = FALSE)
sh_audit <- rbind(sh_audit, data.frame(step = 10, quantity = "TOTAL direct medical (point est.)",
  formula = "steps 6+7+8+9", value = round(B_hosp*mpc_h + B_nonh*mpc_a +
    B_nonh*B_psub*mpc_s + B_nonh*B_pchr*mpc_c, 2), stringsAsFactors = FALSE))

# ------------------------------------------------------------
# 7. Fixed-R0 representative outbreak (R0_FIX): costs without the take-off lottery
# ------------------------------------------------------------
fx_base <- cost_of(G$fixed_base_pd)$cost
fx_vac  <- cost_of(G$fixed_vac_pd)$cost
fx_av   <- fx_base - fx_vac
sh_fixed <- rbind(
  data.frame(scenario = sprintf("No vaccine (fixed R0 = %.1f)", G$R0_FIX),
             hosp_inpatient = fmt(q3(fx_base[, "hosp_inpatient"])),
             out_acute      = fmt(q3(fx_base[, "out_acute"])),
             out_subacute   = fmt(q3(fx_base[, "out_subacute"])),
             out_chronic    = fmt(q3(fx_base[, "out_chronic"])),
             TOTAL_direct_medical = fmt(q3(fx_base[, "total_direct_medical"])),
             stringsAsFactors = FALSE),
  data.frame(scenario = sprintf("Disease-blocking vaccine (fixed R0 = %.1f)", G$R0_FIX),
             hosp_inpatient = fmt(q3(fx_vac[, "hosp_inpatient"])),
             out_acute      = fmt(q3(fx_vac[, "out_acute"])),
             out_subacute   = fmt(q3(fx_vac[, "out_subacute"])),
             out_chronic    = fmt(q3(fx_vac[, "out_chronic"])),
             TOTAL_direct_medical = fmt(q3(fx_vac[, "total_direct_medical"])),
             stringsAsFactors = FALSE),
  data.frame(scenario = "COST AVERTED (fixed R0)",
             hosp_inpatient = fmt(q3(fx_av[, "hosp_inpatient"])),
             out_acute      = fmt(q3(fx_av[, "out_acute"])),
             out_subacute   = fmt(q3(fx_av[, "out_subacute"])),
             out_chronic    = fmt(q3(fx_av[, "out_chronic"])),
             TOTAL_direct_medical = fmt(q3(fx_av[, "total_direct_medical"])),
             stringsAsFactors = FALSE))

# ------------------------------------------------------------
# 8. Validation: the BORROWED formulas reproduce Goncalves' published Rio 2019 totals
# ------------------------------------------------------------
N0 <- 38830; N1g <- 0.537*N0; N2g <- 0.52*N0
g_ac <- 2*23.17*N0 + 0.47*N0*(4.33+1.13) + 0.53*N0*((4.33+1.13)/2 + (4.24+20.95+89.33)/3)
g_sb <- 3*23.17*N1g + 0.755*N1g*16.28 + 0.245*N1g*(((1.06+9.16)/2 + 16.10)/2)
g_cr <- 3*23.17*N2g + 0.404*N2g*244.78 + 0.596*N2g*(109.86+0.87)
g_hp <- 256*347.37
sh_check <- data.frame(
  phase = c("Acute", "Post-acute", "Chronic", "Hospitalisations"),
  published_BRL = c(2740818, 1759900, 4732986, 88927),
  reproduced_BRL = round(c(g_ac, g_sb, g_cr, g_hp), 0),
  diff = round(c(g_ac-2740818, g_sb-1759900, g_cr-4732986, g_hp-88927), 0),
  pct_diff = sprintf("%.4f%%", 100*c(g_ac/2740818-1, g_sb/1759900-1, g_cr/4732986-1, g_hp/88927-1)),
  stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 9. Notes + write
# ------------------------------------------------------------
sh_notes <- data.frame(item = c(
 "Disease", "Currency", "Unit cost source", "BORROWED-COST CAVEAT", "Treatment regimen source",
 "Epidemic source", "Evaluation window", "Conditioning", "Draws",
 "Denominator (outpatient)", "Denominator (inpatient)", "Phase counts",
 "Deaths", "Hospitalisation cost", "Unused input",
 "Why medians do not add up", "Averted costs", "Fixed-R0 sheet", "Scope"),
 detail = c(
 "Mayaro virus (MAYV), Caldas Novas, hypothetical outbreak.",
 "2019 Brazilian reais (BRL). No inflation or PPP adjustment applied.",
 "Goncalves et al. 2024, Rev Bras Epidemiol 27:e240026 (Rio de Janeiro, 2019) -- a CHIKUNGUNYA study.",
 "There is NO published MAYV cost-of-illness study. CHIKV unit costs and regimens are BORROWED as the closest alphavirus proxy (same acute-arthralgia care pathway). MAYV is generally milder, so these are a CHIKV-equivalent UPPER bound, not measured MAYV costs.",
 "Brazilian MoH / SES-RJ chikungunya flowchart, 10/01/19.",
 "MAYV_ca_engine_results.rds -- unified Monte Carlo off the MAYV LHS ensemble.",
 sprintf("52 weeks, index %d-%d (2025-W24 -> 2026-W22; hybrid CHIKV-beta + dry-season envelope).", min(G$EVAL_WIN), max(G$EVAL_WIN)),
 sprintf("CONDITIONAL ON AN OUTBREAK: %d of %d draws took off (P = %.1f%%, attack > %.1f%% of susceptibles). Unconditional means would be dominated by non-take-off draws.", length(ok), ND, 100*G$p_outbreak, G$OUTBREAK_ATTACK_THRESH),
 sprintf("%d. Cost parameters are drawn by Latin hypercube and paired row-wise with the epidemic draws, so cost and epidemic uncertainty propagate jointly.", ND),
 "NON-HOSPITALISED modelled symptomatic cases (hospitalised excluded to avoid double counting).",
 "Hospitalised cases from the engine (symptomatic x hospitalisation rate).",
 "NESTED/cumulative: a chronic patient also consumed acute and sub-acute care and is charged at each phase. (The DALY layer uses a disjoint exit funnel -- different question, different structure.)",
 "MAYV has no confirmed attributable death, so the engine sets CFR = 0. This layer is direct MEDICAL cost only and never costed deaths, so the zero has no effect here.",
 "Per ADMISSION, not per day. Blended as public_share x public + (1-public_share) x private.",
 "'Average hospitalisation LOS, days' is NOT used -- the stay cost is already per admission. Retained in costs.xlsx for reference only.",
 "Every figure is the median (and 2.5-97.5th percentile) of the per-draw distribution. The median of a sum is not the sum of medians, so component medians will not add exactly to the TOTAL median. The per-draw totals ARE internally consistent.",
 "Computed per draw as (baseline - scenario) before summarising, so baseline and scenario stay paired and the UI reflects the correlated uncertainty.",
 sprintf("costs_fixedR0 removes the take-off lottery: every draw is a genuine outbreak at R0 = %.1f, so it answers 'what does an outbreak of this size cost?' rather than 'what does the average season cost?'.", G$R0_FIX),
 "DIRECT MEDICAL costs only. Excludes indirect/productivity losses, which were 97% of total costs in Goncalves."),
 stringsAsFactors = FALSE)

write_xlsx(list(notes = sh_notes, unit_costs = sh_units, cost_per_case = sh_percase,
                case_counts = sh_counts, costs_by_scenario = sh_costs,
                cost_averted = sh_averted, costs_fixedR0 = sh_fixed,
                audit_point_estimate = sh_audit, goncalves_check = sh_check),
           "MAYV_ca_costs.xlsx")
saveRDS(list(cost_pd = cost_pd, count_pd = count_pd, outbreak = ok,
             pc = list(acute = pc_acute, sub = pc_sub, chr = pc_chr, hosp = pc_hosp)),
        "MAYV_ca_costs.rds")

cat("Wrote MAYV_ca_costs.xlsx (9 sheets) + MAYV_ca_costs.rds\n\n")
cat("Borrowed-formula replication check (Goncalves CHIKV Rio 2019):\n"); print(sh_check)
cat(sprintf("\nBaseline MAYV direct medical cost | outbreak (BRL 2019, median [95%% UI]); %d/%d draws:\n",
            length(ok), ND))
bm <- cost_pd[["No vaccine (baseline)"]][ok, , drop = FALSE]
for (k in COMP) cat(sprintf("  %-22s %s\n", k, fmt(q3(bm[, k]))))
cat("\nCost averted by the disease-blocking vaccine | outbreak:\n")
am <- (base - cost_pd[[vac[1]]])[ok, , drop = FALSE]
for (k in COMP) cat(sprintf("  %-22s %s\n", k, fmt(q3(am[, k]))))
