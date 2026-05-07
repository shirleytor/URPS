# WLS Analysis: Medicaid Expansion and Mortality
#
# MODEL:
#
#   mort_wndr_ac_i = alpha_s(i)
#                  + beta1 * treat_i
#                  + beta2 * smw_i
#                  + beta3 * (treat_i x smw_i)
#                  + epsilon_i
#
#   weighted by w_i = pop_total_i  (county population size)
#   fixed effects:  strata (absorbs k-1 dummy variables)
#   SE:             cluster-robust, clustered at stratum level
#
#   Parameter of interest: beta3 (treat x smw interaction)
#
#   mort_wndr_ac is an aggregate count over population n_i, 
#   Var(epsilon_i) = sigma^2 / n_i. WLS with w_i = n_i corrects this.
#
# PAPER BENCHMARKS (Mann et al. 2024, AOAS):
#   WONDER — SMW counties: delta=0.97 (0.94, 0.99)  p~0.03
#   WONDER — Overall:      delta=1.02 (1.00, 1.06)  n.s.
#
#
# FILES NEEDED:
#   analysis_dataset.csv


library(tidyverse)
library(fixest)

# LOAD AND PREPARE DATA

dat <- read_csv("analysis_dataset.csv", show_col_types = FALSE) |>
  mutate(FIPS = str_pad(FIPS, 5, "left", "0")) |>
  filter(!is.na(strata)) |>
  filter(!is.na(mort_wndr_ac)) |>
  mutate(
    treat  = as.factor(mdcdExp),
    smw    = as.factor(maj_white),
    strata = as.factor(strata)
  )

cat("Dataset Summary\n")
cat(sprintf("In-match (non-suppressed): %d counties\n", nrow(dat)))
cat(sprintf("  Treated: %d,  Control: %d\n",
            sum(dat$mdcdExp == 1), sum(dat$mdcdExp == 0)))
cat(sprintf("  Strata:  %d\n", nlevels(dat$strata)))
cat(sprintf("  SMW counties: %d  (T=%d, C=%d)\n",
            sum(dat$maj_white == 1),
            sum(dat$maj_white == 1 & dat$mdcdExp == 1),
            sum(dat$maj_white == 1 & dat$mdcdExp == 0)))
cat("\n")


# 1. PRIMARY MODEL
#
#   feols(mort_wndr_ac ~ treat + smw + treat:smw | strata,
#         weights = ~pop_total, cluster = ~strata)
#
#   | strata   = strata fixed effects (absorbed, no dummies in output)
#   weights    = population size
#   cluster    = stratum-level cluster-robust SE
cat("1. PRIMARY MODEL\n")
cat("   mort_wndr_ac ~ treat + smw + treat:smw | strata\n")
cat("   weights = pop_total,  cluster = strata\n")
cat("   Parameter of interest: treat1:smw1\n")

fit_primary <- feols(
  mort_wndr_ac ~ treat + smw + treat:smw | strata,
  weights = ~pop_total,
  cluster = ~strata,
  data    = dat
)

summary(fit_primary)

cat("\n Interpretation of treat1:smw1 \n")
cat("  Negative = Medicaid expansion REDUCED mort_wndr_ac more in SMW counties\n")
cat("             than in non-SMW counties  (consistent with paper)\n")
cat("  Positive = expansion increased mortality in SMW relative to non-SMW\n\n")


# 2. SECONDARY OUTCOMES
# Same model applied to all four mortality outcomes


cat("2. SECONDARY OUTCOMES  (treat:smw interaction coefficient)\n")

outcomes <- c(
  mort_wndr_ac        = "AC death count (all)",
  mort_wndr_ac_white  = "AC death count (White)",
  mort_wndr_hca       = "HCA death count (all)",
  mort_wndr_hca_white = "HCA death count (White)"
)

cat(sprintf("  %-30s  %10s  %10s  %8s\n",
            "Outcome", "Estimate", "SE", "p-value"))
cat("  ", paste(rep("-", 62), collapse = ""), "\n", sep = "")

results_secondary <- map_dfr(names(outcomes), function(oc) {
  d <- dat |> filter(!is.na(.data[[oc]]))
  
  fit <- feols(
    as.formula(paste0(oc, " ~ treat + smw + treat:smw | strata")),
    weights = ~pop_total,
    cluster = ~strata,
    data    = d
  )
  
  ct  <- coeftable(fit)
  key <- "treat1:smw1"
  
  if (key %in% rownames(ct)) {
    est <- ct[key, "Estimate"]
    se  <- ct[key, "Std. Error"]
    pv  <- ct[key, "Pr(>|t|)"]
    sig <- ifelse(pv < .001, "***", ifelse(pv < .01, "** ",
                                           ifelse(pv < .05, "*  ", ifelse(pv < .10, ".  ", "   "))))
    
    cat(sprintf("  %-30s  %10.3f  %10.3f  %8.4f  %s\n",
                outcomes[oc], est, se, pv, sig))
    
    tibble(outcome = oc, label = outcomes[oc],
           estimate = est, se = se, p = pv, n = nrow(d))
  }
})


# 3. SMW THRESHOLD SENSITIVITY
#    Replace maj_white (97.5th pctile) with broader definitions


cat("3. SMW THRESHOLD SENSITIVITY\n")
cat("   Outcome: mort_wndr_ac\n")

thresholds <- c(
  maj_white.975 = "97.5th pctile [paper]",
  maj_white.95  = "95th pctile",
  maj_white.9   = "90th pctile",
  maj_white.86  = "86th pctile"
)

cat(sprintf("  %-25s  %5s  %10s  %10s  %8s\n",
            "Threshold", "N_SMW", "Estimate", "SE", "p-value"))
cat("  ", paste(rep("-", 64), collapse = ""), "\n", sep = "")

results_threshold <- map_dfr(names(thresholds), function(thr) {
  d <- dat |>
    mutate(smw_thr = as.factor(.data[[thr]])) |>
    filter(!is.na(mort_wndr_ac))
  
  n_smw <- sum(d[[thr]] == 1, na.rm = TRUE)
  
  fit <- feols(
    mort_wndr_ac ~ treat + smw_thr + treat:smw_thr | strata,
    weights = ~pop_total,
    cluster = ~strata,
    data    = d
  )
  
  ct  <- coeftable(fit)
  key <- "treat1:smw_thr1"
  
  if (key %in% rownames(ct)) {
    est <- ct[key, "Estimate"]
    se  <- ct[key, "Std. Error"]
    pv  <- ct[key, "Pr(>|t|)"]
    sig <- ifelse(pv < .001, "***", ifelse(pv < .01, "** ",
                                           ifelse(pv < .05, "*  ", ifelse(pv < .10, ".  ", "   "))))
    
    cat(sprintf("  %-25s  %5d  %10.3f  %10.3f  %8.4f  %s\n",
                thresholds[thr], n_smw, est, se, pv, sig))
    
    tibble(threshold = thr, label = thresholds[thr], n_smw = n_smw,
           estimate = est, se = se, p = pv)
  }
})


# 4  LINEAR COMBINATION: treat + treat:smw (Additional)
#
# Question: Is the TOTAL effect of expansion on SMW counties negative and statistically significant?
#
# For SMW counties (smw=1), the predicted treatment effect is:
#   E[Y | treat=1, smw=1] - E[Y | treat=0, smw=1] = beta1 + beta3
#
# This is the sum of the main effect of treatment (beta1) and the
# interaction (beta3). We test H0: beta1 + beta3 = 0 using a Wald test.
#
# Var(beta1 + beta3) = Var(beta1) + Var(beta3) + 2*Cov(beta1, beta3)
# The covariance term is critical — cannot just sum the two SEs.

cat("LINEAR COMBINATION TEST: treat + treat:smw\n")
cat("     H0: beta1 + beta3 = 0  (no total treatment effect in SMW counties)\n")

# Manual Wald test using vcov
ct  <- coeftable(fit_primary)
# beta1 (treat main effect)
b1  <- ct["treat1", "Estimate"]
# beta3 (interaction)
b3  <- ct["treat1:smw1", "Estimate"]      
sum_est <- b1 + b3

# cluster-robust variance-covariance matrix
V   <- vcov(fit_primary)                  
se_sum <- sqrt(V["treat1", "treat1"] +
                 V["treat1:smw1", "treat1:smw1"] +
                 2 * V["treat1", "treat1:smw1"])

# Satterthwaite or fixest default df
df  <- degrees_freedom(fit_primary, type = "t")  
t_stat <- sum_est / se_sum
p_val  <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)

cat(sprintf("  beta1 (treat)       = %10.3f\n", b1))
cat(sprintf("  beta3 (treat:smw)   = %10.3f\n", b3))
cat(sprintf("  Sum (beta1 + beta3) = %10.3f\n", sum_est))
cat(sprintf("  SE of sum           = %10.3f\n", se_sum))
cat(sprintf("  t-statistic         = %10.4f\n", t_stat))
cat(sprintf("  df                  = %10.0f\n", df))
cat(sprintf("  p-value             = %10.6f\n", p_val))
cat(sprintf("  95%% CI              = [%.3f, %.3f]\n",
            sum_est - qt(0.975, df) * se_sum,
            sum_est + qt(0.975, df) * se_sum))

sig <- ifelse(p_val < .001, "***", ifelse(p_val < .01, "**",
                                          ifelse(p_val < .05, "*",   ifelse(p_val < .10, ".",   "n.s."))))
cat(sprintf("\n  Conclusion: Sum = %.3f, p = %.4f  (%s)\n", sum_est, p_val, sig))
cat(sprintf("  Direction:  %s\n\n",
            ifelse(sum_est < 0,
                   "NEGATIVE — total treatment effect reduces mortality in SMW counties",
                   "POSITIVE — total treatment effect increases mortality in SMW counties")))

# Cross-check with lincom() if available in fixest version
cat("  --- Cross-check via lincom (if available) ---\n")
tryCatch({
  lc <- lincom(fit_primary, c("treat1 + treat1:smw1" = 1), cluster = ~strata)
  print(lc)
}, error = function(e) {
  cat("  lincom() not available in this fixest version; manual Wald test above is valid.\n")
})



# 5. SUMMARY

cat("SUMMARY\n")

cat("  PAPER BENCHMARK (Mann et al. 2024, WONDER data):\n")
cat("    SMW subgroup:  delta=0.97 (0.94, 0.99)  p~0.03\n")
cat("    Overall:       delta=1.02 (1.00, 1.06)  n.s.\n\n")

cat("  OUR MODEL (treat:smw interaction, mort_wndr_ac):\n")
ct_primary <- coeftable(fit_primary)
if ("treat1:smw1" %in% rownames(ct_primary)) {
  est <- ct_primary["treat1:smw1", "Estimate"]
  se  <- ct_primary["treat1:smw1", "Std. Error"]
  pv  <- ct_primary["treat1:smw1", "Pr(>|t|)"]
  cat(sprintf("    treat:smw  Estimate=%.3f  SE=%.3f  p=%.4f\n", est, se, pv))
  cat(sprintf("    Direction: %s\n\n",
              ifelse(est < 0,
                     "NEGATIVE — consistent with paper (expansion reduced mortality in SMW)",
                     "POSITIVE — opposite to paper")))
}

cat("  METHODOLOGICAL DIFFERENCES vs. PAPER (relevant to RQ3):\n")
cat("    1. Effect scale:    additive count difference vs. multiplicative delta\n")
cat("    2. Test statistic:  WLS t-test vs. aligned rank-sum (ART)\n")
cat("    3. Weights:         pop_total vs. age-standardized m_ki\n")
cat("    4. Suppression:     drop NA rows vs. partial ordering\n")
cat("    5. Covariate adj.:  none (matching handles it) vs. NB regression\n")