# =============================================================================
# WLS Analysis: Medicaid Expansion and Mortality
# uses fixest::feols()
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
# WHY WLS:
#   mort_wndr_ac is an aggregate count over population n_i,
#   so Var(epsilon_i) = sigma^2 / n_i. WLS with w_i = n_i corrects this.
#
# WHY NO COVARIATE ADJUSTMENT:
#   Matching structure (strata FE) already handles covariate balance.
#
# PAPER BENCHMARKS (Mann et al. 2024, AOAS):
#   WONDER — SMW counties: delta=0.97 (0.94, 0.99)  p~0.03
#   WONDER — Overall:      delta=1.02 (1.00, 1.06)  n.s.
#
#
# FILES NEEDED:
#   analysis_dataset.csv
# =============================================================================

library(tidyverse)
library(fixest)

# =============================================================================
# 0. LOAD AND PREPARE DATA
# =============================================================================

dat <- read_csv("analysis_dataset.csv", show_col_types = FALSE) |>
  mutate(FIPS = str_pad(FIPS, 5, "left", "0")) |>
  filter(!is.na(strata)) |>
  filter(!is.na(mort_wndr_ac)) |>
  mutate(
    treat  = as.factor(mdcdExp),
    smw    = as.factor(maj_white),
    strata = as.factor(strata)
  )

cat("=== Dataset Summary ===\n")
cat(sprintf("In-match (non-suppressed): %d counties\n", nrow(dat)))
cat(sprintf("  Treated: %d,  Control: %d\n",
            sum(dat$mdcdExp == 1), sum(dat$mdcdExp == 0)))
cat(sprintf("  Strata:  %d\n", nlevels(dat$strata)))
cat(sprintf("  SMW counties: %d  (T=%d, C=%d)\n",
            sum(dat$maj_white == 1),
            sum(dat$maj_white == 1 & dat$mdcdExp == 1),
            sum(dat$maj_white == 1 & dat$mdcdExp == 0)))
cat("\n")


# =============================================================================
# 1. PRIMARY MODEL
#
#   feols(mort_wndr_ac ~ treat + smw + treat:smw | strata,
#         weights = ~pop_total, cluster = ~strata)
#
#   | strata   = strata fixed effects (absorbed, no dummies in output)
#   weights    = population size
#   cluster    = stratum-level cluster-robust SE
# =============================================================================

cat("======================================================================\n")
cat("1. PRIMARY MODEL\n")
cat("   mort_wndr_ac ~ treat + smw + treat:smw | strata\n")
cat("   weights = pop_total,  cluster = strata\n")
cat("   Parameter of interest: treat1:smw1\n")
cat("======================================================================\n\n")

fit_primary <- feols(
  mort_wndr_ac ~ treat + smw + treat:smw | strata,
  weights = ~pop_total,
  cluster = ~strata,
  data    = dat
)

summary(fit_primary)

cat("\n--- Interpretation of treat1:smw1 ---\n")
cat("  Negative = Medicaid expansion REDUCED mort_wndr_ac more in SMW counties\n")
cat("             than in non-SMW counties  (consistent with paper)\n")
cat("  Positive = expansion increased mortality in SMW relative to non-SMW\n\n")


# =============================================================================
# 2. SECONDARY OUTCOMES
#    Same model applied to all four mortality outcomes
# =============================================================================

cat("======================================================================\n")
cat("2. SECONDARY OUTCOMES  (treat:smw interaction coefficient)\n")
cat("======================================================================\n\n")

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


# =============================================================================
# 3. SMW THRESHOLD SENSITIVITY
#    Replace maj_white (97.5th pctile) with broader definitions
# =============================================================================

cat("\n======================================================================\n")
cat("3. SMW THRESHOLD SENSITIVITY\n")
cat("   Outcome: mort_wndr_ac\n")
cat("======================================================================\n\n")

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


# =============================================================================
# 4. SUMMARY
# =============================================================================

cat("\n======================================================================\n")
cat("4. SUMMARY\n")
cat("======================================================================\n\n")

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
