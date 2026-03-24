# =============================================================================
# Visualization: Medicaid Expansion & Mortality WLS Results
# =============================================================================

library(tidyverse)
library(fixest)
library(ggplot2)
library(patchwork)

# Theme
theme_clean <- function() {
  theme_minimal(base_family = "Arial", base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 15, color = "#1F3864", margin = margin(b = 6)),
      plot.subtitle    = element_text(size = 11, color = "#595959", margin = margin(b = 12)),
      plot.caption     = element_text(size = 9,  color = "#888888", margin = margin(t = 10)),
      axis.title       = element_text(size = 11, color = "#333333"),
      axis.text        = element_text(size = 10, color = "#444444"),
      panel.grid.major = element_line(color = "#EEEEEE"),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position  = "bottom",
      legend.title     = element_text(size = 10),
      legend.text      = element_text(size = 10)
    )
}

BLUE   <- "#1F3864"
RED    <- "#C00000"
GRAY   <- "#888888"
GREEN  <- "#375623"

# Load Data
dat <- read_csv("analysis_dataset.csv", show_col_types = FALSE) |>
  mutate(FIPS = str_pad(FIPS, 5, "left", "0")) |>
  filter(!is.na(strata)) |>
  filter(!is.na(mort_wndr_ac)) |>
  mutate(
    treat  = as.factor(mdcdExp),
    smw    = as.factor(maj_white),
    strata = as.factor(strata)
  )

cat("Data loaded:", nrow(dat), "counties\n")

# Helper: run model and extract treat:smw coefficient
run_model <- function(outcome, data = dat) {
  d <- data |> filter(!is.na(.data[[outcome]]))
  fit <- feols(
    as.formula(paste0(outcome, " ~ treat + smw + treat:smw | strata")),
    weights = ~pop_total,
    cluster = ~strata,
    data    = d
  )
  ct  <- coeftable(fit)
  key <- "treat1:smw1"
  if (!key %in% rownames(ct)) return(NULL)
  tibble(
    outcome  = outcome,
    estimate = ct[key, "Estimate"],
    se       = ct[key, "Std. Error"],
    p        = ct[key, "Pr(>|t|)"],
    ci_lo    = ct[key, "Estimate"] - 1.96 * ct[key, "Std. Error"],
    ci_hi    = ct[key, "Estimate"] + 1.96 * ct[key, "Std. Error"],
    n        = nrow(d)
  )
}


# =============================================================================
# PLOT 1: Coefficient Plot — treat:smw across four outcomes
# FIX: geom_errorbarh → geom_errorbar(orientation="y")
#      ggsave: explicit width=9, height=5 (landscape)
# =============================================================================
cat("Building Plot 1: Coefficient plot...\n")

outcomes_labels <- c(
  mort_wndr_ac        = "All-cause\n(all race)",
  mort_wndr_ac_white  = "All-cause\n(White)",
  mort_wndr_hca       = "Healthcare-amenable\n(all race)",
  mort_wndr_hca_white = "Healthcare-amenable\n(White)"
)

coef_data <- map_dfr(names(outcomes_labels), run_model) |>
  mutate(
    label = outcomes_labels[outcome],
    label = factor(label, levels = rev(outcomes_labels)),
    sig   = ifelse(p < 0.05, "p < 0.05", "p \u2265 0.05")
  )

p1 <- ggplot(coef_data, aes(x = estimate, y = label, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = GRAY, linewidth = 0.7) +
  # FIX 1: geom_errorbar with orientation = "y" instead of geom_errorbarh
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                orientation = "y",
                width = 0.2, linewidth = 0.9) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c("p < 0.05" = RED, "p \u2265 0.05" = BLUE),
    name   = "Significance"
  ) +
  geom_text(aes(label = sprintf("%.0f\n(p=%.2f)", estimate, p)),
            hjust = -0.15, size = 3.2, color = "#333333", lineheight = 0.9) +
  scale_x_continuous(
    limits = c(-800, 200),
    labels = scales::comma
  ) +
  labs(
    title    = "Effect of Medicaid Expansion in SMW Counties",
    subtitle = "treat \u00d7 smw interaction coefficient (WLS with strata FE)",
    x        = "Estimated difference in age-adjusted death count\n(SMW treated vs. SMW control, relative to non-SMW)",
    y        = NULL,
    caption  = "Error bars = 95% CI. Negative = expansion reduced mortality in SMW counties."
  ) +
  theme_clean() +
  theme(legend.position = "right")

# FIX 2: explicit landscape dimensions, no auto-translation
ggsave("plot1_coef.png", p1, width = 9, height = 5, dpi = 150, bg = "white")
cat("  Saved: plot1_coef.png\n")


# =============================================================================
# PLOT 2a & 2b: Bar Charts — split into two separate figures
#   2a: SMW Counties only
#   2b: Non-SMW Counties only
# Each has its own y-axis scale to avoid 30x magnitude distortion in facets.
# =============================================================================
cat("Building Plot 2a/2b: Bar charts (split by subgroup)...\n")

make_bar_data <- function(smw_value) {
  dat |>
    filter(maj_white == smw_value) |>
    mutate(
      treat_lab = ifelse(mdcdExp == 1,
                         "Treated\n(Medicaid expanded)",
                         "Control\n(No expansion)")
    ) |>
    group_by(treat_lab) |>
    summarise(
      mean_mort = weighted.mean(mort_wndr_ac, pop_total, na.rm = TRUE),
      n         = n(),
      .groups   = "drop"
    ) |>
    mutate(
      treat_lab = factor(treat_lab,
                         levels = c("Control\n(No expansion)",
                                    "Treated\n(Medicaid expanded)"))
    )
}

make_bar_plot <- function(data, title_suffix, subtitle_extra = "") {
  ggplot(data, aes(x = treat_lab, y = mean_mort, fill = treat_lab)) +
    geom_col(width = 0.55, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.0f\nn=%d", mean_mort, n)),
              vjust = -0.4, size = 4.0, fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(values = c(
      "Treated\n(Medicaid expanded)" = BLUE,
      "Control\n(No expansion)"      = "#A8C4E0"
    ), guide = "none") +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.18))
    ) +
    labs(
      title    = paste0("Weighted Mean Age-Adjusted Death Count\n", title_suffix),
      subtitle = paste0("All-cause mortality (mort_wndr_ac), weighted by county population",
                        subtitle_extra),
      x        = NULL,
      y        = "Weighted mean age-adjusted death count",
      caption  = "Weights = pop_total. SMW = supermajority White counties (97.5th percentile). n = county count."
    ) +
    theme_clean()
}

# Plot 2a — SMW Counties
bar_smw    <- make_bar_data(smw_value = 1)
p2a <- make_bar_plot(bar_smw,
                     title_suffix  = "— SMW Counties (97.5th percentile)",
                     subtitle_extra = "")

ggsave("plot2a_bar_smw.png", p2a, width = 6, height = 5.5, dpi = 150, bg = "white")
cat("  Saved: plot2a_bar_smw.png\n")

# Plot 2b — Non-SMW Counties
bar_nonsmw <- make_bar_data(smw_value = 0)
p2b <- make_bar_plot(bar_nonsmw,
                     title_suffix  = "— Non-SMW Counties",
                     subtitle_extra = "")

ggsave("plot2b_bar_nonsmw.png", p2b, width = 6, height = 5.5, dpi = 150, bg = "white")
cat("  Saved: plot2b_bar_nonsmw.png\n")


# =============================================================================
# PLOT 3: Threshold Sensitivity — estimate + CI across four SMW thresholds
# (no geometry changes needed here — already uses geom_errorbar vertically)
# =============================================================================
cat("Building Plot 3: Threshold sensitivity...\n")

thresholds <- c(
  maj_white.975 = "97.5th\n(paper, N=306)",
  maj_white.95  = "95th\n(N=499)",
  maj_white.9   = "90th\n(N=832)",
  maj_white.86  = "86th\n(N=997)"
)

thr_data <- map_dfr(names(thresholds), function(thr) {
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
  if (!key %in% rownames(ct)) return(NULL)
  tibble(
    threshold = thresholds[thr],
    n_smw     = n_smw,
    estimate  = ct[key, "Estimate"],
    se        = ct[key, "Std. Error"],
    p         = ct[key, "Pr(>|t|)"],
    ci_lo     = ct[key, "Estimate"] - 1.96 * ct[key, "Std. Error"],
    ci_hi     = ct[key, "Estimate"] + 1.96 * ct[key, "Std. Error"]
  )
}) |>
  mutate(
    threshold = factor(threshold, levels = thresholds),
    sig       = ifelse(p < 0.05, "p < 0.05", "p \u2265 0.05")
  )

p3 <- ggplot(thr_data, aes(x = threshold, y = estimate, color = sig, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = GRAY, linewidth = 0.7) +
  geom_line(color = "#AAAAAA", linewidth = 0.8) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.15, linewidth = 1.0) +
  geom_point(size = 5) +
  scale_color_manual(
    values = c("p < 0.05" = RED, "p \u2265 0.05" = BLUE),
    name   = "Significance"
  ) +
  geom_text(aes(label = sprintf("%.0f\n(p=%.3f)", estimate, p)),
            vjust = -1.1, size = 3.2, color = "#333333", lineheight = 0.85) +
  scale_y_continuous(
    limits = c(-1200, 300),
    labels = scales::comma
  ) +
  labs(
    title    = "SMW Threshold Sensitivity Analysis",
    subtitle = "treat \u00d7 smw interaction estimate as SMW definition broadens",
    x        = "SMW percentile threshold",
    y        = "Estimated treat \u00d7 smw coefficient",
    caption  = "Error bars = 95% CI. Outcome: mort_wndr_ac. Negative = expansion reduced mortality in SMW counties."
  ) +
  theme_clean() +
  theme(legend.position = "right")

ggsave("plot3_threshold.png", p3, width = 9, height = 5.5, dpi = 150, bg = "white")
cat("  Saved: plot3_threshold.png\n")


# =============================================================================
# PLOT 4: Forest Plot — all results combined
# FIX 1: geom_errorbarh → geom_errorbar(orientation = "y")
# FIX 2: ggsave explicit landscape dimensions
# =============================================================================
cat("Building Plot 4: Forest plot...\n")

# Primary results (four outcomes, 97.5th pctile)
primary_results <- map_dfr(names(outcomes_labels), run_model) |>
  mutate(
    row_label = paste0(outcomes_labels[outcome]),
    group     = "Primary (97.5th pctile SMW)"
  )

# Threshold results (mort_wndr_ac, four thresholds)
thr_labels <- c(
  maj_white.975 = "AC (all)  \u2014  97.5th pctile",
  maj_white.95  = "AC (all)  \u2014  95th pctile",
  maj_white.9   = "AC (all)  \u2014  90th pctile",
  maj_white.86  = "AC (all)  \u2014  86th pctile"
)

thr_results <- map_dfr(names(thr_labels), function(thr) {
  d <- dat |>
    mutate(smw_thr = as.factor(.data[[thr]])) |>
    filter(!is.na(mort_wndr_ac))
  fit <- feols(
    mort_wndr_ac ~ treat + smw_thr + treat:smw_thr | strata,
    weights = ~pop_total,
    cluster = ~strata,
    data    = d
  )
  ct  <- coeftable(fit)
  key <- "treat1:smw_thr1"
  if (!key %in% rownames(ct)) return(NULL)
  tibble(
    outcome   = "mort_wndr_ac",
    row_label = thr_labels[thr],
    group     = "Threshold sensitivity",
    estimate  = ct[key, "Estimate"],
    se        = ct[key, "Std. Error"],
    p         = ct[key, "Pr(>|t|)"],
    ci_lo     = ct[key, "Estimate"] - 1.96 * ct[key, "Std. Error"],
    ci_hi     = ct[key, "Estimate"] + 1.96 * ct[key, "Std. Error"],
    n         = nrow(d)
  )
})

forest_data <- bind_rows(primary_results, thr_results) |>
  mutate(
    sig       = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                          p < 0.05  ~ "*",   p < 0.10 ~ ".",   TRUE ~ ""),
    row_label = factor(row_label, levels = rev(unique(row_label))),
    color_grp = ifelse(p < 0.05, "Significant", "Not significant")
  )

p4 <- ggplot(forest_data, aes(x = estimate, y = row_label, color = color_grp)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = GRAY, linewidth = 0.8) +
  # FIX 1: geom_errorbar with orientation = "y"
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                orientation = "y",
                width = 0.3, linewidth = 0.9) +
  geom_point(aes(shape = group), size = 4) +
  geom_text(aes(x = max(ci_hi) + 80,
                label = sprintf("%.0f [%.0f, %.0f] %s", estimate, ci_lo, ci_hi, sig)),
            hjust = 0, size = 3.0, color = "#333333", family = "mono") +
  scale_color_manual(
    values = c("Significant" = RED, "Not significant" = BLUE),
    name   = NULL
  ) +
  scale_shape_manual(
    values = c("Primary (97.5th pctile SMW)" = 16, "Threshold sensitivity" = 15),
    name   = NULL
  ) +
  scale_x_continuous(
    limits = c(-1400, 800),
    labels = scales::comma
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "Forest Plot: All WLS Results",
    subtitle = "treat \u00d7 smw interaction coefficient | mort_wndr_ac unless noted",
    x        = "Estimate (age-adjusted death count difference)",
    y        = NULL,
    caption  = "Error bars = 95% CI.  * p<0.05  ** p<0.01  *** p<0.001"
  ) +
  theme_clean() +
  theme(
    strip.text       = element_text(face = "bold", size = 11, color = BLUE),
    strip.background = element_rect(fill = "#EEF2F8", color = NA),
    legend.position  = "bottom"
  )

# FIX 2: explicit landscape dimensions
ggsave("plot4_forest.png", p4, width = 11, height = 7, dpi = 150, bg = "white")
cat("  Saved: plot4_forest.png\n")

cat("\nAll plots saved:\n")
cat("  plot1_coef.png         — Coefficient plot (4 outcomes)\n")
cat("  plot2a_bar_smw.png     — Bar chart: SMW counties\n")
cat("  plot2b_bar_nonsmw.png  — Bar chart: Non-SMW counties\n")
cat("  plot3_threshold.png    — Threshold sensitivity\n")
cat("  plot4_forest.png       — Forest plot (all results)\n")