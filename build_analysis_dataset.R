# =============================================================================
# Build Analysis Dataset
# Output: analysis_dataset.csv
#
# INPUT FILES NEEDED=:
#   1. public_mod_dat.Rdata               — matching structure (from author GitHub)
#   2. supermajority_white_xwalk.Rdata    — SMW indicator (from author)
#   3. supermajority_white_xwalk_05_26.Rdata — SMW thresholds (from author)
#   4. wonder_2014_ac_age_adj.txt         — CDC WONDER: all-cause, all race
#   5. wonder_2014_ac_white_age_adj.txt   — CDC WONDER: all-cause, White
#   6. wonder_2014_hca_age_adj.txt        — CDC WONDER: HCA, all race
#   7. wonder_2014_hca_white_age_adj.txt  — CDC WONDER: HCA, White
#
# This Script:
#   1. Loads matching structure from public_mod_dat.Rdata
#   2. Loads SMW indicators from xwalk files
#   3. Loads and cleans four CDC WONDER mortality files
#   4. Constructs age-adjusted death counts (mort_wndr) per Mann et al. (2024)
#   5. Merges all sources by FIPS
#   6. Saves analysis_dataset.csv
#
# NOTE ON VARIABLES:
#   mort_wndr_X = Deaths x (Age Adjusted Rate / Crude Rate)
#     This converts raw death counts to age-standardized counts,
#     correcting for differences in county age structure.
#     Source: Mann et al. (2024) Section 5.2
#
#   pop_total = working-age (25-64) population from WONDER all-cause query
#   pop_white = working-age (25-64) White population from WONDER White query
#
#   weight for WLS = pop_total (county population size)
#     Rationale: mort_wndr is an aggregate over population n_i,
#     so Var(epsilon_i) = sigma^2 / n_i; WLS with w_i = n_i corrects this.
#     Note: pscore2 in public_mod_dat is derived from restricted data
#     (confirmed by Dr. Mann), so overlap weighting is not appropriate.
# =============================================================================

library(tidyverse)

# =============================================================================
# 1. LOAD MATCHING STRUCTURE
# =============================================================================

load("public.mod.dat.Rdata")   # loads: public.dat

# Key variables:
#   FIPS         — 5-digit county code (stored as integer, will zero-pad)
#   mdcdExp      — treatment indicator (1 = Medicaid expansion by mid-2014)
#   pscore2      — log-odds propensity score (from restricted data, not used)
#   matches.final — matched stratum ID (NA = trimmed/excluded from match)

match_dat <- public.dat |>
  select(FIPS, mdcdExp, pscore2, strata = matches.final) |>
  mutate(FIPS = str_pad(as.character(FIPS), 5, "left", "0"))

cat(sprintf("Matching structure: %d counties, %d in-match, %d strata\n",
            nrow(match_dat),
            sum(!is.na(match_dat$strata)),
            n_distinct(match_dat$strata, na.rm = TRUE)))


# =============================================================================
# 2. LOAD SMW INDICATORS
# =============================================================================

load("supermajority.white.xwalk.Rdata")       # loads: maj.white.xwalk
load("supermajority.white.xwalk_05.26.Rdata") # loads: maj.white.xwalk.upd

# maj.white.xwalk:     FIPS + maj_white (0/1, 97.5th percentile threshold)
# maj.white.xwalk.upd: FIPS + four threshold columns

smw_dat <- maj.white.xwalk |>
  mutate(FIPS = str_pad(as.character(FIPS), 5, "left", "0")) |>
  left_join(
    maj.white.xwalk.upd |>
      mutate(FIPS = str_pad(as.character(FIPS), 5, "left", "0")),
    by = "FIPS"
  )

cat(sprintf("SMW xwalk: %d counties, %d SMW (97.5th pctile)\n",
            nrow(smw_dat), sum(smw_dat$maj_white, na.rm = TRUE)))


# =============================================================================
# 3. LOAD AND CLEAN CDC WONDER FILES
#
#   Each file is tab-delimited with a Notes footer (rows where County is NA).
#   Suppressed values (Deaths <= 10) appear as "Suppressed" — set to NA.
#   "Unreliable" Age Adjusted Rate values — set to NA.
# =============================================================================

read_wonder <- function(path) {
  read_tsv(path, col_types = cols(.default = "c"), show_col_types = FALSE) |>
    filter(!is.na(County)) |>                          # drop Notes footer rows
    filter(`County Code` != "") |>
    rename(
      FIPS            = `County Code`,
      deaths_raw      = Deaths,
      population      = Population,
      crude_rate      = `Crude Rate`,
      age_adj_rate    = `Age Adjusted Rate`
    ) |>
    mutate(
      FIPS         = str_pad(FIPS, 5, "left", "0"),
      deaths_raw   = suppressWarnings(as.numeric(deaths_raw)),   # "Suppressed" -> NA
      population   = suppressWarnings(as.numeric(population)),
      crude_rate   = suppressWarnings(as.numeric(
                       str_replace(crude_rate, "Unreliable", NA_character_))),
      age_adj_rate = suppressWarnings(as.numeric(
                       str_replace(age_adj_rate, "Unreliable", NA_character_))),
      # Age-adjusted death count (Mann et al. 2024, Section 5.2)
      mort_wndr    = deaths_raw * (age_adj_rate / crude_rate)
    ) |>
    select(FIPS, deaths_raw, population, crude_rate, age_adj_rate, mort_wndr)
}

wonder_ac       <- read_wonder("wonder_2014_ac_age_adj.txt")
wonder_ac_white <- read_wonder("wonder_2014_ac_white_age_adj.txt")
wonder_hca      <- read_wonder("wonder_2014_hca_age_adj.txt")
wonder_hca_white<- read_wonder("wonder_2014_hca_white_age_adj.txt")

cat(sprintf("WONDER ac:        %d counties, %d suppressed\n",
            nrow(wonder_ac), sum(is.na(wonder_ac$mort_wndr))))
cat(sprintf("WONDER ac_white:  %d counties, %d suppressed\n",
            nrow(wonder_ac_white), sum(is.na(wonder_ac_white$mort_wndr))))
cat(sprintf("WONDER hca:       %d counties, %d suppressed\n",
            nrow(wonder_hca), sum(is.na(wonder_hca$mort_wndr))))
cat(sprintf("WONDER hca_white: %d counties, %d suppressed\n",
            nrow(wonder_hca_white), sum(is.na(wonder_hca_white$mort_wndr))))


# =============================================================================
# 4. MERGE ALL SOURCES
# =============================================================================

dat <- match_dat |>
  left_join(smw_dat, by = "FIPS") |>
  left_join(
    wonder_ac |> rename(
      pop_total      = population,
      mort_wndr_ac   = mort_wndr,
      rate_ac        = age_adj_rate
    ) |> select(FIPS, pop_total, mort_wndr_ac, rate_ac),
    by = "FIPS"
  ) |>
  left_join(
    wonder_ac_white |> rename(
      pop_white          = population,
      mort_wndr_ac_white = mort_wndr,
      rate_ac_white      = age_adj_rate
    ) |> select(FIPS, pop_white, mort_wndr_ac_white, rate_ac_white),
    by = "FIPS"
  ) |>
  left_join(
    wonder_hca |> rename(
      mort_wndr_hca = mort_wndr,
      rate_hca      = age_adj_rate
    ) |> select(FIPS, mort_wndr_hca, rate_hca),
    by = "FIPS"
  ) |>
  left_join(
    wonder_hca_white |> rename(
      mort_wndr_hca_white = mort_wndr,
      rate_hca_white      = age_adj_rate
    ) |> select(FIPS, mort_wndr_hca_white, rate_hca_white),
    by = "FIPS"
  )


# =============================================================================
# 5. VALIDATION CHECKS
# =============================================================================

cat("\n=== Validation ===\n")
cat(sprintf("Total counties:     %d\n", nrow(dat)))
cat(sprintf("In-match:           %d\n", sum(!is.na(dat$strata))))
cat(sprintf("SMW (97.5th):       %d total, %d in-match\n",
            sum(dat$maj_white == 1, na.rm = TRUE),
            sum(dat$maj_white == 1 & !is.na(dat$strata), na.rm = TRUE)))
cat(sprintf("mort_wndr_ac NAs:   %d (suppressed WONDER outcomes)\n",
            sum(is.na(dat$mort_wndr_ac))))
cat(sprintf("pop_total range:    %d to %d\n",
            min(dat$pop_total, na.rm = TRUE),
            max(dat$pop_total, na.rm = TRUE)))

# Check FIPS coverage
fips_missing <- sum(!dat$FIPS %in% wonder_ac$FIPS)
cat(sprintf("FIPS not in WONDER: %d\n", fips_missing))


# =============================================================================
# 6. SAVE
# =============================================================================

write_csv(dat, "analysis_dataset.csv")

cat("\nSaved: analysis_dataset.csv\n")
cat(sprintf("  %d rows x %d columns\n", nrow(dat), ncol(dat)))
cat("  Columns:", paste(names(dat), collapse = ", "), "\n")
