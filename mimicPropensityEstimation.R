#Script to calculate the Hajeck estimator for treatment effect based on approximated
#propensity scores

#Load the data, assumes the base directory is URPS/analysis
pub_data <- local(get(load("../data/public.mod.dat.Rdata")))

#Clean the data to get rid of the unwanted values, also makes mort_wndr which we will
#renormailize later, but I didn't want to mess any of the calculation up
age_adj <- read.delim(file = "../data-raw/wonder/wonder_2014_ac_age_adj.txt")
age_adj <- age_adj %>%
  mutate(across(Deaths:Age.Adjusted.Rate, ~ as.numeric(case_when(.x %in% c("Suppressed", "Unreliable", "Missing") ~ "",
                                                                 TRUE ~ .x))), 
         Population = as.numeric(Population)) %>%
  mutate(FIPS = sprintf('%05d', County.Code),
         adj = Age.Adjusted.Rate/Crude.Rate,
         mort1 = Age.Adjusted.Rate/100000*Population,
         mort = Deaths*adj) %>%
  mutate(mort_wndr = case_when(is.na(Crude.Rate) & !is.na(Deaths) ~ Deaths,
                               TRUE ~ round(mort))) %>%
  select(FIPS, mort_wndr, Population) %>%
  drop_na(mort_wndr, Population)

join_data <- age_adj %>%
  inner_join(pub_data, by = "FIPS") %>%
  mutate(matches.final = matches.final) %>%
  drop_na(matches.final)

#Load white populations from WNDR
white_age_adj <- read.delim(file = "../data-raw/wonder/wonder_2014_ac_white_age_adj.txt")
white_age_adj <- white_age_adj %>%
  mutate(across(Deaths:Age.Adjusted.Rate, ~ as.numeric(case_when(.x %in% c("Suppressed", "Unreliable", "Missing") ~ "",
                                                                 TRUE ~ .x)))) %>%
  mutate(FIPS = sprintf('%05d', County.Code)) %>%
  select(FIPS, WPopulation = Population)

join_data <- join_data %>%
  left_join(white_age_adj, by = "FIPS")

#Calculate the super majority white counties by taking a high percentile of the 
#proportions of white people in each county.
join_data$WhitePercent <- join_data$WPopulation / join_data$Population
join_data$SMW <- join_data$WhitePercent >= quantile(join_data$WhitePercent, 0.985)

temp_matches.final <- join_data[["matches.final"]]

#Finds the number of treated groups in each match.
treatNumInGroup <- join_data %>%
  filter(mdcdExp == 1) %>%
  count(matches.final) %>%
  mutate(matches.final = as.character(matches.final)) %>%
  rename(treated_n = n)
#Finds number of total elements of each match
numsPerGroup <- join_data %>%
  count(matches.final) %>%
  rename(Freq = n) %>%
  mutate(matches.final = as.character(matches.final)) %>%
  left_join(treatNumInGroup, by = "matches.final") %>%
  mutate(treated_n = replace_na(treated_n, 0))

#Calculates the propensity score estimates
fauxPScoresExp <- numsPerGroup$treated_n / numsPerGroup$Freq
fauxPScoresNoExp <- 1 - fauxPScoresExp

fauxScoresDat <- data.frame(
  scoresExp = fauxPScoresExp,
  scoresNoExp = fauxPScoresNoExp,
  matches.final = numsPerGroup[["matches.final"]]
)

fauxScoresDat <- fauxScoresDat %>% 
  left_join(join_data %>% 
              mutate(matches.final = as.character(matches.final)) %>%
              select(matches.final, SMW, Population, mort_wndr, mdcdExp, FIPS), 
            by = "matches.final")

fauxScoresDat <- fauxScoresDat %>% drop_na(mort_wndr)

fauxScoresDat <- fauxScoresDat %>% mutate(mort_wndr = mort_wndr / Population * 100000 )

fauxScoresDat <- fauxScoresDat %>%
  filter(scoresExp > 0, scoresNoExp > 0)

#fauxWeights <- (fauxScoresDat$Population) / (fauxScoresDat$mdcdExp * fauxScoresDat$scoresExp + (1 - fauxScoresDat$mdcdExp) * fauxScoresDat$scoresExp)

Hajeck_Est <- fauxScoresDat %>%
  filter(SMW == TRUE) %>%
    summarise(
      treated = sum((mdcdExp * mort_wndr) / scoresExp) / sum(mdcdExp / scoresExp),
      control = sum(((1-mdcdExp) * mort_wndr) / scoresNoExp) / sum((1-mdcdExp) / scoresNoExp),
      ATE = treated - control
  )


