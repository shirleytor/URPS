pub_data <- local(get(load("../data/public.mod.dat.Rdata")))

age_adj <- read.delim(file = "../data-raw/wonder/wonder_2014_ac_age_adj.txt")
age_adj <- age_adj %>%
  mutate(across(Deaths:Age.Adjusted.Rate, ~ as.numeric(case_when(.x %in% c("Suppressed", "Unreliable", "Missing") ~ "",
                                                                 TRUE ~ .x)))) %>%
  mutate(FIPS = sprintf('%05d', County.Code),
         adj = Age.Adjusted.Rate/Crude.Rate,
         mort1 = Age.Adjusted.Rate/100000*Population,
         mort = Deaths*adj) %>%
  mutate(mort_wndr = case_when(is.na(Crude.Rate) & !is.na(Deaths) ~ Deaths,
                               TRUE ~ round(mort))) %>%
  select(FIPS, mort_wndr, Population)

join_data <- pub_data %>%
  left_join(age_adj, by = "FIPS") %>%
  mutate(matches.final = matches.final) %>%
  drop_na(matches.final)

temp_matches.final <- as.data.frame(join_data[["matches.final"]])

numsPerGroup <- as.data.frame(table(temp_matches.final))
treatNumInGroup <- join_data %>%
  filter(mdcdExp == 1) %>%
  count(matches.final)

fauxPScoresExp <- treatNumInGroup[["n"]] / numsPerGroup[["Freq"]]
fauxPScoresNoExp <- 1 - (treatNumInGroup[["n"]] / numsPerGroup[["Freq"]])

fauxScoresDat <- data.frame(
   scoresExp = fauxPScoresExp,
   scoresNoExp = fauxPScoresNoExp,
   matches.final = treatNumInGroup[["matches.final"]]
)
fauxScoresDat <- fauxScoresDat %>% 
  left_join(join_data %>% select(matches.final,Population,mort_wndr,mdcdExp,FIPS), by = join_by(matches.final))

fauxWeights <- (fauxScoresDat$Population) / (fauxScoresDat$mdcdExp * fauxScoresDat$scoresExp + (1 - fauxScoresDat$mdcdExp) * fauxScoresDat$scoresExp)




