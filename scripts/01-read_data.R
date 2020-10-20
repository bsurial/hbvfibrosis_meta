library(tidyverse)
library(readxl)
library(here)


df <- read_excel("data/fib_hbv_extraced_dominik_update_0820.xlsx", .name_repair = "universal",
                 na = c("x", "NA"))

df <- df %>% 
  filter(!is.na(author)) %>% 
  select(-starts_with("empty"))


# Function to clean labels from excel

clean_string <- function(string) {
  string <- str_remove_all(string,"[:digit:]+:? ")
  string <- str_remove_all(string, " {2,}")
  string <- str_remove(string, "^ ")
  string <- str_replace_all(string, "; ", ";")
  str_split(string, ";")
}

l_studypop2 <- c("Referral/Teaching hospital", "Primary care/General population")
l_setting <- c("urban", "rural", "both", "unable to tell")
l_institution <- c("teaching (university) hospital", 
                   "other hospital", "health center", "community",
                   "not specified", "prison and referral to hospital", 
                   "combination of hospital and other settings")
l_popluation_type <- c(" hospitalized", "oupatient clinic", "VCT center", "blood donors", 
                       "prisoners", "prison officers", " pregnant women (PMTCT ANC)", 
                       "household or other general population surveys", "female sex workers", 
                       "drug users", "STI clinic or other outpatient clinic", 
                       "HIV-infected children's home", 
                       "MSM ", "medicolegal specimens (bodies)", "occupational group")

l_study_type <- "1: randomized controlled trial (RCT); 2: prospective cohort study; 3: case-control study; 4: cross-sectional study; 5: chart review"
l_study_type <- clean_string(l_study_type)[[1]]

l_age_group <- "1: adults; 2: children; 3: adolescents; 4: adults and children; 5: adults and adolescents; 6: children and adolescents"
l_age_group <- clean_string(l_age_group)[[1]]

l_sex <- "1: both; 2: men only; 3: women only"
l_sex <- clean_string(l_sex)[[1]]

df <- df %>% 
  mutate(setting = factor(setting, levels = 1:4,  
                          labels = l_setting),
         institution = factor(institution, levels = 1:7,
                              labels = l_institution),
         population_type = factor(population_type, levels = 1:15, 
                             labels = l_popluation_type),
         study_type = factor(study_type, levels = 1:length(l_study_type), 
                             labels = l_study_type),
         age_group = factor(age_group, levels = 1:length(l_age_group), 
                            labels = l_age_group),
         sex = factor(sex, levels = 1:length(l_sex), 
                      labels = l_sex),
         studypop2 = factor(studypop2, levels = 1:length(l_studypop2), 
                            labels = l_studypop2))



df %>% 
  select(setting, institution, population_type, study_type, age_group, sex, studypop2) %>% 
  summary()


## Regroup categories

df <- df %>% 
  mutate(institution_simple = fct_collapse(institution,
    "teaching hospital" = "teaching (university) hospital",
    "other hospital or health center" = c("other hospital", "health center"),
    "prison" = "prison and referral to hospital")) %>% 
  mutate(population_type_simple = fct_collapse(population_type, 
    "outpatient clinic" = "oupatient clinic", 
    "population screening (incl. blood donors and occ groups)" = c("blood donors", 
                                                                  "household or other general population surveys",
                                                                   "occupational group"),
    "prisoners" = "prisoners")) %>% 
  mutate(setting_simple = fct_collapse(setting, 
    "urban" = c("urban", "both"),
    "rural" = "rural"))


df %>% 
  count(setting_simple, setting)

Hmisc::describe(df %>% select(1:prop_fibrosis_assessed, prop_chb_female, prop_hbe_ag_tested_positive, cutoff_fibrosis, cutoff_cirrhosis))

median(df$prop_hbe_ag_tested_positive, na.rm = TRUE)

df %>% 
  filter(prop_chb_female > 1) %>% 
  select(author, n_chb, n_chb_female)


write_rds(df, here("processed", "01-hbv_fibrosis.rds"))


