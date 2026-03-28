###########################
## Veterinary Case Study ##
###########################


library(tidyverse)
library(janitor)
library(lubridate)
library(survival)
library(survminer)

sea_snakes <- read_csv("vet_data.csv") %>%
  clean_names()

#Standardise and clean dataset
sea_snakes_clean <- sea_snakes %>%
  mutate(
    
    status = status %>% str_trim() %>% str_to_lower(),
    species = species %>% str_trim() %>% str_squish(),
    sex = sex %>% str_trim() %>% str_to_lower() %>% str_squish(),
    maturity = maturity %>% str_trim() %>% str_to_lower() %>% str_squish(),
    suburb = suburb %>% str_trim() %>% str_to_upper(),
    
    
    across(where(is.character), ~na_if(.x, "")),
    
    
    date_admitted = dmy(date_admit),
    status_date = dmy(status_date2),
    
    
    outcome = case_when(
      status %in% c("released", "wildlife released") ~ "survived",
      status %in% c("died in hospital", "euthanised", "doa", "unassisted death", "died in care") ~ "died",
      status %in% c("in hospital", "sent to carer", "under vet care") ~ "ongoing",
      TRUE ~ NA_character_
    ),
    
    mortality = case_when(
      outcome == "survived" ~ 0,
      outcome == "died" ~ 1,
      TRUE ~ NA_real_
    ),
    
    
    maturity = case_when(
      maturity %in% c("sub adult", "subadult") ~ "subadult",
      TRUE ~ maturity
    ),
    
    
    sex = case_when(
      sex %in% c("unknown", "unknown sex") ~ "unknown",
      TRUE ~ sex
    ),
    
    
    species_scientific = case_when(
      species == "Elegant Sea Snake" ~ "Hydrophis elegans",
      species == "Olive-Headed Sea Snake" ~ "Hydrophis major",
      species %in% c("Yellow-Bellied Sea Snake", "Yellow-bellied Sea Snake") ~ "Hydrophis platurus",
      species %in% c("REPTILE - DUBOIS SEA SNAKE", "Dubois Sea Snake", "Dubosis sea snake") ~ "Aipysurus duboisii",
      species == "Stokes Sea Snake" ~ "Hydrophis stokesii",
      species == "Turtle Headed Sea Snake" ~ "Emydocephalus annulatus",
      species == "Small-Headed Sea Snake" ~ "Hydrophis macdowelli",
      species == "Ornate Sea Snake" ~ "Hydrophis ornatus",
      species == "Spectacled Sea Snake" ~ "Hydrophis kingii",
      species == "Horned Sea Snake" ~ "Hydrophis peronii",
      species %in% c("Ocellated Sea Snake", "Spotted sea snake") ~ "Hydrophis ocellatus",
      TRUE ~ NA_character_
    ),
    
    
    categories = str_squish(categories),
    cause_group = case_when(
      categories %in% c("infection/disease", "Infection/Disease") ~ "Infection/Disease",
      categories %in% c("trauma", "Trauma") ~ "Trauma",
      categories %in% c("infection/disease - trauma",
                        "Trauma Infection/Disease",
                        "Trauma  Infection/Disease",
                        "Trauma / Infection/Disease",
                        "Infection/Disease Other",
                        "Infection/Disease  Other",
                        "Trauma Other",
                        "Trauma  Other",
                        "infection/disease - other",
                        "trauma - infection/disease") ~ "Mixed",
      categories %in% c("nad", "NAD", "Other NAD", "Other  NAD") ~ "Unknown",
      categories %in% c("Other", "other") ~ "Other",
      categories == "Orphan/dependant young / NAD" ~ "Dependant Young",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(species))


sea_snakes_outcomes <- sea_snakes_clean %>%
  filter(outcome %in% c("survived", "died"))

#Mortality (%)
overall_mortality <- sea_snakes_outcomes %>%
  summarise(
    total_cases = n(),
    total_died = sum(mortality, na.rm = TRUE),
    total_survived = sum(outcome == "survived", na.rm = TRUE),
    mortality_rate = total_died / total_cases,
    mortality_percent = mortality_rate * 100,
    survival_percent = 100 - mortality_percent
  )

overall_mortality

#Mortality by species (%)
mortality_species <- sea_snakes_outcomes %>%
  group_by(species_scientific) %>%
  summarise(
    total_cases = n(),
    deaths = sum(mortality, na.rm = TRUE),
    survivors = sum(outcome == "survived", na.rm = TRUE),
    mortality_percent = deaths / total_cases * 100,
    survival_percent = survivors / total_cases * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(total_cases))

mortality_species

#Survival Curve
sea_snakes_outcomes <- sea_snakes_outcomes %>%
  mutate(
    days_in_care = as.numeric(status_date - date_admitted),
    outcome = factor(outcome, levels = c("died", "survived"))
  )

time_in_care_summary <- sea_snakes_outcomes %>%
  group_by(outcome) %>%
  summarise(
    cases = n(),
    mean_days = mean(days_in_care, na.rm = TRUE),
    median_days = median(days_in_care, na.rm = TRUE),
    max_days = max(days_in_care, na.rm = TRUE),
    .groups = "drop"
  )

time_in_care_summary

fit_overall <- survfit(
  Surv(days_in_care, mortality) ~ 1,
  data = sea_snakes_outcomes
)

ggsurvplot(
  fit_overall,
  data = sea_snakes_outcomes,
  xlab = "Days in Care",
  ylab = "Survival Probability",
  title = "Survival Probability of Stranded Sea Snakes During Rehabilitation",
  conf.int = TRUE,
  risk.table = TRUE
)

#Cause of affliction
mortality_cause <- sea_snakes_outcomes %>%
  group_by(cause_group) %>%
  summarise(
    total_cases = n(),
    deaths = sum(mortality, na.rm = TRUE),
    survivors = sum(outcome == "survived", na.rm = TRUE),
    mortality_percent = deaths / total_cases * 100,
    survival_percent = survivors / total_cases * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(total_cases))

mortality_cause

#Demographic Analysis
mortality_maturity <- sea_snakes_outcomes %>%
  group_by(maturity) %>%
  summarise(
    total_cases = n(),
    deaths = sum(mortality, na.rm = TRUE),
    survivors = sum(outcome == "survived", na.rm = TRUE),
    mortality_percent = deaths / total_cases * 100,
    survival_percent = survivors / total_cases * 100,
    .groups = "drop"
  )

mortality_maturity

chisq_maturity <- chisq.test(table(sea_snakes_outcomes$maturity,
                                   sea_snakes_outcomes$outcome))
chisq_maturity


mortality_sex <- sea_snakes_outcomes %>%
  group_by(sex) %>%
  summarise(
    total_cases = n(),
    deaths = sum(mortality, na.rm = TRUE),
    survivors = sum(outcome == "survived", na.rm = TRUE),
    mortality_percent = deaths / total_cases * 100,
    survival_percent = survivors / total_cases * 100,
    .groups = "drop"
  )

mortality_sex

chisq_sex <- chisq.test(table(sea_snakes_outcomes$sex,
                              sea_snakes_outcomes$outcome))
chisq_sex

#Seasonal-sex variation
sea_snakes_clean <- sea_snakes_clean %>%
  mutate(
    admit_mo = factor(
      admit_mo,
      levels = c("Jan","Feb","Mar","Apr","May","Jun",
                 "Jul","Aug","Sep","Oct","Nov","Dec"),
      ordered = TRUE
    )
  )

sea_snakes_sex <- sea_snakes_clean

sex_month_counts <- sea_snakes_sex %>%
  count(admit_mo, sex) %>%
  tidyr::complete(
    admit_mo = levels(sea_snakes_sex$admit_mo),
    sex = c("male", "female", "unknown"),
    fill = list(n = 0)
  ) %>%
  mutate(
    admit_mo = factor(
      admit_mo,
      levels = c("Jan","Feb","Mar","Apr","May","Jun",
                 "Jul","Aug","Sep","Oct","Nov","Dec"),
      ordered = TRUE
    )
  )

ggplot(sex_month_counts,
       aes(x = admit_mo, y = n, group = sex, colour = sex)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_manual(
    values = c(
      "male" = "#0072B2",
      "female" = "#D55E00",
      "unknown" = "grey50"
    )
  ) +
  labs(
    title = "Seasonal Pattern of Sea Snake Strandings by Sex",
    x = "Month",
    y = "Number of Strandings",
    colour = "Sex"
  ) +
  theme_classic()

