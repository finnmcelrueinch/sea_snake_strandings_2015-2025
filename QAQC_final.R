setwd("C:/Users/finnm/OneDrive/Desktop/Uni/Year 4/Dissertation/Data")
library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)
library(janitor)
library(sf)
library(ggplot2)
library(patchwork)
library(forcats)
library(ggtext) 
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(cowplot)
library(spatstat)
library(spatstat.geom)
library(viridis)
library(stars)
library(purrr)
library(dplyr)

master_data <- read_excel("master_data.xlsx")
View(master_data)


#####################
# Quality Assurance #
#####################

#Combine dates so date_found > date_posted. Date_posted only used if date_found NA

clean_date_text <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x == ""] <- NA
  x
}

master_data <- master_data %>%
  mutate(
    date_found  = clean_date_text(date_found),
    date_posted = clean_date_text(date_posted)
  )

parse_text_date <- function(x) {
  
  x <- as.character(x)
  x <- str_trim(x)
  
  out_date <- as.Date(rep(NA, length(x)))
  precision <- rep(NA_character_, length(x))
  
  for (i in seq_along(x)) {
    
    if (is.na(x[i]) || x[i] == "") next
    
    
    xi <- str_replace_all(x[i], "/", "-")
    
    
    if (str_detect(xi, "^\\d{4}-\\d{2}-\\?$")) {
      out_date[i] <- as.Date(paste0(substr(xi, 1, 7), "-01"))
      precision[i] <- "month"
      
      
    } else if (str_detect(xi, "^\\d{4}-\\?-\\?$")) {
      out_date[i] <- as.Date(paste0(substr(xi, 1, 4), "-01-01"))
      precision[i] <- "year"
      
      
    } else if (str_detect(xi, "^\\d{4}-\\d{2}-\\d{2}$")) {
      out_date[i] <- as.Date(xi)
      precision[i] <- "day"
      
      
    } else if (str_detect(xi, "^\\d{1,2}-\\d{1,2}-\\d{4}$")) {
      out_date[i] <- suppressWarnings(dmy(xi))
      precision[i] <- "day"
    }
  }
  
  tibble(date = out_date, precision = precision)
}

found_parsed  <- parse_text_date(master_data$date_found)
posted_parsed <- parse_text_date(master_data$date_posted)

master_data <- master_data %>%
  mutate(
    date_found_parsed  = found_parsed$date,
    date_found_prec    = found_parsed$precision,
    date_posted_parsed = posted_parsed$date,
    date_posted_prec   = posted_parsed$precision
  )

master_data <- master_data %>%
  mutate(
    event_date = coalesce(date_found_parsed, date_posted_parsed),
    date_source = case_when(
      !is.na(date_found_parsed) ~ "found",
      is.na(date_found_parsed) & !is.na(date_posted_parsed) ~ "posted",
      TRUE ~ NA_character_
    ),
    date_precision = case_when(
      date_source == "found"  ~ date_found_prec,
      date_source == "posted" ~ date_posted_prec,
      TRUE ~ NA_character_
    )
  )

summary(master_data$event_date)
table(master_data$date_source)
table(master_data$date_precision)

#Ensure separation of records with date_found  and total strandings 
master_data <- master_data %>%
  mutate(
    found_day_precise =
      date_source == "found" & date_precision == "day",
    
    any_date_available =
      !is.na(event_date)
  )
table(master_data$found_day_precise)
table(master_data$date_precision)
table(master_data$date_source)

#Identify duplicate records (some are as a result of mating or mass strandings, but aren't true duplicates)

n_records <- nrow(master_data)
dup_rows <- master_data %>%
  duplicated() %>%
  sum()
dup_by_core <- master_data %>%
  count(genus, species, event_date, lat, lon) %>%
  filter(n > 1)
list(
  n_records = n_records,
  duplicate_rows = dup_rows,
  duplicate_core_records = nrow(dup_by_core)
)
dup_full <- master_data %>%
  group_by(across(everything())) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  arrange(event_date, stated_location)

nrow(dup_full)
dup_full %>%
  select(event_date, stated_location, genus, species, poster, url)
print(
  dup_full %>% 
    select(event_date, stated_location, genus, species, poster, url),
  n = Inf
)

#Check Taxonomy Consistency (NAs flagged but not incorrect - assigned records with no species as '[genus] sp.' and those without any identified as 'unidentified')
master_data <- master_data %>%
  mutate(
    genus   = str_trim(genus),
    species = str_trim(species),
    
    genus   = str_to_lower(genus),
    species = str_to_lower(species),
    
    
    genus = if_else(
      genus %in% c("na", "n/a", "n.a.", "unknown", "unk", "?"),
      NA_character_,
      genus
    ),
    species = if_else(
      species %in% c("na", "n/a", "n.a.", "unknown", "unk", "?"),
      NA_character_,
      species
    )
  )

master_data <- master_data %>%
  mutate(
    genus   = if_else(!is.na(genus), str_to_sentence(genus), NA_character_),
    species = if_else(!is.na(species), str_to_lower(species), NA_character_)
  )

master_data <- master_data %>%
  mutate(
    genus_missing   = is.na(genus),
    species_missing = is.na(species)
  )

master_data <- master_data %>%
  mutate(
    genus_format_ok = case_when(
      genus_missing ~ NA,
      TRUE ~ str_detect(genus, "^[A-Z][a-z]+$")
    ),
    species_format_ok = case_when(
      species_missing ~ NA,
      TRUE ~ str_detect(species, "^[a-z]+$")
    )
  )
master_data <- master_data %>%
  mutate(
    taxon_label = case_when(
      genus_missing & species_missing ~ "unidentified",
      !genus_missing & species_missing ~ paste0(genus, " sp."),
      !genus_missing & !species_missing ~ paste(genus, species),
      TRUE ~ "unidentified"
    )
  )

table(master_data$taxon_label)


master_data %>%
  filter(
    !is.na(genus),
    !is.na(species),
    species != "sp.",
    species != "sp"
  ) %>%
  distinct(genus, species) %>%
  nrow()
length(unique(master_data$species))
length(unique(master_data$genus))


#Missing coordinates
master_data <- master_data %>%
  mutate(
    lat = as.numeric(lat),
    lon = as.numeric(lon),
    
    coord_missing = is.na(lat) | is.na(lon),
    coord_zero = lat == 0 | lon == 0
  )

#Coordinate range checks
master_data <- master_data %>%
  mutate(
    coord_out_of_range = lat < -45 | lat > -5 | lon < 105 | lon > 165
  )

#QA summaries
table(master_data$coord_missing)
table(master_data$coord_zero)
table(master_data$coord_out_of_range)

missing_coords <- master_data %>%
  filter(coord_missing) %>%
  select(event_date, stated_location, state, lat, lon, poster, url)

print(missing_coords, n = Inf)


out_of_range_coords <- master_data %>%
  filter(coord_out_of_range) %>%
  select(event_date, stated_location, state, lat, lon, poster, url)

out_of_range_coords

#Temporal QA (any records out of date bounds?)
master_data <- master_data %>%
  mutate(
    event_future = event_date > Sys.Date(),
    event_too_old = event_date < as.Date("2014-12-31")
  )

table(master_data$event_future)
table(master_data$event_too_old)

#Logical consistency (i.e. identify non-strandings (0) associated with mortalities (1))
master_data <- master_data %>%
  mutate(
    mortality = as.integer(mortality),
    stranded = as.integer(stranded),
    
    stranded_mortality_conflict =
      stranded == 0 & mortality == 1
  )

table(master_data$stranded_mortality_conflict)
conflicts <- master_data %>%
  filter(stranded_mortality_conflict) %>%
  select(
    event_date,
    genus,
    species,
    stated_location,
    stranded,
    mortality,
    stranding_reason,
    mortality_reason,
    poster,
    url
  )

print(conflicts, n = Inf)

## Save Cleaned dataset, so dont have to run it every time 
write_csv(master_data, "2026-02-23_cleaned_data.csv")
