#!/usr/bin Rscript

## This script uses the metadata downloaded from ATB and categorizes the host species based on categories. The categories are:
##Meat producing poultry, Adult cattle, layers, veal calves, pets, small ruminants, pigs, surface water and water birds, wild birds, non-water environment and human
##
## usage:
## bin/organise_ATB-metadata
## Run in the base campylobacter-crisprscape folder!
##
## TODO: Does not account for if the date of the metadata changes.
## 
## 


print("Enabling packages...")
library(tidyverse)

library(naniar)


print("Importing...")
#ENA metadata is downloaded from ATB and was filtered for original analysis in descriptive_statistics
colnames <- read_table("doc/colnames_metadata_filtered.txt", col_names = c("one", "two"), show_col_types = F)
metadata_source <-
  read_tsv("data/ATB/ena_metadata.20240801-filtered.tsv.gz",
           col_names = colnames$two, show_col_types = FALSE)

#further filtering metadata as not all are needed for this analysis
metadata <- select(metadata_source, c(sample_accession, country, scientific_name, tax_id, isolation_source, host))

#in the metadata there is a distinction between unclear data and data that is simply missing/unspecified. missing data is turned into NA and data that is present in isolation source but not in host is copied to host.
missingno <- metadata %>%
  mutate(across(c(isolation_source, host), ~str_to_lower(.x)))%>%
  replace_with_na(replace = list(isolation_source = c(
    "missing", "other", "not collected", "no source specified", "mising"
  ), host = c(
    "missing", "other", "not collected", "no source specified", "mising"
  ))) %>% mutate(host = ifelse(is.na(host), isolation_source, host))

print("constructing categories...")
#consolidates the many differing inputs of host into more managable categories: Meat producing poultry, Adult cattle, layers, veal calves, pets, small ruminants, pigs, surface water and water birds, wild birds, non-water environment and human. 
#the last category (other/undetermined) consists of all species not part of these initial categories, or unclear enough to not be useful.
definedmeta <- missingno %>%
  mutate(category = host) %>%
  mutate(
    category = str_replace(
      category,
      ".*(field|environment|pasture|sediment|soil|surfaces|crates).*",
      "dry environment"
    )
  ) %>%
  mutate(
    category = str_replace(
      category,
      ".*(calf|veal).*",
      "calves"
    )
  ) %>%
  mutate(
    category = str_replace(
      category,
      ".*(taurus|cattle|milk|boi?vine|cow|heifer|steer|beef|dairy|indicus).*",
      "adult cattle"
    )
  ) %>% 
  mutate(
    category = ifelse(str_detect(
      isolation_source,
      ".*((?<!v)egg|layer).*"), yes = "egg layer", no = host)
    ) %>%
  mutate(category = str_replace(category, ".*(water|lagoon|sewage|river|wetland|fulica\\satra|platyrhynchos|duck).*", "surface water and water birds")) %>%
  mutate(
    category = str_replace(
      category,
      ".*(pheasant|phasianus|gallopavo|hen|chi[ec]ken|turkey|broiler|gallus|cb-|drumsticks|poultry).*",
      "meat producing poultry"
    )
  ) %>%
  mutate(category = str_replace(
    category,
    ".*(sus\\sscrofa|porcine|pig|swine|sow|pork).*",
    "swine/pig"
  )) %>%
  mutate(
    category = str_replace(
      category,
      ".*(puppy|canii?ne|cat|feline|dog|kitten|pet|familiaris)(?!tle).*",
      "pet animals"
    )
  ) %>%
  mutate(category = str_replace(category, ".*(human|clinical|guillain|sapiens).*", "human")) %>%
  mutate(
    category = str_replace(
      category,
      ".*(avian|cloaca|(?<!water\\s)bird|columba livia|crow|corvus|dove).*",
      "wild avian"
    )
  ) %>%
  mutate(category = str_replace(category, ".*(sheep|aries|ovine|goat|hircus).*", "small ruminants")) %>%
  mutate(category = str_replace(category, ".*(human|sapiens).*", "human")) %>%
  mutate(
    category = str_replace(
      category,
      "^(?!(dry environment|calves|adult cattle|egg layer|meat producing poultry|swine/pig|pet animals|wild avian|small ruminants|food|feces|human|surface water and water birds)).*",
      "other/undetermined"
    ))





write_tsv(definedmeta, "doc/filtered_ena_metadata.tsv")
print("done!")