#-----------------------------------------------#
### Strigosella africana Minion 16S           ###
### Cleaning the raw data                     ###
### Date: 02/10/2023                          ###
### Author: Josh Leon                         ###
### Software Versions: R v 4.1.1              ###
###                    tidyverse v 1.3.2      ###
###                    taxize v 0.9.100       ###
#-----------------------------------------------#

# IMPORTS ####
library(tidyverse)
library(taxize)

# DATA CLEANUP ####
## Sample Data ####
sam_data_capreef <- readxl::read_xlsx("./data/Strigosella_sam_data.xlsx") %>% 
  filter(Location == "Capitol Reef")

## Minion Data ####
minion <- read_csv("./data/1st_MINION_files/classification_16s_barcode-v1 (1).csv")
# Getting Just the Capitol Reef Sequences
cap_reffbarcodes <- 
  c("barcode01",
  "barcode02",
  "barcode03",
  "barcode04",
  "barcode05",
  "barcode06",
  "barcode07",
  "barcode08",
  "barcode09",
  "barcode10",
  "barcode11",
  "barcode12")

minion_capreef <- 
  minion %>% 
  filter(barcode %in% cap_reffbarcodes)

# Removing the Word "barcode" in the minion data to allow for joining
minion_capreef <- 
  minion_capreef %>% 
  mutate_at("barcode", str_replace, "barcode", "")
  
# Removing any uneeded columns
minion_capreef <- 
  minion_capreef %>% 
  select(c("barcode", "species"))

# Remove any NAs
minion_capreef <- 
  minion_capreef %>% 
  filter(!is.na(species))

# Seperating the species and genus into two columns
minion_capreef <- 
  minion_capreef %>% 
  separate(col = species,
           into = c("genus", "species"),
           sep = " ")

# FINALIZING THE DATA ####
# Join the Minion Data with the Sample Data
full <- full_join(minion_capreef, sam_data_capreef)

# Getting Phylum, Order, and Family for each sequence based of the genus
unique_genera <- full$genus %>% unique()
t <- tax_name(sci = unique_genera,
              get = c("phylum","order", "family", "genus"), 
              db = "ncbi")

# Adding The Ranks to the Full Dataframe
full <- 
  merge(full, t, by="genus")

# Re-ordering columns and removing any unneeded columns
full$db <- NULL
full$query <- NULL

cleaned_strigosella <- 
  full[, c(2,4,5,6,7,8,9,10,11,12,1,3)]

saveRDS(object = cleaned_strigosella, file = "./data/cleaned_strigosella",)










