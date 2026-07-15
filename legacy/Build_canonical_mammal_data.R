###############################################################################
# Build Canonical Mammal Diet Data
#
# Converts the original cat, tegu and rat dietary spreadsheets into the
# canonical long format used throughout the analyses.
#
# This script does not have to be run. 
# It is recorded for legacy and transparency!
###############################################################################

###############################################################################
## Packages
###############################################################################

Require::Require("data.table")
Require::Require("readxl")
Require::Require("openxlsx")

###############################################################################
## Helper functions
###############################################################################

source("functions/helpers.R")

###############################################################################
## Cats
###############################################################################

cats_raw <- extract_diet_data(
  workbook = "data/original/Psiri_Gato.xlsx",
  abundance_sheet = "cats_abund",
  volume_sheet = "cats_vol",
  predator = "Felis_catus"
)

write.csv(cats_raw, "data/intermediate/cats_raw.csv")

###############################################################################
## Tegu
###############################################################################

# I have two sets of samples!

tegu_raw1 <- extract_diet_data(
  workbook = "data/original/Psiri_Teiu.xlsx",
  abundance_sheet = "tegu_abund",
  volume_sheet = "tegu_vol",
  predator = "Salvator_merianae"
)

tegu_raw <- tegu_raw1
# GAVE UP FOR NOW: Volume is giving weird results, not what we need it seems
# Potentially revisit later on
# tegu_raw2 <- extract_new_tegu_data(
#   workbook = "data/original/Psiri_Teiu2.xlsx"
# )
# 
# tegu_raw <- c(tegu_raw1, tegu_raw2)
write.csv(tegu_raw, "data/intermediate/tegu_raw.csv")


###############################################################################
## Rats
###############################################################################

rats_raw <- extract_diet_data(
  workbook = "data/original/Psiri_Rato.xlsx",
  abundance_sheet = "rat_abund",
  volume_sheet = "rat_vol",
  predator = "Rattus_sp"
)

write.csv(rats_raw, "data/intermediate/rats_raw.csv")

