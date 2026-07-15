
###############################################################################
## PART I. Taxonomic reconstruction
###############################################################################

###############################################################################
## 1. Load packages
###############################################################################

Require::Require("readxl")
Require::Require("data.table")

###############################################################################
## 2. Read input workbook
###############################################################################

input_file <- "data/foodweb_input.xlsx"
sheet_names <- excel_sheets(input_file)

sheets <- lapply(sheet_names, function(x)
  as.data.table(read_excel(input_file, sheet = x))
)
names(sheets) <- sheet_names

diet_raw <- copy(sheets$geral)
setnames(diet_raw, "...1", "prey")
taxon_map <- copy(sheets$taxon_map)
psiri_map <- copy(sheets$psiri_map)
psiri_taxon_map <- copy(sheets$psiri_taxon_map)
final_group_map <- copy(sheets$final_group_map)

# ###############################################################################
# Convert all upper case to lower case
# ###############################################################################
colnames(diet_raw) <- tolower(colnames(diet_raw))
diet_raw$prey <- tolower(diet_raw$prey)

taxon_map$standardized_taxon <- tolower(taxon_map$standardized_taxon)
taxon_map$original_taxon <- tolower(taxon_map$original_taxon)

psiri_map$psiri_taxon <- tolower(psiri_map$psiri_taxon)
psiri_map$standardized_taxon <- tolower(psiri_map$standardized_taxon)

psiri_taxon_map$standard_name <- tolower(psiri_taxon_map$standard_name)
psiri_taxon_map$original_name <- tolower(psiri_taxon_map$original_name)

final_group_map$final_group <- tolower(final_group_map$final_group)
final_group_map$prey <- tolower(final_group_map$prey)

# ###############################################################################
# ## Legacy sheets (used only for validation while rebuilding the pipeline)
# ###############################################################################
# 
# legacy_psiri <- copy(sheets$PSi)
# diet_psiri_mod <- copy(sheets$psiri_mod)
# network_psiri <- copy(sheets$`agora vai`)
# network_final <- copy(sheets$Planilha2)
# diet_harmonized <- copy(sheets$ordem)


###############################################################################
## 3. Source Functions
###############################################################################

source("functions/helpers.R")

###############################################################################
## 4. Taxonomic standardization
###############################################################################
preds <-   c("salvator_merianae", "felis_catus", "rattus_rattus", "rhinella")

diet_collapsed <- standardize_taxa(diet_raw, taxon_map)

###############################################################################
## 5. Aggregate taxa for PSIRI
###############################################################################

stopifnot(!anyDuplicated(diet_collapsed$prey))
stopifnot(all(!is.na(diet_collapsed$prey)))

diet_psiri_input <- aggregate_psiri_taxa(diet_collapsed, psiri_map)

  
###############################################################################
## 6. Build abundance matrix
###############################################################################

abundance_matrix <- copy(diet_psiri_input)

setnames(
  abundance_matrix,
  old = c(
    "teiu",
    "gato",
    "rato",
    "sapo"
  ),
  new = preds
)


###############################################################################
## PART II. PSIRI reconstruction and network analysis
###############################################################################

###############################################################################
## 7. Read original stomach summary tables
###############################################################################

psiri_file <- "data/PSIRI_article.xlsx"
psiri_sheet_names <- excel_sheets(psiri_file)
psiri_sheets <- lapply(
  psiri_sheet_names,
  function(x)
    as.data.table(
      read_excel(
        psiri_file,
        sheet = x
      )
    )
)
names(psiri_sheets) <- psiri_sheet_names

###############################################################################
## Cats
###############################################################################

cat_summary <- read_psiri_block(
  sheet = psiri_sheets$`Psiri gato,teiu, rato`,
  start_col = 1,
  predator = "Felis_catus"
)

###############################################################################
## Tegu
###############################################################################

tegu_summary <- read_psiri_block(
  sheet = psiri_sheets$`Psiri gato,teiu, rato`,
  start_col = 10,
  predator = "Salvator_merianae"
)

###############################################################################
## Rat
###############################################################################

rat_summary <- read_psiri_block(
  sheet = psiri_sheets$`Psiri gato,teiu, rato`,
  start_col = 19,
  predator = "Rattus_rattus"
)

###############################################################################
## Toad
###############################################################################

## NOTE: PN and PV for Rhinella come directly from the table from Felipe.
## They are already rounded to two decimal places by the original authors.

toad_summary <- read_psiri_summary(
  sheet = psiri_sheets$`psiri sapos`,
  predator = "Rhinella"
)

###############################################################################
## Validate imported values tables
###############################################################################

cat_tegu_rat <- rbindlist(
  list(
    cat_summary,
    tegu_summary,
    rat_summary
  ),
  use.names = TRUE
)

## Expected predators
stopifnot(
  identical(
    sort(unique(cat_tegu_rat$predator)),
    sort(c(
      "Felis_catus",
      "Salvator_merianae",
      "Rattus_rattus"
    ))
  )
)

## No missing prey names
stopifnot(
  !any(is.na(cat_tegu_rat$prey))
)

## No duplicated prey within predators
stopifnot(
  !anyDuplicated(
    cat_tegu_rat[
      ,
      .(predator, prey)
    ]
  )
)

## Counts and volumes must be non-negative
stopifnot(all(cat_tegu_rat$count >= 0))
stopifnot(all(cat_tegu_rat$volume >= 0))
stopifnot(all(cat_tegu_rat$fo >= 0))

## Quick summary
cat_tegu_rat[
  ,
  .(
    n_prey = .N,
    total_FO = sum(fo),
    total_count = sum(count),
    total_volume = sum(volume)
  ),
  by = predator
]

stopifnot(
  unique(toad_summary$predator) == "Rhinella"
)

stopifnot(
  !anyDuplicated(toad_summary$prey)
)

stopifnot(
  all(
    c("prey", "predator", "fo", "pn", "pv", "psiri") %in%
      names(toad_summary)
  )
)

###############################################################################
## 8. Combine stomach summary tables
###############################################################################

###############################################################################
## Additional tegu stomachs
##
## When available, add a new worksheet named "tegu_added".
##
## Expected format:
##
## Same layout as the original
## "Psiri gato,teiu, rato" worksheet.
##
## i.e. exactly the same format returned by read_psiri_sheet().
###############################################################################

cat_tegu_rat$prey <- tolower(cat_tegu_rat$prey)
cat_tegu_rat$predator <- tolower(cat_tegu_rat$predator)

toad_summary$prey <- tolower(toad_summary$prey)
toad_summary$predator <- tolower(toad_summary$predator)

diet_summary <- combine_diet_summaries(
  vertebrate_summary = cat_tegu_rat,
  published_summary = toad_summary
)

###############################################################################
## 9. Calculate PSIRI where missing
###############################################################################

diet_summary[
  is.na(psiri),
  psiri := fo * (pn + pv) / 200
]

###############################################################################
## 10. Harmonize interaction matrices
###############################################################################

psiri_matrix <- dcast(
  diet_summary,
  prey ~ predator,
  value.var = "psiri",
  fill = 0
)

setcolorder(
  psiri_matrix,
  c(
    "prey",
    preds
  )
)

###############################################################################
## Validate PSIRI matrix
###############################################################################

stopifnot(!anyDuplicated(psiri_matrix$prey))
stopifnot(
  all(
    names(psiri_matrix) ==
      c(
        "prey",
        preds
      )
  )
)

stopifnot(!anyDuplicated(abundance_matrix$prey))
stopifnot(!anyDuplicated(psiri_matrix$prey))

stopifnot(all(diet_summary$pn >= 0))
stopifnot(all(diet_summary$pv >= 0))
stopifnot(all(diet_summary$psiri >= 0))

diet_summary[
  ,
  .(
    pn = round(sum(pn), 5),
    pv = round(sum(pv), 5)
  ),
  by = predator
]

###############################################################################
## 11. Standardize prey names in PSIRI matrix
###############################################################################

psiri_matrix_std <- merge(
  psiri_matrix,
  psiri_taxon_map,
  by.x = "prey",
  by.y = "original_name",
  all.x = TRUE,
  sort = FALSE
)

## Every prey must have a mapping

stopifnot(!any(is.na(psiri_matrix_std$standard_name)))
psiri_matrix_std[,prey := standard_name]
psiri_matrix_std[,standard_name := NULL]

## Collapse duplicated prey created by harmonization
psiri_matrix_std <- psiri_matrix_std[
  ,
  lapply(.SD, sum, na.rm = TRUE),
  by = prey,
  .SDcols = preds
]

setorder(psiri_matrix_std, prey)
stopifnot(!anyDuplicated(psiri_matrix_std$prey))

###############################################################################
## 12. CANONICAL PSIRI MATRIX
###############################################################################

psiri_matrix <- copy(psiri_matrix_std)
rm(psiri_matrix_std)

###############################################################################
## 13. Aggregate abundance matrix to final prey groups
###############################################################################

abundance_final <- merge(
  abundance_matrix,
  final_group_map,
  by.x = "prey",
  by.y = "prey",
  all.x = TRUE,
  sort = FALSE
)

stopifnot(!any(is.na(abundance_final$final_group)))

abundance_final <- abundance_final[,lapply(.SD, sum, na.rm = TRUE),
  by = final_group,
  .SDcols = preds
]

setorder(
  abundance_final,
  final_group
)

stopifnot(!anyDuplicated(abundance_final$final_group))

###############################################################################
## 14. Build final PSIRI table
###############################################################################

psiri_final <- merge(
  psiri_matrix,
  final_group_map,
  by.x = "prey",
  by.y = "prey",
  all.x = TRUE,
  sort = FALSE
)

stopifnot(!any(is.na(psiri_final$final_group)))

psiri_final <- psiri_final[  ,
  lapply(.SD, sum, na.rm = TRUE),
  by = final_group,
  .SDcols = preds
]

setorder(psiri_final, final_group)

###############################################################################
## Validate PSIRI aggregation
###############################################################################

old_totals <- psiri_matrix[
  ,
  lapply(.SD, sum),
  .SDcols = preds
]

new_totals <- psiri_final[
  ,
  lapply(.SD, sum),
  .SDcols = preds
]

stopifnot(
  isTRUE(all.equal(
    unlist(old_totals),
    unlist(new_totals),
    tolerance = 1e-10
  ))
)
###############################################################################
## Validate ABUNDANCE aggregation
###############################################################################

old_totals <- abundance_matrix[,
  lapply(.SD, sum, na.rm = TRUE),
  .SDcols = preds
]

new_totals <- abundance_final[,
  lapply(.SD, sum),
  .SDcols = preds
]

stopifnot(
  isTRUE(all.equal(
    unlist(old_totals),
    unlist(new_totals),
    tolerance = 1e-10
  ))
)

###############################################################################
## 15. Setting proper column names
###############################################################################

setnames(
  psiri_final,
  "final_group",
  "prey"
)
setorder(psiri_final, prey)

setnames(
  abundance_final,
  "final_group",
  "prey"
)
setorder(abundance_final, prey)

stopifnot(
  identical(
    abundance_final$prey,
    psiri_final$prey
  )
)

###############################################################################
## 16. Calculate density-weighted predator consumption
###############################################################################

samples <- c("salvator_merianae" = 22, # tegu (n=22) Gaiotto et al., 2020
             "felis_catus" = 78, # cats (n=78) Gaiotto et al., 2020
             "rattus_rattus" = 10, # rats (n=10) Gaiotto et al., 2020
             "rhinella" = 66) # ???? toads? (n=143) Tolledo & Toledo (2015) Not matching the text!

densities <- c("salvator_merianae" = 3.98, # tegu 3.98 ind/ha Abrahão et al. 2019
               "felis_catus" = 0.71, # feral cats, 0.71 ind/ha , Dias et al. 2017
               "rattus_rattus" = 37, # rats 37 ind/ha, Russell et al. 2018
               "rhinella" = 10.35) # Toads, extrapolation from Solomon Islands, Pikacha et al. 2015

predator_consumption <- copy(abundance_final)
predators <- names(samples)

predator_consumption[,
  (predators) := Map(`/`, .SD, samples),
  .SDcols = predators
]

predator_consumption[,
  (predators) := Map(`*`, .SD, densities),
  .SDcols = predators
]

###############################################################################
## 17. Define food web nodes
###############################################################################

resource_nodes <- abundance_final$prey

predator_nodes <- colnames(abundance_final)[!colnames(abundance_final)=="prey"]

## EXCLUDE RATS as they are also prey
resource_nodes <- resource_nodes[!resource_nodes == "rattus_rattus"]

food_web_nodes <- c(
  resource_nodes,
  predator_nodes
)

###############################################################################
## 18. Create food web matrix
###############################################################################

food_web_matrix <- matrix(
  0,
  nrow = length(food_web_nodes),
  ncol = length(food_web_nodes),
  dimnames = list(
    food_web_nodes,
    food_web_nodes
  )
)

food_web_matrix[
  psiri_final$prey,
  c(
    "salvator_merianae",
    "felis_catus",
    "rattus_rattus",
    "rhinella"
  )
] <- as.matrix(
  psiri_final[, .(
    salvator_merianae,
    felis_catus,
    rattus_rattus,
    rhinella
  )]
)

# ###############################################################################
# ## 18. Build interaction matrix
# ###############################################################################

