
### EXPLORATION
## Concluded: Só pra explicar um pouco melhor o pedido da Ju: estou reorganizando 
# toda a análise da rede trófica para deixá-la completamente reproduzível 
# (desde os dados brutos até a matriz final utilizada no modelo). Estou fazendo 
# isso porque encontrei alguns errinhos na tabela final e quero ter certeza que 
# o paper esta correto. Consegui reconstruir praticamente toda a pipeline: 
# padronização da taxonomia, agregação das presas, construção da matriz de 
# interações e todas as etapas posteriores. O único ponto que ainda não consegui 
# reproduzir é o cálculo do PSIRI (por isso o pedido da Ju). Revisei todas as
# planilhas do arquivo original e elas contêm apenas os valores finais de PSIRI 
# (já calculados), além dos dados agregados de FO, NO e Volume. Pelo que 
# consegui verificar, como a Ju falou, esses dados agregados não são suficientes 
# para reproduzir exatamente o cálculo do PSIRI, pois esse índice depende das
# informações por estômago (ou por amostra individual). Vocês sabem se ainda 
# existem as planilhas originais com os dados de cada estômago (ou de cada 
# amostra), ou o arquivo utilizado para calcular o PSIRI antes desses valores 
# serem copiados para a planilha atual? Se esses arquivos ainda existirem, isso
# permitirá reconstruir toda a análise de forma totalmente reproduzível. Obrigada!


###############################################################################
## Validate PSIRI Mapping
###############################################################################

setdiff(diet_psiri_input$prey, legacy_psiri$PSIRI)
setdiff(legacy_psiri$PSIRI, diet_psiri_input$prey)

# NOTE: Reptiles_(scale) was retained as Squamata during taxonomic
# standardization. The original Excel workflow removed these records,
# but no rationale was documented.

# Reconstruct PSIRI
setnames(
  legacy_psiri,
  c(
    "PSIRI",
    "Felis_catus...2",
    "Rattus_rattus...3",
    "Salvator_merianae...4",
    "Rhinella...5",
    "...6",
    "...7",
    "ABUN",
    "Salvator_merianae...9",
    "Felis_catus...10",
    "Rattus_rattus...11",
    "Rhinella...12"
  ),
  c(
    "prey",
    "psiri_cat",
    "psiri_rat",
    "psiri_tegu",
    "psiri_toad",
    "blank1",
    "blank2",
    "prey_count_name",
    "count_tegu",
    "count_cat",
    "count_rat",
    "count_toad"
  )
)

legacy_psiri[, c("blank1", "blank2") := NULL]
legacy_psiri

raw_psiri <- fread("data/raw_data_NO_FO.csv")

str(raw_psiri)

unique(raw_psiri$`Species Evaluated`)

unique(raw_psiri$`Species Eaten`)
setnames(
  raw_psiri,
  old = c(
    "Species Eaten",
    "Species Evaluated",
    "FO",
    "NO",
    "Volume"
  ),
  new = c(
    "prey",
    "predator",
    "fo",
    "count",
    "volume"
  )
)

raw_psiri <- raw_psiri[
  !prey %in% c("Teiu", "Gato", "Rato", "Sapo")
]

raw_psiri[is.na(fo), fo := 0]
raw_psiri[is.na(count), count := 0]
raw_psiri[is.na(volume), volume := 0]

## TEST PSIRI CALC
tegu <- raw_psiri[predator == "Teiu"]
sum(tegu$count)
sum(tegu$volume)
tegu[
  ,
  `:=`(
    PN = 100 * count / sum(count),
    PV = 100 * volume / sum(volume)
  )
]
tegu[
  ,
  PSIRI_test := fo * (PN + PV) / 200
]
tegu

tegu <- copy(raw_psiri[predator == "Teiu"])

## Replace missing values
tegu[is.na(fo), fo := 0]
tegu[is.na(count), count := 0]
tegu[is.na(volume), volume := 0]

## Totals
N_total <- sum(tegu$count)
V_total <- sum(tegu$volume)

## Candidate metrics
tegu[, `:=`(
  FO_percent = fo,
  FO_prop    = fo / 100,
  NO_percent = 100 * count / N_total,
  NO_prop    = count / N_total,
  VOL_percent = 100 * volume / V_total,
  VOL_prop    = volume / V_total
)]
tegu[, `:=`(
  
  cand1 = FO_percent,
  
  cand2 = NO_percent,
  
  cand3 = VOL_percent,
  
  cand4 = FO_percent * NO_percent / 100,
  
  cand5 = FO_percent * VOL_percent / 100,
  
  cand6 = FO_percent * (NO_percent + VOL_percent) / 200,
  
  cand7 = FO_prop * (NO_percent + VOL_percent) / 2,
  
  cand8 = (NO_percent + VOL_percent) / 2
  
)]
legacy_tegu <- legacy_psiri[
  ,
  .(
    prey,
    psiri = psiri_tegu
  )
]

compare <- merge(
  tegu,
  legacy_tegu,
  by = "prey",
  all.x = TRUE
)

compare
for(i in paste0("cand",1:8)){
  
  cat(
    i, " : ",
    cor(compare[[i]], compare$psiri,
        use = "complete.obs"),
    "\n"
  )
  
}

###############################################################################
## Compare abundance matrix with quadratica2
###############################################################################

legacy_freq <- copy(sheets$quadratica2)
legacy_freq
freq_matrix <- copy(diet_collapsed)

setnames(
  freq_matrix,
  old = c("Teiu","Gato","Rato","Sapo"),
  new = c(
    "Salvator_merianae",
    "Felis_catus",
    "Rattus_rattus",
    "Rhinella"
  )
)

freq_matrix
setdiff(freq_matrix$prey, legacy_psiri$prey)
setdiff(legacy_psiri$prey, freq_matrix$prey)
freq_matrix[prey == "Trachylepis_atlantica"]
freq_matrix[prey == "Seeds"]

tegu_compare <- merge(
  tegu[, .(
    prey,
    fo,
    count,
    volume
  )],
  legacy_tegu,
  by = "prey",
  all = TRUE
)

tegu_compare[order(-psiri)]

n_tegu <- 22

tegu[, FO_prop := fo / n_tegu]

tegu[, PN := count / sum(count)]

tegu[, PV := volume / sum(volume)]

tegu[, PSIRI_manuscript := 100 * FO_prop * (PN + PV) / 2]
raw_psiri[fo > 22]

tegu[, PSIRI_test2 := (fo / 100) * (PN + PV) / 2 * 100]

merge(
  tegu[, .(prey, PSIRI_test2)],
  legacy_tegu,
  by = "prey",
  all = TRUE
)[order(-psiri)]
tegu[, PSIRI_test := fo * (PN + PV) / 2]
compare <- merge(
  tegu[, .(prey, PSIRI_test)],
  legacy_tegu,
  by = "prey"
)

compare


###############################################################################
## 6. Prepare legacy PSIRI table
###############################################################################

# The original workbook stores the final PSIRI values but not the calculations
# used to obtain them. Until the original stomach-level data are recovered,
# these values are treated as validated legacy input.

setnames(
  legacy_psiri,
  old = c(
    "PSIRI",
    "Felis_catus...2",
    "Rattus_rattus...3",
    "Salvator_merianae...4",
    "Rhinella...5",
    "...6",
    "...7",
    "ABUN",
    "Salvator_merianae...9",
    "Felis_catus...10",
    "Rattus_rattus...11",
    "Rhinella...12"
  ),
  new = c(
    "prey",
    "psiri_cat",
    "psiri_rat",
    "psiri_tegu",
    "psiri_toad",
    "blank1",
    "blank2",
    "prey_name",
    "count_tegu",
    "count_cat",
    "count_rat",
    "count_toad"
  )
)

legacy_psiri[, c("blank1", "blank2") := NULL]

###############################################################################
## 7. Correct documented legacy typos
###############################################################################

legacy_psiri[
  prey == "Chepalopoda",
  prey := "Cephalopoda"
]

legacy_psiri[
  prey == "Decapoda",
  prey := "Brachyura"
]

legacy_psiri[
  prey == "Rattus_spp._",
  prey := "Rattus_rattus"
]

###############################################################################
## 8. Validate taxonomy
###############################################################################

missing_taxa <- setdiff(
  legacy_psiri$prey,
  diet_psiri_input$prey
)

extra_taxa <- setdiff(
  diet_psiri_input$prey,
  legacy_psiri$prey
)

cat("Taxa only in legacy PSIRI:\n")
print(missing_taxa)

cat("\nTaxa only in reconstructed table:\n")
print(extra_taxa)

###############################################################################
## 11. Validate reconstructed abundance matrix
###############################################################################

legacy_abundance <- copy(
  legacy_psiri[
    ,
    .(
      prey,
      Felis_catus       = count_cat,
      Rattus_rattus     = count_rat,
      Salvator_merianae = count_tegu,
      Rhinella          = count_toad
    )
  ]
)

comparison <- merge(
  freq_matrix,
  legacy_abundance,
  by = "prey",
  suffixes = c("_reconstructed", "_legacy"),
  all = TRUE
)

comparison[
  Salvator_merianae_reconstructed != Salvator_merianae_legacy |
    Felis_catus_reconstructed != Felis_catus_legacy |
    Rattus_rattus_reconstructed != Rattus_rattus_legacy |
    Rhinella_reconstructed != Rhinella_legacy
]

## Expected discrepancies
## --------------------
## - Decapoda: legacy copy/paste error
## - Diptera: reconstructed table collapses larvae + pupae
## - Coleoptera: reconstructed table includes Hyperaulax_ridleyi
## - Squamata: retained pending clarification

expected_differences <- c(
  "Brachyura",
  "Cephalopoda",
  "Coleoptera",
  "Diptera",
  "Gastropoda",
  "Rattus_rattus",
  "Squamata (temporary)"
)

unexpected <- setdiff(
  abundance_matrix$prey,
  c(legacy_abundance$prey, expected_differences)
)

stopifnot(length(unexpected) == 0)

###############################################################################
## Investigate legacy PSIRI calculation
###############################################################################

diet_summary_compare <- copy(diet_summary)

diet_summary_compare[
  predator == "Teiu",
  predator := "Salvator_merianae"
]

diet_summary_compare[
  predator == "Gato",
  predator := "Felis_catus"
]

diet_summary_compare[
  predator == "Rato",
  predator := "Rattus_rattus"
]

diet_summary_compare[
  predator == "Sapo",
  predator := "Rhinella"
]

legacy_long <- melt(
  psiri_matrix,
  id.vars = "prey",
  variable.name = "predator",
  value.name = "psiri_legacy"
)

comparison <- merge(
  diet_summary_compare[
    ,
    .(
      prey,
      predator,
      fo,
      count,
      volume,
      PN,
      PV,
      PSIRI
    )
  ],
  legacy_long,
  by = c("prey", "predator")
)

comparison[
  order(predator, -psiri_legacy)
]

comparison[
  ,
  .(
    correlation = cor(
      PSIRI,
      psiri_legacy
    )
  ),
  by = predator
]

###############################################################################
## Validate reconstructed PSIRI
###############################################################################

summary(diet_summary$PSIRI)

diet_summary[order(predator, -PSIRI)]

network_final
str(network_final)
names(network_final)
dim(network_final)


###############################################################################
## Candidate mapping between prey categories and network compartments
###############################################################################

###############################################################################
## Candidate mapping between prey categories and network compartments
###############################################################################

network_map <- data.table(
  
  prey = c(
    
    ## -----------------------------------------------------------------------
    ## r1 - Terrestrial invertebrates
    ## -----------------------------------------------------------------------
    "Acari",
    "Araneae",
    "Blattaria",
    "Chilopoda",
    "Coleoptera",
    "Dermaptera",
    "Diplopoda",
    "Diptera",
    "Formicidae",
    "Gastropoda",
    "Hemiptera",
    "Hymenoptera",
    "Isoptera",
    "Ixodidae",
    "Lepidoptera_larvae",
    "Opiliones",
    "Orthoptera",
    "Phasmida",
    "Unidentified_insects",
    
    ## -----------------------------------------------------------------------
    ## r2 - Aquatic invertebrates
    ## -----------------------------------------------------------------------
    "Amphipoda",
    "Cephalopoda",
    "Decapoda",
    "Isopoda",
    
    ## -----------------------------------------------------------------------
    ## Vertebrates
    ## -----------------------------------------------------------------------
    "Rodentia",
    "Rattus_rattus",
    "Kerodon_rupestris",
    "Birds",
    "Amphisbaena_ridleyi",
    "Trachylepis_atlantica",
    
    ## -----------------------------------------------------------------------
    ## Tentative assignment
    ## -----------------------------------------------------------------------
    "Squamata (temporary)",
    
    ## -----------------------------------------------------------------------
    ## Still unresolved
    ## -----------------------------------------------------------------------
    "Mus_musculus",
    "Mammal",
    "Fish",
    "Animal_tissue",
    "Seeds",
    "Fruits",
    "Vegetal",
    "Unidentified_material"
    
  ),
  
  group = c(
    
    ## r1
    rep("r1", 19),
    
    ## r2
    rep("r2", 4),
    
    ## Vertebrates
    "r3",   # Rodentia
    "m2",   # Rattus_rattus
    "r4",   # Kerodon_rupestris
    "r5",   # Birds
    "r6",   # Amphisbaena_ridleyi
    "m1",   # Trachylepis_atlantica
    
    ## Tentative
    "m1",   # Squamata (temporary)
    
    ## Unknown
    NA,     # Mus_musculus
    NA,     # Mammal
    NA,     # Fish
    NA,     # Animal_tissue
    NA,     # Seeds
    NA,     # Fruits
    NA,     # Vegetal
    NA      # Unidentified_material
    
  )
  
)

network_map[is.na(group)]
setdiff(
  abundance_matrix$prey,
  network_map$prey
)
setdiff(
  network_map$prey,
  abundance_matrix$prey
)

###############################################################################
## Aggregate abundance by network compartment
###############################################################################

abundance_long <- melt(
  abundance_matrix,
  id.vars = "prey",
  variable.name = "predator",
  value.name = "count"
)

abundance_long <- merge(
  abundance_long,
  network_map,
  by = "prey",
  all.x = TRUE
)

candidate_rawFreq <- abundance_long[
  !is.na(group),
  .(
    count = sum(count, na.rm = TRUE)
  ),
  by = .(group, predator)
]

candidate_rawFreq

rawFreq_long <- melt(
  rawFreq,
  id.vars = "Column1",
  variable.name = "predator",
  value.name = "legacy"
)

setnames(rawFreq_long, "Column1", "group")

###############################################################################
## Convert network abbreviations to species names
###############################################################################

code_map <- data.table(
  predator = c(
    "p1",
    "m2",
    "p2",
    "p3"
  ),
  predator_name = c(
    "Felis_catus",
    "Rattus_rattus",
    "Salvator_merianae",
    "Rhinella"
  )
)

rawFreq_long <- merge(
  rawFreq_long,
  code_map,
  by = "predator"
)

rawFreq_long[, predator := predator_name]
rawFreq_long[, predator_name := NULL]

comparison <- merge(
  candidate_rawFreq,
  rawFreq_long,
  by = c("group", "predator"),
  all = TRUE
)

comparison[
  order(group, predator)
]

###############################################################################
## Investigate unmapped prey categories
###############################################################################

unmapped <- abundance_long[
  is.na(group),
  .(
    prey,
    predator,
    count
  )
]

dcast(
  unmapped,
  prey ~ predator,
  value.var = "count",
  fill = 0
)[order(prey)]

unmapped[
  ,
  .(
    total = sum(count, na.rm = TRUE)
  ),
  by = prey
][order(-total)]

unmapped[
  ,
  .(
    total = sum(count, na.rm = TRUE)
  ),
  by = predator
]


###############################################################################
## 13. Reverse-engineer aggregation rules
###############################################################################

evaluate_network <- function(exclude = character()) {
  
  tmp <- copy(abundance_long)
  
  ## Remove candidate taxa
  tmp <- tmp[!prey %in% exclude]
  
  ## Aggregate
  candidate <- tmp[
    !is.na(group),
    .(
      count = sum(count, na.rm = TRUE)
    ),
    by = .(group, predator)
  ]
  
  ## Complete missing combinations
  candidate <- CJ(
    group = unique(rawFreq_long$group),
    predator = unique(rawFreq_long$predator)
  )[
    candidate,
    on = .(group, predator)
  ]
  
  candidate[is.na(count), count := 0]
  
  comparison <- merge(
    candidate,
    rawFreq_long,
    by = c("group","predator")
  )
  
  comparison[
    ,
    difference := count - legacy
  ]
  
  list(
    
    total_error =
      sum(abs(comparison$difference)),
    
    max_error =
      max(abs(comparison$difference)),
    
    comparison = comparison
    
  )
  
}

baseline <- evaluate_network()
baseline$total_error
baseline$comparison[
  order(group, predator)
]

evaluate_network(
  exclude = "Diptera_larvae"
)$total_error

evaluate_network(
  exclude = c(
    "Diptera_larvae",
    "Diptera_pupa"
  )
)$total_error

candidate_taxa <- c(
  "Diptera_larvae",
  "Diptera_pupa",
  "Coleoptera_larvae",
  "Coleoptera_pupa",
  "Unidentified_insects"
)

results <- rbindlist(
  lapply(candidate_taxa, function(x){
    out <- evaluate_network(
      exclude = x
    )
    data.table(
      excluded = x,
      total_error = out$total_error,
      max_error = out$max_error
    )
  })
)

results[
  order(total_error)
]

pairs <- combn(
  candidate_taxa,
  2,
  simplify = FALSE
)

pair_results <- rbindlist(
  
  lapply(pairs, function(x){
    
    out <- evaluate_network(
      exclude = x
    )
    
    data.table(
      
      excluded =
        paste(x, collapse = " + "),
      
      total_error =
        out$total_error,
      
      max_error =
        out$max_error
      
    )
  })
)
pair_results[
  order(total_error)
]

###############################################################################
## Candidate aggregation function
###############################################################################

aggregate_candidate <- function(
    diet,
    merge_diptera_larvae = TRUE,
    merge_diptera_pupae = TRUE,
    merge_coleoptera_larvae = TRUE,
    merge_coleoptera_pupae = TRUE,
    squamata_to_trachylepis = TRUE,
    mus_to_rodentia = FALSE
){
  
  dt <- copy(diet)
  
  ## ----------------------------------------------------------
  ## Diptera
  ## ----------------------------------------------------------
  
  if(merge_diptera_larvae){
    
    dt[
      prey == "Diptera_larvae",
      prey := "Diptera"
    ]
    
  }
  
  if(merge_diptera_pupae){
    
    dt[
      prey == "Diptera_pupa",
      prey := "Diptera"
    ]
    
  }
  
  ## ----------------------------------------------------------
  ## Coleoptera
  ## ----------------------------------------------------------
  
  if(merge_coleoptera_larvae){
    
    dt[
      prey == "Coleoptera_larvae",
      prey := "Coleoptera"
    ]
    
  }
  
  if(merge_coleoptera_pupae){
    
    dt[
      prey == "Coleoptera_pupa",
      prey := "Coleoptera"
    ]
    
  }
  
  ## ----------------------------------------------------------
  ## Squamata
  ## ----------------------------------------------------------
  
  if(squamata_to_trachylepis){
    
    dt[
      prey == "Squamata",
      prey := "Trachylepis_atlantica"
    ]
    
  }else{
    
    dt[
      prey == "Squamata",
      prey := "Squamata (temporary)"
    ]
    
  }
  
  ## ----------------------------------------------------------
  ## Mus musculus
  ## ----------------------------------------------------------
  
  if(mus_to_rodentia){
    
    dt[
      prey == "Mus_musculus",
      prey := "Rodentia"
    ]
    
  }
  
  ## Aggregate
  
  dt <- dt[
    ,
    lapply(.SD, sum, na.rm = TRUE),
    by = prey
  ]
  
  return(dt)
  
}

candidate <- aggregate_candidate(
  diet_collapsed
)

candidate

###############################################################################
## 11. Build broad taxonomic abundance table
###############################################################################

broad_abundance <- merge(
  abundance_matrix,
  broad_group_map,
  by.x = "prey",
  by.y = "psiri_taxon",
  all.x = TRUE
)

stopifnot(!any(is.na(broad_abundance$broad_group)))

broad_abundance <- broad_abundance[
  ,
  lapply(.SD, sum, na.rm = TRUE),
  by = broad_taxon,
  .SDcols = c(
    "Salvator_merianae",
    "Felis_catus",
    "Rattus_rattus",
    "Rhinella"
  )
]

setorder(
  broad_abundance,
  broad_taxon
)

###############################################################################
## 12. Validate broad taxonomic abundance table
###############################################################################

stopifnot(!anyDuplicated(broad_abundance$broad_group))

stopifnot(
  sum(broad_abundance$Salvator_merianae, na.rm = TRUE) ==
    sum(abundance_matrix$Salvator_merianae, na.rm = TRUE)
)

stopifnot(
  sum(broad_abundance$Felis_catus, na.rm = TRUE) ==
    sum(abundance_matrix$Felis_catus, na.rm = TRUE)
)

stopifnot(
  sum(broad_abundance$Rattus_rattus, na.rm = TRUE) ==
    sum(abundance_matrix$Rattus_rattus, na.rm = TRUE)
)

stopifnot(
  sum(broad_abundance$Rhinella, na.rm = TRUE) ==
    sum(abundance_matrix$Rhinella, na.rm = TRUE)
)

###############################################################################
## 14. Construct network abundance table
###############################################################################

network_abundance <- merge(
  broad_abundance,
  network_group_map,
  by = "broad_taxon",
  all.x = TRUE
  
)

network_abundance <- network_abundance[
  !is.na(network_group),
  lapply(.SD, sum, na.rm = TRUE),
  by = network_group,
  .SDcols = c(
    "Salvator_merianae",
    "Felis_catus",
    "Rattus_rattus",
    "Rhinella"
  )
]

setorder(
  network_abundance,
  network_group
)

network_abundance

###############################################################################
## 15. Read legacy network matrix
###############################################################################

rawFreq <- fread("data/rawFreq.csv")

setnames(
  rawFreq,
  "Column1",
  "network_group"
)

rawFreq
 
