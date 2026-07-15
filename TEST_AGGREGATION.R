
# DISCOVERED THAT THE AGGREGATION WAS DONE AT A BROADER SCALE. THIS SCRIPT IS USELESS NOW.
# Use only runme.R. Exploration.R is also useful to understand what was done, but a bit messy!

###############################################################################
## STEP 1. Candidate aggregation
###############################################################################
if (FALSE){
  candidate <- copy(diet_collapsed)
  
  ## ---- Diptera ----
  
  candidate[
    prey == "Diptera_larvae",
    prey := "Diptera"
  ]
  
  candidate[
    prey == "Diptera_pupa",
    prey := "Diptera"
  ]
  
  ## ---- Coleoptera ----
  
  candidate[
    prey == "Coleoptera_larvae",
    prey := "Coleoptera"
  ]
  
  candidate[
    prey == "Coleoptera_pupa",
    prey := "Coleoptera"
  ]
  
  ## ---- Decapoda ----
  
  candidate[
    prey == "Brachyura",
    prey := "Decapoda"
  ]
  
  ## ---- Squamata ----
  
  candidate[
    prey == "Squamata",
    prey := "Trachylepis_atlantica"
  ]
  
  ## Aggregate
  
  candidate <- candidate[
    ,
    lapply(.SD, sum, na.rm = TRUE),
    by = prey
  ]
  
  setnames(
    
    candidate,
    
    c(
      "Teiu",
      "Gato",
      "Rato",
      "Sapo"
    ),
    
    c(
      "Salvator_merianae",
      "Felis_catus",
      "Rattus_rattus",
      "Rhinella"
    )
    
  )
  
  candidate_long <- melt(
    
    candidate,
    
    id.vars = "prey",
    
    variable.name = "predator",
    
    value.name = "count"
    
  )
  
  candidate_long <- merge(
    
    candidate_long,
    
    network_map,
    
    by = "prey",
    
    all.x = TRUE
    
  )
  
  candidate_rawFreq <- candidate_long[
    
    !is.na(group),
    
    .(
      count = sum(count, na.rm = TRUE)
    ),
    
    by = .(
      group,
      predator
    )
    
  ]
  comparison <- merge(
    
    candidate_rawFreq,
    
    rawFreq_long,
    
    by = c(
      "group",
      "predator"
    ),
    
    all = TRUE
    
  )
  
  comparison[is.na(count), count := 0]
  
  comparison[is.na(legacy), legacy := 0]
  
  comparison[
    ,
    difference := count - legacy
  ]
  
  comparison[
    order(
      group,
      predator
    )
  ]
  
  sum(abs(comparison$difference))
  
  
  ###############################################################################
  ## Apply candidate aggregation rules
  ###############################################################################
  
  build_candidate <- function(
    merge_diptera_larvae = TRUE,
    merge_diptera_pupae = TRUE,
    merge_coleoptera_larvae = TRUE,
    merge_coleoptera_pupae = TRUE,
    merge_brachyura = TRUE,
    squamata_target = c(
      "Trachylepis_atlantica",
      "Squamata (temporary)"
    )
  ){
    
    squamata_target <- match.arg(squamata_target)
    
    candidate <- copy(diet_collapsed)
    
    ## Diptera
    if(merge_diptera_larvae)
      candidate[prey == "Diptera_larvae", prey := "Diptera"]
    
    if(merge_diptera_pupae)
      candidate[prey == "Diptera_pupa", prey := "Diptera"]
    
    ## Coleoptera
    if(merge_coleoptera_larvae)
      candidate[prey == "Coleoptera_larvae", prey := "Coleoptera"]
    
    if(merge_coleoptera_pupae)
      candidate[prey == "Coleoptera_pupa", prey := "Coleoptera"]
    
    ## Crustacea
    if(merge_brachyura)
      candidate[prey == "Brachyura", prey := "Decapoda"]
    
    ## Squamata
    candidate[
      prey == "Squamata",
      prey := squamata_target
    ]
    
    candidate <- candidate[
      ,
      lapply(.SD, sum, na.rm = TRUE),
      by = prey
    ]
    
    candidate
  }
  
  ###############################################################################
  ## Evaluate one aggregation hypothesis
  ###############################################################################
  
  evaluate_candidate <- function(candidate){
    
    candidate <- copy(candidate)
    
    setnames(
      candidate,
      c("Teiu","Gato","Rato","Sapo"),
      c(
        "Salvator_merianae",
        "Felis_catus",
        "Rattus_rattus",
        "Rhinella"
      )
    )
    
    candidate_long <- melt(
      candidate,
      id.vars = "prey",
      variable.name = "predator",
      value.name = "count"
    )
    
    candidate_long <- merge(
      candidate_long,
      network_map,
      by = "prey",
      all.x = TRUE
    )
    
    candidate_rawFreq <- candidate_long[
      !is.na(group),
      .(
        count = sum(count, na.rm = TRUE)
      ),
      by = .(
        group,
        predator
      )
    ]
    
    comparison <- merge(
      candidate_rawFreq,
      rawFreq_long,
      by = c("group","predator"),
      all = TRUE
    )
    
    comparison[is.na(count), count := 0]
    comparison[is.na(legacy), legacy := 0]
    
    comparison[
      ,
      difference := count - legacy
    ]
    
    list(
      comparison = comparison,
      total_error = sum(abs(comparison$difference))
    )
    
  }
  
  test <- evaluate_candidate(
    build_candidate()
  )
  
  test$total_error
  
  ###############################################################################
  ## Search all aggregation combinations
  ###############################################################################
  
  settings <- expand.grid(
    
    merge_diptera_larvae    = c(TRUE, FALSE),
    merge_diptera_pupae     = c(TRUE, FALSE),
    
    merge_coleoptera_larvae = c(TRUE, FALSE),
    merge_coleoptera_pupae  = c(TRUE, FALSE),
    
    merge_brachyura         = c(TRUE, FALSE),
    
    squamata_target = c(
      "Trachylepis_atlantica",
      "Squamata (temporary)"
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  results <- rbindlist(
    
    lapply(seq_len(nrow(settings)), function(i){
      
      out <- evaluate_candidate(
        
        build_candidate(
          
          merge_diptera_larvae =
            settings$merge_diptera_larvae[[i]],
          
          merge_diptera_pupae =
            settings$merge_diptera_pupae[[i]],
          
          merge_coleoptera_larvae =
            settings$merge_coleoptera_larvae[[i]],
          
          merge_coleoptera_pupae =
            settings$merge_coleoptera_pupae[[i]],
          
          merge_brachyura =
            settings$merge_brachyura[[i]],
          
          squamata_target =
            settings$squamata_target[[i]]
          
        )
        
      )
      
      data.table(
        
        merge_diptera_larvae =
          settings$merge_diptera_larvae[[i]],
        
        merge_diptera_pupae =
          settings$merge_diptera_pupae[[i]],
        
        merge_coleoptera_larvae =
          settings$merge_coleoptera_larvae[[i]],
        
        merge_coleoptera_pupae =
          settings$merge_coleoptera_pupae[[i]],
        
        merge_brachyura =
          settings$merge_brachyura[[i]],
        
        squamata_target =
          settings$squamata_target[[i]],
        
        total_error =
          out$total_error
        
      )
      
    })
    
  )
  
  results[
    order(total_error)
  ]
  
  
  best <- evaluate_candidate(
    
    build_candidate(
      
      merge_diptera_larvae    = FALSE,
      merge_diptera_pupae     = FALSE,
      
      merge_coleoptera_larvae = FALSE,
      merge_coleoptera_pupae  = FALSE,
      
      merge_brachyura         = TRUE,
      
      squamata_target = "Trachylepis_atlantica"
      
    )
    
  )
  
  best$comparison[
    order(-abs(difference))
  ]
  
  candidate <- build_candidate(
    
    merge_diptera_larvae = FALSE,
    merge_diptera_pupae = FALSE,
    merge_coleoptera_larvae = FALSE,
    merge_coleoptera_pupae = FALSE,
    merge_brachyura = TRUE,
    squamata_target = "Trachylepis_atlantica"
    
  )
  
  candidate[
    ,
    .(
      prey,
      Sapo
    )
  ][order(-Sapo)]
  
  network_map[group == "r1"][order(prey)]
  network_map[group == "r2"][order(prey)]
  setdiff(candidate$prey, network_map$prey)
  
  network_map[
    prey == "Gastropoda",
    group := "r2"
  ]
  
  rhinella_raw <- raw_psiri[
    predator == "Sapo",
    .(
      prey,
      raw_count = count
    )
  ]
  
  rhinella_reconstructed <- abundance_matrix[
    ,
    .(
      prey,
      reconstructed = Rhinella
    )
  ]
  
  merge(
    rhinella_raw,
    rhinella_reconstructed,
    by = "prey",
    all = TRUE
  )
  
  hexapoda <- candidate[
    prey %in% c(
      "Blattaria",
      "Coleoptera",
      "Dermaptera",
      "Diptera",
      "Formicidae",
      "Hemiptera",
      "Hymenoptera",
      "Isoptera",
      "Lepidoptera_larvae",
      "Orthoptera",
      "Phasmida"
    ),
    lapply(.SD, sum),
    .SDcols = c("Teiu", "Gato", "Rato", "Sapo")
  ]
  
  candidate[
    prey %in% c(
      "Acari",
      "Araneae",
      "Blattaria",
      "Chilopoda",
      "Coleoptera",
      "Coleoptera_larvae",
      "Coleoptera_pupa",
      "Dermaptera",
      "Diplopoda",
      "Diptera",
      "Diptera_larvae",
      "Diptera_pupa",
      "Formicidae",
      "Hemiptera",
      "Hymenoptera",
      "Isoptera",
      "Ixodidae",
      "Lepidoptera_larvae",
      "Opiliones",
      "Orthoptera",
      "Phasmida",
      "Unidentified_insects"
    ),
    .(prey, Sapo)
  ][order(-Sapo)]
  
}
