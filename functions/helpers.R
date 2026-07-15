### HELPER FUNCTIONS

###############################################################################
## Harmonize prey taxonomy
###############################################################################

harmonize_taxonomy <- function(diet_raw, taxon_map, 
                               finalNaming) {
  ## Check required columns
  stopifnot("prey" %in% names(diet_raw))
  stopifnot(all(c("Original_Name", "Standardized_Name", 
                  "Ecological_Group", "Keep") %in% names(taxon_map)))
  
  ## Rename prey column temporarily
  out <- copy(diet_raw)
  
  ## Join taxonomy
  out <- merge(
    out,
    taxon_map,
    by.x = "prey",
    by.y = "Original_Name",
    all.x = TRUE,
    sort = FALSE
  )
  
  # Remove categories that need exclusion
  out <- out[Keep==TRUE, ]
  
  ## Check for taxa without mapping
  missing <- out[is.na(Standardized_Name), prey]
  
  if(length(missing) > 0){
    stop(
      "The following taxa have no mapping:\n",
      paste(missing, collapse = ", ")
    )
  }

  ## Replace original names
  if (finalNaming == "Ecological"){
    out[, prey := Ecological_Group]
  } else {
    if (finalNaming == "Standardized"){
      out[, prey := Standardized_Name]
    } else {
      stop("finalNaming needs to be provided as either 'Ecological' or 'Standardized'")
    }
    
  }
  ## Remove helper columns
  out[, Standardized_Name := NULL]
  out[, Ecological_Group := NULL]
  out[, Keep := NULL]
  
  return(out)
}

###############################################################################
## Aggregate Preys
###############################################################################

aggregate_raw_diet <- function(dt, predator){
  if (!predator %in% c("cats", "rats", "tegu", "toads"))
    stop("predator needs to be one of: cats, rats, tegu, toads")
  if (predator %in% c("cats", "rats", "tegu")){
    dt <- dt[,.(count = sum(count),
                volume_ml = sum(volume_ml)),
             by = .(sample_id, predator, prey)]
  }
  if (predator == c("toads")){
    dt <- dt[,.(count_N = sum(count_N),
                occurances_NO = sum(occurances_NO),
                volume_ml = sum(volume_ml)),
             by = prey]
  }
  return(dt)
}

###############################################################################
## Calculate IRI 
###############################################################################

calculate_iri <- function(
    diet_data,
    n_stomachs = NULL,
    rawData = TRUE){
  
  if(rawData){
    
    ## Cats, rats and tegu
    summary <- diet_data[
      ,
      .(
        N  = sum(count),
        NO = uniqueN(sample_id),
        total_volume_ml  = sum(volume_ml)
      ),
      by = prey
    ]
    
    n_stomachs <- uniqueN(diet_data$sample_id)
    
  } else {
    
    ## Aggregated datasets (e.g. toads)
    
    if(is.null(n_stomachs))
      stop("When rawData = FALSE, n_stomachs must be provided.")
    
    summary <- diet_data[
      ,
      .(
        N  = sum(count_N),
        NO = sum(occurances_NO),
        total_volume_ml  = sum(volume_ml)
      ),
      by = prey
    ]
    
  }
  
  ## Relative importance metrics
  summary[
    ,
    `:=`(
      percent_N  = 100 * N / sum(N),
      FO         = 100 * NO / n_stomachs,
      percent_V  = 100 * total_volume_ml / sum(total_volume_ml)
    )
  ]
  
  summary[
    ,
    IRI := FO * (percent_N + percent_V)
  ]
  
  summary[
    ,
    percent_IRI := 100 * IRI / sum(IRI)
  ]
  
  setcolorder(
    summary,
    c(
      "prey",
      "N",
      "NO",
      "total_volume_ml",
      "percent_N",
      "FO",
      "percent_V",
      "IRI",
      "percent_IRI"
    )
  )
  return(summary)
}

###############################################################################
# calculate_effect_matrix(): Computes the T matrix of total (direct + 
#                            indirect) effects.
#
# Note: upper triangle of W is sign-flipped to encode antagonistic interactions.
#       Diagonal (self-effects) set to NA before summarising contributions.
#
###############################################################################

calculate_effect_matrix <- function(mat, R) {
  # NOTE: Consumers_table was removed due to changes in the code. 
  # Here just kept as legacy
  
  ## -------------------------------------------------------------------------
  ## Input checks
  ## -------------------------------------------------------------------------
  
  stopifnot(is.matrix(mat))
  stopifnot(nrow(mat) == ncol(mat))
  
  stopifnot(is.numeric(R))
  stopifnot(length(R) == nrow(mat))
  
  stopifnot(!anyNA(mat))
  stopifnot(!anyNA(R))
  
  S <- nrow(mat)
  
  ## -------------------------------------------------------------------------
  ## Normalize interaction strengths
  ##
  ## Original mathematics:
  ## Each row is divided by its total interaction strength.
  ## -------------------------------------------------------------------------
  
  row_totals <- rowSums(mat)
  stopifnot(all(row_totals > 0))
  
  W <- mat / row_totals
  W[which(upper.tri(mat))] <- (-1)*W[which(upper.tri(mat))] #only for antagonisms
  
  
  ## -------------------------------------------------------------------------
  ## Interaction matrix
  ## -------------------------------------------------------------------------
  
  P <- diag(R)
  
  ## -------------------------------------------------------------------------
  ## Total effects matrix
  ##
  ## T = (I - P W)^(-1)
  ## -------------------------------------------------------------------------
  
  I <- diag(S)
  
  cat("Maximum |W|:", max(abs(W)), "\n")
  cat("Maximum |PW|:", max(abs(P %*% W)), "\n")
  
  Tmat <- solve(I - P %*% W)
  rownames(Tmat) <- rownames(mat)
  colnames(Tmat) <- colnames(mat)
  
  ## Ignore self-effects when computing summaries
  
  T_no_diag <- Tmat
  diag(T_no_diag) <- NA
  
  ## -------------------------------------------------------------------------
  ## Direct Effects
  ## -------------------------------------------------------------------------
  
  direct <- P %*% W
  rownames(direct) <- rownames(mat)
  colnames(direct) <- colnames(mat)
  
  ## -------------------------------------------------------------------------
  ## Mean contribution (Tout)
  ## -------------------------------------------------------------------------
  
  mean_contribution <- colMeans(T_no_diag, na.rm = TRUE)
  
  ## -------------------------------------------------------------------------
  ## Binary interaction matrix
  ## -------------------------------------------------------------------------
  
  B <- mat
  B[B>0] <- 1
  
  ## -------------------------------------------------------------------------
  ## Proportion of indirect effects
  ## -------------------------------------------------------------------------
  
  indirect_only <- T_no_diag * (1 - B)
  
  prop_ind_eff <-
    sum(indirect_only, na.rm = TRUE) /
    sum(T_no_diag, na.rm = TRUE)
  
  prop_ind_eff_tin <-
    rowSums(indirect_only, na.rm = TRUE) /
    rowSums(T_no_diag, na.rm = TRUE)
  
  prop_ind_eff_tout <-
    colSums(indirect_only, na.rm = TRUE) /
    colSums(T_no_diag, na.rm = TRUE)
  
  ## -------------------------------------------------------------------------
  ## Return
  ## -------------------------------------------------------------------------
  
  finalList <- list(
    ## Mean total effect exerted by each species
    mean_total_effect = mean_contribution,
    ## Network-wide proportion of indirect effects
    proportion_indirect_effects = prop_ind_eff,
    ## Proportion of effects received by each species that are indirect
    proportion_indirect_effects_received = prop_ind_eff_tin,
    ## Proportion of effects exerted by each species that are indirect
    proportion_indirect_effects_exerted = prop_ind_eff_tout,
    ## Direct-effect matrix (P × W)
    direct_effect_matrix = direct,
    ## Total-effect matrix (direct + indirect)
    total_effect_matrix = Tmat 
  )
  return(finalList)
}

###############################################################################
## Generic scenario function
###############################################################################

simulate_management_scenario <- function(
    effect_matrix,
    response_strength,
    modifications,
    redistribution = c("none", "predator_diet")){
  # "none" = legacy implementation: when a species disappears, it disappears
  # "predator_diet" = redistribute the removed prey pressure among the 
  #                   predator's remaining prey.  
  redistribution <-
    match.arg(redistribution)
  scenario_matrix <- effect_matrix
  
  ## -------------------------------------------------------------------------
  ## Apply all requested modifications
  ##
  ## modifications must contain:
  ##
  ## source     = species exerting the effect (column)
  ## target     = species receiving the effect (row)
  ## reduction  = proportion removed (0–1)
  ##
  ## If target = NA, all outgoing interactions from source are modified.
  ## -------------------------------------------------------------------------
  
  for(i in seq_len(nrow(modifications))){
    
    source <-
      modifications$source[i]
    
    target <-
      modifications$target[i]
    
    reduction <-
      modifications$reduction[i]
    
    if(is.na(target)){
      
      ############################################################
      ## Store original effects BEFORE removing the species
      ############################################################
      
      original_column <-
        scenario_matrix[, source]
      
      ############################################################
      ## Remove (or reduce) the managed species
      ############################################################
      
      scenario_matrix[, source] <-
        original_column *
        (1 - reduction)
      
      ############################################################
      ## Optional predator diet redistribution
      ############################################################
      
      if(redistribution == "predator_diet"){
        
        ## Which predators consume this prey?
        # predators <-
        #   names(original_column)[abs(original_column) > 0]
        predators <-
          colnames(scenario_matrix)[
            abs(original_column) > 0
          ]
        
        ## Don't redistribute to the removed species itself
        predators <-
          setdiff(predators, source)
        
        for(predator in predators){
          
          ## Current diet of this predator
          diet <-
            abs(scenario_matrix[, predator])
          
          ## Amount of diet removed
          removed_amount <-
            abs(original_column[predator]) * reduction
          
          if(removed_amount == 0)
            next
          
          ## Remaining prey
          diet[source] <-
            diet[source] * (1 - reduction)
          
          remaining_prey <-
            names(diet)[diet > 0]
          
          remaining_prey <-
            setdiff(remaining_prey, source)
          
          if(length(remaining_prey) == 0)
            next
          
          ## Redistribute proportionally
          weights <-
            diet[remaining_prey] /
            sum(diet[remaining_prey])
          
          diet[remaining_prey] <-
            diet[remaining_prey] +
            removed_amount * weights
          
          ## Restore original signs
          scenario_matrix[, predator] <-
            diet * sign(scenario_matrix[, predator])
          
        }
        
      }
      
    } else {
      
      ############################################################
      ## Modify only one interaction
      ############################################################
      
      scenario_matrix[target, source] <-
        scenario_matrix[target, source] *
        (1 - reduction)
      
    }
    # OLDER IMPLEMENTATION 
    # if(is.na(target)){
    #   
    #   scenario_matrix[, source] <-
    #     scenario_matrix[, source] *
    #     (1 - reduction)
    #   
    # }else{
    #   
    #   scenario_matrix[target, source] <-
    #     scenario_matrix[target, source] *
    #     (1 - reduction)
    #   }
    }
  
  ## -------------------------------------------------------------------------
  ## Recalculate network
  ## -------------------------------------------------------------------------
  
  scenario_results <-
    calculate_effect_matrix(
      mat = scenario_matrix,
      R   = response_strength
    )
  
  delta_total_effect_matrix <-
    scenario_results$total_effect_matrix -
    baseline_total_effect_matrix
  
  delta_mean_total_effect <-
    scenario_results$mean_total_effect -
    baseline_mean_total_effect
  
  sign_changes <-
    sum(
      sign(baseline_total_effect_matrix) !=
        sign(scenario_results$total_effect_matrix),
      na.rm = TRUE
    )
  
  list(
    
    scenario_matrix = scenario_matrix,
    
    results = scenario_results,
    
    delta_total_effect_matrix =
      delta_total_effect_matrix,
    
    delta_mean_total_effect =
      delta_mean_total_effect,
    
    sign_changes = sign_changes,
    
    winners =
      sort(
        delta_mean_total_effect,
        decreasing = TRUE
      ),
    
    losers =
      sort(
        delta_mean_total_effect
      )
    
  )
  
}