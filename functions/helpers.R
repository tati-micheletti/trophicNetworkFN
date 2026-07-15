### HELPER FUNCTIONS

###############################################################################
## Harmonize prey taxonomy
###############################################################################

### PSIRI
# harmonize_taxonomy <- function(diet_raw, taxon_map) {
#   
#   ## Check required columns
#   stopifnot("prey" %in% names(diet_raw))
#   stopifnot(all(c("original_taxon", "standardized_taxon") %in% names(taxon_map)))
#   
#   ## Rename prey column temporarily
#   out <- copy(diet_raw)
#   
#   ## Join taxonomy
#   out <- merge(
#     out,
#     taxon_map,
#     by.x = "prey",
#     by.y = "original_taxon",
#     all.x = TRUE,
#     sort = FALSE
#   )
#   
#   ## Check for taxa without mapping
#   missing <- out[is.na(standardized_taxon), prey]
#   
#   if(length(missing) > 0){
#     stop(
#       "The following taxa have no mapping:\n",
#       paste(missing, collapse = ", ")
#     )
#   }
#   ## Replace original names
#   out[, prey := standardized_taxon]
#   ## Remove helper column
#   out[, standardized_taxon := NULL]
#   return(out)
# } 

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
  
  ## Check for taxa without mapping
  missing <- out[is.na(Standardized_Name), prey]
  
  if(length(missing) > 0){
    stop(
      "The following taxa have no mapping:\n",
      paste(missing, collapse = ", ")
    )
  }
  out <- out[Keep==TRUE,]
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
## Collapse duplicated taxa
###############################################################################

collapse_taxa <- function(diet_standardized) {
  
  stopifnot("prey" %in% names(diet_standardized))
  
  numeric_cols <- setdiff(names(diet_standardized), "prey")
  
  collapsed <- diet_standardized[
    ,
    lapply(.SD, function(x) sum(x, na.rm = TRUE)),
    by = prey,
    .SDcols = numeric_cols
  ]
  
  ## Replace zeros with NA to match the original workbook
  collapsed[
    ,
    (numeric_cols) := lapply(.SD, function(x) fifelse(x == 0, NA_real_, x)),
    .SDcols = numeric_cols
  ]
  
  setorder(collapsed, prey)
  
  return(collapsed)
}

###############################################################################
## Standardize taxa
###############################################################################

standardize_taxa <- function(diet_raw, taxon_map){
  diet_standardized <- harmonize_taxonomy(diet_raw,taxon_map)
  collapse_taxa(diet_standardized)
}

###############################################################################
# indirect(): Computes the T matrix of total (direct + indirect) effects.
#
# Arguments:
#   mat  - Square interaction matrix (unipartite) or biadjacency matrix
#          (bipartite, not yet implemented). Entry [i,j] = effect of j on i.
#   R    - Numeric vector of interaction-dependency parameters (one per species).
#          Length = nrow(mat) for unipartite; nrow+ncol for bipartite.
#   type - "unipartite" (default) or "bipartite".
#
# Returns a list:
#   mean_contribution  - Mean total effect exerted by each species (col means of T)
#   prop_ind_eff       - Network-wide proportion of effects that are indirect
#   prop_ind_eff_tin   - Per-species proportion of received effects that are indirect
#   prop_ind_eff_tout  - Per-species proportion of exerted effects that are indirect
#   direct             - Direct-effect matrix (P %*% W)
#   Tmat               - Full S x S total-effects matrix T = (I - PW)^{-1}
#
# Note: upper triangle of W is sign-flipped to encode antagonistic interactions.
#       Diagonal (self-effects) set to NA before summarising contributions.
#
# Computes indirect effects from a matrix (mat) and a vector of interaction dependencies (R)
# R is a single vector with length = sum(nrow(mat), ncol(mat)) corresponding to the R value of the ROWS and COLUMNS
# Example: 
# indirect(mat=matrix(rnbinom(Na*Np,1,0.2),nrow=Na,ncol=Np),R=0.5)
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
## Prepare PSIRI
###############################################################################

aggregate_psiri_taxa <- function(diet_collapsed, psiri_map) {
  
  ## check columns
  stopifnot("prey" %in% names(diet_collapsed))
  
  stopifnot(
    all(
      c("standardized_taxon",
        "psiri_taxon") %in% names(psiri_map)
    )
  )
  
  out <- copy(diet_collapsed)
  
  out <- merge(
    out,
    psiri_map,
    by.x = "prey",
    by.y = "standardized_taxon",
    all.x = TRUE,
    sort = FALSE
  )
  
  missing <- sort(unique(out[is.na(psiri_taxon), prey]))
  
  if (length(missing) > 0) {
    stop(
      sprintf(
        "Missing PSIRI mappings for %d taxa:\n%s",
        length(missing),
        paste(missing, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  
  out[, prey := psiri_taxon]
  out[, psiri_taxon := NULL]
  
  collapse_taxa(out)
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
## Read raw data for psiri -- work from Ju, Felipe and further, the 
## new toad additions
###############################################################################

## For Tegu, Cats, Rats
read_psiri_block <- function(
    sheet,
    start_col,
    predator,
    prey_col   = 1,
    fo_col     = 2,
    count_col  = 4,
    volume_col = 6
){
  cols <- start_col + c(
    prey_col,
    fo_col,
    count_col,
    volume_col
  ) - 1
  
  out <- copy(
    sheet[
      ,
      .(
        prey   = trimws(as.character(.SD[[1]])),
        fo     = as.numeric(.SD[[2]]),
        count  = as.numeric(.SD[[3]]),
        volume = as.numeric(.SD[[4]])
      ),
      .SDcols = cols
    ]
  )
  
  out[, predator := predator]
  
  ## Remove empty rows
  out <- out[
    !is.na(prey) &
      prey != ""
  ]
  ## Remove title rows
  out <- out[
    !grepl("\\(n =", prey)
  ]
  ## Remove footer rows
  out <- out[
    !grepl("^Fonte|^%PSIRI", prey)
  ]
  ## Arrange columns
  setcolorder(
    out,
    c(
      "prey",
      "predator",
      "fo",
      "count",
      "volume"
    )
  )
  out
}

## For Toads
read_psiri_summary <- function(
    sheet,
    predator
){
  
  out <- copy(
    sheet[
      ,
      .(
        prey  = trimws(as.character(ITENS)),
        fo    = as.numeric(`FO%`),
        pn    = as.numeric(`PN %`),
        pv    = as.numeric(`PV%`),
        psiri = as.numeric(`/100`)
      )
    ]
  )
  
  out[, predator := predator]
  
  ## Remove empty rows
  out <- out[
    !is.na(prey) &
      prey != ""
  ]
  
  ## Arrange columns
  setcolorder(
    out,
    c(
      "prey",
      "predator",
      "fo",
      "pn",
      "pv",
      "psiri"
    )
  )
  
  out
}

combine_diet_summaries <- function(
    vertebrate_summary,
    published_summary
){
  
  vertebrate_summary <- copy(vertebrate_summary)
  
  ## Calculate PN and PV
  vertebrate_summary[
    ,
    pn := 100 * count / sum(count),
    by = predator
  ]
  
  vertebrate_summary[
    ,
    pv := 100 * volume / sum(volume),
    by = predator
  ]
  
  vertebrate_summary[
    is.nan(pn),
    pn := 0
  ]
  
  vertebrate_summary[
    is.nan(pv),
    pv := 0
  ]
  
  ## PSIRI will be calculated later
  vertebrate_summary[
    ,
    psiri := NA_real_
  ]
  
  setcolorder(
    vertebrate_summary,
    c(
      "prey",
      "predator",
      "fo",
      "count",
      "volume",
      "pn",
      "pv",
      "psiri"
    )
  )
  
  published_summary <- copy(published_summary)
  
  published_summary[
    ,
    `:=`(
      count = NA_real_,
      volume = NA_real_
    )
  ]
  
  setcolorder(
    published_summary,
    c(
      "prey",
      "predator",
      "fo",
      "count",
      "volume",
      "pn",
      "pv",
      "psiri"
    )
  )
  
  rbindlist(
    list(
      vertebrate_summary,
      published_summary
    ),
    use.names = TRUE,
    fill = TRUE
  )
  
}

###############################################################################
#' Extract dietary observations from original spreadsheets
#'
#' Converts the original wide-format dietary spreadsheets (abundance and
#' volume) into a canonical long-format table suitable for downstream
#' food-web analyses.
#'
#' Each row of the returned table represents a single prey item recorded in
#' a single sample, containing the sample identifier, the original prey
#' description, the observed number of individuals and the corresponding
#' mean prey volume.
#'
#' Only prey items with abundance > 0 are retained. Missing volume values
#' are replaced with zero.
#'
#' @param workbook Character. Path to the original Excel workbook.
#' @param abundance_sheet Character. Name of the worksheet containing prey
#'   abundance data.
#' @param volume_sheet Character. Name of the worksheet containing prey
#'   volume data.
#'
#' @return A data.table with the columns:
#' \describe{
#'   \item{sample_id}{Unique sample identifier.}
#'   \item{prey_original}{Original prey description exactly as recorded in the source data.}
#'   \item{count}{Number of prey individuals observed in the sample.}
#'   \item{mean_volume_ml}{Mean volume (mL) of the prey item.}
#' }
#'
#' @author Tatiane Micheletti
#'
###############################################################################

extract_diet_data <- function(workbook,
                              abundance_sheet,
                              volume_sheet,
                              predator){
      abundance <- as.data.table(
      read_excel(
        workbook,
        sheet = abundance_sheet
      )
    )
    
    volume <- as.data.table(
      read_excel(
        workbook,
        sheet = volume_sheet
      )
    )
    # Had to manually make sure the names matched
    stopifnot(
      identical(
        names(abundance),
        names(volume)
      )
    )
    
    sample_ids <- abundance[[1]]
    prey_names <- names(abundance)[-1]
    
    diet_data <- rbindlist(
      lapply(
        seq_along(prey_names),
        function(i){
          data.table(
            sample_id = sample_ids,
            prey = prey_names[i],
            count = abundance[[i + 1]],
            volume_ml = volume[[i + 1]]
          )
        }
      )
    )
    
    diet_data <- diet_data[
      !is.na(count) &
        count > 0
    ]
    
    diet_data[
      is.na(volume_ml),
      volume_ml := 0
    ]

    diet_data[, predator := predator]
    setcolorder(
      diet_data,
      c(
        "sample_id",
        "predator",
        "prey",
        "count",
        "volume_ml"
      )
    )
    return(diet_data)
  }

extract_new_tegu_data <- function(workbook){
  
  sheets <- grep(
    "^SACO",
    excel_sheets(workbook),
    value = TRUE
  )
  
  tegu_raw <- rbindlist(
    lapply(
      sheets,
      function(sheet){
        
        x <- as.data.table(
          read_excel(
            workbook,
            sheet = sheet,
            col_names = FALSE
          )
        )
        
        ## Keep only stomach-content table
        x <- x[, 1:9]
        
        setnames(
          x,
          c(
            "sample_id",
            "identification",
            "prey",
            "count",
            "individual",
            "length_mm",
            "width_mm",
            "volume_ml",
            "total_volume_ml"
          )
        )
        
        ## Remove title rows
        x <- x[-c(1, 2)]
        
        ## Fill merged cells
        fill_cols <- c(
          "sample_id",
          "identification",
          "prey",
          "count"
        )
        for(col in fill_cols){
          last <- NA
          for(i in seq_len(nrow(x))){
            if(!is.na(x[[col]][i])){
              last <- x[[col]][i]
            }else{
              x[[col]][i] <- last
            }
          }
        }
        ## Safe numeric conversion
        clean_numeric <- function(z){
          z <- trimws(as.character(z))
          z[z %in% c("", "-", "NA")] <- NA
          as.numeric(z)
        }
        
        x[
          ,
          `:=`(
            count = clean_numeric(count),
            individual = clean_numeric(individual),
            volume_ml = clean_numeric(volume_ml)
          )
        ]
        
        ## Keep only measured prey individuals
        x <- x[
          !is.na(individual) &
            !is.na(volume_ml)
        ]
        
        ## Number measured should never exceed number counted
        stopifnot(all(x$n_measured <= x$count))
        
        ## Estimate mean prey volume from measured individuals
        x <- x[
          ,
          .(
            count = first(count),
            mean_volume_ml = mean(volume_ml)
          ),
          by = .(
            sample_id,
            prey
          )
        ]
        
        x[, predator := "Salvator_merianae"]
        
        setcolorder(
          x,
          c(
            "sample_id",
            "predator",
            "prey",
            "count",
            "mean_volume_ml"
          )
        )
        
        x
        
      }
    ),
    use.names = TRUE
  )
  
  tegu_raw[]
}

###############################################################################
## Calculate PSIRI 
###############################################################################

calculate_psiri <- function(
    diet_data,
    n_stomachs = NULL,
    rawData = TRUE){
  if(!rawData){
    if (is.null(n_stomachs)) 
      stop("When rawData == FALSE, n_stomachs must be provided")
    ## Calculate totals from individual data
    # i.e., toads
    ## Aggregated datasets (e.g. toads) lack stomach-level composition,
    ## therefore prey-specific PN and PV cannot be reconstructed.
    ## Aggregated percentages are used instead.
    summary <- diet_data[,`:=`(
        PN = 100 * count_N / sum(count_N),
        PV = 100 * volume_ml / sum(volume_ml),
        FO = 100 * occurances_NO / n_stomachs
      )
    ]
    summary[,PSIRI := FO * (PN + PV) / 2]
    } else {
      ## Calculate totals from stomach-level data
      # i.e., cats, tegu, rats
      summary <- diet_data[
        ,
        `:=`(
          PN_stomach = 100 * count / sum(count),
          PV_stomach = 100 * volume_ml / sum(volume_ml)
        ),
        by = sample_id
      ]
      summary <- diet_data[
        ,
        .(
          N  = sum(count),
          NO = uniqueN(sample_id),
          V  = sum(volume_ml),
          PN = mean(PN_stomach),
          PV = mean(PV_stomach)
        ),
        by = prey
      ]
      
      summary[
        ,
        `:=`(
          FO = 100 * NO / uniqueN(diet_data$sample_id),
          PSIRI = (FO / 100) * (PN + PV) / 2
        )
      ]
    }
  # summary[
  #   ,
  #   `:=`(
  #     PN = 100 * N / sum(N),
  #     FO = 100 * NO / n_stomachs,
  #     PV = 100 * V / sum(V)
  #   )
  # ]
  # summary[,PSIRI := (FO/100) * (PN + PV) / 2]
  return(summary)
}

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
