###############################################################################
## PART I. Taxonomic reconstruction
###############################################################################

###############################################################################
## 1. Load packages
###############################################################################

if(!require("Require")){
  install.packages("Require")
}
library("Require")
Require::Require("readxl")
Require::Require("data.table")
Require::Require("igraph")
Require::Require("ggraph")
Require::Require("ggplot2")
Require::Require("reshape2")

###############################################################################
## 2. Read input workbook and create output dir
###############################################################################

input_file <- "data/foodwebData.xlsx"
sheet_names <- excel_sheets(input_file)

sheets <- lapply(sheet_names, function(x)
  as.data.table(read_excel(input_file, sheet = x))
)
names(sheets) <- sheet_names

mapping <- copy(sheets$mapping)
species_parameters <- copy(sheets$species_parameters)
weighting <- copy(sheets$weighting)
cats <- copy(sheets$cats)
rats <- copy(sheets$rats)
tegu <- copy(sheets$tegu)
toads <- copy(sheets$toads)
toads[, volume_ml := mean_volume_mL_V * count_N] # Toads needs one more calculation

dir.create("outputs", showWarnings = FALSE)

###############################################################################
## 3. Source Functions
###############################################################################

source("functions/helpers.R")

###############################################################################
## 4. Harmonize prey names
## Note: The last column in the Workbook FoodwebData called `mapping` (column 
## Ecological_Group) is currently the one used for the final groupping of diet 
## items. 
## TO REMOVE AN ITEM: Simply change in the column `Keep` from TRUE to FALSE, 
##                    save the notebook and re-run the analysis from the 
##                    beginning
## TO CHANGE THE GROUPPING: I strongly suggest creating a new column and using 
##                          the columns name below as argument for `finalNaming`
## /!\ ATTENTION /!\ Several names in different groups match! If you need to 
## exclude  an item for one group, make sure that the same item is NOT in other
## groups. One "easy and clever" way to do this is to slightly change the name 
## of the item in the plans you do NOT want to change and add them as a new row
## on the table, making sure the columns Keep is set to TRUE.
###############################################################################

cats  <- harmonize_taxonomy(cats,  mapping, finalNaming = "Ecological")
rats  <- harmonize_taxonomy(rats,  mapping, finalNaming = "Ecological")
tegu  <- harmonize_taxonomy(tegu,  mapping, finalNaming = "Ecological")
toads <- harmonize_taxonomy(toads, mapping, finalNaming = "Ecological")

###############################################################################
# 5. Aggregate prey types (done PER SAMPLE when possible)
###############################################################################

cats  <- aggregate_raw_diet(cats, predator = "cats")
rats  <- aggregate_raw_diet(rats, predator = "rats")
tegu  <- aggregate_raw_diet(tegu, predator = "tegu")
toads <- aggregate_raw_diet(toads, predator = "toads")

###############################################################################
## 6. Calculate IRI for all species
###############################################################################

toads_iri <- calculate_iri(diet_data = toads, 
                           n_stomachs = 137, 
                           rawData = FALSE)

cats_iri <- calculate_iri(diet_data = cats, 
                          rawData = TRUE)

rats_iri <- calculate_iri(diet_data = rats, 
                          rawData = TRUE)

tegu_iri <- calculate_iri(diet_data = tegu, 
                          rawData = TRUE)

###############################################################################
# 7. Step Combine predator diets
###############################################################################

cats_iri[, predator := "Felis_catus"]
rats_iri[, predator := "Rattus_spp"]
tegu_iri[, predator := "Salvator_merianae"]
toads_iri[, predator := "Rhinella"]

diet_network <- rbindlist(
  list(
    cats_iri,
    rats_iri,
    tegu_iri,
    toads_iri
  ),
  use.names = TRUE
)

###############################################################################
## PART III. Final interaction table
###############################################################################

diet_network <- diet_network[,.(predator, 
                                prey, 
                                interaction_strength = percent_IRI)]

setorder(diet_network, predator, -interaction_strength)

###############################################################################
## 8. Calculate density-weighted predator consumption
###############################################################################

# IRI → diet composition (used for prey → predator effects).
# Density-weighted N → predator pressure on prey (used for predator → prey effects).

## Number of stomachs analysed
samples <- c(
  Felis_catus       = 78, # cats (n=78) Gaiotto et al., 2020
  Salvator_merianae = 22, # tegu (n=22) Gaiotto et al., 2020
  Rattus_spp        = 10, # rats (n=10) Gaiotto et al., 2020
  Rhinella          = 137 # Number of individual samples collected for this study
)

## Predator densities (individuals / ha)
densities <- c(
  Felis_catus       = 0.71, # feral cats, 0.71 ind/ha , Dias et al. 2017
  Salvator_merianae = 3.98, # tegu 3.98 ind/ha Abrahão et al. 2019
  Rattus_spp        = 37, # rats 37 ind/ha, Russell et al. 2018
  Rhinella          = 10.35 # Toads, extrapolation from Solomon Islands, Pikacha et al. 2015
)

predator_consumption <- rbindlist(
  list(
    cats_iri[, .(
      predator = "Felis_catus",
      prey,
      N
    )],
    rats_iri[, .(
      predator = "Rattus_spp",
      prey,
      N
    )],
    tegu_iri[, .(
      predator = "Salvator_merianae",
      prey,
      N
    )],
    toads_iri[, .(
      predator = "Rhinella",
      prey,
      N
    )]
  ))
  
## Build prey × predator consumption matrix
predator_consumption <- as.data.table(dcast(
  predator_consumption,
  prey ~ predator,
  value.var = "N",
  fill = 0
))

## Convert counts to per-predator consumption
predator_consumption[, names(samples) := Map(`/`, .SD, samples),
                     .SDcols = names(samples)]

## Weight by predator density
predator_consumption[,names(densities) := Map(`*`,.SD,densities), 
                     .SDcols = names(densities)]

## Keep prey as the first column
setcolorder(
  predator_consumption,
  c("prey", names(samples)))

###############################################################################
## 8. Build predator-prey interaction matrix
###############################################################################

## All network nodes
network_species <- sort(unique(c(
  predator_consumption$prey,
  predator_consumption |> names() |> setdiff("prey")
)))

## Predator -> prey interaction matrix
predator_prey_matrix <- matrix(
  0,
  nrow = length(network_species),
  ncol = length(network_species),
  dimnames = list(
    network_species,
    network_species
  )
)

## Fill predator -> prey interactions
for(i in seq_len(nrow(predator_consumption))){
  prey <- predator_consumption$prey[i]
  for(predator in names(samples)){
    predator_prey_matrix[prey, predator] <- predator_consumption[i,get(predator)]
  }
}

stopifnot(all(predator_consumption[,prey] %in% rownames(predator_prey_matrix)))
stopifnot(names(samples) %in% colnames(predator_prey_matrix))

###############################################################################
## 9.  Build normalized resource dependency matrix
###############################################################################

resource_dependency_matrix <- matrix(
  0,
  nrow = length(network_species),
  ncol = length(network_species),
  dimnames = list(
    network_species,
    network_species
  )
)

## Calculate relative diet composition
diet_network[
  ,
  dependency :=
    interaction_strength / sum(interaction_strength),
  by = predator
]
## Fill matrix
## Rows = consumers (predators)
## Columns = resources (prey)
## This reproduces the orientation of t(psiri) in the original model.
for(i in seq_len(nrow(diet_network))){
  resource_dependency_matrix[
    diet_network$predator[i],
    diet_network$prey[i]
  ] <- diet_network$dependency[i]
  
}

###############################################################################
## 10. Read species parameters
###############################################################################

## Ensure every network species has parameters
missing_species <- setdiff(
  network_species,
  species_parameters$species
)

if(length(missing_species) > 0){
  stop(
    "Missing parameters for: ",
    paste(missing_species, collapse = ", ")
  )
  
}

## Order parameters to match the network
setkey(species_parameters, species)
species_parameters <- species_parameters[network_species]

## Response strength vector
response_strength <-
  species_parameters$response_strength

names(response_strength) <-
  species_parameters$species

###############################################################################
## 11. Prepare food-web model inputs
###############################################################################

foodweb_model <- list(
  ## Species in network order
  species = network_species,
  ## Top-down interactions
  predator_prey_matrix = predator_prey_matrix,
  ## Bottom-up interactions
  resource_dependency_matrix = resource_dependency_matrix,
  ## Species-specific response strengths
  response_strength = response_strength,
  ## Predator densities
  predator_densities = densities,
  ## Number of diet samples
  predator_samples = samples
  )

###############################################################################
## 14. Visualize the food web
###############################################################################
#TODO: Save the plot!
nodes <- data.table(
  name = network_species
)

nodes <- merge(
  nodes,
  species_parameters[, .(
    species,
    label,
    trophic_group
  )],
  by.x = "name",
  by.y = "species",
  all.x = TRUE
)

edges <- copy(diet_network)
setnames(
  edges,
  c("prey", "predator"),
  c("from", "to")
)

foodweb_graph <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = TRUE
)

foodweb_plot <- ggraph(
  foodweb_graph,
  layout = "stress"
) +
  
  geom_edge_link(
    aes(width = interaction_strength),
    alpha = 0.7,
    colour = "grey40",
    end_cap = circle(5, "mm"),
    start_cap = circle(5, "mm"),
    arrow = arrow(
      length = unit(3.5, "mm"),
      type = "closed"
    )
  ) +
  
  geom_node_point(
    aes(fill = trophic_group),
    shape = 21,
    size = 6,
    colour = "black",
    stroke = 0.4
  ) +
  
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = 3.8,
    point.padding = unit(0.5, "lines"),
    box.padding = unit(0.7, "lines")
  ) +
  
  scale_edge_width(range = c(0.2, 2.8)) +
  
  theme_graph(base_family = "sans")
  
###############################################################################
## 15. Build the direct interaction matrices
###############################################################################

## Top-down effects
predation_matrix <- copy(predator_prey_matrix)
## Bottom-up effects
resource_dep_matrix <- copy(resource_dependency_matrix)
## Species response
response_strength <- response_strength

stopifnot(
  identical(
    rownames(predation_matrix),
    rownames(resource_dep_matrix)
  )
)

stopifnot(
  identical(
    colnames(predation_matrix),
    colnames(resource_dep_matrix)
  )
)

stopifnot(
  identical(
    rownames(predation_matrix),
    names(response_strength)
  )
)

###############################################################################
## 16. Build direct-effect matrix
###############################################################################

## Top-down direct effects (Analogous to Bw)
predator_effect_matrix <- copy(predator_prey_matrix) # Analogous to Bw

# Because psiri had the order in the columns as prey-first, predators at the end,
# we reset the order of resource_dependency_matrix

network_order <- c(
  "Amphisbaena_ridleyi",
  "Aquatic_invertebrates",
  "Birds",
  "Fish",
  "Kerodon_rupestris",
  "Mammals",
  "Mus_musculus",
  "Plants",
  "Reptiles",
  "Terrestrial_invertebrates",
  "Trachylepis_atlantica",
  "Rhinella",
  "Rattus_spp",
  "Salvator_merianae",
  "Felis_catus"
)

stopifnot(setequal(network_order, rownames(predator_prey_matrix)))
stopifnot(setequal(network_order, rownames(resource_dependency_matrix)))

predator_effect_matrix <-
  predator_effect_matrix[
    network_order,
    network_order
  ]

resource_dependency_matrix <-
  resource_dependency_matrix[
    network_order,
    network_order
  ]
response_strength <- response_strength[network_order]

## Bottom-up direct effects (Analogous to M = 0.1*t(psiri), as psiri was in %)
# IMPORTANT: resource_dependency_matrix is ALREADY transposed. This means 
# that t(psiri) ~ resource_dependency_matrix/100.  Therefore, we multiply by 100,
# before multiplying by 0.1 and do not transpose.
resource_effect_matrix <- 0.1 * 100 * resource_dependency_matrix

###############################################################################
## Eventually, we can potentially assign weights for the effects. 
#  Originally, there were no weighting of effects
###############################################################################

# ## Relative importance of each process (VALUES FROM EXCEL FILE, TAB 'weighting')
# predation_weight <- weighting[parameter == "predation_weight", value]
# resource_weight  <- weighting[parameter == "resource_weight", value]
predation_weight <- 1
resource_weight <- 1

# Now we need to: M <- M + Bw, so 
effects_Matrix <- (predation_weight * predator_effect_matrix) +
                  (resource_weight  * resource_effect_matrix)

###############################################################################
## 17. Calculate total direct and indirect effects
###############################################################################

effect_results <- calculate_effect_matrix(
  mat = effects_Matrix,
  R = response_strength
)

###############################################################################
## 18. Visualize network effects
###############################################################################

## Total-effect matrix
total_effect_matrix <- effect_results$total_effect_matrix
diag(total_effect_matrix) <- NA

## Direct-effect matrix
direct_effect_matrix <- effect_results$direct_effect_matrix
diag(direct_effect_matrix) <- NA

###############################################################################
## Positive and negative effects
###############################################################################

positive_effect_matrix <- pmax(total_effect_matrix, 0)
negative_effect_matrix <- pmin(total_effect_matrix, 0)

positive_effect_out <-
  colSums(
    positive_effect_matrix,
    na.rm = TRUE
  )

negative_effect_out <-
  colSums(
    negative_effect_matrix,
    na.rm = TRUE
  )

###############################################################################
## Species rankings
###############################################################################

positive_effect_ranking <-
  data.frame(
    species = names(sort(positive_effect_out, decreasing = TRUE)),
    total_positive_effect =
      sort(positive_effect_out, decreasing = TRUE)
  )

negative_effect_ranking <-
  data.frame(
    species = names(sort(negative_effect_out)),
    total_negative_effect =
      sort(negative_effect_out)
  )

###############################################################################
## Prepare heatmaps
###############################################################################

total_effect_heatmap <-
  ggplot(
    melt(total_effect_matrix),
    aes(
      Var2,
      factor(
        Var1,
        levels = rev(rownames(total_effect_matrix))
      ),
      fill = value
    )
  ) +
  geom_tile(color = "grey90") +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    na.value = "grey90"
  ) +
  coord_fixed() +
  labs(
    x = "Acting species",
    y = "Affected species",
    fill = "Total effect"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x =
      element_text(angle = 45, hjust = 1)
  )

direct_effect_heatmap <-
  ggplot(
    melt(direct_effect_matrix),
    aes(
      Var2,
      factor(
        Var1,
        levels = rev(rownames(direct_effect_matrix))
      ),
      fill = value
    )
  ) +
  geom_tile(color = "grey90") +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    na.value = "grey90"
  ) +
  coord_fixed() +
  labs(
    x = "Acting species",
    y = "Affected species",
    fill = "Direct effect"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x =
      element_text(angle = 45, hjust = 1)
  )


###############################################################################
## 19. Sensitivity analysis
##
## Evaluate how robust the network results are to uncertainty in
## species response strengths (R).
###############################################################################

n_simulations <- 5000

###############################################################################
## Baseline ranking
###############################################################################

baseline_effects <- effect_results$mean_total_effect

baseline_rank <- names(
  sort(baseline_effects,
       decreasing = TRUE)
)

###############################################################################
## Objects to store results
###############################################################################

spearman_correlations <- numeric(n_simulations)
top1_species <- character(n_simulations)
top5_species <- vector(
  mode = "list",
  length = n_simulations)
ranking_identical <- numeric(n_simulations)

###############################################################################
## Run simulations
###############################################################################

baseline_ranking <-
  names(sort(baseline_effects, decreasing = TRUE))

for(i in seq_len(n_simulations)){
  
  ## -------------------------------------------------------------
  ## Perturb response strengths (±20%)
  ## -------------------------------------------------------------
  
  response_strength_temp <-
    response_strength *
    runif(length(response_strength),
          min = 0.8,
          max = 1.2)
  
  ## Keep values within ecological bounds
  
  response_strength_temp <-
    pmin(1,
         pmax(0,
              response_strength_temp))
  
  names(response_strength_temp) <-
    names(response_strength)
  
  response_strength_temp <-
    response_strength_temp[
      rownames(direct_effect_matrix)
    ]
  
  ## -------------------------------------------------------------
  ## Recalculate effects
  ## -------------------------------------------------------------
  
  effect_results_temp <-
    calculate_effect_matrix(
      mat = effects_Matrix,
      R   = response_strength_temp
    )
 
  temp_effects <-
    effect_results_temp$mean_total_effect
  
  ## -------------------------------------------------------------
  ## Ranking stability
  ## -------------------------------------------------------------
  
  spearman_correlations[i] <-
    cor(
      baseline_effects,
      temp_effects,
      method = "spearman"
    )
  
  ranking <-
    names(
      sort(temp_effects,
           decreasing = TRUE)
    )
  
  top1_species[i] <- ranking[1]
  top5_species[[i]] <- ranking[1:5]
  
  ranking_identical[i] <-
    identical(ranking, baseline_ranking)  
}

###############################################################################
## Summary statistics
###############################################################################

sensitivity_summary <-
  data.frame(
    statistic = c(
      "Mean Spearman correlation",
      "Median Spearman correlation",
      "Minimum Spearman correlation",
      "Maximum Spearman correlation",
      "Proportion of identical rankings"
    ),
    value = c(
      mean(spearman_correlations),
      median(spearman_correlations),
      min(spearman_correlations),
      max(spearman_correlations),
      mean(ranking_identical)
    )
  )

###############################################################################
## Frequency each species is ranked first
###############################################################################

top_rank_stability <-
  sort(
    prop.table(table(top1_species)),
    decreasing = TRUE
  )

top_rank_stability <-
  data.frame(
    species = names(top_rank_stability),
    proportion_top1 =
      as.numeric(top_rank_stability)
  )

###############################################################################
## Frequency each species appears in the Top 5
###############################################################################

top5_stability <-
  sapply(names(baseline_effects), function(sp){
    
    mean(
      sapply(
        top5_species,
        function(x) sp %in% x
      )
    )
    
  })

top5_stability <-
  sort(top5_stability,
       decreasing = TRUE)

top5_stability <-
  data.frame(
    species = names(top5_stability),
    proportion_top5 = top5_stability
  )

###############################################################################
## Plot
###############################################################################

response_strength_sensitivity_plot <-
  ggplot(
    data.frame(
      correlation = spearman_correlations
    ),
    aes(correlation)
  ) +
  
  geom_histogram(
    bins = 30,
    fill = "grey70",
    colour = "black"
  ) +
  
  theme_classic() +
  
  xlab("Spearman correlation with baseline species ranking") +
  ylab("Number of simulations") +
  ggtitle("Robustness of species rankings to uncertainty in response strengths")

plot_top5_stability <-
  ggplot(
    top5_stability,
    aes(
      x = reorder(species, proportion_top5),
      y = proportion_top5
    )
  ) +
  
  geom_col(fill = "steelblue") +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Top-5 ranking stability",
    subtitle = "Proportion of simulations in which each species ranked among the five most influential",
    x = NULL,
    y = "Proportion of simulations"
  ) +
  
  theme_bw()
  
###############################################################################
## 20. Scenarios Analysis
###############################################################################

## ============================================================================
## Baseline
## ============================================================================

baseline_results <- effect_results

baseline_total_effect_matrix <-
  baseline_results$total_effect_matrix

baseline_mean_total_effect <-
  baseline_results$mean_total_effect

###############################################################################
## Define scenarios
###############################################################################

management_scenarios <- list()

###############################################################################
## Predator reductions
###############################################################################

managed_predators <-
  c(
    "Felis_catus",
    "Rattus_spp",
    "Salvator_merianae",
    "Rhinella"
  )

management_levels <-
  c(
    0.25,
    0.50,
    0.75,
    0.90,
    0.9999
  )

for(predator in managed_predators){
  for(level in management_levels){
    scenario_name <-
      paste0(
        predator,
        "_",
        round(level * 100),
        "percent"
      )
    
    management_scenarios[[scenario_name]] <-
      data.frame(
        source = predator,
        target = NA,
        reduction = level,
        stringsAsFactors = FALSE
      )
  }
  }

###############################################################################
## Combined predator eradications
###############################################################################

management_scenarios$Cats_and_Rats <-
  data.frame(
    source = c(
      "Felis_catus",
      "Rattus_spp"
    ),
    target = NA,
    reduction = 0.99999
  )

management_scenarios$Cats_and_Rhinella <-
  data.frame(
    source = c(
      "Felis_catus",
      "Rhinella"
    ),
    target = NA,
    reduction = 0.99999
  )

management_scenarios$Cats_and_Tegu <-
  data.frame(
    source = c(
      "Felis_catus",
      "Salvator_merianae"
    ),
    target = NA,
    reduction = 0.99999
  )

management_scenarios$Rats_and_Tegu <-
  data.frame(
    source = c(
      "Rattus_spp",
      "Salvator_merianae"
    ),
    target = NA,
    reduction = 0.99999
  )

management_scenarios$Rats_and_Rhinella <-
  data.frame(
    source = c(
      "Rattus_spp",
      "Rhinella"
    ),
    target = NA,
    reduction =0.99999
  )

management_scenarios$Tegu_and_Rhinella <-
  data.frame(
    source = c(
      "Salvator_merianae",
      "Rhinella"
    ),
    target = NA,
    reduction = 0.99999
  )

management_scenarios$All_Invasive <-
  data.frame(
    source = c(
      "Felis_catus",
      "Rattus_spp",
      "Salvator_merianae",
      "Rhinella"
    ),
    target = NA,
    reduction = 0.99999
  )

###############################################################################
## Run all scenarios
###############################################################################

######################## TO USE FOR THE PLOTS
management_results <- list()

for(name in names(management_scenarios)){
  management_results[[name]] <- simulate_management_scenario(
    effect_matrix = effects_Matrix,
    response_strength = response_strength,
    modifications = management_scenarios[[name]],
    redistribution = "predator_diet")
  
  management_results[[name]]$delta_received_effect <-
    rowSums(
      management_results[[name]]$
        delta_total_effect_matrix,
      na.rm = TRUE
    )
  }


########################


######################## TO USE FOR COMPARISON ###########################################################
management_results_none <- list()

for(name in names(management_scenarios)){
  management_results_none[[name]] <-
    simulate_management_scenario(
      effect_matrix = effects_Matrix,
      response_strength = response_strength,
      modifications = management_scenarios[[name]],
      redistribution = "none")
  
  management_results_none[[name]]$delta_received_effect <-
    rowSums(
      management_results_none[[name]]$
        delta_total_effect_matrix,
      na.rm = TRUE
    )
}

management_results_diet <- list()

for(name in names(management_scenarios)){
  management_results_diet[[name]] <-
    simulate_management_scenario(
      effect_matrix = effects_Matrix,
      response_strength = response_strength,
      modifications = management_scenarios[[name]],
      redistribution = "predator_diet")
  
  management_results_diet[[name]]$delta_received_effect <-
    rowSums(
      management_results_diet[[name]]$
        delta_total_effect_matrix,
      na.rm = TRUE
    )
}

comparison_trachylepis <-
  data.frame(
    scenario = names(management_scenarios),
    
    no_redistribution =
      sapply(
        management_results_none,
        function(x)
          x$delta_mean_total_effect[
            "Trachylepis_atlantica"
          ]
      ),
    
    predator_redistribution =
      sapply(
        management_results_diet,
        function(x)
          x$delta_mean_total_effect[
            "Trachylepis_atlantica"
          ]
      )
  )

comparison_trachylepis$difference <-
  comparison_trachylepis$predator_redistribution -
  comparison_trachylepis$no_redistribution
comparison_summary <-
  data.frame(
    scenario = names(management_scenarios),
    
    rank_correlation =
      sapply(
        names(management_scenarios),
        function(name){
          
          cor(
            management_results_none[[name]]$
              delta_mean_total_effect,
            
            management_results_diet[[name]]$
              delta_mean_total_effect,
            
            method = "spearman"
          )
          
        }
      ),
    
    mean_absolute_difference =
      sapply(
        names(management_scenarios),
        function(name){
          
          mean(
            abs(
              management_results_none[[name]]$
                delta_mean_total_effect -
                
                management_results_diet[[name]]$
                delta_mean_total_effect
            )
          )
          
        }
      ),
    
    max_absolute_difference =
      sapply(
        names(management_scenarios),
        function(name){
          
          max(
            abs(
              management_results_none[[name]]$
                delta_mean_total_effect -
                
                management_results_diet[[name]]$
                delta_mean_total_effect
            )
          )
          
        }
      )
  )

comparison_summary
# comparison_summary <-
#   data.frame(
#     scenario = names(management_scenarios),
#     
#     correlation =
#       sapply(
#         names(management_scenarios),
#         function(name){
#           cor(
#             management_results_none[[name]]$results$mean_total_effect,
#             management_results_diet[[name]]$results$mean_total_effect,
#             method = "spearman"
#           )
#         }
#       ),
#     
#     mean_absolute_difference =
#       sapply(
#         names(management_scenarios),
#         function(name){
#           
#           mean(
#             abs(
#               management_results_none[[name]]$
#                 delta_mean_total_effect -
#                 
#                 management_results_diet[[name]]$
#                 delta_mean_total_effect
#             )
#           )
#           
#         }
#       )
#   )

comparison_plot <-
  ggplot(
    comparison_trachylepis,
    aes(
      no_redistribution,
      predator_redistribution
    )
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  
  geom_point(size = 3) +
  
  geom_text(
    aes(label = scenario),
    hjust = -0.1,
    size = 3
  ) +
  
  theme_bw() +
  
  labs(
    x = "Without predator redistribution",
    y = "With predator redistribution",
    title = "Sensitivity of Trachylepis response to behavioural assumptions"
  )

###########################################################


###############################################################################
## Common limits for management plots
###############################################################################

## Largest species response across ALL scenarios

max_species_change <-
  max(
    sapply(
      management_results,
      function(x)
        max(abs(x$delta_mean_total_effect),
            na.rm = TRUE)
    )
  )

## Largest matrix change across ALL scenarios

max_matrix_change <-
  max(
    sapply(
      management_results,
      function(x)
        max(abs(x$delta_total_effect_matrix),
            na.rm = TRUE)
    )
  )

###############################################################################
## Scenario × Species response matrix
###############################################################################

management_species_matrix <-
  do.call(
    rbind,
    lapply(
      management_results,
      function(x)
        x$delta_mean_total_effect
    )
  )

rownames(management_species_matrix) <-
  names(management_results)

management_species_matrix <-
  management_species_matrix[
    ,
    names(sort(baseline_effects,
               decreasing = TRUE)),
    drop = FALSE
  ]

###############################################################################
## Summary table
###############################################################################

management_summary <-
  do.call(
    rbind,
    lapply(
      names(management_results),
      function(name){
        
        data.frame(
          scenario = name,
          species = names(
            management_results[[name]]$delta_mean_total_effect
          ),
          delta_mean_total_effect =
            as.numeric(
              management_results[[name]]$delta_mean_total_effect
            ),
          stringsAsFactors = FALSE
        )
        
      }
    )
  )

###############################################################################
## Number of sign changes
###############################################################################

management_sign_changes <-
  data.frame(
    scenario = names(management_results),
    sign_changes =
      sapply(
        management_results,
        function(x) x$sign_changes
      )
  )

###############################################################################
## Responses of endemic species
###############################################################################

endemic_species <-
  c(
    "Trachylepis_atlantica",
    "Amphisbaena_ridleyi",
    "Kerodon_rupestris"
  )

management_endemic_species <-
  do.call(
    rbind,
    lapply(
      names(management_results),
      function(name){
        
        data.frame(
          scenario = name,
          species = endemic_species,
          delta_mean_total_effect =
            management_results[[name]]$
            delta_mean_total_effect[endemic_species],
          stringsAsFactors = FALSE
        )
        
      }
    )
  )

###############################################################################
## Store matrices for every scenario
###############################################################################

management_direct_effect_matrices <- list()
management_total_effect_matrices <- list()
management_delta_effect_matrices <- list()

###############################################################################
## Common colour scale across all management scenarios
###############################################################################

for(name in names(management_results)){
  
  management_direct_effect_matrices[[name]] <-
    management_results[[name]]$direct_effect_matrix
  
  management_total_effect_matrices[[name]] <-
    management_results[[name]]$total_effect_matrix
  
  management_delta_effect_matrices[[name]] <-
    management_results[[name]]$delta_total_effect_matrix
}

max_delta_effect <-
  max(
    abs(
      unlist(management_delta_effect_matrices)
    ),
    na.rm = TRUE
  )
###############################################################################
## Overall management comparison
###############################################################################

management_summary2 <-
  data.frame(
    scenario =
      names(management_results),
    
    mean_change =
      sapply(
        management_results,
        function(x)
          mean(
            x$delta_mean_total_effect
          )
      ),
    
    positive_species =
      sapply(
        management_results,
        function(x)
          sum(
            x$delta_mean_total_effect > 0
          )
      ),
    
    negative_species =
      sapply(
        names(management_results),
        function(name){
          
          delta <-
            management_results[[name]]$
            delta_mean_total_effect
          
          managed_species <-
            management_scenarios[[name]]$source
          
          delta <-
            delta[
              !(names(delta) %in% managed_species)
            ]
          sum(delta < 0)
        }
      ),
    
    sign_changes =
      sapply(
        management_results,
        function(x)
          x$sign_changes
      )
  )

management_summary2$Trachylepis <-
  sapply(
    management_results,
    function(x)
      x$delta_mean_total_effect["Trachylepis_atlantica"]
  )

management_summary2$Amphisbaena <-
  sapply(
    management_results,
    function(x)
      x$delta_mean_total_effect["Amphisbaena_ridleyi"]
  )

management_summary2$Kerodon <-
  sapply(
    management_results,
    function(x)
      x$delta_mean_total_effect["Kerodon_rupestris"]
  )
###############################################################################
## Heatmaps
###############################################################################

management_heatmaps <- list()

for(name in names(management_results)){
  
  delta_matrix <-
    management_delta_effect_matrices[[name]]
  
  ## Remove managed species from the plot
  managed_species <-
    unique(
      management_scenarios[[name]]$source
    )
  
  delta_matrix <-
    delta_matrix[
      !(rownames(delta_matrix) %in% managed_species),
      !(colnames(delta_matrix) %in% managed_species),
      drop = FALSE
    ]
  
  diag(delta_matrix) <- NA
  
  df <-
    reshape2::melt(delta_matrix)
  
  colnames(df) <-
    c(
      "Receiver",
      "Source",
      "Effect"
    )
  
  management_heatmaps[[name]] <-
    
    ggplot(
      df,
      aes(
        Source,
        Receiver,
        fill = Effect
      )
    ) +
    
    geom_tile() +
    
    scale_fill_gradient2(
      low = "#762A83",
      mid = "white",
      high = "#1B7837",
      midpoint = 0,
      limits = c(-max_delta_effect,
                 max_delta_effect)
    ) +
    
    coord_fixed() +
    
    theme_bw() +
    
    theme(
      axis.text.x =
        element_text(
          angle = 90,
          hjust = 1,
          size = 8
        ),
      axis.text.y =
        element_text(size = 8)
    ) +
    labs(
      title = paste("Change relative to baseline:", name),
      x = "Effect source",
      y = "Effect receiver",
      fill = expression(Delta*" total effect")
    )
}

###############################################################################
## Ecological responses to management
###############################################################################

max_ecological_change <-
  max(
    sapply(
      management_results,
      function(x)
        max(
          abs(x$delta_received_effect),
          na.rm = TRUE
        )
    )
  )

############################################################
## Only plot complete eradication scenarios
############################################################

plot_scenarios <-
  c(
    "Felis_catus_100percent",
    "Rattus_spp_100percent",
    "Salvator_merianae_100percent",
    "Rhinella_100percent",
    "All_Invasive",
    "Cats_and_Rats",
    "Cats_and_Rhinella",
    "Cats_and_Tegu",
    "Rats_and_Tegu",
    "Rats_and_Rhinella",
    "Tegu_and_Rhinella"
  )

management_ecological_plots <- list()
released_pressure_plots <- list()

for(name in plot_scenarios){
  df <-
    data.frame(
      species =
        names(
          management_results[[name]]$
            delta_received_effect
        ),
      delta =
        as.numeric(
          management_results[[name]]$
            delta_received_effect
        )
    )
  
  managed_species <-
    management_scenarios[[name]]$source
  
  df <-
    df[
      !(df$species %in% managed_species),
    ]
  
  df <-
    df[
      order(df$delta),
    ]
  
  cat(name, "\n")
  print(dim(df))
  print(head(df))
  
  management_ecological_plots[[name]] <-
    
    ggplot(
      df,
      aes(
        x = reorder(species, delta),
        y = delta,
        fill = delta > 0
      )
    ) +
    
    geom_col() +
    
    coord_flip() +
    
    scale_y_continuous(
      limits = c(
        -max_ecological_change,
        max_ecological_change
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "#B2182B",
        "#1B7837"
      ),
      guide = "none"
    ) +
    
    geom_hline(
      yintercept = 0,
      colour = "black"
    ) +
    
    labs(
      title = paste(
        "Ecological response:",
        name
      ),
      subtitle =
        "Positive values indicate species experiencing less total pressure",
      x = NULL,
      y = expression(Delta*" received pressure")
    ) +
    
    theme_bw()
  

  ############################################################
  ## Direct predation released
  ############################################################
  
  managed_species <-
    unique(
      management_scenarios[[name]]$source
    )
  
  ############################################################
  ## Build matrix containing ONLY removed predators
  ############################################################
  
  released_matrix <-
    effect_results$direct_effect_matrix[
      ,
      managed_species,
      drop = FALSE
    ]
  
  ############################################################
  ## Convert predation pressure into released pressure
  ############################################################
  
  released_matrix <-
    -released_matrix
  
  released_matrix[
    released_matrix <= 0
  ] <- NA
  
  ############################################################
  ## Removed predators should not appear as receivers
  ############################################################
  
  released_matrix <-
    released_matrix[
      !(rownames(released_matrix) %in% managed_species),
      ,
      drop = FALSE
    ]
  
  ############################################################
  ## Long format
  ############################################################
  
  df <-
    reshape2::melt(
      released_matrix,
      na.rm = TRUE
    )
  
  colnames(df) <-
    c(
      "Receiver",
      "Source",
      "Released"
    )
  
  ############################################################
  ## Plot
  ############################################################
  
  released_pressure_plots[[name]] <-
    
    ggplot(
      df,
      aes(
        Source,
        Receiver,
        fill = Released
      )
    ) +
    
    geom_tile() +
    
    scale_fill_gradient(
      low = "white",
      high = "#2166AC"
    ) +
    
    coord_fixed() +
    
    theme_bw() +
    
    theme(
      axis.text.x =
        element_text(
          angle = 90,
          hjust = 1,
          size = 8
        ),
      axis.text.y =
        element_text(size = 8)
    ) +
    
    labs(
      title = paste(
        "Direct predation released:",
        name
      ),
      x = "Removed predator",
      y = "Released prey",
      fill = "Released\npressure"
    )
}

###############################################################################
## Mean ecosystem improvement
###############################################################################

plot_management_summary <-
  
  ggplot(
    management_summary2,
    aes(
      x = reorder(
        scenario,
        mean_change
      ),
      y = mean_change,
      fill = mean_change > 0
    )
  ) +
  
  geom_col() +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(
      -max(abs(management_summary2$mean_change)),
      max(abs(management_summary2$mean_change))
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "#B2182B",
      "#1B7837"
    ),
    guide = "none"
  ) +
  geom_hline(
    yintercept = 0
  ) +
  labs(
    title = "Overall ecosystem response",
    x = NULL,
    y = expression(Delta*" mean total effect")
  ) +
  theme_bw()

###############################################################################
## 20. Save outputs
###############################################################################

## ---------------------------------------------------------------------------
## Save plots
## ---------------------------------------------------------------------------

ggsave(
  filename = "outputs/Figure_01_Foodweb.png",
  plot = foodweb_plot,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = "outputs/Figure_02_DirectEffectsHeatmap.png",
  plot = direct_effect_heatmap,
  width = 10,
  height = 9,
  dpi = 300
)

ggsave(
  filename = "outputs/Figure_03_TotalEffectsHeatmap.png",
  plot = total_effect_heatmap,
  width = 10,
  height = 9,
  dpi = 300
)

ggsave(
  filename = "outputs/Figure_04_SensitivityHistogram.png",
  plot = response_strength_sensitivity_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "outputs/Figure_09_TopSpeciesStability.png",
  plot = plot_top5_stability,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "outputs/Figure_10_Management_summary.png",
  plot_management_summary,
  width = 8,
  height = 6,
  dpi = 300
)

###############################################################################
## Save management scenario heatmaps
###############################################################################

for(name in names(management_heatmaps)){
  
  ggsave(
    filename =
      file.path(
        "outputs",
        paste0(
          "Management_",
          gsub("[^A-Za-z0-9]", "_", name),
          ".png"
        )
      ),
    plot = management_heatmaps[[name]],
    width = 9,
    height = 8,
    dpi = 300
  )
  
}

############################################################
## Save direct-release plots
############################################################

for(name in names(released_pressure_plots)){
  
  ggsave(
    filename = file.path(
      "outputs",
      paste0(
        "DirectPredationReleased_",
        name,
        ".png"
      )
    ),
    plot = released_pressure_plots[[name]],
    width = 8,
    height = 7,
    dpi = 300
  )
  
}
###############################################################################
## Save ecological-response plots
###############################################################################

for(name in names(management_ecological_plots)){
  ggsave(
    filename =
      file.path(
        "outputs",
        paste0(
          "EcologicalResponse_",
          gsub("[^A-Za-z0-9]", "_", name),
          ".png"
        )
      ),
    plot = management_ecological_plots[[name]],
    width = 8,
    height = 6,
    dpi = 300
  )
}

## ---------------------------------------------------------------------------
## Save tables
## ---------------------------------------------------------------------------

write.csv(
  direct_effect_matrix,
  "outputs/direct_effect_matrix.csv",
  row.names = TRUE
)

write.csv(
  total_effect_matrix,
  "outputs/total_effect_matrix.csv",
  row.names = TRUE
)

write.csv(
  effects_Matrix,
  "outputs/effects_input_matrix.csv",
  row.names = TRUE
)

write.csv(
  sensitivity_summary,
  "outputs/sensitivity_summary.csv",
  row.names = FALSE
)

write.csv(
  management_summary,
  "outputs/management_summary.csv",
  row.names = FALSE
)

###############################################################################
## Save scenario summaries
###############################################################################

write.csv(
  management_sign_changes,
  "outputs/management_sign_changes.csv",
  row.names = FALSE
)

write.csv(
  management_endemic_species,
  "outputs/management_endemic_species.csv",
  row.names = FALSE
)

write.csv(
  management_summary2,
  "outputs/management_summary2.csv",
  row.names = FALSE
)

## ---------------------------------------------------------------------------
## Save complete objects
## ---------------------------------------------------------------------------

saveRDS(
  effect_results,
  "outputs/effect_results.rds"
)

saveRDS(
  management_results,
  "outputs/management_results.rds"
)

saveRDS(
  management_summary,
  "outputs/management_summary.rds"
)

saveRDS(
  sensitivity_summary,
  "outputs/sensitivity_summary.rds"
)

###############################################################################
## Save all scenario matrices
###############################################################################

saveRDS(
  management_direct_effect_matrices,
  "outputs/management_direct_effect_matrices.rds"
)

saveRDS(
  management_total_effect_matrices,
  "outputs/management_total_effect_matrices.rds"
)

saveRDS(
  management_delta_effect_matrices,
  "outputs/management_delta_effect_matrices.rds"
)

cat("\nOutputs successfully saved to 'outputs/'\n")
  