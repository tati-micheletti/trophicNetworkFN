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
## 2. Read input workbook
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

###############################################################################
## 3. Source Functions
###############################################################################

source("functions/helpers.R")

###############################################################################
## 4. Harmonize prey names
###############################################################################

cats  <- harmonize_taxonomy(cats,  mapping, finalNaming = "Ecological")
rats  <- harmonize_taxonomy(rats,  mapping, finalNaming = "Ecological")
tegu  <- harmonize_taxonomy(tegu,  mapping, finalNaming = "Ecological")
toads <- harmonize_taxonomy(toads, mapping, finalNaming = "Ecological")

###############################################################################
# 5. Aggregate prey types
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
## 19. Save outputs
###############################################################################

ggsave(
  "outputs/FoodWeb.png",
  foodweb_plot,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  "outputs/TotalEffectMatrix.png",
  total_effect_heatmap,
  width = 10,
  height = 9,
  dpi = 300
)

ggsave(
  "outputs/DirectEffectMatrix.png",
  direct_effect_heatmap,
  width = 10,
  height = 9,
  dpi = 300
)

###############################################################################
## Save matrices
###############################################################################

write.csv(
  effects_Matrix,
  "outputs/EffectsMatrix.csv"
)

write.csv(
  direct_effect_matrix,
  "outputs/DirectEffectMatrix.csv"
)

write.csv(
  total_effect_matrix,
  "outputs/TotalEffectMatrix.csv"
)

write.csv(
  positive_effect_matrix,
  "outputs/PositiveEffectMatrix.csv"
)

write.csv(
  negative_effect_matrix,
  "outputs/NegativeEffectMatrix.csv"
)

###############################################################################
## Save rankings
###############################################################################

write.csv(
  positive_effect_ranking,
  "outputs/PositiveEffectRanking.csv",
  row.names = FALSE
)

write.csv(
  negative_effect_ranking,
  "outputs/NegativeEffectRanking.csv",
  row.names = FALSE
)

###############################################################################
## Save summary statistics
###############################################################################

summary_results <- data.frame(
  metric = c(
    "Mean_total_effect",
    "Proportion_indirect_effects"
  ),
  value = c(
    mean(effect_results$mean_total_effect),
    effect_results$proportion_indirect_effects
  )
)

write.csv(
  summary_results,
  "outputs/NetworkSummary.csv",
  row.names = FALSE
)

###############################################################################
## Save complete R object
###############################################################################

saveRDS(
  effect_results,
  "outputs/EffectResults.rds"
)

saveRDS(
  foodweb_model,
  "outputs/FoodWebModel.rds"
)

