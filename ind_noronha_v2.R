# =============================================================================
# Total Effect Analysis — Fernando de Noronha Food Web
# Computes the T matrix of total (direct + indirect) effects using
# PSIRI-based interaction strengths and density-weighted predator consumption.
# Includes sensitivity analysis on R and predator removal simulations.
# =============================================================================
#Total effect analysis based on dietary data for the Noronha food web
Require::Require("data.table")
Require::Require("reshape")
Require::Require("ggplot2")
Require::Require("cowplot")
# Require::Require("dplyr")

psiri <- as.matrix(read.csv("data/temp/psiri.csv", header = TRUE, row.names = 1)) #PSRI for predators
rawFreq <- as.matrix(read.csv("data/temp/rawFreq.csv", header = TRUE, row.names = 1)) #Raw frequencies of resource use computed from dietary studies
spp <- read.csv("data/temp/species.csv", header = TRUE)

ord <- c('r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'm1', 'm2', 'p1','p2', 'p3')
psiri <- psiri[ord,ord]
rawFreq <- rawFreq[ord,ord]
spp <- as.data.table(spp)
spp[, abr := factor(abr, levels = ord)]
setkey(spp, "abr")

lookup <- setNames(spp$name, as.character(spp$abr))
rownames(psiri) <- rownames(rawFreq) <- lookup[rownames(psiri)]
colnames(psiri) <- colnames(rawFreq) <- lookup[colnames(psiri)]

# =============================================================================
# Build interaction matrix M
# Columns 8-11 are the four invasive predators: cats, tegu, rats, toads.
# Per-capita consumption = raw frequency / sample size, then scaled by density.
# =============================================================================
Bw <- rawFreq
samplesInvaders <- c("felis_catus" = 78, # cats (n=78) Gaiotto et al., 2020
                     "salvator_merianae" = 22, # tegu (n=22) Gaiotto et al., 2020
                     "rattus_rattus" = 10, # rats (n=10) Gaiotto et al., 2020
                     "rhinella" = 66) # ???? toads? (n=143) Tolledo & Toledo (2015) Not matching the text!

Bw[,names(samplesInvaders)] <- Bw[,names(samplesInvaders)]/matrix(samplesInvaders, 
                                                                  ncol = 4, 
                                                                  nrow = nrow(rawFreq), 
                                                                  byrow = TRUE) # Average individual consumption (Frequency/sample size)

densities <- c("felis_catus" = 0.71, # feral cats, 0.71 ind/ha , Dias et al. 2017
               "salvator_merianae" = 3.98, # tegu 3.98 ind/ha Abrahão et al. 2019
               "rattus_rattus" = 37, # rats 37 ind/ha, Russell et al. 2018
               "rhinella" = 10.35) # Toads, extrapolation from Solomon Islands, Pikacha et al. 2015

Bw[,names(samplesInvaders)] <- Bw[,names(samplesInvaders)]*matrix(densities, 
                                                                  ncol = 4, 
                                                                  nrow = nrow(rawFreq), 
                                                                  byrow = TRUE) # weighting consumption by estimated density of the predator

# M matrix - effects of columns over rows
M <- 0.1 * t(psiri) # effect of prey over predators - energetic efficiency of predators is low
M2 <- M  + Bw # adding the effect of predators over prey


# R.vec: species-specific interaction dependency / sensitivity parameters
# (one value per species, controls how strongly each species responds to net effects)
R.vec <- c(
  invertebrados_ter     = 0.1,
  invertebrados_aq      = 0.1,
  rodentia               = 0.5,
  kerodon_rupestris      = 0.5,
  birds                  = 0.1,
  amphisbaena_ridleyi    = 0.2,
  trachylepis_atlantica  = 0.5,
  rattus_rattus          = 0.8,
  felis_catus            = 0.8,
  salvator_merianae      = 0.8,
  rhinella               = 0.8
)

stopifnot(setequal(names(R.vec), colnames(M)))  # catches typos/missing species
R.vec <- R.vec[colnames(M)]  # reorders defensively to match M regardless of how R.vec was typed

#######################
#                     #
# FUNCTION indirect   #
#                     #
#######################

# indirect()
# Computes the T matrix of total (direct + indirect) effects.
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
#==========================================================================

# indirect <- function(mat, R) {
#   
#   ## -------------------------------------------------------------------------
#   ## Input checks
#   ## -------------------------------------------------------------------------
#   
#   stopifnot(is.matrix(mat))
#   stopifnot(nrow(mat) == ncol(mat))
#   
#   stopifnot(is.numeric(R))
#   stopifnot(length(R) == nrow(mat))
#   
#   stopifnot(!anyNA(mat))
#   stopifnot(!anyNA(R))
#   
#   S <- nrow(mat)
#   
#   ## -------------------------------------------------------------------------
#   ## Build dependency matrix (W)
#   ##
#   ## Each row is divided by its total interaction strength so rows sum to 1.
#   ## -------------------------------------------------------------------------
#   
#   row_totals <- rowSums(mat)
#   stopifnot(all(row_totals > 0))
#   W <- mat / row_totals
#   
#   ## -------------------------------------------------------------------------
#   ## Sign convention
#   ##
#   ## Upper triangle = antagonistic effects
#   ##
#   ## (We'll probably improve this later because it currently depends on species
#   ## ordering.)
#   ## -------------------------------------------------------------------------
#   
#   consumers <- spp$name[spp$is_consumer]
#   prey <- setdiff(rownames(W), consumers)
#   
#   W[prey, consumers] <- -W[prey, consumers]
#   
#   ## -------------------------------------------------------------------------
#   ## Interaction matrix
#   ## -------------------------------------------------------------------------
#   
#   P <- diag(R)
#   
#   ## -------------------------------------------------------------------------
#   ## Total effects matrix
#   ##
#   ## T = (I - P W)^(-1)
#   ## -------------------------------------------------------------------------
#   
#   I <- diag(S)
#   
#   Tmat <- solve(I - P %*% W)
#   
#   ## Ignore self-effects when computing summaries
#   
#   T_no_diag <- Tmat
#   diag(T_no_diag) <- NA
#   
#   ## -------------------------------------------------------------------------
#   ## Mean contribution (Tout)
#   ## -------------------------------------------------------------------------
#   
#   mean_contribution <- colMeans(T_no_diag, na.rm = TRUE)
#   
#   ## -------------------------------------------------------------------------
#   ## Binary adjacency matrix indicating the presence (1) or absence (0)
#   ## of direct interactions. Used to distinguish species pairs connected
#   ## only through indirect pathways from those with direct interactions.
#   ## -------------------------------------------------------------------------
#   
#   B <- (mat > 0) * 1
#   
#   ## -------------------------------------------------------------------------
#   ## Proportion of indirect effects
#   ## -------------------------------------------------------------------------
#   
#   indirect_only <- T_no_diag * (1 - B)
#   
#   prop_ind_eff <-
#     sum(indirect_only, na.rm = TRUE) /
#     sum(T_no_diag, na.rm = TRUE)
#   
#   prop_ind_eff_tin <-
#     rowSums(indirect_only, na.rm = TRUE) /
#     rowSums(T_no_diag, na.rm = TRUE)
#   
#   prop_ind_eff_tout <-
#     colSums(indirect_only, na.rm = TRUE) /
#     colSums(T_no_diag, na.rm = TRUE)
#   
#   ## -------------------------------------------------------------------------
#   ## Return
#   ## -------------------------------------------------------------------------
#   
#   list(
#     mean_contribution = mean_contribution,
#     prop_ind_eff      = prop_ind_eff,
#     prop_ind_eff_tin  = prop_ind_eff_tin,
#     prop_ind_eff_tout = prop_ind_eff_tout,
#     direct            = P %*% W,
#     Tmat              = Tmat
#   )
# }

#######################

# T.matrix <- indirect(M2, R = R.vec, type = "unipartite") 
T.matrix <- calculate_effect_matrix(M2, R = R.vec)

Tm <- T.matrix$Tmat

#Negative Tout
NT <- pmin(Tm, 0)

NTout <- colSums(NT)
sort(NTout)
hist(as.numeric(Tm), col = 'darkgrey')

#Positive Tout
PT <- (Tm > 0)*Tm
PTout <- colSums(PT)
sort(PTout, decreasing = TRUE)

mean(NT)
hist(as.numeric(PTout), col = 'darkgrey')

# =============================================================================
# Visualisation — heatmap and histograms of T matrix values
# =============================================================================
##### plotting ####
####ggplot

# Transform the matrix in long format
Tm_plot <- Tm
diag(Tm_plot) <- NA
df <- as.data.table(as.table(Tm_plot))
colnames(df) <- c("x", "y", "value")

p1 <- ggplot(df, aes(x = x, y = ordered(y, levels = rev(sort(unique(y)))), fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#FF0000",
                       mid = "#FFFFFF",
                       high = "#87CEFA", n.breaks=15, limits = c(-0.5,1), na.value = "lightgrey") +
  coord_fixed()+
  xlab("")+ ylab("")

#################### hist #############

df<-data.frame(as.numeric(PT))
colnames(df)<-"value"

a<-ggplot(df, aes(x=value)) +
  geom_histogram(fill="darkgrey", color="black", alpha=1) +
  ylim(c(0,80))+
  theme_classic()+ xlab("Positive T matrix values")+ ylab("Count")

df2<-data.frame(as.numeric(NT))
colnames(df)<-"value"

b<-ggplot(df,aes(x=value)) +
  geom_histogram(fill="darkgrey", color="black", alpha=1) +
  ylim(c(0,80))+
  theme_classic()+ xlab("Negative T matrix values")+ ylab("Count")

plot_grid(b, a, labels = c('psiri', 'rawFreq'))

# Shouldn't I remove below?
# pdf("hist.pdf", width = 14, height = 8)
# data %>%
#   ggplot(aes(x=value, fill=type)) +
#   geom_histogram(  position = 'identity') +
#   scale_fill_manual(name="Type",labels=c("Negative","Positive", "Zero"),values=c("skyblue1", "tomato", "grey")) +
#   theme_classic(base_size = 15) + xlab("T matrix values")+ ylab("Count")
# dev.off()

## PLOTTING



#============================================
# Sensitivity Analysis: robustness of T_out species rankings to choice of R.
# Draws 1000 random R vectors; Spearman-correlates resulting T_out ranks with
# baseline to check whether influence rankings are stable regardless of R.
#Sensitivity test to R values

#1. Change R
#2. correlation between Tout values

r.neg <- numeric(1000)
r.pos <- numeric(1000)
for(i in seq_len(1000)){
  R.temp <- runif(nrow(psiri)) #randomly assigning R
  Tm.temp <- indirect(M, R = R.temp) 
  Tm.temp <- Tm.temp$Tmat
  
  NT.temp <- (Tm.temp < 0)*Tm.temp
  NTout.temp <- colSums(NT.temp)
  r.neg[i] <- cor(NTout, NTout.temp, method = 'spearman')
  
  PT.temp <- (Tm.temp > 0)*Tm.temp
  PTout.temp <- colSums(PT.temp)
  r.pos[i] <- cor(PTout, PTout.temp, method = 'spearman')
}
  par(mfrow = c(1,2))
  hist(r.neg , col = "darkgrey", border = "white",
       xlab = "correlation coefficient (r)", main = "Negative Tout")
  hist(r.pos , col = "darkgrey", border = "white",
       xlab = "correlation coefficient (r)", main = "Positive Tout")
  
  
  mean(r.neg)
  mean(c(r.pos,r.neg))
#even when randomly assigning R values, the relative rank of species influence is consistent
  

  
  
  
# #=========================================================
# #Removal Simulations
# #how removing the predators change the net effects
# #1. Number of signal shifts
# #2. delta Tout
# #3. Changes in the effect over endemic species (Trachileps)
#   
#   #TODO: simulate several eradication scenarios
#   species_removed <- c("felis_catus", "salvator_merianae")
#   M.rem <- M
#   M.rem[, species_removed] <- 0.1 * M.rem[, species_removed]
#   
#   Tm.rem <- indirect(M.rem, R = R.vec) 
#   
#   Tm.rem <- Tm.rem$Tmat
#   colnames(Tm.rem) <- colnames(M)
#   
#   #Negative Tout
#   NT <- (Tm.rem < 0)*Tm.rem
#   NTout <- colSums(NT)
#   order(abs(NTout),decreasing = T)
#   
#   #Positive Tout
#   PT <- (Tm > 0)*Tm
#   PTout <- colSums(PT)
#   order(abs(PTout),decreasing = T)
#   
#   
#   Tm[,1],Tm.rem[,1]
