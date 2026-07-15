#Total effect analysis based on dietary data for the Noronha food web

#setwd("D:/Dropbox/Backup/Colaboracoes/MS-Noronha/data")
A <- as.matrix(read.table("psiri.txt", header = T, row.names = 1)) #PSRI for predators
B <- as.matrix(read.table("geral.txt", header = T, row.names = 1)) #Raw frequencies of resource use computed from dietary studies
spp <- read.table("sp_list.txt", header = T)

ord <- c("r1", 'r2', 'r3', 'r4', 'r5', 'r6', 'm1', 'm2', 'p1','p2', 'p3')

A <- A[ord,ord]
B <- B[ord,ord]

Bw <- B
Bw[,8:11] <- B[,8:11]/matrix(c(78,22,10,66), ncol = 4, nrow = nrow(B), byrow = T ) #Average individual consumption (Frequency/sample size)

densities <- c(0.71, 3.98, 37, 10.35)
Bw[,8:11] <- Bw[,8:11]*matrix(densities, ncol = 4, nrow = nrow(B), byrow = T ) #weighting consumption by estimated density of the predator

rowSums(A)
rowSums(Bw)


#M matrix - effects of columns over rows
M <- 0.1*t(A) # effect of prey over predators - energetic efficiency of predators is low
M <- M  + Bw#adding the effect of predators over prey

R.vec <- c(0.1,0.1,0.5,0.5,0.5,0.1,0.2,0.8,0.8,0.8,0.8)

T.matrix <- indirect(M, R = R.vec, type = "unipartite") 

Tm <- T.matrix$Tmat
colnames(Tm) <- colnames(M)

#Pred.eff <- Tm[,1:4]
#rownames(Pred.eff) <- spp[,1]

#Negative Tout
NT <- (Tm < 0)*Tm
NTout <- colSums(NT)
sort(NTout)
hist(as.numeric(Tm), col = 'darkgrey')


#Positive Tout
PT <- (Tm > 0)*Tm
PTout <- colSums(PT)
sort(PTout, decreasing = T)

mean(NT)

hist(as.numeric(PTout), col = 'darkgrey')

##### plotting ####
####ggplot
install.packages("reshape")
library(reshape)
library(ggplot2)

# Transform the matrix in long format
diag(Tm)=NA

df <- melt(Tm)
colnames(df) <- c("x", "y", "value")


pdf("graph2.pdf", width = 14, height = 8)
ggplot(df, aes(x = x, y = ordered(y, levels = rev(sort(unique(y)))), fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#FF0000",
                       mid = "#FFFFFF",
                       high = "#87CEFA", n.breaks=15, limits = c(-0.5,1), na.value = "lightgrey") +
  coord_fixed()+
  xlab("")+ ylab("")
dev.off()

#################### hist #############

df<-data.frame(as.numeric(PT))
colnames(df)<-"value"

library(dplyr)
library(cowplot)

pdf("hist_pos.pdf", width = 14, height = 8)
a<-ggplot(df, aes(x=value)) +
  geom_histogram(fill="darkgrey", color="black", alpha=1) +
  ylim(c(0,80))+
  theme_classic()+ xlab("Positive T matrix values")+ ylab("Count")
dev.off()

df2<-data.frame(as.numeric(NT))
colnames(df)<-"value"

pdf("hist_neg.pdf", width = 14, height = 8)
b<-ggplot(df,aes(x=value)) +
  geom_histogram(fill="darkgrey", color="black", alpha=1) +
  ylim(c(0,80))+
  theme_classic()+ xlab("Negative T matrix values")+ ylab("Count")
dev.off()

plot_grid(b, a, labels = c('A', 'B'))


data<-read.table("clipboard", header=T)
data

pdf("hist.pdf", width = 14, height = 8)
data %>%
  ggplot(aes(x=value, fill=type)) +
  geom_histogram(  position = 'identity') +
  scale_fill_manual(name="Type",labels=c("Negative","Positive", "Zero"),values=c("skyblue1", "tomato", "grey")) +
  theme_classic(base_size = 15) + xlab("T matrix values")+ ylab("Count")
dev.off()
p



library(plotrix)

# generate colors that show negative values in red to brown
# and positive in blue-green to green
cellcol<-matrix(rep("#000000",100),nrow=10)
cellcol[Tm<0]<-color.scale(Tm[Tm<0],c(0,1),c(0,0),0.2)
cellcol[Tm>0]<-color.scale(Tm[Tm>0],0,c(0.2,0),c(0,1))


cellcol=color.scale(Tm,extremes=c("#FF7F50","#ADD8E6"))
color2D.matplot(Tm,cellcolors=cellcol,xlab="Columns",ylab="Rows",
                do.hex=FALSE,border="#F5F5F5")


pal = colorRampPalette(c("steelblue", "white", "darkorange"))
val.grad <- c(do.breaks(c(-min_max[1], 0),49), do.breaks(c(0, min_max[2]),50))

heatmap(as.matrix(Tm), Colv = NA,Rowv = NA,col = hcl.colors(200, palette = "RdBu"))


color2D.matplot(Tm,extremes=c("#FF6347","#FFA07A","white","#E0FFFF","#B0E0E6"),border="white", axes=FALSE, xlab="", ylab="")

color2D.matplot(Tm,extremes=c("#FF6347","white","#B0E0E6"),border="white", axes=FALSE, xlab="", ylab="")

color2D.matplot(Tm,extremes=c("white","#2B83BA","#ABDDA4","#FFFFBF","#FDAE61","#D7191C"),border=NA, axes=FALSE, xlab="", ylab="")


#============================================
#Sensitivity test to R values

#1. Change R
#2. correlation between Tout values

r.neg <- c()
r.pos <- c()
for(i in 1:1000){
  R.temp <- runif(nrow(A)) #randomly assigning R
  Tm.temp <- indirect(M, R = R.temp, type = "unipartite") 
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
  

#=========================================================
#Removal Simulations
#how removing the predators change the net effects
#1. Number of signal shifts
#2. delta Tout
#3. Changes in the effect over endemic species (Trachileps)

  M.rem <- M
  M.rem[,1] <- 0.1*M.rem[,1]
  Tm.rem <- indirect(M.rem, R = R.vec, type = "unipartite") 
  
  Tm.rem <- Tm.rem$Tmat
  colnames(Tm.rem) <- colnames(M)
  
  #Negative Tout
  NT <- (Tm.rem < 0)*Tm.rem
  NTout <- colSums(NT)
  order(abs(NTout),decreasing = T)
  
  #Positive Tout
  PT <- (Tm > 0)*Tm
  PTout <- colSums(PT)
  order(abs(PTout),decreasing = T)
  
  
Tm[,1],Tm.rem[,1]

# Computes indirect effects from a matrix (mat) and a vector of interaction dependencies (R)
# R is a single vector with length = sum(nrow(mat), ncol(mat)) corresponding to the R value of the ROWS and COLUMNS
# Example: 
# indirect(mat=matrix(rnbinom(Na*Np,1,0.2),nrow=Na,ncol=Np),R=0.5)
#==========================================================================
indirect <- function(mat, R, type = "bipartite"){
  if(!is.matrix(mat)){
    stop('mat must be a matrix')
  }
  if(!is.numeric(R)){
    stop('R must be a numeric vector')
  }
  
  if(type == "bipartite"){
    stop('R must be a numeric vector')
    
    Na <- nrow(mat)
    Np <- ncol(mat)
    S  <- Na + Np # species richness
    
    if(length(R) != S){
      stop('R must be a vector of length = nrow(mat) + ncol(mat)')
    }
    
    # adjacency matrix
    full <- matrix(0, S, S)
    full[1:Na, (Na+1):S] <- mat
    full[(Na+1):S,1:Na]  <- t(mat)
  }else{
    full <- mat 
    
    S <- nrow(mat)
  }
  # Bidependence matrix - dependence of rows over columns
  W <- full/rowSums(full)
  
  #W[which(lower.tri(M))] <- (-1)*W[which(lower.tri(M))] #only for antagonisms
  W[which(upper.tri(M))] <- (-1)*W[which(upper.tri(M))] #only for antagonisms
  
  # Identity matrix
  I <- diag(S)
  
  # interaction effect matrix
  P <- diag(R, S, S)
  
  # T matrix (effects of the columns on the rows )
  MT <- I-(P%*%W)
  MT <- solve(MT)
  MT. <- matrix(0, S ,S)
  
  #Nakajima and Higashi 1995
  #for (i in 1:S){
  # for (j in 1:S){
  #  MT.[i,j] <- MT[i,j]/((MT[i,i]*MT[j,j])-(MT[i,j]*MT[j,i])) 
  #} 
  #} 
  #diag(MT.) <- NA
  
  # MT <- MT-I
  
  #Disconsidering auto-indirect effects
  MT_zdiag <- MT
  diag(MT_zdiag) <- NA
  
  # In_ind<-apply(MT,1,sum) # effect suscetibility --> 1/(1-P)
  #sum_Out.ind <- colSums(MT_zdiag, na.rm = T) #indirect effect contribution
  #Out.ind <- Out.ind/sum(MT_zdiag)
  Out.ind <- apply(MT_zdiag,2,mean, na.rm = TRUE)
  # Computing ratio of indirect/direct effects
  # Binary matrix
  B <- full
  B[B>0] <- 1
  
  # Proportion of indirect effects in the whole network
  Ind.net <- sum(MT_zdiag*(1-B),na.rm=T)/sum(MT_zdiag, na.rm=TRUE) 
  
  # Proportion of indirect effects over each species
  Ind.sp.in <- apply(MT_zdiag*(1-B),1,sum,na.rm=T)/apply(MT_zdiag,1,sum,na.rm=TRUE) 
  # Proportion of indirect effects of each species
  Ind.sp.out <- apply(MT_zdiag*(1-B),2,sum,na.rm=T)/apply(MT_zdiag,2,sum,na.rm=TRUE) 
  
  results <- list(
    #summed_contribution = sum_Out.ind, #total species contribution
    mean_contribution = Out.ind, #mean species contribution
    prop_ind_eff = Ind.net, #proportion of indirect effects
    prop_ind_eff_tin = Ind.sp.in, 
    prop_ind_eff_tout = Ind.sp.out,
    direct = (P%*%W),
    Tmat = MT #Total effects matrix
    
  )
  return(results)
  
}
