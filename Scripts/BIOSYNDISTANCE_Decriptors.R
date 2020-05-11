####Chemoinformatics

#goal: extraction of matrix of physiochemical properties from list of PUBCHEMIDs

#packages used:
#if (!requireNamespace("BiocManager", quietly=TRUE))
# install.packages("BiocManager")
#BiocManager::install("ChemmineR")

#install & library require CRAN packages
requiredpackages <- c("rJava", "rcdklibs", "iterators", "rcdk","ade4","philentropy","vegan","reshape2","dendextend")

for (pkg in requiredpackages) {
  if (pkg %in% rownames(installed.packages()) == FALSE)
  {install.packages(pkg)}
  if (pkg %in% rownames(.packages()) == FALSE)
  {library(pkg, character.only=T)}
}



#SETWORKING DIRECTORY for the data folder 
setwd("Data/")


#import data with compounds as rows and predictors as columns
#IMPORT ENZYME DATA
ENZYMEDATABASE <- read.csv("Global_Compound_by_Enzyme_Table.csv")

ENZYMES <- ENZYMEDATABASE
ENZYMES <- ENZYMES[!duplicated(ENZYMES$PUBCHEMID), ]


#convert rownames to pubchem ID
ENZYMES <- data.frame(ENZYMES[, -1], row.names = ENZYMES[, 1])
#chemistry
CHEMISTRY <- read.csv("3D_Global_Compound_by_Enzyme_Table_DESCRIPTORS.csv")
CHEMISTRY <- CHEMISTRY[,-1]
CHEMISTRY <- CHEMISTRY[!duplicated(CHEMISTRY), ]

#convert rownames to pubchem ID
CHEMISTRY <- data.frame(CHEMISTRY[, -1], row.names = CHEMISTRY[, 1])

# in order to incorporate negative values into the bray curtis 
#each column had the absolute value of the minimum value per column added to each cell to presereve inter intra compound variation

i<-1
for (i in i:dim(CHEMISTRY)[2]) {
  columnminimum<-min(CHEMISTRY[,i])
  CHEMISTRY[,i]<- CHEMISTRY[,i] + abs(columnminimum)
}



#read in chemical class to color tanglegram chemical class
CLASS <- read.csv("Compound_Class.2.8.19.csv")

CLASS <- CLASS[!duplicated(CLASS$PUBchemID), ]
#convert rownames to pubchem ID
CLASS <- data.frame(CLASS, row.names = CLASS[, 1])



#full tanglegram


#create distance matricies
CHEMISTRY.DIST<-vegdist(CHEMISTRY, method = "bray")
ENZYMES.DIST<-vegdist(ENZYMES,method = "bray")


#mantel correlations
ENZYCHEM.mantel <-
  mantel.randtest(ENZYMES.DIST, CHEMISTRY.DIST, nrepet = 99999)
ENZYCHEM.mantel

#create a directory to store Global outputs
dir.create("../Output_Tables/Global")


#create dataframe to store outputs 
GlobalOutput<-data.frame(Dataset=0,MantelR=0, SimPval=0, HypothesisTest=0, Std.obs=0, Expectation=0, Variance=0)
GlobalOutput[1,]<-c("Global",ENZYCHEM.mantel$obs,ENZYCHEM.mantel$pvalue ,ENZYCHEM.mantel$alter,ENZYCHEM.mantel$expvar[1],ENZYCHEM.mantel$expvar[2],ENZYCHEM.mantel$expvar[3])
#store output
write.csv(GlobalOutput,paste("../Output_Tables/Global/","Globalmantel.csv",sep=""))

#create a directory to store Global outputs
dir.create("../Output_Plots/Global")

#store Plots of simulations 
pdf(paste("../Output_Plots/Global/","Global_mantel_Chemistry_Enzyme.pdf",sep="")) # writing a PDF to file
plot(ENZYCHEM.mantel)                   # makes the actual plot
dev.off()


#begin to create dendrograms  

#create dendrograms
Chemistry.dend <-
  as.dendrogram(hclust(CHEMISTRY.DIST, method = "complete"))
ENZYMES.dend <- as.dendrogram(hclust(ENZYMES.DIST, method = "complete"))



#CREATE COLOR LIST
colors_to_use <- as.numeric(CLASS$ALLCLASS)

CHEMCOLORS <- colors_to_use[order.dendrogram(Chemistry.dend)]
ENZYCOLORS <- colors_to_use[order.dendrogram(ENZYMES.dend)]



#produce single Chemistry dendrogram
#color nodes by chemical class

pdf(paste("../Output_Plots/Global/","Global_Chemistry_dendrogram.pdf",sep="")) # writing a PDF to file
Chemistry.dend %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.4) %>%
  set("labels_cex", 0.25) %>%
  set("leaves_col", CHEMCOLORS) %>% plot
dev.off()
#produce single enzyme dendrogram
#color nodes by chemical class
pdf(paste("../Output_Plots/Global/","Global_Enzyme_dendrogram.pdf",sep="")) # writing a PDF to file
ENZYMES.dend %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.4) %>%
  set("labels_cex", 0.25) %>%
  set("leaves_col", ENZYCOLORS) %>% plot
dev.off()


#create a tangle gram 
# Extract labels from dendrogram on the left
labels <- ENZYMES.dend %>% set("labels_to_char") %>% labels
#Using a metadata table with colours create a vector of colours
labels <- as.data.frame(labels)
labels2 <-
  merge(labels,
        CLASS,
        by.x = "labels",
        by.y = "PUBchemID",
        sort = F)
labeler <- as.character(labels2$ALLCLASS)



#plot tangle gram
dend.chemenzy <- dendlist(ENZYMES.dend,Chemistry.dend)

# color legend
#black = 1 red = 2 green = 3 blue = 4 light blue = 5
#make tangle gram...remove labels apply node color based on group


pdf(paste("../Output_Plots/Global/","Global_Enzyme_Chemistry_Tanglegram.pdf",sep="")) # writing a PDF to file
dend.chemenzy %>% 
  #untangle(method = "random", R = 10) %>%
  untangle(method = "step2side") %>%
  tanglegram(
    highlight_distinct_edges = FALSE,
    highlight_branches_lwd =FALSE,
    common_subtrees_color_branches = F,
    color_lines = labeler,
    margin_inner = .5
  )

dev.off()




# #colors are recovered subtrees
# dend.chemenzy %>% untangle(method = "random", R = 10) %>%
#   untangle(method = "step1") %>%
#   tanglegram(common_subtrees_color_branches = T,
#              margin_inner = .5)





#permutation test of sig baker gamma
set.seed(23801)
the_cor <- cor_bakers_gamma(Chemistry.dend, Chemistry.dend)
the_cor2 <- cor_bakers_gamma(Chemistry.dend, ENZYMES.dend)
R <- 100
cor_bakers_gamma_results <- numeric(R)
dend_mixed <- Chemistry.dend
for (i in 1:R) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
  cor_bakers_gamma_results[i] <-
    cor_bakers_gamma(Chemistry.dend, dend_mixed)
}

pdf(paste("../Output_Plots/Global/","Global_Tanglegram_bakersgamma.pdf",sep="")) # writing a PDF to file

plot(density(cor_bakers_gamma_results),
     main = "Baker's gamma distribution under H0",
     xlim = c(-1, 1))
abline(v = 0, lty = 2)
abline(v = the_cor, lty = 2, col = 2)
abline(v = the_cor2, lty = 2, col = 4)
legend("topleft",cex = 0.75,
       legend = c("Identity", "Correlation"),
       fill = c(2, 4))
(sum(the_cor2 < cor_bakers_gamma_results) / R)
#p. value unnecessary if distributions do not touch
title(sub = paste(
  "One sided p-value:",
  "cor =",
  round(sum(the_cor < cor_bakers_gamma_results) / R, 4),
  " ; cor2 =",
  round(sum(the_cor2 < cor_bakers_gamma_results) / R, 4)
))

dev.off()




