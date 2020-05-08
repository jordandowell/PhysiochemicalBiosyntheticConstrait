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



#SETWORKING DIRECTORY
setwd("/Data")


#import data with compounds as rows and predictors as columns
#IMPORT ENZYME DATA
ENZYMEDATABASE <- read.csv("Data/Global_Compound_by_Enzyme_Table.csv")

ENZYMES <- ENZYMEDATABASE
ENZYMES <- ENZYMES[!duplicated(ENZYMES$PUBCHEMID), ]

#total 137 compounds
#convert rownames to pubchem ID
ENZYMES <- data.frame(ENZYMES[, -1], row.names = ENZYMES[, 1])
#chemistry
CHEMISTRY <- read.csv("Data/3D_Global_Compound_by_Enzyme_Table_DESCRIPTORS.csv")
CHEMISTRY <- CHEMISTRY[,-1]
CHEMISTRY <- CHEMISTRY[!duplicated(CHEMISTRY), ]

#convert rownames to pubchem ID
CHEMISTRY <- data.frame(CHEMISTRY[, -1], row.names = CHEMISTRY[, 1])

# inorder to incorporate negative values into the bray curtis 
#each column had the absolute value of the minimum value per column added to each cell to presereve inter intra compound variation

i<-1
for (i in i:dim(CHEMISTRY)[2]) {
  columnminimum<-min(CHEMISTRY[,i])
  CHEMISTRY[,i]<- CHEMISTRY[,i] + abs(columnminimum)
}



#chemical class
CLASS <- read.csv("Data/Compound_Class.2.8.19.csv")

CLASS <- CLASS[!duplicated(CLASS$PUBchemID), ]
#convert rownames to pubchem ID
CLASS <- data.frame(CLASS, row.names = CLASS[, 1])



#full tanglegram

#create dendro grams
#gowers distance used for mixed data types
#For each variable, a particular distance metric that works
#well for that data type and is used to scale between 0-1
#Then a linear combination of those user specied weights
#(most simply an average) is calculated to create the final distance matrix
#for quantitative data = range normalzed Manhattan distance
#ordinal = variable is first ranked then Manhattan with adjustment for ties
#nominal = variables of k categories are first converted into k binary columns and then a Dice coefficient is used

Chemistry.dend <-
  as.dendrogram(hclust(vegdist(CHEMISTRY, method = "bray"), method = "complete"))
ENZYMES.dend <- (as.dendrogram(hclust(vegdist(ENZYMES,
                                              method = "bray"), method = "complete")))



#CREATE COLOR LIST
colors_to_use <- as.numeric(CLASS$ALLCLASS)

CHEMCOLORS <- colors_to_use[order.dendrogram(Chemistry.dend)]
ENZYCOLORS <- colors_to_use[order.dendrogram(ENZYMES.dend)]


#mantel test betoween chemisty and enzyme dataset
mantel(vegdist(CHEMISTRY, method = "bray"),
       vegdist(ENZYMES,
               method = "bray"))


#produce single Chemistry dendrogram
#color nodes by chemical class
Chemistry.dend %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.7) %>%
  set("labels_cex", 0.5) %>%
  set("leaves_col", CHEMCOLORS) %>% plot

#produce single enzyme dendrogram
#color nodes by chemical class
ENZYMES.dend %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.7) %>%
  set("labels_cex", 0.5) %>%
  set("leaves_col", ENZYCOLORS) %>% plot

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







Bk_plot(Chemistry.dend, ENZYMES.dend, main = "CORRECT Bk plot \n(based on dendrograms)")






#start

#STUDYSPECIFIC
#setworking directory for study
setwd("~/Documents/UCF DIssertation/QUALS/QUALS DATA/STUDY016")

#import data with samples as rows and compounds as predictors
ALLSAMPLES <-
  read.csv(
    "barbosa2016_C.csv",
    header = T,
    row.names = 1,
    check.names = F
  )
View(ALLSAMPLES)
#segment samples of interest
SAMPLES <- ALLSAMPLES[]

#remove invariant columns
SAMPLES <- SAMPLES[,
                   sapply(SAMPLES,
                          function(col)
                            length(unique(col))) > 1]


#SAMPLES
nrow(SAMPLES)
#compounds
ncol(SAMPLES)
#create a numeric compound list
COMPOUNDLIST <- as.list(colnames(SAMPLES))

#prune enzyme and chemistry set for compounds found in study
STUDY.ENZYMES <-
  subset(ENZYMES, rownames(ENZYMES) %in% COMPOUNDLIST)
STUDY.CHEMISTRY <-
  subset(CHEMISTRY, rownames(CHEMISTRY) %in% COMPOUNDLIST)
#remove any columns that are invariant
STUDY.ENZYMES <- STUDY.ENZYMES[,
                               sapply(STUDY.ENZYMES,
                                      function(col)
                                        length(unique(col))) > 1]
#remove any columns that are invariant
STUDY.CHEMISTRY <- STUDY.CHEMISTRY[,
                                   sapply(STUDY.CHEMISTRY,
                                          function(col)
                                            length(unique(col))) > 1]


dim(STUDY.CHEMISTRY)

#create 3 distance matricies
#gower for enzymes* binary data reduces to a sorenson index
ENZYMES.DIST <- vegdist(STUDY.ENZYMES, method = "gower")

#gower for chemistry
CHEMISTRY.DIST <- vegdist(STUDY.CHEMISTRY, method = "gower")



#Sample 1-Correlation matrix
SAMPLES.DIST <- as.dist(1 - cor(as.matrix((SAMPLES))))
#covariance because we are interested in flux which has scale variation
#are two trees significantly different

plot(CHEMISTRY.DIST, ENZYMES.DIST)
cor(CHEMISTRY.DIST, ENZYMES.DIST)


#mantel correlations
ENZYCHEM.mantel <-
  mantel.randtest(ENZYMES.DIST, CHEMISTRY.DIST, nrepet = 99999)
ENZYCHEM.mantel
pdf("016_C_mantel_c_e.pdf") # writing a PDF to file
plot(ENZYCHEM.mantel)                   # makes the actual plot
dev.off()

#if collinear cannont assess difference.

Chemistry.mantel <-
  mantel.randtest(SAMPLES.DIST, CHEMISTRY.DIST, nrepet = 99999)
Chemistry.mantel
pdf("016_C_mantel_c_s.pdf") # writing a PDF to file
plot(Chemistry.mantel)                  # makes the actual plot
dev.off()


ENZYMES.mantel <-
  mantel.randtest(SAMPLES.DIST, ENZYMES.DIST, nrepet = 99999)
ENZYMES.mantel

pdf("016_C_mantel_e_s.pdf") # writing a PDF to file
plot(ENZYMES.mantel)                 # makes the actual plot
dev.off()


#average distance in profile
CHEMISTRY.pairwise <- subset(melt(as.matrix(CHEMISTRY.DIST),
                                  varnames = c("row", "col")), value != 0)
#View(CHEMISTRY.pairwise)
mean(CHEMISTRY.pairwise[, 3])

ENZYMES.pairwise <- subset(melt(as.matrix(ENZYMES.DIST),
                                varnames = c("row", "col")), value != 0)
mean(ENZYMES.pairwise[, 3])
#variance partitioning


ENZYMECORR <- subset(melt((as.matrix(ENZYMES.DIST))))
SAMPLECORR <- subset(melt((as.matrix(SAMPLES.DIST))))
CHEMCORR <- subset(melt((as.matrix(CHEMISTRY.DIST))))


hello <- varpart(SAMPLECORR$value,  ~ ENZYMECORR$value,  ~ CHEMCORR$value)
(hello)
# fraction [a]:
rda.ENZYME.CHEM <-
  rda (SAMPLECORR$value ~ ENZYMECORR$value + Condition (CHEMCORR$value))
# fraction [c]:
rda.CHEM.ENZYME <-
  rda (SAMPLECORR$value ~ CHEMCORR$value + Condition (ENZYMECORR$value))
#fractions [a+b+c]:
rda.all <-
  rda (SAMPLECORR$value ~ ENZYMECORR$value + CHEMCORR$value)
# fractions [a+b]:
rda.ENZYME <- rda (SAMPLECORR$value ~ ENZYMECORR$value)
# fractions [b+c]:
rda.CHEM <- rda (SAMPLECORR$value ~ CHEMCORR$value)


#fractions [a+b+c]:
anova(rda.all)
# fraction [a]:
anova(rda.ENZYME.CHEM)
# fraction [c]:
anova(rda.CHEM.ENZYME)
# fractions [a+b]:
anova(rda.ENZYME)
# fractions [b+c]:
anova(rda.CHEM)


plot(
  hello,
  digits = 2,
  Xnames = c('E', 'T'),
  bg = c('navy', 'tomato'),
  id.size = 3
)


pdf("016_C_varpart.pdf") # writing a PDF to file
plot(
  hello,
  digits = 2,
  Xnames = c('E', 'T'),
  bg = c('navy', 'tomato'),
  id.size = 3
)                    # makes the actual plot
dev.off()                     # closes the PDF file





####making tangle gram
#######################################################
######   Biosynthetically informed Distances   ########
#######################################################
###############   Source Functions   ##################
#######################################################
############################# Robert R. Junker ########
##############################edit Jordan A. Dowell####
library(vegan)
library(cluster)
library(GUniFrac)

BioSynDist <- function(biosynthesisInfo, SecMetComp) {
  clus_Comp <-
    hclust(vegdist(biosynthesisInfo, method = "gower"), method = "average")
  unifracs <-
    GUniFrac(SecMetComp, as.phylo(clus_Comp), alpha = c(0, 0.5, 1))$unifracs
  BSD <- unifracs[, , "d_1"] # Weighted UniFrac
  #BioSynDist <- unifracs[, , "d_UW"] # Unweighted UniFrac
  #BSD <- unifracs[, , "d_VAW"] # Variance adjusted weighted UniFrac
  #BioSynDist <- unifracs[, , "d_0"]      # GUniFrac with alpha 0
  #BioSynDist <- unifracs[, , "d_0.5"]    # GUniFrac with alpha 0.5
  
  r <- list("BSDist" = BSD,
            "dend_Comp" = clus_Comp)
  
  #return(r)
  invisible(r)
  
}




#Calculate biosynthetically informed distances d(A,B) between samples

#CHEMICALLY INFORMED
CHEMr <-
  (BioSynDist(as.data.frame(STUDY.CHEMISTRY), as.data.frame(SAMPLES)))
CSDist <- as.dist(CHEMr$BSDist)

#ENZYME INFORMED
ENZYr <- BioSynDist(STUDY.ENZYMES, SAMPLES)
ESDist <- as.dist(ENZYr$BSDist)
#mantel of new distance metrics
mantel.randtest(CSDist, ESDist)
pdf("016_C_mantel_informed.pdf") 
plot(mantel.randtest(CSDist, ESDist))
dev.off()
Chemistry.informed.dend <-
  as.dendrogram(hclust(CSDist, method = "average"))
ENZYMES.informed.dend <-
  as.dendrogram(hclust(ESDist, method = "average"))


dend.chemenzy <- dendlist(as.dendrogram(Chemistry.informed.dend),
                          as.dendrogram(ENZYMES.informed.dend))

#color  tangle gram by treatments.
#import metadata
#METADATA
META <- read.csv("META_016_C.csv")
View(META)
#CLASS <- CLASS[!duplicated(CLASS$PUBchemID),]
#convert rownames to SAMPLE ID
META <- data.frame(META, row.names = META[, 1])
# Extract labels from dendrogram on the left
labels <-
  Chemistry.informed.dend %>% set("labels_to_char") %>% labels
#Using a metadata table with colours create a vector of colours
labels <- as.data.frame(labels)
labels2 <-
  merge(labels,
        META,
        by.x = "labels",
        by.y = "Sample",
        sort = F)
#View(labels2)
#CHANGE BASED ON TREATMENT CHOICE
METALABEL <- as.numeric(labels2[,2])
View(labels2)
#make tangle gram...remove labels apply node color based on group
set.seed(23801)
pdf("016_C_tanglegram_META_1.pdf") 
dend.chemenzy %>% untangle(method = "random", R = 10) %>%
  untangle(method = "step1") %>%
  tanglegram(
    common_subtrees_color_branches = F,
    color_lines = METALABEL,
    margin_inner = .5
  )
dev.off()
pdf("016_C_Bkplot_META_1.pdf")
Bk_plot(Chemistry.informed.dend, ENZYMES.informed.dend, main = "Bk plot \n(based on dendrograms)")
dev.off()
# #permutation test of sig baker gamma
# set.seed(23801)
# the_cor <- cor_bakers_gamma(Chemistry.informed.dend, ENZYMES.informed.dend)
# the_cor2 <- cor_bakers_gamma(Chemistry.informed.dend, ENZYMES.informed.dend)
# R <- 1000
# cor_bakers_gamma_results <- numeric(R)
# dend_mixed <- as.dendrogram(Chemistry.informed.dend)
# is.dendrogram(dend_mixed)
# for (i in 1:R) {
#   dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
#   cor_bakers_gamma_results[i] <-
#     cor_bakers_gamma(Chemistry.informed.dend, ENZYMES.informed.dend)
# }
# plot(density(cor_bakers_gamma_results),
#      main = "Baker's gamma distribution under H0",
#      xlim = c(-1, 1))
# abline(v = 0, lty = 2)
# abline(v = the_cor, lty = 2, col = 2)
# abline(v = the_cor2, lty = 2, col = 4)
# legend("topleft",
#        legend = c("cor", "cor2"),
#        fill = c(2, 4))
# (sum(the_cor2 < cor_bakers_gamma_results) / R)
# #p. value unnecessary if distributions do not touch
# title(sub = paste(
#   "One sided p-value:",
#   "cor =",
#   round(sum(the_cor < cor_bakers_gamma_results) / R, 4),
#   " ; cor2 =",
#   round(sum(the_cor2 < cor_bakers_gamma_results) / R, 4)
# ))


Bk_plot(Chemistry.informed.dend, ENZYMES.informed.dend, main = "Bk plot \n(based on dendrograms)")



