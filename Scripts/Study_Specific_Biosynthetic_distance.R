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
#start
#Get a list of study directories
StudyPaths<-list.dirs(path = ".", full.names = F, recursive = FALSE)

#create dataframes for outputs
GlobalMantel<-data.frame()
GlobalVariantionpartitioning<-data.frame()
GlobalCorrelation<-data.frame()
i<-1
j<-1
for( i in 1:length(StudyPaths)){
  print(StudyPaths[i])

  
  #get file names 
  Studyfiles<-list.files(path =StudyPaths[i], pattern = "Sampl" )

  
  for (j in 1:length(Studyfiles)) {
    StudyName<-tools::file_path_sans_ext(Studyfiles[j])
    
    #create output directory 
    dir.create(paste("../Output_Tables/",StudyName,sep=""))
    dir.create(paste("../Output_Plots/",StudyName,sep=""))
    #save output directory for later 
    currentTableDirectory<- paste("../Output_Tables/",StudyName,sep="")
    currentPlotDirectory<-paste("../Output_Plots/",StudyName,sep="")
    
    #import data with samples as rows and compounds as predictors
    ALLSAMPLES <-
      read.csv(paste(StudyPaths[i],"/",Studyfiles[j],sep = ""),
        header = T,
        row.names = 1,
        check.names = F
      )
    
    #segment samples of interest
    SAMPLES <- ALLSAMPLES[]
    
    #remove invariant columns
    SAMPLES <- SAMPLES[,
                       sapply(SAMPLES,
                              function(col)
                                length(unique(col))) > 1]
    
    
    #SAMPLES
    numberofsamples<-nrow(SAMPLES)
    #compounds
    numberofcompounds<-ncol(SAMPLES)
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
    
    
    #create 3 distance matricies
    #bray for enzymes
    ENZYMES.DIST <- vegdist(STUDY.ENZYMES, method = "bray")
    
    #bray for chemistry
    CHEMISTRY.DIST <- vegdist(STUDY.CHEMISTRY, method = "bray")
    
    
    
    #Sample 1-Correlation matrix
    SAMPLES.DIST <- as.dist(1 - cor(as.matrix((SAMPLES))))
    
    
    #create dataframe to store outputs 
    MantelOutput<-data.frame(Dataset=0,Comparisons=0,MantelR=0, SimPval=0, HypothesisTest=0, Std.obs=0, Expectation=0, Variance=0)
    
    #mantel correlations
    ENZYCHEM.mantel <-
      mantel.randtest(ENZYMES.DIST, CHEMISTRY.DIST, nrepet = 9999)
    
    MantelOutput[1,]<-c(StudyName,"EC",ENZYCHEM.mantel$obs,ENZYCHEM.mantel$pvalue ,ENZYCHEM.mantel$alter,ENZYCHEM.mantel$expvar[1],ENZYCHEM.mantel$expvar[2],ENZYCHEM.mantel$expvar[3])
    
    
    pdf(paste(currentPlotDirectory,"/_Mantel_Chemisty_Enzyme.pdf", sep="")) # writing a PDF to file
    plot(ENZYCHEM.mantel)                   # makes the actual plot
    dev.off()
    
    
    
    
    #if collinear cannont assess difference.
    
    Chemistry.mantel <-
      mantel.randtest(SAMPLES.DIST, CHEMISTRY.DIST, nrepet = 9999)
    Chemistry.mantel
    
    MantelOutput[2,]<-c(StudyName,"CS",Chemistry.mantel$obs,Chemistry.mantel$pvalue ,Chemistry.mantel$alter,Chemistry.mantel$expvar[1],Chemistry.mantel$expvar[2],Chemistry.mantel$expvar[3])
    
    pdf(paste(currentPlotDirectory,"/_Mantel_Chemistry_Samples.pdf", sep="")) # writing a PDF to file
    plot(Chemistry.mantel)                  # makes the actual plot
    dev.off()
    
    
    ENZYMES.mantel <-
      mantel.randtest(SAMPLES.DIST, ENZYMES.DIST, nrepet = 9999)
    
    MantelOutput[3,]<-c(StudyName,"ES",ENZYMES.mantel$obs,ENZYMES.mantel$pvalue ,ENZYMES.mantel$alter,ENZYMES.mantel$expvar[1],ENZYMES.mantel$expvar[2],ENZYMES.mantel$expvar[3])
    
    pdf(paste(currentPlotDirectory,"/_Mantel_Enzyme_Samples.pdf", sep="")) # writing a PDF to file
    plot(ENZYMES.mantel)                 # makes the actual plot
    dev.off()
    
    #bind to Global output dataset
    
    GlobalMantel<-rbind(GlobalMantel,MantelOutput)
    
    #average distance in profile
    CHEMISTRY.pairwise <- subset(melt(as.matrix(CHEMISTRY.DIST),
                                      varnames = c("row", "col")), value != 0)
    
    
    ENZYMES.pairwise <- subset(melt(as.matrix(ENZYMES.DIST),
                                    varnames = c("row", "col")), value != 0)
    
    #variance partitioning
    Meandistances<-data.frame(Study=StudyName, ChemistryDist=mean(CHEMISTRY.pairwise[, 3]),Chemistry.SD=sd(CHEMISTRY.pairwise[, 3]),Chemistry.SE=(sd(CHEMISTRY.pairwise[, 3]) / sqrt(length(CHEMISTRY.pairwise[, 3]))),
                              EnzymeDist=mean(ENZYMES.pairwise[, 3]),
                              ENZYME.SD=sd(ENZYMES.pairwise[, 3]),
                              ENZYME.SE=(sd(ENZYMES.pairwise[, 3]) / sqrt(length(ENZYMES.pairwise[, 3]))))               
    
    GlobalCorrelation<-rbind(GlobalCorrelation,Meandistances)
    
    
    
    
    #variance partitioning
    
    ENZYMECORR <- subset(melt((as.matrix(ENZYMES.DIST))))
    SAMPLECORR <- subset(melt((as.matrix(SAMPLES.DIST))))
    CHEMCORR <- subset(melt((as.matrix(CHEMISTRY.DIST))))
    
    
    varpartfullmodel <- varpart(SAMPLECORR$value,  ~ ENZYMECORR$value,  ~ CHEMCORR$value)
  
    
    #create data frame of variance partiton componenets
    Varpartition<-varpartfullmodel$part$fract
    Varpartition<-rbind(Varpartition,varpartfullmodel$part$indfract)
    Varpartition[,(length(Varpartition)+1)]<- row.names(Varpartition)
    colnames(Varpartition)[5] <- "Component"
 
    
    # fraction [E|C]:
    rda.ENZYME.CHEM <-
      rda (SAMPLECORR$value ~ ENZYMECORR$value + Condition (CHEMCORR$value))
    # fraction [C|E]:
    rda.CHEM.ENZYME <-
      rda (SAMPLECORR$value ~ CHEMCORR$value + Condition (ENZYMECORR$value))
    #fractions [E|C+E:C+C|E]:
    rda.all <-
      rda (SAMPLECORR$value ~ ENZYMECORR$value + CHEMCORR$value)
    # fractions [E|C+E:C]:
    rda.ENZYME <- rda (SAMPLECORR$value ~ ENZYMECORR$value)
    # fractions [E:C+C|E]:
    rda.CHEM <- rda (SAMPLECORR$value ~ CHEMCORR$value)
    
    
    #add column to store pvalue
    
    Varpartition$Pvalue<-"NA"
    Varpartition$STudy<- StudyName
    
    #add pvalue to all components 
    #assess significance via one-way permutation
    
    #fractions [a+b+c]:
    fraction.abc<-anova(rda.all, permutations = 10000)
    Varpartition[3,6]<-fraction.abc$`Pr(>F)`[1]
    # fraction [a]:
    fraction.a<-anova(rda.ENZYME.CHEM, permutations = 10000)
    Varpartition[4,6]<-fraction.a$`Pr(>F)`[1]
    
    # fraction [c]:
    fraction.c<-anova(rda.CHEM.ENZYME, permutations = 10000)
    Varpartition[6,6]<-fraction.c$`Pr(>F)`[1]
    
    # fractions [a+b]:
    fraction.ab<-anova(rda.ENZYME, permutations = 10000)
    Varpartition[1,6]<-fraction.ab$`Pr(>F)`[1]
    
    # fractions [b+c]:
    fraction.bc<-anova(rda.CHEM, permutations = 10000)
    Varpartition[2,6]<-fraction.bc$`Pr(>F)`[1]
    
    
    
    GlobalVariantionpartitioning<-rbind(GlobalVariantionpartitioning,Varpartition)
    
    pdf(paste(currentPlotDirectory,"/_VariancePartitioning.pdf", sep="")) # writing a PDF to file
    plot(
      varpartfullmodel,
      digits = 2,
      cutoff = 0.01,
      Xnames = c('E', 'C'),
      bg = c('navy', 'tomato'),
      id.size = 3
    )                    # makes the actual plot
    dev.off()                     # closes the PDF file
    
    
    ########################    
  }
 
  }



#write outputs to CSv

write.csv(GlobalMantel,"../Output_Tables/Global/Global_Mantel.csv")
write.csv(GlobalCorrelation,"../Output_Tables/Global/Global_Correlation.csv")


#rename partiitions before saving 
GlobalVariantionpartitioning$Component<-as.factor(GlobalVariantionpartitioning$Component)
levels(GlobalVariantionpartitioning$Component)[1:7]<- c("[E|C]","[E|C+E:C]","[E|C+E:C+C|E]","[E:C]","[E:C+C|E]","[C|E]", "Residuals" )




write.csv(GlobalVariantionpartitioning,"../Output_Tables/Global/Global_VariancePartition.csv")














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
    hclust(vegdist(biosynthesisInfo, method = "bray "), method = "average")
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


