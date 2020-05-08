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

#warning open babel must be installed first: http://openbabel.org/wiki/Category:Installation
#install & library require Bioconductor packages
requiredpackages <- c("ChemmineOB","ChemmineR")

for (pkg in requiredpackages) {
  if (pkg %in% rownames(installed.packages()) == FALSE)
  {BiocManager::install(pkg)}
  if (pkg %in% rownames(.packages()) == FALSE)
  {library(pkg, character.only=T)}
}


#SETWORKING DIRECTORY
setwd("/Data")

#IMPORT ENZYME DATA
ENZYMEDATABASE <- read.csv("Global_Compound_by_Enzyme_Table.csv")


#use online chemmine web service to get sdf from PUBchemIDs
#https://chemminetools.ucr.edu
# this can be accomplished internally through R but I have not added this yet 


#IMPORT SDF LIST
#rcdk list for getting descriptors
MOLS <- load.molecules("Global_Compound_by_Enzyme_Table.csv.sdf")
#chemminerR list if you want to visualize any compound
sdfset <- read.SDFset("Global_Compound_by_Enzyme_Table.csv.sdf")
#generate3D coordinates and add this to the sdf file
sdfset3d <- generate3DCoords(sdfset)

#write Sdf file
write.SDF(sdfset3d, file = "3D_Global_Compound_by_Enzyme_Table.sdf")
#rcdk list for getting descriptors
MOLS <- load.molecules("3D_Global_Compound_by_Enzyme_Table.sdf")
#gather all molecular descriptors
#create a vector of uniqure variables from all categorites of descriptors
descNames <-
  unique(unlist(sapply(get.desc.categories(), get.desc.names)))


compounddescriptors <-  eval.desc(MOLS, descNames)


#remove any variables that are all NA

compounddescriptors<-compounddescriptors[,colSums(is.na(compounddescriptors))<nrow(compounddescriptors)]



compounddescriptors <-
  data.frame(ENZYMEDATABASE[, 1], compounddescriptors[,-1])

#write csv of chemical descriptors
write.csv(compounddescriptors,
          "3D_Global_Compound_by_Enzyme_Table_DESCRIPTORS.csv")

