#assessment of adj R2ss
if(!require(coin)){install.packages("coin")}
if(!require(FSA)){install.packages("FSA")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(multcompView)){install.packages("multcompView")}
library(tidyverse)
library(ggsignif)
library(broom)
library(ggpubr)
library(rcompanion)



####fulldata set Wilcox & violin plots


#read in data
GlobalVarPart<-read.csv("../Output_Tables/Global/Global_VariancePartition.csv")

View(GlobalVarPart)


# your list of combinations you want to compare
GlobalVarPart.Combinations <- combn(levels(GlobalVarPart$Component), 2, simplify = FALSE)
View(GlobalVarPart.Combinations)



PT<-pairwisePermutationTest(Adj.R.squared ~ Component,g =STudy,
                        data = GlobalVarPart,
                       # distribution = (nresample = 10000),
                        method="bonferroni")

PVALUES<-cldList(comparison = PT$Comparison,
        p.value    = PT$p.value,
        swap.colon = FALSE,
        threshold  = 0.05)


#order pvalue letters by proposed order of graphing data

target<-c("[E|C]","[E|C+E:C]","[E:C]","[E:C+C|E]","[C|E]","[E|C+E:C+C|E]", "Residuals" )


PVALUES<-PVALUES[match(target, PVALUES$Group),]



# for wilcoxnsigned rank test or pariwise T test 
# pv <-
#   tidy(with(
#     GlobalVarPart[GlobalVarPart$Component,],
#     pairwise.wilcox.test(
#       GlobalVarPart$Adj.R.squared,
#       GlobalVarPart$Component,
#       paired = T,
#       p.adjust.method = "bonferroni"
#     )
#   ))
# View(pv)
# 
# pv.matrix<- pairwise.wilcox.test(
#   GlobalVarPart$Adj.R.squared,
#   GlobalVarPart$Component,
#   paired = T,
#   p.adjust.method = "bonferroni"
# )
# View(round(pv.matrix$p.value, digits = 2))

# ?pairwise.wilcox.test
# #  data preparation
# GlobalVarPart.Combinations2 <- do.call(rbind.data.frame, GlobalVarPart.Combinations)
# colnames(GlobalVarPart.Combinations2) <- colnames(pv)[-3]
# # subset the pvalues, by merging the CN list
# pv_final <-
#   merge(GlobalVarPart.Combinations2,
#         pv,
#         by.x = c("group2", "group1"),
#         by.y = c("group1", "group2"))
# 
# # fix ordering
# pv_final <- pv_final[order(pv_final$group1),]
# # set signif level
# pv_final$map_signif <-
#   ifelse(pv_final$p.value > 0.05,
#          "", "*")
# 
# 
# # ifelse(pv_final$p.value > 0.01, "**", "*"))
# View(pv_final)
# # subset
# gr <- pv_final$p.value <= 0.05
# the plot

pdf("../Output_Plots/Global/Global_Variaiance_Partition.pdf")
ggviolin(
  GlobalVarPart,
  x = "Component",
  y = "Adj.R.squared",
  color = "Component",
  order = c("[E|C]","[E|C+E:C]","[E:C]","[E:C+C|E]","[C|E]","[E|C+E:C+C|E]", "Residuals" ),
  ylab = expression(paste("Adjusted R" ^ "2")),
  xlab = "Variance Components",
  add = c("boxplot","jitter"),
  add.params = list(fill = "white")
) +
  ylim(0,1) + 
  stat_summary( geom = "text", label = PVALUES$Letter,fun.y = max,vjust=-2,hjust=2)+
  theme(panel.grid.major = element_blank(),text = element_text(size=10), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()



####################################################################

#across treatment combinations 
#combinations


TreatmentVarPart<-read.csv("../Output_Tables/Global/Treatments_VariancePartition.csv")

Treatments<-c("BIOTIC","CHEMICAL","CONTROL")

i<-3

for (i in 1:length(Treatments)) {
  

VarPart_OI<- filter(TreatmentVarPart, TreatmentVarPart$TREATMENT %in% Treatments[i])

VarPart.Combinations <- combn(levels(VarPart_OI$Component), 2, simplify = FALSE)




PT<-pairwisePermutationTest(Adj.R.squared ~ Component,g =STudy,
                            data = VarPart_OI,
                            # distribution = (nresample = 10000),
                            method="bonferroni")

PVALUES<-cldList(comparison = PT$Comparison,
                 p.value    = PT$p.value,
                 swap.colon = FALSE,
                 threshold  = 0.05)


#order pvalue letters by proposed order of graphing data

target<-c("[E|C]","[E|C+E:C]","[E:C]","[E:C+C|E]","[C|E]","[E|C+E:C+C|E]", "Residuals" )


PVALUES<-PVALUES[match(target, PVALUES$Group),]






pdf(paste("../Output_Plots/Global/",Treatments[i],"_Variaiance_Partition.pdf",sep = ""))

ggviolin(
  VarPart_OI,
  x = "Component",
  y = "Adj.R.squared",
  color = "Component",
  order = c("[E|C]","[E|C+E:C]","[E:C]","[E:C+C|E]","[C|E]","[E|C+E:C+C|E]", "Residuals" ),
  ylab = expression(paste("Adjusted R" ^ "2")),
  xlab = "Variance Components",
  add = c("boxplot","jitter"),
  add.params = list(fill = "white")
) +
  ylim(0,1) + 
  stat_summary( geom = "text", label = PVALUES$Letter,fun.y = max,vjust=-2,hjust=2)+
  theme(panel.grid.major = element_blank(),text = element_text(size=10), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))




dev.off()

}
















