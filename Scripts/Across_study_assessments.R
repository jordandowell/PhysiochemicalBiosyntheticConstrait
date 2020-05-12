#assessment of adj R2ss
if(!require(coin)){install.packages("coin")}
if(!require(FSA)){install.packages("FSA")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(multcompView)){install.packages("multcompView")}
library(tidyverse)
library(ggsignif)
library(broom)
library(ggpubr)

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
  order = c("[E|C+E:C+C|E]","[E|C]","[E|C+E:C]","[E:C]","[E:C+C|E]","[C|E]", "Residuals" ),
  ylab = expression(paste("Adjusted R" ^ "2")),
  xlab = "Variance Components",
  add = c("boxplot","jitter"),
  add.params = list(fill = "white")
) +
  ylim(0,1) + 
  stat_summary( geom = "text", label = PVALUES$Letter,fun.y = max,vjust=-2,hjust=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()



####################################################################
#combinations

BIOTIC <- filter(dataset, dataset$Treatment %in% c("BIOTIC"))
HORMONE <- filter(dataset, dataset$Treatment %in% c("HORMONE"))
CONTROL <- filter(dataset, dataset$Treatment %in% c("CONTROL"))



View(CONTROL)



# your list of combinations you want to compare
CN <- combn(levels(dataset$RDA), 2, simplify = FALSE)
# the pvalues. I use broom and tidy to get a nice formatted dataframe. Note, I turned off the adjustment of the pvalues.
View(HORMONE)

pv <-
  tidy(with(
    dataset[dataset$RDA,],
    pairwise.wilcox.test(
      dataset$adjR2,
      dataset$RDA,
      paired = T,
      p.adjust.method = "bonferroni"
    )
  ))
View(pv)
#  data preparation
CN2 <- do.call(rbind.data.frame, CN)
colnames(CN2) <- colnames(pv)[-3]
# subset the pvalues, by merging the CN list
pv_final <-
  merge(CN2,
        pv,
        by.x = c("group2", "group1"),
        by.y = c("group1", "group2"))

# fix ordering
pv_final <- pv_final[order(pv_final$group1),]
# set signif level
pv_final$map_signif <-
  ifelse(pv_final$p.value > 0.05,
         "", "*")


# ifelse(pv_final$p.value > 0.01, "**", "*"))
View(pv_final)
# subset
gr <- pv_final$p.value <= 0.05
# the plot
ggviolin(
  dataset,
  x = "RDA",
  y = "adjR2",
  color = "RDA",
  order = c("A+B+C", "A", "A+B", "B", "B+C", "C"),
  ylab = expression(paste("Adjusted R" ^ "2")),
  xlab = "Explained Variance Component",
  add = "boxplot",
  add.params = list(fill = "white")
) +
  ylim(0,2.5) +
  geom_signif(
    comparisons = CN[gr],
    # y_position = seq(0.75,1.75,by=0.1),
    step_increase = 0.12,
    annotation = pv_final$map_signif[gr]
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

?ggviolin

CN[gr]





# subset
gr <- pv_final$p.value <= 0.05


ggviolin(
  dataset,
  x = "RDA",
  y = "AdjR2",
  color = "RDA",
  order = c("A+B+C", "A", "A+B", "B", "B+C", "C"),
  ylab = expression(paste("Adjusted R" ^ "2")),
  xlab = "Explained Variance Component",
  add = "boxplot",
  add.params = list(fill = "white")
) +
  stat_compare_means(
    comparisons = CN[gr],
    method = "wilcox.test",
    label = "p.signif",
    color = "red"
  )#+ # Add pairwise comparisons p-value
#           stat_compare_means(label.y = 50)     # Add global p-value

library(bestNormalize)
dataset2<- filter(dataset, dataset$RDA %in% c("B+C"))
kruskal.test(AdjR2 ~ Treatment, data = datatset2)


View(dataset2)
results <-
  pairwise.wilcox.test(dataset2$adjR2,
                       dataset2$Treatment,
                       paired = F,
                       p.adjust.method = "bonferroni")
results


newdata <- bestNormalize(dataset$adjR2)
newdata


dataset$x <-    newdata$x.t

kruskal.test(dataset$adjR2 ~ Treatment, data = dataset)
results <- pairwise.wilcox.test(dataset$adjR2,
                                dataset$Treatment,
                                paired = F,
                                p.adjust.method = "bonferroni")
results

View(dataset)
anova(dataset$adjR2 ~ Treatment, data = dataset)

ggdensity(dataset$x)
View(result$p.value)
? pairwise.wilcox.test
? kruskal.test
? ggviolin
? p.adjust

install.packages('devtools') #assuming it is not already installed

library(devtools)

install_github('andreacirilloac/updateR')

library(updateR)

updateR(admin_password = 'Unlv1991!')
mean(BIOTIC$adjR2)
mean(BIOTIC$adjR2)


library(coin)

independence_test(BIOTIC$adjR2 ~ BIOTIC$RDA, 
                  data = BIOTIC)

View(BIOTIC)
library(rcompanion)

PT = pairwisePermutationTest(adjR2 ~ RDA,
                             data = BIOTIC,
                             distribution = approximate(nresample = 10000),
                             method="bonferroni")

?pairwiseper
PT

?pairwisePermutationTest
