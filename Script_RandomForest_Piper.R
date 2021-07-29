library(randomForest)
library(varSelRF)
library(VSURF)
library(Boruta)
library(psych)
library(dplyr)


piper <- read.csv("allsp_presabs_4tiss_boruta.csv", header =TRUE)
x1 <- piper[,c(3:1313)]
x1[is.na(x1)] <- 0
c1 <- as.factor(paste(piper$tissue))
piper_randfor <- randomForest(x1, c1, importance=TRUE, proximity=TRUE, oob.prox=TRUE, ntree=2000)
plot(piper_randfor)
importance(piper_randfor)
varImpPlot(piper_randfor)

piper_boruta <- Boruta(x1, c1)
plot(piper_boruta)
plotImpHistory(piper_boruta)
piper_boruta_table <- getImpRfZ(x1, c1)
boruta_sigfeatures <- getSelectedAttributes(piper_boruta)
write.csv(boruta_sigfeatures, "piper_boruta_sigfeats.csv")
write.csv(piper_boruta_table, "piper_boruta_table.csv")

boruta_classes <- read.csv("piper_boruta_table_classes_curated.csv", header = TRUE)
boruta_class_means <- aov(boruta_importance ~ chem_class, data = boruta_classes)
summary(boruta_class_means)
varImpPlot(piper_boruta)
TukeyHSD(boruta_class_means)

mean_by_class <- tapply(boruta_classes$boruta_importance, boruta_classes$chem_class, mean)
mean_by_class
sd_by_class <- tapply(boruta_classes$boruta_importance, boruta_classes$chem_class, sd)
sd_by_class

