library(randomForest)
library(varSelRF)
library(VSURF)
library(Boruta)
library(psych)
library(dplyr)

#-------------------------------------------------------------
#  Machine Learning to explore organ-level affinities of compounds
#--------------------------------------------------

########## Random Forest analysis - can sample identity be predicted from chemical composition?

#### Remove "sp" and "PlantID" columns from Data_Peak_Table.csv
#### Convert ion intensities to presence/absence
#### Save as new file "Data_Peak_Table_randfor.csv"

piper <- read.csv("Data_Peak_Table_randfor.csv", header =TRUE) 
x1 <- piper[,c(3:1313)]
x1[is.na(x1)] <- 0
c1 <- as.factor(paste(piper$tissue))
piper_randfor <- randomForest(x1, c1, importance=TRUE, proximity=TRUE, oob.prox=TRUE, ntree=2000)
plot(piper_randfor)
importance(piper_randfor)
varImpPlot(piper_randfor)

########## Boruta analysis - identifying which compounds were individually informative in Random Forest

piper_boruta <- Boruta(x1, c1)
plot(piper_boruta)
plotImpHistory(piper_boruta)
piper_boruta_table <- getImpRfZ(x1, c1)
boruta_sigfeatures <- getSelectedAttributes(piper_boruta)
write.csv(boruta_sigfeatures, "piper_boruta_sigfeats.csv") 
#this gives a list of the compounds that had significantly higher variable importance scores than shadow variables
#for more information on these compounds, look them up in the molecular networking classification annotation 



