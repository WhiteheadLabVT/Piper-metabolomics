################### PERMANOVA with pairwise tests using lmPerm and Tukey's HSD


####### convert chemical structural dissimilarity matrix to columnar pairwise format for lmPerm
### "csd" = dissimilarity matrix, "csd_col" = converted matrix

csd_col <- csd[lower.tri(csd)] 
who.vs.who <- expand.grid(rownames(csd), rownames(csd)) 
who <- who.vs.who[lower.tri(csd),] 
names(csd_col) <- paste(who[,1], who[,2], sep=".vs.")



library(lmPerm)

######## Example: comparison between tissues with species pooled. TIC data included.
mod.piper <- aovp(dissim ~ tissue*species, data = csd_col, maxIter = 10000, perm = "Prob")
summary(mod.piper)
TukeyHSD(mod.piper)
