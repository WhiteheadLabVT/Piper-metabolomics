if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(sciplot)
library(ggplot2)
library(psych)
library(vegan)
library(PMCMR)
library(dunn.test)
library(agricolae)
library(lmPerm)
library(limma)
library(venneuler)




########################NMDS w/ Bray-Curtis
div.ticuni <- read.csv("allsp_diversity_ticuni_4tiss_v2.csv", header = F)
rnames <- read.csv("allsp_diversity_rownames_4tiss_v2.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div.ticuni) = rnames$row
colnames(div.ticuni) = cnames$col

div.ticuni.mdsexp <- read.csv("allsp_diversity_ticuni_4tiss_nmdsexp.csv", header = T)
row.names(div.ticuni.mdsexp) = rnames$row
div.ticuni.mds<-metaMDS(div.ticuni, distance = "bray", k = 6, trymax = 100, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
cor(vegdist(div.ticuni), dist(div.ticuni.mds$points))^2 
plot(div.ticuni.mds, choices = c(1, 2), type="n", ylim = c(-0.8, 0.8), xlim = c(-1.2,1.4)) #plots the ordination axes
points(div.ticuni.mds, display = "sites", pch=as.numeric(div.ticuni.mdsexp$Symbol),
       col=as.character(div.ticuni.mdsexp$Color), cex = 2)     

text(div.ticuni.mds, pos = 4, cex = 0.3, display = "sites")



div.css.allsamps <- read.csv("css_p12_ticuni_dissim_nohead.csv", header = F)
names.allsamps <- read.csv("cscs_matrix_headernames.csv", header = T)
row.names(div.css.allsamps) = names.allsamps$sample1
colnames(div.css.allsamps) = names.allsamps$sample1
div.css.mdsexp <- read.csv("css_piper12spp_ticuni_nmdsexp.csv", header = T)
row.names(div.css.mdsexp) = names.allsamps$sample1
div.css.allsamps.mds <-metaMDS(div.css.allsamps, k = 6, maxit = 10000, trymax = 1000)# sratmax = 0.9999999999, sfgrmin = 1e-10)
plot(div.css.allsamps.mds, choices = c(1, 2), type = "n") #xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4)) #plots the ordination axes
points(div.css.allsamps.mds, display = "sites", pch=as.numeric(div.css.mdsexp$Symbol),
       col=as.character(div.css.mdsexp$Color), cex = 2)




############# Bray-Curtis dissimilarity scores

bc.ticuni<-vegdist(div.ticuni, method="bray", binary=FALSE, na.rm = TRUE, upper = FALSE) 
bc.ticuni.matrix <- as.matrix(bc.ticuni)
dist.ticuni <- bc.ticuni.matrix[lower.tri(bc.ticuni.matrix)] 
who.vs.who <- expand.grid(rownames(bc.ticuni.matrix), rownames(bc.ticuni.matrix)) 
who <- who.vs.who[lower.tri(bc.ticuni.matrix),] 
names(dist.ticuni) <- paste(who[,1], who[,2], sep=".vs.") 
write.csv(dist.ticuni,"filename.csv")













