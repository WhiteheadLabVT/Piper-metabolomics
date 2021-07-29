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

######################################### abundance and diversity raw numbers
div.raw <- read.csv("all_pos_chemdiv_raw_4tiss.csv", header = T)
div.raw.3tiss <- read.csv("all_pos_chemdiv_raw_3tiss_nolf.csv", header = T)
hist(div.raw.3tiss$n_cpd_tiss_uniq)
#right-skewed, do ln transform
hist(div.raw.3tiss$ln_nctu)
#looks good
shapiro.test(div.raw$ln_nctu)
#passes
hist(div.raw$n_cpd_sp)
#good
shapiro.test(div.raw$n_cpd_sp)
#fails; attempts to transform failed
shapiro.test(div.raw.3tiss$ln_pct_uniq_tiss)

hist(div.raw.3tiss$n_cpd_tiss)
#good
shapiro.test(div.raw$n_cpd_tiss)
#passes
hist(div.raw$pct_cpd_tiss_persp)
#slight left skew
shapiro.test(div.raw$pct_cpd_tiss_persp)
#fails; not transformable

hist(div.raw$asqrt_pct_uniq_tiss)
#right-skew, try ln transform
hist(div.raw$ln_pct_uniq_tiss)
#good
shapiro.test(div.raw$ln_pct_uniq_tiss)
#passes

boxplot(n_cpd_tiss~tissue, data = div.raw.3tiss,ylim = c(0,700), xlab = "Piper Tissue", ylab = "Compound Richness")

mod.rich <- aov(ln_nctu~tissue, data = div.raw)
summary(mod.rich)
Df Sum Sq Mean Sq F value Pr(>F)
tissue       3  34690   11563    1.45  0.241
Residuals   44 350842    7974    

mod.rich.3tiss <- aov(n_cpd_tiss~tissue, data = div.raw.3tiss)
summary(mod.rich.3tiss)

mod.uniq.3tiss <- aov(ln_nctu~tissue, data = div.raw.3tiss)
#Df Sum Sq Mean Sq F value Pr(>F)  
#tissue       2  1.069  0.5344   3.374 0.0464 *
#  Residuals   33  5.228  0.1584 

TukeyHSD(mod.uniq.3tiss)
#diff         lwr       upr     p adj
#2_rf-1_uf -0.2666667 -0.66537253 0.1320392 0.2430174
#3_sd-1_uf  0.1500000 -0.24870586 0.5487059 0.6297380
#3_sd-2_rf  0.4166667  0.01796081 0.8153725 0.0390061
boxplot(n_cpd_tiss_uniq~tissue, data = div.raw.3tiss, ylim = c(0,80), xlab = "Piper Tissue", ylab = "# tissue-specific compounds per species", par(cex.lab=1.5), par(cex.axis=1.5))



mod.pctuniq.3tiss <- aov(ln_pct_uniq_tiss~tissue, data = div.raw.3tiss)
summary(mod.pctuniq.3tiss)
TukeyHSD(mod.pctuniq.3tiss)

boxplot(pct_uniq_cpd_tiss~tissue, data = div.raw.3tiss, ylim = c(0,20), xlab = "Piper Tissue", ylab = "% tissue-specific compounds per species", par(cex.lab=1.5), par(cex.axis=1.5))

kruskal.test(pct_cpd_tiss_persp ~ tissue, data = div.raw)
posthoc.kruskal.dunn.test(pct_cpd_tiss_persp ~ tissue, data = div.raw, p.adjust.method = "bonferroni")
boxplot(pct_cpd_tiss_persp ~ tissue, data = div.raw)

mod.uniq <- aov(n_cpd_tiss_uniq ~ tissue, data = div.raw)
summary(mod.uniq)
TukeyHSD(mod.uniq)
boxplot(pct_uniq_cpd_tiss ~ tissue, ylim = c(0,25), data = div.raw)

mod.uniq.frt <- aov(n_cpd_tiss_uniq ~ tissue, data = div.raw.3tiss)
summary(mod.uniq.frt)
TukeyHSD(mod.uniq.frt)
boxplot(n_cpd_tiss_uniq ~ tissue, ylim = c(0,40), data = div.raw.3tiss)


########################### venn diagram of compound richness
venn.comps <- read.csv("ncomps_venn.csv", header = T)

fvl_subfrt <- cbind(venn.comps$lf, venn.comps$frt_randsub)
venn.fvl_subfrt <- vennCounts(fvl_subfrt)
vennDiagram(venn.fvl_subfrt, names = c("Leaf", "Fruit"), circle.col = c("darkblue", "limegreen"), cex = c(3,1))

fvl_all <- cbind(venn.comps$lf, venn.comps$frt_all)
venn.fvl_all <- vennCounts(fvl_all)
vennDiagram(venn.fvl_all)

all4tiss <- cbind(venn.comps$lf, venn.comps$uf, venn.comps$rf, venn.comps$sd)
venn.4tiss <- vennCounts(all4tiss)
vennDiagram(venn.4tiss, names = c("Leaf", "Unripe pulp", "Ripe pulp", "Ripe seed"), circle.col = c("darkblue", "deepskyblue", "seagreen", "limegreen"), cex = c(2,1.5,1.3))

frt_5sppa <- cbind(venn.comps$adu_frt,venn.comps$aur_frt,venn.comps$bio_frt,venn.comps$col_frt,venn.comps$gen_frt)
frt_5sppb <- cbind(venn.comps$gla_frt,venn.comps$mul_frt,venn.comps$pel_frt,venn.comps$ret_frt,venn.comps$san_frt)
#venn.comps$sil_frt,venn.comps$umb_frt)
venn.frt_sppa <- vennCounts(frt_5sppa)
vennDiagram(venn.frt_sppa)
venn.frt_sppb <- vennCounts(frt_5sppb)
vennDiagram(venn.frt_sppb)

lf_5sppa <- cbind(venn.comps$adu_lf,venn.comps$aur_lf,venn.comps$bio_lf,venn.comps$col_lf,venn.comps$gen_lf)
lf_5sppb <- cbind(venn.comps$gla_lf,venn.comps$mul_lf,venn.comps$pel_lf,venn.comps$ret_lf,venn.comps$san_lf)
#venn.comps$sil_lf,venn.comps$umb_lf)
venn.lf_sppa <- vennCounts(lf_5sppa)
vennDiagram(venn.lf_sppa)
venn.lf_sppb <- vennCounts(lf_5sppb)
vennDiagram(venn.lf_sppb)

adu <- cbind(venn.comps$adu_lf, venn.comps$adu_frt)
venn.adu <- vennCounts(adu)
vennDiagram(venn.adu, names = c("adu_leaf", "adu_fruit"))

aur <- cbind(venn.comps$aur_lf, venn.comps$aur_frt)
venn.aur <- vennCounts(aur)
vennDiagram(venn.aur, names = c("aur_leaf", "aur_fruit"))

bio <- cbind(venn.comps$bio_lf, venn.comps$bio_frt)
venn.bio <- vennCounts(bio)
vennDiagram(venn.bio, names = c("bio_leaf", "bio_fruit"))

col <- cbind(venn.comps$col_lf, venn.comps$col_frt)
venn.col <- vennCounts(col)
vennDiagram(venn.col, names = c("col_leaf", "col_fruit"))
########################### fruit vs. leaf

div.raw.fvl <- read.csv("all_pos_chemdiv_raw_fvl_v3.csv", header = T)

hist(div.raw.fvl$n_cpd_tiss_uniq)
#right-skewed, ln transform insufficient, do rrt transform
hist(div.raw.fvl$rrt_nctu)
#looks good
shapiro.test(div.raw.fvl$rrt_nctu)
#passes
hist(div.raw.fvl$n_cpd_sp)
#good
shapiro.test(div.raw.fvl$n_cpd_sp)
#passes
shapiro.test(div.raw.fvl$rrt_nctu)

hist(div.raw.fvl$n_cpd_tiss)
#good
shapiro.test(div.raw.fvl$n_cpd_tiss)
#passes
hist(div.raw.fvl$cu_pct_cpd_persp)
#good
shapiro.test(div.raw.fvl$cu_pct_cpd_persp)
#passes

boxplot(pct_uniq_tiss_sptot~tissue, data = div.raw.fvl, ylim = c(0,20), xlab = "Piper Tissue", ylab = "% unique compounds in tissue/ species richness")

boxplot(n_cpd_tiss~tissue, data = div.raw.fvl, xlab = "Piper Tissue", ylab = "Compound Richness")
t.test(n_cpd_tiss~tissue, data =div.raw.fvl, paired=TRUE, conf.level=0.95)
#
data:  n_cpd_tiss by tissue
#t = 0.65089, df = 11, p-value = 0.5285
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -47.23290  86.89957
#sample estimates:
#  mean of the differences 
#19.83333

boxplot(n_cpd_tiss_uniq~tissue, data = div.raw.fvl, ylim = c(0,250), xlab = "Piper Tissue", ylab = "# tissue-specific compounds per species", par(cex.lab=1.5), par(cex.axis=1.5))
t.test(rrt_nctu~tissue, data =div.raw.fvl, paired=TRUE, conf.level=0.95)

t = -1.4565, df = 11, p-value = 0.1732
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.036230230  0.007374443
sample estimates:
  mean of the differences 
-0.01442789 
boxplot(pct_uniq_cpd_tiss~tissue, data = div.raw.fvl, ylim = c(0,50), xlab = "Piper Tissue", ylab = "% tissue-specific compounds")

t.test(rrt_pctu~tissue, data =div.raw.fvl, paired=TRUE, conf.level=0.95)
#NS
#t = -1.8293, df = 11, p-value = 0.09457
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.04772747  0.00440162
#sample estimates:
#  mean of the differences 
#-0.02166293 

boxplot(pct_cpd_tiss_persp~tissue, data = div.raw.fvl, xlab = "Piper Tissue", ylab = "% species richness in tissue")
wilcox.test(pct_cpd_tiss_persp ~ tissue, data = div.raw.fvl, paired=TRUE)


mod.uniq <- aov(ln_pct_uniq_tiss ~ tissue, data = div.raw)
summary(mod.uniq)
TukeyHSD(mod.uniq)
boxplot(pct_uniq_cpd_tiss ~ tissue, ylim = c(0,25), data = div.raw)



div_abd <- read.csv("all_pos_div_abd.csv", header = T)
boxplot(num_feats ~ sp_tissue, data = div_abd)
boxplot(num_feats ~ tissue, data = div_abd)
boxplot(num_feats ~ sp, data = div_abd)
boxplot(sum_TIC ~ sp_tissue, data = div_abd)
boxplot(sum_TIC ~ tissue, data = div_abd)
boxplot(sum_TIC ~ sp, data = div_abd)

div_abd <- read.csv("div_abd_summary_sptissue.csv", header = T)

describeBy(div_abd, div_abd$sp_tissue)

###################### randomly subsampled pooled fruit (uf+rf+sd) vs. leaves

frt_lf <- read.csv("div_abd_frt_vs_lf.csv", header = T)
#reciprocal root transformation -> passes shapiro test
shapiro.test(frt_lf$rec_frt)
shapiro.test(frt_lf$rec_lf)

t.test(frt_lf$rec_frt, 
       frt_lf$rec_lf, 
       paired=TRUE, 
       conf.level=0.95)
#Paired t-test

#data:  frt_lf$rec_frt and frt_lf$rec_lf
#t = -2.4696, df = 11, p-value = 0.03115
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.063348596 -0.003642689
#sample estimates:
#  mean of the differences 
#-0.03349564 

boxplot(uniq_cpd ~ tissue, data = frt_lf)

############################### diversity indexes
div <- read.csv("allsp_diversity_v2_100f.csv", header = F)
rnames <- read.csv("allsp_diversity_rownames_4tiss_v2.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div) = rnames$row
colnames(div) = cnames$col

div.fvl <- read.csv("allsp_diversity_100f_fvl.csv", header = F)
rnames.fvl <- read.csv("allsp_diversity_rownames_fvl_v2.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div.fvl) = rnames.fvl$row
colnames(div.fvl) = cnames$col

div.ticuni <- read.csv("allsp_diversity_ticuni_4tiss_v2.csv", header = F)
rnames <- read.csv("allsp_diversity_rownames_4tiss_v2.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div.ticuni) = rnames$row
colnames(div.ticuni) = cnames$col

div.ticuni.3t <- read.csv("allsp_diversity_ticuni_3tiss.csv", header = F)
rnames <- read.csv("allsp_diversity_rownames_3tiss.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div.ticuni.3t) = rnames$row
colnames(div.ticuni.3t) = cnames$col

div.fvl.ticuni <- read.csv("allsp_diversity_ticuni_fvl_v2.csv", header = F)
rnames.fvl <- read.csv("allsp_diversity_rownames_fvl_v2.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div.fvl.ticuni) = rnames.fvl$row
colnames(div.fvl.ticuni) = cnames$col


div.fvl.means <- read.csv("allsp_diversity_100f_fvl_means.csv", header = F)
rnames.fvl <- read.csv("allsp_diversity_rownames_fvl_means.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div.fvl.means) = rnames.fvl$row
colnames(div.fvl.means) = cnames$col

div.means <- read.csv("allsp_diversity_v2_100f_means.csv", header = F)
rnames <- read.csv("allsp_diversity_rownames_4tiss_means.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div.means) = rnames$row
colnames(div.means) = cnames$col


############################# Shannon H
H <- diversity(div.fvl)
write.csv(H, "allsp_div_H_100f_fvl.csv")

###parse sp, tissue, etc into columns as necessary then re-import

allsp_H <- read.csv("allsp_div_H_100f_v2.csv", header = T) # leaf, unripe fruit, ripe fruit, seed
allsp_H_fvl <- read.csv("allsp_div_H_v3.csv", header = T) # just bootstrapped fruit "samples" vs. leaves

H_spmeans <- aggregate(allsp_H[, 7], list(allsp_H$sp), mean)
write.csv(H_spmeans, "allsp_H_spmeans.csv")
H_tissmeans <- aggregate(allsp_H[, 7], list(allsp_H$tissue), mean)
write.csv(H_tissmeans, "allsp_H_tissmeans.csv")
H_groform_means <- aggregate(allsp_H[, 7], list(allsp_H$groform), mean)
write.csv(H_groform_means, "allsp_H_groform_means.csv")
H_hab_means <- aggregate(allsp_H[, 7], list(allsp_H$habitat), mean)
write.csv(H_hab_means, "allsp_H_habitat_means.csv")
H_tissmeans_fvl <- aggregate(allsp_H_fvl[, 7], list(allsp_H_fvl$tissue), mean)
write.csv(H_tissmeans_fvl, "allsp_H_tissmeans_fvl.csv")

######################## Simpson

simp <- diversity(div, index = "simpson")
write.csv(simp, "allsp_div_simp.csv")

###parse sp, tissue, etc into columns as necessary then re-import

allsp_simp <- read.csv("allsp_div_simp_v2.csv", header = T) # leaf, unripe fruit, ripe fruit, seed
allsp_simp_fvl <- read.csv("allsp_div_simp_v3.csv", header = T) # just bootstrapped fruit "samples" vs. leaves




########################NMDS w/ Bray-Curtis

div <- read.csv("allsp_diversity.csv", header = F)
rnames <- read.csv("allsp_diversity_rownames.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames.csv", header = T)
row.names(div) = rnames$row
colnames(div) = cnames$col
div.mds<-metaMDS(div, distance = "bray", k = 10, trymax = 100, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)

plot(div.mds, choices = c(1, 2), type="n", ylim = c(-0.8,0.8), xlim = c(-1, 1)) #plots the ordination axes
points(div.mds, display = "sites", pch = 21, col = "black", bg = c("dark green", "mediumblue", "yellow", "cyan"), cex = c(1.5))#displays both sites and species on the same plot.  Try choosing just "sites" to reduce clutter
text(div.mds, pos = 4, cex = 0.3, display = "sites")

div.fvl.mds <- metaMDS(div.fvl, distance = "bray", k = 10, trymax = 100, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
plot(div.fvl.mds, choices = c(1, 2), type="n", ylim = c(-0.7,1), xlim = c(-0.5, 0.5)) #plots the ordination axes
points(div.fvl.mds, display = "sites", pch = 21:21, col = "black", bg = c("green", "light blue"), cex = c(1.5)) #displays both sites and species on the same plot.  Try choosing just "sites" to reduce clutter
text(div.fvl.mds, pos = 4, cex = 0.3, display = "sites")

div.ticuni.mdsexp <- read.csv("allsp_diversity_ticuni_4tiss_nmdsexp.csv", header = T)
row.names(div.ticuni.mdsexp) = rnames$row
div.ticuni.mds<-metaMDS(div.ticuni, distance = "bray", k = 6, trymax = 100, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
cor(vegdist(div.ticuni), dist(div.ticuni.mds$points))^2 
plot(div.ticuni.mds, choices = c(1, 2), type="n", ylim = c(-0.8, 0.8), xlim = c(-1.2,1.4)) #plots the ordination axes
points(div.ticuni.mds, display = "sites", pch=as.numeric(div.ticuni.mdsexp$Symbol),
       col=as.character(div.ticuni.mdsexp$Color), cex = 2)     

text(div.ticuni.mds, pos = 4, cex = 0.3, display = "sites")

div.css <- read.csv("css_p12_means_matrix_mirrored_nohead.csv", header = F)
names <- read.csv("cscs_p12ticuni_means_matrix_headnames.csv", header = T)
row.names(div.css) = names$name
colnames(div.css) = names$name
div.css.mdsexp <- read.csv("css_p12_means_matrix_nmdsexp.csv", header = T)
row.names(div.css.mdsexp) = names$name
div.css.mds <-metaMDS(div.css, k = 10, maxit = 100000, trymax = 1000, sratmax = 0.9999999999, sfgrmin = 1e-10)
plot(div.css.mds, choices = c(1, 2), type = "n", ylim = c(-0.4,0.4), xlim = c(-0.1, 0.1)) #plots the ordination axes
points(div.css.mds, display = "sites", pch=as.numeric(div.css.mdsexp$Symbol),
       col=as.character(div.css.mdsexp$Color), cex = 3.5)

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
##### PCoA - not used
#PCOA.css <- pcoa(div.css, correction = "lingoes")
#biplot.pcoa(PCOA.css)


div.css.test <- read.csv("css_p12_means_matrix_mirrored_testmini.csv", header = F)
names.test <- read.csv("css_p12_means_matrix_testmini_headnames.csv", header = T)
row.names(div.css.test) = names.test$name
colnames(div.css.test) = names.test$name
div.css.mds.test <- metaMDS(div.css.test, k = 10, autotransform = FALSE, trymax = 100)

NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

NMDS.scree(div.css)

div.fvl.ticuni.mds<-metaMDS(div.fvl.ticuni, distance = "bray", k = 10, trymax = 100, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
plot(div.fvl.ticuni.mds, choices = c(1, 2), type="n", ylim = c(-0.5,0.5), xlim = c(-0.5, 0.7)) #plots the ordination axes
points(div.fvl.ticuni.mds, display = "sites", pch = 21:21, col = "black", bg = c("green", "light blue"), cex = c(1.5)) #displays both sites and species on the same plot.  Try choosing just "sites" to reduce clutter
text(div.fvl.ticuni.mds, pos = 4, cex = 0.3, display = "sites")

div.lf <- read.csv("allsp_diversity_v2_100f_lf.csv", header = F)
rnames <- read.csv("allsp_diversity_rownames_lf.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div.lf) = rnames$row
colnames(div.lf) = cnames$col
div.lf.mds<-metaMDS(div.lf, distance = "bray", k = 10, trymax = 100, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
plot(div.lf.mds, choices = c(1, 2), type="n", xlim = c(-0.2,0.2), ylim = c(-0.5,1)) #plots the ordination axes
points(div.lf.mds, display = "sites", pch = 21:21, col = "black", bg = "gold", cex = c(1.5)) #displays both sites and species on the same plot.  Try choosing just "sites" to reduce clutter
text(div.lf.mds, pos = 4, cex = 0.7, display = "sites")

div.uf <- read.csv("allsp_diversity_v2_100f_uf.csv", header = F)
rnames <- read.csv("allsp_diversity_rownames_uf.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div.uf) = rnames$row
colnames(div.uf) = cnames$col
div.uf.mds<-metaMDS(div.uf, distance = "bray", k = 8, trymax = 100, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
plot(div.uf.mds, choices = c(1, 2), type="n", xlim = c(-0.5,0.5), ylim = c(-1,1)) #plots the ordination axes
points(div.uf.mds, display = "sites", pch = 21:21, col = "black", bg = "gold", cex = c(1.5)) #displays both sites and species on the same plot.  Try choosing just "sites" to reduce clutter
text(div.uf.mds, pos = 4, cex = 0.7, display = "sites")


############# Bray-Curtis dissimilarity scores

bc<-vegdist(div.fvl, method="bray", binary=FALSE, na.rm = TRUE, upper = FALSE) 
bc.matrix <- as.matrix(bc)
dist <- bc.matrix[lower.tri(bc.matrix)] 
who.vs.who <- expand.grid(rownames(bc.matrix), rownames(bc.matrix)) 
who <- who.vs.who[lower.tri(bc.matrix),] 
names(dist) <- paste(who[,1], who[,2], sep=".vs.") 
write.csv(dist,"allsp_100f_braycurtis_fvl_v4.csv")

bc.ticuni<-vegdist(div.ticuni, method="bray", binary=FALSE, na.rm = TRUE, upper = FALSE) 
bc.ticuni.matrix <- as.matrix(bc.ticuni)
dist.ticuni <- bc.ticuni.matrix[lower.tri(bc.ticuni.matrix)] 
who.vs.who <- expand.grid(rownames(bc.ticuni.matrix), rownames(bc.ticuni.matrix)) 
who <- who.vs.who[lower.tri(bc.ticuni.matrix),] 
names(dist.ticuni) <- paste(who[,1], who[,2], sep=".vs.") 
write.csv(dist.ticuni,"allsp_braycurtis_4tiss_ticuni.csv")

bc.ticuni.fvl<-vegdist(div.fvl.ticuni, method="bray", binary=FALSE, na.rm = TRUE, upper = FALSE) 
bc.ticuni.fvl.matrix <- as.matrix(bc.ticuni.fvl)
dist.ticuni.fvl <- bc.ticuni.fvl.matrix[lower.tri(bc.ticuni.fvl.matrix)] 
who.vs.who <- expand.grid(rownames(bc.ticuni.fvl.matrix), rownames(bc.ticuni.fvl.matrix)) 
who <- who.vs.who[lower.tri(bc.ticuni.fvl.matrix),] 
names(dist.ticuni.fvl) <- paste(who[,1], who[,2], sep=".vs.") 
write.csv(dist.ticuni.fvl,"allsp_braycurtis_fvl_ticuni.csv")



########################## Bray-Curtis group comparisons

###between tissues, species pooled

bc_itbsp <- read.csv("allsp_100f_braycurtis_4tiss_itbsp.csv", header = TRUE)


mod.bc.itbsp <- aovp(bc_dissim ~ tiss1, data = bc_itbsp, maxIter = 10000, perm = "Prob")
summary(mod.bc.itbsp)
# Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
# tiss1          3   13.193    4.3978 10000 < 2.2e-16 ***
# Residuals   2552   60.571    0.0237  

TukeyHSD(mod.bc.itbsp)
#diff          lwr         upr     p adj
#rf-lf  0.04476344  0.022752361  0.06677452 0.0000011
#sd-lf -0.12364862 -0.145963313 -0.10133392 0.0000000
#uf-lf  0.06000304  0.037688343  0.08231774 0.0000000
#sd-rf -0.16841206 -0.190423140 -0.14640098 0.0000000
#uf-rf  0.01523960 -0.006771484  0.03725068 0.2832517
#uf-sd  0.18365166  0.161336959  0.20596635 0.0000000

########## all pairwise comparisons

bc_4tiss <- read.csv("allsp_100f_braycurtis_4tiss_v2.csv", header = T)
mod.bc_4tiss <- aovp(bc_dissim~comp, data = bc_4tiss, maxIter = 10000, perm = "Prob")
summary(mod.bc_4tiss)

TukeyHSD(mod.bc_4tiss)
# results in word doc

bargraph.CI(x.factor = comp, response = bc_dissim, data = bc_4tiss, err.width = 0.1,las = 2, ylim = c(0,1), ylab = "Bray-Curtis Dissimilarity")

bc_fvl <- read.csv("allsp_100f_braycurtis_fvl_v5.csv", header = T)
mod.bc_fvl<- aovp(bc_dissim~comp, data = bc_fvl, maxIter = 10000, perm = "Prob")
summary(mod.bc_fvl)
#               Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
# comp          16   16.102   1.00638 10000 < 2.2e-16 ***
# Residuals   2539   12.799   0.00504

TukeyHSD(mod.bc_fvl)
# results in word doc

bargraph.CI(x.factor = comp, response = bc_dissim, data = bc_fvl, err.width = 0.1,las = 2, ylim = c(0,1), ylab = "Bray-Curtis Dissimilarity")


bc_4tiss.ticuni <- read.csv("allsp_braycurtis_4tiss_ticuni_v2.csv", header = T)
mod.bc_4tiss.ticuni <- aovp(bc_dissim~comp, data = bc_4tiss.ticuni, maxIter = 10000, perm = "Prob")
summary(mod.bc_4tiss.ticuni)
#               Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
# comp           20   35.671   1.78354 10000 < 2.2e-16 ***
# Residuals   10419   55.083   0.00529    

TukeyHSD(mod.bc_4tiss.ticuni)
# results in word doc

bargraph.CI(x.factor = comp, response = bc_dissim, data = bc_4tiss.ticuni, err.width = 0.1,las = 2, ylim = c(0,0.6), ylab = "Bray-Curtis Dissimilarity")

####### 3 fruit tissues
bc_3tiss.ticuni <- read.csv("allsp_braycurtis_3tiss_ticuni.csv", header = T)
mod.bc_3tiss.ticuni <- aovp(bc_sim~comp, data = bc_3tiss.ticuni, maxIter = 10000, perm = "Prob")
summary(mod.bc_3tiss.ticuni)
#               Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
# comp          

TukeyHSD(mod.bc_3tiss.ticuni)
# results in word doc

bargraph.CI(x.factor = comp, response = bc_sim, data = bc_3tiss.ticuni, err.width = 0.1,las = 2, ylim = c(0.4,0.6), ylab = "Bray-Curtis Similarity")


######################## subsampled fruit vs. leaves
bc.ticuni.fvl <- read.csv("allsp_braycurtis_fvl_ticuni.csv", header = T)
mod.bc.fvl <- aovp(bc_dissim~comp, data = bc.ticuni.fvl, maxIter = 10000, perm = "Prob")
summary(mod.bc.fvl)
#              Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
#comp           1   2.1259   2.12594 10000 < 2.2e-16 ***
#Residuals   1186   4.3127   0.00364 

bargraph.CI(x.factor = comp, response = bc_sim, data = bc.ticuni.fvl, err.width = 0.1,las = 2, ylim = c(0,0.6), ylab = "Bray-Curtis Dissimilarity")

############### species-level comparisons

bc.ticuni.fvl.comps <- read.csv("allsp_braycurtis_fvl_ticuni_v3.csv", header = T)
mod.bc.fvl.comps <- aovp(bc_dissim~comp, data = bc.ticuni.fvl.comps, maxIter = 10000, perm = "Prob")
summary(mod.bc.fvl.comps)
#               Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
# comp          16   9.2106   0.57566 10000 < 2.2e-16 ***
# Residuals   2539   9.5144   0.00375 

TukeyHSD(mod.bc.fvl.comps)
# results in txt doc

bargraph.CI(x.factor = comp, response = bc_sim, data = bc.ticuni.fvl.comps, err.width = 0.1,las = 2, ylim = c(0.4,0.8), ylab = "Bray-Curtis similarity")

## just species vs. species, no group summaries 

bc.ticuni.fvl.comps <- read.csv("allsp_braycurtis_fvl_ticuni_spcomps.csv", header = T)
hist(bc.ticuni.fvl.comps$sq_bcsim)
shapiro.test(bc.ticuni.fvl.comps$logit_bcsim)
#data:  bc.ticuni.fvl.comps$logit_bcsim
#W = 0.95651, p-value = 0.001388
mod.bc.fvl.comps <- aovp(bc_dissim~comp, data = bc.ticuni.fvl.comps, maxIter = 10000, perm = "Prob")
summary(mod.bc.fvl.comps)
#               Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
# comp          16   9.2106   0.57566 10000 < 2.2e-16 ***
# Residuals   2539   9.5144   0.00375 

TukeyHSD(mod.bc.fvl.comps)
# results in txt doc

############### intra-species ripe pulp vs. leaf
bc_rfvl <- read.csv("allsp_braycurtis_ticuni_rfvl.csv", header = T)
mod.bc_rfvl <- aovp(bc_sim~comp, data = bc_rfvl, maxIter = 10000, perm = "Prob")
summary(mod.bc_rfvl)
#               Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
#comp  13   15.543   1.19561 10000 < 2.2e-16 ***
#Residuals   7472   38.416   0.00514 

TukeyHSD(mod.bc_rfvl)
# results in txt doc

bargraph.CI(x.factor = comp, response = bc_sim, data = bc_rfvl, err.width = 0.1,las = 2, ylim = c(0,1), ylab = "Bray-Curtis Similarity")

############### intra-species ripe seed vs. leaf
bc_sdvl <- read.csv("allsp_braycurtis_ticuni_sdvl.csv", header = T)
mod.bc_sdvl <- aovp(bc_sim~comp, data = bc_sdvl, maxIter = 10000, perm = "Prob")
summary(mod.bc_sdvl)
#               Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
#comp          13   12.965   0.99731 10000 < 2.2e-16 ***
#Residuals   7469   38.356   0.00514   

TukeyHSD(mod.bc_sdvl)
# results in txt doc

bargraph.CI(x.factor = comp, response = bc_sim, data = bc_sdvl, err.width = 0.1,las = 2, ylim = c(0,1), ylab = "Bray-Curtis Similarity")

############### intra-species unripe pulp vs. leaf
bc_ufvl <- read.csv("allsp_braycurtis_ticuni_ufvl.csv", header = T)
mod.bc_ufvl <- aovp(bc_sim~comp, data = bc_ufvl, maxIter = 10000, perm = "Prob")
summary(mod.bc_ufvl)
#               Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
# 

TukeyHSD(mod.bc_ufvl)
# results in word doc

bargraph.CI(x.factor = comp, response = bc_sim, data = bc_ufvl, err.width = 0.1,las = 2, ylim = c(0,1), ylab = "Bray-Curtis Similarity")


########################## Beta diversity 
div.sp.tiss <- read.csv("div_spp+tissue_groups_alltiss.csv", header = T)
z <- betadiver(div, "z")
mod.sp <- with(div.sp.tiss, betadisper(z, sp))
mod.sp

#No. of Positive Eigenvalues: 114
#No. of Negative Eigenvalues: 66

#Average distance to median:
#  adu    aur    bio    col    gen    gla    mul    pel    ret    san    sil    umb 
#0.3331 0.3249 0.3272 0.3583 0.3962 0.3001 0.3376 0.3793 0.3527 0.3354 0.3187 0.2506 

TukeyHSD(mod.sp)

boxplot(mod.sp)


z <- betadiver(div.ticuni, "z")
mod.ticuni.tiss <- with(div.sp.tiss, betadisper(z, tiss))
mod.ticuni.tiss

#No. of Positive Eigenvalues: 114
#No. of Negative Eigenvalues: 66

#Average distance to median:
#  frt     lf     rf     sd     uf 
#0.4088 0.3616 0.3996 0.4118 0.4023 

TukeyHSD(mod.ticuni.tiss)


w <- betadiver(div.ticuni, "w")
modw.ticuni.tiss <- with(div.sp.tiss, betadisper(w, tiss))
modw.ticuni.tiss


TukeyHSD(modw.ticuni.tiss)

boxplot(modw.ticuni.tiss, xlab = "Piper Tissue", ylab = "Beta Diversity", par(cex.lab=1.7), par(cex.axis=1.7))

simt <- betadiver(div.ticuni, "sim")

t <- betadiver(div.ticuni, "t")
modt.ticuni.tiss <- with(div.sp.tiss, betadisper(t, tiss))
modt.ticuni.tiss

TukeyHSD(modt.ticuni.tiss)

boxplot(modt.ticuni.tiss, xlab = "Piper Tissue", ylab = "Beta Diversity", par(cex.lab=1.7), par(cex.axis=1.7))

################## 3 fruit tissues
div.sp.3t <- read.csv("div_spp+tissue_groups_3t.csv", header = T)
w <- betadiver(div.ticuni.3t, "w")
modw.ticuni.3t <- with(div.sp.3t, betadisper(w, tiss))
modw.ticuni.3t

TukeyHSD(modw.ticuni.3t)

boxplot(modw.ticuni.3t, xlab = "Piper Tissue", ylab = "Beta Diversity (w)", par(cex.lab=1.7), par(cex.axis=1.7))

t <- betadiver(div.ticuni.3t, "t")
modt.ticuni.3t <- with(div.sp.3t, betadisper(t, tiss))
modt.ticuni.3t

wb <- betadisper(w, div.sp.3t$tiss)
anova(wb)

sorb <- betadisper(sor, div.sp.3t$tiss)
anova(sorb)


boxplot(modt.ticuni.3t, xlab = "Piper Tissue", ylab = "Beta Diversity (t)", par(cex.lab=1.7), par(cex.axis=1.7))


sor <- betadiver(div.ticuni.3t, "sor")
modsor.ticuni.3t <- with(div.sp.3t, betadisper(sor, tiss))
modsor.ticuni.3t

TukeyHSD(modsor.ticuni.3t)

boxplot(modsor.ticuni.3t, xlab = "Piper Tissue", ylab = "Beta Diversity (t)", par(cex.lab=1.7), par(cex.axis=1.7))

################################### fruit vs leaf
div.groups.fvl <- read.csv("div_spp+tissue_groups_fvl.csv", header = T)
z <- betadiver(div.fvl.ticuni, "z")
w <- betadiver(div.fvl.ticuni, "w")
mod.sp <- with(div.groups.fvl, betadisper(z, sp))
mod.sp

#No. of Positive Eigenvalues: 58
#No. of Negative Eigenvalues: 13

#Average distance to median:
#  adu    aur    bio    col    gen    gla    mul    pel    ret    san    sil    umb 
#0.3267 0.2834 0.3612 0.3131 0.3277 0.2687 0.3013 0.3438 0.3403 0.3392 0.3676 0.3108 

boxplot(mod.sp)


mod.tiss <- with(div.groups.fvl, betadisper(z, tiss))
mod.tiss


#No. of Positive Eigenvalues: 67
#No. of Negative Eigenvalues: 4

#Average distance to median:
#  frt     lf 
#0.4240 0.3705 

#Eigenvalues for PCoA axes:
#  (Showing 8 of 71 eigenvalues)
#PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
#1.7253 1.2708 1.0976 0.9435 0.7790 0.6676 0.5568 0.4818 

mod.tiss_w <- with(div.groups.fvl, betadisper(w, tiss))
mod.tiss_w

wb <- betadisper(w, div.groups.fvl$tiss)
anova(wb)
#Df   Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     1 0.055077 0.055077  45.204 3.973e-09 ***
#Residuals 70 0.085289 0.001218 


TukeyHSD(mod.tiss)
#           diff         lwr         upr      p adj
#lf-frt -0.0535295 -0.06940768 -0.03765133     0

TukeyHSD(mod.tiss_w)
#           diff         lwr         upr      p adj
#lf-frt -0.0535295 -0.06940768 -0.03765133     0

boxplot(mod.tiss, xlab = "Piper Tissue", ylab = "Beta Diversity", par(cex.lab=1.7), par(cex.axis=1.7))

sor <- betadiver(div.fvl.ticuni, "sor")
sorb <- betadisper(sor, div.groups.fvl$tiss)
anova(sorb)

wb <- betadisper(w, div.groups.fvl$tiss)
anova(wb)

################### means
div.groups.fvl <- read.csv("div_spp+tissue_groups_fvl_means.csv", header = T)
z <- betadiver(div.fvl.means, "z")
mod.sp <- with(div.groups.fvl, betadisper(z, sp))
mod.sp

mod.tiss <- with(div.groups.fvl, betadisper(z, tiss))
mod.tiss

TukeyHSD(mod.tiss)
#diff         lwr         upr     p adj
#lf-frt -0.03155748 -0.06503373 0.001918783 0.0634053

######################################## all tissues separate
div.groups.alltiss <- read.csv("div_spp+tissue_groups_alltiss.csv", header = T)
z <- betadiver(div.alltiss, "z")
w <- betadiver(div.alltiss, "w")
mod.sp <- with(div.groups.alltiss, betadisper(z, sp))
mod.sp

#No. of Positive Eigenvalues: 95
#No. of Negative Eigenvalues: 49

#Average distance to median:
#  adu    aur    bio    col    gen    gla    mul    pel    ret    san    sil    umb 
#0.3099 0.3084 0.3520 0.3801 0.3973 0.3101 0.3491 0.3884 0.3655 0.3568 0.3364 0.2688 

TukeyHSD(mod.sp)

boxplot(mod.sp)


mod.tiss <- with(div.groups.alltiss, betadisper(z, tiss))
mod.tiss

#No. of Positive Eigenvalues: 95
#No. of Negative Eigenvalues: 49

#Average distance to median:
#  lf     rf     sd     uf 
#0.3616 0.3996 0.4118 0.4023 

TukeyHSD(mod.tiss)
#diff          lwr        upr     p adj
#rf-lf  0.037990620  0.006064565 0.06991667 0.0126075
#sd-lf  0.050200241  0.018056258 0.08234422 0.0004657
#uf-lf  0.040660596  0.008516614 0.07280458 0.0068887
#sd-rf  0.012209621 -0.019716433 0.04413568 0.7529044
#uf-rf  0.002669977 -0.029256078 0.03459603 0.9963547
#uf-sd -0.009539644 -0.041683627 0.02260434 0.867149

boxplot(mod.tiss)
plot(mod.tiss)

################### means

div.groups.alltiss <- read.csv("div_spp+tissue_groups_4tiss_means.csv", header = T)
z <- betadiver(div.means, "z")
mod.sp <- with(div.groups.alltiss, betadisper(z, sp))
mod.sp

TukeyHSD(mod.sp)
#NS

mod.tiss <- with(div.groups.alltiss, betadisper(z, tiss))
mod.tiss

# No. of Positive Eigenvalues: 46
#No. of Negative Eigenvalues: 1

#Average distance to median:
#  lf     rf     sd     uf 
#0.3027 0.3326 0.3328 0.3327 

#Eigenvalues for PCoA axes:
#  (Showing 8 of 47 eigenvalues)
#PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
#0.9429 0.6506 0.5008 0.4808 0.3902 0.3091 0.3007 0.2488 

TukeyHSD(mod.tiss)
#NS

######################################################### Shannon H all tissues separate
allsp_H <- read.csv("allsp_div_H_100f_v2.csv", header = T) 
hist(allsp_H$Shannon_H)
#somewhat left-skewed, try squaring
hist(allsp_H$H_sq)
#better
shapiro.test(allsp_H$H_sq)
#passes
hist(allsp_H$evenness)
shapiro.test(allsp_H$evenness)
hist(allsp_H$even_sq)
shapiro.test(allsp_H$even_sq)
hist(allsp_H$even_cu)
shapiro.test(allsp_H$even_cu)
#passes
################################## comparing species
bargraph.CI(x.factor = sp, response = Shannon_H, data = allsp_H, err.width = 0.1,las = 1, ylim = c(0,5), xlab = "Piper Species", ylab = "Chem Diversity (Shannon H)")
mod.Hsp <- aov(H_sq ~ sp, data = allsp_H)
summary(mod.Hsp)
#              Df Sum Sq Mean Sq F value   Pr(>F)    
# sp           11  589.7   53.60   6.715 7.03e-09 ***
# Residuals   133 1061.8    7.98  

TukeyHSD(mod.Hsp)
#              diff         lwr        upr     p adj
aur-adu -2.3515589 -6.27560138  1.5724836 0.6967235
bio-adu -1.9365227 -5.86056516  1.9875198 0.8900991
col-adu -1.6023391 -5.45352379  2.2488457 0.9646441
gen-adu -4.4340205 -8.35806298 -0.5099780 0.0130250
gla-adu  1.3954015 -2.52864101  5.3194440 0.9894315
mul-adu -1.0952518 -5.01929428  2.8287907 0.9987045
pel-adu -4.6682327 -8.51941743 -0.8170480 0.0050588
ret-adu  1.1904419 -2.73360062  5.1144844 0.9972585
san-adu -1.3396552 -5.26369765  2.5843873 0.9924419
sil-adu -2.4779781 -6.40202057  1.4460644 0.6229391
umb-adu  1.6972234 -2.22681905  5.6212659 0.9534329
bio-aur  0.4150362 -3.42275304  4.2528255 0.9999999
col-aur  0.7492198 -3.01404216  4.5124818 0.9999499
gen-aur -2.0824616 -5.92025086  1.7553277 0.8122430
gla-aur  3.7469604 -0.09082889  7.5847496 0.0624924
mul-aur  1.2563071 -2.58148216  5.0940964 0.9946938
pel-aur -2.3166738 -6.07993580  1.4465882 0.6597944
ret-aur  3.5420008 -0.29578850  7.3797900 0.1005946
san-aur  1.0119037 -2.82588553  4.8496930 0.9992344
sil-aur -0.1264192 -3.96420845  3.7113701 1.0000000
umb-aur  4.0487823  0.21099307  7.8865716 0.0289915
col-bio  0.3341836 -3.42907838  4.0974456 1.0000000
gen-bio -2.4974978 -6.33528708  1.3402914 0.5774654
gla-bio  3.3319242 -0.50586511  7.1697134 0.1572067
mul-bio  0.8412709 -2.99651838  4.6790601 0.9998697
pel-bio -2.7317100 -6.49497201  1.0315520 0.4034715
ret-bio  3.1269645 -0.71082472  6.9647538 0.2327081
san-bio  0.5968675 -3.24092174  4.4346568 0.9999960
sil-bio -0.5414554 -4.37924467  3.2963339 0.9999985
umb-bio  3.6337461 -0.20404315  7.4715354 0.0816773
gen-col -2.8316814 -6.59494342  0.9315806 0.3470984
gla-col  2.9977405 -0.76552144  6.7610025 0.2632672
mul-col  0.5070873 -3.25617472  4.2703493 0.9999991
pel-col -3.0658936 -6.75312229  0.6213350 0.2068227
ret-col  2.7927809 -0.97048106  6.5560429 0.3685640
san-col  0.2626839 -3.50057808  4.0259459 1.0000000
sil-col -0.8756390 -4.63890101  2.8876230 0.9997669
umb-col  3.2995625 -0.46369949  7.0628245 0.1469902
gla-gen  5.8294220  1.99163271  9.6672112 0.0000875
mul-gen  3.3387687 -0.49902056  7.1765580 0.1550438
pel-gen -0.2342122 -3.99747420  3.5290498 1.0000000
ret-gen  5.6244624  1.78667309  9.4622516 0.0001868
san-gen  3.0943653 -0.74342393  6.9321546 0.2466648
sil-gen  1.9560424 -1.88174686  5.7938317 0.8670760
umb-gen  6.1312439  2.29345467  9.9690332 0.0000276
mul-gla -2.4906533 -6.32844253  1.3471360 0.5816909
pel-gla -6.0636342 -9.82689617 -2.3003722 0.0000225
ret-gla -0.2049596 -4.04274888  3.6328296 1.0000000
san-gla -2.7350566 -6.57284590  1.1027326 0.4328793
sil-gla -3.8733796 -7.71116883 -0.0355903 0.0457271
umb-gla  0.3018220 -3.53596731  4.1396112 1.0000000
pel-mul -3.5729809 -7.33624290  0.1902811 0.0798069
ret-mul  2.2856937 -1.55209561  6.1234829 0.7049011
san-mul -0.2444034 -4.08219263  3.5933859 1.0000000
sil-mul -1.3827263 -5.22051556  2.4550630 0.9882501
umb-mul  2.7924752 -1.04531404  6.6302645 0.3996693
ret-pel  5.8586746  2.09541258  9.6219366 0.0000505
san-pel  3.3285775 -0.43468444  7.0918395 0.1382618
sil-pel  2.1902546 -1.57300737  5.9535166 0.7342750
umb-pel  6.3654561  2.60219415 10.1287181 0.0000066
san-ret -2.5300970 -6.36788629  1.3076922 0.5573222
sil-ret -3.6684199 -7.50620921  0.1693693 0.0753400
umb-ret  0.5067816 -3.33100769  4.3445708 0.9999993
sil-san -1.1383229 -4.97611219  2.6994663 0.9977559
umb-san  3.0368786 -0.80091067  6.8746679 0.2725692
umb-sil  4.1752015  0.33741226  8.0129908 0.0205583

################################## comparing tissues
bargraph.CI(x.factor = tissue, response = Shannon_H, data = allsp_H, err.width = 0.1,las = 1, ylim = c(0,5), xlab = "Piper Tissue", ylab = "Chem Diversity (Shannon H)", cex = 1.5)
mod.Ht <- aov(H_sq~ tissue, data = allsp_H)
summary(mod.Ht)
#              Df Sum Sq Mean Sq F value   Pr(>F)    
# tissue        3  221.1   73.70   7.265 0.000144 ***
# Residuals   141 1430.3   10.14

TukeyHSD(mod.Ht)
#diff        lwr        upr     p adj
rf-lf  0.9016537 -1.0369015  2.8402089 0.6220089
sd-lf -2.4639154 -4.4157032 -0.5121276 0.0070363
uf-lf -0.3639685 -2.3157563  1.5878193 0.9623561
sd-rf -3.3655691 -5.3041243 -1.4270139 0.0000781
uf-rf -1.2656222 -3.2041774  0.6729330 0.3288900
uf-sd  2.0999469  0.1481591  4.0517347 0.0296287


######################################### evenness
bargraph.CI(x.factor = tissue, response = evenness, data = allsp_H, err.width = 0.1,las = 1, ylim = c(0,0.8), xlab = "Piper Tissue", ylab = "Shannon Evenness", cex = 1.5)
mod.et <- aov(even_cu~ tissue, data = allsp_H)
summary(mod.et)

# Df Sum Sq  Mean Sq F value   Pr(>F)    
# tissue        3 0.0919 0.030635   6.913 0.000224 ***
# Residuals   141 0.6249 0.004432   

TukeyHSD(mod.et)
#              diff          lwr          upr     p adj
#rf-lf  0.03038770 -0.01013065  0.0709060457 0.2122163
#sd-lf -0.03970302 -0.08049795  0.0010919063 0.0595921
#uf-lf -0.01105146 -0.05184638  0.0297434725 0.8952322
#sd-rf -0.07009072 -0.11060907 -0.0295723696 0.0000834
#uf-rf -0.04143915 -0.08195750 -0.0009208035 0.0429180
#uf-sd  0.02865157 -0.01214336  0.0694464938 0.2654878

################################# comparing ecological variables
bargraph.CI(x.factor = groform, response = Shannon_H, data = allsp_H, err.width = 0.1,las = 1, ylim = c(0,5), xlab = "Piper Growth Form", ylab = "Chem Diversity (Shannon H)")
mod.Hgroform <- aov(H_sq ~ groform, data = allsp_H)
summary(mod.Hgroform)
#              Df Sum Sq Mean Sq F value   Pr(>F)    
# groform       3  180.3   60.10   5.761 0.000957 ***
# Residuals   141 1471.1   10.43 

TukeyHSD(mod.Hgroform)
#                diff         lwr       upr     p adj
shrub-herb  3.6617108  1.09256807 6.2308535 0.0017021
tree-herb   4.0477457  1.42203588 6.6734555 0.0005672
vine-herb   2.8816178 -0.01037808 5.7736136 0.0512046
tree-shrub  0.3860349 -1.24022419 2.0122940 0.9264612
vine-shrub -0.7800930 -2.80839444 1.2482084 0.7496617
vine-tree  -1.1661279 -3.26561935 0.9333635 0.4741663


bargraph.CI(x.factor = habitat, response = Shannon_H, data = allsp_H, err.width = 0.1,las = 1, ylim = c(0,5), xlab = "Piper Successional Habitat", ylab = "Chem Diversity (Shannon H)", cex = 1.5)
mod.Hhab <- aov(H_sq ~ habitat, data = allsp_H)
summary(mod.Hhab)
#              Df Sum Sq Mean Sq F value Pr(>F)  
# habitat       2    102   50.98   4.672 0.0108 *
# Residuals   142   1550   10.91 

TukeyHSD(mod.Hhab)
#                 diff       lwr        upr     p adj
late-early  3.1418901  0.704731  5.5790493 0.0075698
mid-early   0.3632949 -1.000080  1.7266697 0.8032487
mid-late   -2.7785952 -5.252748 -0.3044426 0.0235395

mod.Hgrohab <- aov(H_sq~ habitat*groform, data = allsp_H)
summary(mod.Hgrohab)
# 

TukeyHSD(mod.Hgrohab)
$`habitat`


$groform


mod.Hght <- aov(H_sq ~ habitat*groform*tissue, data = allsp_H)
summary(mod.Hght)
     

TukeyHSD(mod.Hght)



############################################################# bootstrapped frt vs leaf
allsp_H_fvl <- read.csv("allsp_div_H_100f_fvl_v2.csv", header = T)
hist(allsp_H_fvl$Shannon_H)
#left-skewed, try squaring
hist(allsp_H_fvl$H_sq)
#better
shapiro.test(allsp_H_fvl$H_sq)
#passes
hist(allsp_H_fvl$evenness)
hist(allsp_H_fvl$even_sq)
shapiro.test(allsp_H_fvl$even_sq)
#passes
################################## comparing species
bargraph.CI(x.factor = sp, response = Shannon_H, data = allsp_H_fvl, err.width = 0.1,las = 1, ylim = c(0,5), xlab = "Piper Tissue", ylab = "Chem Diversity (Shannon H)")

mod.Hsp <- aov(H_sq ~ sp, data = allsp_H_fvl)
summary(mod.Hsp)
#             Df Sum Sq Mean Sq F value  Pr(>F)   
# sp          11  188.1  17.096   3.035 0.00278 **
# Residuals   60  337.9   5.632  

TukeyHSD(mod.Hsp)
#           diff        lwr         upr     p adj
aur-adu -2.11869490 -6.7774283  2.54003851 0.9208588
bio-adu -1.94259198 -6.6013254  2.71614144 0.9553345
col-adu -0.91539666 -5.5741301  3.74333675 0.9999356
gen-adu -2.95458031 -7.6133137  1.70415310 0.5861165
gla-adu  0.08510802 -4.5736254  4.74384143 1.0000000
mul-adu -0.15070908 -4.8094425  4.50802433 1.0000000
pel-adu -4.51777483 -9.1765082  0.14095859 0.0653886
ret-adu -0.01689437 -4.6756278  4.64183905 1.0000000
san-adu -0.35862854 -5.0173619  4.30010487 1.0000000
sil-adu -3.42577216 -8.0845056  1.23296125 0.3598033
umb-adu  0.94578810 -3.7129453  5.60452151 0.9999110
bio-aur  0.17610293 -4.4826305  4.83483634 1.0000000
col-aur  1.20329824 -3.4554352  5.86203165 0.9991117
gen-aur -0.83588541 -5.4946188  3.82284800 0.9999741
gla-aur  2.20380292 -2.4549305  6.86253633 0.8991941
mul-aur  1.96798582 -2.6907476  6.62671923 0.9511826
pel-aur -2.39907992 -7.0578133  2.25965349 0.8368416
ret-aur  2.10180054 -2.5569329  6.76053395 0.9247614
san-aur  1.76006637 -2.8986670  6.41879978 0.9780484
sil-aur -1.30707725 -5.9658107  3.35165616 0.9981197
umb-aur  3.06448300 -1.5942504  7.72321642 0.5310271
col-bio  1.02719532 -3.6315381  5.68592873 0.9998007
gen-bio -1.01198834 -5.6707217  3.64674507 0.9998275
gla-bio  2.02769999 -2.6310334  6.68643341 0.9403582
mul-bio  1.79188289 -2.8668505  6.45061630 0.9749155
pel-bio -2.57518285 -7.2339163  2.08355056 0.7666890
ret-bio  1.92569761 -2.7330358  6.58443102 0.9579527
san-bio  1.58396344 -3.0747700  6.24269685 0.9903463
sil-bio -1.48318018 -6.1419136  3.17555323 0.9943813
umb-bio  2.88838008 -1.7703533  7.54711349 0.6192061
gen-col -2.03918365 -6.6979171  2.61954976 0.9381009
gla-col  1.00050468 -3.6582287  5.65923809 0.9998456
mul-col  0.76468758 -3.8940458  5.42342099 0.9999895
pel-col -3.60237817 -8.2611116  1.05635525 0.2873453
ret-col  0.89850230 -3.7602311  5.55723571 0.9999465
san-col  0.55676812 -4.1019653  5.21550153 0.9999996
sil-col -2.51037550 -7.1691089  2.14835791 0.7938843
umb-col  1.86118476 -2.7975486  6.51991817 0.9669361
gla-gen  3.03968833 -1.6190451  7.69842174 0.5434293
mul-gen  2.80387123 -1.8548622  7.46260464 0.6608819
pel-gen -1.56319451 -6.2219279  3.09553890 0.9913257
ret-gen  2.93768595 -1.7210475  7.59641936 0.5945823
san-gen  2.59595178 -2.0627816  7.25468519 0.7576686
sil-gen -0.47119184 -5.1299253  4.18754157 0.9999999
umb-gen  3.90036842 -0.7583650  8.55910183 0.1874631
mul-gla -0.23581710 -4.8945505  4.42291631 1.0000000
pel-gla -4.60288285 -9.2616163  0.05585057 0.0556714
ret-gla -0.10200238 -4.7607358  4.55673103 1.0000000
san-gla -0.44373656 -5.1024700  4.21499685 1.0000000
sil-gla -3.51088018 -8.1696136  1.14785323 0.3237332
umb-gla  0.86068008 -3.7980533  5.51941349 0.9999652
pel-mul -4.36706574 -9.0257992  0.29166767 0.0861829
ret-mul  0.13381472 -4.5249187  4.79254813 1.0000000
san-mul -0.20791945 -4.8666529  4.45081396 1.0000000
sil-mul -3.27506307 -7.9337965  1.38367034 0.4282844
umb-mul  1.09649719 -3.5622362  5.75523060 0.9996276
ret-pel  4.50088046 -0.1578530  9.15961387 0.0674824
san-pel  4.15914629 -0.4995871  8.81787970 0.1237070
sil-pel  1.09200267 -3.5667307  5.75073608 0.9996418
umb-pel  5.46356293  0.8048295 10.12229634 0.0091958
san-ret -0.34173417 -5.0004676  4.31699924 1.0000000
sil-ret -3.40887779 -8.0676112  1.24985562 0.3672019
umb-ret  0.96268247 -3.6960509  5.62141588 0.9998941
sil-san -3.06714362 -7.7258770  1.59158979 0.5296983
umb-san  1.30441664 -3.3543168  5.96315005 0.9981535
umb-sil  4.37156026 -0.2871732  9.03029367 0.0854904


################################## comparing tissues
bargraph.CI(x.factor = tissue, response = Shannon_H, data = allsp_H_fvl, err.width = 0.1,las = 1, ylim = c(0,5), xlab = "Piper Tissue", ylab = "Chem Diversity (Shannon H)", cex = 1.3)
t.test(H_sq ~ tissue, data = allsp_H_fvl, paired = TRUE, conf.level = 0.95)

# t = -1.0799, df = 35, p-value = 0.2876
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.9201155  0.5866931
#sample estimates:
#  mean of the differences 
#-0.6667112 

####################### evenness
bargraph.CI(x.factor = tissue, response = evenness, data = allsp_H_fvl, err.width = 0.1,las = 1, ylim = c(0,1), xlab = "Piper Tissue", ylab = "Shannon Evenness")
t.test(even_sq ~ tissue, data = allsp_H_fvl, paired = TRUE, conf.level = 0.95)

#t = -0.7869, df = 35, p-value = 0.4366
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.03606903  0.01591800
#sample estimates:
#  mean of the differences 
#-0.01007552 

################################# comparing ecological variables
bargraph.CI(x.factor = groform, response = Shannon_H, data = allsp_H_fvl, err.width = 0.1,las = 1, ylim = c(0,5), xlab = "Piper Growth Form", ylab = "Chem Diversity (Shannon H)")
mod.Hgroform <- aov(H_sq ~ groform, data = allsp_H_fvl)
summary(mod.Hgroform)
#             Df Sum Sq Mean Sq F value Pr(>F)  
# groform      3   78.4  26.118   3.968 0.0115 *
# Residuals   68  447.6   6.583 

TukeyHSD(mod.Hgroform)
#                  diff        lwr      upr     p adj
shrub-herb  3.6375733  0.6156116 6.659535 0.0119626
tree-herb   3.7990541  0.7147773 6.883331 0.0096644
vine-herb   2.7295342 -0.6491217 6.108190 0.1547406
tree-shrub  0.1614808 -1.6890853 2.012047 0.9956824
vine-shrub -0.9080391 -3.2161005 1.400022 0.7289096
vine-tree  -1.0695199 -3.4585904 1.319551 0.6420505

mod.Hhab <- aov(H_sq ~ habitat, data = allsp_H_fvl)
summary(mod.Hhab)
#             Df Sum Sq Mean Sq F value Pr(>F)
# habitat      2   15.4   7.679   1.038   0.36
# Residuals   69  510.6   7.400 

bargraph.CI(x.factor = habitat, response = Shannon_H, data = allsp_H_fvl, err.width = 0.1,las = 1, ylim = c(0,5), xlab = "Piper Successional Habitat", ylab = "Chem Diversity (Shannon H)")





###########################################################Simpson diversity 
div <- read.csv("allsp_diversity_v2_100f.csv", header = F)
rnames <- read.csv("allsp_diversity_rownames_4tiss_v2.csv", header = T)
cnames <- read.csv("allsp_diversity_colnames_v2.csv", header = T)
row.names(div) = rnames$row
colnames(div) = cnames$col
simp <- diversity(div, index = "simpson")
write.csv(simp, "allsp_div_simp_100f.csv")

###parse sp, tissue, etc into columns as necessary then re-import

allsp_simp <- read.csv("allsp_div_simp_100f_v2.csv", header = T) #lf, uf, rf, sd

hist(allsp_simp$simpson)
#very left-skewed, power and arcsin-sqrt transformations didn't work

################################## comparing species
boxplot(simpson ~ sp, data = allsp_simp, ylim = c(0.9,1.0), xlab = "Piper Species", ylab = "Chem Diversity (Simpson)")
kruskal.test(simpson ~ sp, data = allsp_simp)
# ruskal-Wallis chi-squared = 30.081, df = 11, p-value = 0.001539

posthoc.kruskal.dunn.test(simpson ~ sp, data = allsp_simp, p.adjust.method = "bonferroni")

adu   aur   bio   col   gen   gla   mul   pel   ret   san   sil  
aur 1.000 -     -     -     -     -     -     -     -     -     -    
  bio 1.000 1.000 -     -     -     -     -     -     -     -     -    
  col 1.000 1.000 1.000 -     -     -     -     -     -     -     -    
  gen 1.000 1.000 0.824 1.000 -     -     -     -     -     -     -    
  gla 1.000 1.000 1.000 1.000 0.323 -     -     -     -     -     -    
  mul 1.000 1.000 1.000 1.000 0.150 1.000 -     -     -     -     -    
  pel 0.810 1.000 0.521 1.000 1.000 0.191 0.084 -     -     -     -    
  ret 1.000 1.000 1.000 1.000 0.036 1.000 1.000 0.019 -     -     -    
  san 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 -     -    
  sil 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 -    
  umb 1.000 1.000 1.000 1.000 0.859 1.000 1.000 0.545 1.000 1.000 1.000

################################## comparing tissues
boxplot(simpson ~ tissue, data = allsp_simp, ylim = c(0.9,1.0), xlab = "Piper Tissue", ylab = "Chem Diversity (Simpson)")
kruskal.test(simpson ~ tissue, data = allsp_simp)
# Kruskal-Wallis chi-squared = 29.938, df = 3, p-value = 1.422e-06

posthoc.kruskal.dunn.test(simpson ~ tissue, data = allsp_simp, p.adjust.method = "bonferroni")
# 
lf      rf      sd     
rf 0.42693 -       -      
sd 0.00391 9.8e-07 -      
uf 1.00000 1.00000 0.00045

################################# comparing ecological variables
boxplot(simpson ~ groform, data = allsp_simp, ylim = c(0.9,1.0), xlab = "Piper Growth Form", ylab = "Chem Diversity (Simpson)")
kruskal.test(simpson ~ groform, data = allsp_simp)
# Kruskal-Wallis chi-squared = 11.431, df = 3, p-value = 0.009611

posthoc.kruskal.dunn.test(simpson ~ groform, data = allsp_simp, p.adjust.method = "bonferroni")

  herb   shrub  tree  
  shrub 0.1126 -      -     
  tree  0.0058 0.6416 -     
  vine  0.0601 1.0000 1.0000

boxplot(simpson ~ habitat, data = allsp_simp, ylim = c(0.9,1.0), xlab = "Piper Successional Habitat", ylab = "Chem Diversity (Simpson)")

kruskal.test(simpson ~ habitat, data = allsp_simp)
# Kruskal-Wallis chi-squared = 2.9487, df = 2, p-value = 0.2289





# 







############################################################### Simpson frt vs. leaves
allsp_simp <- read.csv("allsp_div_simp_100f_fvl.csv", header = T) 

################################## comparing species
bargraph.CI(x.factor = sp, response = simp, data = allsp_simp, err.width = 0.1,las = 1, ylim = c(0,1), xlab = "Piper Species", ylab = "Chem Diversity (Simpson)")
mod.simpsp <- aov(simp~ sp, data = allsp_simp)
summary(mod.simpsp)
#  


################################## comparing tissues
bargraph.CI(x.factor = tissue, response = simp, data = allsp_simp, err.width = 0.1,las = 1, ylim = c(0,1), xlab = "Piper Tissue", ylab = "Chem Diversity (Simpson)")
mod.simpt <- aov(simp ~ tissue, data = allsp_simp)
summary(mod.simpt)
#       




################################# comparing ecological variables
bargraph.CI(x.factor = groform, response = simp, data = allsp_simp, err.width = 0.1,las = 1, ylim = c(0,1), xlab = "Piper Growth Form", ylab = "Chem Diversity (Simpson)")
mod.simpgrof <- aov(simp ~ groform, data = allsp_simp)
summary(mod.simpgrof)
#



bargraph.CI(x.factor = habitat, response = simp, data = allsp_simp, err.width = 0.1,las = 1, ylim = c(0,1), xlab = "Piper Successional Habitat", ylab = "Chem Diversity (Simpson)")
mod.simphab <- aov(simp ~ habitat, data = allsp_simp)
summary(mod.simphab)
#  

mod.simpgrohab <- aov(simp ~ habitat*groform, data = allsp_simp)
summary(mod.simpgrohab)
#  


mod.simpght <- aov(simp ~ habitat*groform*tissue, data = allsp_simp)
summary(mod.simpght)
#
TukeyHSD(mod.simpght)


############################### correlations of organ chemical richness across species
rich_corr <- read.csv("chemrich_corr.csv", header = T)
mod.richcorr_fvl <- lm(chemrich_lf ~ chemrich_frt, data = rich_corr)
summary(mod.richcorr_fvl)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-130.51  -65.04   11.50   34.99  149.79 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  279.4506   111.9928   2.495   0.0317 *
#  chemrich_frt   0.2846     0.3005   0.947   0.3658  

plot(chemrich_lf ~ chemrich_frt, data = rich_corr)
abline(mod.richcorr_fvl, col = "red")

mod.richcorr_rfvl <- lm(chemrich_lf ~ chemrich_rf, data=rich_corr)
summary(mod.richcorr_rfvl)
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-132.250  -68.371    9.573   40.491  149.053 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 264.2227   107.7383   2.452   0.0341 *
#  chemrich_rf   0.3176     0.2808   1.131   0.2843  

plot(chemrich_lf ~ chemrich_rf, data = rich_corr)
abline(mod.richcorr_rfvl, col = "red")

mod.richcorr_rfvuf <- lm(chemrich_rf ~ chemrich_uf, data=rich_corr)
summary(mod.richcorr_rfvuf)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-31.144 -19.607   0.045  12.451  39.300 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   6.6312    29.1716   0.227    0.825    
#chemrich_uf   0.9465     0.0733  12.912 1.46e-07 ***
#Residual standard error: 23.29 on 10 degrees of freedom
#Multiple R-squared:  0.9434,	Adjusted R-squared:  0.9378 
#F-statistic: 166.7 on 1 and 10 DF,  p-value: 1.462e-07

plot(chemrich_rf ~ chemrich_uf, data = rich_corr)
abline(mod.richcorr_rfvuf, col = "red")

mod.richcorr_rfvsd <- lm(chemrich_rf ~ chemrich_sd, data=rich_corr)
summary(mod.richcorr_rfvsd)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-51.769 -38.952   1.503  26.291  58.949 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  30.5304    52.7327   0.579    0.575    
#chemrich_sd   1.0988     0.1646   6.676 5.53e-05 ***
#Residual standard error: 41.91 on 10 degrees of freedom
#Multiple R-squared:  0.8167,	Adjusted R-squared:  0.7984 
#F-statistic: 44.57 on 1 and 10 DF,  p-value: 5.529e-05

plot(chemrich_rf ~ chemrich_sd, data = rich_corr)
abline(mod.richcorr_rfvsd, col = "red")



###############################

class_rich <- read.csv("compound_class_richness.csv", header = T)
bargraph.CI(x.factor = chem_class, response = n_compounds, data = class_rich, ylim = c(0,500), las = 2)
ggplot() +
  geom_bar(aes(y = n_compounds, x = chem_class, fill = chem_class), data = class_rich, stat='identity', las = 1)
           