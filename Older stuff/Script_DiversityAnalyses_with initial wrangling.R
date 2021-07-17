#Note, newest vegan won't work with current adonis code. Updating would be really arduous because
#we are using the pairwise.adonis.2 function and that won't work either, so I used this code to install a previous version
library(devtools)
install_version("vegan", version = "2.5-6", repos = "http://cran.us.r-project.org")



library(dplyr)
library(lme4)
library(ggplot2)
library(gridExtra)
library(multcomp)
library(vegan)
library(lmPerm)
library(viridis)
library(scales)
library(VennDiagram)




#-----------------------------------------------------
#   Data wrangling and functions--run this first and save
#   workspace, then all other sections can be run independently
#-----------------------------------------------------

peaks <- read.csv("peaktable_allsamples_allcompounds.csv")
rownames(peaks) <- peaks[,1]
peaks <- peaks[,-c(1:2)]
peaks <- as.data.frame(t(peaks))
peaks$SampleID <- rownames(peaks)
peaks <- mutate(peaks, richness=rowSums(peaks[1:1311]>100), richness_all=rowSums(peaks[1:1311] != 0))
peaks <- mutate(peaks, extraction=gsub("_.*", "",  gsub("[a-z]*_[a-z]*_", "", SampleID)))
peaks <- peaks[,c(1312:1315, 1:1311)]
peaks[1:5]

IDs <- read.csv("extraction_data_plant_ID.csv")
IDs <- IDs[1:149,2:5]
colnames(IDs)[1] <- "extraction"
IDs$extraction <- as.character(IDs$extraction)
IDs <- mutate(IDs, PlantID=gsub("#", "", gsub(" .*", "", sample_code)))

#re-coding some weird plantIDs for silvivagum and generalense, were just labeled with date or trail marker
IDs$PlantID[which(IDs$PlantID=="19/7/11")] <- "2"
IDs$PlantID[which(IDs$PlantID=="14/5/10")] <- "3"
IDs$PlantID[which(IDs$PlantID=="SSO570")] <- "1"
IDs$PlantID[which(IDs$PlantID=="SSO630")] <- "2"

IDs$PlantID[which(IDs$sample_code=="#7 5/7/12")] <- "7b"



div <- full_join(peaks, IDs, by="extraction")
div <- div[,c(1,4,1316:1319, 2,3,5:1315)]
div[1:9]

#re-coding level for "ripe seed" as just seed for ease of reference
div$tissue <- recode(div$tissue, "ripe seed"="seed")



#correcting spelling
div$sp <- recode(div$sp, multiplinervium="multiplinervum")

#making sp a factor
div$sp <- as.factor(div$sp)

#making sure all plant IDs are unique
div$PlantID <- paste(div$sp, div$PlantID, sep="_")

#making plant ID a factor
div$PlantID <- as.factor(div$PlantID)

#saving current div as the div_all version
div_all <- div

#write.csv(div_all, "Data_Peak_Table.csv")

#now eliminating blank rows and replacing all values < 100 with zero.
#This is per recommendation by Jerry that many of these values are just noise and not
#real peaks in the chromatograms
div <- div[which(!is.na(div$richness)),]
for (i in 9:1319){
  d <- div[i]
  d[which(d<100),] <- 0
  div[i] <- d
}

div <- droplevels(div)

#also droping extra samples for colonense for ripe fruit and peltatum unripe
#where we don't have the other tissues
#could also keep this in for some analyses, but it messes with model fits where we are trying to
#block for individual and also the constrained rarefaction


div <- div[which(div$PlantID !="colonense_12"),]
div <- div[which(div$PlantID !="peltatum_1"),]

#also will create a new variable with all fruit tissues recoded as just fruit,
#so we can examine diversity in fruit as a whole organ in some cases
div$tissue2 <- recode(div$tissue, ripe="frt", unripe="frt", seed="frt")
div <- div[c(1:8, 1320, 9:1319)]

#structural diversity data, intrasample data (all compounds in a sample compared
# to each other)

SD <- read.csv("chem_structural_similarity_intrasample.csv")

#metric shows the "chemical similarity between samples, but I would like to look at 
#diversity, so I will do 1-similarity
SD$SD <- 1-SD$chem_similarity_internal

#adding sample ID info to SD dataset
d.temp <- IDs[,c(1,5)]
colnames(SD)[3] <- "extraction"
SD$extraction <- as.character(SD$extraction)

SD <- full_join(SD, d.temp, by="extraction")

SD <- SD[which(!is.na(SD$run_date)),]

SD$species <- as.factor(SD$species)
SD$species <- recode(SD$species, adu="aduncum", aur="auritum", bio= "biolleyi",
                     col="colonense", gen="generalense", gla="glabrescens",
                     mul="multiplinervum", pel="peltatum", ret="reticulatum",
                     san="sancti-felicis", sil="silvivagum", umb="umbricola")

SD$tissue <- as.factor(SD$tissue)
SD$tissue <- factor(SD$tissue, levels=c("lf","sd", "uf", "rf"))

SD$tissue <- recode(SD$tissue, lf = "leaf", sd = "seed", 
                                      uf="unripe", rf = "ripe")


#making sure all plant IDs are unique
SD$PlantID <- paste(SD$species, SD$PlantID, sep="_")

SID <- div[1:2]

#adding sample ID
SD <- full_join(SD, SID, by="extraction")

SD_write <- SD[c(8,1:2,7,5)]

#write.csv(SD_write, "Data_Intrasample_Structural_Similarity.csv")


###Structural diversity data, intersample--overall structural composition of each sample
##compared to all other samples in a distance matrix

SD_inter <- read.csv("chem_struct_sim_piper12spp_presabs.csv")
SD_inter <- SD_inter[,1:146] #seems like whole matrix was in there twice? Just taking the first one for now--need to check with Jerry

#data are similarity scores, and I want diversity, so I will take 1-sim
SD_inter[2:146] <- 1-SD_inter[2:146]


#Need to extract explanatory variables for this

#species
SD_inter <- mutate(SD_inter, species=gsub("_.*", "", X))# gsub("[a-z]*_", "", X))) 
SD_inter$species <- recode(SD_inter$species, adu="aduncum", aur="auritum", bio= "biolleyi",
                     col="colonense", gen="generalense", gla="glabrescens",
                     mul="multiplinervum", pel="peltatum", ret="reticulatum",
                     san="sancti-felicis", sil="silvivagum", umb="umbricola")

#tissue
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
SD_inter <- mutate(SD_inter, tissue=substrRight(gsub("_[0-9]*", "", X), 2))
SD_inter$tissue <- recode(SD_inter$tissue, lf = "leaf", sd = "seed", 
                    uf="unripe", rf = "ripe")

#plantID

SD_inter <- mutate(SD_inter, extraction=gsub("[a-z]*_", "",substr(X, 1, nchar(X)-7)))

#adding sample ID info to SD_inter dataset
d.temp <- IDs[,c(1,5)]
SD_inter <- full_join(SD_inter, d.temp, by="extraction")
SD_inter <- SD_inter[-c(146:149),]

SD_inter$PlantID <- paste(SD_inter$species, SD_inter$PlantID, sep="_")


write.csv(SD_inter, "Data_Intersample_Structural_Similarity.csv")
#note when I exported this I hadn't run the one line where I did 1-sim to get diversity--I want the uploaded
#final data to match the similarity scores we got as output

#removing col_12 and pel_1 with no other tissues
SD_inter <- SD_inter[which(SD_inter$PlantID !="colonense_12"),]
SD_inter <- SD_inter[which(SD_inter$PlantID !="peltatum_1"),]
SD_inter <- SD_inter[,-which(colnames(SD_inter)=="col_rf_29_180510")]
SD_inter <- SD_inter[,-which(colnames(SD_inter)=="pel_uf_52_180511")]



SD_inter <- SD_inter[, c(145:148, 1:144)]




##setting color palette for all figures
##viridis colors
show_col(viridis_pal()(4))
show_col(viridis_pal()(60))

#would like to have
#leaves=green, seeds=yellow, ripe=purple, unripe=blue

#ordering factor levels so it plots in this order on graphs
div$tissue <- factor(div$tissue, levels=c("leaf","seed", "unripe", "ripe"))

pal <- viridis_pal()(4)[c(3,4,2,1)]
levels(div$tissue)
show_col(pal)

#palette for just leaves and fruit
pal2 <- pal[c(1,4)]


save.image("Workspace_PiperChem")


#subsetting data for Oikos review

d.temp <- div[which(div$sp=="reticulatum"),] #c(1,2,3,5,6, 7)]
d.temp.SD <- SD[which(SD$species=="reticulatum"),c(3,6)]

d.Pret <- full_join(d.temp, d.temp.SD, by="extraction")
plot(d.Pret$tissue, d.Pret$SD)
plot(d.Pret$tissue, d.Pret$richness)

d.Pret <- d.Pret[-c(2,7)]
write.csv(d.Pret, file="Oikos_Fig1_Pret_richness.csv")


#-------------------------------------------------------------
#  Overall PERMANOVAs to assess differences in composition across samples
#--------------------------------------------------

load("Workspace_PiperChem_9dec")

d.temp <- div[-c(1:9)]
d.temp.expl <- div[c(3,5:6)]
m1 <- adonis2(d.temp ~ tissue*sp, strata = 'PlantID', data = d.temp.expl, method = "bray", binary=TRUE, permutations = 999)
m1
#strong effects of tissue, species, and their interaction


#Trying a pairwise adonis to assess pairwise differences among factor levels
#using the function from Pedro Martinez Arbizu on github
#https://github.com/pmartinezarbizu/pairwiseAdonis

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- ad[1]
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 


m1.pw <- pairwise.adonis2(d.temp ~ tissue*sp, strata='PlantID', data = d.temp.expl, 
                          method = "bray", binary=TRUE, nperm=999)
m1.pw
#differences for all pairwise contrasts, but there are sometimes
#strong interactions with species


#checking that the pairwise.adonis function is doing what I think it is doing
d.temp <- div[div$tissue %in% c("seed", "unripe"),-c(1:9)]
d.temp.expl <- div[div$tissue %in% c("seed", "unripe"), c(3,5:6)]
m <- adonis2(d.temp ~ tissue*sp, strata='PlantID', data=d.temp.expl, 
             method="bray", binary=TRUE)
m

#The F-stats are the same as the pairwise.adonis but the p-values are 
#different so I think there is definitely some p-value correction here, 
#but I am not sure what it is in the pairwise.adonis2 function. 
#With pairwise.adonis you can specify
#the p-value correction term but it does not accept interactions or strata



#now trying to split by species and look for tissue level differences
adonis.byspecies <- data.frame(sp=character(), Fstat=numeric(), P=numeric(), 
                        leaf.vs.seed.F=numeric(), leaf.vs.seed.P=numeric(),
                        leaf.vs.unripe.F=numeric(), leaf.vs.unripe.P=numeric(),
                        leaf.vs.ripe.F=numeric(), leaf.vs.ripe.P=numeric(),
                        seed.vs.unripe.F=numeric(), seed.vs.unripe.P=numeric(),
                        seed.vs.ripe.F=numeric(), seed.vs.ripe.P=numeric(),
                        unripe.vs.ripe.F=numeric(), unripe.vs.ripe.P=numeric())


for(i in 1: length(levels(div$sp))){
  sp <- levels(div$sp)[i]
  d <- div[which(div$sp==sp),]
  d.temp <- d[-c(1:9)]
  d.temp.expl <- d[1:9]
  m <- adonis2(d.temp ~ tissue, data=d.temp.expl, 
               method="bray", binary=TRUE)
  Fstat <- m$F[1]
  P <- m$`Pr(>F)`[1]
  m.pw <- pairwise.adonis2(d.temp ~ tissue, data = d.temp.expl, 
                           method = "bray", binary=TRUE)
  F.l.s <- m.pw$leaf_vs_seed$F.Model[1]
  P.l.s <- m.pw$leaf_vs_seed$`Pr(>F)`[1]
  F.l.u <- m.pw$leaf_vs_unripe$F.Model[1]
  P.l.u <- m.pw$leaf_vs_unripe$`Pr(>F)`[1]
  F.l.r <- m.pw$leaf_vs_ripe$F.Model[1]
  P.l.r <- m.pw$leaf_vs_ripe$`Pr(>F)`[1]
  F.s.u <- m.pw$unripe_vs_seed$F.Model[1]
  P.s.u <- m.pw$unripe_vs_seed$`Pr(>F)`[1]
  F.s.r <- m.pw$ripe_vs_seed$F.Model[1]
  P.s.r <- m.pw$ripe_vs_seed$`Pr(>F)`[1]
  F.u.r <- m.pw$ripe_vs_unripe$F.Model[1]
  P.u.r <- m.pw$ripe_vs_unripe$`Pr(>F)`[1]
  newrow <- data.frame(sp=as.character(sp), Fstat=as.numeric(Fstat), P=as.numeric(P), 
                       leaf.vs.seed.F=as.numeric(F.l.s), leaf.vs.seed.P=as.numeric(P.l.s),
                       leaf.vs.unripe.F=as.numeric(F.l.u), leaf.vs.unripe.P=as.numeric(P.l.u),
                       leaf.vs.ripe.F=as.numeric(F.l.r), leaf.vs.ripe.P=as.numeric(P.l.r),
                       seed.vs.unripe.F=as.numeric(F.s.u), seed.vs.unripe.P=as.numeric(P.s.u),
                       seed.vs.ripe.F=as.numeric(F.s.r), seed.vs.ripe.P=as.numeric(P.s.r),
                       unripe.vs.ripe.F=as.numeric(F.u.r), unripe.vs.ripe.P=as.numeric(P.u.r))
  adonis.byspecies <- rbind(adonis.byspecies, newrow)
}



#Basically, we always see an overall effect of tissue, but 
#never any significant pairwise contrasts. Once we get down to that level
#we just have N=3, so I just think that we do not have the sample size to do 
#PERMANOVAs on these very small subsets...we are getting warnings that 
#Set of permutations < 'minperm'. Generating entire set.
#'nperm' >= set of all permutations: complete enumeration.

#The minimum p-value for a PERMANOVA is determined by the number
#of permutations, e.g. for nperm=999 the minimum is p=0.001, which
#means that the differences between groups could not be replicated
#in any of the 999 permutations

#So with N=3 and strata = plantID, I think we only have 9 permutations,
#so the minimum p is 0.1

#I think we can just report the overall effects of tissue for the individual
#species



#can also try this with the lmPerm package

#calculate the distance matrix
w <- vegdist(d.temp, binary=TRUE, method="bray")
plot(w)
w <- as.matrix(w)


d.temp.dist <- data.frame(w=w[lower.tri(w)], tissue=NA, sp=NA )
t.vs.t <- expand.grid(d.temp.expl$tissue, d.temp.expl$tissue)
t <- t.vs.t[lower.tri(w),]
d.temp.dist$tissue <- paste(t[,1], t[,2], sep=".vs.")
sp.vs.sp <- expand.grid(d.temp.expl$sp, d.temp.expl$sp)
sp <- sp.vs.sp[lower.tri(w),]
d.temp.dist$sp <- paste(sp[,1], sp[,2], sep=".vs.")
d.temp.dist$tissue <- as.factor(d.temp.dist$tissue)
d.temp.dist$sp <- as.factor(d.temp.dist$sp)


#PERMANOVA model
m1 <- aovp(w ~ tissue*sp, data = d.temp.dist, maxIter = 10000, perm = "Prob")
summary(m1)
TukeyHSD(m1)





##Now trying a PERMANOVA on the structural data

#in this case, we feed adonis2 our own distance matrix rather than the dataframe
d.temp <- as.dist(SD_inter[,6:148])
d.temp.expl <- SD_inter[,c(1:5)]
m2 <- adonis2(d.temp ~ tissue*species, strata='PlantID', data = d.temp.expl, permutations = 999)
m2
plot(d.temp)

m2.pw <- pairwise.adonis2(d.temp ~ tissue*species, strata='PlantID', data = d.temp.expl, 
                          nperm=999)
m2.pw


#now trying to split by species and look for tissue level differences
adonis.byspecies.SD <- data.frame(sp=character(), Fstat=numeric(), P=numeric(), 
                               leaf.vs.seed.F=numeric(), leaf.vs.seed.P=numeric(),
                               leaf.vs.unripe.F=numeric(), leaf.vs.unripe.P=numeric(),
                               leaf.vs.ripe.F=numeric(), leaf.vs.ripe.P=numeric(),
                               seed.vs.unripe.F=numeric(), seed.vs.unripe.P=numeric(),
                               seed.vs.ripe.F=numeric(), seed.vs.ripe.P=numeric(),
                               unripe.vs.ripe.F=numeric(), unripe.vs.ripe.P=numeric())


for(i in 1: length(levels(div$sp))){
  sp <- levels(div$sp)[i]
  d <- div[which(div$sp==sp),]
  d.temp <- d[-c(1:9)]
  d.temp.expl <- d[1:9]
  m <- adonis2(d.temp ~ tissue, data=d.temp.expl, 
               method="bray", binary=TRUE)
  Fstat <- m$F[1]
  P <- m$`Pr(>F)`[1]
  m.pw <- pairwise.adonis2(d.temp ~ tissue, data = d.temp.expl, 
                           method = "bray", binary=TRUE)
  F.l.s <- m.pw$leaf_vs_seed$F.Model[1]
  P.l.s <- m.pw$leaf_vs_seed$`Pr(>F)`[1]
  F.l.u <- m.pw$leaf_vs_unripe$F.Model[1]
  P.l.u <- m.pw$leaf_vs_unripe$`Pr(>F)`[1]
  F.l.r <- m.pw$leaf_vs_ripe$F.Model[1]
  P.l.r <- m.pw$leaf_vs_ripe$`Pr(>F)`[1]
  F.s.u <- m.pw$unripe_vs_seed$F.Model[1]
  P.s.u <- m.pw$unripe_vs_seed$`Pr(>F)`[1]
  F.s.r <- m.pw$ripe_vs_seed$F.Model[1]
  P.s.r <- m.pw$ripe_vs_seed$`Pr(>F)`[1]
  F.u.r <- m.pw$ripe_vs_unripe$F.Model[1]
  P.u.r <- m.pw$ripe_vs_unripe$`Pr(>F)`[1]
  newrow <- data.frame(sp=as.character(sp), Fstat=as.numeric(Fstat), P=as.numeric(P), 
                       leaf.vs.seed.F=as.numeric(F.l.s), leaf.vs.seed.P=as.numeric(P.l.s),
                       leaf.vs.unripe.F=as.numeric(F.l.u), leaf.vs.unripe.P=as.numeric(P.l.u),
                       leaf.vs.ripe.F=as.numeric(F.l.r), leaf.vs.ripe.P=as.numeric(P.l.r),
                       seed.vs.unripe.F=as.numeric(F.s.u), seed.vs.unripe.P=as.numeric(P.s.u),
                       seed.vs.ripe.F=as.numeric(F.s.r), seed.vs.ripe.P=as.numeric(P.s.r),
                       unripe.vs.ripe.F=as.numeric(F.u.r), unripe.vs.ripe.P=as.numeric(P.u.r))
  adonis.byspecies.SD <- rbind(adonis.byspecies.SD, newrow)
}


#combined by species data for table

tbl <- cbind(adonis.byspecies[1:3], adonis.byspecies.SD[2:3])

write.csv(tbl, file="adonis.byspecies.csv")







#-------------------------------------------------------------
#  Gamma diversity
#-----------------------------------------------------------

load("Workspace_PiperChem")
#load("Workspace_PiperChem_Gamma")  #last rarefaction run with 5000 reps

#-----------Venn------------------
#First, will show a venn diagram with the total count of compounds detected in each
#tissue type and shared among tissues


levels(div$tissue)
d <- div[which(div$tissue=="leaf"),10:1320]
lf <- colnames(d)[which(colSums(d) !=0)]

d <- div[which(div$tissue=="ripe"),10:1320]
r <- colnames(d)[which(colSums(d) !=0)]

d <- div[which(div$tissue=="seed"),10:1320]
s <- colnames(d)[which(colSums(d) !=0)]

d <- div[which(div$tissue=="unripe"),10:1320]
u <- colnames(d)[which(colSums(d) !=0)]

x <- list(lf, s, u, r)


#make a list of compounds shared in each category of overlap, e.g. in everything, in 
#seeds and leaves, etc.
overlap <- calculate.overlap(x)
 

#use this to double check groupings...can see data for a particular compund and which
#tissues it was found in
div[,c("tissue","577.1393_3.1992")]

all <- length(overlap$a6)
l.u.s <- length(overlap$a12)
l.r.s <- length(overlap$a11)
l.u.r <- length(overlap$a5)
r.u.s <- length(overlap$a7)
l.s <- length(overlap$a15)
l.u <- length(overlap$a4)
l.r <- length(overlap$a10)
u.s <- length(overlap$a13)
r.s <- length(overlap$a8)
u.r <- length(overlap$a2)
l <- length(overlap$a9)



#set1=lf
#set2=r
#set3=s
#set4=u

palVenn <- pal[c(1,4,2,3)]
png(height=500, width=500, filename="Venn.tiff", type="cairo")
#tiff(height=500, width=500, filename="Venn.tiff", type="cairo")
#cairo_pdf(height=7, width=7, filename="Venn.pdf")
#cairo_ps(height=7, width=7, filename="Venn.eps")
venn <- draw.quad.venn(length(lf), length(r), length(s), length(u), length(intersect(lf, r)),
                       length(intersect(lf, s)),length(intersect(lf, u)),length(intersect(r,s)),
                       length(intersect(r,u)),length(intersect(s,u)), 
                       length(Reduce(intersect, list(lf,r,s))), 
                       length(Reduce(intersect, list(lf,r,u))),
                       length(Reduce(intersect, list(lf,s,u))),
                       length(Reduce(intersect, list(r,s,u))),
                       length(Reduce(intersect, list(lf,r,u,s))),
                       category=c("leaves", "ripe\npulp", "seeds","unripe\npulp"), 
                       fill=palVenn, alpha=c(0.6,0.6,0.6,0.6), cex=2, cat.cex=2,
                       cat.just=list(c(0.5,0),c(0.5,0),c(0.5,0),c(0.5,0)))
dev.off()

#-----------Rarefaction------------

#Our data contains replicated samples (n=3) for each plant species (n=12)
#So we want to know about total gamma diversity across all species
#but we want to keep the information about intraspecific variation 
#as well. So, I am trying a constrained version of rarefaction, basically the same principle as
#spatial rarefaction, where the species accumulation curve always moves to the 
#nearest neighbor. Here instead we will force it to go to the replicate
#within a species, rather than just choosing the next sample at random
#in the species accumulation curve


#Building a function 'CPR' for the actual rarefaction part
#CPR for constrained phytochemical rarefaction

#Inputs for function are 'compounds', 'meta', 'reps', and 'mod' 

#compounds should be a dataframe with samples as rows and compounds as columns
#observations can be presence/absence or peak areas/quantities, but this function uses
#presence/absence only as in sample-based rarefaction

#meta should be a dataframe with metadata. The first column should be a unique sample ID code
#given in the same order as the samples in 'compounds', and the second column
#should be the species identity or other blocking factor in a nested sampling
#design (i.e. any clear grouping factor)

#reps is the number of random starting points for the rarefaction, default is 100

#mod is model or list of models that are accepted by fitspecaccum in vegan; default is to
#try a list c("arrhenius","gitay","lomolino","asymp", "gompertz", "logis")

#CPR returns a list of dataframes: 

#mean_coef contains the mean, SE, and CIs for model
#parameters specified in mod

#mean_AIC contains the mean, SE, and CIs for AIC values for model fit (can be used
#to compare fit among these models)

#spAcumCurves contains the mean and CI for the species accumulation curves (can be
#used for plotting)

#Other notes for myself:
#note that specaccum can also take weights based on the sampling effort
#at each site. What if we included SM total abundance as the weight??, basically controlling
#for "sample effort"; currently this is not implemented in the function

#Use for debugging function:
compounds <- d.temp
meta <- m.temp
reps <- 100
mod=c("arrhenius","gitay","lomolino","asymp", 
           "gompertz", "logis")
      

CPR <- function(compounds, meta, reps=100, mod=c("arrhenius","gitay","lomolino","asymp", 
                                                       "gompertz", "logis")) {
  spAcumCurves <- data.frame(samples=1:nrow(compounds))
  library(vegan)
  library(dplyr)
  library(tidyr)
  N <- reps
  #making a list of dataframes that will have the mod coefficients and AIC values
  data(BCI)
  all_coef <- list()
  for(i in 1:length(mod)){
    m.temp <- mod[i]
    coef_names <- names(coef(fitspecaccum(specaccum(BCI), m.temp)))
    coef <- matrix(nrow=N, ncol=length(coef_names))
    coef <- as.data.frame(coef)
    names(coef) <- c(coef_names)
    coef$AIC <- NA
    all_coef[[i]] <- coef
    names(all_coef)[[i]] <- m.temp
  }
  pb <- txtProgressBar(min = 0, max = N, initial = 0, style=3)
  for(i in 1:N){
    r <- sample(1:nrow(compounds))[1]
    f <- meta[r,1] #define the focal plant
    sp <- meta[r,2] #define the species
    names(meta) <- c("SampleID", "sp")
     #define an order for the species, starting with the focal species and random after that
    spOrd <- data.frame(sp=c(as.character(sp),sample(levels(meta$sp)[which(levels(meta$sp) !=sp)])), 
                        ord=1:length(levels(meta$sp)))
    d.temp <- cbind(meta, compounds)
    d.temp2 <- inner_join(d.temp, spOrd, by="sp") 
    #randomly order samples
    rows <- sample(nrow(d.temp2))
    d.temp2 <- d.temp2[rows,]
    #now sort according to species order, put focal plant first
    d.temp2$ord[which(d.temp2$SampleID==f)] <- 0
    d.temp2 <- d.temp2[order(d.temp2$ord),]
    c <- specaccum(d.temp2[,-c(1,2,ncol(d.temp2))], method="collector")
    spAcumCurves <- cbind(spAcumCurves, c$richness)
    names(spAcumCurves)[i+1] <- f
    for(j in 1:length(mod)){
      fit <- try(fitspecaccum(c, mod[j], control = list(maxiter = 500)))
      if(class(fit)[1]=="try-error"){
        c.temp=rep(NA, ncol(coef))}else{
          c.temp <- c(coef(fit), AIC(fit))}
      all_coef[[j]][i,] <- c.temp
    }
    setTxtProgressBar(pb, i)
  }
  result <- spAcumCurves[-1]
  r.mean <- rowMeans(result)
  r.sd <- apply(result, 1, sd)
  r.n <- ncol(result)
  CI.high <-  r.mean + (1.96*(r.sd/sqrt(r.n)))
  CI.low <- r.mean - (1.96*(r.sd/sqrt(r.n)))
  curves <- data.frame(samples=spAcumCurves$samples, mean=r.mean, 
                       CI_high=CI.high, CI_low=CI.low)
  se <- function(x, na.rm) sd(x, na.rm=na.rm)/sqrt(length(x))
  ci.high <- function(x, na.rm) mean(x, na.rm=na.rm) + (1.96*sd(x, na.rm=na.rm)/sqrt(length(x)))
  ci.low <- function(x, na.rm) mean(x, na.rm=na.rm) - (1.96*sd(x, na.rm=na.rm)/sqrt(length(x)))
  mean_coef <- data.frame(model_type=character(),parameter=character(), mean=numeric(),
                                       SE=numeric(), CI.high=numeric(), CI.low=numeric())
  mean_AIC <- data.frame(model_type=mod, meanAIC=NA, SE=NA, CI.high=NA, CI.low=NA)
  for(i in 1:length(mod)){
    d.temp <- all_coef[[i]]
    d.temp2 <- d.temp[1:ncol(d.temp)-1]
    m.temp <- d.temp2 %>%
      summarise_all(list(mean=mean, SE=se, CI.high=ci.high, CI.low=ci.low), na.rm=TRUE) %>%
      pivot_longer(cols=ends_with(c("mean", "SE", "high", "low"))) %>%
      separate (name, into = c("parameter", "sum_stat"), sep = "_" ) %>%
      pivot_wider(names_from=sum_stat, values_from=value)
    m.temp <- as.data.frame(m.temp)
    m.temp$model_type <- mod[i]
    m.temp <- m.temp[c(ncol(m.temp), 1:ncol(m.temp)-1)] 
    mean_coef <- rbind(mean_coef, m.temp)
    #now summarise AIC values
    mean_AIC[i,2] <- mean(d.temp$AIC, na.rm=TRUE)
    mean_AIC[i,3] <- se(d.temp$AIC, na.rm=TRUE)
    mean_AIC[i,4] <- ci.high(d.temp$AIC, na.rm=TRUE)
    mean_AIC[i,5] <- ci.low(d.temp$AIC, na.rm=TRUE)
  }
  CPR <- list(mean_coef=mean_coef, mean_AIC=mean_AIC, spAcumCurves=curves)
  return(CPR)
}
  



#now applying this across all tissues  
  
allCurves <- data.frame(tiss=factor(levels=levels(div$tissue)), samples=numeric(), mean=numeric(), 
                        CI_high=numeric(), CI_low=numeric())
allCoefEst <- data.frame(tiss=factor(levels=levels(div$tissue)),
                         model_type=character(), parameter=character(), mean=numeric(),
                         SE=numeric(), CI.high=numeric(), CI.low=numeric())
allAIC <- data.frame(tiss=factor(levels=levels(div$tissue)),
                         meanAIC=numeric(), SE=numeric(), CI.high=numeric(), CI.low=numeric())
for(i in 1:length(levels(div$tissue))){
  tiss <- levels(div$tissue)[i]
  d.temp <- div[which(div$tissue==tiss), -c(1:9)]
  m.temp <- div[which(div$tissue==tiss), c(1,3)]
  CPR_out <- CPR(d.temp, m.temp, reps=500) 
  curves <- CPR_out[[3]]
  curves$tiss <- tiss
  curves <- curves[c(5,1:4)]
  allCurves <- rbind(allCurves, curves)
  coefs <- CPR_out[[1]]
  coefs$tiss <- tiss
  coefs <- coefs[c(7,1:6)]
  allCoefEst <- rbind(allCoefEst, coefs)
  AICs <- CPR_out[[2]]
  AICs$tiss <- tiss
  AICs <- AICs[c(6,1:5)]
  allAIC <- rbind(allAIC, AICs) 
}

allCurves
allCoefEst
allAIC

#varies a bit between runs due to random sampling, but asymptote is overall the best fit
#lomolino is close and occasionally better for leaf or seed, but lomolino estimates
#for asymtote are really variable and sometimes WAY higher that observed (e.g. 35,000 compounds) 

asymp.finalest <- allCoefEst[which(allCoefEst$model_type=="asymp" & allCoefEst$parameter=="Asym"),]
asymp.finalest
#save these for table
write.csv(asymp.finalest, file="gammaestimates.csv")

#leaf is lower than all other fruit tissues, 95% CIs do not cross. 

# Plot rarefaction

#ordering factor levels so it plots in this order on graphs
allCurves$tiss <- factor(allCurves$tiss, levels=c("leaf","seed", "unripe", "ripe"))
p <- ggplot(data=allCurves) +
  geom_line(aes(x=samples, y=mean, group=tiss, color=tiss), size=1) +
  #geom_ribbon(aes(x=samples, ymin=CI_low, ymax=CI_high, group=tiss, fill=tiss), alpha=0.2) +
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  xlab('No. of samples') + ylab('Cumulative richness') +
  labs(subtitle="(B)")+
  theme_classic() +
  theme(legend.position=c(0.8,0.4), legend.title=element_blank(),
        legend.key = element_rect(colour = "transparent", fill="transparent")) + 
  #plot.margin = margin(0.25,0.25,0.25,0.75, unit='lines')) +
  xlim(0,37) +
  ylim(0, 1400)
p

#ggsave("Rarefaction.png", width = 8, height = 6.5, units = "cm")
#ggsave("Rarefaction.eps", width = 8, height = 6.5, units = "cm", device=cairo_ps)
ggsave("Rarefaction.pdf", width = 8, height = 6.5, units = "cm", device = cairo_pdf)





#-------------------------------------------------------------------
#  Alpha diversity--looking at average richness in each sample type
#------------------------------------------------------------------------


load("Workspace_PiperChem")


hist(div$richness)
hist(div$richness_all)
div <- div[which(!is.na(div$richness)),]
div <- droplevels(div)

plot(richness ~ tissue, data=div)
plot(richness ~ sp, data=div)

m1 <- lmer(richness ~ tissue * sp + (1|PlantID), data=div) 
summary(m1)
anova(m1, test="F")
drop1(m1, test="Chisq")

interaction.plot(div$tissue, div$sp, div$richness)

p <- ggplot(div, aes(tissue, richness)) + 
  geom_boxplot() +
  facet_wrap(vars(sp))
p

#clear significant interaction between species and tissue, splitting by species

lmm.byspecies <- data.frame(sp=character(), X=numeric(), P=numeric(), leaf=character(), 
                            seed=character(), unripe=character(), ripe=character())

for(i in 1: length(levels(div$sp))){
  sp <- levels(div$sp)[i]
  d <- div[which(div$sp==sp),]
  m <- lmer(richness~tissue + (1|PlantID), data=d)
  #m <- lm(richness~tissue, data=d)
  X <- drop1(m, test="Chisq")[2,3]
  P <- drop1(m, test="Chisq")[2,4]
  #X <- anova(m)[1,4]
  #P <- anova(m)[1,5]
  t <- cld(glht(m, mcp(tissue="Tukey")))
  l <- t$mcletters$Letters
  newrow <- data.frame(sp=as.character(sp), X=as.numeric(X), P=as.numeric(P), 
                          leaf=as.character(l[1]), seed=as.character(l[2]), 
                          unripe=as.character(l[3]), ripe=as.character(l[4]))
  lmm.byspecies <- rbind(lmm.byspecies, newrow)
}

#Some boundary, singular fit errors. These all go away when you run this as a simple lm
#without the plantID random effect. However, plantID explains a large portion of the 
#variance for many of the models

#in most cases, singular fit errors occur for random effects only models used for LRTs
#There are only a few cases where there are singular fit errors for main model 
#(with tissue as a factor) and this is usually when tissue is not significant, as in 
#biolleyi, peltatum, but also occurs for silvivagum

#in all cases, results are qualitatively similar with or without random effect
#so I will keep plantID in models (for consistency and because it explains a large part of
#the variance for many of the species and is an important part of the design)


#pretty plot

TukeyLab <- data.frame(sp=character(), tissue=character(),x=character(), max=numeric())
for (i in 1:length(levels(div$sp))) {
  sp <- levels(div$sp)[i]
  for(j in 1:length(levels(div$tissue))){
    t <- levels(div$tissue)[j]
    newrow <- data.frame(sp=sp, tissue=t, lab=lmm.byspecies[which(lmm.byspecies$sp==sp),
                                                            which(colnames(lmm.byspecies)==t)],
                         max=max(div$richness[which(div$sp==sp & div$tissue==t)]))
    TukeyLab <- rbind(TukeyLab, newrow)
  }
}

TukeyLab2 <- TukeyLab[-c(9:12,29:32),]

ypos_lab <- TukeyLab2$max + 60 
xpos_lab <- c(0.7, 0.9, 1.1, 1.3, 
              1.7, 1.9, 2.1, 2.3,
              #2.7, 2.9, 3.1, 3.3,
              3.7, 3.9, 4.1, 4.3,
              4.7, 4.9, 5.1, 5.3,
              5.7, 5.9, 6.1, 6.3,
              6.7, 6.9, 7.1, 7.3,
              #7.7, 7.9, 8.1, 8.3,
              8.7, 8.9, 9.1, 9.3,
              9.7, 9.9, 10.1, 10.3,
              10.7, 10.9, 11.1, 11.3,
              11.7, 11.9, 12.1, 12.3)


TukeyLab2$xpos <- rep(c(1,2,3,4), 10)

Plabs <- character()
for (i in 1:length(lmm.byspecies$sp)){
  sp <- lmm.byspecies$sp[i]
  if (lmm.byspecies$P[which(lmm.byspecies$sp==sp)] < 0.001){
    lab <- "P < 0.001"
  }else{
    lab <- paste("P ==", round(lmm.byspecies$P[i], 3)) 
  }
  Plabs <- c(Plabs, lab)
}

LRT_Lab <- data.frame(sp=lmm.byspecies$sp, 
                      tissue=factor("leaf", levels=c("leaf", "seed", "unripe", "ripe")),
                      labs=sprintf("italic(X^2) == %.2f~%s", lmm.byspecies$X, Plabs),
                      xpos=rep(2.5, 12), ypos=c(rep(175,4), 700, 175, 175, 700, rep(175, 4)))


p <- ggplot(div, aes(tissue, richness, fill=tissue)) + 
  geom_boxplot() +
  scale_fill_manual(values=pal, aesthetics = "fill") +
  labs(x=element_blank(), y="Compound Richness")+
  facet_wrap(vars(sp), ncol=4)+
  geom_text(data=TukeyLab2, mapping=aes(x=xpos, y=max+50, label=lab)) +
  geom_text(data=LRT_Lab, mapping=aes(x=xpos, y=ypos, label=labs), 
            parse=TRUE, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background =element_rect(fill=pal[4]))+
  theme(strip.text = element_text(colour = 'white')) +
  theme(legend.title=element_blank(), legend.position = "none") +
  theme (axis.text.x = element_text(angle = 45, vjust=0.65))
p


ggsave("Alphadiv.png", width = 15, height = 15, units = "cm")
ggsave("Alphadiv.eps", width = 15, height = 15, units = "cm", device=cairo_ps)
ggsave("Alphadiv.pdf", width = 15, height = 15, units = "cm", device = cairo_pdf)




#alternative plot

lines <- data.frame(lines=c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5))

xlab <- substitute(paste(italic('Piper'), " species"))

TukeyLab <- data.frame(sp=character(), tissue=character(),x=character(), max=numeric())
for (i in 1:length(levels(div$sp))) {
  sp <- levels(div$sp)[i]
  for(j in 1:length(levels(div$tissue))){
    t <- levels(div$tissue)[j]
    newrow <- data.frame(sp=sp, tissue=t, lab=lmm.byspecies[which(lmm.byspecies$sp==sp),
                                                            which(colnames(lmm.byspecies)==t)],
                         max=max(div$richness[which(div$sp==sp & div$tissue==t)]))
    TukeyLab <- rbind(TukeyLab, newrow)
  }
}

TukeyLab2 <- TukeyLab[-c(9:12,29:32),]

ypos_lab <- TukeyLab2$max + 60 
xpos_lab <- c(0.7, 0.9, 1.1, 1.3, 
              1.7, 1.9, 2.1, 2.3,
              #2.7, 2.9, 3.1, 3.3,
              3.7, 3.9, 4.1, 4.3,
              4.7, 4.9, 5.1, 5.3,
              5.7, 5.9, 6.1, 6.3,
              6.7, 6.9, 7.1, 7.3,
              #7.7, 7.9, 8.1, 8.3,
              8.7, 8.9, 9.1, 9.3,
              9.7, 9.9, 10.1, 10.3,
              10.7, 10.9, 11.1, 11.3,
              11.7, 11.9, 12.1, 12.3)

p <- ggplot(div, aes(sp, richness, fill=tissue)) + 
  geom_boxplot() +
  labs(x=xlab, y="Compound Richness") +
  scale_fill_manual(values=pal, aesthetics = "fill") +
  geom_vline(data = lines, aes(xintercept = as.numeric(lines)), lty="dotted") +
  annotate("text", x = xpos_lab, y = ypos_lab, label = TukeyLab2$lab, size=3) +
  theme_classic() +
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, vjust=0.65))
p
ggsave("Alphadiv2.png", width = 25, height = 10, units = "cm")
#ggsave("Alphadiv2.eps", width = 30, height = 7, units = "cm", device=cairo_ps)
#ggsave("Alphadiv2.pdf", width = 30, height = 7, units = "cm", device = cairo_pdf)





#-----------------------------------------------------------------
#  Beta-diversity
#------------------------------------------------------------------

load("Workspace_PiperChem")

d.temp <- div[-c(1:9)]

#first calculating the pairwise distances for all samples
#this would be the same as betadiver(d.temp, "w")
#it is equivalent to the Sorensen dissimilarity index
#β_w = (b+c)/(2 a + b + c).  where a are shared species, b and c are unique species to each site
w <- vegdist(d.temp, binary=TRUE, method="bray")
plot(w)
hist(w)

#then we can look at the distances between each sample and the group centroid
#but we need to define a specific group (e.g. tissue)
wb <- betadisper(w, div$tissue)
wb
anova(wb)
permutest(wb) 
TukeyHSD(wb)
boxplot(wb)



#note that the results are very different if we use the quantitative data
w2 <- vegdist(d.temp, method="bray")
plot(w2)
hist(w2)
wb2 <- betadisper(w2, div$tissue)
wb2
permutest(wb2)
TukeyHSD(wb2)
boxplot(wb2)




#Data for plot plot
beta.div <- data.frame(tissue=wb$group, beta.div=wb$distances)
beta.div$tissue <- factor(beta.div$tissue, levels=c("leaf","seed", "unripe", "ripe"))



#can also look at differences for species
wb_sp <- betadisper(w, div$sp)
anova(wb_sp)
TukeyHSD(wb_sp)
boxplot(wb_sp)
#main difference is that generalense seems to have higher beta-diversity
#than other species. Note all the tissues are in there, so this probably is driven
#by major differences among tissues


##in the dataset above, all samples are together, so the beta-div
#in fruits might be due to within species variance in addition
#to across species variance. Two things to do: 1) look at beta
#diversity on species averages (for across species variance); 2) look at beta div for each
#species separately (to see if there is within species variance)


#with the species averaged first
div_SpAvg <- div %>%
  group_by(sp, tissue, tissue2) %>%
  summarize_if(is.numeric, max)

d.temp <- div_SpAvg[-c(1:5)]
w <- vegdist(d.temp, binary=TRUE, method="bray")
plot(w)

#then we can look at the distances between each sample and the group centroid
#but we need to define a specific group (e.g. tissue)
wb <- betadisper(w, div_SpAvg$tissue)
anova(wb)
TukeyHSD(wb)
boxplot(wb)
#same trend but not significant, likely because our sample size is only 12 now


#now trying to split by species
beta.byspecies <- data.frame(sp=character(), Fstat=numeric(), P=numeric(), leaf=character(), 
                             seed=character(), unripe=character(), ripe=character())

for(i in 1: length(levels(div$sp))){
  sp <- levels(div$sp)[i]
  d <- div[which(div$sp==sp),]
  d2 <- d[-c(1:9)]
  w <- vegdist(d2, binary=TRUE, method="bray")
  wb <- betadisper(w, d$tissue)
  d3 <- data.frame(tissue=wb$group, beta.div=wb$distances)
  m <- aov(beta.div ~ tissue, data=d3)
  Fstat <- summary(m)[[1]][1,4]
  P <- summary(m)[[1]][1,5]
  t <- cld(glht(m, mcp(tissue="Tukey")))
  l <- t$mcletters$Letters
  newrow <- data.frame(sp=as.character(sp), Fstat=as.numeric(Fstat), P=as.numeric(P), 
                       leaf=as.character(l[1]), seed=as.character(l[2]), 
                       unripe=as.character(l[3]), ripe=as.character(l[4]))
  beta.byspecies <- rbind(beta.byspecies, newrow)
}

#there is only one case where there are significant differences in 
#beta diversity across tissues within a species, that is for colonense
#Is this due to very low sample size (n=3 samples per tissue per species)? Or is it that 
#the main tissue beta-diversity effects are driven by species
#level differences (i.e. fruits of species A are very different from species B, etc)

#I am still struggling with the fact that just doing tissue ignores
#this important factor in the dataset (species) and maybe would
#be like psuedoreplication??  and you can't do two factors
#simultaneously in betadisper

#solution suggested here is to do the two-factor adonis with interaction, and if 
#there is no interaction, to just interpret each factor individually in
#betadisper
#https://stat.ethz.ch/pipermail/r-sig-ecology/2010-September/001524.html

d.temp <- div[-c(1:9, 50:1311)]
d.temp.expl <- div[c(3,5)]
w <- vegdist(d.temp, binary=TRUE, method="bray")
wb <- adonis2(d.temp ~ tissue*sp, data=d.temp.expl)
wb

#There is a strong interaction here between tissue and species, so that doesn't really work. 



##Also just found this paper:
#Anderson, Marti J. (2014) "Permutational multivariate analysis of variance (PERMANOVA)."
#See Table 2 and Fig. 4 in that paper...they show an example where PERMANOVA is used to
#partition the variance in multivariate composition among different spatial scales
#something like this would be perfect, to show both the within species and across species
#components of variation

## Using adonis to quantify the contribution of species to tissue-level variance
# subset data by tissue and run adonis with species as the independent variable
# Use Sum of Squares of residuals = SSE to calculate sigma^2 = SSE/(n-v)

d.temp <- div[which(div$tissue=="leaf"), -c(1:9)]
d.temp.expl <- div[which(div$tissue=="leaf"), c(3,5)]
adonis_leaf <- adonis2(d.temp ~ sp, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_leaf


d.temp <- div[which(div$tissue=="ripe"), -c(1:9)]
d.temp.expl <- div[which(div$tissue=="ripe"), c(3,5)]
adonis_ripe <- adonis2(d.temp ~ sp, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_ripe

d.temp <- div[which(div$tissue=="unripe"), -c(1:9)]
d.temp.expl <- div[which(div$tissue=="unripe"), c(3,5)]
adonis_unripe <- adonis2(d.temp ~ sp, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_unripe

d.temp <- div[which(div$tissue=="seed"), -c(1:9)]
d.temp.expl <- div[which(div$tissue=="seed"), c(3,5)]
adonis_seed <- adonis2(d.temp ~ sp, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_seed


#For each tissue, these values can tell us the proportion of variance in composition that is due to 
#interspecific vs intraspecific variation


prop.variance <- data.frame(
  tissue="leaf", Fstat=adonis_leaf$F[1], P=adonis_leaf$`Pr(>F)`[1],
  prop.var.sp=adonis_leaf$SumOfSqs[1]/adonis_leaf$SumOfSqs[3],
  prop.var.resid=adonis_leaf$SumOfSqs[2]/adonis_leaf$SumOfSqs[3]
)

prop.variance <- rbind(prop.variance, data.frame(
  tissue="seed", Fstat=adonis_seed$F[1], P=adonis_seed$`Pr(>F)`[1],
  prop.var.sp=adonis_seed$SumOfSqs[1]/adonis_seed$SumOfSqs[3],
  prop.var.resid=adonis_seed$SumOfSqs[2]/adonis_seed$SumOfSqs[3]
))

prop.variance <- rbind(prop.variance, data.frame(
  tissue="unripe", Fstat=adonis_unripe$F[1], P=adonis_unripe$`Pr(>F)`[1],
  prop.var.sp=adonis_unripe$SumOfSqs[1]/adonis_unripe$SumOfSqs[3],
  prop.var.resid=adonis_unripe$SumOfSqs[2]/adonis_unripe$SumOfSqs[3]
))

prop.variance <- rbind(prop.variance, data.frame(
  tissue="ripe", Fstat=adonis_ripe$F[1], P=adonis_ripe$`Pr(>F)`[1],
  prop.var.sp=adonis_ripe$SumOfSqs[1]/adonis_ripe$SumOfSqs[3],
  prop.var.resid=adonis_ripe$SumOfSqs[2]/adonis_ripe$SumOfSqs[3]
))





##Now trying a beta-diversity analysis on the structural data--so this is something 
#like "structural beta diversity"

d.temp <- as.dist(SD_inter[,6:148])
plot(d.temp)

#then we can look at the distances between each sample and the group centroid
#but we need to define a specific group (e.g. tissue)
wb <- betadisper(d.temp, SD_inter$tissue)
permutest(wb)
TukeyHSD(wb)
boxplot(wb)

##structural beta-diversity higher for fruits (at least seeds and ripe) than in leaves

#Data for plot
s.beta.div <- data.frame(tissue=wb$group, beta.div=wb$distances)
s.beta.div$tissue <- factor(s.beta.div$tissue, levels=c("leaf","seed", "unripe", "ripe"))


## Using adonis to quantify the contribution of species to tissue-level variance
# subset data by tissue and run adonis with species as the independent variable
# Use Sum of Squares of residuals = SSE to calculate sigma^2 = SSE/(n-v)

d.temp <- as.dist(SD_inter[which(SD_inter$tissue=="leaf"), 
                           which(grepl("lf" , colnames(SD_inter)))  ])
d.temp.expl <- SD_inter[which(SD_inter$tissue=="leaf"), c(1:5)]
adonis_leaf <- adonis2(d.temp ~ species, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_leaf


d.temp <- as.dist(SD_inter[which(SD_inter$tissue=="seed"), 
                           which(grepl("sd" , colnames(SD_inter)))  ])
d.temp.expl <- SD_inter[which(SD_inter$tissue=="seed"), c(1:5)]
adonis_seed <- adonis2(d.temp ~ species, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_seed

d.temp <- as.dist(SD_inter[which(SD_inter$tissue=="unripe"), 
                           which(grepl("uf" , colnames(SD_inter)))  ])
d.temp.expl <- SD_inter[which(SD_inter$tissue=="unripe"), c(1:5)]
adonis_unripe <- adonis2(d.temp ~ species, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_unripe

d.temp <- as.dist(SD_inter[which(SD_inter$tissue=="ripe"), 
                           which(grepl("rf" , colnames(SD_inter)))  ])
d.temp.expl <- SD_inter[which(SD_inter$tissue=="ripe"), c(1:5)]
adonis_ripe <- adonis2(d.temp ~ species, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_ripe

#For each tissue, these values can tell us the proportion of variance in composition that is due to 
#interspecific vs intraspecific variation


prop.variance.SD <- data.frame(
  tissue="leaf", Fstat=adonis_leaf$F[1], P=adonis_leaf$`Pr(>F)`[1],
  prop.var.sp=adonis_leaf$SumOfSqs[1]/adonis_leaf$SumOfSqs[3],
  prop.var.resid=adonis_leaf$SumOfSqs[2]/adonis_leaf$SumOfSqs[3]
)

prop.variance.SD <- rbind(prop.variance.SD, data.frame(
  tissue="seed", Fstat=adonis_seed$F[1], P=adonis_seed$`Pr(>F)`[1],
  prop.var.sp=adonis_seed$SumOfSqs[1]/adonis_seed$SumOfSqs[3],
  prop.var.resid=adonis_seed$SumOfSqs[2]/adonis_seed$SumOfSqs[3]
))

prop.variance.SD <- rbind(prop.variance.SD, data.frame(
  tissue="unripe", Fstat=adonis_unripe$F[1], P=adonis_unripe$`Pr(>F)`[1],
  prop.var.sp=adonis_unripe$SumOfSqs[1]/adonis_unripe$SumOfSqs[3],
  prop.var.resid=adonis_unripe$SumOfSqs[2]/adonis_unripe$SumOfSqs[3]
))

prop.variance.SD <- rbind(prop.variance.SD, data.frame(
  tissue="ripe", Fstat=adonis_ripe$F[1], P=adonis_ripe$`Pr(>F)`[1],
  prop.var.sp=adonis_ripe$SumOfSqs[1]/adonis_ripe$SumOfSqs[3],
  prop.var.resid=adonis_ripe$SumOfSqs[2]/adonis_ripe$SumOfSqs[3]
))


prop.variance.all <- rbind(prop.variance, prop.variance.SD)

write.csv(prop.variance.all, file="prop_variance_all.csv")



#Nice composite plot of compositional and structural beta diversity

ylab <- expression(paste("Distance to centroid"))
p <- ggplot(beta.div, aes(tissue, beta.div, fill=tissue)) + 
  geom_boxplot(show.legend = FALSE) +
  labs(x="", y=ylab) +
  labs(subtitle="(A)") +
  scale_fill_manual(values=pal) +
  scale_x_discrete(labels=c("leaf", "seed", "unripe\npulp", "ripe\npulp")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_blank()) +
  ylim(0.25,0.57) +
  annotate("text", x = 1:4, y = c(0.43, 0.56, 0.50, 0.48), label = c("A", "B", "B", "B")) 
p


p2 <- ggplot(s.beta.div, aes(tissue, beta.div, fill=tissue)) + 
  geom_boxplot(show.legend = FALSE) +
  labs(x="", y=ylab) +
  labs(subtitle="(B)") +
  scale_fill_manual(values=pal) +
  scale_x_discrete(labels=c("leaf", "seed", "unripe\npulp", "ripe\npulp")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylim(0.2,0.67) +
  annotate("text", x = 1:4, y = c(0.48, 0.64, 0.52, 0.51), label = c("A", "B", "AB", "B")) 
p2

cairo_pdf(height=5, width=3, family="sans",filename="betadiv.pdf")
grid.arrange(p, p2)
dev.off()

tiff(height=375, width=225, filename="betadiv.tiff", type="cairo") #0.8*widths
grid.arrange(p, p2)
dev.off()

postscript(height=5, width=3,family="sans", file="betadiv.eps")
grid.arrange(p, p2)
dev.off()


#------------------------------------------------------------------------
# Structural diversity
#--------------------------------------------------------------------

load("Workspace_PiperChem")

hist(SD$SD)
plot(SD$SD ~ SD$tissue)
plot(SD$SD ~ SD$species)


m1 <- lmer(SD ~ tissue * species + (1|PlantID), data=SD) 
#m1 <- lm(SD ~ tissue * species, data=SD) 
summary(m1)
#anova(m1, test="F")
drop1(m1, test="Chisq")

#singular fit errors go away when you drop random effect, but overall
#result of a strong tissue x species interaction is clear either way

interaction.plot(SD$tissue, SD$species, SD$SD)

p <- ggplot(SD, aes(tissue, SD)) + 
  geom_boxplot() +
  facet_wrap(vars(species))
p

#clear significant interaction between species and tissue, splitting by species

lmm.byspecies <- data.frame(species=character(), X=numeric(), P=numeric(), lf=character(), 
                             sd=character(), uf=character(), rf=character())
for(i in 1: length(levels(SD$species))){
  species <- levels(SD$species)[i]
  d <- SD[which(SD$species==species),]
  m <- lmer(SD~tissue + (1|PlantID), data=d)
  #m <- lm(SD~tissue, data=d)
  X <- drop1(m, test="Chisq")[2,3]
  P <- drop1(m, test="Chisq")[2,4]
  #X <- anova(m)[1,4]
  #P <- anova(m)[1,5]
  t <- cld(glht(m, mcp(tissue="Tukey")))
  l <- t$mcletters$Letters
  newrow <- data.frame(species=as.character(species), X=as.numeric(X), P=as.numeric(P), 
                       lf=as.character(l[1]), sd=as.character(l[2]), uf=as.character(l[3]),
                       rf=as.character(l[4]))
  lmm.byspecies <- rbind(lmm.byspecies, newrow)
}

#some boundary, singular fit errors here all go away when you run this as a simple lm
#without the plantID random effect. However, results are qualitatively similar
#so I will keep plant in

#pretty plot


lmm.byspecies$lf <- as.character(lmm.byspecies$lf)
lmm.byspecies$rf <- as.character(lmm.byspecies$rf)
lmm.byspecies$uf <- as.character(lmm.byspecies$uf)
lmm.byspecies$sd <- as.character(lmm.byspecies$sd)

colnames(lmm.byspecies)[4:7] <- c("leaf", "seed", "unripe", "ripe")

TukeyLab <- data.frame(species=character(), tissue=character(),x=character(), max=numeric())
for (i in 1:length(levels(SD$species))) {
  sp <- levels(SD$species)[i]
  for(j in 1:length(levels(SD$tissue))){
    t <- levels(SD$tissue)[j]
    newrow <- data.frame(species=sp, tissue=t, x=lmm.byspecies[which(lmm.byspecies$species==sp),
                                                          which(colnames(lmm.byspecies)==t)],
                         max=max(SD$SD[which(SD$species==sp & SD$tissue==t)]))
    TukeyLab <- rbind(TukeyLab, newrow)
  }
}

TukeyLab2 <- TukeyLab[-c(5:16,29:36),]

#ypos_lab <- TukeyLab2$max + 60 

TukeyLab2$xpos <- rep(c(1,2,3,4), 7)
TukeyLab2$ypos <- TukeyLab2$max + 0.1 

Plabs <- character()
for (i in 1:length(lmm.byspecies$sp)){
  sp <- lmm.byspecies$sp[i]
  if (lmm.byspecies$P[which(lmm.byspecies$sp==sp)] < 0.001){
    lab <- "P < 0.001"
  }else{
    lab <- paste("P ==", round(lmm.byspecies$P[i], 3)) 
  }
  Plabs <- c(Plabs, lab)
}

LRT_Lab <- data.frame(species=lmm.byspecies$sp, 
                      tissue=factor("leaf", levels=c("leaf", "seed", "unripe", "ripe")),
                      labs=sprintf("italic(X^2) == %.2f~%s", lmm.byspecies$X, Plabs),
                      xpos=c(rep(2.5, 4), 3, rep(2.5,7)), 
                      ypos=c(rep(0.9965, 4), 0.9998, rep(0.9965, 7)))

#see https://r-graphics.org/recipe-annotate-facet 
#for good examples on how to do this


p <- ggplot(SD, aes(tissue, SD, fill=tissue)) + 
  geom_boxplot() +
  ylim(0.9955, 1) +
  scale_fill_manual(values=pal, aesthetics = "fill") +
  labs(x=element_blank(), y="Structural Complexity")+
  facet_wrap(vars(species), ncol=4)+
  geom_text(data=TukeyLab2, mapping=aes(x=xpos, y=max+0.0007, label=x)) +
  geom_text(data=LRT_Lab, mapping=aes(x=xpos, y=ypos, label=labs), 
            parse=TRUE, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background =element_rect(fill=pal[4]))+
  theme(strip.text = element_text(colour = 'white')) +
  theme(legend.title=element_blank(), legend.position = "none") +
  theme (axis.text.x = element_text(angle = 45, vjust=0.65))
p


ggsave("Strucdiv.png", width = 18, height = 15, units = "cm")
ggsave("Strucdiv.eps", width = 18, height = 15, units = "cm", device=cairo_ps)
ggsave("Strucdiv.pdf", width = 18, height = 15, units = "cm", device = cairo_pdf)




#alternative plot--might have some bugs...changed factor level
#labels for SD and also sp to species

SD$tissue <- factor(SD$tissue, levels=c("lf","sd", "uf", "rf"))
lines <- data.frame(lines=c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5))

xlab <- substitute(paste(italic('Piper'), " species"))

lmm.byspecies$lf <- as.character(lmm.byspecies$lf)
lmm.byspecies$rf <- as.character(lmm.byspecies$rf)
lmm.byspecies$uf <- as.character(lmm.byspecies$uf)
lmm.byspecies$sd <- as.character(lmm.byspecies$sd)

TukeyLab <- data.frame(sp=character(), tissue=character(),x=character(), max=numeric())
for (i in 1:length(levels(SD$species))) {
  sp <- levels(SD$species)[i]
  for(j in 1:length(levels(SD$tissue))){
    t <- levels(SD$tissue)[j]
    newrow <- data.frame(sp=sp, tissue=t, x=lmm.byspecies[which(lmm.byspecies$species==sp),
                            which(colnames(lmm.byspecies)==t)],
                         max=max(SD$SD[which(SD$species==sp & SD$tissue==t)]))
    TukeyLab <- rbind(TukeyLab, newrow)
  }
}

TukeyLab2 <- TukeyLab[-c(5:16,29:36),]

ypos_lab <- TukeyLab2$max + 60 
xpos_lab <- c(0.7, 0.9, 1.1, 1.3, 
              #1.7, 1.9, 2.1, 2.3,
              #2.7, 2.9, 3.1, 3.3,
              #3.7, 3.9, 4.1, 4.3,
              4.7, 4.9, 5.1, 5.3,
              5.7, 5.9, 6.1, 6.3,
              6.7, 6.9, 7.1, 7.3,
              #7.7, 7.9, 8.1, 8.3,
              #8.7, 8.9, 9.1, 9.3,
              9.7, 9.9, 10.1, 10.3,
              10.7, 10.9, 11.1, 11.3,
              11.7, 11.9, 12.1, 12.3)
p <- ggplot(SD, aes(species, SD, fill=tissue)) + 
  geom_boxplot() +
  labs(x=xlab, y="Structural Diversity") +
  scale_fill_manual(values=pal, aesthetics = "fill") +
  geom_vline(data = lines, aes(xintercept = as.numeric(lines)), lty="dotted") +
  #geom_text(aes(y = ypos_lab, label = richness, group = tissue)) + 
  annotate("text", x = xpos_lab, y = ypos_lab, label = TukeyLab2$x) +
  theme(legend.title=element_blank()) +
  theme_bw()
p
ggsave("SD.tiff", width = 25, height = 6.5, units = "cm")
#ggsave("SD.eps", width = 8, height = 6.5, units = "cm", device=cairo_ps)
#ggsave("SD.pdf", width = 8, height = 6.5, units = "cm", device = cairo_pdf)



#-----------------------
###### OLDER STUFF###
#---------------------------


###Rarefaction function

# Rarefaction function
# Written by Will Wetzel for Wetzel and Whitehead (2020) Ecology Letters

library(compiler)
rare = function(mat, reps=50, maxSamples = NA) {
  N = nrow(mat)
  if(is.na(maxSamples) | maxSamples > nrow(mat)) 
    MS = nrow(mat) else MS = maxSamples
    
    dat = data.frame(samples = 1:MS, shan=NA, shanS=NA, shanLower=NA, shanUpper=NA, 
                     rich=NA, richS=NA, richLower10=NA, richUpper90=NA,
                     even=NA, evenLower10=NA, evenUpper90=NA)
    
    pb = txtProgressBar(min = 0, max = MS, initial = 0, style=3)
    
    for(n in 1:MS){
      # Sample units
      temp = replicate(reps, sample(N, n, replace = FALSE))
      if(n == 1) lumpedSamples = mat[temp,] else {
        lumpedSamples = t(apply(temp, 2, function(z) colSums(mat[z,]) ))
      }
      shannons = diversity(lumpedSamples)
      richnesses = specnumber(lumpedSamples) #apply(lumpedSamples, 1, function(z) sum(z > 0))
      evennesses = shannons/log(richnesses)
      # Save means
      dat[n,'shan'] = mean(shannons, na.rm=TRUE)
      dat[n,'shanS'] = sd(shannons, na.rm=TRUE)
      dat[n,'shanLower'] = quantile(shannons, probs=0.025, na.rm=TRUE)
      dat[n,'shanUpper'] = quantile(shannons, probs=0.975, na.rm=TRUE)
      dat[n,'rich'] = mean(richnesses)
      dat[n,'richS'] = sd(richnesses, na.rm=TRUE)
      dat[n,'richLower10'] = quantile(richnesses, probs=0.10, na.rm=TRUE)
      dat[n,'richUpper90'] = quantile(richnesses, probs=0.90, na.rm=TRUE)
      dat[n, 'even'] = mean(evennesses)
      dat[n,'evenLower10'] = quantile(evennesses, probs=0.10, na.rm=TRUE)
      dat[n,'evenUpper90'] = quantile(evennesses, probs=0.90, na.rm=TRUE)
      
      setTxtProgressBar(pb, n)
    }
    return(dat)
    
}

rare = cmpfun(rare)



#ALTERNATIVE to this constrained rarefaction is just to average across individuals first
#for each species. Then we only have 12 samples for each tissue, but each sample is independent.
#also, we can use more traditional rarefaction and get gamma estimates (e.g. jack1) for the asymptotes

#combining species data--this counts a compound as present in a species if it was found in 
#any of the three samples
div_SpAvg <- div %>%
  group_by(sp, tissue, tissue2) %>%
  summarize_at(-c(1:6), max)


# Peak areas only
comps_SpAvg <- div_SpAvg[,-c(1:3)]


out_SpAvg = data.frame(tiss = sort(div_SpAvg$tissue), rich=NA, richLower10=NA, richUpper90=NA)

for(i in sort(unique(div_SpAvg$tissue))) {
  tempd = subset(comps_SpAvg, div_SpAvg$tissue == i)
  tempo = rare(tempd, reps=50)    #increase reps to 5000 if we keep in paper--I reduced for speed
  
  out_SpAvg[out_SpAvg$tiss == i, 'samples'] = tempo$samples
  out_SpAvg[out_SpAvg$tiss == i, 'rich'] = tempo$rich
  out_SpAvg[out_SpAvg$tiss == i, 'richLower10'] = tempo$richLower10
  out_SpAvg[out_SpAvg$tiss == i, 'richUpper90'] = tempo$richUpper90
  out_SpAvg[out_SpAvg$tiss == i, 'rich1'] = tempo$rich[tempo$samples==1]
}


# Plot rarefaction
p <- ggplot(data=out_SpAvg) +
  geom_line(aes(x=samples, y=rich, group=tiss, color=tiss), size=1) +
  geom_ribbon(aes(x=samples, ymin=richLower10, ymax=richUpper90, group=tiss, fill=tiss), alpha=0.2)+
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  xlab('No. of samples') + ylab('Cumulative richness') +
  theme_classic() +
  theme(legend.position=c(0.8,0.6), legend.title=element_blank()) + 
  #plot.margin = margin(0.25,0.25,0.25,0.75, unit='lines')) +
  xlim(0,13) +
  ylim(0, 1400)
p


d.temp <- as.data.frame(div_SpAvg)
estRich_all <- data.frame(tiss=factor(levels=levels(d.temp$tissue)), Species=numeric(), 
                          chao=numeric(), chao.se=numeric(), jack1=numeric(),
                          jack1.se=numeric(), jack2=numeric(), boot=numeric(), 
                          boot.se=numeric(), n=numeric())
for(i in 1:length(levels(d.temp$tissue))) {
  tiss <- levels(d.temp$tissue)[i]
  d.temp2 <- d.temp[which(d.temp$tissue==tiss),]
  c <- specaccum(d.temp2[,-c(1:3)], method="random")
  estRich <- specpool(d.temp2[,-c(1:3)])
  row <- as.data.frame(c(tiss, estRich))
  colnames(row)[1] <- "tiss" 
  estRich_all <- rbind(estRich_all, row)
}
estRich_all

#fruits have higher diversity by all estimates


#One more option...to compare leaves generally to "fruit" we can use a random sampling
#of fruit samples. We have 3x as many samples from fruit, but we cannot use classic rarefaction 
#to correct for uneven sampling effort because the fact that the three fruit samples came from 
#the same fruit means that there is a correlated structure in the fruit data that would not
#be accounted for in classic rarefaction. This leads to an artificial pattern where we reach
#the fruit asymptote quickly in the accumulation curve (because once you have sampled most fruits
#once you are just repeating samples from the same fruit and the accumulation of new compounds
#slows way down).


##TO DO: modify sampling below so that when we take the random sample of fruit, we take one sample
#from each species--randomly choose a tissue, then choose that tissue for a species
#Also do this same sampling scheme for species accumulation curves; prob this needs a function
#that will give both the species accumulation curves and the estimates as an output


d.temp <- as.data.frame(div_SpAvg)
estRich_all <- data.frame(tiss=factor(levels=levels(d.temp$tissue2)), Species=numeric(), 
                          chao=numeric(), chao.se=numeric(), jack1=numeric(),
                          jack1.se=numeric(), jack2=numeric(), boot=numeric(), 
                          boot.se=numeric(), n=numeric())
s <- nrow(d.temp[which(d.temp$tissue2=="leaf"),])
for(i in 1:length(levels(d.temp$tissue2))) {
  tiss <- levels(d.temp$tissue2)[i]
  if(tiss=="frt"){
    d.temp2 <- d.temp[which(d.temp$tissue2=="frt"),]
    rows <- sample(nrow(d.temp2), s)
    d.temp2 <- d.temp2[rows,]
  } else {
    d.temp2 <- d.temp[which(d.temp$tissue2=="leaf"),]
  }
  c <- specaccum(d.temp2[,-c(1:3)], method="random")
  estRich <- specpool(d.temp2[,-c(1:3)])
  row <- as.data.frame(c(tiss, estRich))
  colnames(row)[1] <- "tiss" 
  estRich_all <- rbind(estRich_all, row)
}
estRich_all





#OLDER STUFF


#spatial rarefaction original code across all tissues

allCurves <- list()
means <- data.frame(samples=1:37)
for(i in 1:length(levels(div$tissue))){
  tiss <- levels(div$tissue)[i]
  d.temp <- div[which(div$tissue==tiss), -c(2, 4:9)]
  N <- nrow(d.temp)
  spAcumCurves <- data.frame(samples=1:N)
  for(j in 1:N){
    f <- d.temp$SampleID[j]
    sp <- as.character(d.temp$sp[which(d.temp$SampleID==f)])
    d.temp2 <- d.temp
    spOrd <- data.frame(sp=c(sp,sample(levels(d.temp2$sp)[which(levels(d.temp2$sp) !=sp)])), 
                        ord=1:length(levels(d.temp2$sp)))
    d.temp2 <- inner_join(d.temp2, spOrd, by="sp") 
    rows <- sample(nrow(d.temp2))
    d.temp2 <- d.temp2[rows,]
    d.temp2$ord[which(d.temp2$SampleID==f)] <- 0
    d.temp2 <- d.temp2[order(d.temp2$ord),]
    rich <- vector()
    for (k in 1:nrow(d.temp2)){
      d.temp3 <- d.temp2[1:k,]
      count <- length(which(colSums(d.temp3[-c(1:2,1314)]) !=0))
      rich[k] <- count
    }
    spAcumCurves <- cbind(spAcumCurves, rich) 
  }
  spAcumCurves$mean <- rowMeans(spAcumCurves[-1])
  allCurves[[i]] <- spAcumCurves
  names(allCurves)[i] <- tiss
  if(length(spAcumCurves$mean)==37){
    means <- cbind(means, spAcumCurves$mean)} else {
      means <- cbind(means, c(spAcumCurves$mean, NA))
    }
  colnames(means)[i+1] <- tiss
}

pal

p <- ggplot(data=means) +
  geom_line(aes(x=samples, y=leaf), col=pal[1], size=1) +
  geom_line(aes(x=samples, y=seed), col=pal[2], size=1) +
  geom_line(aes(x=samples, y=unripe), col=pal[3], size=1) +
  geom_line(aes(x=samples, y=ripe), col=pal[4], size=1) +
  labs(x='No. samples', y = 'Cumulative richness')
p

##spatial rarefaction with just lvs vs frt

allCurves2 <- list()
means2 <- data.frame(samples=1:109)
for(i in 1:length(levels(div$tissue2))){
  tiss <- levels(div$tissue2)[i]
  d.temp <- div[which(div$tissue2==tiss), -c(2, 4:9)]
  N <- nrow(d.temp)
  spAcumCurves <- data.frame(samples=1:N)
  for(j in 1:N){
    f <- d.temp$SampleID[j]
    sp <- as.character(d.temp$sp[which(d.temp$SampleID==f)])
    d.temp2 <- d.temp
    spOrd <- data.frame(sp=c(sp,sample(levels(d.temp2$sp)[which(levels(d.temp2$sp) !=sp)])), 
                        ord=1:length(levels(d.temp2$sp)))
    d.temp2 <- inner_join(d.temp2, spOrd, by="sp") 
    rows <- sample(nrow(d.temp2))
    d.temp2 <- d.temp2[rows,]
    d.temp2$ord[which(d.temp2$SampleID==f)] <- 0
    d.temp2 <- d.temp2[order(d.temp2$ord),]
    rich <- vector()
    for (k in 1:nrow(d.temp2)){
      d.temp3 <- d.temp2[1:k,]
      count <- length(which(colSums(d.temp3[-c(1:2,1314)]) !=0))
      rich[k] <- count
    }
    spAcumCurves <- cbind(spAcumCurves, rich) 
  }
  spAcumCurves$mean <- rowMeans(spAcumCurves[-1])
  allCurves2[[i]] <- spAcumCurves
  names(allCurves2)[i] <- tiss
  if(length(spAcumCurves$mean)==109){
    means2 <- cbind(means2, spAcumCurves$mean)} else {
      means2 <- cbind(means2, c(spAcumCurves$mean, rep(NA, 109-length(spAcumCurves$mean))))
    }
  colnames(means2)[i+1] <- tiss
}


p <- ggplot(data=means2) +
  geom_line(aes(x=samples, y=leaf), col=pal[1], size=1) +
  geom_line(aes(x=samples, y=frt), col=pal[4], size=1) +
  labs(x='No. samples', y = 'Cumulative richness')
p

###Hmnnnn...I don't think this works because the shape of the curves
#will be influenced by the "spatial" constraints. For the fruits, we do 
#a set of 9 related samples, then move on to the next species,
#this makes the curve more flat
#whereas for leaves, we do a set of 3 related samples, then move on
#to next species, so it is going to rise a bit faster I think
#probably best just to do each tissue individually


#Next using rarefaction, i.e. diversity accumulation curves for
#different tissue types, to examine gamma diversity across tissues


#first will examine diversity in fruit as a whole organ, across all fruit tissues

# Peak areas only
comps <- div[,-c(1:9)]

# Rarefaction

#takes a long time with 5000 reps (pretty much all day), can load
#previous workspace from run on 04/21/20
load("Workspace_PiperChem_Gamma")

out = data.frame(tiss = sort(div$tissue2), rich=NA, richLower10=NA, richUpper90=NA)

for(i in sort(unique(div$tissue2))) {
  tempd = subset(comps, div$tissue2 == i)
  tempo = rare(tempd, reps=50)    #increase reps to 5000 if we keep in paper--I reduced for speed
  
  out[out$tiss == i, 'samples'] = tempo$samples
  out[out$tiss == i, 'rich'] = tempo$rich
  out[out$tiss == i, 'richLower10'] = tempo$richLower10
  out[out$tiss == i, 'richUpper90'] = tempo$richUpper90
  out[out$tiss == i, 'rich1'] = tempo$rich[tempo$samples==1]
}

pal2 <- pal[1,4]

# Plot rarefaction
p <- ggplot(data=out) +
  geom_line(aes(x=samples, y=rich, group=tiss, color=tiss), size=1) +
  geom_ribbon(aes(x=samples, ymin=richLower10, ymax=richUpper90, group=tiss, fill=tiss), alpha=0.2)+
  scale_colour_manual(values=pal2) +
  scale_fill_manual(values=pal2) +
  xlab('No. of samples') + ylab('Cumulative richness') +
  theme_classic() +
  theme(legend.position=c(0.8,0.6), legend.title=element_blank()) + 
  #plot.margin = margin(0.25,0.25,0.25,0.75, unit='lines')) +
  xlim(0,110) +
  ylim(0, 1400)
p

#ggsave("Rarefaction_lfvsfrt.tiff", width = 8, height = 6.5, units = "cm")
ggsave("Rarefaction_lfvsfrt.png", width = 8, height = 6.5, units = "cm")
#ggsave("Rarefaction_lfvsfrt.eps", width = 8, height = 6.5, units = "cm", device=cairo_ps)
#ggsave("Rarefaction_lfvsfrt.pdf", width = 8, height = 6.5, units = "cm", device = cairo_pdf)




#Repeating rarefaction, but including all tissues individually 

# Rarefaction by tissue type
out.allt = data.frame(tiss = sort(div$tissue), rich=NA, richLower10=NA, richUpper90=NA)


for(i in sort(unique(div$tissue))) {
  tempd = subset(comps, div$tissue == i)
  tempo = rare(tempd, reps=5000)    #increase reps to 5000 if we keep in paper--I reduced for speed
  
  out.allt[out.allt$tiss == i, 'samples'] = tempo$samples
  out.allt[out.allt$tiss == i, 'rich'] = tempo$rich
  out.allt[out.allt$tiss == i, 'richLower10'] = tempo$richLower10
  out.allt[out.allt$tiss == i, 'richUpper90'] = tempo$richUpper90
  out.allt[out.allt$tiss == i, 'rich1'] = tempo$rich[tempo$samples==1]
}


# Plot rarefaction by tissue

#re-ordering factor levels--won't need this if I re-run rarefaction with
#new data wrangling code above as I did this in the start but I haven't 
#re-run the 5000 rep rarefaction with it
out.allt2 <- out.allt 
out.allt2$tiss <- factor(out.allt2$tiss, levels=c("leaf","seed", "unripe", "ripe"))


p <- ggplot(data=out.allt2) +
  geom_line(aes(x=samples, y=rich, group=tiss, color=tiss), size=1) +
  geom_ribbon(aes(x=samples, ymin=richLower10, ymax=richUpper90, group=tiss, fill=tiss), alpha=0.2)+
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  xlab('No. of samples') + ylab('Cumulative richness') +
  theme_classic() +
  theme(legend.position=c(0.8,0.5), legend.title=element_blank()) +  
  #plot.margin = margin(0.25,0.25,0.25,0.75, unit='lines')) +
  xlim(1,40) +
  ylim(0, 1300)
p

ggsave("Rarefaction_alltiss.tiff", width = 8, height = 6.5, units = "cm")
#ggsave("Rarefaction_alltiss.eps", width = 8, height = 6.5, units = "cm", device=cairo_ps)
#ggsave("Rarefaction_alltiss.pdf", width = 8, height = 6.5, units = "cm", device = cairo_pdf)


#to save workspace after a rarefaction run, last save on 04/21/20 with
#5000 reps
#save.image("Workspace_PiperChem_Gamma")




load("Workspace_PiperChem")

#can also try this with the species averaged first
div_SpAvg <- div %>%
  group_by(sp, tissue, tissue2) %>%
  summarize_at(-c(1:6), max)


#another option for species "averaging", could we count the number
#of samples where we detected the compound as the "count"

div_SpAvg <- div %>%
  group_by(sp, tissue, tissue2) %>%
  summarize_at(-c(1:6), ~sum(. != 0))


# Peak areas only
comps_SpAvg <- div_SpAvg[,-c(1:3)]


out_SpAvg = data.frame(tiss = sort(div_SpAvg$tissue2), rich=NA, richLower10=NA, richUpper90=NA)

for(i in sort(unique(div_SpAvg$tissue2))) {
  tempd = subset(comps_SpAvg, div_SpAvg$tissue2 == i)
  tempo = rare(tempd, reps=50)    #increase reps to 5000 if we keep in paper--I reduced for speed
  
  out_SpAvg[out_SpAvg$tiss == i, 'samples'] = tempo$samples
  out_SpAvg[out_SpAvg$tiss == i, 'rich'] = tempo$rich
  out_SpAvg[out_SpAvg$tiss == i, 'richLower10'] = tempo$richLower10
  out_SpAvg[out_SpAvg$tiss == i, 'richUpper90'] = tempo$richUpper90
  out_SpAvg[out_SpAvg$tiss == i, 'rich1'] = tempo$rich[tempo$samples==1]
}


# Plot rarefaction
p <- ggplot(data=out_SpAvg) +
  geom_line(aes(x=samples, y=rich, group=tiss, color=tiss), size=1) +
  geom_ribbon(aes(x=samples, ymin=richLower10, ymax=richUpper90, group=tiss, fill=tiss), alpha=0.2)+
  scale_colour_manual(values=pal2) +
  scale_fill_manual(values=pal2) +
  xlab('No. of samples') + ylab('Cumulative richness') +
  theme_classic() +
  theme(legend.position=c(0.8,0.6), legend.title=element_blank()) + 
  #plot.margin = margin(0.25,0.25,0.25,0.75, unit='lines')) +
  xlim(0,50) +
  ylim(0, 1400)
p


#This does not really make sense because we have different sampling dependency
#the addtional fruit samples are not independent from one another (3 tissues from
#same species, so it is going to make that curve artificially flatten
#out once each species is sampled)


##and with species averaged but tissues separate

out_SpAvg_allt = data.frame(tiss = sort(div_SpAvg$tissue), rich=NA, richLower10=NA, richUpper90=NA)

for(i in sort(unique(div_SpAvg$tissue))) {
  tempd = subset(comps_SpAvg, div_SpAvg$tissue == i)
  tempo = rare(tempd, reps=50)    #increase reps to 5000 if we keep in paper--I reduced for speed
  
  out_SpAvg_allt[out_SpAvg_allt$tiss == i, 'samples'] = tempo$samples
  out_SpAvg_allt[out_SpAvg_allt$tiss == i, 'rich'] = tempo$rich
  out_SpAvg_allt[out_SpAvg_allt$tiss == i, 'richLower10'] = tempo$richLower10
  out_SpAvg_allt[out_SpAvg_allt$tiss == i, 'richUpper90'] = tempo$richUpper90
  out_SpAvg_allt[out_SpAvg_allt$tiss == i, 'rich1'] = tempo$rich[tempo$samples==1]
}


# Plot rarefaction by tissue


p <- ggplot(data=out_SpAvg_allt) +
  geom_line(aes(x=samples, y=rich, group=tiss, color=tiss), size=1) +
  geom_ribbon(aes(x=samples, ymin=richLower10, ymax=richUpper90, group=tiss, fill=tiss), alpha=0.2)+
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  xlab('No. of samples') + ylab('Cumulative richness') +
  theme_classic() +
  theme(legend.position=c(0.8,0.5), legend.title=element_blank()) +  
  #plot.margin = margin(0.25,0.25,0.25,0.75, unit='lines')) +
  xlim(1,15) +
  ylim(300, 1300)
p




