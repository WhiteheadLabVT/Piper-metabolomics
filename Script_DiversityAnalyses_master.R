
library(dplyr) #version 0.8.3
library(lme4) #version 1.1-21
library(ggplot2) #version 3.3.5
library(gridExtra) #version 2.3
library(multcomp) #version 1.4-10
library(vegan)  #version 2.5-5 functions are buggy with later versions
library(lmPerm) #version 2.1.0
library(viridis) #version 0.5.1
library(scales) #version 1.0.0
library(VennDiagram) #version 1.6.20


#-----------------------------------------------------
#   Data wrangling--run this first and save
#   workspace, then all other sections can be run independently
#-----------------------------------------------------


#reading in peak tables
div_all <- read.csv("Data_Peak_Table.csv")


#replacing all values < 100 with zero.
#This is based on examining noise levels in chromatograms. 
div <- div_all
for (i in 5:1315){
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


#calculating richness
div <- mutate(div, richness=rowSums(div[5:1315]>100))

#ordering columns
div <- div[c(1:4, 1316:1317, 5:1315)]

div$sp <- as.factor(div$sp)
div$tissue2 <- as.factor(div$tissue2)
div$PlantID <- as.factor(div$PlantID)

#read in structural similarity data, intrasample data (all compounds in a sample compared
# to each other)

SD <- read.csv("Data_Intrasample_Structural_Similarity.csv")
SD <- SD[-c(146:149),]


#similarity metric shows the chemical similarity between samples, but we would like to look at 
#diversity, so I will do 1-similarity
SD$SD <- 1-SD$chem_similarity_internal



###read in structural diversity data, intersample--overall structural composition of each sample
##compared to all other samples in a distance matrix

SD_inter <- read.csv("Data_Intersample_Structural_Similarity.csv", row.names=1)

#data are similarity scores, and I want diversity, so I will take 1-sim
SD_inter <- 1-SD_inter


#removing two samples, col_12 (ripe fruit) and pel_1 (unripe fruits) that are from plants with no other tissues
SD_inter <- SD_inter[,-which(colnames(SD_inter)=="col_rf_29_180510")]
SD_inter <- SD_inter[,-which(colnames(SD_inter)=="pel_uf_52_180511")]
SD_inter <- SD_inter[-which(rownames(SD_inter)=="col_rf_29_180510"),]
SD_inter <- SD_inter[-which(rownames(SD_inter)=="pel_uf_52_180511"),]


#this also needs to have a set of explanatory variables for some analyses

SD_inter_expl <- data.frame(SampleID=rownames(SD_inter))
SD_inter_expl <- left_join(SD_inter_expl, SD[1:4], by='SampleID')



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



#-------------------------------------------------------------
#  Overall PERMANOVAs to assess differences in composition across samples
#--------------------------------------------------

load("Workspace_PiperChem")

d.temp <- div[-c(1:6)]
d.temp.expl <- div[c(2:5)]

m1 <- adonis2(d.temp ~ tissue*sp, strata = 'PlantID', data = d.temp.expl, method = "bray", binary=TRUE, permutations = 999)
m1
#strong effects of tissue, species, and their interaction

## Note, in newest version of vegan strata is deprecated, so you need to set up using permutations argument
#perm <- how(nperm = 999)
#setBlocks(perm) <- with(d.temp.expl, PlantID)
#m1 <- adonis2(d.temp ~ tissue*sp, data = d.temp.expl, method = "bray", binary=TRUE, permutations = perm)



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



#now we will split by species and look for tissue level differences for each species
#for these we cannot do the pairwise contrasts between each pair of tissue types because once we 
#get to that level with N=3 per species and strata = plantID, there are only 9 possible permutations,
#so the minimum p is 0.1

adonis.byspecies <- data.frame(sp=character(), Fstat=numeric(), P=numeric())

for(i in 1: length(levels(div$sp))){
  sp <- levels(div$sp)[i]
  d <- div[which(div$sp==sp),]
  d.temp <- d[-c(1:9)]
  d.temp.expl <- d[1:9]
  m <- adonis2(d.temp ~ tissue, data=d.temp.expl, 
               method="bray", binary=TRUE)
  Fstat <- m$F[1]
  P <- m$`Pr(>F)`[1]
  newrow <- data.frame(sp=as.character(sp), Fstat=as.numeric(Fstat), P=as.numeric(P))
  adonis.byspecies <- rbind(adonis.byspecies, newrow)
}

adonis.byspecies
#always see an overall effect of tissue





##Now trying a PERMANOVA on the structural data

#in this case, we feed adonis2 our own distance matrix rather than the dataframe
d.temp <- as.dist(SD_inter)
d.temp.expl <- SD_inter_expl
m2 <- adonis2(d.temp ~ tissue*species, strata='PlantID', data = d.temp.expl, permutations = 999)
m2
plot(d.temp)

m2.pw <- pairwise.adonis2(d.temp ~ tissue*species, strata='PlantID', data = d.temp.expl, 
                          nperm=999)
m2.pw


#now trying to split by species and look for tissue level differences
adonis.byspecies.SD <- data.frame(sp=character(), Fstat=numeric(), P=numeric())


for(i in 1: length(levels(div$sp))){
  sp <- levels(div$sp)[i]
  d <- div[which(div$sp==sp),]
  d.temp <- d[-c(1:9)]
  d.temp.expl <- d[1:9]
  m <- adonis2(d.temp ~ tissue, data=d.temp.expl, 
               method="bray", binary=TRUE)
  Fstat <- m$F[1]
  P <- m$`Pr(>F)`[1]
  newrow <- data.frame(sp=as.character(sp), Fstat=as.numeric(Fstat), P=as.numeric(P))
  adonis.byspecies.SD <- rbind(adonis.byspecies.SD, newrow)
}


#combined by species data for table

tbl <- cbind(adonis.byspecies[1:3], adonis.byspecies.SD[2:3])

write.csv(tbl, file="./Outputs/Table2_adonis.byspecies.csv")







#-------------------------------------------------------------
#  Gamma diversity
#-----------------------------------------------------------

load("Workspace_PiperChem")
#load("Workspace_PiperChem_Gamma")  #last rarefaction run with 5000 reps

#-----------Venn------------------
#First, will show a venn diagram with the total count of compounds detected in each
#tissue type and shared among tissues


levels(div$tissue)
d <- div[which(div$tissue=="leaf"),7:1317]
lf <- colnames(d)[which(colSums(d) !=0)]

d <- div[which(div$tissue=="ripe"),7:1317]
r <- colnames(d)[which(colSums(d) !=0)]

d <- div[which(div$tissue=="seed"),7:1317]
s <- colnames(d)[which(colSums(d) !=0)]

d <- div[which(div$tissue=="unripe"),7:1317]
u <- colnames(d)[which(colSums(d) !=0)]

x <- list(lf, s, u, r)


#make a list of compounds shared in each category of overlap, e.g. in everything, in 
#seeds and leaves, etc.
overlap <- calculate.overlap(x)
 

#use this to double check groupings...can see data for a particular compund and which
#tissues it was found in
div[,c("tissue","X577.1393_3.1992")]

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
#png(height=500, width=500, filename="./Outputs/Fig1_Venn.tiff", type="cairo")
cairo_pdf(height=7, width=7, filename="./Outputs/Fig1_Venn.pdf")
#cairo_ps(height=7, width=7, filename="./Outputs/Fig1_Venn.eps")
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
#compounds <- div[which(div$tissue=="leaf"), -c(1:6)]
#meta <- div[which(div$tissue=="leaf"), c(1,2)]
#reps <- 100
#mod=c("arrhenius","gitay","lomolino","asymp", 
 #          "gompertz", "logis")
      

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
      pivot_longer(cols=everything()) %>% 
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
  d.temp <- div[which(div$tissue==tiss), -c(1:6)]
  m.temp <- div[which(div$tissue==tiss), c(1,2)]
  CPR_out <- CPR(d.temp, m.temp, reps=5000) 
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
write.csv(asymp.finalest, file="./Outputs/Table3_gammaestimates.csv")

#save.image("Workspace_PiperChem_Gamma")



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

ggsave("./Outputs/Figure1_Rarefaction.png", width = 8, height = 6.5, units = "cm")
ggsave("./Outputs/Figure1_Rarefaction.eps", width = 8, height = 6.5, units = "cm", device=cairo_ps)
ggsave("./Outputs/Figure1_Rarefaction.pdf", width = 8, height = 6.5, units = "cm", device = cairo_pdf)


#-------------------------------------------------------------------
#  Alpha diversity--looking at average richness in each sample type
#------------------------------------------------------------------------


load("Workspace_PiperChem")


hist(div$richness)

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

#Some boundary, singular fit warnings. These all go away when you run this as a simple lm
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


ggsave("./Outputs/Figure4_Alphadiv.png", width = 15, height = 15, units = "cm")
ggsave("./Outputs/Figure4_Alphadiv.eps", width = 15, height = 15, units = "cm", device=cairo_ps)
ggsave("./Outputs/Figure4_Alphadiv.pdf", width = 15, height = 15, units = "cm", device = cairo_pdf)


#-----------------------------------------------------------------
#  Beta-diversity
#------------------------------------------------------------------

load("Workspace_PiperChem")

#------------Chemical composition-----------
#First looking at beta-diversity (i.e. sample-to-sample variance) 
#in chemical composition, based on compound presence/absence

d.temp <- div[-c(1:6)]

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


#Data for plot
beta.div <- data.frame(tissue=wb$group, beta.div=wb$distances)
beta.div$tissue <- factor(beta.div$tissue, levels=c("leaf","seed", "unripe", "ripe"))


#can also look at differences for species
wb_sp <- betadisper(w, div$sp)
anova(wb_sp)
TukeyHSD(wb_sp)
boxplot(wb_sp)
#generalense seems to have higher beta-diversity
#than other species. Note all the tissues are in there, so this probably is driven
#by major differences among tissues


##there are effects of tissue, and effects of species, but we can't include
#two factors simultaneously in betadisper. So for the differences in 
#beta-diversity across tissues above, this might be due to within species 
#variance in addition to across species variance. 

#We could look at beta-diversity across tissues with the species averaged first

div_SpAvg <- div %>%
  group_by(sp, tissue, tissue2) %>%
  summarize_if(is.numeric, max)

d.temp <- div_SpAvg[-c(1:5)]
w <- vegdist(d.temp, binary=TRUE, method="bray")
plot(w)
wb <- betadisper(w, div_SpAvg$tissue)
anova(wb)
TukeyHSD(wb)
boxplot(wb)
#same trend but not significant, likely because our sample size is only 12 now
#I also do not like this because it is ignoring an important scale of 
#variation (within-species) for beta diversity


##In this paper:
#Anderson, Marti J. (2014) "Permutational multivariate analysis of variance (PERMANOVA)."
#See Table 2 and Fig. 4...they show an example where PERMANOVA is used to
#partition the variance in multivariate composition among different spatial 
#scales. Following this model, we can show both the within species and across species
#components of variation

## Using adonis to quantify the contribution of species to tissue-level variance
# subset data by tissue and run adonis with species as the independent variable
# Use Sum of Squares of residuals = SSE to calculate sigma^2 = SSE/(n-v)

d.temp <- div[which(div$tissue=="leaf"), -c(1:6)]
d.temp.expl <- div[which(div$tissue=="leaf"), c(2,4)]
adonis_leaf <- adonis2(d.temp ~ sp, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_leaf


d.temp <- div[which(div$tissue=="ripe"), -c(1:6)]
d.temp.expl <- div[which(div$tissue=="ripe"), c(2,4)]
adonis_ripe <- adonis2(d.temp ~ sp, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_ripe

d.temp <- div[which(div$tissue=="unripe"), -c(1:6)]
d.temp.expl <- div[which(div$tissue=="unripe"), c(2,4)]
adonis_unripe <- adonis2(d.temp ~ sp, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_unripe

d.temp <- div[which(div$tissue=="seed"), -c(1:6)]
d.temp.expl <- div[which(div$tissue=="seed"), c(2,4)]
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

prop.variance
#"prop.var.sp" shows the proportion of variance explained by species
#"prop.var.resid" shows the proportion of variance explained by
#individual sample within species

#------------Structural composition-----------
#Now looking at structural beta-diversity (i.e. sample-to-sample variance) 
#in structural composition, based on intersample distances in structural 
#similarity cosine scores

d.temp <- as.dist(SD_inter)
plot(d.temp)

#then we can look at the distances between each sample and the group centroid
#but we need to define a specific group (e.g. tissue)
wb <- betadisper(d.temp, SD_inter_expl$tissue)
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

d.temp <- as.dist(SD_inter[which(SD_inter_expl$tissue=="leaf"), 
                           which(SD_inter_expl$tissue=="leaf")])
d.temp.expl <- SD_inter_expl[which(SD_inter_expl$tissue=="leaf"),]
adonis_leaf <- adonis2(d.temp ~ species, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_leaf


d.temp <- as.dist(SD_inter[which(SD_inter_expl$tissue=="seed"), 
                           which(SD_inter_expl$tissue=="seed")])
d.temp.expl <- SD_inter_expl[which(SD_inter_expl$tissue=="seed"),]
adonis_seed <- adonis2(d.temp ~ species, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_seed

d.temp <- as.dist(SD_inter[which(SD_inter_expl$tissue=="unripe"), 
                           which(SD_inter_expl$tissue=="unripe")])
d.temp.expl <- SD_inter_expl[which(SD_inter_expl$tissue=="unripe"),]
adonis_unripe <- adonis2(d.temp ~ species, strata='PlantID', data = d.temp.expl, method = "bray", permutations = 999)
adonis_unripe

d.temp <- as.dist(SD_inter[which(SD_inter_expl$tissue=="ripe"), 
                           which(SD_inter_expl$tissue=="ripe")])
d.temp.expl <- SD_inter_expl[which(SD_inter_expl$tissue=="ripe"),]
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

write.csv(prop.variance.all, file="./Outputs/Table4_prop_variance_all.csv")



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

cairo_pdf(height=5, width=3, family="sans",filename="./Outputs/Figure6_betadiv.pdf")
grid.arrange(p, p2)
dev.off()

tiff(height=375, width=225, filename="./Outputs/Figure6_betadiv.tiff", type="cairo") #0.8*widths
grid.arrange(p, p2)
dev.off()

postscript(height=5, width=3,family="sans", file="./Outputs/Figure6_betadiv.eps")
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


ggsave("./Outputs/Figure5_Strucdiv.png", width = 18, height = 15, units = "cm")
ggsave("./Outputs/Figure5_Strucdiv.eps", width = 18, height = 15, units = "cm", device=cairo_ps)
ggsave("./Outputs/Figure5_Strucdiv.pdf", width = 18, height = 15, units = "cm", device = cairo_pdf)



