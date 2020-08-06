library(dplyr)
library(lme4)
library(ggplot2)
library(multcomp)
library(vegan)
library(viridis)
library(scales)
library(VennDiagram)
install.packages("remotes")
library(remotes)
remotes::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

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

#making plant ID a factor
div$PlantID <- as.factor(div$PlantID)

#correcting spelling
div$sp <- recode(div$sp, multiplinervium="multiplinervum")

#making sure all plant IDs are unique
div$PlantID <- paste(div$sp, div$PlantID, sep="_")

#saving current div as the div_all version
div_all <- div

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

#also droping an extra sample for colonense for ripe fruit where we don't have the other tissues
#could also keep this in for some analyses, but it messes with model fits where we are trying to
#block for individual and also the constrained rarefaction

div <- div[which(div$PlantID !="colonense_12"),]

#also will create a new variable with all fruit tissues recoded as just fruit,
#so we can examine diversity in fruit as a whole organ in some cases
div$tissue2 <- recode(div$tissue, ripe="frt", unripe="frt", seed="frt")
div <- div[c(1:8, 1320, 9:1319)]



#structural diversity data, intrasample data (all compounds in a sample compared
# to each other)

SD <- read.csv("chem_structural_similarity_intrasample.csv")

#metric shows the "chemical similarity between samples, but I would like to look at 
#diversity, so I will take the inverse
SD$SD <- 1/SD$chem_similarity_internal

#adding sample ID info to SD dataset
d.temp <- IDs[,c(1,5)]
colnames(SD)[3] <- "extraction"
SD$extraction <- as.character(SD$extraction)

SD <- full_join(SD, d.temp, by="extraction")

SD <- SD[which(!is.na(SD$run_date)),]



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



##setting color palette for all figures
##viridis colors
show_col(viridis_pal()(4))

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
  CPR_out <- CPR(d.temp, m.temp, reps=50) 
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

allCoefEst[which(allCoefEst$model_type=="asymp" & allCoefEst$parameter=="Asym"),]

#leaf is lower than all other fruit tissues, 95% CIs do not cross. 

# Plot rarefaction

#ordering factor levels so it plots in this order on graphs
allCurves$tiss <- factor(allCurves$tiss, levels=c("leaf","seed", "unripe", "ripe"))
p <- ggplot(data=allCurves) +
  geom_line(aes(x=samples, y=mean, group=tiss, color=tiss), size=1) +
  geom_ribbon(aes(x=samples, ymin=CI_low, ymax=CI_high, group=tiss, fill=tiss), alpha=0.2) +
  scale_colour_manual(values=pal) +
  scale_fill_manual(values=pal) +
  xlab('No. of samples') + ylab('Cumulative richness') +
  theme_classic() +
  theme(legend.position=c(0.8,0.5), legend.title=element_blank()) + 
  #plot.margin = margin(0.25,0.25,0.25,0.75, unit='lines')) +
  xlim(0,37) +
  ylim(0, 1400)
p


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

#some boundary, singular fit errors here all go away when you run this as a simple lm
#without the plantID random effect. However, results are qualitatively similar
#so I will keep plant in


#pretty plot

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
ggsave("Alphadiv.png", width = 25, height = 10, units = "cm")
#ggsave("Alphadiv.eps", width = 30, height = 7, units = "cm", device=cairo_ps)
#ggsave("BetadivAlphadiv.pdf", width = 30, height = 7, units = "cm", device = cairo_pdf)



##Alternative plot

TukeyLab2$xpos <- rep(c(1,2,3,4), 10)

p <- ggplot(div, aes(tissue, richness, fill=tissue)) + 
  geom_boxplot() +
  scale_fill_manual(values=pal, aesthetics = "fill") +
  labs(x=element_blank(), y="Compound Richness")+
  facet_wrap(vars(sp), ncol=6)+
  geom_text(data=TukeyLab2, mapping=aes(x=xpos, y=max+50, label=lab)) +
  theme_bw() +
  theme(legend.title=element_blank(), legend.position = "none")
p


ggsave("Alphadiv2.png", width = 25, height = 15, units = "cm")






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

#then we can look at the distances between each sample and the group centroid
#but we need to define a specific group (e.g. tissue)
wb <- betadisper(w, div$tissue)
anova(wb)
TukeyHSD(wb)
boxplot(wb)


#note that the results are very different if we use the quantitative data
w2 <- vegdist(d.temp, method="bray")
plot(w2)

#an alternative to the anova is to perform a permutation test
permutest(wb) 


#can also look at differences for species
wb_sp <- betadisper(w, div$sp)
anova(wb_sp)
TukeyHSD(wb_sp)
boxplot(wb_sp)
#main difference is that generalense seems to have higher beta-diversity
#than other species. But this is probably not that useful
#considering all the tissues are in there


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
#same trend but not significant


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
#So it seems the main tissue beta-diversity effects are driven by species
#level differences

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

#There is a strong interaction here. 

w


##Also just found this paper:
#Anderson, Marti J. (2014) "Permutational multivariate analysis of variance (PERMANOVA)."
#See Table 2 and Fig. 4 in that paper...they show an example where PERMANOVA is used to
#partition the variance in multivariate composition among different spatial scales
#something like this would be perfect, to show both the within species and across species
#components of variation

## Using adonis to quantify the contribution of species to tissue-level variance
# subset mz matrix and IDs matrix by tissue and run adonis with species as the independent variable
# Use Sum of Squares of residuals = SSE to calculate sigma^2 = SSE/(n-v)

div_matrix_bytissue <- div_matrix
div_matrix_bytissue$tissue <- IDs_matrix$tissue
div_matrix_leaf <- filter(div_matrix_bytissue, tissue == "leaf")
div_matrix_leaf <- select(div_matrix_leaf, -tissue)
IDs_matrix_leaf <- filter(IDs_matrix, tissue == "leaf")

adonis_leaf <- adonis2(div_matrix_leaf ~ sp, data = IDs_matrix_leaf, method = "bray", permutations = 9999)

div_matrix_bytissue$tissue <- IDs_matrix$tissue
div_matrix_ripe <- filter(div_matrix_bytissue, tissue == "ripe")
div_matrix_ripe <- select(div_matrix_ripe, -tissue)
IDs_matrix_ripe <- filter(IDs_matrix, tissue == "ripe")

adonis_ripe <- adonis2(div_matrix_ripe ~ sp, data = IDs_matrix_ripe, method = "bray", permutations = 9999)

div_matrix_bytissue$tissue <- IDs_matrix$tissue
div_matrix_unripe <- filter(div_matrix_bytissue, tissue == "unripe")
div_matrix_unripe <- select(div_matrix_unripe, -tissue)
IDs_matrix_unripe <- filter(IDs_matrix, tissue == "unripe")

adonis_unripe <- adonis2(div_matrix_unripe ~ sp, data = IDs_matrix_unripe, method = "bray", permutations = 9999)

div_matrix_bytissue$tissue <- IDs_matrix$tissue
div_matrix_seed <- filter(div_matrix_bytissue, tissue == "seed")
div_matrix_seed <- select(div_matrix_seed, -tissue)
IDs_matrix_seed <- filter(IDs_matrix, tissue == "seed")

adonis_seed <- adonis2(div_matrix_seed ~ sp, data = IDs_matrix_seed, method = "bray", permutations = 9999)



beta.div <- data.frame(tissue=wb$group, beta.div=wb$distances)
beta.div$tissue <- factor(beta.div$tissue, levels=c("leaf","ripe seed", "unripe", "ripe"))

ylab <- expression(paste("Beta diversity (", omega, ")"))
p <- ggplot(beta.div, aes(tissue, beta.div, fill=tissue)) + 
  geom_boxplot(show.legend = FALSE) +
  labs(x="", y=ylab) +
  scale_fill_manual(values=pal) +
  scale_x_discrete(labels=c("leaf", "seed", "unripe\npulp", "ripe\npulp")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylim(0.25,0.57) +
  annotate("text", x = 1:4, y = c(0.43, 0.56, 0.50, 0.48), label = c("A", "B", "B", "B")) 
p
ggsave("Betadiv.tiff", width = 8, height = 6.5, units = "cm")
ggsave("Betadiv.eps", width = 8, height = 6.5, units = "cm", device=cairo_ps)
ggsave("Betadiv.pdf", width = 8, height = 6.5, units = "cm", device = cairo_pdf)



#------------------------------------------------------------------------
# Structural diversity
#--------------------------------------------------------------------

hist(SD$SD)
plot(SD$SD ~ SD$tissue)
plot(SD$SD ~ SD$species)


m1 <- lmer(SD ~ tissue * species + (1|PlantID), data=SD) 
#m1 <- lm(SD ~ tissue * species, data=SD) 
summary(m1)
anova(m1, test="F")
drop1(m1, test="Chisq")

interaction.plot(SD$tissue, SD$species, SD$SD)

p <- ggplot(SD, aes(tissue, SD)) + 
  geom_boxplot() +
  facet_wrap(vars(species))
p

#clear significant interaction between speciesecies and tissue, splitting by species

lmm.byspecies <- data.frame(species=character(), X=numeric(), P=numeric(), lf=character(), 
                             sd=character(), uf=character(), rf=character())
SD$tissue <- factor(SD$tissue, levels=c("lf","sd", "uf", "rf"))
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

##viridis colors
show_col(viridis_pal()(4))
#leaves=green, seeds=yellow, ripe=purple, unripe=blue 

pal <- viridis_pal()(4)[c(3,4,2,1)]

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
