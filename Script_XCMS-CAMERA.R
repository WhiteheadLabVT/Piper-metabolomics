######################## Pre-processing, deconvolution, and annotation of LC-MS data with XCMS and CAMERA 

### on first use of XCMS
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("xcms", "MSnbase", "limma", "CAMERA"))

### on subsequent uses, to check for updates
BiocManager::install()

library(multtest)
library(MSnbase)
library(xcms)
library(snow)
library(muma)
library(pvclust)
library(data.table)
library(vegan)
library(splitstackshape)
library(CAMERA)

imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)
                     
# the following code assumes samples are located in directory listed below, in separate folders labeled "Blank" and "Sample"
setwd("path to folder containing CDF files to be analyzed")

xcmsSet(polarity = "positive")
files_1 <- list.files(recursive = TRUE, full.names = TRUE)
files <- files_1[endsWith(files_1, "CDF")]
s_groups <- sapply(1:length(files), function(x) unlist(strsplit(files[x], split = "/"))[2])
s_names <- sapply(1:length(files), function(x) sub(unlist(strsplit(files[x], split = "/"))[3], pattern = ".CDF", replacement = "", fixed = TRUE))

cwparam <- CentWaveParam(ppm=15, peakwidth=c(4,12), snthresh=20, prefilter=c(10,1000)) #p.149 in manual
#increased snthresh to 20 and prefilter to 10,1000 from previous version
dparam1 <- PeakDensityParam(sampleGroups = s_groups, bw=10, binSize=0.05, minSamples = 1, minFraction = 0.01) #p.126
oparam <- ObiwarpParam(binSize=1) # bin size in m/z for alignment profiling
dparam2 <- PeakDensityParam(sampleGroups = s_groups, bw = 3, binSize = 0.025, minSamples = 1, minFraction = 0.01)
fpparam <- FillChromPeaksParam(expandMz = 0, expandRt = 0, ppm = 0) #p.64
pheno <- data.frame(sample_name = s_names, sample_group = s_groups, stringsAsFactors = FALSE)
raw_data <- readMSData(files, msLevel. = 1, mode="onDisk")

xcmsexp <- findChromPeaks(object = raw_data, param = cwparam)
xcmsexp <- groupChromPeaks(object = xcmsexp, param = dparam1) # first grouping = broader scale, high bandwidth and m/z bin size
xcmsadjustrtime <- adjustRtime(object = xcmsexp, param = oparam)
xcmsexp <- groupChromPeaks(object = xcmsexp, param = dparam2) # second grouping = fine scale, low bandwidth and m/z bin size
xcmsexp <- fillChromPeaks(xcmsexp, fpparam)
xset <- as(xcmsexp, "xcmsSet")
sampclass(xset) <- sapply(1:length(files), function(x) unlist(strsplit(files[x], split = "/"))[2])
## applying CAMERA functions
xset1 <- xsAnnotate(xs=xset,polarity="positive")
xset2 <- groupFWHM(xset1, perfwhm=0.7) # CAMERA p.30
xset3 <- findIsotopes(xset2, ppm=20, mzabs=0.015,intval="intb")
xset4 <- groupCorr(xset3,cor_eic_th=0.1, pval=1.0, graphMethod="lpc", calcIso = TRUE, calcCiS = TRUE, calcCaS = ifelse(
 sum(sampclass(xset) == "sample") > 3, TRUE, FALSE)) #p. 28
xsetFA <- findAdducts(xset4, polarity="positive", rules = NULL) 
xset5<- getPeaklist(xsetFA)
xset5$rt_in_min<- (xset5$rt)/60	
xset5$mz_round <- round((xset5$mz),4)
xset5$rt_round <- round((xset5$rt_in_min),4)
xset5$mz_rt<- paste(xset5$mz_round,xset5$rt_round, sep="_")

# Delete rows found in blank and rows where ratio of blank to sample is greater than 0.33
xset5[is.na(xset5)] <- 0
if(length(which(sampclass(xset) == "Blank"))>1) {
  xset5$Blank_Average <- rowMeans(xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Blank"))],na.rm=T)
} else {
  xset5$Blank_Average <- xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Blank"))]
}
xset5$TIC_Average <- rowMeans(xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Sample")), drop=FALSE],na.rm=T)
xset5$TIC_ratio <- xset5$Blank_Average/xset5$TIC_Average
xset6 <- xset5[xset5$Blank==0 & xset5$TIC_ratio < 0.33, ]

write.csv(xset6, "aduncum_v5_pos.csv")
