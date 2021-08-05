###### Part 1 of functions to calculate Chemical Structural Similarity (CSS), adapted from Sedio et al.Applications in Plant Sciences 2018 


## Set working directory to location of Cytoscape network files (downloaded from GNPS)
setwd(".../ProteoSAFe-METABOLOMICS-SNETS-b7ed685e-view_network")

## the 8-character unique ID for the GNPS network, e.g. "d1f7f083"
code = "b7ed685e"

## the date the script is run, to append to the output data file
date = 20181204

## filename under which to save the results
outfile = paste("MolecNetsChemTraits_piper12sp", date, ".RData",sep="")


molecNetsTraits = function(code, date, clone = TRUE, outfile = paste("MolecNetsChemTraits", date, ".RData",sep="")){
	wkdir = getwd()
	network = read.table(".../ProteoSAFe-METABOLOMICS-SNETS-b7ed685e-view_network/networkedges_selfloop/3e498e23668742138671bd51c08d2971.pairsinfo", header = TRUE)
	network$CLUSTERID1 = as.factor(network$CLUSTERID1)
	network$CLUSTERID2 = as.factor(network$CLUSTERID2)
	netshort = read.table(".../ProteoSAFe-METABOLOMICS-SNETS-b7ed685e-view_network/networkedges/4a515a24d4db41d29edadb6dd63bdc1c.pairsinfo", header = TRUE)
	names(netshort) = names(network)
	netshort$CLUSTERID1 = as.factor(netshort$CLUSTERID1)
	netshort$CLUSTERID2 = as.factor(netshort$CLUSTERID2)
		
	### enter sample gravimetric variables (extract mass or leaf mass)
	sampmass = read.csv("FreshMass.csv", header = T)
	sampmass$sampmass = as.numeric(as.character(sampmass$sampmass))
	atts = read.table(".../ProteoSAFe-METABOLOMICS-SNETS-b7ed685e-view_network/clusterinfosummarygroup_attributes/24ae481144eb43f5a628db324ba1ce82..out", header = TRUE, comment.char = "")
	atts$cluster.index = as.factor(atts$cluster.index)
	blanks = sort(union(grep("blank", atts$AllFiles), grep("blank", atts$AllFiles))) #remove blanks
	defs = atts[-blanks,]
	defs$cluster.index = as.factor(defs$cluster.index)
	defs = droplevels(defs)
		
	write.table(defs, file = paste("AttributesBlanksRemoved", date, ".txt",sep=""), quote = F, row.names = F, sep = "	")
	netshortdefs = netshort[which(netshort$CLUSTERID1 %in% defs$cluster.index),]
	netshortdefs = netshortdefs[which(netshortdefs$CLUSTERID2 %in% defs$cluster.index),]
	netdefs = network[which(network$CLUSTERID1 %in% defs$cluster.index),]
	netdefs = netdefs[which(netdefs$CLUSTERID2 %in% defs$cluster.index),]
	write.table(netshortdefs, file = paste("NetworkBlanksRemovedNoSingletons", date, ".txt",sep=""), quote = F, row.names = F, sep = "	")
	write.table(netdefs, file = paste("NetworkBlanksRemoved", date, ".txt",sep=""), quote = F, row.names = F, sep = "	")
	
	speccodes = read.table(paste(wkdir,"/SampleSpecMap.txt",sep=""), header = T, sep = "	")
	ints = read.table(file = paste(wkdir,"/clusterinfogroup/",list.files(paste(wkdir,"/clusterinfogroup/",sep = "")),sep = ""), header = T, sep = "	", comment.char = "")
	ints$X.ClusterIdx = as.factor(ints$X.ClusterIdx)
	ints$Conc = NA
	totalnum = nrow(speccodes)
	cat("Calculating ion intensity of compounds in each sample", "\n")
	for(i in 1:totalnum){
		records = grep(as.character(speccodes$SpecCode[i]),ints$X.Filename)
		sampname = as.character(speccodes$Sample[i])
		if(sampname != "blank"){
			if(sampname %in% sampmass$Sample){
				sampmass = sampmass$sampmass[which(sampmass$Sample == sampname)]
				ints$Conc[records] = (ints$X.PrecIntensity[records])/sampmass
			}
		}
		cat("Finished sample", i, "of ", totalnum, "\n")
	}
	leafsamps = as.character(levels(as.factor(speccodes$Sample)))
	sampsByCompounds = as.data.frame(matrix(0,nrow = nrow(speccodes), ncol = length(levels(defs$cluster.index))))
	names(sampsByCompounds) = as.character(levels(defs$cluster.index))
	for (i in 1:nrow(speccodes)){
		cat("Analyzing", as.character(speccodes$Species[i]), "\n")
		specsample = ints[grep(speccodes$SpecCode[i], ints$X.Filename),] # selects a sample, then selects the MS that came from it
		origfile = as.character(speccodes$OrigFilename[i])
		# for each individual MS, with a unique SpecIdx, measure the amount per leaf mass (i.e. concentration or Conc)
		for (j in 1:nrow(specsample)){
			specIndex = specsample$X.SpecIdx[j]
			if(length(grep(pattern = paste(origfile,":",specIndex, sep = ""), defs$AllFiles)) > 0){
				defrows = grep(pattern = paste(origfile,":",specIndex, sep = ""), defs$AllFiles)
				clusterIndexes = as.character(defs$cluster.index[defrows])
				sampsByCompounds[i,clusterIndexes] = specsample$Conc[j]
			}
		}
		cat("Completed sample #", i, "out of", nrow(speccodes), "\n")
	}
	row.names(sampsByCompounds) = speccodes$OrigFilename
	sampsByCompounds = sampsByCompounds[,which(colSums(sampsByCompounds)>0)]
	save(sampsByCompounds, ints, defs, network, netshort, speccodes, sampmass, file = outfile)
	write.csv(sampsByCompounds, "sampsByCompounds.csv")
	spp = levels(speccodes$Species)
	sppByCompounds = as.data.frame(matrix(0,nrow = length(spp),ncol = ncol(sampsByCompounds)))
	row.names(sppByCompounds) = spp
	names(sppByCompounds) = names(sampsByCompounds)
	for(i in 1:length(spp)){
		if(length(grep(as.character(spp[i]), speccodes$OrigFilename))==1){
			sppByCompounds[i,] = sampsByCompounds[grep(as.character(spp[i]), speccodes$OrigFilename),]
		}
		if(length(grep(as.character(spp[i]), speccodes$OrigFilename))>1){
			sampsSpp = sampsByCompounds[grep(as.character(spp[i]), speccodes$OrigFilename),]
			sppByCompounds[i,] = colMeans(sampsSpp)
		}
	}
	save(sppByCompounds, sampsByCompounds, ints, defs, network, netshort, speccodes, sampmass, file = outfile)
	cat("Analysis Complete","\n")
}

write.csv(network, "network.csv") #network data in columnar format, containing all pairwise similarity scores for all molecular features.
	
