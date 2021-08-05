###### Part 2 of functions to Chemical Structural Similarity (CSS), adapted from Sedio et al.Applications in Plant Sciences 2018 

## Requires the following from Part 1:
## sampsByCompounds: rows are samples, columns are molecular features, entries are ion intensity (we use 1/0 to indicate presence/absence)
## row names can be sample or species names, column names must be feature ID numbers ('CLUSTERID') from the GNPS output (to match the first two columns of 'network')
## network: "network.csv" or your curated version. 

## Set working directory to location of Cytoscape network files 
setwd("../ProteoSAFe-METABOLOMICS-SNETS-b7ed685e-view_network")

## if network and/or compound list were curated or otherwise adjusted after MolecNetsChemTraits, load adjusted files
network <- read.csv("network_curated.csv", header = TRUE)
network$CLUSTERID1 = as.factor(network$CLUSTERID1)
network$CLUSTERID2 = as.factor(network$CLUSTERID2)
sampsByCompounds <- read.csv("sampsByCompounds_apf_ticuni.csv", header = FALSE)
sampnames <- read.csv("sampsByCompounds_all_pos_full_rownames.csv", header = TRUE)
clustnames <- read.csv("sampsByCompounds_all_pos_full_colnames.csv", header = TRUE)
row.names(sampsByCompounds) = sampnames$sample
colnames(sampsByCompounds) = clustnames$colnames

## date that MolecNetsTraits was run
date = 20181204

## if calculating CSCS for species pairs rather than sample pairs, set to TRUE
species = FALSE

outfile = paste("CSCS", date, ".RData",sep="")


calcCSCS = function(date = 20181204, species = FALSE, outfile)
  
	load(paste("MolecNetsChemTraits_piper12sp20181204.Rdata",sep="")) #load workspace from Part 1
	rm(list = c("defs", "ints", "speccodes", "freshmass"))
	if(species == TRUE){
		sampsByCompoundsOrig = sampsByCompounds
		sampsByCompounds = sppByCompounds
	}
	nspp = nrow(sampsByCompounds)
	ncomps = ncol(sampsByCompounds)
	net.comps = c(levels(network$CLUSTERID1), levels(network$CLUSTERID2))
	nspp = nrow(sampsByCompounds)
	ncomps = ncol(sampsByCompounds)
	pairwise.spp = as.data.frame(matrix(0,nrow = nspp, ncol = nspp))
	names(pairwise.spp) = row.names(sampsByCompounds)
	row.names(pairwise.spp) = row.names(sampsByCompounds)
	sampsCompsStand = sampsByCompounds
	for(i in 1:nrow(sampsByCompounds)){	
		sampsCompsStand[i,] = sampsByCompounds[i,]/sum(sampsByCompounds[i,])
	}
	diags = pairwise.spp
	for (k in 1:nspp){
		sppX = as.character(row.names(sampsCompsStand)[k])
		cat("Comparing ", sppX, " to itself", "\n", sep = "")
		sppXonly = sampsCompsStand[k,which(sampsCompsStand[k,]>0)]
		ncomps = length(sppXonly)
		pairwise.comps = as.data.frame(matrix(0, ncol = ncomps, nrow = ncomps))
		names(pairwise.comps) = names(sppXonly)
		row.names(pairwise.comps) = names(sppXonly)
		for (y in 1:ncomps){
			comp1 = names(sppXonly)[y]
			links.comp1 = network[which((network$CLUSTERID1 == comp1)|(network$CLUSTERID2 == comp1)),]
			for(z in 1:ncomps){
				comp2 = names(sppXonly)[z]
				if(comp2 %in% c(as.character(links.comp1$CLUSTERID1),as.character(links.comp1$CLUSTERID2))){
					if(comp1 == comp2){
						pairwise.comps[y,z] = 1
					}
					if(comp1 != comp2){				
						if((comp1 %in% links.comp1$CLUSTERID1)&(comp2 %in% links.comp1$CLUSTERID2)){
							pairwise.comps[y,z] = network$Cosine[which((network$CLUSTERID1 == comp1) & (network$CLUSTERID2 == comp2))]
						}
						if((comp1 %in% links.comp1$CLUSTERID2)&(comp2 %in% links.comp1$CLUSTERID1)){
							pairwise.comps[y,z] = network$Cosine[which((network$CLUSTERID2 == comp1) & (network$CLUSTERID1 == comp2))]
						}
					}				
				}
			}
		}
		diags[k,k] = sum(((outer(as.numeric(sppXonly), as.numeric(sppXonly)))*pairwise.comps),na.rm = T)	
	}
	save(sampsByCompounds, pairwise.comps, diags, file = outfile)
	
	write.csv(diags, "Data_Intrasample_Structural_Similarity_matrix.csv") # within-sample similarity scores for structural complexity analysis
	
	for (i in 1:nspp){
		spp1 = as.character(row.names(sampsCompsStand)[i])
		for (j in i:nspp){
			spp2 = as.character(row.names(sampsCompsStand)[j])
			cat("Comparing ", spp1, " to ", spp2, "\n", sep = "")
			#identify which compounds are in each species
			spp1comps = sampsCompsStand[spp1,]
			spp2comps = sampsCompsStand[spp2,]
			spp_pair = rbind(spp1comps,spp2comps)
			paircomps = spp_pair[,which(colSums(spp_pair)>0)]
			#make a pairwise.comps matrix for only those compounds found in either species in the species pair
			ncomps = ncol(paircomps)
			pairwise.comps = as.data.frame(matrix(0, ncol = ncomps, nrow = ncomps))
			names(pairwise.comps) = names(paircomps)
			row.names(pairwise.comps) = names(paircomps)
			cat("Beginning pairwise comparison of compounds","\n", sep = "")
			for (y in 1:ncomps){
				comp1 = names(paircomps)[y]
				links.comp1 = network[which((network$CLUSTERID1 == comp1)|(network$CLUSTERID2 == comp1)),]
				for(z in 1:ncomps){
					comp2 = names(paircomps)[z]
					if(comp2 %in% c(as.character(links.comp1$CLUSTERID1),as.character(links.comp1$CLUSTERID2))){
						if(comp1 == comp2){
							pairwise.comps[y,z] = 1
						}
						if(comp1 != comp2){				
							if((comp1 %in% links.comp1$CLUSTERID1)&(comp2 %in% links.comp1$CLUSTERID2)){
								pairwise.comps[y,z] = network$Cosine[which((network$CLUSTERID1 == comp1) & (network$CLUSTERID2 == comp2))]
							}
							if((comp1 %in% links.comp1$CLUSTERID2)&(comp2 %in% links.comp1$CLUSTERID1)){
								pairwise.comps[y,z] = network$Cosine[which((network$CLUSTERID2 == comp1) & (network$CLUSTERID1 == comp2))]
							}
						}				
					}
				}
			}
			pairwise.spp[i,j] = pairwise.spp[j,i] = sum(((outer(as.numeric(paircomps[1,]), as.numeric(paircomps[2,])))*pairwise.comps), na.rm = T)/max(diags[i,i], diags[j,j])
		}
	}
	cscs = pairwise.spp
	sampsByCompounds = sampsByCompoundsOrig
	save(sampsByCompounds, sppByCompounds, sampsCompsStand, diags, cscs, file = outfile)
	cat("Completed cacluation of CSCS for all sample pairs","\n")
}
write.csv(cscs, "Data_Intersample_Structural_similarity.csv")
   # this is the sample x sample matrix

################# for converting triangular matrix to columnar format - used for Data_Intersample_Structural_similarity and Data_Intrasample_Structural_similarity
#move row and column headers to separate files
cscs <- read.csv("Data_Intersample_Structural_similarity_noheaders.csv", header = F)
names <- read.csv("css_matrix_header_names.csv", header = TRUE)
row.names(cscs) = names$sample1
colnames(cscs) = names$sample2

dist <- cscs[lower.tri(cscs)] 
who.vs.who <- expand.grid(rownames(cscs), rownames(cscs)) 
who <- who.vs.who[lower.tri(cscs),] 
names(dist) <- paste(who[,1], who[,2], sep=".vs.") 
dist 

write.csv(dist, "columnar.csv")
