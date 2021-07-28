

# Piper-metabolomics

These data and code are to accompany the manuscript: “ Comparative metabolomics of fruits and leaves in a hyperdiverse lineage suggests fruits are a key incubator of phytochemical diversification”
by Gerald F. Schneider, Diego Salazar Amoretti, Sherry B. Hildreth, Richard F. Helm, Susan R. Whitehead

Publication information: Accepted for publication in Frontiers in Plant Science, *official date TBD*

DOI: *TBD*

Our overall objective in this study is to test the hypothesis that fruits can act as incubators of phytochemical diversification in plants. First, we describe the occurrence patterns of secondary metabolites across leaves, fruit pulp, and seeds in 12 *Piper* species, providing baseline data for understanding *Piper* secondary metabolite function. We use untargeted mass spectrometry-based metabolomics, molecular networking, and in-silico fragmentation modeling to characterize undescribed metabolites at the class level, followed by machine learning and distance-based methods to compare composition across organs and species. Second, we use these data to test predictions of high relative diversity in fruits derived from our hypothesis of fruit-driven phytochemical diversification. We compare multiple dimensions of phytochemical diversity across leaves and fruits, including the richness at multiple scales (alpha and gamma diversity), variability (beta diversity), and structural complexity of secondary metabolites.

## Data files

This repository contains 3 data files:

### 1)  Data_Intersample_Structural_Similarity.csv
This is a raw output file from the chemical structural similarity analysis. It is a matrix displaying the mean cosine-scored pairwise chemical structural similarity between each sample. Each row and each column represent one sample, as indicated by the sample name in the row and column headings. Each cell shows the similarity of the sample indicated in the row heading to the sample indicated in the column heading.
   
### 2) Data_Intersample_Structural_Similarity.csv
This table is derived from the chemical structural similarity analysis. It provides the mean cosine-scored pairwise chemical structural similarity between each compound in a given sample, shown for all samples and with the species, tissue type/organ, and individual plant ID listed alongside each sample ID. 

Variables include:

SampleID: the identification code for each individual sample as run on the UPLC-MS instrument. The format is "species_organ_sampleNumber_dateRun".

Species: the epithet of the *Piper* species to which the sample belongs.

Tissue: the tissue or organ from which the sample was extracted. Leaf = mature leaf, ripe = ripe pulp, unripe = unripe pulp, seed = ripe seed. 

PlantID: the identification number of the individual plant from which the sample was collected. 

chem_similarity_internal: the mean cosine-scored pairwise chemical structural similarity between each compound in a given sample.

### 3) Data_Peak_Table.csv
This is a table containing the curated output of all molecular features and ion abundances yielded by XCMS-CAMERA processing of raw data from low ionization energy UPLC-MS experiments. All putative molecular ions and their TIC abundances are included for each sample. Features annotated by CAMERA as adducts, isotopologues, or contaminants, have been removed.  

Variables include:

SampleID: the identification code for each individual sample as run on the UPLC-MS instrument. The format is "species_organ_sampleNumber_dateRun".

sp: the species epithet of the *Piper* species to which the sample belongs

tissue: the tissue or organ from which the sample was extracted. Leaf = mature leaf, ripe = ripe pulp, unripe = unripe pulp, seed = ripe seed. 

PlantID: the identification number of the individual plant from which the sample was collected. 

All other variables are molecular ions, represented by their mass:charge ratio and retention time in minutes, written as "m/z_retention time"



## Analysis Scripts

This repository contains 4 R scripts.

### 1) Script_DiversityAnalyses_master.R,
This script was used for Bray-Curtis, NMDS, machine learning, and all diversity analyses.

### 2) molec_net_piper12sp_pt1_v3.R
This script was used to prepare fragmentation spectra for chemical structural similarity analyses.

### 3) molec_net_piper12sp_pt2_v4.R
This script was used to conduct chemical structural similarity analyses.

### 4) xcms_params_pos_CDF_2019_07_17.R
This script was used to apply XCMS and CAMERA to extract and annotate molecular ion peak tables from raw UPLC-MS files.
