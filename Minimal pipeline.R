# Note this process is extremely resource intensive (especially the total amount of RAM required).
# Run on a workstation with dual Xeon X5690's (12C/24T@3.5GHz) with 80gb of DDR3 (1333MHz) memory.
# You've been warned. 

##### Back End ####

rm(list=ls()) # Clear the workspace if required

# Install 
install_packages <- function(){
  cran_pkg <- c("RColorBrewer", "devtools", "fitdistrplus", "plotly", "reshape", "lars", "magrittr", "factoextra", 
                "Rserve", "RColorBrewer", "xtable", "som", "ROCR", "RJSONIO", "gplots", "e1071", 
                "caTools", "igraph", "randomForest", "Cairo", "pls", "pheatmap", "lattice", "rmarkdown",
                "knitr", "data.table", "pROC", "Rcpp", "caret", "ellipse", "scatterplot3d")
  bioconductor_pkg <- c("xcms", "IPO", "CAMERA", "impute", "pcaMethods", "siggenes", "globaltest",
                        "GlobalAncova", "Rgraphviz", "KEGGgraph", "preprocessCore", "genefilter",
                        "SSPA", "sva")
  list_installed <- installed.packages()
  new_cran <- subset(cran_pkg, !(cran_pkg %in% list_installed[, "Package"]))
  if(length(new_cran)!=0){
  install.packages(new_cran, dependencies = TRUE)
  print(c(new_cran, " packages added..."))
  }
  new_bio <- subset(bioconductor_pkg, !(bioconductor_pkg %in% list_installed[, "Package"]))
  if(length(new_bio)!=0){
    source("https://bioconductor.org/biocLite.R")
    biocLite(new_bio, dependencies = TRUE, ask = FALSE)
    print(c(new_bio, " packages added..."))
  }
  print("No new packages added...")
}
install_packages()
install.packages("backports")
library("devtools")
devtools::install_github("xia-lab/MetaboAnalystR", build_vignettes=TRUE)
install.packages("MALDIquantForeign")
# Libraries (All will call xcms)
library(IPO) 
library(CAMERA)
library(MetaboAnalystR)
library(RColorBrewer)
library(MALDIquant)
library(MALDIquantForeign)

##### 0.1 Preparation (Can run all - Update? What's needed?) #####
# Prior to analysis, first define the variables of the experimental design in the following 
# section. These paramters will be called at various points of the analysis, so you don't 
# have to. This analysis is currently tailored to a four group analysis, but will be modified
# to enable a two group analysis in future. 

## Debug function - Give index of all non-numeric data in the df. 
## Will show up as integer(0) if there is none. 
which.nonnum <- function(x) {
  badNum <- is.na(suppressWarnings(as.numeric(as.character(x))))
  which(badNum & !is.na(x))
}
which.nonnum(diffreport[,17:88])

  # Using the mzXML's due to increased  information regarding the metadata. May change back to 
  # CDF for ease (no special conversion software required) if no difference is made. 
setwd("D:/Metabolomics")

# [A] - DoE parameters:

  # Total number of technical replicates:
  Treps <- 72
  
  # Total number of groups analysed: (placeholder, just 4 for now)
  Tgroups = 4

  # Name the conditions for the experiment. 
  # (in the alphabetical order that they appear in the working directory)
  Condition1 = "2>2"
  Condition2 = "2>8"
  Condition3 = "8>2"
  Condition4 = "8>8"
  
  Treatment1 = "Low Protein"
  Treatment2 = "High Protein"
  
# [B] - Subsequent grouping information:
  
  # Create grouping for technical replicates (TR)
  TR_groups <- data.frame(matrix(NA, nrow=1, ncol=Treps))
  TR_groups[1:(Treps*0.25)]<-Condition1
  TR_groups[((Treps*0.25)+1):(Treps*0.50)]<-Condition2
  TR_groups[((Treps*0.50)+1):(Treps*0.75)]<-Condition3
  TR_groups[((Treps*0.75)+1):(Treps)]<-Condition4
  TR_groups <- t(TR_groups)
  
  # Create grouping for averaged results (AV) - Provided that there is three technical replicates for each sample.
  AV_groups <- data.frame(matrix(NA, nrow=1, ncol=(Treps/3)))
  AV_groups[1:((Treps/3)*0.25)]<-Condition1
  AV_groups[(((Treps/3)*0.25)+1):((Treps/3)*0.50)]<-Condition2
  AV_groups[(((Treps/3)*0.50)+1):((Treps/3)*0.75)]<-Condition3
  AV_groups[(((Treps/3)*0.75)+1):((Treps/3))]<-Condition4
  AV_groups <- t(AV_groups)
  
  # Generate a numeric list of the technical reps in the working directory #
  # This needs to be changed each time, or input a number list yourself (currently set to mzXML). 
  FileNames <- list.files(pattern = ".mzXML", recursive = TRUE)
  FileNames <- as.data.frame(FileNames)
  rownames(FileNames) <- FileNames[,1]
  FileNames <- FileNames[,-1]
  
  library(tidyr)
  
  FileNames <- rownames(FileNames) %>%
    strsplit( "mzXML/2-2/ac-090218-0" ) 
  FileNames <- as.data.frame(FileNames)
  FileNames <- FileNames[-1,]
  colnames(FileNames) <- FileNames[2,]
  FileNames <- t(FileNames)
  FileNames <- as.data.frame(FileNames)
  rownames(FileNames) <- FileNames[,1]
  FileNames <- FileNames[,-1]
  
  FileNames <- rownames(FileNames)%>%
    strsplit( "mzXML/2-8/ac-090218-0" )
  FileNames <- as.data.frame(FileNames)
  FileNames <- FileNames[-1,]
  colnames(FileNames) <- FileNames[2,]
  FileNames <- t(FileNames)
  FileNames <- as.data.frame(FileNames)
  rownames(FileNames) <- FileNames[,1]
  FileNames <- FileNames[,-1]
  
  FileNames <- rownames(FileNames)%>%
    strsplit( "mzXML/8-2/ac-090218-0" )
  FileNames <- as.data.frame(FileNames)
  FileNames <- FileNames[-1,]
  colnames(FileNames) <- FileNames[2,]
  FileNames <- t(FileNames)
  FileNames <- as.data.frame(FileNames)
  rownames(FileNames) <- FileNames[,1]
  FileNames <- FileNames[,-1]
  
  FileNames <- rownames(FileNames)%>%
    strsplit( "mzXML/8-8/ac-090218-0" )
  FileNames <- as.data.frame(FileNames)
  FileNames <- FileNames[-1,]
  colnames(FileNames) <- FileNames[2,]
  FileNames <- t(FileNames)
  FileNames <- as.data.frame(FileNames)
  rownames(FileNames) <- FileNames[,1]
  FileNames <- FileNames[,-1]
  
  rownames(FileNames) <- substr(rownames(FileNames),1,nchar(rownames(FileNames))-6) # delete the last 4 characters (.CDF)
  
  AVFileNames = FileNames[seq(3, nrow(FileNames), 3), ]
  AVFileNames[,1] <- rownames(AVFileNames)
  AVFileNames <- as.numeric(AVFileNames[,1])/3
  AVFileNames <- as.data.frame(AVFileNames)

##### 1.1 Data Import: #####
# Use following for two factor analysis example: https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html
setwd("D:/Metabolomics")
  ## IMPORTANT ##
  # To provide correct grouping information for analysis, be sure to place the raw data files 
  # into their own subfolders according to group! (e.g. CDF>'2-2')
  # Use either 'CDFs' or 'mzXMLs' depending on YOUR data format. 
  # Note ProteoWizard can convert most outputs to usable formats if conversion software isn't provided.
# CDFs <- list.files(path="D:/Metabolomics", pattern=".CDF", full.names = TRUE, recursive = TRUE)
mzXMLs <- list.files(path="D:/Metabolomics", pattern=".mzXML", full.names = TRUE, recursive = TRUE)

##### 1.2 Establish approximate baseline parameters #####
# Before using IPO to optimise parameters automatically (step 1.3) it is best to define as many
# parameters before using IPO to optimise the rest. Some parameters can be closely approximated 
# (such as ppm accuracy) from the methodology or analyser used - ask the technician for details!
# We'll be using MALDIquant to establish some of these parameters. Note these may not correspond
# exactly to the deconvolution process performed by XCMS. (both 8513!! :] )

# [A] - Identify mass to charge (m.z.) and intensity ranges of the data
spectra <- import(mzXMLs) # Process data with MALDIquant (~45 mins)
length(spectra) # Number of MassSpectrum objects generated from the data. If the same as the number 
                # of elements in the 'raw_data' from earlier then that's perfect!
spectra[1:2] # Can check the spectra for their properties. (M.Z. Range)
intensity.range <- range(list(spectra[[1]]@intensity)) # Intensity range (for later use)
mz.range <- range(list(spectra[[1]]@mass)) # Extract the m.z range for spectra 1 and save for later. 
                                          # Same for all spectra, but getting the range from all nested
                                          # list obejcts of the S4 object is just wasting too much time.


# [B] - Determine spectra quality
any(sapply(spectra, isEmpty)) # Test whether all spectra contain the same number of data points 
                              # and are not empty. FALSE means they aren't empty, good.(Just a quick QC step)

# [C] - Set noise threshold (snthresh) estimation
# To eliminate the background 'noise' from actual peaks, it is best to visualise a number of example
# spectra with the thresholds plotted to ascertain background noise. Being strict would retain only 
# the strongest peaks (fewer results) but with greater confidence, but being too relaxed (e.g. snthr =1)
# could introduce noise into the final dataset. A balance between both is desired for a metabolic
# overview, but be aware of the potential bias that could be introduced into the data if the bounds are
# set too close to either extreme.

samples <- factor(sapply(spectra, function(x)metaData(x)$sampleName))
spectra <- transformIntensity(spectra,method="sqrt")
avgSpectra <- averageMassSpectra(spectra, labels=samples, method="mean")
noise <- estimateNoise(avgSpectra[[1]])

# Use the function below to quickly conduct multiple plots
snthresh.plot <- function(a,b,x,y) {
  plot(avgSpectra[[1]], xlim=c(a, b), ylim=c(x, y))
  lines(noise, col="red")
  lines(noise[,1], noise[, 2]*2, col="blue")
  lines(noise[,1], noise[, 2]*3, col="green")  
  lines(noise[,1], noise[, 2]*4, col="purple")
  lines(noise[,1], noise[, 2]*5, col="cyan")
  lines(noise[,1], noise[, 2]*6, col="orange")
  lines(noise[,1], noise[, 2]*7, col="grey")
  lines(noise[,1], noise[, 2]*8, col="maroon")
  lines(noise[,1], noise[, 2]*9, col="yellow")
  lines(noise[,1], noise[, 2]*10, col="cornflowerblue")
  lines(noise[,1], noise[, 2]*11, col="salmon")
  lines(noise[,1], noise[, 2]*12, col="brown")
  lines(noise[,1], noise[, 2]*13, col="pink")
  lines(noise[,1], noise[, 2]*14, col="darkblue")
  lines(noise[,1], noise[, 2]*15, col="magenta")
}

# Plots - 4 values, (a,b) = m.z range, (x,y) = Intensity range.
        # Refer to the above colours for snthr value deemed most suitable.
plot(avgSpectra[[1]]) # Whole overview (the intensity ranges correspond to the non-transformed data...)
#snthresh.plot(mz.range[1],mz.range[2],intensity.range[1],intensity.range[2]) # Whole overview FIX!!
snthresh.plot(mz.range[1],mz.range[1]+50,intensity.range[1],intensity.range[1]+50) # Front-end zoom
snthresh.plot(mz.range[2]-50,mz.range[2],intensity.range[1],intensity.range[1]+50) # End-end zoom
snthresh.plot(mz.range[2]/2,(mz.range[2]/2)+50,intensity.range[1],intensity.range[1]+50) # Mid-end zoom

# We're looking to set a level which cut's out the noise from the data, so for the End-end zoom plot, 
# This is fairly obvious, choose the snthr line which sits ABOVE the wavelets. For the other plots, use 
# your best judgement to determine the correct snthr boundary to eliminate the noise. I suggest
# (especially for close calls) to save a two value range.
#
# In this case, snthr =3 seems to be appropriate for the majority of data, whilst snthr = 5 is safe,
# but may inadvertently cut out some real peaks.
#
# In any case, save two values for IPO to optimise later on, replace the two values below (low,high):

snthr.range <- c(3,5)

# [D] Retention Time
# The data was collected using MALDI (matrix-assisted laser desorption-ionisation) WITHOUT any prior
# separation with LC/HPLC/GC step. Common retention times of say 2000-3000 seconds DO NOT APPLY. Use 
# The following functions to estimate common times 
head(rtime(raw_data)) # rtime for some of the points
mean(rtime(raw_data)) # mean rtime
median(rtime(raw_data)) # median rtime
range(rtime(raw_data)) # min/max rtime
plot(rtime(raw_data)) # check rtime consistency (gaps / outliers)

rtime <- range(rtime(raw_data)) # Set this variable for later use

# [C] - Chromatogram Plots - FYI - Raw data overview.

  ## Set colour and grouping information
  group_colors <- paste0(brewer.pal(4, "Set1")[1:4], Treps)
  names(group_colors) <- c(Condition1, Condition2, Condition3, Condition4)
 
  ## Plot the total ion chromatogram (TIC).
  # A TIC (Total Ion Chromatogram) is a chromatogram created by summing up intensities of all mass 
  # spectral peaks belonging to the same scan.
  tic_raw <- chromatogram(raw_data, aggregationFun = "sum")
  plot(tic_raw, col = group_colors[raw_data$sample_group], 
       xlab = "Retention Time (seconds)", 
       ylab = "Intensity (CPS)")
 
  ## Plot the base peak chromatogram (BPC).
  # The base peak chromatogram is similar to the TIC chromatogram, however it monitors only the most
  # intense peak in each spectrum. This means that the base peak chromatogram represents the  intensity
  # of the most intense peak at every point in the analysis. Base peak chromatograms often have a
  # cleaner look and thus are more informative than TIC chromatograms because the background  is
  # reduced by focusing on a single analyte at every point.
  #
  # TL;DR: The BPC will often look similar to the TIC, but less noisy. However, chromatographic peaks
  # consisting of signal at many masses can be understated by the BPC.
  bpc_raw <- chromatogram(raw_data, aggregationFun = "max")
  plot(bpc_raw, col = group_colors[raw_data$sample_group],
       xlab = "Retention Time (seconds)", 
       ylab = "Intensity (CPS)")
    
  
  
## Need to look at the individual parameters for peakpickingParameters and see how many I can optimise.

##### 1.3 Optimise XCMS parameters with IPO #####

  # [A] - Define the IPO parameters for testing

  # If the data is centroided, use 'centWave' in the first arguement to load the default parameters.
  # If the data is profile, use 'matchedFilter' instead. - Need matchedFilter for our use.
  # Note not all of the parameters listed may be equivalent / used.
peakpickingParameters <- getDefaultXcmsSetStartingParams('matchedFilter')
peakpickingParameters
peakpickingParameters$profStep <- c(1) # step size (in m/z) to use for profile generation from the raw data files
peakpickingParameters$bw <- c(5) # allowable retention time deviation (seconds) ''
peakpickingParameters$minfrac <- c(0.5) # minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group.
peakpickingParameters$minsamp <- c(1) # minimum number of samples necessary in at least one of the sample groups for it to be a valid group.
peakpickingParameters$max <- c(10) # representing the maximum number of peaks that are expected/will
                                       # be identified per slice.
peakpickingParameters$step <- c(0.1, 0.2) # step size to use for profile generation
peakpickingParameters$fwhm <- c(40, 50) # specifying the full width at half maximum of matched filtration
                                        # gaussian model peak. Only used to calculate the actual sigma,
peakpickingParameters$ppm <- c(5,10) # ppm tolerance for identification
peakpickingParameters$snthresh <- snthr.range # defines the signal to noise ratio cutoff (from earlier).
peakpickingParameters$mzdiff <- c(0.01) # Default from centWave which allows overlap, '0' was previous value... let's see. 
peakpickingParameters ## Check the parameters if desired. 

# [B] - Run the defined IPO parameters to determine XCMS 

time.xcmsSet <- system.time({ # measuring time
  resultPeakpicking <- 
    optimizeXcmsSet(files = mzXMLs, 
                    params = peakpickingParameters, 
                    BPPARAM = SnowParam(), #linux; MulicoreParam() #windows; SnowParam()
                    subdir = "IPO",
                    plot = TRUE)
})


# [C] - Optimise peak picking result
resultPeakpicking$best_settings$result
optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset

# [D] - Optimise retention time correction and grouping parameters
retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- 1
retcorGroupParameters$gapExtend <- 2.7
time.RetGroup <- system.time({ # measuring time
  resultRetcorGroup <-
    optimizeRetGroup(xset = optimizedXcmsSetObject, 
                     params = retcorGroupParameters, 
                     nSlaves = 24, 
                     subdir = NULL,
                     plot = TRUE)
})


# [D] - Running times and session info (If interested)
time.xcmsSet # time for optimizing peak picking parameters, SECONDS (13 HOURS)
time.RetGroup # time for optimizing retention time correction and grouping parameters (2 HOURS)
sessionInfo() # Platform information

##### 1.4 Write/Run optimised R script #####
# -- Copy the response from the terminal in the section below! 
writeRScript(resultPeakpicking$best_settings$parameters, 
             resultRetcorGroup$best_settings)

# New R script:
library(xcms)
library(Rmpi)
library(CAMERA)

rm(xset)

xset <- xcmsSet( 
  method   = "matchedFilter",
  fwhm     = 25,
  snthresh = 1,
  step     = 0.04,
  steps    = 1,
  sigma    = 10.6166128758281,
  max      = 10,
  mzdiff   = 0.01,
  index    = FALSE)
xset <- retcor( 
  xset,         
  method         = "obiwarp",
  plottype       = "none",
  distFunc       = "cor_opt",
  profStep       = 1,
  center         = 69,
  response       = 1,
  gapInit        = 0.4,
  gapExtend      = 2.7,
  factorDiag     = 2,
  factorGap      = 1,
  localAlignment = 0)
xset <- group( ## 153307 mz slices...
  xset,
  method  = "density",
  bw      = 38,
  mzwid   = 0.015,
  minfrac = 0.7,
  minsamp = 1,
  max     = 50)

xset <- fillPeaks(xset, method="chrom")
xset # (info on dataset) # Should have the sample classes (grouping) according to the file structure. 

##### 1.5 Metabolome Report #####

# 22711 variables

# The 'diffreport' is an xcms generated report of the identified metabolomic features.
# There are many columns contained in the report, such as name, statistics and peak features
# Brief descriptions of columns are provided below, see table 1 for reference (DOI: 10.1021/ac300698c)

#   Column name	 -     Explanation
# name	         -  Feature name (arbitrary), formed by nominal mass and retention time, e.g., M120T7
# fold change	   -  Fold change (ratio of the mean intensities)
# p-value	       -  p-value (Welch t-test, unequal variances)
# m/z	           -  m/z value (median value for the aligned features)
#                -  (Min/Max also provided)
# retention time -	retention time (median value for the aligned features)
#                -  (Min/Max also provided)
# - My definitions
# x2.2 (e.g)     -  Number of peaks identifying this feature (prior to the fillpeaks/imputation function)

## Generate the XCMS diffreport of total identified features 
diffreport <- diffreport(xset)

## Adjust the column headers for downstream compatibility
names(diffreport)[names(diffreport) == 'mzmed'] <- 'm.z'
names(diffreport)[names(diffreport) == 'pvalue'] <- 'p.value'
names(diffreport)[names(diffreport) == 'tstat'] <- 't.score'

## Save the diffreport
write.table(diffreport, file="diffreport.txt", sep="\t", row.names = FALSE)

## Reload the generated diffreport if required
setwd("D:/Metabolomics")
diffreport <-read.table("diffreport.txt", sep="\t", header=TRUE)

##### 1.6 Average TR's into their respective groups #####

## Pull the transformed data back into the diffreport to keep track. 
diffreport[,17:(16+Treps)] <- data_eRD_tmm # Pull transformed data
colnames(diffreport)[17:(16+Treps)] <- row.names(FileNames) # Shorten names
View(colnames(diffreport)) # Check shortening was successful

## Save the diffreport (transformed/normalised)
write.table(diffreport, file="t_diffreport.txt", sep="\t", row.names = FALSE)

## Group the technical reps according to sample. 
AVdata <- t(diffreport[,17:(16+Treps)]) # Just using the data, no need to pull back from diffreprt
dim(AVdata)
AVdata <- rowsum(AVdata, rep(1:(Treps/3), each=3))/3 
dim(AVdata) # should = Treps / 3
rownames(AVdata) <- AVFileNames$AVFileNames
View(row.names(AVdata)) # Correct sample numbering

## Create & Save the diffreport (averaged & transformed/normalised)
AVdata <- t(AVdata)
diffreport[, c(17:(16+Treps))] <- list(NULL) # ()! - Remove old intensity data.
diffreport <- cbind(diffreport,AVdata)
write.table(diffreport, file="AVdiffreport.txt", sep="\t", row.names = FALSE) # Save AVeraged & Transformed diffreport

##### 1.7 Annotation and Pathway Enrichment Analysis <- not useful overall, but good for separate groups ##### 

library(MetaboAnalystR)

# Annotation
PeakListProfile <- diffreport[,c(6,4,3)]
dim(diffreport) # == 40 columns (AV intensities, =/= TRs).
colnames(diffreport[,c(6,4,3)]) # == 'm.z' 'p.value' 't.score' (change made in 2.0)

## Kegg
AV_mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
AV_mSet <- Read.PeakListData(AV_mSet, "AVdiffreport.txt");
AV_mSet <- UpdateMummichogParameters(AV_mSet, "6.0", "positive", 0.05); # PPM, ion mode, significance value:
  # Numeric, specify the p-value cutoff to define significant m/z features from reference m/z features
AV_mSet <- SanityCheckMummichogData(AV_mSet)
AV_mSet <- PerformMummichog(AV_mSet, "dme_kegg", "fisher", "gamma") # organism_database, test, distribution.
AV_mSet$matches.res$Common.Name <- doKEGG2NameMapping(AV_mSet$matches.res$Matched.Compound)
names(AV_mSet$matches.res)[names(AV_mSet$matches.res) == 'Query.Mass'] <- 'm.z'
names(AV_mSet$matches.res)[names(AV_mSet$matches.res) == 'Matched.Compound'] <- 'kegg.Matched.Compund'
names(AV_mSet$matches.res)[names(AV_mSet$matches.res) == 'Matched.Form'] <- 'kegg.Matched.Form'
names(AV_mSet$matches.res)[names(AV_mSet$matches.res) == 'Mass.Diff'] <- 'kegg.Mass.Diff'
names(AV_mSet$matches.res)[names(AV_mSet$matches.res) == 'Common.Name'] <- 'kegg.Common.Name'
write.table(AV_mSet$matches.res, file="AVKeggMatchedCompoundTable.txt", sep="\t", row.names = FALSE)
AV_MCT<-AV_mSet$matches.res
AV_MCT$m.z <- as.numeric(AV_MCT$m.z) # Fix output

# Natural merge the diffreport to the  MCT to provide annotated results
AnnotatedDiffReport <- merge.data.frame(diffreport, AV_MCT, by = "m.z")
write.table(AnnotatedDiffReport, file="ANAVdiffReport.txt", sep="\t", row.names = FALSE)

##### 1.8 Create a final df for MetabR, doi: 10.1186/1756-0500-5-596 #####
# The first three columns MUST be called 'ID', 'Subject' and 'Group'. Put whatever else you want in.

metabolites <- t(AnnotatedDiffReport[, c(17:40)])
colnames(metabolites) <- AnnotatedDiffReport$kegg.Matched.Compund

metainfo <- matrix(nrow = 24, ncol = 0 )
metainfo <- as.data.frame(metainfo)
metainfo$ID <- AVFileNames$AVFileNames
metainfo$Subject <- paste(AV_groups,"_Sample",AVFileNames$AVFileNames)
metainfo$Group <- c(rep("2>2", 6), rep("2>8", 6), rep("8>2", 6), rep("8>8", 6))
metainfo$Past.Treatment <- c(rep("Low.Protein", 12), rep("High.Protein", 12))
metainfo$Current.Treatment <- c(rep("Low.Protein", 6), rep("High.Protein", 6), rep("Low.Protein", 6), rep("High.Protein", 6))
metainfo$Type.Treatment <- c(rep("Control", 6), rep("Switch", 6), rep("Switch", 6), rep("Control", 6))

MetabRdata <- cbind(metainfo, metabolites)
write.csv(MetabRdata, "MetabRdata.csv", row.names = FALSE)


##### 1.9 MetabR Preparation #####

# Note, this 'Program' does tend to break, especially at startup. Warnings are fine (honestly, see the manual),
# but errors are well... not. Make sure you've actually selected a correct experimental design as otherwise this
# can break the 'program'
# Once loaded, follow the steps in the GUI. See manual for reference!!
# ,,,,,,,
# TO RUN # - After setup (below) just run ALL of the code in the next section. DO NOT do anything until the
# '''''''    console has finished loading, and be sure to check the manual if stuck. 
# http://metabr.r-forge.r-project.org/ #

install.packages(pkgs=c("lme4", "gplots", "R.oo", "lawstat", "gWidgetsRGtk2", "gWidgets"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue", version = "3.8")

install.packages("lme4")
library(lme4)
library(gWidgets)
library(qvalue)
library(Matrix)
setwd("D:/Metabolomics") # set the wd to wherever the data is stored!!
rm(list=ls()) # Clear the workspace if required
##### 2.0 MetabR - Run all, follow the GUI ####
# set variables
data.pos="NULL"
data.neg="NULL"

# set gtk layout
d = ggroup(container=gwindow("MetabR--Data input"), horizontal=FALSE)
file.type.label=glabel("File type", container=d)
file.type=gradio(c("csv", "txt"), selected=1, horizontal=TRUE, container=d) 

file.label=glabel("Data files",container=d)

data.pos.group=ggroup(container=d, horizontal=TRUE)
btn.pos <- gbutton(
  text 		= "Data 1", 
  container 	= data.pos.group,
  handler	= function (...) data.pos<<-gfile("Select a file",type="open")
)

data.neg.group=ggroup(container=d, horizontal=TRUE)
btn.neg <- gbutton(
  text 		= "Data 2", 
  container 	= data.neg.group,
  handler	= function (...) data.neg<<-gfile("Select a file",type="open")
)
btn.neg.clear=gbutton(
  text		= "Clear Data 2",
  container	= data.neg.group,
  handler	= function(...) data.neg<<-"NULL"
)


col.start.label <- glabel("Column of first data measurement", container = d)
col.start.group=ggroup(container=d, horizontal=TRUE)
col.start <- gedit("7", container = d, width=10, coerce.with=as.numeric)

ok = gbutton(text="OK", container=d,handler = function (...)
{ 
  col.start=svalue(col.start)	
  #options(show.error.messages=FALSE)
  #try.neg	= try(data.neg)
  #options(show.error.messages=TRUE)
  file.type	= svalue(file.type)
  if (file.type=="csv")
  {
    data.pos <<- read.csv(file=data.pos, header=TRUE, row.names=1, sep=",", check.names=FALSE)
    #if (class(try.neg)!="try-error")
    if (data.neg!="NULL")
    {
      data.neg <<- read.csv(file=data.neg, header=TRUE, row.names=1, sep=",", check.names=FALSE)
    }
  }else	{
    data.pos <<- read.table(file=data.pos, header=TRUE, row.names=1, check.names=FALSE)
    #if (class(try.neg)!="try-error")
    if (data.neg!="NULL")
    {
      data.neg <<- read.table(file=data.neg, header=TRUE, row.names=1, check.names=FALSE)
    }
  }
  #if (class(try.neg)=="try-error"){data.neg<<-NULL}
  #dispose(d)
  
  ###############################################################################################################################
  ###############################################################################################################################
  ##Begining of MetabR GUI and function
  
  
  w = gwindow("MetabR", width=1200, height=600)
  nb = gnotebook(tab.pos=2, container=w)
  
  fixed.continuous = 0
  random.factors = 0
  
  ##Required parameters##
  
  required.tab = ggroup(container=nb, label="Required", use.scrollwindow = TRUE, horizontal=TRUE)  
  
  if (as.numeric(col.start) >= 5)
  {
    required.1 = gframe(container=required.tab, horizontal=FALSE)
    fixed.continuous.label = glabel("Fixed-effect variables", container=required.1)
    fixed.continuous = gcheckboxgroup(colnames(data.pos)[3:(col.start-2)], container=required.1, 
                                      use.table=TRUE, coerce.with=svalue)
    
    random.factors.label = glabel("Random-effect variables", container = required.1)
    
    random.factors = gcheckboxgroup(colnames(data.pos)[3:(col.start-2)], container=required.1, 
                                    use.table=TRUE, coerce.with=svalue)
  }
  
  
  required.2 = gframe(container=required.tab, horizontal=FALSE)
  criterion.label = glabel("Criterion", container=required.2)
  criterion = gradio(c("Pval", "Qval", "Fold-change", "Pval+Fold-change", "Qval+Fold-change"), 
                     selected=1, horizontal=FALSE, container=required.2)
  plot.choice.label = glabel("Mean plot choice", container=required.2)
  plot.choice = gradio(c("Pval", "Qval"), selected=1, horizontal=FALSE, container=required.2)
  
  required.3 = gframe(container=required.tab, horizontal=FALSE)
  log.base.label <- glabel("Log base", container = required.3)
  log.base <- gedit("2", container = required.3, width=10, coerce.with=as.numeric)
  f.cutoff.label <- glabel("Fold-change cutoff", container = required.3)
  f.cutoff<- gedit("1.5", container = required.3, width=10, coerce.with=as.numeric)
  p.cutoff.label <- glabel("P or Q cutoff", container = required.3)
  p.cutoff<- gedit("0.05", container = required.3, width=10, coerce.with=as.numeric)
  output.label <- glabel("Output", container = required.3)
  output<- gedit("Exp1", container = required.3)
  
  
  #}) #stop here to just check correct parameter input through GUI
  
  ##Mean plot tab (optional)##
  
  meanplot.tab = ggroup(container=nb, label="Mean plot", use.scrollwindow = FALSE, horizontal=TRUE)   
  
  meanplot.1 = gframe(container=meanplot.tab, horizontal=FALSE)
  las.mean.label = glabel("las", container = meanplot.1)
  las.mean = gedit("3", container = meanplot.1, width=1, coerce.with=as.numeric)
  cex.mean.label = glabel("cex (par)", container = meanplot.1)
  cex.mean = gedit("1", container = meanplot.1, width=1, coerce.with=as.numeric)
  cex.axis.mean.label = glabel("cex.axis", container = meanplot.1)
  cex.axis.mean = gedit("1.2", container = meanplot.1, width=1, coerce.with=as.numeric)
  cex.lab.mean.label = glabel("cex.lab", container = meanplot.1)
  cex.lab.mean = gedit("1.2", container = meanplot.1, width=1, coerce.with=as.numeric)
  
  meanplot.2 = gframe(container=meanplot.tab, horizontal=FALSE)
  cex.main.mean.label = glabel("cex.main", container = meanplot.2)
  cex.main.mean = gedit("1.2", container = meanplot.2, coerce.with=as.numeric, width=1)
  cex.sub.mean.label = glabel("cex.sub", container = meanplot.2)
  cex.sub.mean = gedit("1.2", container = meanplot.2, coerce.with=as.numeric, width=1)
  mex.mean.label = glabel("mex", container = meanplot.2)
  mex.mean = gedit("1.2", container = meanplot.2, coerce.with=as.numeric, width=1)
  cex.mean.charsize.label = glabel("cex (text)", container = meanplot.2)
  cex.mean.charsize = gedit("1.2", container = meanplot.2, coerce.with=as.numeric, width=1)
  
  meanplot.3 = gframe(container=meanplot.tab, horizontal=FALSE)
  height.mean.label = glabel("height", container = meanplot.3)
  height.mean = gedit("7", container = meanplot.3, coerce.with=as.numeric, width=10)
  width.mean.label = glabel("width", container = meanplot.3)
  width.mean = gedit("7", container = meanplot.3, coerce.with=as.numeric, width=10)
  mfrow.mean.label = glabel("mfrow", container = meanplot.3)
  mfrow.mean = gedit("c(1,1)", container = meanplot.3, coerce.with=svalue, width=10)
  mar.mean.label = glabel("mar", container = meanplot.3)
  mar.mean = gedit("c(6,4,3,1)", container = meanplot.3, coerce.with=svalue, width=10)
  
  
  ##Residual plot##	
  
  residual.tab = ggroup(container=nb, label="Residual plot, unlabeled", horizontal=TRUE)   
  
  residual.1 = gframe(container=residual.tab, horizontal=FALSE)
  las.summaryplot.label = glabel("las", container = residual.1)
  las.summaryplot = gedit("3", container = residual.1, coerce.with=as.numeric, width=10)
  cex.summaryplot.label = glabel("cex", container = residual.1)
  cex.summaryplot = gedit("1", container = residual.1, coerce.with=as.numeric, width=10)
  cex.axis.summaryplot.label = glabel("cex.axis", container = residual.1)
  cex.axis.summaryplot<- gedit("1.2", container = residual.1, coerce.with=as.numeric, width=10)
  cex.lab.summaryplot.label <- glabel("cex.lab", container = residual.1)
  cex.lab.summaryplot<- gedit("1.2", container = residual.1, coerce.with=as.numeric, width=10)
  cex.main.summaryplot.label <- glabel("cex.main", container = residual.1)
  cex.main.summaryplot<- gedit("1.2", container = residual.1, coerce.with=as.numeric, width=10)
  
  residual.2=gframe(container=residual.tab, horizontal=FALSE)
  cex.sub.summaryplot.label <- glabel("cex.sub", container = residual.2)
  cex.sub.summaryplot<- gedit("1.2", container = residual.2, coerce.with=as.numeric, width=10)
  mex.summaryplot.label <- glabel("mex", container = residual.2)
  mex.summaryplot<- gedit("1", container = residual.2, coerce.with=as.numeric, width=10)
  height.summaryplot.label <- glabel("height", container = residual.2)
  height.summaryplot<- gedit("7", container = residual.2, coerce.with=as.numeric, width=10)
  width.summaryplot.label <- glabel("width", container = residual.2)
  width.summaryplot<- gedit("20", container = residual.2, coerce.with=as.numeric, width=10)
  mar.summaryplot.label <- glabel("mar", container = residual.2)
  mar.summaryplot<- gedit("c(6,4,5,1)", container = residual.2, coerce.with=svalue, width=10)
  
  
  ##Summary Plot with metabolite names##
  
  residual.labeled.tab = ggroup(container=nb, label="Residual plot, labeled", horizontal=TRUE)   
  
  residual.labeled.1=gframe(container=residual.labeled.tab, horizontal=FALSE)
  las.summaryplot.labeled.label <- glabel("las", container = residual.labeled.1)
  las.summaryplot.labeled<- gedit("3", container = residual.labeled.1, coerce.with=as.numeric, width=10)
  cex.summaryplot.labeled.label <- glabel("cex", container = residual.labeled.1)
  cex.summaryplot.labeled<- gedit("1", container = residual.labeled.1, coerce.with=as.numeric, width=10)
  cex.axis.summaryplot.labeled.label <- glabel("cex.axis", container = residual.labeled.1)
  cex.axis.summaryplot.labeled<- gedit("1.2", container = residual.labeled.1, coerce.with=as.numeric, width=10)
  cex.lab.summaryplot.labeled.label <- glabel("cex.lab", container = residual.labeled.1)
  cex.lab.summaryplot.labeled<- gedit("1.2", container = residual.labeled.1, coerce.with=as.numeric, width=10)
  cex.main.summaryplot.labeled.label <- glabel("cex.main", container = residual.labeled.1)
  cex.main.summaryplot.labeled<- gedit("1.2", container = residual.labeled.1, coerce.with=as.numeric, width=10)
  
  residual.labeled.2=gframe(container=residual.labeled.tab, horizontal=FALSE)
  cex.sub.summaryplot.labeled.label <- glabel("cex.sub", container = residual.labeled.2)
  cex.sub.summaryplot.labeled<- gedit("1", container = residual.labeled.2, coerce.with=as.numeric, width=10)
  mex.summaryplot.labeled.label <- glabel("mex", container = residual.labeled.2)
  mex.summaryplot.labeled<- gedit("1", container = residual.labeled.2, coerce.with=as.numeric, width=10)
  height.summaryplot.labeled.label <- glabel("height", container = residual.labeled.2)
  height.summaryplot.labeled<- gedit("15", container = residual.labeled.2, coerce.with=as.numeric, width=10)
  width.summaryplot.labeled.label <- glabel("width", container = residual.labeled.2)
  width.summaryplot.labeled<- gedit("200", container = residual.labeled.2, coerce.with=as.numeric, width=10)
  mar.summaryplot.labeled.label <- glabel("mar", container = residual.labeled.2)
  mar.summaryplot.labeled<- gedit("c(20,4,5,1)", container = residual.labeled.2, coerce.with=svalue, width=10)
  
  
  ##Pre- vs. post-normalization plots##
  
  prevpost.tab = ggroup(container=nb, label="Pre- vs. post-normalization plots", horizontal=TRUE)   
  
  prevpost.1=gframe(container=prevpost.tab, horizontal=FALSE)
  las.prevpost.label <- glabel("las", container = prevpost.1)
  las.prevpost <- gedit("3", container = prevpost.1, coerce.with=as.numeric, width=10)
  cex.prevpost.label <- glabel("cex", container = prevpost.1)
  cex.prevpost <- gedit("1", container = prevpost.1, coerce.with=as.numeric, width=10)
  cex.axis.prevpost.label <- glabel("cex.axis", container = prevpost.1)
  cex.axis.prevpost <- gedit("1.2", container = prevpost.1, coerce.with=as.numeric, width=10)
  cex.lab.prevpost.label <- glabel("cex.lab", container = prevpost.1)
  cex.lab.prevpost <- gedit("1.2", container = prevpost.1, coerce.with=as.numeric, width=10)
  cex.main.prevpost.label <- glabel("cex.main", container = prevpost.1)
  cex.main.prevpost <- gedit("1.2", container = prevpost.1, coerce.with=as.numeric, width=10)
  cex.sub.prevpost.label <- glabel("cex.sub", container = prevpost.1)
  cex.sub.prevpost <- gedit("1", container = prevpost.1, coerce.with=as.numeric, width=10)
  
  
  prevpost.2=gframe(container=prevpost.tab, horizontal=FALSE)
  mex.prevpost.label <- glabel("mex", container = prevpost.2)
  mex.prevpost <- gedit("1", container = prevpost.2, coerce.with=as.numeric, width=10)
  height.prevpost.label <- glabel("height", container = prevpost.2)
  height.prevpost <- gedit("3", container = prevpost.2, coerce.with=as.numeric, width=10)
  width.prevpost.label <- glabel("width", container = prevpost.2)
  width.prevpost <- gedit("6", container = prevpost.2, coerce.with=as.numeric, width=10)
  mfrow.prevpost.label <- glabel("mfrow", container = prevpost.2)
  mfrow.prevpost <- gedit("c(1,2)", container = prevpost.2, coerce.with=svalue, width=10)
  mar.prevpost.label <- glabel("mar", container = prevpost.2)
  mar.prevpost<- gedit("c(6,4,5,1)", container = prevpost.2, coerce.with=svalue, width=10)
  
  
  ##Levene test parameters##	
  levene.tab = ggroup(container=nb, label="Levene test", horizontal=TRUE)   
  
  levene.1=gframe(container=levene.tab, horizontal=FALSE)
  location.choice.label <- glabel("location", container = levene.1)
  location.choice <- gedit("median", container = levene.1, coerce.with=as.character)
  trim.alpha.choice.label <- glabel("trim.alpha", container = levene.1)
  trim.alpha.choice <- gedit("0", container = levene.1, coerce.with=as.numeric)
  bootstrap.choice.label <- glabel("bootstrap", container = levene.1)
  bootstrap.choice <- gedit("FALSE", container = levene.1, coerce.with=as.character)
  
  levene.2=gframe(container=levene.tab, horizontal=FALSE)
  num.bootstrap.choice.label <- glabel("num.bootstrap", container = levene.2)
  num.bootstrap.choice <- gedit("1000", container = levene.2, coerce.with=as.numeric)
  kruskal.test.choice.label <- glabel("kruskal.test", container = levene.2)
  kruskal.test.choice <- gedit("FALSE", container = levene.2, coerce.with=as.character)
  correction.method.choice.label <- glabel("correction.method", container = levene.2)
  correction.method.choice <- gedit("none", container = levene.2, coerce.with=as.character)
  
  ##Heatmap parameters##
  heatmap.tab = ggroup(container=nb, label="Heatmap", horizontal=TRUE)
  
  heatmap.1=gframe(container=heatmap.tab, horizontal=FALSE)
  las.heatmap.label <- glabel("las", container = heatmap.1)
  las.heatmap<- gedit("0", container = heatmap.1, coerce.with=as.numeric, width=10)
  cex.heatmap.label <- glabel("cex", container = heatmap.1)
  cex.heatmap<- gedit("1", container = heatmap.1, coerce.with=as.numeric, width=10)
  cex.axis.heatmap.label <- glabel("cex.axis", container = heatmap.1)
  cex.axis.heatmap<- gedit("1", container = heatmap.1, coerce.with=as.numeric, width=10)
  cex.lab.heatmap.label <- glabel("cex.lab", container = heatmap.1)
  cex.lab.heatmap<- gedit("1", container = heatmap.1, coerce.with=as.numeric, width=10)
  cex.main.heatmap.label <- glabel("cex.main", container = heatmap.1)
  cex.main.heatmap<- gedit("1", container = heatmap.1, coerce.with=as.numeric, width=10)
  
  heatmap.2=gframe(container=heatmap.tab, horizontal=FALSE)
  cex.sub.heatmap.label <- glabel("cex.sub", container = heatmap.2)
  cex.sub.heatmap<- gedit("1", container = heatmap.2, coerce.with=as.numeric, width=10)
  mex.heatmap.label <- glabel("mex", container = heatmap.2)
  mex.heatmap<- gedit("1", container = heatmap.2, coerce.with=as.numeric, width=10)
  height.heatmap.label <- glabel("height", container = heatmap.2)
  height.heatmap<- gedit("7", container = heatmap.2, coerce.with=as.numeric, width=10)
  width.heatmap.label <- glabel("width", container = heatmap.2)
  width.heatmap<- gedit("7", container = heatmap.2, coerce.with=as.numeric, width=10)
  margins.heatmap.label <- glabel("margins", container = heatmap.2)
  margins.heatmap<- gedit("c(10,10)", container = heatmap.2, coerce.with=svalue, width=10)
  
  ##Pathway Projector file parameters##
  pp.tab = ggroup(container=nb, label="Pathway Projector", horizontal=TRUE)
  
  pp.1 = gframe(container = pp.tab, horizontal = FALSE)
  treatment.label = glabel("Treatment Group", container = pp.1)
  treatment = gradio(levels(as.factor(data.pos$Group)), container = pp.1, selected = 1, use.table = TRUE, coerce.with=svalue )
  control.label = glabel("Control Group", container = pp.1)
  control = gradio(levels(as.factor(data.pos$Group)), container = pp.1, selected = 2, use.table = TRUE, coerce.with=svalue )
  dot.size.label = glabel("Dot Size", container = pp.1)
  dot.size = gedit("80", container = pp.1, coerce.with = svalue, width = 10)
  criterionPP.label = glabel("Criterion", container = pp.1, horizontal = FALSE)
  criterionPP = gradio(c("Pval", "Qval", "Fold-change", "Pval+Fold-change", "Qval+Fold-change"), 
                       selected=1, horizontal=FALSE, container = pp.1, coerce.with = svalue)
  PPchoice.label = glabel("Plot label", container = pp.1)
  PPchoice = gradio(c("Pval", "Qval"), container = pp.1, selected = 1, horizontal = TRUE, coerce.with = svalue) 
  
  
  pp.2 = gframe(container = pp.tab, horizontal = FALSE)
  
  t1.frame = gframe(container = pp.2, horizontal = FALSE)
  t1.label = glabel("Threshold 1", container = t1.frame)
  t1pq.frame = gframe(container = t1.frame, horizontal = TRUE)
  t1pq.label = glabel("P or Q", container = t1pq.frame)
  t1pq = gedit("0.05", container = t1pq.frame, coerce.with = as.numeric, width = 10)
  t1f.frame = gframe(container = t1.frame, horizontal = TRUE)
  t1f.label = glabel("Fold-change", container = t1f.frame)
  t1f = gedit("1.5", container = t1f.frame, coerce.with = as.numeric, width = 10)
  
  t2.frame = gframe(container = pp.2, horizontal = FALSE)
  t2.label = glabel("Threshold 2", container = t2.frame)
  t2pq.frame = gframe(container = t2.frame, horizontal = TRUE)
  t2pq.label = glabel("P or Q", container = t2pq.frame)
  t2pq = gedit("0.05", container = t2pq.frame, coerce.with = as.numeric, width = 10)
  t2f.frame = gframe(container = t2.frame, horizontal = TRUE)
  t2f.label = glabel("Fold-change", container = t2f.frame)
  t2f = gedit("2.0", container = t2f.frame, coerce.with = as.numeric, width = 10)
  
  t3.frame = gframe(container = pp.2, horizontal = FALSE)
  t3.label = glabel("Threshold 3", container = t3.frame)
  t3pq.frame = gframe(container = t3.frame, horizontal = TRUE)
  t3pq.label = glabel("P or Q", container = t3pq.frame)
  t3pq = gedit("0.05", container = t3pq.frame, coerce.with = as.numeric, width = 10)
  t3f.frame = gframe(container = t3.frame, horizontal = TRUE)
  t3f.label = glabel("Fold-change", container = t3f.frame)
  t3f = gedit("3.0", container = t3f.frame, coerce.with = as.numeric, width = 10)
  
  fontsizePP.label = glabel("Font Size", container = pp.2)
  fontsizePP = gedit("80", container = pp.2, coerce.with = svalue, width = 10)
  
  svalue(nb)=1
  
  btn.run <- gbutton(text = "Run", container = required.tab,
                     handler = function(...)	
                     { 	
                       criterion		= svalue(criterion)
                       plot.choice		= svalue(plot.choice)
                       log.base		<<- svalue(log.base)
                       f.cutoff		= svalue(f.cutoff)
                       p.cutoff		= svalue(p.cutoff)
                       
                       if (as.numeric(col.start) >= 5)
                       {	
                         random.factors	<<- svalue(random.factors)
                         fixed.continuous	<<- svalue(fixed.continuous)
                       }
                       output		= svalue(output)
                       
                       las.mean		= svalue(las.mean)	
                       cex.mean		= svalue(cex.mean)			
                       cex.axis.mean	= svalue(cex.axis.mean)			
                       cex.lab.mean	= svalue(cex.lab.mean)			 
                       cex.main.mean	= svalue(cex.main.mean)			
                       cex.sub.mean	= svalue(cex.sub.mean)			
                       mex.mean		= svalue(mex.mean)			
                       cex.mean.charsize = svalue(cex.mean.charsize)	
                       height.mean		= svalue(height.mean)
                       width.mean		= svalue(width.mean)				 
                       
                       las.prevpost	= svalue(las.prevpost)			
                       cex.prevpost	= svalue(cex.prevpost)
                       cex.axis.prevpost	= svalue(cex.axis.prevpost)
                       cex.lab.prevpost	= svalue(cex.lab.prevpost)			 
                       cex.main.prevpost	= svalue(cex.main.prevpost)
                       cex.sub.prevpost	= svalue(cex.sub.prevpost)
                       mex.prevpost	= svalue(mex.prevpost)
                       height.prevpost	= svalue(height.prevpost)
                       width.prevpost	= svalue(width.prevpost)
                       
                       las.summaryplot		= svalue(las.summaryplot)			
                       cex.summaryplot		= svalue(cex.summaryplot)
                       cex.axis.summaryplot	= svalue(cex.axis.summaryplot)
                       cex.lab.summaryplot	= svalue(cex.lab.summaryplot)
                       cex.main.summaryplot	= svalue(cex.main.summaryplot)
                       cex.sub.summaryplot	= svalue(cex.sub.summaryplot)
                       mex.summaryplot		= svalue(mex.summaryplot)		
                       height.summaryplot	= svalue(height.summaryplot)
                       width.summaryplot		= svalue(width.summaryplot)
                       
                       las.summaryplot.labeled		= svalue(las.summaryplot.labeled)			
                       cex.summaryplot.labeled 	= svalue(cex.summaryplot.labeled)
                       cex.axis.summaryplot.labeled	= svalue(cex.axis.summaryplot.labeled)		
                       cex.lab.summaryplot.labeled	= svalue(cex.lab.summaryplot.labeled)
                       cex.main.summaryplot.labeled	= svalue(cex.main.summaryplot.labeled)
                       cex.sub.summaryplot.labeled	= svalue(cex.sub.summaryplot.labeled)			
                       mex.summaryplot.labeled		= svalue(mex.summaryplot.labeled)			
                       height.summaryplot.labeled 	= svalue(height.summaryplot.labeled) 		 
                       width.summaryplot.labeled 	= svalue(width.summaryplot.labeled)
                       
                       mfrow.mean			= svalue(mfrow.mean)		
                       mar.mean			= svalue(mar.mean)	 
                       mfrow.prevpost		= svalue(mfrow.prevpost)	
                       mar.prevpost		= svalue(mar.prevpost)
                       mar.summaryplot		= svalue(mar.summaryplot)	
                       mar.summaryplot.labeled = svalue(mar.summaryplot.labeled)
                       
                       location.choice 			= svalue(location.choice) 
                       trim.alpha.choice 		= svalue(trim.alpha.choice)
                       bootstrap.choice 			= svalue(bootstrap.choice)
                       num.bootstrap.choice 		= svalue(num.bootstrap.choice)
                       kruskal.test.choice 		= svalue(kruskal.test.choice)
                       correction.method.choice 	= svalue(correction.method.choice)
                       
                       margins.heatmap	= svalue(margins.heatmap)
                       las.heatmap		= svalue(las.heatmap)			
                       cex.heatmap		= svalue(cex.heatmap)
                       cex.axis.heatmap	= svalue(cex.axis.heatmap)
                       cex.lab.heatmap	= svalue(cex.lab.heatmap)
                       cex.main.heatmap	= svalue(cex.main.heatmap)
                       cex.sub.heatmap	= svalue(cex.sub.heatmap)
                       mex.heatmap		= svalue(mex.heatmap)		
                       height.heatmap	= svalue(height.heatmap)
                       width.heatmap	= svalue(width.heatmap)
                       
                       treatment		= svalue(treatment)
                       control		= svalue(control)
                       dot.size		= svalue(dot.size)
                       fontsizePP		= svalue(fontsizePP)
                       criterionPP		= svalue(criterionPP)
                       PPchoice		= svalue(PPchoice)
                       t1pq 			= svalue(t1pq)
                       t1f			= svalue(t1f)
                       t2pq 			= svalue(t2pq)
                       t2f			= svalue(t2f)
                       t3pq 			= svalue(t3pq)
                       t3f			= svalue(t3f)
                       ###############################################################################################################################
                       ###############################################################################################################################
                       
                       time=paste(paste(Sys.Date(), "_", format(Sys.time(),"%H%M"), sep="" ) )
                       
                       
                       #######
                       ## 1 ##
                       #######
                       ##Log-transform ion counts if specified by the user.
                       ##########################################################################
                       
                       if (mode(log.base)=="numeric" && !is.na(log.base))
                       {
                         data.pos=cbind(data.pos[,1:(col.start-2)], log( as.matrix(data.pos[,(col.start-1):ncol(data.pos)]), log.base ))
                         if (mode(data.neg)=="list")
                         {data.neg=cbind(data.neg[,1:(col.start-2)], log( as.matrix(data.neg[,(col.start-1):ncol(data.neg)]), log.base )) }
                       } 
                       
                       #######
                       ## 2 ##
                       #######
                       ##Use user-specified random- and fixed-effect model terms to construct either a lmer or lm model.  Use 
                       ##function "lmer" from package "lme4", or "lm" for each metabolite, predicting metabolite level using 
                       ##fixed-effect and random-effect variables; then add residuals to each group mean.
                       ################################################################################################ 
                       ##normalize negative mode metabolites 
                       
                       library(lme4)
                       data.adj.neg<-c()
                       
                       if (mode(data.neg)=="list")
                       {
                         group.test.1=c()
                         names.neg.2=c()
                         names.neg=colnames(data.neg[(col.start-1):ncol(data.neg)])
                         neg.good=c()
                         neg.bad=c()
                         for (i in (col.start-1):ncol(data.neg))
                         {
                           lmer.data<-as.data.frame(data.neg[,c(1:(col.start-2),i)])
                           colnames(lmer.data) = c(colnames(data.neg[1:(col.start-2)]), "Metabolite")
                           
                           #######################################################################################################################
                           
                           if (length(fixed.continuous)!=0 && as.numeric(col.start) >= 5)
                           {
                             fixed.continuous.model="+"
                             for (yy in 1:length(fixed.continuous))
                             {
                               if (yy<length(fixed.continuous))
                               {
                                 fixed.continuous.model=paste(fixed.continuous.model, "as.numeric(as.vector(", fixed.continuous[yy], "))", "+", sep="")
                               }else
                               {
                                 fixed.continuous.model=paste(fixed.continuous.model, "as.numeric(as.vector(", fixed.continuous[yy], "))", sep="")
                               }
                             }	
                           }else
                           {
                             fixed.continuous.model=c()	
                           }
                           
                           #if ( mode(fixed.categorical)=="character" )
                           #{
                           #	fixed.categorical.model=c("+")
                           #	for (yy in 1:length(fixed.categorical))
                           #	{
                           #		if (yy<length(fixed.categorical))
                           #			{
                           #				fixed.categorical.model=paste(fixed.categorical.model, "as.factor(", fixed.categorical[yy], ")", "+", sep="")
                           #			}else
                           #			{
                           #				fixed.categorical.model=paste(fixed.categorical.model, "as.factor(", fixed.categorical[yy], ")", sep="")
                           #			}
                           #	}	
                           #}else
                           #{
                           #	fixed.categorical.model=c()
                           #}
                           
                           if ( length(random.factors)!=0 && as.numeric(col.start) >= 5)
                           {
                             r.include=c()
                             for (r in random.factors)
                             {
                               r.levels.good=c()
                               for (rr in levels(lmer.data[,r]))
                               {
                                 if ( all(is.na(lmer.data[which(lmer.data[,r]==rr),"Metabolite"]))=="FALSE" )	
                                 {r.levels.good=c(r.levels.good, rr)}
                                 
                               }
                               if (length(r.levels.good)>=2)
                               {r.include=c(r.include, r)}
                               
                             }
                             if (mode(r.include)=="character")
                             {
                               random.model="+"
                               for (yy in 1:length(r.include))
                               {
                                 if (yy<length(r.include))
                                 {
                                   random.model=paste(random.model, "(1|", r.include[yy], ")", "+", sep="")
                                 }else
                                 {
                                   random.model=paste(random.model, "(1|", r.include[yy], ")", sep="")
                                 }
                                 
                               }
                               
                               model.type="LMER"
                             }else	{
                               random.model=c()
                               model.type="LM"
                             }	
                             
                             
                             
                             
                           }else{
                             random.model=c()
                             model.type="LM"
                           }
                           linear.model<<-paste("as.numeric(as.vector(Metabolite)) ~ 1", "+ as.factor(Group)", fixed.continuous.model, random.model, sep="")
                           model.type<<-model.type
                           #######################################################################################################################
                           options(show.error.messages=FALSE)
                           if (model.type=="LMER") 
                           {
                             fit = try( lmer(formula=linear.model, data=lmer.data, contrasts=FALSE, model=TRUE) )
                           }else { fit = try( lm(formula=linear.model, data=lmer.data) ) }
                           
                           options(show.error.messages=TRUE)
                           
                           if (mode(fit)=="character")
                           {	
                             neg.bad=c(neg.bad, i)
                           }else	{
                             check.result="yes"
                             if (model.type=="LMER")
                             {
                               if ( is.na(fit@deviance[2]) ) {check.result="no" } 
                             }else	{ 
                               if ( all(residuals(fit)==0) )  { check.result="no" } 
                             }	
                             n=c()
                             for (j in levels(lmer.data$Group))
                             {
                               n = c(n,length(which(!is.na(lmer.data[which(lmer.data$Group==j),"Metabolite"]))) )
                             }
                             
                             
                             
                             if (check.result=="no" )
                             {
                               neg.bad=c(neg.bad,i)
                             }else	{
                               neg.good=c(neg.good, i)				
                               fit.res=vector(mode="character", nrow(lmer.data))
                               fit.res[which(!is.na(data.neg[,i]))]=residuals(fit)
                               fit.res=as.numeric(as.vector(replace(fit.res, fit.res=="", NA)))
                               group.one.non.na=levels(as.factor(as.character(lmer.data$Group[which(!is.na(lmer.data$Metabolite))])))   
                               
                               if (model.type=="LMER") 
                               {
                                 treat.means=as.vector( c(fit@fixef[1], fit@fixef[2:length(group.one.non.na)]+fit@fixef[1]) )
                               }else	{
                                 treat.means=as.vector( c(fit$coefficients[1], fit$coefficients[2:length(group.one.non.na)]+fit$coefficients[1]) )
                               }
                               treat.means.all=vector(mode="character", length=nrow(lmer.data))
                               
                               for (hh in 1:length(group.one.non.na))
                               {
                                 treat.means.all[which(lmer.data$Group==group.one.non.na[hh])] = treat.means[hh] 
                               }
                               treat.means.all = as.numeric(replace(treat.means.all, treat.means.all=="", NA))
                               
                               data.adj.neg<-cbind( data.adj.neg, fit.res+treat.means.all )
                             }
                           }		
                           
                         }
                         colnames(data.adj.neg)=colnames(data.neg)[neg.good]
                         
                       }
                       
                       ##normalize positive mode metabolites 
                       
                       group.test.1=c()
                       names.pos.2=c()
                       data.adj.pos<-c()
                       names.pos=colnames(data.pos[(col.start-1):ncol(data.pos)])
                       pos.good=c()
                       pos.bad=c()
                       
                       for (i in (col.start-1):ncol(data.pos))
                       {
                         lmer.data<-as.data.frame(data.pos[,c(1:(col.start-2),i)])
                         colnames(lmer.data) = c(colnames(data.pos[1:(col.start-2)]), "Metabolite")
                         
                         #######################################################################################################################
                         
                         if (length(fixed.continuous)!=0 && as.numeric(col.start) >= 5)
                         {
                           fixed.continuous.model="+"
                           for (yy in 1:length(fixed.continuous))
                           {
                             if (yy<length(fixed.continuous))
                             {
                               fixed.continuous.model=paste(fixed.continuous.model, "as.numeric(as.vector(", fixed.continuous[yy], "))", "+", sep="")
                             }else
                             {
                               fixed.continuous.model=paste(fixed.continuous.model, "as.numeric(as.vector(", fixed.continuous[yy], "))", sep="")
                             }
                           }	
                         }else
                         {
                           fixed.continuous.model=c()	
                         }
                         
                         #if ( mode(fixed.categorical)=="character" )
                         #{
                         #	fixed.categorical.model=c("+")
                         #	for (yy in 1:length(fixed.categorical))
                         #	{
                         #		if (yy<length(fixed.categorical))
                         #			{
                         #				fixed.categorical.model=paste(fixed.categorical.model, "as.factor(", fixed.categorical[yy], ")", "+", sep="")
                         #			}else
                         #			{
                         #				fixed.categorical.model=paste(fixed.categorical.model, "as.factor(", fixed.categorical[yy], ")", sep="")
                         #			}
                         #	}	
                         #}else
                         #{
                         #	fixed.categorical.model=c()
                         #}
                         
                         if ( length(random.factors)!=0 && as.numeric(col.start) >= 5)
                         {
                           r.include=c()
                           for (r in random.factors)
                           {
                             r.levels.good=c()
                             for (rr in levels(lmer.data[,r]))
                             {
                               if ( all(is.na(lmer.data[which(lmer.data[,r]==rr),"Metabolite"]))=="FALSE" )	
                               {r.levels.good=c(r.levels.good, rr)}
                               
                             }
                             if (length(r.levels.good)>=2)
                             {r.include=c(r.include, r)}
                             
                           }
                           if (mode(r.include)=="character")
                           {
                             random.model="+"
                             for (yy in 1:length(r.include))
                             {
                               if (yy<length(r.include))
                               {
                                 random.model=paste(random.model, "(1|", r.include[yy], ")", "+", sep="")
                               }else
                               {
                                 random.model=paste(random.model, "(1|", r.include[yy], ")", sep="")
                               }
                               
                             }
                             
                             model.type="LMER"
                           }else	{
                             random.model=c()
                             model.type="LM"
                           }	
                           
                           
                           
                           
                         }else{
                           random.model=c()
                           model.type="LM"
                         }
                         linear.model<<-paste("as.numeric(as.vector(Metabolite)) ~ 1", "+ as.factor(Group)", fixed.continuous.model, random.model, sep="")
                         
                         #######################################################################################################################
                         options(show.error.messages=FALSE)
                         if (model.type=="LMER") 
                         {
                           fit = try( lmer(formula=linear.model, data=lmer.data, contrasts=FALSE, model=TRUE) )
                         }else { fit = try( lm(formula=linear.model, data=lmer.data) ) }
                         
                         options(show.error.messages=TRUE)
                         
                         if (mode(fit)=="character")
                         {	
                           pos.bad=c(pos.bad, i)
                         }else	{
                           check.result="yes"
                           if (model.type=="LMER")
                           {
                             if ( is.na(fit@deviance[2])) {check.result="no" } 
                           }else	{ 
                             if ( all(residuals(fit)==0) )  { check.result="no" } 
                           }	
                           n=c()
                           for (j in levels(lmer.data$Group))
                           {
                             n = c(n,length(which(!is.na(lmer.data[which(lmer.data$Group==j),"Metabolite"]))) )
                           }
                           
                           
                           
                           if (check.result=="no" )
                           {
                             pos.bad=c(pos.bad,i)
                           }else	{	pos.good=c(pos.good, i)				
                           fit.res=vector(mode="character", nrow(lmer.data))
                           fit.res[which(!is.na(data.pos[,i]))]=residuals(fit)
                           fit.res=as.numeric(as.vector(replace(fit.res, fit.res=="", NA)))
                           group.one.non.na=levels(as.factor(as.character(lmer.data$Group[which(!is.na(lmer.data$Metabolite))])))   
                           
                           if (model.type=="LMER") 
                           {
                             treat.means=as.vector( c(fit@fixef[1], fit@fixef[2:length(group.one.non.na)]+fit@fixef[1]) )
                           }else	{
                             treat.means=as.vector( c(fit$coefficients[1], fit$coefficients[2:length(group.one.non.na)]+fit$coefficients[1]) )
                           }
                           treat.means.all=vector(mode="character", length=nrow(lmer.data))
                           
                           for (hh in 1:length(group.one.non.na))
                           {
                             treat.means.all[which(lmer.data$Group==group.one.non.na[hh])] = treat.means[hh] 
                           }
                           treat.means.all = as.numeric(replace(treat.means.all, treat.means.all=="", NA))
                           
                           data.adj.pos<-cbind( data.adj.pos, fit.res+treat.means.all )
                           }
                         }		
                         
                       }
                       colnames(data.adj.pos)=colnames(data.pos)[pos.good]
                       
                       if (mode(data.neg)=="list")	
                       { 
                         data.adj<-cbind(as.character(data.neg$Subject), data.adj.neg, data.adj.pos) 
                         colnames(data.adj)=c( "Subject", colnames(data.adj[,2:ncol(data.adj)]) )
                         row.names(data.adj)=row.names(data.pos)
                       } else 
                       {
                         #data.adj=data.adj.pos
                         #data.adj<-cbind(as.character(data.pos$Subject), data.adj.neg, data.adj.pos) 
                         data.adj<-cbind(as.character(data.pos$Subject), data.adj.pos) 
                         
                         colnames(data.adj)=c("Subject", colnames(data.adj[,2:ncol(data.adj)]))
                         row.names(data.adj)=row.names(data.pos)
                       }
                       
                       write.csv(data.adj, file=paste(output, "_normalized data_", time, ".csv", sep=""))
                       
                       #######
                       ## 3 ##
                       #######
                       ##Average the pairs of rows that contained replicates, combine data into one object, then write to csv file.
                       ############################################################################################################
                       compile.vector=c()
                       data.adj=as.data.frame(data.adj)
                       subjects=levels(as.factor(data.adj$Subject))
                       groups=c()
                       for (i in 1:length(subjects))
                       {
                         group.cats=as.matrix( data.pos[which(data.adj$Subject==subjects[i]), 2] )
                         group.cats=as.vector(group.cats[1,])
                         groups=rbind(groups, group.cats)
                         
                         for (j in 2:ncol(data.adj))
                         {
                           rep.vals=as.vector( data.adj[which(data.adj$Subject==subjects[i]),j] )
                           compile.vector=c( compile.vector, mean(as.numeric(rep.vals), na.rm=TRUE) )
                         }
                       }
                       
                       data.compile=matrix(compile.vector, length(subjects), ncol(data.adj)-1, byrow=TRUE )
                       data.compile=as.data.frame( cbind(groups, data.compile) )
                       rownames(data.compile)=subjects
                       colnames(data.compile)=c( colnames(data.pos[2]), colnames(data.adj[2:ncol(data.adj)]) )
                       
                       ##write normalized data to csv file
                       write.csv(data.compile,file=paste(output, "_normalized and compiled data_", time, ".csv", sep=""))
                       
                       #######
                       ## 4 ##
                       #######
                       ##Perform Levene's test for equality of variances among treatment groups.
                       #########################################################################
                       library("lawstat")
                       
                       levene.vector=c()
                       
                       levene.matrix=c()
                       
                       missing.metabolites=c()
                       complete.metabolites=c()
                       levene.na=c()
                       levene.good=c()
                       for (i in 2:ncol(data.compile))
                       {
                         
                         
                         complete.metabolites=c(complete.metabolites, colnames(data.compile[i]))
                         options(show.error.messages=FALSE)
                         levene.met = try ( levene.test( as.numeric(as.vector(data.compile[,i])), data.compile$Group, 
                                                         location=location.choice, trim.alpha=trim.alpha.choice, bootstrap = bootstrap.choice, 
                                                         num.bootstrap=num.bootstrap.choice, kruskal.test=kruskal.test.choice, 
                                                         correction.method=correction.method.choice ) )
                         options(show.error.messages=TRUE)
                         
                         if (mode(levene.met)=="character" || is.na(levene.met$statistic)) 
                         {	levene.na=c(levene.na, i) 
                         }else	{
                           levene.matrix=rbind(levene.matrix, c(levene.met$statistic, levene.met$p.value))
                           levene.good=c(levene.good,i)	
                         }						
                       }	
                       
                       #write.csv(data.compile[,c(1,levene.na)], file="levene.na.csv")
                       
                       levene.dataframe=as.data.frame(levene.matrix)
                       row.names(levene.dataframe)=colnames(data.compile)[levene.good]
                       colnames(levene.dataframe)=c("Test Statistic", "P-value")
                       
                       non.equal.var=c()
                       non.equal.var_names=c()
                       for (i in 1:nrow(levene.dataframe))
                       {
                         
                         
                         if (levene.dataframe[i,2]<0.05)
                         {
                           non.equal.var=rbind(non.equal.var, levene.dataframe[i,])
                           non.equal.var_names=c( non.equal.var_names, row.names(levene.dataframe)[i] )
                         }		
                       }
                       
                       if (mode(non.equal.var)=="list")
                       {
                         colnames(non.equal.var)=c("Test Statistic", "P-value")
                         non.equal.var$"P-value"=format(non.equal.var$"P-value", scientific=FALSE)
                       }
                       
                       detach("package:lawstat")
                       detach("package:VGAM") 
                       detach("package:stats4")
                       
                       #######
                       ## 5 ##
                       #######
                       ##Construct summary plot of mean metabolite ("x") vs. residual error ("y").
                       ###########################################################################
                       ##create vector of metabolite means and matrix of residuals from data.compile
                       mean.vector=c()
                       resid.matrix=c()
                       for (i in 2:ncol(data.compile))
                       {
                         mean.vector=c(mean.vector, mean(as.numeric(as.vector(data.compile[,i])), na.rm=TRUE) )
                         
                         resid.matrix=cbind(resid.matrix, as.numeric(as.vector(data.compile[,i]))-
                                              mean(as.numeric(as.vector(data.compile[,i])),na.rm=TRUE) )
                       }
                       
                       ##create matrix with metabolite means repeated in every row
                       mean.matrix=c()
                       for (i in 1:nrow(data.compile))
                       {
                         mean.matrix=rbind(mean.matrix, mean.vector)	
                       }
                       
                       
                       colnames(resid.matrix)=colnames(data.compile[2:ncol(data.compile)])
                       colnames(mean.matrix)=colnames(data.compile[2:ncol(data.compile)])
                       
                       
                       colnames(resid.matrix)=colnames(data.compile[2:ncol(data.compile)])
                       overall.mean = mean.matrix[1,]
                       overall.mean = sort(overall.mean, decreasing=FALSE)
                       mean.names = names(overall.mean)
                       
                       
                       mean.matrix.sort = mean.matrix[, order(mean.matrix[1,], decreasing=FALSE)]
                       resid.matrix.sort = resid.matrix[, order(mean.matrix[1,], decreasing=FALSE)]
                       
                       plot.c=vector(length=length(data.compile$Group))
                       plot.s=vector(length=length(data.compile$Group))
                       rainbow.vector=rainbow(length(levels(as.factor(data.compile$Group))))
                       
                       for (i in 1:length(levels(as.factor(data.compile$Group))))
                       {
                         plot.s[which(data.compile$Group==levels(as.factor(data.compile$Group))[i] ) ] = i
                         plot.c[which(data.compile$Group==levels(as.factor(data.compile$Group))[i] ) ] = rainbow.vector[i]
                       }
                       
                       
                       pdf(file=paste(output, "_residual plot_", time, ".pdf", sep=""), width=width.summaryplot, height=height.summaryplot)
                       par(	#mar=mar.summaryplot, 
                         las=las.summaryplot, 
                         cex=cex.summaryplot, 
                         cex.axis=cex.axis.summaryplot, 
                         cex.lab=cex.lab.summaryplot, 
                         cex.main=cex.main.summaryplot, 
                         cex.sub=cex.sub.summaryplot, 
                         mex=mex.summaryplot)
                       
                       layout(matrix(c(1,2), nrow=1), widths=c(1/2, 1/2))
                       
                       
                       plot(mean.matrix.sort, resid.matrix.sort,
                            pch=plot.s,
                            col=plot.c, 
                            cex.main=cex.main.summaryplot,
                            cex.lab=cex.lab.summaryplot,
                            cex.axis=cex.axis.summaryplot,
                            main="Metabolite Mean vs. Residual Error", 
                            xlab="Metabolite Mean", 
                            ylab="Residual Error" )
                       
                       abline(a=0,b=0)
                       
                       plot(mean.matrix.sort, resid.matrix.sort,
                            type="n", axes=FALSE, ann=FALSE, xpd=NA )
                       
                       legend(	"topleft", 
                               legend=levels(as.factor(data.compile$Group)), 
                               text.col=rainbow.vector,
                               pch=1:length(levels(as.factor(data.compile$Group))),
                               col=rainbow.vector )
                       
                       dev.off()
                       
                       pdf(file=paste(output, "_residual plot labeled_", time, ".pdf", sep=""), width=width.summaryplot.labeled, height=height.summaryplot.labeled)
                       par(	#mar=mar.summaryplot.labeled, 
                         las=las.summaryplot.labeled, 
                         cex=cex.summaryplot.labeled, 
                         cex.axis=cex.axis.summaryplot.labeled, 
                         cex.lab=cex.lab.summaryplot.labeled, 
                         cex.main=cex.main.summaryplot.labeled, 
                         cex.sub=cex.sub.summaryplot.labeled,
                         mex=mex.summaryplot.labeled)
                       
                       layout(matrix(c(1,2), nrow=1), widths=c(1/2, 1/2))
                       
                       
                       plot(	mean.matrix.sort, 
                             resid.matrix.sort, 
                             pch=plot.s,
                             col=plot.c, 
                             cex.main=cex.main.summaryplot.labeled,
                             cex.lab=cex.lab.summaryplot.labeled,
                             cex.axis=cex.axis.summaryplot.labeled,
                             main="Metabolite Mean vs. Residual Error", 
                             xlab="Metabolite Mean", 
                             ylab="Residual Error" )
                       
                       axis(side=1, at=overall.mean, labels=mean.names)
                       abline(a=0,b=0)
                       
                       plot(mean.matrix.sort, resid.matrix.sort,
                            type="n", axes=FALSE, ann=FALSE, xpd=NA )
                       
                       legend(	"topleft", 
                               legend=levels(as.factor(data.compile$Group)), 
                               text.col=rainbow.vector,
                               pch=1:length(levels(as.factor(data.compile$Group))),
                               col=rainbow.vector )
                       
                       dev.off()
                       
                       
                       #########################################################################################################################
                       
                       #######
                       ## 6 ##
                       #######
                       ##Perform Shapiro-Wilk test of normality.
                       #########################################
                       shapiro.matrix=c()
                       shapiro.bad=c()
                       shapiro.good=c()
                       options(show.error.messages=FALSE)
                       
                       for (i in 2:ncol(data.compile))
                       {
                         shapiro.met=try( shapiro.test(as.numeric(as.vector(data.compile[,i]))) )
                         if (mode(shapiro.met)=="list")
                         {
                           shapiro.matrix=rbind(shapiro.matrix, c(shapiro.met$statistic, shapiro.met$p.value))
                           shapiro.good=c(shapiro.good, colnames(data.compile)[i])
                         }else	{shapiro.bad=c(shapiro.bad, i)}
                       }
                       
                       options(show.error.messages=TRUE)
                       
                       shapiro.dataframe=as.data.frame(shapiro.matrix)
                       row.names(shapiro.dataframe)=colnames(data.compile[shapiro.good])
                       colnames(shapiro.dataframe)=c("W", "P-value")
                       
                       non.normal.dataframe=c()
                       non.normal.matrix=c()
                       for (i in 1:nrow(shapiro.dataframe))
                       {
                         if (as.numeric(shapiro.dataframe$'P-value'[i])<0.05)
                         {
                           non.normal.matrix=rbind(non.normal.matrix, shapiro.dataframe[i,])
                         }
                       }
                       if (mode(non.normal.matrix)=="list")
                       {
                         non.normal.dataframe=as.data.frame(non.normal.matrix)
                         non.normal.dataframe$"P-value"=format(c(non.normal.dataframe$"P-value"), scientific=FALSE)
                       }
                       ##############################################################################################
                       
                       #######
                       ## 7 ##
                       #######
                       ##Loop through metabolites, performing ANOVA in each loop with Tukey HSD post-hoc comparison.  The user specifies a p-value 
                       ##or fold-change cutoff value in parameter #8.  Any metabolite that had at least one between-group Tukey post-hoc comparison 
                       ##with a p-value below the cutoff or a fold-change above the cutoff will be printed in a list.   
                       ############################################################################################################################  
                       library(gplots)
                       
                       metabolites=colnames(data.compile[2:ncol(data.compile)])
                       
                       amodel=as.numeric(as.vector(Metabolite)) ~ as.factor(Group) 
                       
                       tukey.sig.metabolite=c()
                       tukey.sig.pval=c()
                       tukey.sig.qval=c()
                       tukey.sig.foldch=c()
                       tukey.all.pval=c()
                       
                       foldchange.sig.metabolite=c()
                       foldchange.sig.pval=c()
                       foldchange.sig.qval=c()
                       foldchange.sig.foldch=c()
                       foldchange.all.foldch=c()
                       
                       qvalue.sig.metabolite=c()
                       qvalue.sig.qval=c()
                       qvalue.sig.pval=c()
                       qvalue.sig.foldch=c()
                       
                       pf.sig.metabolite=c()
                       pf.sig.pval=c()
                       pf.sig.qval=c()
                       pf.sig.foldch=c()
                       
                       qf.sig.metabolite=c()
                       qf.sig.qval=c()
                       qf.sig.pval=c()
                       qf.sig.foldch=c()
                       
                       mean.plot.data=c()
                       
                       aov.good=c()
                       aov.bad=c()
                       
                       for (i in 2:ncol(data.compile)) 
                       {
                         aov.data=data.compile[,c(1,i)]
                         colnames(aov.data)=c("Group", "Metabolite")
                         row.names.vector=c()
                         
                         ##replace Nan and Inf with NA for simplicity
                         NA.vector=as.numeric(as.vector(aov.data$Metabolite))
                         NA.vector=replace(NA.vector, which(is.nan(NA.vector)), NA )
                         NA.vector=replace( NA.vector, which(is.infinite(NA.vector)), NA)
                         
                         ##if any treatment groups have fewer than 3 data points, change all 
                         ##measurements to NA and anova will not be done
                         n=c()
                         for (j in levels(as.factor(aov.data$Group)))
                         {
                           level.na = length( which(!is.na(NA.vector[which(aov.data$Group==j)])) )
                           n = c(n, as.numeric(as.vector(level.na)) )
                           if (level.na <= 2)
                           {
                             NA.vector=replace(NA.vector, which(as.character(aov.data$Group)==j), NA )
                           }
                         }
                         
                         aov.data=cbind(aov.data[1:(ncol(aov.data)-1)], NA.vector )
                         colnames(aov.data)[ncol(aov.data)]="Metabolite"
                         
                         ##anova
                         options(show.error.messages=FALSE)
                         fit=try( aov(formula=amodel, data=aov.data, projections=TRUE) ) 
                         options(show.error.messages=TRUE)
                         
                         ##if anova didn't return error and at least 2 treatment groups have 3+ measurements,
                         ##then determine which of these metabolites have significant fold-changes and p- or q-values
                         p=NaN
                         if (mode(fit)=="list")
                         {	tukey.matrix<-as.matrix(TukeyHSD(fit)[[1]]) 
                         p=tukey.matrix[,4]
                         }
                         
                         if (any(!is.nan(p))  && length(which(n>=3))>=2)
                         {
                           ##threecheck--make sure both treatment groups being compared have at least 3 data points
                           threecheck=c()
                           aov.good=c(aov.good, i)
                           tukey.matrix<-as.matrix(TukeyHSD(fit)[[1]])
                           
                           ##if user selected log-transformation:
                           if (mode(log.base)=="numeric" && !is.na(log.base)) 
                           {
                             foldch.dataframe=c()
                             aov.data.backtransformed=as.data.frame( cbind( as.character(aov.data$Group), 
                                                                            log.base^as.numeric(as.vector(aov.data$Metabolite)) ) )
                             colnames(aov.data.backtransformed)=c("Group", "Metabolite")
                             row.names(aov.data.backtransformed)=row.names(data.compile)		
                             
                             ##two loops to calculate pairwise mean fold-changes 
                             for (gg in 1:(length(levels(as.factor(aov.data.backtransformed$Group)))-1) )
                             {
                               for (hh in (gg+1):length(levels(as.factor(aov.data.backtransformed$Group))))
                               {
                                 row.names.vector=c(row.names.vector, paste(levels(as.factor(aov.data.backtransformed$Group))[hh], "-", levels(as.factor(aov.data.backtransformed$Group))[gg], sep=""))
                                 vector.h = as.numeric(as.vector(aov.data.backtransformed[which(aov.data.backtransformed$Group==levels(as.factor(aov.data.backtransformed$Group))[hh]),2]))
                                 mean.h=mean( vector.h, na.rm=TRUE )
                                 n.h = length (which(!is.na(vector.h)) )
                                 vector.g = as.numeric(as.vector(aov.data.backtransformed[which(aov.data.backtransformed$Group==levels(as.factor(aov.data.backtransformed$Group))[gg]),2]))
                                 mean.g=mean( vector.g, na.rm=TRUE )
                                 n.g = length (which(!is.na(vector.g)) )
                                 
                                 if (any(is.na(c(mean.h,mean.g)))) 
                                 {	foldch.dataframe=c(foldch.dataframe, NA)
                                 } else {foldch.dataframe=c(foldch.dataframe, mean.h/mean.g)}
                                 if (n.h<=2 || n.g<=2){threecheck=c(threecheck,"bad") }else{threecheck=c(threecheck,"good")}
                               }
                             }
                             foldch.dataframe=as.data.frame(foldch.dataframe)
                             
                             ##if user did not select log-transformation:	
                           } else {
                             foldch.dataframe=c()
                             for (gg in 1:(length(levels(as.factor(aov.data$Group)))-1) )
                             {
                               for (hh in (gg+1):length(levels(as.factor(aov.data$Group))))
                               {
                                 row.names.vector=c(row.names.vector, paste(levels(as.factor(aov.data$Group))[hh], "-", levels(as.factor(aov.data$Group))[gg], sep=""))
                                 vector.h = as.numeric(as.vector(aov.data[which(aov.data$Group==levels(as.factor(aov.data$Group))[hh]),2]))
                                 mean.h=mean( vector.h, na.rm=TRUE )
                                 n.h = length (which(!is.na(vector.h)) )
                                 vector.g = as.numeric(as.vector(aov.data[which(aov.data$Group==levels(as.factor(aov.data$Group))[gg]),2]))
                                 mean.g=mean( vector.g, na.rm=TRUE )
                                 n.g = length (which(!is.na(vector.g)) )
                                 
                                 
                                 if (any(is.na(c(mean.h,mean.g)))) 
                                 {	foldch.dataframe=c(foldch.dataframe, NA)
                                 } else {foldch.dataframe=c(foldch.dataframe, mean.h/mean.g)}
                                 if (n.h<=2 || n.g<=2){threecheck=c(threecheck,"bad") }else{threecheck=c(threecheck,"good")}
                                 
                                 
                               }
                             }
                             foldch.dataframe=as.data.frame(foldch.dataframe)
                           }
                           
                           colnames(foldch.dataframe)="Fold-change"
                           row.names(foldch.dataframe)=row.names.vector
                           
                           tukey.matrix.withNA=matrix(data="", nrow=nrow(foldch.dataframe), ncol=ncol(tukey.matrix), dimnames=list(row.names(foldch.dataframe), colnames(tukey.matrix)))
                           tukey.matrix.withNA[which(!is.na(foldch.dataframe$"Fold-change")), ] = tukey.matrix	
                           tukey.matrix=tukey.matrix.withNA
                           
                           tukey.matrix[which(as.character(threecheck)=="bad"), 4] = NA
                           
                           tukey.all.pval=rbind(tukey.all.pval, as.numeric(as.vector(tukey.matrix[,ncol(tukey.matrix)])))
                           foldchange.all.foldch=rbind(foldchange.all.foldch, as.numeric(as.vector(foldch.dataframe[,1])))
                           
                         }else	{ aov.bad=c(aov.bad, colnames(data.compile)[i]) }
                       }
                       
                       
                       row.names(tukey.all.pval)=colnames(data.compile)[aov.good]
                       colnames(tukey.all.pval)=row.names(tukey.matrix)
                       tukey.all.pval.dataframe=as.data.frame(tukey.all.pval)
                       tukey.all.pval.dataframe=format(tukey.all.pval.dataframe, scientific=FALSE)
                       
                       row.names(foldchange.all.foldch)=colnames(data.compile)[aov.good]
                       colnames(foldchange.all.foldch)=row.names(foldch.dataframe)
                       foldchange.all.foldch.dataframe=as.data.frame(foldchange.all.foldch)
                       foldchange.all.foldch.dataframe=format(foldchange.all.foldch.dataframe, scientific=FALSE)
                       
                       library(qvalue)
                       
                       qgood=c()
                       qvalues = matrix(NA, nrow(tukey.all.pval.dataframe), ncol(tukey.all.pval.dataframe))
                       colnames(qvalues) = colnames(tukey.all.pval.dataframe)
                       row.names(qvalues) = row.names(tukey.all.pval.dataframe)
                       
                       for (i in 1:ncol(tukey.all.pval.dataframe))
                       {
                         pcol = tukey.all.pval.dataframe[,i]
                         ppos = which(!is.na(as.numeric(as.vector(pcol))))
                         p = as.numeric(as.vector(pcol[ppos]))
                         
                         ###########################################################################################################################
                         lambda = seq(0, 0.9, 0.05)
                         pi0.method = "smoother"
                         fdr.level = NULL
                         robust = FALSE
                         gui = FALSE
                         smooth.df = 3
                         smooth.log.pi0 = FALSE 
                         
                         m <- length(p)
                         
                         pi0 <- rep(0, length(lambda))
                         for (j in 1:length(lambda)) {
                           pi0[j] <- mean(p >= lambda[j])/(1 - lambda[j])
                         }
                         if (smooth.log.pi0) 
                           pi0 <- log(pi0)
                         spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
                         pi0 <- predict(spi0, x = max(lambda))$y
                         if (smooth.log.pi0) 
                           pi0 <- exp(pi0)
                         pi0 <- min(pi0, 1)
                         
                         if (pi0 > 0) 
                         {
                           u <- order(p)
                           qvalue.rank <- function(x) {
                             idx <- sort.list(x)
                             fc <- factor(x)
                             nl <- length(levels(fc))
                             bin <- as.integer(fc)
                             tbl <- tabulate(bin)
                             cs <- cumsum(tbl)
                             tbl <- rep(cs, tbl)
                             tbl[idx] <- tbl
                             return(tbl)
                           }
                           
                           v <- qvalue.rank(p)
                           qvalue <- pi0 * m * p/v
                           if (robust) {
                             qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
                           }
                           
                           qvalue[u[m]] <- min(qvalue[u[m]], 1)
                           for (ii in (m - 1):1) {
                             qvalue[u[ii]] <- min(qvalue[u[ii]], qvalue[u[ii + 1]], 1)
                           }
                           
                           retval <- list(pi0 = pi0, qvalues = qvalue, 
                                          pvalues = p, lambda = lambda)
                           
                           class(retval) <- "qvalue"
                           
                           qvalues[ppos, i] = retval$qvalues 
                           qgood = c(qgood, i)
                         }
                         ###########################################################################################################################
                         
                       }
                       qvalues.dataframe = format(qvalues, scientific = FALSE)
                       qvalues.dataframe = as.data.frame(qvalues.dataframe)
                       
                       
                       #######
                       ## 8 ##
                       #######
                       ##Create a csv file that can be uploaded into Pathway Projector and map pathways and label metabolites with 
                       ##p- or q-values. 
                       ########################################################################################################
                       
                       
                       if ( any( colnames(foldchange.all.foldch.dataframe) == paste(treatment, "-", control, sep="") )  )
                       {
                         pfcol = which( colnames(foldchange.all.foldch.dataframe) == paste(treatment, "-", control, sep="") )
                         p.fch = cbind( as.numeric(as.vector(foldchange.all.foldch.dataframe[,pfcol])), 
                                        as.numeric(as.vector(tukey.all.pval.dataframe[,pfcol])) )
                         q.fch = cbind( as.numeric(as.vector(foldchange.all.foldch.dataframe[,pfcol])), 
                                        as.numeric(as.vector(qvalues.dataframe[,pfcol])) )
                       }
                       
                       if ( any( colnames(foldchange.all.foldch.dataframe) == paste(control, "-", treatment, sep="") )  )
                       {
                         pfcol = which( colnames(foldchange.all.foldch.dataframe) == paste(control, "-", treatment, sep="") )
                         p.fch = cbind( as.numeric(as.vector(foldchange.all.foldch.dataframe[,pfcol])), 
                                        as.numeric(as.vector(tukey.all.pval.dataframe[,pfcol])) )
                         q.fch = cbind( as.numeric(as.vector(foldchange.all.foldch.dataframe[,pfcol])), 
                                        as.numeric(as.vector(qvalues.dataframe[,pfcol])) )
                       }
                       
                       row.names(p.fch) = row.names(foldchange.all.foldch.dataframe)
                       colnames(p.fch) = c("Fold-change", "P-value")
                       row.names(q.fch) = row.names(foldchange.all.foldch.dataframe)
                       colnames(q.fch) = c("Fold-change", "Q-value")
                       
                       p.fch = format(p.fch, scientific = FALSE) 
                       q.fch = format(q.fch, scientific = FALSE) 
                       
                       fcolor = c()
                       fsig = c()
                       
                       if (criterionPP == "Pval" || criterion == "Qval")
                       {##1
                         if (criterionPP == "Pval")
                         {
                           PPdata = p.fch
                         }
                         
                         if (criterionPP == "Qval")
                         {
                           PPdata = q.fch
                         }
                         
                         PPdata = as.data.frame(PPdata)
                         
                         if ( any( colnames(foldchange.all.foldch.dataframe) == paste(control, "-", treatment, sep="") )  )
                         { PPdata[,1] = 1/as.numeric(as.vector(PPdata[,1])) }
                         
                         for (i in 1:nrow(PPdata))
                         {
                           if (!is.na(as.numeric(as.matrix(PPdata)[i,1])))
                           {
                             if (as.numeric(PPdata[i,1]) >= as.numeric(1.0))
                             {
                               if (as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t1pq))
                               {
                                 fsig=c(fsig, i)
                                 fcolor[i]="#800000"
                                 if (as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t2pq))
                                 {
                                   fcolor[i]="#FF0000"
                                   if (as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t3pq))	
                                   {
                                     fcolor[i]="#FF00FF"
                                   }				
                                 }
                               }
                             }
                             
                             if (as.numeric(PPdata[i,1]) < as.numeric(1.0))
                             {
                               if (as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t1pq))
                               {
                                 fsig=c(fsig, i)
                                 fcolor[i]="#000080"
                                 if (as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t2pq))
                                 {
                                   fcolor[i]="#0000FF"
                                   if (as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t3pq))	
                                   {
                                     fcolor[i]="#00FFFF"
                                   }				
                                 }
                               }
                             }
                           }
                         }
                       }
                       
                       if (criterionPP == "Fold-change")
                       {
                         if (PPchoice == "Pval")
                         {
                           PPdata = p.fch
                         }
                         if (PPchoice == "Qval")
                         {
                           PPdata = q.fch
                         }
                         
                         PPdata = as.data.frame(PPdata)
                         
                         if ( any( colnames(foldchange.all.foldch.dataframe) == paste(control, "-", treatment, sep="") )  )
                         { PPdata[,1] = 1/as.numeric(as.vector(PPdata[,1])) }
                         
                         for (i in 1:nrow(PPdata))
                         {
                           if (!is.na(as.numeric(as.matrix(PPdata)[i,1])))
                           {
                             if (as.numeric(PPdata[i,1]) >= as.numeric(1.0))
                             {
                               if (as.numeric(PPdata[i,1]) > as.numeric(t1f))
                               {
                                 fsig=c(fsig, i)
                                 fcolor[i]="#800000"
                                 if (as.numeric(PPdata[i,1]) > as.numeric(t2f))
                                 {
                                   fcolor[i]="#FF0000"
                                   if (as.numeric(PPdata[i,1]) > as.numeric(t3f))	
                                   {
                                     fcolor[i]="#FF00FF"
                                   }				
                                 }
                               }
                             }
                             
                             if (as.numeric(PPdata[i,1]) < as.numeric(1.0))
                             {
                               if (as.numeric(PPdata[i,1]) < 1/as.numeric(t1f))
                               {
                                 fsig=c(fsig, i)
                                 fcolor[i]="#000080"
                                 if (as.numeric(PPdata[i,1]) < 1/as.numeric(t2f))
                                 {
                                   fcolor[i]="#FF0000"
                                   if (as.numeric(PPdata[i,1]) < 1/as.numeric(t3f))	
                                   {
                                     fcolor[i]="#00FFFF"
                                   }				
                                 }
                               }
                             }
                           }
                         }
                       }
                       
                       if (criterionPP == "Pval+Fold-change" || criterionPP == "Qval+Fold-change")
                       {
                         if (PPchoice == "Pval")
                         {
                           PPdata = p.fch
                         }
                         if (PPchoice == "Qval")
                         {
                           PPdata = q.fch
                         }
                         
                         PPdata = as.data.frame(PPdata)
                         
                         if ( any( colnames(foldchange.all.foldch.dataframe) == paste(control, "-", treatment, sep="") )  )
                         { PPdata[,1] = 1/as.numeric(as.vector(PPdata[,1])) }
                         
                         for (i in 1:nrow(PPdata))
                         {
                           if (!is.na(as.numeric(as.matrix(PPdata)[i,1])))
                           {
                             if (as.numeric(PPdata[i,1]) >= as.numeric(1.0))
                             {
                               if (as.numeric(PPdata[i,1]) > as.numeric(t1f) && as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t1pq))
                               {
                                 fsig=c(fsig, i)
                                 fcolor[i]="#800000"
                                 if (as.numeric(PPdata[i,1]) > as.numeric(t2f) && as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t2pq))
                                 {
                                   fcolor[i]="#FF0000"
                                   if (as.numeric(PPdata[i,1]) > as.numeric(t3f) && as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t3pq))	
                                   {
                                     fcolor[i]="#FF00FF"
                                   }				
                                 }
                               }
                             }
                             
                             if (as.numeric(PPdata[i,1]) < as.numeric(1.0))
                             {
                               if (as.numeric(PPdata[i,1]) < 1/as.numeric(t1f) && as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t1pq))
                               {
                                 fsig=c(fsig, i)
                                 fcolor[i]="#000080"
                                 if (as.numeric(PPdata[i,1]) < 1/as.numeric(t2f) && as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t2pq))
                                 {
                                   fcolor[i]="#FF0000"
                                   if (as.numeric(PPdata[i,1]) < 1/as.numeric(t3f) && as.numeric(as.matrix(PPdata)[i,2]) < as.numeric(t3pq))	
                                   {
                                     fcolor[i]="#00FFFF"
                                   }				
                                 }
                               }
                             }
                           }
                         }
                       }
                       
                       pq1 = as.numeric(as.vector(PPdata[,2]))
                       pq1 = formatC(pq1, digits = 2, format = "fg")
                       
                       pqvector=c()
                       if (PPchoice == "Pval")
                       { 
                         for ( i in which(!is.na(fcolor))) 
                         {pqvector=c(pqvector, paste("p=", pq1[i], sep="") ) }
                       }
                       
                       if (PPchoice == "Qval")
                       {
                         for ( i in which(!is.na(fcolor))) 
                         {pqvector=c(pqvector, paste("q=", pq1[i], sep="") ) }
                       }
                       
                       if (mode(fsig) != "NULL")
                       {
                         PPmatrix = cbind(fcolor[which(!is.na(fcolor))], dot.size, pqvector, fontsizePP)
                         PPdataframe = as.data.frame(PPmatrix)
                         row.names(PPdataframe) = row.names(PPdata)[fsig]
                         
                         file.name=paste(output, "_colors_", treatment, "-", control, "_", time, ".csv", sep="")
                         write.table(PPdataframe, file=file.name, col.names = FALSE, sep = ",")
                       }
                       
                       
                       
                       for (i in 1:length(aov.good))
                       {
                         if (min(as.numeric(as.vector(tukey.all.pval.dataframe[aov.good[i]-1,])), na.rm=TRUE) < p.cutoff )
                         {
                           tukey.sig.metabolite = cbind(tukey.sig.metabolite, as.numeric(as.vector(data.compile[,aov.good[i]])) )
                           colnames(tukey.sig.metabolite)[ncol(tukey.sig.metabolite)] = colnames(data.compile)[aov.good[i]]
                           
                           tukey.sig.pval = rbind(tukey.sig.pval, as.numeric(as.vector(tukey.all.pval.dataframe[aov.good[i]-1,])))
                           row.names(tukey.sig.pval)[nrow(tukey.sig.pval)] = row.names(tukey.all.pval.dataframe)[aov.good[i]-1]
                           
                           tukey.sig.qval = rbind(tukey.sig.qval, as.numeric(as.vector(as.matrix(qvalues.dataframe[aov.good[i]-1,]))))
                           row.names(tukey.sig.qval)[nrow(tukey.sig.qval)] = row.names(qvalues.dataframe)[aov.good[i]-1]
                           
                           tukey.sig.foldch = rbind(tukey.sig.foldch, as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,])))
                           row.names(tukey.sig.foldch)[nrow(tukey.sig.foldch)] = row.names(foldchange.all.foldch.dataframe)[aov.good[i]-1]
                         }
                         
                         ##I don't know why, but using "as.numeric" in the following line coerced the vector to a factor.  So
                         ##"as.numeric" was omitted.	
                         if (min(as.numeric(as.vector(as.matrix(qvalues.dataframe[(aov.good[i]-1),]))), na.rm=TRUE) < p.cutoff )
                         {
                           qvalue.sig.metabolite = cbind(qvalue.sig.metabolite, as.numeric(as.vector(data.compile[,aov.good[i]])) )
                           colnames(qvalue.sig.metabolite)[ncol(qvalue.sig.metabolite)] = colnames(data.compile)[aov.good[i]]
                           
                           qvalue.sig.pval = rbind(qvalue.sig.pval, as.numeric(as.vector(tukey.all.pval.dataframe[aov.good[i]-1,])))
                           row.names(qvalue.sig.pval)[nrow(qvalue.sig.pval)] = row.names(tukey.all.pval.dataframe)[aov.good[i]-1]
                           
                           qvalue.sig.qval = rbind(qvalue.sig.qval, as.numeric(as.vector(as.matrix(qvalues.dataframe[aov.good[i]-1,]))))
                           row.names(qvalue.sig.qval)[nrow(qvalue.sig.qval)] = row.names(qvalues.dataframe)[aov.good[i]-1]
                           
                           qvalue.sig.foldch = rbind(qvalue.sig.foldch, as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,])))
                           row.names(qvalue.sig.foldch)[nrow(qvalue.sig.foldch)] = row.names(foldchange.all.foldch.dataframe)[i]
                         }
                         
                         
                         if ( max( c( as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,])), 
                                      1/as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,]))), na.rm=TRUE ) >= f.cutoff)
                         {
                           foldchange.sig.metabolite = cbind(foldchange.sig.metabolite, as.numeric(as.vector(data.compile[,aov.good[i]])) )
                           colnames(foldchange.sig.metabolite)[ncol(foldchange.sig.metabolite)] = colnames(data.compile)[aov.good[i]]
                           
                           foldchange.sig.pval = rbind(foldchange.sig.pval, as.numeric(as.vector(tukey.all.pval.dataframe[aov.good[i]-1,])))
                           row.names(foldchange.sig.pval)[nrow(foldchange.sig.pval)] = row.names(tukey.all.pval.dataframe)[aov.good[i]-1]
                           
                           foldchange.sig.qval = rbind(foldchange.sig.qval, as.numeric(as.vector(as.matrix(qvalues.dataframe[aov.good[i]-1,]))))
                           row.names(foldchange.sig.qval)[nrow(foldchange.sig.qval)] = row.names(qvalues.dataframe)[i]
                           
                           foldchange.sig.foldch = rbind(foldchange.sig.foldch, as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,])))
                           row.names(foldchange.sig.foldch)[nrow(foldchange.sig.foldch)] = row.names(foldchange.all.foldch.dataframe)[aov.good[i]-1]
                         }
                         
                         if (min(as.numeric(as.vector(tukey.all.pval.dataframe[aov.good[i]-1,])), na.rm=TRUE) < p.cutoff )
                         {
                           if ( max( c( as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,])), 
                                        1/as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,]))), na.rm=TRUE ) >= f.cutoff)
                           {
                             pf.sig.metabolite = cbind(pf.sig.metabolite, as.numeric(as.vector(data.compile[,aov.good[i]])) )
                             colnames(pf.sig.metabolite)[ncol(pf.sig.metabolite)] = colnames(data.compile)[aov.good[i]]
                             
                             pf.sig.pval = rbind(pf.sig.pval, as.numeric(as.vector(tukey.all.pval.dataframe[aov.good[i]-1,])))
                             row.names(pf.sig.pval)[nrow(pf.sig.pval)] = row.names(tukey.all.pval.dataframe)[aov.good[i]-1]
                             
                             pf.sig.qval = rbind(pf.sig.qval, as.numeric(as.vector(as.matrix(qvalues.dataframe[aov.good[i]-1,]))))
                             row.names(pf.sig.qval)[nrow(pf.sig.qval)] = row.names(qvalues.dataframe)[aov.good[i]-1]
                             
                             pf.sig.foldch = rbind(pf.sig.foldch, as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,])))
                             row.names(pf.sig.foldch)[nrow(pf.sig.foldch)] = row.names(foldchange.all.foldch.dataframe)[aov.good[i]-1]
                           }
                         }
                         
                         if (min(as.numeric(as.vector(as.matrix(qvalues.dataframe[i,]))), na.rm=TRUE) < p.cutoff )
                         {
                           if ( max( c( as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,])), 
                                        1/as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,]))), na.rm=TRUE ) >= f.cutoff)
                           {
                             qf.sig.metabolite = cbind(qf.sig.metabolite, as.numeric(as.vector(data.compile[,aov.good[i]])) )
                             colnames(qf.sig.metabolite)[ncol(qf.sig.metabolite)] = colnames(data.compile)[aov.good[i]]
                             
                             qf.sig.pval = rbind(qf.sig.pval, as.numeric(as.vector(tukey.all.pval.dataframe[aov.good[i]-1,])))
                             row.names(qf.sig.pval)[nrow(qf.sig.pval)] = row.names(tukey.all.pval.dataframe)[aov.good[i]-1]
                             
                             qf.sig.qval = rbind(qf.sig.qval, as.numeric(as.vector(as.matrix(qvalues.dataframe[aov.good[i]-1,]))))
                             row.names(qf.sig.qval)[nrow(qf.sig.qval)] = row.names(qvalues.dataframe)[aov.good[i]-1]
                             
                             qf.sig.foldch = rbind(qf.sig.foldch, as.numeric(as.vector(foldchange.all.foldch.dataframe[aov.good[i]-1,])))
                             row.names(qf.sig.foldch)[nrow(qf.sig.foldch)] = row.names(foldchange.all.foldch.dataframe)[aov.good[i]-1]
                           }
                         }
                       }
                       
                       if (mode(tukey.sig.metabolite) != "NULL") { row.names(tukey.sig.metabolite) = row.names(data.compile) }
                       if (mode(tukey.sig.pval) != "NULL") { colnames(tukey.sig.pval) = colnames(tukey.all.pval.dataframe) }
                       if (mode(tukey.sig.qval) != "NULL") { colnames(tukey.sig.qval) = colnames(qvalues.dataframe) }
                       if (mode(tukey.sig.foldch) != "NULL") { colnames(tukey.sig.foldch) = colnames(foldchange.all.foldch.dataframe) }
                       
                       if (mode(qvalue.sig.metabolite) != "NULL") { row.names(qvalue.sig.metabolite) = row.names(data.compile) }
                       if (mode(qvalue.sig.pval) != "NULL") { colnames(qvalue.sig.pval) = colnames(tukey.all.pval.dataframe) }
                       if (mode(qvalue.sig.qval) != "NULL") { colnames(qvalue.sig.qval) = colnames(qvalues.dataframe) }
                       if (mode(qvalue.sig.foldch) != "NULL") { colnames(qvalue.sig.foldch) = colnames(foldchange.all.foldch.dataframe) }
                       
                       if (mode(foldchange.sig.metabolite) != "NULL") { row.names(foldchange.sig.metabolite) = row.names(data.compile) }
                       if (mode(foldchange.sig.pval) != "NULL") { colnames(foldchange.sig.pval) = colnames(tukey.all.pval.dataframe) }
                       if (mode(foldchange.sig.qval) != "NULL") { colnames(foldchange.sig.qval) = colnames(qvalues.dataframe) }
                       if (mode(foldchange.sig.foldch) != "NULL") { colnames(foldchange.sig.foldch) = colnames(foldchange.all.foldch.dataframe) }
                       
                       if (mode(pf.sig.metabolite) != "NULL") { row.names(pf.sig.metabolite) = row.names(data.compile) }
                       if (mode(pf.sig.pval) != "NULL") { colnames(pf.sig.pval) = colnames(tukey.all.pval.dataframe) }
                       if (mode(pf.sig.qval) != "NULL") { colnames(pf.sig.qval) = colnames(qvalues.dataframe) }
                       if (mode(pf.sig.foldch) != "NULL") { colnames(pf.sig.foldch) = colnames(foldchange.all.foldch.dataframe) }
                       
                       if (mode(qf.sig.metabolite) != "NULL") { row.names(qf.sig.metabolite) = row.names(data.compile) }
                       if (mode(qf.sig.pval) != "NULL") { colnames(qf.sig.pval) = colnames(tukey.all.pval.dataframe) }
                       if (mode(qf.sig.qval) != "NULL") { colnames(qf.sig.qval) = colnames(qvalues.dataframe) }
                       if (mode(qf.sig.foldch) != "NULL") { colnames(qf.sig.foldch) = colnames(foldchange.all.foldch.dataframe) }
                       
                       if (criterion == "Pval") 
                       {
                         mean.plot.data = tukey.sig.metabolite
                         mean.plot.foldch = tukey.sig.foldch
                         mean.plot.foldch = format(mean.plot.foldch, scientific=FALSE)
                         
                         if (plot.choice == "Pval") {mean.plot.pval = tukey.sig.pval} else {mean.plot.pval = tukey.sig.qval}
                         mean.plot.pval = format(mean.plot.pval, scientific=FALSE)
                       }
                       
                       if (criterion == "Qval") 
                       {
                         mean.plot.data = qvalue.sig.metabolite
                         mean.plot.foldch = qvalue.sig.foldch
                         mean.plot.foldch = format(mean.plot.foldch, scientific=FALSE)
                         
                         if (plot.choice == "Pval") {mean.plot.pval = qvalue.sig.pval} else {mean.plot.pval = qvalue.sig.qval}
                         mean.plot.pval = format(mean.plot.pval, scientific=FALSE)
                       }
                       
                       if (criterion == "Fold-change") 
                       {
                         mean.plot.data = foldchange.sig.metabolite
                         mean.plot.foldch = foldchange.sig.foldch
                         mean.plot.foldch = format(mean.plot.foldch, scientific=FALSE)
                         
                         if (plot.choice == "Pval") {mean.plot.pval = foldchange.sig.pval} else {mean.plot.pval = foldchange.sig.qval}
                         mean.plot.pval = format(mean.plot.pval, scientific=FALSE)
                       }
                       
                       if (criterion == "Pval+Fold-change") 
                       {
                         mean.plot.data = pf.sig.metabolite
                         mean.plot.foldch = pf.sig.foldch
                         mean.plot.foldch = format(mean.plot.foldch, scientific=FALSE)
                         
                         if (plot.choice == "Pval") {mean.plot.pval = pf.sig.pval} else {mean.plot.pval = pf.sig.qval}
                         mean.plot.pval = format(mean.plot.pval, scientific=FALSE)
                       }
                       
                       if (criterion == "Qval+Fold-change") 
                       {
                         mean.plot.data = qf.sig.metabolite
                         mean.plot.foldch = qf.sig.foldch
                         mean.plot.foldch = format(mean.plot.foldch, scientific=FALSE)
                         
                         if (plot.choice == "Pval") {mean.plot.pval = qf.sig.pval} else {mean.plot.pval = qf.sig.qval}
                         mean.plot.pval = format(mean.plot.pval, scientific=FALSE)
                       }
                       
                       
                       library(R.oo)
                       
                       if (mode(mean.plot.data) != "NULL")
                       {#1
                         mean.plot.data = cbind(as.character(data.compile$Group), mean.plot.data)
                         colnames(mean.plot.data)[1] = "Group"
                         
                         mean.plot.data = as.data.frame(mean.plot.data)
                         mean.plot.foldch = as.data.frame(mean.plot.foldch)
                         mean.plot.pval = as.data.frame(mean.plot.pval)
                         
                         
                         
                         pdf(file=paste(output, "_meanplots_", time, ".pdf", sep=""), width=width.mean, height=height.mean)
                         par(mfrow=mfrow.mean, mar=mar.mean, las=las.mean, cex=cex.mean, cex.axis=cex.axis.mean, cex.lab=cex.lab.mean, 
                             cex.main=cex.main.mean, cex.sub=cex.sub.mean, mex=mex.mean)
                         
                         for (i in 2:ncol(mean.plot.data))
                         {#2
                           fold.ch.plot=as.numeric(as.vector(as.matrix(mean.plot.foldch[i-1,])))	
                           fold.ch.plot=as.numeric(fold.ch.plot[which(!is.na(fold.ch.plot))])
                           fold.ch.plot=format(fold.ch.plot, digits=3)
                           plot.data=mean.plot.data[,c(1,i)]
                           colnames(plot.data)=c("Group", "Metabolite")
                           
                           ciw.matrix=matrix(nrow=length(levels(as.factor(plot.data$Group))), ncol=3)
                           row.names(ciw.matrix)=levels(as.factor(plot.data$Group))
                           colnames(ciw.matrix)=c("ciw", "ci.max", "ci.min")
                           which.row=1
                           for (j in levels(as.factor(plot.data$Group)))
                           {
                             n = as.numeric(as.vector( plot.data$Metabolite[which(plot.data$Group==j)] ))
                             n = length(which(!is.na(n)))
                             group.stdev = sqrt(var(as.numeric(as.vector(plot.data$Metabolite[which(plot.data$Group==j)])), na.rm=TRUE))
                             group.ciw=qt(0.975, n-1) * group.stdev / sqrt(n)
                             group.mean = mean(as.numeric(as.vector(plot.data$Metabolite[which(plot.data$Group==j)])), na.rm=TRUE)
                             max.ciw = group.mean + group.ciw
                             min.ciw = group.mean - group.ciw
                             ciw.matrix[which.row, c(1,2,3)] = c(group.ciw, max.ciw, min.ciw)
                             which.row=which.row+1
                           }
                           ciw.dataframe=as.data.frame(ciw.matrix)
                           ciw.min.min = min(as.numeric(as.vector(ciw.dataframe$ci.min)),na.rm=TRUE)
                           ciw.max.max = max(as.numeric(as.vector(ciw.dataframe$ci.max)),na.rm=TRUE)
                           
                           y.step = (cex.mean.charsize*ciw.max.max-ciw.min.min)/20 
                           
                           
                           y.spots = 2*y.step + ciw.max.max 
                           for (qq in 2:length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min))))))
                           {
                             pos2add = vector(length=qq)
                             pos2add = replace(pos2add, pos2add=="FALSE", max(as.numeric(as.vector(y.spots)))+y.step )
                             y.spots = c(y.spots, pos2add)
                           }
                           y.spots=sort(y.spots, decreasing=TRUE)
                           
                           y.lim.max = y.step + max(y.spots)
                           y.lim.min = ciw.min.min - y.step
                           
                           test=plotmeans(formula=amodel, data=plot.data, connect=FALSE, main=colnames(mean.plot.data[i]), xlab="", 
                                          ylab="Metabolite", digits=3, ylim =	c(ciw.min.min-4*y.step, y.lim.max), 
                                          mean.labels=FALSE, n.label=FALSE, ci.label=FALSE)
                           
                           plot.data.mean=c()
                           
                           for (xx in 1:length(levels(as.factor(plot.data$Group))))
                           {#3
                             mean.xx = mean(as.numeric(as.vector(plot.data[which(plot.data$Group==levels(as.factor(plot.data$Group))[xx]),2])), na.rm=TRUE)
                             plot.data.mean=c(plot.data.mean, mean.xx  )	
                           }#3end
                           plot.data.mean.digits=format(plot.data.mean, digits=3)
                           text( 1:length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min))))),  ciw.min.min-y.step, 
                                 plot.data.mean.digits[which(!is.na(as.numeric(as.vector(ciw.dataframe$ciw))))], cex=cex.mean.charsize, pos=1)
                           
                           mean.plot.pval = as.matrix(mean.plot.pval)
                           if (min(mean.plot.pval[i-1,], na.rm=TRUE)<p.cutoff)
                           { #4
                             if (ncol(mean.plot.pval)>2)
                             { #5
                               if (max(as.numeric(as.vector(mean.plot.pval[i-1,])), na.rm=TRUE)>=p.cutoff)
                               {#6
                                 pval.row=as.numeric(as.vector(mean.plot.pval[i-1,])) + 0.0000000001
                                 pval.row=replace(pval.row, is.na(pval.row), 100)
                                 pval.row = pval.row[which(as.numeric(as.vector(pval.row))!=100)]
                                 tri.matrix=matrix(0,length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min))))), 
                                                   length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min))))) )
                                 for (ss in 1:length(pval.row))
                                   #for(ss in 1:17)
                                 {	
                                   for (tt in 1:(ncol(tri.matrix)-1))
                                   {
                                     if (max(tri.matrix[,tt], na.rm=TRUE) == 0 ) 
                                       #if (all(tri.matrix[,tt] == "") ) 
                                     {
                                       tri.matrix[tt+1,tt]=pval.row[ss]	
                                       break	
                                     }else{
                                       if ( max(which(tri.matrix[,tt]!=0)) < nrow(tri.matrix) ) 
                                         
                                       {
                                         tri.matrix[(max(which(tri.matrix[,tt]!=0))+1),tt]=pval.row[ss]
                                         break
                                       }
                                       
                                     }
                                   }
                                   
                                   
                                 }
                                 #tri.matrix=replace(tri.matrix, tri.matrix==100, NA)
                                 
                                 pp=apply(tri.matrix, c(1,2), as.numeric)
                                 pp1=pp
                                 dim=nrow(pp)
                                 pp=as.matrix(pp)
                                 alpha=p.cutoff
                                 group = matrix(data=0,nrow=dim, ncol=dim*(dim-1)/2) 
                                 members=matrix(data=0,nrow=dim,ncol=1) 
                                 gcode=1;  ngroup=1 ;
                                 ii=1; 
                                 while (ii>0) 
                                 {#7
                                   kk=0
                                   flag=0; 
                                   for (jj in (ii+1) : dim)
                                   {   #8 # go down row, find group members  
                                     if (pp[jj,ii] > alpha)  
                                     {    #9 # jj and ii are the same  
                                       # check jj against members  
                                       mm=1; 
                                       while(mm <= kk  )
                                       {
                                         ll=members[mm,1] 
                                         if (jj>ll) { test1=pp[jj,ll] } else  {  test1=pp[ll,jj] }
                                         if (test1<0) { test1=-test1 }
                                         if(test1 < alpha) {break}  # need new group  
                                         mm=mm+1
                                       }
                                       if (mm== (kk+1)) 
                                       {#10
                                         mm=ii+1; 
                                         while(mm <= dim) 
                                         { #11
                                           if (mm==jj) { mm=mm+1}  #skip jj (on diagonal) 
                                           if (mm>dim) { break }
                                           if (jj>mm) { test1=pp[jj,mm] } else { test1=pp[mm,jj] }
                                           if ( (test1 > alpha) && (-pp[mm,ii] > alpha) ) 
                                           { 
                                             # previous grouped mean mm may belong in this group  
                                             # so check if already in and current members 
                                             # don't conflict  
                                             ll=1; 
                                             while (ll <= kk )
                                             {
                                               nn=members[ll,1] 
                                               if (nn==mm | nn==0) { break }
                                               if (nn<mm) { test1=pp[mm,nn] } else  { test1=pp[nn,mm] }
                                               if(test1<0.0) { test1=-test1 }
                                               if(test1<alpha) { break }
                                               ll=ll+1
                                             } 
                                             if(ll== (kk+1) )
                                             {	
                                               group[mm,ngroup]=gcode 
                                               kk=kk+1;  members[ll,1]=mm 
                                             } 
                                           } 
                                           mm=mm+1
                                         } #11end
                                         pp[jj,ii]=-pp[jj,ii]   # set so not put in next group  
                                         if(kk>0)
                                         {
                                           for (mm in 1 : kk )
                                           {
                                             ll=members[mm,1] 
                                             # set so not used again  
                                             if (pp[jj,ll]>0) { pp[jj,ll]=-pp[jj,ll] }
                                             if (pp[ll,jj]>0) { pp[ll,jj]=-pp[ll,jj] } 
                                           }
                                         } 
                                         group[jj,ngroup]=gcode 
                                         kk=kk+1;   members[kk,1]=jj 
                                       }else {flag=1 } #10end
                                     } #9end
                                   } #8 end
                                   if(kk==0) 
                                   {   # no members  
                                     for (jj in 1 : ngroup){ if (group[ii,jj] != 0) {break}   }
                                     # not in a group yet, so set flag  
                                     if(jj== (ngroup+1)) {   kk=kk+1 }
                                   } 
                                   if(kk!=0) 
                                   {    # need to set current mean  
                                     group[ii,ngroup]=gcode 
                                     ngroup=ngroup+1;  gcode=gcode+1 
                                   } 
                                   
                                   if(flag!=0) { ii=ii-1}  #  need another group for this mean 
                                   ii=ii+1
                                   if(ii==dim) {break}
                                 } #7end
                                 
                                 ngroup=ngroup-1 
                                 group=group[,1:ngroup] 
                                 letters.1=matrix(data="A",nrow=dim, ncol=ngroup) 
                                 avalue=charToInt('A')
                                 group=as.matrix(group)
                                 for (ii in 1 : dim) for (jj in 1:ngroup)
                                 {
                                   if (group[ii,jj]==0) letters.1[ii,jj]=""
                                   else letters.1[ii,jj]= intToChar(avalue-1+group[ii,jj])
                                 }
                                 letters.plot=c()
                                 for (vv in 1:nrow(letters.1))
                                 {
                                   letters.vector=c()
                                   for (ww in 1:ncol(letters.1))
                                   {
                                     letters.vector=paste(letters.vector,letters.1[vv,ww], sep="")
                                   }
                                   if (letters.vector=="") letters.vector= intToChar(avalue+ngroup)
                                   letters.plot=c(letters.plot,letters.vector)
                                 }
                               }else{#6end
                                 letters.plot=LETTERS[1:length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min)))))]
                               }	
                             }#5end
                             if (ncol(mean.plot.pval)==1) {letters.plot=c("A","B") }
                             
                           }#4end
                           
                           if (min(as.numeric(as.vector(mean.plot.pval[i-1,])), na.rm=TRUE)>=p.cutoff)
                           {
                             letters.plot=vector(length=length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min))))) )
                             letters.plot[1:length(letters.plot)]=""
                           }
                           ##		x.start=1
                           x.spots=c()
                           for (nn in 1:length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min))))))
                           {
                             x.spots=c(x.spots,nn:length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min))))))
                             ##			x.start=x.start+1
                           }
                           ##		if (length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min)))))>2)
                           ##		{
                           ##			y.spots=factor(10.25)
                           ##			for (qq in 1:(length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min)))))-2))
                           ##			{
                           ##				y.spots=c(as.numeric(as.vector(y.spots)), as.numeric(as.vector(levels(as.factor(y.spots)))), 
                           ##					max(as.numeric(as.vector(y.spots)))+cex.mean.charsize/40)  ##orig = -1
                           ##				y.spots=factor(y.spots)
                           ##			}
                           ##			y.spots=as.numeric(as.vector(y.spots))
                           ##			y.spots=sort(y.spots, decreasing=TRUE)/10*y.lim.max
                           ##			mean.spots=vector(length=length(which(!is.na(as.numeric(as.vector(ciw.dataframe$ci.min))))))
                           ##			mean.spots=replace(mean.spots, mean.spots=="FALSE", max(y.spots)+cex.mean.charsize/40)	##orig = y.lim.max
                           ##			y.spots=c(mean.spots, y.spots)
                           ##		}else{y.spots=c(0.9*y.lim.max, 0.9*y.lim.max, 0.8*y.lim.max) }
                           
                           text(x.spots, y.spots, labels=c(letters.plot,fold.ch.plot),cex=cex.mean.charsize )	##orig = 5
                           
                         }#2end
                         dev.off()
                       }#1end
                       
                       mean.plot.pval = as.data.frame(mean.plot.pval)
                       
                       #options(stringsAsFactors = FALSE)
                       
                       #######
                       ## 9 ##
                       #######
                       ##Compare pre- vs. post-normalization data.  Scatterplots are made plotting X (non-interesting variable) 
                       ##vs. Y (Metabolite across groups).  The purpose is to show that normalization removes the relationship
                       ##between the non-interesting variable and metabolite level.
                       ########################################################################################################
                       if (as.numeric(col.start) >= 5)
                       {
                         pdf(file=paste(output, "_pre.v.post_", time, ".pdf", sep=""), width=width.prevpost, height=height.prevpost)
                         par(mfrow=mfrow.prevpost, #mar=mar.prevpost, 
                             las=las.prevpost, cex=cex.prevpost, cex.axis=cex.axis.prevpost, 
                             cex.lab=cex.lab.prevpost, 	cex.main=cex.main.prevpost, cex.sub=cex.sub.prevpost, mex=mex.prevpost)
                         
                         options(scipen=-3)
                         if (mode(data.neg)=="list")	
                         {
                           data.postnorm.neg=cbind( data.neg[,1:(col.start-2)], data.adj.neg )
                           data.prenorm.neg=cbind(data.neg[,1:(col.start-2)],data.neg[neg.good])
                           
                           color.vector=c("red", "darkblue", "black", "chartreuse4", "darkorange3", "darkmagenta", 
                                          "deeppink", "coral4")
                           
                           for (i in 3:(col.start-2))
                           {
                             for (j in (col.start-1):ncol(data.prenorm.neg))
                             {
                               plot.frame.pre.neg = cbind(data.prenorm.neg[,c(2,i,j)], NA, NA )
                               colnames(plot.frame.pre.neg)[4:5]=c("Color", "Style")
                               rainbow.vector=rainbow(length(levels(as.factor(plot.frame.pre.neg$Group))))
                               for (ll in 1:length ( levels(as.factor(plot.frame.pre.neg$Group))) )
                               {
                                 plot.frame.pre.neg$Color[which(as.character(plot.frame.pre.neg$Group) ==  levels(as.factor(plot.frame.pre.neg$Group))[ll]) ] = rainbow.vector[ll]
                                 plot.frame.pre.neg$Style[which(as.character(plot.frame.pre.neg$Group) ==  levels(as.factor(plot.frame.pre.neg$Group))[ll]) ] = ll
                                 
                               }
                               
                               plot.frame.pre.neg = plot.frame.pre.neg[which(!is.na(plot.frame.pre.neg[,3])), ]
                               y.max.pre=max(as.numeric(as.vector(plot.frame.pre.neg[,3])), na.rm=TRUE)
                               y.min.pre=min(as.numeric(as.vector(plot.frame.pre.neg[,3])), na.rm=TRUE)
                               y.lim.max.pre=y.max.pre + 0.1*(y.max.pre - y.min.pre)
                               y.plus.pre=y.max.pre + 0.05*(y.max.pre - y.min.pre)
                               
                               plot.frame.post.neg = cbind(data.postnorm.neg[,c(2,i,j)], NA, NA )
                               colnames(plot.frame.post.neg)[4:5]=c("Color", "Style")
                               rainbow.vector=rainbow(length(levels(as.factor(plot.frame.post.neg$Group))))
                               for (ll in 1:length ( levels(as.factor(plot.frame.post.neg$Group))) )
                               {
                                 plot.frame.post.neg$Color[which(as.character(plot.frame.post.neg$Group) ==  levels(as.factor(plot.frame.post.neg$Group))[ll]) ] = rainbow.vector[ll]
                                 plot.frame.post.neg$Style[which(as.character(plot.frame.post.neg$Group) ==  levels(as.factor(plot.frame.post.neg$Group))[ll]) ] = ll
                                 
                               }
                               
                               plot.frame.post.neg = plot.frame.post.neg[which(!is.na(plot.frame.post.neg[,3])), ]
                               y.max.post=max(as.numeric(as.vector(plot.frame.post.neg[,3])), na.rm=TRUE)
                               y.min.post=min(as.numeric(as.vector(plot.frame.post.neg[,3])), na.rm=TRUE)
                               y.lim.max.post=y.max.post + 0.1*(y.max.post - y.min.post)
                               y.plus.post=y.max.post + 0.05*(y.max.post - y.min.post)
                               
                               plot.c = c()
                               plot.s=c()	
                               for (cc in levels(as.factor(as.vector(plot.frame.pre.neg$Group))))
                               {
                                 plot.c = c( plot.c, levels(as.factor(plot.frame.pre.neg$Color[which(plot.frame.pre.neg$Group==cc)])) )					
                                 plot.s = c( plot.s, levels(as.factor(plot.frame.pre.neg$Style[which(plot.frame.pre.neg$Group==cc)])) )					
                                 
                               }		 
                               
                               layout(matrix(c(1,2,3), nrow=1), widths=c(1/3, 1/3, 1/3))
                               
                               if (length(levels(plot.frame.pre.neg[,2]))==0)
                               {
                                 x.lim.min=min(as.numeric(as.vector(plot.frame.pre.neg[,2])), na.rm=TRUE)
                                 x.lim.max=max(as.numeric(as.vector(plot.frame.pre.neg[,2])), na.rm=TRUE)
                                 x.lim.min.postnorm=min(as.numeric(as.vector(plot.frame.post.neg[,2])), na.rm=TRUE)
                                 x.lim.max.postnorm=max(as.numeric(as.vector(plot.frame.post.neg[,2])), na.rm=TRUE)
                                 
                                 
                                 plot(as.numeric(as.vector(plot.frame.pre.neg[,2])), 
                                      as.numeric(as.vector(plot.frame.pre.neg[,3])),
                                      xlim=c(x.lim.min-0.05*(x.lim.max-x.lim.min), x.lim.max+0.05*(x.lim.max-x.lim.min)), 
                                      ylim=c(y.min.pre - 0.05*(as.numeric(y.plus.pre) - as.numeric(y.min.pre)), 
                                             as.numeric(y.plus.pre)), 
                                      xlab="", 
                                      ylab=colnames(plot.frame.pre.neg)[3], 
                                      cex.lab=cex.lab.prevpost,
                                      cex.axis=cex.axis.prevpost,
                                      cex.main=cex.main.prevpost,
                                      
                                      main=c(colnames(plot.frame.pre.neg)[2], "vs", colnames(plot.frame.pre.neg)[3], "Un-normalized"), 
                                      col=plot.frame.pre.neg$Color, 
                                      pch=plot.frame.pre.neg$Style,
                                      type="p")
                                 
                                 abline(lsfit(as.numeric(as.vector(plot.frame.pre.neg[,2])), 
                                              as.numeric(as.vector(plot.frame.pre.neg[,3])) ))
                                 
                                 plot(as.numeric(as.vector(plot.frame.post.neg[,2])), 
                                      as.numeric(as.vector(plot.frame.post.neg[,3])), 
                                      xlim=c(x.lim.min.postnorm-0.05*(x.lim.max.postnorm-x.lim.min.postnorm), 
                                             x.lim.max.postnorm+0.05*(x.lim.max.postnorm-x.lim.min.postnorm)), 
                                      ylim=c(y.min.post - 0.05*(as.numeric(y.plus.post) - as.numeric(y.min.post)), 
                                             as.numeric(y.plus.post)), 
                                      xlab="", 
                                      ylab=colnames(plot.frame.post.neg)[3], 
                                      cex.lab=cex.lab.prevpost,
                                      cex.axis=cex.axis.prevpost,
                                      cex.main=cex.main.prevpost,
                                      
                                      main=c(colnames(plot.frame.post.neg)[2], "vs", colnames(plot.frame.post.neg)[3], "Normalized"), 
                                      col=plot.frame.post.neg$Color, 
                                      pch=plot.frame.post.neg$Style,
                                      type="p")
                                 
                                 abline(lsfit(as.numeric(as.vector(plot.frame.post.neg[,2])), 
                                              as.numeric(as.vector(plot.frame.post.neg[,3])) ))
                                 
                                 plot(as.numeric(as.vector(plot.frame.post.neg[,2])), 
                                      as.numeric(as.vector(plot.frame.post.neg[,3])), 
                                      xlim=c(x.lim.min.postnorm-0.05*(x.lim.max.postnorm-x.lim.min.postnorm), 
                                             x.lim.max.postnorm+0.05*(x.lim.max.postnorm-x.lim.min.postnorm)), 
                                      ylim=c(y.min.post - 0.05*(as.numeric(y.plus.post) - as.numeric(y.min.post)), 
                                             as.numeric(y.plus.post)), 
                                      type="n", axes=FALSE, ann=FALSE, xpd=NA)
                                 
                                 test=legend( "topleft", 
                                              legend=levels(as.factor(as.vector(plot.frame.pre.neg$Group))),   
                                              plot=TRUE,
                                              text.col=plot.c, pch=as.numeric(plot.s), col=plot.c, xpd=NA )
                                 
                               }else
                               {
                                 plot(plot.frame.pre.neg[,2], as.numeric(as.vector(plot.frame.pre.neg[,3])), 
                                      ylim=c(y.min.pre - 0.05*(as.numeric(y.plus.pre) - as.numeric(y.min.pre)), 
                                             as.numeric(y.plus.pre)), 						
                                      xlab="", 
                                      ylab=colnames(plot.frame.pre.neg[3]), 
                                      cex.lab=cex.lab.prevpost,
                                      cex.axis=cex.axis.prevpost,
                                      cex.main=cex.main.prevpost,
                                      
                                      main=c(colnames(plot.frame.pre.neg[2]), "vs", colnames(plot.frame.pre.neg[3]), 
                                             "Un-normalized"))
                                 
                                 points(plot.frame.pre.neg[,2], as.numeric(as.vector(plot.frame.pre.neg[,3])), 
                                        pch=plot.frame.pre.neg$Style, col=plot.frame.pre.neg$Color)
                                 
                                 plot(plot.frame.post.neg[,2], as.numeric(as.vector(plot.frame.post.neg[,3])), 
                                      ylim=c(y.min.post - 0.05*(as.numeric(y.plus.post) - as.numeric(y.min.post)), 
                                             as.numeric(y.plus.post)), 	
                                      xlab="", 
                                      ylab=colnames(plot.frame.post.neg)[3],
                                      cex.lab=cex.lab.prevpost,
                                      cex.axis=cex.axis.prevpost,
                                      cex.main=cex.main.prevpost,
                                      
                                      
                                      main=c(colnames(plot.frame.post.neg)[2], "vs", colnames(plot.frame.post.neg)[3], 
                                             "Normalized"))
                                 
                                 points(plot.frame.post.neg[,2], as.numeric(as.vector(plot.frame.post.neg[,3])), 
                                        pch=plot.frame.post.neg$Style, col=plot.frame.post.neg$Color)
                                 
                                 plot(as.numeric(as.vector(1:length(plot.frame.post.neg[,3]))), as.numeric(as.vector(plot.frame.post.neg[,3])), 
                                      ylim=c(y.min.post - 0.05*(as.numeric(y.plus.post) - as.numeric(y.min.post)), 
                                             as.numeric(y.plus.post)), 	
                                      type="n", axes=FALSE, ann=FALSE, xpd=NA)
                                 
                                 test=legend( "topleft", 
                                              legend=levels(as.factor(as.vector(plot.frame.pre.neg$Group))), 
                                              plot=TRUE,
                                              text.col=plot.c, 
                                              pch=as.numeric(plot.s), 
                                              col=plot.c, 
                                              xpd=NA )
                                 
                               }
                             }
                           }
                         }
                         
                         data.postnorm.pos=cbind( data.pos[,1:(col.start-2)], data.adj.pos )
                         data.prenorm.pos=cbind(data.pos[,1:(col.start-2)],data.pos[pos.good])
                         
                         color.vector=c("red", "darkblue", "black", "chartreuse4", "darkorange3", "darkmagenta", 
                                        "deeppink", "coral4")
                         
                         for (i in 3:(col.start-2))
                         {
                           for (j in (col.start-1):ncol(data.prenorm.pos))
                           {
                             plot.frame.pre.pos = cbind(data.prenorm.pos[,c(2,i,j)], NA, NA )
                             colnames(plot.frame.pre.pos)[4:5]=c("Color", "Style")
                             rainbow.vector=rainbow(length(levels(as.factor(plot.frame.pre.pos$Group))))
                             for (ll in 1:length ( levels(as.factor(plot.frame.pre.pos$Group))) )
                             {
                               plot.frame.pre.pos$Color[which(as.character(plot.frame.pre.pos$Group) ==  levels(as.factor(plot.frame.pre.pos$Group))[ll]) ] = rainbow.vector[ll]
                               plot.frame.pre.pos$Style[which(as.character(plot.frame.pre.pos$Group) ==  levels(as.factor(plot.frame.pre.pos$Group))[ll]) ] = ll
                               
                             }
                             
                             plot.frame.pre.pos = plot.frame.pre.pos[which(!is.na(plot.frame.pre.pos[,3])), ]
                             y.max.pre=max(as.numeric(as.vector(plot.frame.pre.pos[,3])), na.rm=TRUE)
                             y.min.pre=min(as.numeric(as.vector(plot.frame.pre.pos[,3])), na.rm=TRUE)
                             y.lim.max.pre=y.max.pre + 0.1*(y.max.pre - y.min.pre)
                             y.plus.pre=y.max.pre + 0.05*(y.max.pre - y.min.pre)
                             
                             plot.frame.post.pos = cbind(data.postnorm.pos[,c(2,i,j)], NA, NA )
                             colnames(plot.frame.post.pos)[4:5]=c("Color", "Style")
                             rainbow.vector=rainbow(length(levels(as.factor(plot.frame.post.pos$Group))))
                             for (ll in 1:length ( levels(as.factor(plot.frame.post.pos$Group))) )
                             {
                               plot.frame.post.pos$Color[which(as.character(plot.frame.post.pos$Group) ==  levels(as.factor(plot.frame.post.pos$Group))[ll]) ] = rainbow.vector[ll]
                               plot.frame.post.pos$Style[which(as.character(plot.frame.post.pos$Group) ==  levels(as.factor(plot.frame.post.pos$Group))[ll]) ] = ll
                               
                             }
                             
                             plot.frame.post.pos = plot.frame.post.pos[which(!is.na(plot.frame.post.pos[,3])), ]
                             y.max.post=max(as.numeric(as.vector(plot.frame.post.pos[,3])), na.rm=TRUE)
                             y.min.post=min(as.numeric(as.vector(plot.frame.post.pos[,3])), na.rm=TRUE)
                             y.lim.max.post=y.max.post + 0.1*(y.max.post - y.min.post)
                             y.plus.post=y.max.post + 0.05*(y.max.post - y.min.post)
                             
                             plot.c = c()
                             plot.s=c()	
                             for (cc in levels(as.factor(as.vector(plot.frame.pre.pos$Group))))
                             {
                               plot.c = c( plot.c, levels(as.factor(plot.frame.pre.pos$Color[which(plot.frame.pre.pos$Group==cc)])) )					
                               plot.s = c( plot.s, levels(as.factor(plot.frame.pre.pos$Style[which(plot.frame.pre.pos$Group==cc)])) )					
                               
                             }		 
                             
                             layout(matrix(c(1,2,3), nrow=1), widths=c(1/3, 1/3, 1/3))
                             
                             if (length(levels(plot.frame.pre.pos[,2]))==0)
                             {
                               x.lim.min=min(as.numeric(as.vector(plot.frame.pre.pos[,2])), na.rm=TRUE)
                               x.lim.max=max(as.numeric(as.vector(plot.frame.pre.pos[,2])), na.rm=TRUE)
                               x.lim.min.postnorm=min(as.numeric(as.vector(plot.frame.post.pos[,2])), na.rm=TRUE)
                               x.lim.max.postnorm=max(as.numeric(as.vector(plot.frame.post.pos[,2])), na.rm=TRUE)
                               
                               
                               plot(as.numeric(as.vector(plot.frame.pre.pos[,2])), 
                                    as.numeric(as.vector(plot.frame.pre.pos[,3])),
                                    xlim=c(x.lim.min-0.05*(x.lim.max-x.lim.min), x.lim.max+0.05*(x.lim.max-x.lim.min)), 
                                    ylim=c(y.min.pre - 0.05*(as.numeric(y.plus.pre) - as.numeric(y.min.pre)), 
                                           as.numeric(y.plus.pre)), 
                                    xlab="", 
                                    ylab=colnames(plot.frame.pre.pos)[3],
                                    cex.lab=cex.lab.prevpost,
                                    cex.axis=cex.axis.prevpost,
                                    cex.main=cex.main.prevpost,
                                    main=c(colnames(plot.frame.pre.pos)[2], "vs", colnames(plot.frame.pre.pos)[3], "Un-normalized"), 
                                    col=plot.frame.pre.pos$Color, 
                                    pch=plot.frame.pre.pos$Style,
                                    type="p")
                               
                               abline(lsfit(as.numeric(as.vector(plot.frame.pre.pos[,2])), 
                                            as.numeric(as.vector(plot.frame.pre.pos[,3])) ))
                               
                               plot(as.numeric(as.vector(plot.frame.post.pos[,2])), 
                                    as.numeric(as.vector(plot.frame.post.pos[,3])), 
                                    xlim=c(x.lim.min.postnorm-0.05*(x.lim.max.postnorm-x.lim.min.postnorm), 
                                           x.lim.max.postnorm+0.05*(x.lim.max.postnorm-x.lim.min.postnorm)), 
                                    ylim=c(y.min.post - 0.05*(as.numeric(y.plus.post) - as.numeric(y.min.post)), 
                                           as.numeric(y.plus.post)), 
                                    xlab="", 
                                    ylab=colnames(plot.frame.post.pos)[3], 
                                    main=c(colnames(plot.frame.post.pos)[2], "vs", colnames(plot.frame.post.pos)[3], "Normalized"), 
                                    col=plot.frame.post.pos$Color, 
                                    pch=plot.frame.post.pos$Style,
                                    type="p")
                               
                               abline(lsfit(as.numeric(as.vector(plot.frame.post.pos[,2])), 
                                            as.numeric(as.vector(plot.frame.post.pos[,3])) ))
                               
                               plot(as.numeric(as.vector(plot.frame.post.pos[,2])), 
                                    as.numeric(as.vector(plot.frame.post.pos[,3])), 
                                    xlim=c(x.lim.min.postnorm-0.05*(x.lim.max.postnorm-x.lim.min.postnorm), 
                                           x.lim.max.postnorm+0.05*(x.lim.max.postnorm-x.lim.min.postnorm)), 
                                    ylim=c(y.min.post - 0.05*(as.numeric(y.plus.post) - as.numeric(y.min.post)), 
                                           as.numeric(y.plus.post)), 
                                    type="n", axes=FALSE, ann=FALSE, xpd=NA)
                               
                               test=legend( "topleft", 
                                            legend=levels(as.factor(as.vector(plot.frame.pre.pos$Group))),   
                                            plot=TRUE,
                                            text.col=plot.c, pch=as.numeric(plot.s), col=plot.c, xpd=NA )
                               
                             }else
                             {
                               plot(plot.frame.pre.pos[,2], as.numeric(as.vector(plot.frame.pre.pos[,3])), 
                                    ylim=c(y.min.pre - 0.05*(as.numeric(y.plus.pre) - as.numeric(y.min.pre)), 
                                           as.numeric(y.plus.pre)), 						
                                    xlab="", 
                                    ylab=colnames(plot.frame.pre.pos[3]), 
                                    cex.lab=cex.lab.prevpost,
                                    cex.axis=cex.axis.prevpost,
                                    cex.main=cex.main.prevpost,
                                    main=c(colnames(plot.frame.pre.pos[2]), "vs", colnames(plot.frame.pre.pos[3]), 
                                           "Un-normalized"))
                               
                               points(plot.frame.pre.pos[,2], as.numeric(as.vector(plot.frame.pre.pos[,3])), 
                                      pch=plot.frame.pre.pos$Style, col=plot.frame.pre.pos$Color)
                               
                               plot(plot.frame.post.pos[,2], as.numeric(as.vector(plot.frame.post.pos[,3])), 
                                    ylim=c(y.min.post - 0.05*(as.numeric(y.plus.post) - as.numeric(y.min.post)), 
                                           as.numeric(y.plus.post)), 	
                                    xlab="", 
                                    ylab=colnames(plot.frame.post.pos)[3], 
                                    cex.lab=cex.lab.prevpost,
                                    cex.axis=cex.axis.prevpost,
                                    cex.main=cex.main.prevpost,
                                    
                                    main=c(colnames(plot.frame.post.pos)[2], "vs", colnames(plot.frame.post.pos)[3], 
                                           "Normalized"))
                               
                               points(plot.frame.post.pos[,2], as.numeric(as.vector(plot.frame.post.pos[,3])), 
                                      pch=plot.frame.post.pos$Style, col=plot.frame.post.pos$Color)
                               
                               plot(as.numeric(as.vector(1:length(plot.frame.post.pos[,3]))), as.numeric(as.vector(plot.frame.post.pos[,3])), 
                                    ylim=c(y.min.post - 0.05*(as.numeric(y.plus.post) - as.numeric(y.min.post)), 
                                           as.numeric(y.plus.post)), 	
                                    type="n", axes=FALSE, ann=FALSE, xpd=NA)
                               
                               test=legend( "topleft", 
                                            legend=levels(as.factor(as.vector(plot.frame.pre.pos$Group))), 
                                            plot=TRUE,
                                            text.col=plot.c, 
                                            pch=as.numeric(plot.s), 
                                            col=plot.c, 
                                            xpd=NA )
                               
                             }
                           }
                         }
                         dev.off()
                       } #end of section 9 pre- vs. post-normalization plotting (if there were variables to normalize to)
                       
                       options(scipen=0)
                       
                       
                       f = foldchange.all.foldch.dataframe
                       p = tukey.all.pval.dataframe
                       qval = qvalues.dataframe
                       levene = levene.dataframe
                       shapiro = shapiro.dataframe
                       
                       write.csv(foldchange.all.foldch.dataframe, file=paste(output, "_fold.changes_", time, ".csv", sep=""))
                       write.csv(tukey.all.pval.dataframe, file=paste(output, "_pvalues_", time, ".csv", sep=""))
                       write.csv(qvalues.dataframe, file=paste(output, "_qvalues_", time, ".csv", sep=""))
                       
                       pdf(paste(output, "_heatmap_", time, ".pdf", sep=""), width=width.heatmap, height=height.heatmap)
                       heatmap.data=data.compile[aov.good]
                       #heatmap.matrix=c()
                       #heatmap.good=c()
                       #for (i in 2:ncol(heatmap.data))
                       #{
                       #	if (any(is.na(as.numeric(as.vector(heatmap.data[,i]))))=="FALSE" )
                       #	{	
                       #		heatmap.good=c(heatmap.good, i)	
                       #		heatmap.matrix=cbind(heatmap.matrix, as.numeric(as.vector(heatmap.data[,i])) - mean(as.numeric(as.vector(heatmap.data[,i])))   )	
                       #	}
                       #	
                       #}
                       
                       #colnames(heatmap.matrix)=colnames(heatmap.data)[heatmap.good]
                       #row.names(heatmap.matrix)=row.names(heatmap.data)
                       
                       heatmap.matrix=c()
                       for (i in 1:ncol(heatmap.data))
                       {
                         #heatmap.matrix = cbind(heatmap.matrix, as.numeric(as.vector(heatmap.data[,i])) )
                         
                         #heatmap.vector=as.numeric(as.vector(heatmap.data[,i])) - 
                         #mean(as.numeric(as.vector(heatmap.data[,i])), na.rm=TRUE) 
                         
                         heatmap.vector=as.numeric(as.vector(heatmap.data[,i]))
                         
                         heatmap.vector=replace(heatmap.vector, which(is.nan(heatmap.vector)), NA) 
                         heatmap.matrix=cbind(heatmap.matrix, heatmap.vector) 	
                         
                       }
                       colnames(heatmap.matrix)=colnames(heatmap.data)
                       row.names(heatmap.matrix)=row.names(heatmap.data)
                       
                       
                       par(las=las.heatmap, cex=cex.heatmap, cex.axis=cex.axis.heatmap, 
                           cex.lab=cex.lab.heatmap, cex.main=cex.main.heatmap, cex.sub=cex.sub.heatmap, mex=mex.heatmap)
                       #heatmap.matrix=t(heatmap.matrix)
                       
                       options(show.error.messages=FALSE)
                       try.heatmap = try (heatmap.2(heatmap.matrix, margins=margins.heatmap, col=colorpanel(20, "blue", "black", "red"), 
                                                    trace="none", density.info="none", scale="column") )
                       options(show.error.messages=TRUE)
                       
                       if (class(try.heatmap)=="try-error")
                       {
                         heatmap.matrix=c()
                         for (i in 1:ncol(heatmap.data))
                         {
                           #heatmap.matrix = cbind(heatmap.matrix, as.numeric(as.vector(heatmap.data[,i])) )
                           
                           #heatmap.vector=as.numeric(as.vector(heatmap.data[,i])) - 
                           #mean(as.numeric(as.vector(heatmap.data[,i])), na.rm=TRUE) 
                           
                           heatmap.vector=as.numeric(as.vector(heatmap.data[,i]))
                           
                           heatmap.vector=replace(heatmap.vector, which(is.nan(heatmap.vector)), 0) 
                           heatmap.matrix=cbind(heatmap.matrix, heatmap.vector) 	
                           
                         }
                         colnames(heatmap.matrix)=colnames(heatmap.data)
                         row.names(heatmap.matrix)=row.names(heatmap.data)
                         
                         
                         par(las=las.heatmap, cex=cex.heatmap, cex.axis=cex.axis.heatmap, 
                             cex.lab=cex.lab.heatmap, cex.main=cex.main.heatmap, cex.sub=cex.sub.heatmap, mex=mex.heatmap)
                         #heatmap.matrix=t(heatmap.matrix)
                         
                         options(show.error.messages=FALSE)
                         try.heatmap = try (heatmap.2(heatmap.matrix, margins=margins.heatmap, col=colorpanel(20, "blue", "black", "red"), 
                                                      trace="none", density.info="none", scale="column") )
                         
                         
                         
                       }
                       
                       dev.off()
                       
                       ########
                       ## 10 ##
                       ########
                       ##Print summary statistics and results from ANOVA and Tukey HSD post-hoc tests of metabolites showing significant between-group 
                       ##mean differences.
                       ###############################################################################################################################
                       
                       for (i in 1:ncol(shapiro)) { shapiro[,i] = as.numeric(as.matrix(shapiro)[,i]) }
                       shapiro = format(shapiro, digits = 3)
                       
                       for (i in 1:ncol(levene)) { levene[,i] = as.numeric(as.matrix(levene)[,i]) }
                       levene = format(levene, digits = 3)
                       
                       for (i in 1:ncol(f)) { f[,i] = as.numeric(as.matrix(f)[,i]) }
                       f = format(f, digits = 3)
                       
                       for (i in 1:ncol(p)) { p[,i] = as.numeric(as.matrix(p)[,i]) }
                       p = format(p, digits = 3)
                       
                       for (i in 1:ncol(qval)) { qval[,i] = as.numeric(as.matrix(qval)[,i]) }
                       qval = format(qval, digits = 3)
                       
                       p <<- as.data.frame(p)
                       f <<- as.data.frame(f)
                       qval <<- as.data.frame(qval)
                       levene <<- as.data.frame(levene)
                       shapiro <<- as.data.frame(shapiro)
                       
                       
                       if (mode(non.normal.dataframe)=="list")
                       {
                         for (i in 1:ncol(non.normal.dataframe))
                         { non.normal.dataframe[,i] = as.numeric(non.normal.dataframe[,i]) }
                         non.normal.dataframe = format(non.normal.dataframe, digits = 3)
                         non.normal.dataframe <<- as.data.frame(non.normal.dataframe)
                         
                         cat("(1) Metabolites with non-normally distributed residuals, based on Shapiro-Wilk
                             test of normality:", "\n");
                         print(non.normal.dataframe); cat("\n"); cat("    Enter 'shapiro' to see Shapiro-Wilk tests for all metabolites.", "\n");
                       }else
                       { 	cat("(1) All metabolites passed Shapiro-Wilk test of normality with p-value of at least 0.05.", "\n"); cat("\n")  
                         cat("    Enter 'shapiro' to see Shapiro-Wilk tests for all metabolites.", "\n") }; cat("\n");
                       
                       
                       if (mode(non.equal.var)=="list") 
                       { 	
                         for (i in 1:ncol(non.normal.dataframe))
                         { non.equal.var[,i] = as.numeric(non.equal.var[,i]) }
                         non.equal.var = format(non.equal.var, digits = 3)
                         non.equal.var <<- as.data.frame(non.equal.var)
                         
                         
                         cat("(2) Metabolites with non-equal variance among groups, based on Levene's test of 
                             equality of varance:", "\n")
                         print(non.equal.var); cat("\n"); cat("    Enter 'levene' to see Levene tests for all metabolites.", "\n"); 
                       }else
                       { 	cat("(2) All metabolites passed Levene's test of equality of variance with p-value of at least 0.05.", "\n"); cat("\n");
                         cat("    Enter 'levene' to see Levene tests for all metabolites.", "\n") }; cat("\n");
                       
                       
                       if (mode(mean.plot.data)=="list")
                       {
                         for (i in 1:ncol(mean.plot.foldch)) { mean.plot.foldch[,i] = as.numeric(as.matrix(mean.plot.foldch)[,i]) }
                         mean.plot.foldch = format(mean.plot.foldch, digits = 3)
                         mean.plot.foldch <<- as.data.frame(mean.plot.foldch)	
                         
                         for (i in 1:ncol(mean.plot.pval)) { mean.plot.pval[,i] = as.numeric(as.matrix(mean.plot.pval)[,i]) }
                         mean.plot.pval = format(mean.plot.pval, digits = 3)
                         mean.plot.pval <<- as.data.frame(mean.plot.pval)
                         
                         cat("(3) Between-group mean fold-changes for metabolites showing significant between-group mean differences, 
                             based on chosen significance threshold:", "\n"); print(mean.plot.foldch); cat("\n")
                         cat("(4) ANOVA Tukey HSD p-values or q-values for metabolites showing significant between-group 
                             mean differences, based on chosen significance threshold:", "\n"); print(mean.plot.pval); cat("\n") 
                         cat("    Enter 'f' to see all between-group mean fold-changes.", "\n")
                         cat("    Enter 'p' to see p-values for all metabolites.", "\n") 
                         cat("    Enter 'qval' to see q-values for all metabolites.", "\n") 
                         
                       }else
                       {cat("(3) No metabolites showed significant between-group mean differences.  No mean-plot pdf file was created.", "\n"); cat("\n");
                         cat("    Enter 'f' to see all between-group mean fold-changes.", "\n")
                         cat("    Enter 'p' to see p-values for all metabolites.", "\n") 
                         cat("    Enter 'qval' to see q-values for all metabolites.", "\n") }
                       
                       ##End of MetabR function									
                       ###############################################################################################################################
                       ###############################################################################################################################
                       
                     }
                     
                         )
  
  
  ##End of d function
  ###############################################################################################################################
  ###############################################################################################################################
  
  
  
  
                       }
)	
