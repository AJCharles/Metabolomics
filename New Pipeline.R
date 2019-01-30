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

##### 2.0 Metabolome Report #####

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

##### 2.1 Raw Data Clustering Analysis #####

# Make a PCA for the non-annotated data (technical replicates) to see assess clustering.
# Technical replicates (e.g. 1,2,3 or 10,11,12) *should* roughly cluster in the PCA. 
# This is the raw data prior to control so be wary of the output. 
# Note that the data is subset from the original diffreport to prevent change. 

## Need to Subselect and transpose the data.
PCA_data_TR <- as.data.frame(t(diffreport[,17:(16+Treps)]))
head(t(PCA_data_TR),4) # Note that this is transposed to ease viewing... 
sum(is.na(PCA_data_TR)) # See if there are NA values. If present, run the next line. 
PCA_data_TR[is.na(PCA_data_TR)] <- 0 # Replace any NA values with 0 for the PCA. 
sum(is.na(PCA_data_TR)) # Check that NA values have been replaced. 
                       # Note that if there were very few, then this isn't an issue, otherwise PCA
                       # May not be best appropriate for zero-inflated data. 

# Shorten the row names for plotting
rownames(PCA_data_TR) <- row.names(FileNames)
head(t(PCA_data_TR),4) # Check that the rename has worked

# Save df for later use
write.table(PCA_data_TR, file="PCA_data_TR.txt", sep="\t", row.names = FALSE)

# Create the PCA on the TR.
library(factoextra)

res.pca.TR <- prcomp(PCA_data_TR, scale = TRUE, center = TRUE) 
fviz_pca_ind(res.pca.TR,
             col.ind= TR_groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07", "green", "purple"), #Change colour codes / names as desired.
             addEllipses = TRUE,  
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
fviz_eig(res.pca.TR) # Generate the scree plot for variance represented by each principal component

## Make a PLS-DA for comparison

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mixOmics", version = "3.8")
library(mixOmics)

## PLS-DA (all features)
PLSDA_data_TR <- as.matrix(PCA_data_TR[1:Treps,1:length(PCA_data_TR)])
str(PLSDA_data_TR) # check data actually has been transferred in!
dim(PLSDA_data_TR) # Check correct dimenstions

Y_Factor <- data.frame(matrix(NA, nrow=Treps, ncol=0))
Y_Factor$Diet_group <- c(rep("2>2", 18), rep("2>8", 18), rep("8>2", 18), rep("8>8", 18))
Y_Factor$Num_group <- c(rep("1", 18), rep("2", 18), rep("3", 18), rep("4", 18))
Y_Factor$Diet_group <- as.factor(Y_Factor$Diet_group)
rownames(Y_Factor) <- row.names(FileNames)
str(Y_Factor) # Check change.

plsda_mod_TR <- plsda(PLSDA_data_TR, # 
                      Y_Factor$Diet_group,# either num/diet, just left for now..  
                      ncomp = 10, # the number of components to include in the model. Default to 2.
                      scale = TRUE, # Boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE)
                      mode = c("regression"),  ##### See stats below (1 of 4 choices)
                      tol = 1e-6, # Convergence stopping value.
                      max.iter = 1000000, # Integer, the maximum number of iterations.
                      near.zero.var = FALSE, # boolean, see the internal nearZeroVar function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE.
                      logratio="none",  # one of "none", "CLR". Specifies the log ratio transformation to deal with compositional values that may arise from specific normalisation in sequencing dadta. Default to 'none'.
                      multilevel=NULL, # Sample information for multilevel decomposition for repeated measurements. A numeric matrix or data frame indicating the repeated measures on each individual, i.e. the individuals ID. See examples in ?splsda.
                      all.outputs = TRUE) # Boolean. Computation can be faster when some specific (and non-essential) outputs are not calculated. Default = TRUE.
plotIndiv(plsda_mod_TR, ind.names = Y_Factor$Diet_group, ellipse = TRUE, abline = TRUE, legend =TRUE)

### mode = The type of algorithm to use is specified with the mode argument. 
### Four PLS algorithms are available: PLS regression ("regression"), PLS canonical analysis 
### ("canonical"), redundancy analysis ("invariant") and the classical PLS algorithm ("classic") 
### (see References). Different modes relate on how the Y matrix is deflated across the iterations 
### of the algorithms - i.e. the different components.
# - Regression mode: the Y matrix is deflated with respect to the information extracted/modelled from 
#   the local regression on X. Here the goal is to predict Y from X (Y and X play an asymmetric role). 
#   Consequently the latent variables computed to predict Y from X are different from those computed to 
#   predict X from Y.
# - Canonical mode: the Y matrix is deflated to the information extracted/modelled from the local 
#   regression on Y. Here X and Y play a symmetric role and the goal is similar to a Canonical 
#   Correlation type of analysis.
# - Invariant mode: the Y matrix is not deflated
# - Classic mode: is similar to a regression mode. It gives identical results for the variates and 
#   loadings associated to the X data set, but differences for the loadings vectors associated to the
#   Y data set (different normalisations are used). Classic mode is the PLS2 model as defined by 
#   Tenenhaus (1998), Chap 9.
#      - From what was mentioned in Recton & Lloyd '14, PLS2 is even more dangerous, easily misused.

# Note that in all cases the results are the same on the first component as deflation only starts after
# component 1.

perf.plsda <- perf(plsda_mod_TR, validation = "Mfold", folds = 5, 
                   progressBar = TRUE, auc = TRUE, nrepeat = 100) 
plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

# We use the function perf to evaluate a PLS-DA model, using 5-fold cross-validation repeated 10 times.
# In addition, the function is useful to decide on the optimal number of components to choose. The 
# Balanced Error Rate (BER) and overall error rate is displayed for all prediction distances.
plsda_var <- as.data.frame(plsda_mod_TR$explained_variance)
plsda_x_var <- as.data.frame((plsda_var$X[1:10])*100)
colnames(plsda_x_var) <- "X"
plsda_x_var$Dimensions <- seq(1, 10, by =1)

library(ggplot2)
# Scree plot equivalent for PLS-DA data.
ggplot(data=plsda_x_var, aes(x=Dimensions, y=X)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_line(data=plsda_x_var, aes(x=Dimensions, y=X), colour="black") +
  geom_point(data=plsda_x_var, aes(x=Dimensions, y=X), colour="black") +
  ylab("Percentage of explained variances") +
  title("Scree plot") +
  ggtitle("PLS-DA X Variances") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal() +
  theme(plot.title = element_text(size=10))


# Ok, if the replicates roughly group then continue to the next part; transformation, averaging and annotating the peaks.
# If they don't, then figure it out. The pipeline is just set up for triplicates for now, so can't remove
# problem queries and average two in a specific group instead of three. Possibly something I'll work on in the 
# future, after I extend the capabilities. 

##### 2.2 Data Transformation / Normalisation #####

BiocManager::install("edgeR", version = "3.8")
library(edgeR)

# Data import and preparation
edgeR_Data <- as.data.frame(t(diffreport[,17:(16+Treps)])) # pull data again
View(row.names(edgeR_Data))
rownames(edgeR_Data) <- row.names(FileNames) # correct naming names
View(row.names(edgeR_Data)) # Check the names have been shortened.
edgeR_Data <- t(edgeR_Data)
View(colnames(edgeR_Data)) # Made the df t() accidentally... just leave like this for now.

# Remove NA values
sum(is.na(edgeR_Data)) # Check if there are any NA values present
edgeR_Data[is.na(edgeR_Data)] <- 0 # Replace any NA values with 0 
sum(is.na(edgeR_Data)) # Removed?

# Remove zero values (intensities usually 100's-1000's, so a value of 1 is negligible)
sum(edgeR_Data == 0) # Check if there are any 0 values present 
edgeR_Data[edgeR_Data == 0] <- 1 # '[]'s not '()' when addressing values! use brackets for functions :) 
sum(edgeR_Data == 0) # Check if there are any 0 values present 

## Plot the raw data. ## - See what it's like, compare before and after.

# Red will always be group 1, group 2 blue, group 3 green, group 4 yellow. (Groups should be set in alphabetical order as shown in the working directory)
# Note: par(mfrow = c(r, c)) can be used to create multipanel plots. Use QuickGraphics.R for ease. 
boxplot(log2(edgeR_Data), las = 2, col = c(rep("red", ((Treps)/Tgroups)), 
                                           rep("blue", ((Treps)/Tgroups)), 
                                           rep("green", ((Treps)/Tgroups)), 
                                           rep("yellow", ((Treps)/Tgroups))),
        xlab = "Sample", ylab = bquote(log[2]~"Intensity"),
        outcex = 0.5, cex.axis = 0.4,
        notch = TRUE, main = "Initial data")

# Create colour grouping information for plotDensitites (in order, 1-72):
# Group order (pre sort): 1 = 2-2, 2 = 2-8, 3 = 8-2, 4 = 8-8 | colour: r,b,g,y *key*
pd_order <- Y_Factor[ order(row.names(Y_Factor)), ] # reorder df accoring to rowname
pd_order <- pd_order[,-1] # Drop unnecessary Diet_groups column
pd_order <- as.data.frame(pd_order) # Convert back to a df

pd_order$pd_order <- as.character(pd_order$pd_order)
pd_order$pd_order[pd_order$pd_order == "1"] <- "red"
pd_order$pd_order[pd_order$pd_order == "2"] <- "blue"
pd_order$pd_order[pd_order$pd_order == "3"] <- "green"
pd_order$pd_order[pd_order$pd_order == "4"] <- "yellow"

plotDensities(log2(edgeR_Data), col = pd_order$pd_order,
              legend = FALSE,
              main = "Initial data")

## Data Transformation & Normalisation ##

format(round(colSums(edgeR_Data), digits = 0), big.mark = ",") # We'll be Normalising these values

# calculate the global scaling values for sample loading normalizations
target_eRD <- mean(colSums(edgeR_Data))

# do the sample loading normalization before other normalizations
# there is a different correction factor for each column
norm_facs_eRD <- target_eRD / colSums(edgeR_Data)
data_eRD_sl <- sweep(edgeR_Data, 2, norm_facs_eRD, FUN = "*")

# Plot after Sample Loading normalisations:
boxplot(log2(data_eRD_sl), las = 2, col = c(rep("red", ((Treps)/Tgroups)), 
                                           rep("blue", ((Treps)/Tgroups)), 
                                           rep("green", ((Treps)/Tgroups)), 
                                           rep("yellow", ((Treps)/Tgroups))),
        xlab = "Sample", ylab = bquote(log[2]~"Intensity"),
        outcex = 0.5, cex.axis = 0.4,
        notch = TRUE, main = "eRD_sl data")

plotDensities(log2(data_eRD_sl), col = pd_order$pd_order,
              legend = FALSE,
              main = "eRD_sl data")

# check the columnn totals (should be more like equal than before)
format(round(colSums(data_eRD_sl), digits = 0), big.mark = ",")

# Normalise for library size differences.
eRD_tmm <- calcNormFactors(data_eRD_sl)
data_eRD_tmm <- sweep(data_eRD_sl, 2, eRD_tmm, FUN = "/")

# Plot after library size (& sample loading) normalisations:
boxplot(log2(data_eRD_tmm), las = 2, col = c(rep("red", ((Treps)/Tgroups)), 
                                            rep("blue", ((Treps)/Tgroups)), 
                                            rep("green", ((Treps)/Tgroups)), 
                                            rep("yellow", ((Treps)/Tgroups))),
        xlab = "Sample", ylab = bquote(log[2]~"Intensity"),
        outcex = 0.5, cex.axis = 0.4,
        notch = TRUE, main = "eRD_tmm data")

plotDensities(log2(data_eRD_tmm), col = pd_order$pd_order,
              legend = FALSE,
              main = "eRD_tmm data")

##### 2.3 Transformed Data Clustering Analysis #####

# Check trandformed data
dim(data_eRD_tmm)
sum(is.na(data_eRD_tmm)) # See if there are NA values. Shouldn't be. 
View(row.names(data_eRD_tmm)) # the 'features' are not named. Could pull across, doesn't matter.
View(colnames(data_eRD_tmm)) # the view transposes the data, but should start with 10, end with 66.

# PCA on transformed data (still TR's)
library(factoextra)

res.pca.T_TR <- prcomp(t(data_eRD_tmm), scale = TRUE, center = TRUE) # just transpose here
fviz_pca_ind(res.pca.T_TR,
             col.ind= TR_groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07", "green", "purple"), #Change colour codes / names as desired.
             addEllipses = TRUE,  
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
fviz_eig(res.pca.T_TR) # Generate the scree plot for variance represented by each principal component

# PLS-DA on transformed data (still TR's)
library(mixOmics)
plsda_mod_T_TR <- plsda(t(data_eRD_tmm), # 
                      Y_Factor$Diet_group,# either num/diet, just left for now..  
                      ncomp = 10, # the number of components to include in the model. Default to 2.
                      scale = TRUE, # Boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE)
                      mode = c("regression"),  ##### See stats below (1 of 4 choices)
                      tol = 1e-6, # Convergence stopping value.
                      max.iter = 1000000, # Integer, the maximum number of iterations.
                      near.zero.var = FALSE, # boolean, see the internal nearZeroVar function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE.
                      logratio="none",  # one of "none", "CLR". Specifies the log ratio transformation to deal with compositional values that may arise from specific normalisation in sequencing dadta. Default to 'none'.
                      multilevel=NULL, # Sample information for multilevel decomposition for repeated measurements. A numeric matrix or data frame indicating the repeated measures on each individual, i.e. the individuals ID. See examples in ?splsda.
                      all.outputs = TRUE) # Boolean. Computation can be faster when some specific (and non-essential) outputs are not calculated. Default = TRUE.
plotIndiv(plsda_mod_T_TR, ind.names = Y_Factor$Diet_group, ellipse = TRUE, abline = TRUE, legend =TRUE)

# Be warned, this takes about 15-20 mins to run. 
perf.plsda_T <- perf(plsda_mod_T_TR, validation = "Mfold", folds = 5, 
                   progressBar = TRUE, auc = TRUE, nrepeat = 100) 
plot(perf.plsda_T, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

# Pull data to make 'scree'-like plot:
plsda_var_T_TR <- as.data.frame(plsda_mod_T_TR$explained_variance)
plsda_x_var_T_TR <- as.data.frame((plsda_var_T_TR$X[1:10])*100)
colnames(plsda_x_var_T_TR) <- "X"
plsda_x_var_T_TR$Dimensions <- seq(1, 10, by =1)

library(ggplot2)
# Scree plot equivalent for PLS-DA data.
ggplot(data=plsda_x_var_T_TR, aes(x=Dimensions, y=X)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_line(data=plsda_x_var_T_TR, aes(x=Dimensions, y=X), colour="black") +
  geom_point(data=plsda_x_var_T_TR, aes(x=Dimensions, y=X), colour="black") +
  ylab("Percentage of explained variances") +
  ggtitle("PLS-DA X Variances") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal() +
  theme(plot.title = element_text(size=10))

##### 3.0 Average TR's into their respective groups #####

## Pull the transformed data back into the diffreport to keep track. 
diffreport[,17:(16+Treps)] <- data_eRD_tmm # Pull transformed data
colnames(diffreport)[17:(16+Treps)] <- row.names(FileNames) # Shorten names
View(colnames(diffreport)) # Check shortening was successful

## Save the diffreport (transformed/normalised)
write.table(diffreport, file="t_diffreport.txt", sep="\t", row.names = FALSE)

## Group the technical reps according to sample. 
AV_T_data <- t(data_eRD_tmm) # Just using the data, no need to pull back from diffreprt
AV_T_data <- rowsum(AV_T_data, rep(1:(Treps/3), each=3))/3 
dim(AV_T_data) # should = Treps / 3
View(row.names(AV_T_data)) # The sample numbering is now incorrect

# Generate the a new list correlating to the correct number formatting after averaging 
# the technical replicates and correct the names
row.names(AV_T_data) <- AVFileNames$AVFileNames
View(row.names(AV_T_data)) # The sample numbering is now correct 


## Create & Save the diffreport (averaged & transformed/normalised)
tmp_AVT <- t(AV_T_data) # Definitely correct...
diffreport[, c(17:(16+Treps))] <- list(NULL) # ()! - Remove old intensity data.
diffreport <- cbind(diffreport,tmp_AVT)
diffreport[,(17:(16+(Treps/3)))] <- AV_T_data
write.table(diffreport, file="av_t_diffreport.txt", sep="\t", row.names = FALSE) # Save AVeraged & Transformed diffreport

##### 3.1 Averaged & Transformed Clustering Analysis #####
# Check Averaged & Transformed data
dim(AV_T_data)
sum(is.na(AV_T_data)) # See if there are NA values. Shouldn't be. 
row.names(AV_T_data) # Actual sample number
colnames(AV_T_data) # 'Feature' number

# PCA on Averaged & Transformed data
library(factoextra)

res.pca.AV_T <- prcomp(AV_T_data, scale = TRUE, center = TRUE)
fviz_pca_ind(res.pca.AV_T,
             col.ind= AV_groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07", "green", "purple"), #Change colour codes / names as desired.
             addEllipses = TRUE,  
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
fviz_eig(res.pca.AV_T) # Generate the scree plot for variance represented by each principal component

# PLS-DA on transformed data (still TR's)
library(mixOmics)
dim(AV_T_data)
plsda_mod_AV_T <- plsda(AV_T_data, # 
                        AV_groups[,1],# either num/diet, just left for now..  
                        ncomp = 10, # the number of components to include in the model. Default to 2.
                        scale = TRUE, # Boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE)
                        mode = c("regression"),  ##### See stats below (1 of 4 choices)
                        tol = 1e-6, # Convergence stopping value.
                        max.iter = 1000000, # Integer, the maximum number of iterations.
                        near.zero.var = FALSE, # boolean, see the internal nearZeroVar function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE.
                        logratio="none",  # one of "none", "CLR". Specifies the log ratio transformation to deal with compositional values that may arise from specific normalisation in sequencing dadta. Default to 'none'.
                        multilevel=NULL, # Sample information for multilevel decomposition for repeated measurements. A numeric matrix or data frame indicating the repeated measures on each individual, i.e. the individuals ID. See examples in ?splsda.
                        all.outputs = TRUE) # Boolean. Computation can be faster when some specific (and non-essential) outputs are not calculated. Default = TRUE.
plotIndiv(plsda_mod_AV_T, ind.names = AV_groups[,1], ellipse = TRUE, abline = TRUE, legend =TRUE)

# Be warned, this takes about 15-20 mins to run. 
perf.plsda_T <- perf(plsda_mod_T_TR, validation = "Mfold", folds = 5, 
                     progressBar = TRUE, auc = TRUE, nrepeat = 100) 
plot(perf.plsda_T, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

# Pull data to make 'scree'-like plot:
plsda_var_AV_T <- as.data.frame(plsda_mod_AV_T$explained_variance)
plsda_x_var_AV_T <- as.data.frame((plsda_var_AV_T$X[1:10])*100)
colnames(plsda_x_var_AV_T) <- "X"
plsda_x_var_AV_T$Dimensions <- seq(1, 10, by =1)

library(ggplot2)
# Scree plot equivalent for PLS-DA data.
ggplot(data=plsda_x_var_AV_T, aes(x=Dimensions, y=X)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_line(data=plsda_x_var_AV_T, aes(x=Dimensions, y=X), colour="black") +
  geom_point(data=plsda_x_var_AV_T, aes(x=Dimensions, y=X), colour="black") +
  ylab("Percentage of explained variances") +
  ggtitle("PLS-DA X Variances") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal() +
  theme(plot.title = element_text(size=10))

##### 3.2 Annotation and Pathway Enrichment Analysis <- not useful overall, but good for separate groups ##### 

setwd("D:/Metabolomics")
diffreport <-read.table("av_t_diffreport.txt", sep="\t", header=TRUE)
dim(diffreport) # == 40 columns (AV intensities, =/= TRs).
colnames(diffreport[,c(6,4,3)]) # == 'm.z' 'p.value' 't.score' (change made in 2.0)
# Annotation
PeakListProfile <- diffreport[,c(6,4,3)]
write.table(PeakListProfile, file="PeakListProfile.txt", sep="\t", row.names = FALSE)

library(MetaboAnalystR)
# Get the annotation results from metaboanalyst - need an active internet connection!

## Kegg
kegg_mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
kegg_mSet <- Read.PeakListData(kegg_mSet, "av_t_diffreport.txt");
kegg_mSet <- UpdateMummichogParameters(kegg_mSet, "6.0", "positive", 0.05); # PPM, ion mode, significance value:
  # Numeric, specify the p-value cutoff to define significant m/z features from reference m/z features
kegg_mSet <- SanityCheckMummichogData(kegg_mSet)
kegg_mSet <- PerformMummichog(kegg_mSet, "dme_kegg", "fisher", "gamma") # organism_database, test, distribution.
kegg_mSet$matches.res$Common.Name <- doKEGG2NameMapping(kegg_mSet$matches.res$Matched.Compound)
names(kegg_mSet$matches.res)[names(kegg_mSet$matches.res) == 'Query.Mass'] <- 'm.z'
names(kegg_mSet$matches.res)[names(kegg_mSet$matches.res) == 'Matched.Compound'] <- 'kegg.Matched.Compund'
names(kegg_mSet$matches.res)[names(kegg_mSet$matches.res) == 'Matched.Form'] <- 'kegg.Matched.Form'
names(kegg_mSet$matches.res)[names(kegg_mSet$matches.res) == 'Mass.Diff'] <- 'kegg.Mass.Diff'
names(kegg_mSet$matches.res)[names(kegg_mSet$matches.res) == 'Common.Name'] <- 'kegg.Common.Name'
write.table(kegg_mSet$matches.res, file="KeggMatchedCompoundTable.txt", sep="\t", row.names = FALSE)
Kegg_MCT<-kegg_mSet$matches.res
Kegg_MCT$m.z <- as.numeric(Kegg_MCT$m.z) # Fix output

## BioCyc
biocyc_mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
biocyc_mSet <- Read.PeakListData(biocyc_mSet, "av_t_diffreport.txt");
biocyc_mSet <- UpdateMummichogParameters(biocyc_mSet, "6.0", "positive", 0.05); # PPM, ion mode, significance value:
  # Numeric, specify the p-value cutoff to define significant m/z features from reference m/z features
biocyc_mSet <- SanityCheckMummichogData(biocyc_mSet)
biocyc_mSet <- PerformMummichog(biocyc_mSet, "dme_biocyc", "fisher", "gamma") # organism_database, test, distribution.
names(biocyc_mSet$matches.res)[names(biocyc_mSet$matches.res) == 'Query.Mass'] <- 'm.z'
names(biocyc_mSet$matches.res)[names(biocyc_mSet$matches.res) == 'Matched.Compound'] <- 'biocyc.Matched.Compund'
names(biocyc_mSet$matches.res)[names(biocyc_mSet$matches.res) == 'Matched.Form'] <- 'biocyc.Matched.Form'
names(biocyc_mSet$matches.res)[names(biocyc_mSet$matches.res) == 'Mass.Diff'] <- 'biocyc.Mass.Diff'
write.table(biocyc_mSet$matches.res, file="BioCycMatchedCompoundTable.txt", sep="\t", row.names = FALSE)
biocyc_MCT<-biocyc_mSet$matches.res
biocyc_MCT$m.z <- as.numeric(biocyc_MCT$m.z) # Fix output

# Conduct a full outer join and remove duplicates (overlap by the natural join) of the MCT's 
dim(Kegg_MCT) # how many annotations from Kegg
dim(biocyc_MCT) # how many annotations from biocyc
sum(dim(Kegg_MCT)[1] + dim(biocyc_MCT)[1]) # how many annotations from both (Total)
Combined_MCT <- merge.data.frame(Kegg_MCT, biocyc_MCT, by = "m.z", all = TRUE)
dim(Combined_MCT) # how many annotations after merge. Value is inflated as shared values
                  # are coerced into duplicates. (Basically, counts the natural join twice!!)
                  # this 'double inner bit' is the part we need to be removed! 
dim(Combined_MCT)[1] - (sum(dim(Kegg_MCT)[1] + dim(biocyc_MCT)[1])) # Number of inflated entries

install.packages("tidyverse")
library(tidyverse)

matched_MCT <- distinct(Combined_MCT, m.z, .keep_all = TRUE) # keep_all keeps the adjoining column info :) 
dim(matched_MCT) # Row totals should be more than kegg but not higher than biocyc. Duplicates from both have been
                 # removed, but so have the dupliates provided WITHIN each database. 

### Demonstrate reasoning for duplication removal:

# test peak list profile ('features' from diffreport) - THEY SHOULD BE ALL UNIQUE
dim(PeakListProfile)[1] # how many m.z entries?
test <- distinct(PeakListProfile, m.z, .keep_all = TRUE) # see if there are duplicates... 
dim(test)[1] # how many after 'removing' duplicates? Should be the same.

# test biocyc_MCT - How many entries are removed? Don't forget, this is results from ONE database, not the combination of the two. 
dim(biocyc_MCT)[1] # how many m.z entries
test <- distinct(biocyc_MCT, m.z, .keep_all = TRUE) # see if there are duplicates... 
dim(test)[1] # how many after removing duplicates. Should be less. 

# Must remove duplicates as otherwise when put into the diffreport, the duplicates would 
# inflate the distribution bias (PLP has one m.z query only). Better to miss out on potential
# names then bias the data; can always loop back to check the most sig. hits. didn't have any other
# possible results. DON'T FORGET THIS FOR THE EDGE_r OUTPUT!!! 

write.table(matched_MCT, file="matched_MCT.txt", sep="\t", row.names = FALSE)

## FINALLY - natural merge the diffreport to the matched MCT to provide annotated results
AnnotatedDiffReport <- merge.data.frame(diffreport, matched_MCT, by = "m.z")
write.table(AnnotatedDiffReport, file="an_av_t_diffReport.txt", sep="\t", row.names = FALSE)

## PET ##
# Pathway Enrichment Table -- Pretty meaningless on the overall data. Plus biocyc or kegg
Kegg_PET<-CreateMummichogAnalTable(mSet) #  Pathway Enrichment Table, Publish-ready, (Follow steps below)
# 1. Copy the output from the above command as displayed in the console. In it's entirity.  
# 2. Once copied, in R studio, File > New File > R Sweave
# 3. Copy the entire output into the default location (row 6) as is.
# 4. Click compile pdf, and give it a file name (not displayed on the table). 
# 5. Note that the table is saved as a pdf and the R sweave file (with raw data) is saved  
#    with the same name, with a .Rnw extention (same directory), which may be re-run at any time.

##### 3.3 Annotated & Averaged & Transformed Clustering Analysis.#####
AnnotatedDiffReport <- read.table("an_av_t_diffReport.txt", header = TRUE) # Load in if needed

## Maybe it'd be worth investigating the principle components further.. 
# PCA on Annotated, Averaged & Transformed data
library(factoextra)

res.pca.AN_AV_T <- prcomp(t(AnnotatedDiffReport[,17:40]), scale = TRUE, center = TRUE)
fviz_pca_ind(res.pca.AN_AV_T,
             col.ind= AV_groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07", "green", "purple"), #Change colour codes / names as desired.
             addEllipses = TRUE,  
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
fviz_eig(res.pca.AN_AV_T) # Generate the scree plot for variance represented by each principal component

# PLS-DA on transformed data (still TR's)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mixOmics", version = "3.8")
library(mixOmics)
plsda_mod_AN_AV_T <- plsda(t(AnnotatedDiffReport[,17:40]), # 
                        AV_groups[,1],# either num/diet, just left for now..  
                        ncomp = 10, # the number of components to include in the model. Default to 2.
                        scale = TRUE, # Boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE)
                        mode = c("regression"),  ##### See stats below (1 of 4 choices)
                        tol = 1e-6, # Convergence stopping value.
                        max.iter = 1000000, # Integer, the maximum number of iterations.
                        near.zero.var = FALSE, # boolean, see the internal nearZeroVar function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE.
                        logratio="none",  # one of "none", "CLR". Specifies the log ratio transformation to deal with compositional values that may arise from specific normalisation in sequencing dadta. Default to 'none'.
                        multilevel=NULL, # Sample information for multilevel decomposition for repeated measurements. A numeric matrix or data frame indicating the repeated measures on each individual, i.e. the individuals ID. See examples in ?splsda.
                        all.outputs = TRUE) # Boolean. Computation can be faster when some specific (and non-essential) outputs are not calculated. Default = TRUE.
plotIndiv(plsda_mod_AN_AV_T, ind.names = AV_groups[,1], ellipse = TRUE, abline = TRUE, legend =TRUE)

# Be warned, this takes about 15-20 mins to run. 
perf.plsda_AN_AV_T <- perf(plsda_mod_AN_AV_T, validation = "Mfold", folds = 5, 
                     progressBar = TRUE, auc = TRUE, nrepeat = 100) 
plot(perf.plsda_AN_AV_T, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

# Pull data to make 'scree'-like plot:
plsda_var_AN_AV_T <- as.data.frame(plsda_mod_AN_AV_T$explained_variance)
plsda_x_var_AN_AV_T <- as.data.frame((plsda_var_AN_AV_T$X[1:10])*100)
colnames(plsda_x_var_AN_AV_T) <- "X"
plsda_x_var_AN_AV_T$Dimensions <- seq(1, 10, by =1)

library(ggplot2)
# Scree plot equivalent for PLS-DA data.
ggplot(data=plsda_x_var_AN_AV_T, aes(x=Dimensions, y=X)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_line(data=plsda_x_var_AN_AV_T, aes(x=Dimensions, y=X), colour="black") +
  geom_point(data=plsda_x_var_AN_AV_T, aes(x=Dimensions, y=X), colour="black") +
  ylab("Percentage of explained variances") +
  ggtitle("PLS-DA X Variances") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal() +
  theme(plot.title = element_text(size=10))


##### 3.4 Individual Group PET testing - TO FIX LATER POSSIBLY?? #####
library(xcms)
library(Rmpi)
library(CAMERA)
# 2>2 - I had to copy the 2-2 data four times for it to work... NO IDEA WHY. Must be a minimum number of raw data files required, but there is no mention of this anywhere :(
#     - Just remove the additional 'copied' columns from the resulting diffreport. No need to average them, they should be identical. IT's the same fucking file afterall. 
#     - I'm getting so cick of this though. Like I like most of this, it's a good challenge, but I'm getting lost with the modelling shit at the end. It's pissing me off.
#     - I can't think of any other way in which I would be able to get predictions though. It is probably right. (Mirre made it after all). But I still don't get what it's used for!!!
#     - Oh, and I've run out of time to ask. Fucking ey. 
setwd("D:/cdf/2-2") 
mzXMLs <- list.files(path="D:/Metabolomics/mzXML/2-2", pattern=".mzXML", full.names = TRUE, recursive = TRUE)
xset <- xcmsSet(method = "matchedFilter", fwhm = 25, snthresh = 1, step = 0.04, steps = 1, 
                  sigma = 10.6166128758281, max = 10, mzdiff = 0.01, index = FALSE)
xset <- retcor(xset, method = "obiwarp", plottype = "none", distFunc = "cor_opt", profStep = 1, 
                  center = 69, response = 1, gapInit = 0.4, gapExtend = 2.7, factorDiag = 2, factorGap = 1,
                  localAlignment = 0)
xset <- group(xset, method  = "density", bw = 38, mzwid = 0.015, minfrac = 0.7, minsamp = 1, max = 50)
xset <- fillPeaks(xset, method="chrom") # 153306...
xset # (info on dataset)  
diffreport <- diffreport(xset)
names(diffreport)[names(diffreport) == 'mzmed'] <- 'm.z'
names(diffreport)[names(diffreport) == 'pvalue'] <- 'p.value'
names(diffreport)[names(diffreport) == 'tstat'] <- 't.score'


# 2>8 
# 8>2 
# 8>8


##### 4.0 Significance testing (on the final DF, An/Av/Tr) #####

# Install packages
install.packages("tidyverse")
install.packages("psych")
install.packages("ggExtra")
install.packages("statmod")

# load libraries
library(tidyverse) 
library(limma) 
library(edgeR)
library(psych) 
library(ggExtra) 
library(scales)
library(statmod) 

# Data notes,
# Best to use the xcms generated 'unique' names (for now). These are arbitrary, but can be traced back easily
# to the information contained in both the biocyc and kegg databases via the av_an_t_diffreport or 
# the MCTs. Note that these are the 'most reduced' lists to prevent over representation (false correlation),
# therefore, it'd be best to check the corresponding m.z of these hits to the values stored in the 
# 'Combined_MCT' object which contains ALL database queries, including isoforms. 
#
# Don't forget, you can always merge the dataframes according to these unique names to generate a full report!#

# Create dataframe for edgeR significance testing
ertesting <- AnnotatedDiffReport[17:40]
head(ertesting)
colnames(ertesting) <- AVFileNames$AVFileNames # Not sure why each sample had an X prefix... oh well fixed now.
rownames(ertesting) <- AnnotatedDiffReport$name
rownames(ertesting)
write.csv(ertesting, "ertesting.csv", row.names = TRUE)

# Check datasets - done previously, but now averaged. Should still be applicable. 
boxplot(log2(ertesting), las = 2, col = c(rep("red", ((Treps/3)/Tgroups)), 
                                           rep("blue", ((Treps/3)/Tgroups)), 
                                           rep("green", ((Treps/3)/Tgroups)), 
                                           rep("yellow", ((Treps/3)/Tgroups))),
        xlab = "Sample", ylab = bquote(log[2]~"Intensity"),
        outcex = 0.5, cex.axis = 1,
        notch = TRUE, main = "AN_AV_T_Data")

# Create colour grouping information for plotDensitites (in order, 1-24):
# Group order (pre sort): 1 = 2-2, 2 = 2-8, 3 = 8-2, 4 = 8-8 | colour: r,b,g,y *key*

pd_av_order <- data.frame(matrix(NA, nrow=24, ncol=0))
pd_av_order$samp <- colnames(ertesting)
pd_av_order$col <- c(rep("1", 6), rep("2", 6), rep("3", 6), rep("4", 6))
pd_av_order <- pd_av_order[order(pd_av_order$samp),] # Re-order the df according to sample number.
pd_av_order$col[pd_av_order$col == "1"] <- "red"
pd_av_order$col[pd_av_order$col == "2"] <- "blue"
pd_av_order$col[pd_av_order$col == "3"] <- "green"
pd_av_order$col[pd_av_order$col == "4"] <- "yellow"

# PlotDensities:
plotDensities(log2(ertesting), col = pd_av_order$col,
              legend = FALSE,
              main = "AN_AV_T_Data")


## Significance testing - what metabolites matter most??
## 4 group analysis ##

# set up the sample mapping
group <- c(rep(Condition1, 6), rep(Condition2, 6), rep(Condition3, 6), rep(Condition4, 6))

# make group into factors and set the order
group <- factor(group, levels = c(Condition1, Condition2, Condition3, Condition4))
str(group)

# create a DGEList object with our data
y_eRD <- DGEList(counts = ertesting, group = group)
y_eRD <- calcNormFactors(y_eRD)
y_eRD <- estimateDisp(y_eRD)

# y_eRD is a list: y_eRD$counts is the data, and y_eRD$samples has interesting content
y_eRD$samples
plotBCV(y_eRD, main = "Biological variation")

# Multigroup significance analysis
edg=DGEList(counts=ertesting, samples = colnames(ertesting), genes = rownames(ertesting))

?DGEList

# Create a new phenotype data frame  
pheno <- matrix(nrow = (Treps/3), ncol = 4 )
colnames(pheno) <- c("Sample", "Past", "Current", "Combined")
pheno[,1] <- colnames(ertesting) # Sample numbering
pheno[,4] <- c(rep(Condition1, 6), rep(Condition2, 6), rep(Condition3, 6), rep(Condition4, 6)) # Current (combined) treatment
pheno[,2] <- c(rep(Treatment1, 12), rep(Treatment2,12)) # Set past treatments (only two before switch)
pheno[,3] <- c(rep(Treatment1, 6), rep(Treatment2,6), rep(Treatment1, 6), rep(Treatment2, 6)) # Set current treatment (still only two as the combined reflects the 4 groups)
pheno <- as.data.frame(pheno)
View(pheno)

# Create the model matrix
design=model.matrix(~1+Past*Current,data=pheno)
rownames(design) <- pheno$Sample
edg <- estimateDisp(edg, design, robust=TRUE) #need to add disp to edg object
edg$common.dispersion # how dispersed(ish) is the data?

# Need to double check linear models and what the coefficients are testing! 
fit <- glmQLFit(edg, design) # fit the model on all of our data
colnames(design)              # Always compared to the first group, in our case Condition1. ! 
lrt_1 <- glmQLFTest(fit,coef=2) # Coef 2 is for water in the past, coef1 would be control but as it's the base for the lm it's just gonna be the differential. 
topTags(lrt_1,n=10) # n = top 10, not sample size! 
lrt_2 <- glmQLFTest(fit,coef=3) ### testing whether the current treatment is significant regardless of past treatment. 
topTags(lrt_2,n=10) # n = top 10, not sample size! 
lrt_3 <- glmQLFTest(fit,coef=4) ### Testing whether the current treatment has any significance taking into account of it's past treatment 
topTags(lrt_3,n=10) # n = top 10, not sample size! 
?glmQLFit

coef2 <- as.data.frame(topTags(lrt_1,n=100))
coef3 <- as.data.frame(topTags(lrt_2,n=100))
coef4 <- as.data.frame(topTags(lrt_3,n=100))

write.csv(coef2, "ncoef2.csv")
write.csv(coef3, "ncoef3.csv")
write.csv(coef4, "ncoef4.csv")




## Significance testing - what metabolites matter most??
## 2 group analysis

lowHIGH <- exactTest(y_eRD, pair = c("2>2", "8>8"))
summary(decideTestsDGE(lowHIGH)) # this counts up, down, and unchanged genes (here it is metabolites)
lowHIGH <- topTags(lowHIGH, n = 10000, sort.by = "none")






## Extra plots for stats-based stuff 

# function computes CVs per time point
make_CVs <- function(df) {
  # separate by samples
  group1 <- df[1:6]
  group2 <- df[7:12]
  group3 <- df[13:18]
  group4 <- df[19:24]
  
  group1$ave <- rowMeans(group1)
  group1$sd <- apply(group1[1:6], 1, sd)
  group1$cv <- 100 * group1$sd / group1$ave
  group2$ave <- rowMeans(group2)
  group2$sd <- apply(group2[1:6], 1, sd)
  group2$cv <- 100 * group2$sd / group2$ave 
  
  group3$ave <- rowMeans(group3)
  group3$sd <- apply(group3[1:6], 1, sd)
  group3$cv <- 100 * group3$sd / group3$ave
  group4$ave <- rowMeans(group4)
  group4$sd <- apply(group4[1:6], 1, sd)
  group4$cv <- 100 * group4$sd / group4$ave 
  
  ave_df <- data.frame(group1$ave, group2$ave, group3$ave, group4$ave)
  sd_df <- data.frame(group1$sd, group2$sd, group3$sd, group4$sd)
  cv_df <- data.frame(group1$cv, group2$cv, group3$cv, group4$cv) 
  return(list(ave_df, sd_df, cv_df))
}

list_ertesting <- make_CVs(ertesting)
boxplot(list_ertesting[[3]], las = 2, notch = TRUE, main = "CVs",
        ylim = c(0, 100), names=c(Condition1,Condition2,Condition3,Condition4))

# print out the average median CVs 
print("eRD (%) (SL then SL/TMM):")
(ertesting_cv <- round(mean(apply(list_ertesting[[3]], 2, median)), 2))

# compare biological replicates to each other by condition
pairs.panels(log2(ertesting[1:6]), lm = TRUE, main = bquote("AN_AV_T Data"~.(Condition1)))
pairs.panels(log2(ertesting[7:12]), lm = TRUE, main = bquote("AN_AV_T Data"~.(Condition2)))
pairs.panels(log2(ertesting[13:18]), lm = TRUE, main = bquote("AN_AV_T Data"~.(Condition3)))
pairs.panels(log2(ertesting[19:24]), lm = TRUE, main = bquote("AN_AV_T Data"~.(Condition4)))

# add marginal distrubution histograms to basic correlation plot (good starting point)
ave_eRD <- data.frame(group1 = rowMeans(ertesting[1:6]), 
                      group2 = rowMeans(ertesting[7:12]), 
                      group3 = rowMeans(ertesting[13:18]), 
                      group4 = rowMeans(ertesting[19:24]))
colnames(ave_eRD) <- c(Condition1,Condition2,Condition3,Condition4) # Fix colnames according to DoE
## If you're interested in different combinations (e.g. 1 v 3) then follow the example below... you get the gist.
# group1 v group2
ggplot()
corr_plot <- ggplot(ave_eRD, aes(x = log10(ave_eRD[,1]), y = log10(ave_eRD[,4]))) +
  geom_point() + ggtitle(bquote(.(Condition1) ~ "vs" ~ .(Condition4) ~ ": eRD"))
ggMarginal(corr_plot, type = "histogram")

























##### Old Code (To delete when done) #####
# Remove the rowname prefixes as well, they're ugly. Do later / ask Mirre :) 
# rownames(PCA_data_TR) <- sub("-0*\\-0", "", rownames(PCA_data_TR))

# Generate CAMERA Annotated Diffreport and export to your working directory (extra annotation plus diffreport)
diffreport <- annotateDiffreport(xset, nSlaves=80) # EIC = Extracted Ion Count
cleanParallel # remove the spawned slave processes. Wonder if it'll actually make things quicker??


## Don't need the next two lines?
pheno_data <- data.frame(sample_name = sub(basename(mzXMLs), pattern = ".mzXMLs",
                                           replacement = "", fixed = TRUE),
                         sample_group = c(rep("2-2", 18), rep("2-8", 18), rep("8-2", 18), rep("8-8", 18)),
                         stringsAsFactors = FALSE)

raw_data <- readMSData(files = mzXMLs, pdata = new("NAnnotatedDataFrame", pheno_data), mode = "onDisk")

# index logical(1) specifying whether indicies should be returned instead of values for m/z and retention times

peakpickingParameters$min_peakwidth <- c(2,5)   # expected approximate peak width in chromatographic 
peakpickingParameters$max_peakwidth <- c(15,20) # space. Given as a range (min, max) in seconds
# As we're working with MALDI data without prior LC/HPLC
# the peakwidth paramaters are a lot smaller as the total
# scan time was 120 seconds.. the default settings could
# span near 50% of the scan causing peak detection issues.
peakpickingParameters$ppm <- c(10,20) # defining the maximal tolerated m/z deviation in consecutive scans
# in parts per million (ppm) for the initial ROI definition. NOT the 
# ppm accuracy of the analyser used.
peakpickingParameters$mzdiff <- c(-0.001, 0.010) # Default
# mzdiff represents the minimum difference in m/z dimension required for peaks with overlapping 
# retention times; can be negative to allow overlap. During peak post-processing, peaks defined 
# to be overlapping are reduced to the one peak with the largest signal.
peakpickingParameters$snthresh <- snthr.range 
peakpickingParameters$noise <- as.numeric((noise[1,2])^2)*snthr.range[1] 
# allowing to set a minimum intensity required for centroids to be
# considered in the first analysis step (centroids with intensity <
# noise are omitted from ROI detection).
peakpickingParameters # Check changes have taken hold.
peakpickingParameters$bw <- c(5,10) # defining the bandwidth (standard deviation of the smoothing kernel)
# to be used
peakpickingParameters$max <- c(1,10)
peakpickingParameters$step <- c(0.2, 0.3)
peakpickingParameters$fwhm <- c(40, 50)

peakpickingParameters$steps <- 2 

peakpickingParameters 

?getDefaultXcmsSetStartingParams


Test <- xcmsRaw(mzXMLs[1], profstep = 1, profmethod = 'bin', profparam = list())
plotScan(Test, 118, mzrange = numeric(), ident = FALSE)
?plotScan
?xcmsRaw
all(sapply(spectra, isRegular)) # TRUE = Profile data, False = Centroided data. Important for XCMS.