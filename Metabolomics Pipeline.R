### Dependencies ### 

# Can copy and paste into the terminal tab below directly in R studio (usually)
# If a dependency is already installed, move on to the next one. All systems are different, and I've covered most linux/ubuntu dependencies.
# Always run the last command (update) and make sure R is up to date (run on R 3.4.4, mmv if otherwise).
# Run on Ubuntu 16.04.3 LTS and 18.04 LTS and Windows 10
# Don't forget your renaming bash script to fix the nomencelature if needed.
sudo apt-get moo
sudo apt-get install fortune
sudo apt-get install cowsay
fortune | cowsay
# Ok, let's get serious:
sudo apt-get install libgtk2.0-dev
sudo apt-get install libcairo2-dev
sudo apt-get install libnetcdf-dev
sudo apt-get install xvfb
sudo apt-get install xauth 
sudo apt-get install xfonts-base 
sudo apt-get install libxt-dev
sudo apt-get install r-base-dev 
sudo apt-get install r-cran-rgl
sudo apt-get install texlive-full
sudo apt-get install cmake
sudo apt-get install libxml2-dev 
sudo apt-get install libssl-dev
sudo apt-get update

## Windows - Not sure every one is required, but install them all for now!
## Make sure to hit install all, users, path, dependencies etc etc in every installer.
https://sourceforge.net/projects/gtk-win/files/latest/download
http://www.unidata.ucar.edu/downloads/netcdf/ftp/netCDF4.6.1-NC4-DAP-64.exe
http://mirror.ctan.org/systems/texlive/tlnet/install-tl-windows.exe ## takes a long time to install.
https://cmake.org/files/v3.12/cmake-3.12.0-win64-x64.msi
https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe

## WINDOWS IMPORTANT - Do NOT compile the following packages from scratch, choose no and it will install correctly.
install.packages("XML")
library(XML)
install.packages("stringi")
library(stringi)


### Finally, All to be run in R ###
### You'll also need an active internet connection to install packages and perform many of the functions with MetaboanalystR ###
### Any sections with [@][@][@] are sections which require manual optimisation with regards to your experiment (sample size, names etc)
### You don't have to save every time a write command is issued, take what you want and use common sense if its needed further on.

rm(list=ls()) # Clear the workspace if required

# Libraries - I've put everything into a function to save time. 
# Install packages via compilation if asked. If libraries fail to load, then rerun the installs for the two packages above without compilation.
install_packages <- function(){
  
  cran_pkg <- c("devtools", "magrittr", "factoextra", "Rserve", "RColorBrewer", "xtable", "som", "ROCR", "RJSONIO", "gplots", "e1071", "caTools", "igraph", "randomForest", "Cairo", "pls", "pheatmap", "lattice", "rmarkdown", "knitr", "data.table", "pROC", "Rcpp", "caret", "ellipse", "scatterplot3d")
  bioconductor_pkg <- c("xcms", "IPO", "CAMERA", "impute", "pcaMethods", "siggenes", "globaltest", "GlobalAncova", "Rgraphviz", "KEGGgraph", "preprocessCore", "genefilter", "SSPA", "sva")
  
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
library(devtools)
devtools::install_github("xia-lab/MetaboAnalystR", build_vignettes=TRUE)
devtools::install_github("romainfrancois/nothing")


library(xcms)
library(CAMERA)
library(IPO)
library(MetaboAnalystR)
library(factoextra)
library(magrittr)


# Load the data - The working directory should be kept the same for each project. [@][@][@]
# The loaded data should be kept in one folder, with a single subfolder containing the files for each of the groups.
setwd("/home/alex/Documents/Metabolomics")
files <- list.files(path="/home/alex/Documents/Metabolomics/cdf", pattern=".CDF", full.names = TRUE, recursive = TRUE)

setwd("D:/Metabolomics")
files <- list.files(path="D:/Metabolomics/cdf", pattern=".CDF", full.names = TRUE, recursive = TRUE)

############################################################################################################
############################################################################################################

# Define your experiment 

# Replace the parameters stored in this section with that of the DoE. The pipeline should then handle the majority of 
# the data wrangling for you. 


# Total number of technical replicates:
Treps <- 72
# Total number of groups analysed: (placeholder, just 4 for now)
Tgroups = 4
# Name the conditions for the experiment in the alphabetical order that they appear in the folder.
Condition1 = "2-2"
Condition2 = "2-8"
Condition3 = "8-2"
Condition4 = "8-8"
# Name of the two treatment types (e.g. control/water, highprotein/lowprotein etc) - a 4 group DoE tends to have
# two treatments, which are then 'switched' for half the samples after a certain time-point.
# Treatment 1 should be the first treatment given to the first alphabetical folder.
Treatment1 = "Low Protein"
Treatment2 = "High Protein"

############################################################################################################
############################################################################################################




# Housekeeping to be conducted ahead of time;

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

############################################################################################################

# Beginning of the pipeline;

# IPO - Optimise the parameters that are to be used by XCMS. 
# Should allow for variability and differences in analyser platforms.
# As a result, weeks of coding were dropped due to this process. Oh well.
# If this doesn't work (no isotopes detected) run the script in the next section. 

## This isn't working perfectly for alternate datasets. I have fixed it so that isotopes are now detected; but the best parameters are still
## those that I have previously used for my dataset. Needs attention at a later point; however, it takes hours to test a single
## combination, hence why I haven't done anything more so far. 

### TL;DR, skip to line 226 :) 
   # ?centWave
   # ?matchedFilter
   # ?getDefaultRetGroupStartingParams
## Optimise Peak Picking Parameters   #matchedFilter
peakpickingParameters <- getDefaultXcmsSetStartingParams('matchedFilter') # Maybe try without specifying which? Would it try both? But then what's the point.. 
  #setting levels for step to 0.2 and 0.3 (hence 0.25 is the center point)
peakpickingParameters$bw <- c(5,10) # only 1 value
peakpickingParameters$max <- c(1,10)
peakpickingParameters$step <- c(0.2, 0.3)
peakpickingParameters$fwhm <- c(40, 50)
peakpickingParameters$ppm <- c(4,10) ## IMPORTANT ## PPM accuracy of the analyser! Can enter a range, but best to specify if known.
peakpickingParameters$steps <- 2 #setting only one value for steps therefore this parameter is not optimized
peakpickingParameters$snthresh <- c(1,10) # max two values
peakpickingParameters ## Check the parameters if desired. 

## I've expanded the parameter selection (rather significantly) in an attempt to improve isotope detection tolerance.
## However, this does negatively impact the performance as an increased number of combinations must be tested. 

# says sigma 0 and index false, so maybe I can give centwave a whirl if I specify these??
# The old settings only gave 2... yes 2! mz slices, so I've tweaked the parameters and provided a range for ppm.
# Taking even longer than before, numbers going higher than 25 at least - hopefully that'll equate to a better result this time :) 
# a lot longer than before!! 

# taking hours on an i7-8550u... maybe it's best for the workstation :/ 
# then it crashed (guessing its kernel panic) was doing other work though :(
# It went up to around 72 each time, and was on the third iteration before it fucked it. 
# At least I was working over teamviewer for the rest of it LOL <3 
#
#Using the windows-workstation for now with 12 cores / 24 thread at 3.6ghz with 80gb ram. Should do.
time.xcmsSet <- system.time({ # measuring time
  resultPeakpicking <- 
    optimizeXcmsSet(files = files, 
                    params = peakpickingParameters, 
                    BPPARAM = SnowParam(), #linux MulicoreParam() #windows SnowParam()
                    subdir = "IPO",
                    plot = TRUE)
})


## Optimise peak picking result
resultPeakpicking$best_settings$result
optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset

## Optimise retention time correction and grouping parameters
retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- 1
retcorGroupParameters$gapExtend <- 2.7
time.RetGroup <- system.time({ # measuring time
  resultRetcorGroup <-
    optimizeRetGroup(xset = optimizedXcmsSetObject, 
                     params = retcorGroupParameters, 
                     nSlaves = 4, 
                     subdir = NULL,
                     plot = TRUE)
})


## Running times and session info (If interested - It takes 5-20 minutes to run)  HEREHEREHEREHEREHERE
time.xcmsSet # time for optimizing peak picking parameters
time.RetGroup # time for optimizing retention time correction and grouping parameters

sessionInfo()

## Write R script -- Copy the response from the terminal in the section below! 
writeRScript(resultPeakpicking$best_settings$parameters, 
             resultRetcorGroup$best_settings)

# New R script:
xset <- xcmsSet( 
  method   = "matchedFilter",
  fwhm     = 9,
  snthresh = 1,
  step     = 0.2014,
  steps    = 2,
  sigma    = 3.82198063529811,
  max      = 13,
  mzdiff   = 0.3972,
  index    = FALSE)
xset <- retcor( 
  xset,
  method         = "obiwarp",
  plottype       = "none",
  distFunc       = "cor_opt",
  profStep       = 1,
  center         = 65,
  response       = 1,
  gapInit        = 0.24,
  gapExtend      = 2.7,
  factorDiag     = 2,
  factorGap      = 1,
  localAlignment = 0)
xset <- group( 
  xset,
  method  = "density",
  bw      = 22,
  mzwid   = 0.035,
  minfrac = 0.3,
  minsamp = 1,
  max     = 50)

xset <- fillPeaks(xset)


## Optimised R Script: [@][@][@] -- Delete the following section and paste the new output in. 
##                               -- Unless there's an error using IPO. Then run this; tailored to the MALDI in Heather's lab. 
#                                -- Found 65690 mz slices... still better despite fixing IPO?! 
#                                -- and doesn't error; has enough info to correctly group the samples. 

xset <- xcmsSet( 
  method   = "matchedFilter",
  fwhm     = 8.88,
  snthresh = 1,
  step     = 0.18,
  steps    = 2,
  sigma    = 3.77102089349414,
  max      = 5,
  mzdiff   = 0.44,
  index    = FALSE)
xset <- retcor( 
  xset,
  method         = "obiwarp",
  plottype       = "none",
  distFunc       = "cor_opt",
  profStep       = 1,
  center         = 5,
  response       = 1,
  gapInit        = 0,
  gapExtend      = 2.7,
  factorDiag     = 2,
  factorGap      = 1,
  localAlignment = 0)
xset <- group( 
  xset,
  method  = "density",
  bw      = 22,
  mzwid   = 0.035,
  minfrac = 0.3,
  minsamp = 1,
  max     = 50)

xset <- fillPeaks(xset)

## Generate CAMERA Annotated Diffreport and export to your working directory (extra annotation plus diffreport)
diffreport <- annotateDiffreport(xset)
names(diffreport)[names(diffreport) == 'mzmed'] <- 'm.z'
names(diffreport)[names(diffreport) == 'pvalue'] <- 'p.value'
names(diffreport)[names(diffreport) == 'tstat'] <- 't.score'

write.table(diffreport, file="diffreport.txt", sep="\t", row.names = FALSE)
# diffreport <- read.table("diffreport.txt", header = TRUE) # Read the data back in at a later date if you'd prefer

################# Data Analysis ################

# Make a PCA for the non-annotated data (technical replicates) to see whether they group correctly;
# Technical replicates (e.g. 1,2,3 or 10,11,12) should reoughly cluster in the PCA. If not, then they may be dodgy. 
# Haven't done anything to deal with / remove them for the moment...

## Need to Subselect and transpose the data.
PCA_data_TR <- t(diffreport[,17:(16+Treps)]) 

# Data frame
PCA_data_TR <- as.data.frame(PCA_data_TR)

# FIX any NA values present if any.
PCA_data_TR[is.na(PCA_data_TR)] <- 0
sum(is.na(PCA_data_TR)) # check that there aren't any more NAs 

# Remove the rowname prefixes as well, they're ugly. Do later / ask Mirre :) 
# rownames(PCA_data_TR) <- sub("-0*\\-0", "", rownames(PCA_data_TR))
 
# Create the PCA on the TR.
res.pca.TR <- prcomp(PCA_data_TR, scale = TRUE, center = TRUE) 
fviz_pca_ind(res.pca.TR,
             col.ind= TR_groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07", "green", "purple"), #Change colour codes / names as desired.
             addEllipses = TRUE,  
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
                  # Note, if the names are getting in the way of the PCA, then rename the file by the simple 2 digit 
                  # number. I have a bash script to do this in linux to save time, just ask.

# Ok, if the repliates roughly group then continue to the next part; averaging and annotating the peaks.
# If they don't, then figure it out. The pipeline is just set up for triplicates for now, so can't remove
# problem queries and average two in a specific group instead of three. Possibly something I'll work on in the 
# future, after I extend the capabilities. 

# Annotation
PeakListProfile <- diffreport[,c(6,4,3)]
names(PeakListProfile)[names(PeakListProfile) == 'mzmed'] <- 'm.z'
names(PeakListProfile)[names(PeakListProfile) == 'pvalue'] <- 'p.value'
names(PeakListProfile)[names(PeakListProfile) == 'tstat'] <- 't.score'
write.table(PeakListProfile, file="PeakListProfile.txt", sep="\t", row.names = FALSE)

# Get the annotation result from metaboanalyst - need an active internet connection!
mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
mSet <- Read.PeakListData(mSet, "diffreport.txt");
mSet <- UpdateMummichogParameters(mSet, "6.0", "positive", 0.05); # PPM, ion mode, significance value threshold << SET APPROPRIATELY!
mSet <- SanityCheckMummichogData(mSet)
mSet <- PerformMummichog(mSet, "dme_kegg", "fisher", "gamma") # organism_database, test, distribution.

# Pathway Enrichment Table 
PET<-CreateMummichogAnalTable(mSet) #  Pathway Enrichment Table, Publish-ready, (Follow steps below)
# 1. Copy the output from the above command as displayed in the console. In it's entirity.  
# 2. Once copied, in R studio, File > New File > R Sweave
# 3. Copy the entire output into the default location (row 6) as is.
# 4. Click compile pdf, and give it a file name (not displayed on the table). 
# 5. Note that the table is saved as a pdf and the R sweave file (with raw data) is saved  
#    with the same name, with a .Rnw extention (same directory), which may be re-run at any time.

# Matched Compound Table - Now with common name.
mSet$matches.res$Common.Name<-doKEGG2NameMapping(mSet$matches.res$Matched.Compound)
names(mSet$matches.res)[names(mSet$matches.res) == 'Query.Mass'] <- 'm.z'
write.table(mSet$matches.res, file="MatchedCompoundTable.txt", sep="\t", row.names = FALSE)
MCT<-mSet$matches.res

# Merge the dataframes back together to give completed report 
AnnotatedDiffReport <- merge.data.frame(diffreport, MCT, by = "m.z")
write.table(AnnotatedDiffReport, file="AnnotatedDiffReport.txt", sep="\t", row.names = FALSE)
AnnotatedDiffReport <- read.table("AnnotatedDiffReport.txt", header = TRUE) # Load in if needed

# Create the PCA data frame - The column numbers will be different if not using a total of 72 samples. 
PCA_names <- AnnotatedDiffReport[,(16+Treps+7)]
PCA_intensities <- AnnotatedDiffReport[,17:(16+Treps)]  
PCA_data_TR_AN <- t(PCA_intensities)

# Convert row 1 to column names: (-> Compound Name)
colnames(PCA_data_TR_AN) <- PCA_names
#PCA_data_TR_AN <- PCA_data_TR_AN[-1,]
write.table(PCA_data_TR_AN, file="PCA_data_TR_AN.txt", sep="\t", row.names = FALSE)
   ## You may want to save this as it has the compound names matched to intensities,
   ## I need to remove the names for now to remove the duplicates. This won't work otherwise as the names aren't unique..
   ## Remember, it's not just the NA columns which are non-unique, overrepresenting a single m.z query. 
   ## ---///--- look into a way of combining the column names for those that are removed? May be an idea to use the compound ID number too. 

# Create a PCA for the annotated Technical Replicates. Data should have been reduced.
# FIX any NA values present if any.
PCA_data_TR_AN[is.na(PCA_data_TR_AN)] <- 0
sum(is.na(PCA_data_TR_AN)) # check that there aren't any more NAs 
PCA_data_TR_AN <- as.data.frame(PCA_data_TR_AN)

# remove duplicates and plot PCA --- may remove??
PCA_data_TR_AN <- PCA_data_TR_AN[!duplicated(as.list(PCA_data_TR_AN))]
write.csv(PCA_data_TR_AN, "PCA_data_TR_AN_nodups.csv", row.names = TRUE)

res.pca.TR_AN <- prcomp(PCA_data_TR_AN, scale = TRUE, center = TRUE) ## New error! Damn it :/ 
fviz_pca_ind(res.pca.TR_AN,
             col.ind= TR_groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07", "green", "purple"), #Change colour codes / names as desired.
             addEllipses = TRUE,  
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)


# Now time to average the technical replicates and repeat the PCA.
PCA_AN_data <- rowsum(PCA_data_TR_AN, rep(1:(Treps/3), each=3))/3 
nrow(PCA_AN_data)

# Generate a numeric list of the technical reps in the working directory #
#||| Ok this needs to be changed each time, or input a number list yourself. 
FileNames <- list.files(pattern = ".CDF", recursive = TRUE)
FileNames <- as.data.frame(FileNames)
rownames(FileNames) <- FileNames[,1]
FileNames <- FileNames[,-1]

FileNames <- rownames(FileNames) %>%
  strsplit( "cdf/2-2/" ) 
FileNames <- as.data.frame(FileNames)
FileNames <- FileNames[-1,]
colnames(FileNames) <- FileNames[2,]
FileNames <- t(FileNames)
FileNames <- as.data.frame(FileNames)
rownames(FileNames) <- FileNames[,1]
FileNames <- FileNames[,-1]

FileNames <- rownames(FileNames)%>%
  strsplit( "cdf/2-8/" )
FileNames <- as.data.frame(FileNames)
FileNames <- FileNames[-1,]
colnames(FileNames) <- FileNames[2,]
FileNames <- t(FileNames)
FileNames <- as.data.frame(FileNames)
rownames(FileNames) <- FileNames[,1]
FileNames <- FileNames[,-1]

FileNames <- rownames(FileNames)%>%
  strsplit( "cdf/8-2/" )
FileNames <- as.data.frame(FileNames)
FileNames <- FileNames[-1,]
colnames(FileNames) <- FileNames[2,]
FileNames <- t(FileNames)
FileNames <- as.data.frame(FileNames)
rownames(FileNames) <- FileNames[,1]
FileNames <- FileNames[,-1]

FileNames <- rownames(FileNames)%>%
  strsplit( "cdf/8-8/" )
FileNames <- as.data.frame(FileNames)
FileNames <- FileNames[-1,]
colnames(FileNames) <- FileNames[2,]
FileNames <- t(FileNames)
FileNames <- as.data.frame(FileNames)
rownames(FileNames) <- FileNames[,1]
FileNames <- FileNames[,-1]

rownames(FileNames) <- substr(rownames(FileNames),1,nchar(rownames(FileNames))-4) # delete the last 4 characters (.CDF)

# Generate the a new list correlating to the correct number formatting after averaging the technical replicates
FileNames = FileNames[seq(3, nrow(FileNames), 3), ]
FileNames[,1] <- rownames(FileNames)
FileNames <- as.numeric(FileNames[,1])/3
FileNames <- as.data.frame(FileNames)

# Correct the Sample ID numbering in the dataframe
rownames(PCA_AN_data) <- FileNames[,1]
write.table(PCA_AN_data, file="PCA_AN_data.txt", sep="\t", row.names = TRUE) 

# PCA - Grouped Annotated Data
res.pca_AN_data <- prcomp(PCA_AN_data, scale = TRUE, center = TRUE) 

fviz_pca_ind(res.pca_AN_data,
             col.ind= AV_groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07", "green", "purple"), #Change to what you want :) 
             addEllipses = TRUE,  
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

# Significance testing // code may be broken
resultsP=matrix(NA,ncol=2,nrow=dim(PCA_AN_data)[2])
p=1
while(p<(dim(PCA_AN_data)[2]+1))
  {
StatsResults <- anova(lm(PCA_AN_data[,p]~factor(AV_groups))) # column 1
PP=StatsResults$'Pr(>F)'[1]
resultsP[p,]=c(names(PCA_AN_data)[p],PP)
p=p+1
}
hist(as.numeric(as.character(resultsP[,2])))
plot(PCA_AN_data[,which(as.numeric(as.character(resultsP[,2]))==min(as.numeric(as.character(resultsP[,2]))))]~factor(AV_groups))


# What are the significant hit(s) if any? (Little manual input required, just repeat the final line if multiple compounds identified).
min(as.numeric(as.character(resultsP[,2]))) # list hit(s)
which(as.numeric(as.character(resultsP[,2]))==min(as.numeric(as.character(resultsP[,2])))) # Row position 
resultsP[471,] # Compound and P value -- does it make biological sense?


# Multi PCA - The hunt for significance
            # Boxplots for each principal component, p.val in title. (PC#, Pvalue)


Groups<-factor(AV_groups) 
p=1
while(p<21){
  plot(res.pca_AN_data$x[, p] ~ Groups)
  title(paste(p,anova(lm(res.pca_AN_data$x[, p] ~ AV_groups))[5]))
  p=p+1}

# Compare any interesting principlal components (alter as required) #18,8
## PC5 v PC11    
fviz_pca_ind(res.pca_AN_data,
             col.ind= AV_groups, # color by groups
             axes = c(8,18),
             palette = c("#00AFBB",  "#FC4E07", "green", "purple"), #Change to what you want :) 
             addEllipses = TRUE,  
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

# Plot histograms for each sample to show the intensity distribution 
# (note, only shows majoirity of hits not all, x axis limited to 1000, higher intensities are present in low frequency)                                               

matrixPCA_AN_data=as.matrix(PCA_AN_data)
dfGroups <- as.data.frame(Groups)
par(mfrow = c(4, 6))
p=1 # It'd be nice if I could export these automatically may need to restructure to ggplot and ggsave in the loop
while(p<25){
  hist(matrixPCA_AN_data[p,],
     xlab = "Intensity (CPS)",
     xlim = c(0,10000),
     breaks = 20000,
     main = bquote("Sample" ~ .(FileNames[p,1]) *"," ~ group == .(paste(dfGroups[p,1])))
     )
  p=p+1
}
par(mfrow = c(1, 1))

## May just remove this next section, not based on a group by group basis
## can't split the data as it's simply based on names... which are identified in all?
## Hopefuly the edgeR stuff works just fine.. 



# Pathway Enrichment Analysis, Calculation - Based on whether a compound is present or absent, not accounting for intensity values.
# Intensity-based analysis in edgeR to follow, but this will give some overview of what's found in the data as a whole.
mSet <- InitDataObjects("conc", "pathora", FALSE)
cmpd.vec <- AnnotatedDiffReport$Common.Name
mSet <- Setup.MapData(mSet, cmpd.vec);
mSet <- CrossReferencing(mSet, "name");
mSet <- CreateMappingResultTable(mSet)
mSet <- SetKEGG.PathLib(mSet, "dme") #dme = drosophila, as before, but can specify others.
mSet <- SetMetabolomeFilter(mSet, F);
mSet <- CalculateOraScore(mSet, "rbc", "hyperg")

# Pathway Enrichment Analysis, Results
## Summary Table
Pathway_Results <- read.csv("pathway_results.csv")
colnames(Pathway_Results)[1] <-c("Pathway Name")
View(Pathway_Results)
## Summary Graphic  (Note the following graphics are saved directly to the working directory and are not previewed)
mSet<-PlotPathSummary(mSet, "< name >", "png", 72, width=NA) # This is a graphical representation of the pathways... need to add labels (and white background?), but this will do for now. 
#                                                      |           X = Pathway impact factor, Y = p value
#                                                      |
## Individual Pathway Graphics (put the corresponding \./ number from the viewer in the command below, where indicated)
mSet<-PlotMetPath(mSet, Pathway_Results$`Pathway Name`[1], 528, 480)
## Examples by anually setting the name... what a chore. 
mSet<-PlotMetPath(mSet, "One carbon pool by folate", 528, 480) 
mSet<-PlotMetPath(mSet, "Nicotinate and nicotinamide metabolism", 528, 480)
mSet<-PlotMetPath(mSet, "Terpenoid backbone biosynthesis", 528, 480)
# About compound colors within the pathway - light blue means those metabolites are not in your data and are used as background 
# for enrichment analysis; grey means the metabolite is not in your data and is also excluded from enrichment analysis (only applicable
# if you have used a custom metabolome profile); other colors (varying from yellow to red) means the metabolites are in the data with
# different levels of significance (red are high).

# Part II - edgeR analysis
# Plug in the edgeR stuff - I'll keep the packages etc separate for now until I can test the windows compatibilty.

# At this point, the number of loaded packages in ridiculous; so much so in fact, that R may error with:
# "maximal number of DLLs reached...". Don't let this happen to you. Save what you need and then hit restart
# R from within R studio (Session>Restart R). Or just update a loaded package; that'll allow rstudio to save your work, 
# update and reload your entire environment.

# Install packages
install.packages("tidyverse")
install.packages("psych")
install.packages("ggExtra")
install.packages("statmod")

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("limma")
biocLite("Rcpp")
biocLite("scales")
biocLite("digest")

# load libraries
library(tidyverse) 
library(limma) 
library(edgeR)
library(psych)
library(ggExtra)
library(scales)
library(statmod)

# Create a new dataframe for use in edgeR
edgeR_Data <- as.data.frame(t(PCA_AN_data)) # Reformat the previous data so it's suitable for edgeR analysis
write.csv(edgeR_Data, "edgeR_initial_data.csv", row.names = TRUE) # Write starting dframe to file if desired

# Load in edgeR_data if saved previously for reuse
setwd("D:/Metabolomics/")
edgeR_Data <- read_csv("edgeR_initial_data.csv") # Load if required
edgeR_Data <- as.data.frame(edgeR_Data)
rownames(edgeR_Data) <- make.names(edgeR_Data[,1],TRUE)
edgeR_Data <- edgeR_Data[,-1]

# finish initial data handling steps 
edgeR_names <- as.data.frame(rownames(edgeR_Data)) # just incase they're needed later on.
dim(edgeR_Data)

# take care of any zero values // Is this right? surely 1 =/ 0... 
edgeR_Data[edgeR_Data <= 1] <- 1

# let's see what the starting data looks like || Also red will always be group 1, group 2 blue, group 3 green, group 4 yellow. (Groups should be set in alphabetical order as shown in the working directory)
# par(mfrow = c(1, 1)) // ok this should be fixed for now (automatic regardless)
boxplot(log2(edgeR_Data), las = 2, col = c(rep("red", ((Treps/3)/Tgroups)), 
                                           rep("blue", ((Treps/3)/Tgroups)), 
                                           rep("green", ((Treps/3)/Tgroups)), 
                                           rep("yellow", ((Treps/3)/Tgroups))),
        xlab = "Sample", ylab = bquote(log[2]~"Intensity"),
        notch = TRUE, main = "RAW edgeR data")
# NOTE: density distributions order samples differently than box plots... Adjusted accordingly.. 
plotDensities(log2(edgeR_Data), col = c(rep("red", 1), 
                                        rep("blue", 3), 
                                        rep("green", 6), 
                                        rep("yellow", 1), 
                                        rep("red", 1), 
                                        rep("yellow", 5), 
                                        rep("red", 4), 
                                        rep("blue", 3)), 
              legend = FALSE,
              main = "Raw edgeR data")

## The legend is being a pain right now (too many entries). 
## Come back at a later point; put legend = FALSE in the plotDensities arguement, then define the legend afterwards. 

# check the column totals to see if they are (roughly) equal
print("edgeR_Data:")
format(round(colSums(edgeR_Data), digits = 0), big.mark = ",")

# figure out the global scaling values for sample loading normalizations
target_eRD <- mean(colSums(edgeR_Data))

# do the sample loading normalization before other normalizations
# there is a different correction factor for each column
norm_facs_eRD <- target_eRD / colSums(edgeR_Data)
data_eRD_sl <- sweep(edgeR_Data, 2, norm_facs_eRD, FUN = "*")


# par(mfrow = c(2, 2))
# see what the SL normalized data look like
boxplot(log2(data_eRD_sl), las = 2, col = c(rep("red", ((Treps/3)/Tgroups)), 
                                           rep("blue", ((Treps/3)/Tgroups)), 
                                           rep("green", ((Treps/3)/Tgroups)), 
                                           rep("yellow", ((Treps/3)/Tgroups))),
        xlab = "Sample", ylab = bquote(log[2]~"Intensity"),
        notch = TRUE, main = "sl edgeR data")


# NOTE: density distributions order samples differently than box plots...
plotDensities(log2(data_eRD_sl), col = c(rep("red", 1), 
                                        rep("blue", 3), 
                                        rep("green", 6), 
                                        rep("yellow", 1), 
                                        rep("red", 1), 
                                        rep("yellow", 5), 
                                        rep("red", 4), 
                                        rep("blue", 3)), 
              legend = FALSE,
              main = "sl edgeR data") 
# check the columnn totals (should be more like equal than before)
format(round(colSums(data_eRD_sl), digits = 0), big.mark = ",")

# do TMM on MQ data
eRD_tmm <- calcNormFactors(data_eRD_sl)
data_eRD_tmm <- sweep(data_eRD_sl, 2, eRD_tmm, FUN = "/") # this is data after SL and TMM on original scale
?calcNormFactors
boxplot(log2(data_eRD_tmm), las = 2, col = c(rep("red", ((Treps/3)/Tgroups)), 
                                            rep("blue", ((Treps/3)/Tgroups)), 
                                            rep("green", ((Treps/3)/Tgroups)), 
                                            rep("yellow", ((Treps/3)/Tgroups))),
        xlab = "Sample", ylab = bquote(log[2]~"Intensity"),
        notch = TRUE, main = "sl/tmm edgeR data")

plotDensities(log2(data_eRD_tmm), col = c(rep("red", 1), 
                                         rep("blue", 3), 
                                         rep("green", 6), 
                                         rep("yellow", 1), 
                                         rep("red", 1), 
                                         rep("yellow", 5), 
                                         rep("red", 4), 
                                         rep("blue", 3)), 
              legend = FALSE,
              main = "sl/tmm edgeR data") 

# Go back and pull the plots together for write up
par(mfrow = c(1, 2))
par(mfrow = c(1, 1))

# check final column totals after TMM
format(round(colSums(data_eRD_tmm), digits = 0), big.mark = ",")

# see how things cluster after we have gotten the boxplots and desity plots looking nice
# (Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples).
# Sample 10 (controlADD) doesn't cluster well, ald some other controlADDs don't either.. oh well gives an idea. 
par(mfrow = c(1, 3))
plotMDS(log2(edgeR_Data), col = c(rep("red", ((Treps/3)/Tgroups)), 
                                    rep("blue", ((Treps/3)/Tgroups)), 
                                    rep("green", ((Treps/3)/Tgroups)), 
                                    rep("yellow", ((Treps/3)/Tgroups))), 
        main = "Raw edgeR Data")
plotMDS(log2(data_eRD_sl), col = c(rep("red", ((Treps/3)/Tgroups)), 
                                    rep("blue", ((Treps/3)/Tgroups)), 
                                    rep("green", ((Treps/3)/Tgroups)), 
                                    rep("yellow", ((Treps/3)/Tgroups))), 
        main = "SL edgeR Data")
plotMDS(log2(data_eRD_tmm), col = c(rep("red", ((Treps/3)/Tgroups)), 
                                    rep("blue", ((Treps/3)/Tgroups)), 
                                    rep("green", ((Treps/3)/Tgroups)), 
                                    rep("yellow", ((Treps/3)/Tgroups))), 
        main = "SL/TMM edgeR Data")
par(mfrow = c(1, 1))
?plotMDS

# function computes CVs per time point ### TRY TO DO all 4 ## Worked :) // need to automate
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

# get CVs and averages
list_raw <- make_CVs(edgeR_Data)
list_eRD_sl <- make_CVs(data_eRD_sl)
list_eRD_tmm <- make_CVs(data_eRD_tmm)

# compare CV distributions - coefficient of variance. 
par(mfrow = c(1, 3))
boxplot(list_raw[[3]], las = 2, notch = TRUE, main = "Raw CVs",
        ylim = c(0, 100), names=c(Condition1,Condition2,Condition3,Condition4))
boxplot(list_eRD_sl[[3]], las = 2, notch = TRUE, main = "eRD/SL CVs",
        ylim = c(0, 100), names=c(Condition1,Condition2,Condition3,Condition4))
boxplot(list_eRD_tmm[[3]], las = 2, notch = TRUE, main = "eRD/SL/TMM CVs",
        ylim = c(0, 100), names=c(Condition1,Condition2,Condition3,Condition4))
par(mfrow = c(1, 1))

# print out the average median CVs ## Maybed a good metric to include with the figures... 
print("eRD (%) (SL then SL/TMM):")
(raw_med_cv <- round(mean(apply(list_raw[[3]], 2, median)), 2))
(eRD_sl_med_cv <- round(mean(apply(list_eRD_sl[[3]], 2, median)), 2))
(eRD_tmm_med_cv <- round(mean(apply(list_eRD_tmm[[3]], 2, median)), 2))

# compare biological replicates to each other by condition
par(mfrow = c(2, 2))
pairs.panels(log2(data_eRD_tmm[1:6]), lm = TRUE, main = bquote("eRD SL/TMM"~.(Condition1)))
pairs.panels(log2(data_eRD_tmm[7:12]), lm = TRUE, main = bquote("eRD SL/TMM"~.(Condition2)))
pairs.panels(log2(data_eRD_tmm[13:18]), lm = TRUE, main = bquote("eRD SL/TMM"~.(Condition3)))
pairs.panels(log2(data_eRD_tmm[19:24]), lm = TRUE, main = bquote("eRD SL/TMM"~.(Condition4)))
par(mfrow = c(1, 1))

# add marginal distrubution histograms to basic correlation plot (good starting point)
ave_eRD <- data.frame(group1 = rowMeans(data_eRD_tmm[1:6]), 
                      group2 = rowMeans(data_eRD_tmm[7:12]), 
                      group3 = rowMeans(data_eRD_tmm[13:18]), 
                      group4 = rowMeans(data_eRD_tmm[19:24]))
colnames(ave_eRD) <- c(Condition1,Condition2,Condition3,Condition4) # Fix colnames according to DoE
## If you're interested in different combinations (e.g. 1 v 3) then follow the example below... you get the gist.
# group1 v group2
ggplot()
corr_plot <- ggplot(ave_eRD, aes(x = log10(ave_eRD[,1]), y = log10(ave_eRD[,2]))) +
  geom_point() + ggtitle(bquote(.(Condition1) ~ "vs" ~ .(Condition2) ~ ": eRD"))
ggMarginal(corr_plot, type = "histogram")
# group3 v group4
ggplot()
corr_plot <- ggplot(ave_eRD, aes(x = log10(ave_eRD[,3]), y = log10(ave_eRD[,4]))) +
  geom_point() + ggtitle(bquote(.(Condition3) ~ "vs" ~ .(Condition4) ~ ": eRD"))
ggMarginal(corr_plot, type = "histogram")

# let's do DE testing with edgeR
# set up the sample mapping
group <- c(rep(Condition1, 6), rep(Condition2, 6), rep(Condition3, 6), rep(Condition4, 6))

# make group into factors and set the order
group <- factor(group, levels = c(Condition1, Condition2, Condition3, Condition4))
str(group)

# create a DGEList object with our data
y_eRD <- DGEList(counts = data_eRD_tmm, group = group)
y_eRD <- calcNormFactors(y_eRD)
y_eRD <- estimateDisp(y_eRD)

# y_eRD is a list: y_eRD$counts is the data, and y_eRD$samples has interesting content
y_eRD$samples
plotBCV(y_eRD, main = "Biological variation eRD")
?plotBCV

# Multigroup significance analysis
edg=DGEList(counts=edgeR_Data, samples = colnames(data_eRD_tmm), genes = rownames(data_eRD_tmm))

# Create a new phenotype data frame // I really need to fix the rep counts to automated.. 
pheno <- matrix(nrow = (Treps/3), ncol = 4 )
colnames(pheno) <- c("Sample", "Past", "Current", "Combined")
pheno[,1] <- colnames(edgeR_Data) # Sample numbering
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

coef2 <- as.data.frame(topTags(lrt_1,n=100))
coef3 <- as.data.frame(topTags(lrt_2,n=100))
coef4 <- as.data.frame(topTags(lrt_3,n=100))

write.csv(coef2, "coef2.csv")
write.csv(coef3, "coef3.csv")
write.csv(coef4, "coef4.csv")





# As the pipeline is stil unfinished, the code snippets below repflect things which are a still a work in progress. 







?glmQLFit
### so nothing much here, but tidy up the end of this pipeline, and see if I have anything interesting. Hand the data off back to Andy sad for him. 
### double check intensity alues, and redo on the raw data or one transformed? may be overfitted. 


# Ok, the rest of the code (below) may be interesting in the future,
# but there's certainly no point trying to implement it (as is) without
# ANY significant candidates! 








####################################################

?exactTest

# the exact test object has columns like fold-change, CPM, and p-values
# test all combinations between the 4 groups. The test is set for parwise comparisons. 
# not elegant, but does the job. 
et_control_controlADD_eRD <- exactTest(y_eRD, pair = c("control", "controlADD"))
et_control_water_eRD <- exactTest(y_eRD, pair = c("control", "water"))
et_control_waterREMOVE_eRD <- exactTest(y_eRD, pair = c("control", "waterREMOVE"))
et_cotrolADD_water_eRD <- exactTest(y_eRD, pair = c("controlADD", "water"))
et_controlADD_waterREMOVE_eRD <- exactTest(y_eRD, pair = c("controlADD", "waterREMOVE"))
et_water_waterREMOVE_eRD <- exactTest(y_eRD, pair = c("water", "waterREMOVE"))

summary(decideTestsDGE(et_control_controlADD_eRD)) # this counts up, down, and unchanged genes (here it is metabolites)
summary(decideTestsDGE(et_control_water_eRD)) # this counts up, down, and unchanged genes (here it is metabolites)
summary(decideTestsDGE(et_control_waterREMOVE_eRD)) # this counts up, down, and unchanged genes (here it is metabolites)
summary(decideTestsDGE(et_cotrolADD_water_eRD)) # this counts up, down, and unchanged genes (here it is metabolites)
summary(decideTestsDGE(et_controlADD_waterREMOVE_eRD)) # this counts up, down, and unchanged genes (here it is metabolites)
summary(decideTestsDGE(et_water_waterREMOVE_eRD)) # this counts up, down, and unchanged genes (here it is metabolites)

# ok but no results here :'( 

# the topTags function adds the BH FDR values to an exactTest data frame. Make sure not to change row order!
tt_control_controlADD_eRD <- topTags(et_control_controlADD_eRD, n = 10000, sort.by = "none")
tt_control_water_eRD <- topTags(et_control_water_eRD, n = 10000, sort.by = "none")
tt_control_waterREMOVE_eRD <- topTags(et_control_waterREMOVE_eRD, n = 10000, sort.by = "none")
tt_controlADD_water_eRD <- topTags(et_cotrolADD_water_eRD, n = 10000, sort.by = "none")
tt_controlADD_waterREMOVE_eRD <- topTags(et_controlADD_waterREMOVE_eRD, n = 10000, sort.by = "none")
tt_water_waterREMOVE_eRD <- topTags(et_water_waterREMOVE_eRD, n = 10000, sort.by = "none")

tt_control_controlADD_eRD <- tt_control_controlADD_eRD$table # tt_sl is a list. We just need the data frame table
tt_control_water_eRD <- tt_control_water_eRD$table # tt_sl is a list. We just need the data frame table
tt_control_waterREMOVE_eRD <- tt_control_waterREMOVE_eRD$table # tt_sl is a list. We just need the data frame table
tt_controlADD_water_eRD <- tt_controlADD_water_eRD$table # tt_sl is a list. We just need the data frame table
tt_controlADD_waterREMOVE_eRD <- tt_controlADD_waterREMOVE_eRD$table # tt_sl is a list. We just need the data frame table
tt_water_waterREMOVE_eRD <- tt_water_waterREMOVE_eRD$table # tt_sl is a list. We just need the data frame table


# add the default value as a new column
tt_control_controlADD_eRD$candidate <- "no"
tt_control_controlADD_eRD[which(tt_control_controlADD_eRD$FDR <= 0.10 & tt_control_controlADD_eRD$FDR > 0.05), dim(tt_control_controlADD_eRD)[2]] <- "low"
tt_control_controlADD_eRD[which(tt_control_controlADD_eRD$FDR <= 0.05 & tt_control_controlADD_eRD$FDR > 0.01), dim(tt_control_controlADD_eRD)[2]] <- "med"
tt_control_controlADD_eRD[which(tt_control_controlADD_eRD$FDR <= 0.01), dim(tt_control_controlADD_eRD)[2]] <- "high"
tt_control_controlADD_eRD$candidate <- factor(tt_control_controlADD_eRD$candidate, levels = c("high", "med",  "low", "no"))

tt_control_water_eRD$candidate <- "no"
tt_control_water_eRD[which(tt_control_water_eRD$FDR <= 0.10 & tt_control_water_eRD$FDR > 0.05), dim(tt_control_water_eRD)[2]] <- "low"
tt_control_water_eRD[which(tt_control_water_eRD$FDR <= 0.05 & tt_control_water_eRD$FDR > 0.01), dim(tt_control_water_eRD)[2]] <- "med"
tt_control_water_eRD[which(tt_control_water_eRD$FDR <= 0.01), dim(tt_control_water_eRD)[2]] <- "high"
tt_control_water_eRD$candidate <- factor(tt_control_water_eRD$candidate, levels = c("high", "med",  "low", "no"))

tt_control_waterREMOVE_eRD$candidate <- "no"
tt_control_waterREMOVE_eRD[which(tt_control_waterREMOVE_eRD$FDR <= 0.10 & tt_control_waterREMOVE_eRD$FDR > 0.05), dim(tt_control_waterREMOVE_eRD)[2]] <- "low"
tt_control_waterREMOVE_eRD[which(tt_control_waterREMOVE_eRD$FDR <= 0.05 & tt_control_waterREMOVE_eRD$FDR > 0.01), dim(tt_control_waterREMOVE_eRD)[2]] <- "med"
tt_control_waterREMOVE_eRD[which(tt_control_waterREMOVE_eRD$FDR <= 0.01), dim(tt_control_waterREMOVE_eRD)[2]] <- "high"
tt_control_waterREMOVE_eRD$candidate <- factor(tt_control_waterREMOVE_eRD$candidate, levels = c("high", "med",  "low", "no"))

tt_controlADD_water_eRD$candidate <- "no"
tt_controlADD_water_eRD[which(tt_controlADD_water_eRD$FDR <= 0.10 & tt_controlADD_water_eRD$FDR > 0.05), dim(tt_controlADD_water_eRD)[2]] <- "low"
tt_controlADD_water_eRD[which(tt_controlADD_water_eRD$FDR <= 0.05 & tt_controlADD_water_eRD$FDR > 0.01), dim(tt_controlADD_water_eRD)[2]] <- "med"
tt_controlADD_water_eRD[which(tt_controlADD_water_eRD$FDR <= 0.01), dim(tt_controlADD_water_eRD)[2]] <- "high"
tt_controlADD_water_eRD$candidate <- factor(tt_controlADD_water_eRD$candidate, levels = c("high", "med",  "low", "no"))

tt_controlADD_waterREMOVE_eRD$candidate <- "no"
tt_controlADD_waterREMOVE_eRD[which(tt_controlADD_waterREMOVE_eRD$FDR <= 0.10 & tt_controlADD_waterREMOVE_eRD$FDR > 0.05), dim(tt_controlADD_waterREMOVE_eRD)[2]] <- "low"
tt_controlADD_waterREMOVE_eRD[which(tt_controlADD_waterREMOVE_eRD$FDR <= 0.05 & tt_controlADD_waterREMOVE_eRD$FDR > 0.01), dim(tt_controlADD_waterREMOVE_eRD)[2]] <- "med"
tt_controlADD_waterREMOVE_eRD[which(tt_controlADD_waterREMOVE_eRD$FDR <= 0.01), dim(tt_controlADD_waterREMOVE_eRD)[2]] <- "high"
tt_controlADD_waterREMOVE_eRD$candidate <- factor(tt_controlADD_waterREMOVE_eRD$candidate, levels = c("high", "med",  "low", "no"))

tt_water_waterREMOVE_eRD$candidate <- "no"
tt_water_waterREMOVE_eRD[which(tt_water_waterREMOVE_eRD$FDR <= 0.10 & tt_water_waterREMOVE_eRD$FDR > 0.05), dim(tt_water_waterREMOVE_eRD)[2]] <- "low"
tt_water_waterREMOVE_eRD[which(tt_water_waterREMOVE_eRD$FDR <= 0.05 & tt_water_waterREMOVE_eRD$FDR > 0.01), dim(tt_water_waterREMOVE_eRD)[2]] <- "med"
tt_water_waterREMOVE_eRD[which(tt_water_waterREMOVE_eRD$FDR <= 0.01), dim(tt_water_waterREMOVE_eRD)[2]] <- "high"
tt_water_waterREMOVE_eRD$candidate <- factor(tt_water_waterREMOVE_eRD$candidate, levels = c("high", "med",  "low", "no"))


# what does tt_eRD's look like?
head(tt_control_controlADD_eRD)
head(tt_control_water_eRD)
head(tt_control_waterREMOVE_eRD)
head(tt_controlADD_water_eRD)
head(tt_controlADD_waterREMOVE_eRD)
head(tt_water_waterREMOVE_eRD)
# logFC = log (base 2) fold change
# logCPM = log counts per million

# what does the test p-value distributions look like?
ggplot(tt_control_controlADD_eRD, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_control_controlADD_eRD$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("control_controlADD/edgeR p-value distribution")
# X axis represents 0.01 increments of probability. 
ggplot((tt_control_water_eRD), aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_control_water_eRD$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("control_water/edgeR p-value distribution")
# X axis represents 0.01 increments of probability.
ggplot(tt_control_waterREMOVE_eRD, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_control_waterREMOVE_eRD$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("control_waterREMOVE/edgeR p-value distribution")
# X axis represents 0.01 increments of probability.
ggplot(tt_controlADD_water_eRD, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_cotrolADD_water_eRD$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("cotrolADD_water/edgeR p-value distribution")
# X axis represents 0.01 increments of probability.
ggplot(tt_controlADD_waterREMOVE_eRD, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_controlADD_waterREMOVE_eRD$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("controlADD_waterREMOVE/edgeR p-value distribution")
# X axis represents 0.01 increments of probability.
ggplot(tt_water_waterREMOVE_eRD, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_water_waterREMOVE_eRD$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("water_waterREMOVE/edgeR p-value distribution")
# X axis represents 0.01 increments of probability.

# function for MA plots
pw_ma_plot <- function(frame, x, y, f, title) {
  # frame = data frame with data
  # x, y are the string names of the x and y columns
  # f is the factor for faceting, title is a string for plot titles
  # make the main MA plot
  temp <- data.frame(log2((frame[x] + frame[y])/2), log2(frame[y] / frame[x]), frame[f])
  colnames(temp) <- c("Ave", "FC", "candidate")
  first  <- ggplot(temp, aes(x = Ave, y = FC)) +
    geom_point(aes(color = candidate, shape = candidate)) +
    scale_y_continuous(paste0("logFC (", x, "/", y, ")")) +
    scale_x_continuous("Ave_intensity") +
    ggtitle(title) + 
    geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
    geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") # 2-fold down
  
  # make separate MA plots
  second <- ggplot(temp, aes(x = Ave, y = FC)) +
    geom_point(aes(color = candidate, shape = candidate)) +
    scale_y_continuous(paste0("logFC (", x, "/", y, ")")) +
    scale_x_continuous("Ave_intensity") +
    geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
    geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") + # 2-fold down
    facet_wrap(~ candidate) +
    ggtitle(paste(title, "(separated)", sep=" "))
  # plots inside functions do not automatically display
  print(first)
  print(second)
}

pw_scatter_plot <- function(frame, X, Y, f, title) {
  # frame = data frame with data
  # X, Y are the string names of the x and y columns
  # f is the factor for faceting, title is a string for plot titles
  # make the combined candidate corelation plot
  first <- ggplot(frame, aes_string(X, Y)) +
    geom_point(aes_string(color = f, shape = f)) +
    scale_y_log10() +
    scale_x_log10() +
    ggtitle(title) + 
    geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
    geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") # 2-fold down
  
  # make separate corelation plots
  second <- ggplot(frame, aes_string(X, Y)) +
    geom_point(aes_string(color = f, shape = f)) +
    scale_y_log10() +
    scale_x_log10() +
    geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
    geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold down
    facet_wrap(~ candidate) +
    ggtitle(paste(title, "(separated)", sep=" ")) 
  
  print(first)
  print(second)
}

# for plotting results, we will use the average intensities for the samples
ave_eRD$control_controlADD_candidate <- tt_control_controlADD_eRD$candidate
ave_eRD$control_water_candidate <- tt_control_water_eRD$candidate
ave_eRD$control_waterREMOVE_candidate <- tt_control_waterREMOVE_eRD$candidate
ave_eRD$controlADD_water_candidate <- tt_controlADD_water_eRD$candidate
ave_eRD$controlADD_waterREMOVE_candidate <- tt_controlADD_waterREMOVE_eRD$candidate
ave_eRD$water_waterREMOVE_candidate <- tt_water_waterREMOVE_eRD$candidate

volcano_control_controlADD_eRD <- data.frame(log2(ave_eRD$controlADD / ave_eRD$control), log10(tt_control_controlADD_eRD$FDR)*(-1), ave_eRD$control_controlADD_candidate)
volcano_control_water_eRD <- data.frame(log2(ave_eRD$water / ave_eRD$control), log10(tt_control_water_eRD$FDR)*(-1), ave_eRD$control_water_candidate)
volcano_control_waterREMOVE_eRD <- data.frame(log2(ave_eRD$waterREMOVE / ave_eRD$control), log10(tt_control_waterREMOVE_eRD$FDR)*(-1), ave_eRD$control_waterREMOVE_candidate)
volcano_controlADD_water_eRD <- data.frame(log2(ave_eRD$water / ave_eRD$controlADD), log10(tt_controlADD_water_eRD$FDR)*(-1), ave_eRD$controlADD_water_candidate)
volcano_controlADD_waterREMOVE_eRD <- data.frame(log2(ave_eRD$waterREMOVE / ave_eRD$controlADD), log10(tt_controlADD_waterREMOVE_eRD$FDR)*(-1), ave_eRD$controlADD_waterREMOVE_candidate)
volcano_water_waterREMOVE_eRD <- data.frame(log2(ave_eRD$waterREMOVE / ave_eRD$water), log10(tt_water_waterREMOVE_eRD$FDR)*(-1), ave_eRD$water_waterREMOVE_candidate)


colnames(volcano_control_controlADD_eRD) <- c("FoldChange", "FDR", "candidate")
colnames(volcano_control_water_eRD) <- c("FoldChange", "FDR", "candidate")
colnames(volcano_control_waterREMOVE_eRD) <- c("FoldChange", "FDR", "candidate")
colnames(volcano_controlADD_water_eRD) <- c("FoldChange", "FDR", "candidate")
colnames(volcano_controlADD_waterREMOVE_eRD) <- c("FoldChange", "FDR", "candidate")
colnames(volcano_water_waterREMOVE_eRD) <- c("FoldChange", "FDR", "candidate")

head(volcano_control_controlADD_eRD)
head(volcano_control_water_eRD)
head(volcano_control_waterREMOVE_eRD)
head(volcano_controlADD_water_eRD)
head(volcano_controlADD_waterREMOVE_eRD)
head(volcano_water_waterREMOVE_eRD)


# start with MA plot
pw_ma_plot(ave_eRD, "control", "controlADD", "control_controlADD_candidate", "Control v ControlADD")
pw_ma_plot(ave_eRD, "control", "water", "control_water_candidate", "Control v Water")
pw_ma_plot(ave_eRD, "control", "waterREMOVE", "control_waterREMOVE_candidate", "Control v Water Remove")
pw_ma_plot(ave_eRD, "controlADD", "water", "controlADD_water_candidate", "ControlADD v Water")
pw_ma_plot(ave_eRD, "controlADD", "waterREMOVE", "controlADD_waterREMOVE_candidate", "ControlADD v Water Remove")
pw_ma_plot(ave_eRD, "water", "waterREMOVE", "water_waterREMOVE_candidate", "Water v Water Remove")


# now the scatter plot 
## Well this is broken...
## Ok, turns out that the column name HAS to be candidate... the MA plot was fine though.. 
## Go to the function and see if there's specificity there?
ave_eRD$candidate <- ave_eRD$control_controlADD_candidate ## Hmmm workaround? Can rename plot then remove...
pw_scatter_plot(ave_eRD, "control", "controlADD", "control_controlADD_candidate", "control_controlADD_candidate")
ave_eRD <- ave_eRD[,-11]
ave_eRD$candidate <- ave_eRD$control_water_candidate ## Hmmm workaround? Can rename plot then remove...
pw_scatter_plot(ave_eRD, "control", "water", "control_water_candidate", "control_water_candidate")
ave_eRD <- ave_eRD[,-11]
ave_eRD$candidate <- ave_eRD$control_waterREMOVE_candidate ## Hmmm workaround? Can rename plot then remove...
pw_scatter_plot(ave_eRD, "control", "waterREMOVE", "control_waterREMOVE_candidate", "control_waterREMOVE_candidate")
ave_eRD <- ave_eRD[,-11]
ave_eRD$candidate <- ave_eRD$controlADD_water_candidate ## Hmmm workaround? Can rename plot then remove...
pw_scatter_plot(ave_eRD, "controlADD", "water", "controlADD_water_candidate", "controlADD_water_candidate")
ave_eRD <- ave_eRD[,-11]
ave_eRD$candidate <- ave_eRD$controlADD_waterREMOVE_candidate ## Hmmm workaround? Can rename plot then remove...
pw_scatter_plot(ave_eRD, "controlADD", "waterREMOVE", "controlADD_waterREMOVE_candidate", "controlADD_waterREMOVE_candidate")
ave_eRD <- ave_eRD[,-11]
ave_eRD$candidate <- ave_eRD$water_waterREMOVE_candidate ## Hmmm workaround? Can rename plot then remove...
pw_scatter_plot(ave_eRD, "water", "waterREMOVE", "water_waterREMOVE_candidate", "water_waterREMOVE_candidate")
ave_eRD <- ave_eRD[,-11]

# make a volcano plot ////////////////////////////////////////////////////
ggplot(volcano_control_controlADD_eRD, aes(x = FoldChange, y = FDR)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ylim(c(NA, 0.5)) +
  xlim(c(-1.25, 1.25)) +
  ggtitle("control_controlADD Volcano Plot")

# see how many metabolites are in each DE category
summary(ave_eRD$control_controlADD_candidate)

# collect up some of the data frames and write to disk
final_MQ_frame <- cbind(anno_MQ, data_eRD_tmm, tt_eRD)
write.csv(final_MQ_frame, file = "final_MQ_frame.csv", row.names = FALSE)
# |Important|:
# Be sure to check the identified candidate accessions to the original MQ proteinGroups.txt file. 
# Some Accessions may have been grouped with others which were removed during initial data processing.
# Note, this is due to the peptide assignment being equally likely to belong to either protein. 
# Can help to understand the biological problem, same proteins in a pathway, instead of completely unrelated processes.

## Should I just leave all the accessions as they are? May break the row.names arguement
## Remember you can autofilter the output to pull out the significant hits 


# compare the edgeR testing to a more basic t-test
# do the t-test on log transformed intensities to be safe 
ttest_MQ <- log2(data_MQ_tmm)
# add average ratio columns (non-logged ratios), fold-change column, and row names
ttest_MQ$ave_control <- rowMeans(data_MQ_tmm[1:5])
ttest_MQ$ave_RU  <- rowMeans(data_MQ_tmm[6:10])
ttest_MQ$logFC <- log2(ttest_MQ$ave_RU / ttest_MQ$ave_control)
row.names(ttest_MQ) <- anno_MQ$Accession

# apply the basic two-sample t-test (we will pool variance)
t.result <- apply(ttest_MQ, 1, function(x) t.test(x[1:5], x[6:10], var.equal = TRUE))
# extract the p-value column from the t-test thingy 
ttest_MQ$p_value <- unlist(lapply(t.result, function(x) x$p.value))
# do a Benjamini-Hochberg multiple testing correction
ttest_MQ$fdr <- p.adjust(ttest_MQ$p_value, method = "BH")

# add a DE candidate status column
ttest_MQ$candidate <- "no"
ttest_MQ[which(ttest_MQ$fdr <= 0.10 & ttest_MQ$fdr > 0.05), dim(ttest_MQ)[2]] <- "low"
ttest_MQ[which(ttest_MQ$fdr <= 0.05 & ttest_MQ$fdr > 0.01), dim(ttest_MQ)[2]] <- "med"
ttest_MQ[which(ttest_MQ$fdr <= 0.01), dim(ttest_MQ)[2]] <- "high"
ttest_MQ$candidate <- factor(ttest_MQ$candidate, levels = c("high", "med",  "low", "no"))
head(ttest_MQ)

# count up, down and the rest (FDR less than 0.05)
all <- dim(ttest_MQ)[1]
up <- dim(ttest_MQ[(ttest_MQ$fdr <= 0.05) & (ttest_MQ$logFC > 0.0), ])[1]
down <- dim(ttest_MQ[(ttest_MQ$fdr <= 0.05) & (ttest_MQ$logFC <= 0.0), ])[1]
print("This is like the decideTest in edgeR - 5% FDR cut:")
up 
all - up - down
down
print("Candidate Counts:")
summary(ttest_MQ$candidate)
# Hmmmm... using a t-test candidates fall by the way... 
# Only significant in the edgeR analysis? 
# Anyway lets carry on for now.. 

# what does the test p-value distribution look like?  ## Why is the count arguement there? Ask Mirre :/
ggplot(ttest_MQ, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_MQ$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("MQ data with t-test p-value distribution")

# start with MA plot
pw_ma_plot(ttest_MQ, "ave_control", "ave_RU", "candidate", "MQ t-test data")
# now the scatter plot
pw_scatter_plot(ttest_MQ, "ave_control", "ave_RU", "candidate", "MQ t-test data")

# and the volcano plot
volcano_tt <- data.frame(log2(ttest_MQ$ave_RU / ttest_MQ$ave_control), log10(ttest_MQ$fdr)*(-1), ttest_MQ$candidate)
colnames(volcano_tt) <- c("FoldChange", "FDR", "candidate")
ggplot(volcano_tt, aes(x = FoldChange, y = FDR)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ylim(c(NA, 0.25)) +
  ggtitle("MQ t-test Volcano Plot")












### Next, look into comparing groups / utilising the intensity data. 
### Correlation analysis etc - more functionality in MetabanalystR.
### Maybe take a look into the ?MALDIquant package, implement some of the graphical outputs?













