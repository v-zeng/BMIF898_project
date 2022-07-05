# load packages
library(data.table)
library(dplyr)
# read in data
df <- fread('KIRC_tumor_samples.csv')
df[1:5,1:5]
dim(df) # 485577    483
# Read in keep probes from MacIntyre and preview
keepProbes <- read.csv("FinalMcIntyreKeepProbesAndInfo.csv")
probesClean <- df %>% # normalsPrepared previously 'probes'
  filter(df[[1]] %in% keepProbes$ID) # index column in data frame with double square bracket
colnames(probesClean)[1] <- "ProbeID"
probesClean <- data.frame(probesClean, row.names = 1)
dim(probesClean) # 284482    484
#change first column to "ProbeID" =============================================
# # check missingness
# probesClean[1:5,1:5]
# dim(probesClean)
# probesClean <- t(probesClean)
# # which probes missing all values
# samples<-nrow(probesClean)
# length(which(colSums(is.na(probesClean))==samples))

missingcells = sum(is.na(probesClean))
print("Missing value cells")
print(missingcells)
dimProbes<-dim(probesClean)
totalcells <- dimProbes[1]*dimProbes[2]
# calculating percentage of missing values
percentage = (missingcells * 100 )/(totalcells)
print("Percentage of missing values' cells")
print (percentage) #5.081865

# MICE imputation==============================================================
library('mice')
set.seed(3)
# run multiple (m=5) imputation
imputed <- mice(probesClean,meth="pmm",m=5,maxit=20)
# create dataset after imputation
imputed <- complete(imputed)
write.csv(imputed, "KIRC_tumors_imputed.csv")
###
# omit df rows with NA
# df_noNA <- na.omit(probesClean)
# dim(df_noNA) # 262151    484
# df_noNA[1:5,1:5]
# # use first column as rownames
# df_noNA <- data.frame(df_noNA, row.names = 1)
# df_noNA[1:5,1:5]
###

# save as ____tumor_noNA_keepProbes.csv
#write.csv(df_noNA,"KIRCtumorsProcessedZeroNA.csv")

# merge clinical and tumor samples processed ===================================
# library(data.table)
# library(dplyr)
df <- fread("KIRC_tumors_imputed.csv")
sum(is.na(df))
dim(df)
df[1:5,1:5]
dfSampleOnly_data <- df[,-1] # drop 'ProbeID'
# replace dots '.' with hyphens '-' in row names
names(dfSampleOnly_data) <- gsub('[.]','-',names(dfSampleOnly_data))
dfSampleOnly_data[1:5,1:5]
# substring names of barcodes to only include first 12 characters
names(dfSampleOnly_data) <- substring(names(dfSampleOnly_data),1,12)
names(df) <- substring(names(df),1,12)
# barcode list for clinical data row filtering
barcodeList <- colnames(dfSampleOnly_data)


# read in clinical data
dfClinical <- fread('cBioPortalclinicalKIRC.csv') # contains columns "patientID", "AGE"
# keep rows that match barcodeList from project; can replace with subset()===========
dfClinicalClean <- dfClinical %>% # 
  filter(dfClinical$patientId %in% barcodeList)
dfClinicalClean <- dfClinicalClean[!duplicated(dfClinicalClean$patientId)] # remove any duplicate elements for patientID
#===============================================================================

# select 'patientId' and 'AGE'
dfAgeID <- dfClinicalClean[,c('patientId','AGE')]
names(dfAgeID)[c(1,2)] <- c('PatientID','Age') # format column names appropriately ==== may not need this line

# rename dfAgeID rows with 'ProbeID'
dfAgeID <- data.frame(dfAgeID, row.names = 1) # column '1' used for row names, column is removed automatically

# rename dfSample_data rows with 'ProbeID' column
df<- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
df[1:5,1:5]

# transpose and rename==========================================================
df_t <- t(df)
df_t[1:5,1:5]

# replace dots '.' with hyphens '-' in row names
rownames(df_t) <- gsub('[.]','-',rownames(df_t))
df_t[1:5,1:5]
# merge dfAgeID and df_t========================================================
df_merged <- merge(dfAgeID,df_t,by = 0) # merge by rows, but adds 'row.names' (see Rdocumentation for merge()) - try to remove in one line
df_merged[1:5,1:5]
# rename dfSample_data rows with 'ProbeID' column
df_mergedClean <- data.frame(df_merged, row.names = 1) # column '1' used for row names, column is removed automatically
dim(df_mergedClean) # 275 214
# omit rows with NA (some ages not available)
df_mergedClean <- na.omit(df_mergedClean)
dim(df_mergedClean) # 272 214
df_mergedClean[1:5,1:5]
# csv of above
write.csv(df_mergedClean, 'KIRCtumors_probesAge_merged.csv')
