# create df with age and probes as columns ================================
library(data.table)
library(dplyr)
# need add comma if using rest of list
projectList <- list(#"LUAD",
                    "COAD" 
                    #"THCA",
                    #"LUSC",
                    #"KIRP",
                    # "UCEC",
                    # "LIHC",
                    # "PRAD",
                    # "HNSC"
                    # "BRCA"
                    #"KIRC"
) 
#===========================================
for (project in projectList){

  df <- fread(paste0(project,'NormsCommonProcessedZeroNA.csv')) #updated to '...commonProcessed...' DEC 19
  df <- subset(df,select=-(SD)) # drop last column, which is the SD column we added; or use 'allNormals_TCGABLCA_fall2021.csv' without this line
  dfSample_data <- df[,-c(2:10)] # drop columns which are not 'ProbeID' and sample barcodes
  dfSampleOnly_data <- dfSample_data[,-1] # drop 'ProbeID'
  
  # substring names of barcodes to only include first 12 characters
  names(dfSample_data) <- substring(names(dfSample_data),1,12)
  
  # barcode list for clinical data row filtering
  barcodeList <- colnames(dfSample_data)[-1] # drop "ProbeID" ('V1' in cBioPortal data?) column
  
  
  # read in clinical data
  dfClinical <- fread(paste0('cBioPortalclinical',project,'.csv')) # contains columns "patientID", "AGE"
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
  dfSample_data <- data.frame(dfSample_data, row.names = 1) # column '1' used for row names, column is removed automatically
  
  
  # transpose and rename==========================================================
  df_t <- transpose(dfSample_data)
  
  #redefine column and row names
  colnames(df_t) <- rownames(dfSample_data)
  rownames(df_t) <- colnames(dfSample_data)
  
  # replace dots '.' with hyphens '-' in row names
  rownames(df_t) <- gsub('[.]','-',rownames(df_t))
  
  
  # merge dfAgeID and df_t========================================================
  df_merged <- merge(dfAgeID,df_t,by = 0) # merge by rows, but adds 'row.names' (see Rdocumentation for merge()) - try to remove in one line
  
  # rename dfSample_data rows with 'ProbeID' column
  df_mergedClean <- data.frame(df_merged, row.names = 1) # column '1' used for row names, column is removed automatically
  # dim(df_mergedClean) # 160 267462
  
  # csv of above
  write.csv(df_mergedClean, paste0(project,'_probesAge_mergedCommon.csv'))
}
