# start fitting each project on chronological clock
# start with KIRC (most samples)
# clinical data - indexed or XML?
library(data.table)
library(cBioPortalData)
library(AnVIL)
library(rapiclient)
# create loop through normals to download clinical data
# files <-list.files(pattern = "NormsProcessedZeroNA")
projectList <- list("COAD", # coadread_tcga ????
                    # "LUAD", # current script test for 2 projects================
                     # need add comma if using rest of list
                    # "THCA",
                    # "LUSC",
                    # "KIRP",
                    "UCEC",
                    "LIHC",
                    "PRAD",
                    "HNSC"
                    # "BRCA"
                    #"KIRC"
                    ) 
for (project in projectList){

  # df <- fread(paste0(project,'NormsProcessedZeroNA.csv'))
  # df <- subset(df,select=-(SD)) # drop last column, which is the SD column we added; or use 'allNormals_TCGABLCA_fall2021.csv' without this line
  # df <- df[,-c(2:10)] # drop columns which are not 'ProbeID' and sample barcodes
  # # substring names of barcodes to only include first 12 characters
  # names(df) <- substring(names(df),1,12)
  # # check first five column names
  # head(names(df))
  
  # create list of barcodes for GDCquery clinical data retrieval, which does not include first column ('ProbeID')
  # barcodeList <- colnames(df)[-1]

  # cBioPortalData API ===========================================================
  # check documentation again for arguments
  
  (cbio <- cBioPortal())
  ClinicalData <- clinicalData(cbio, (paste0(tolower(project),"_tcga")))
  write.csv(ClinicalData, paste0("cBioPortalclinical",project,".csv"))
}
