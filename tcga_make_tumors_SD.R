# download tcga tumors
# subset McIntyre probes and no NA probes, sort both by increasing SD

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(data.table)

# list of projects =============================================================
projectList <- list(#"LUAD",
  #"COAD",
  #"THCA",
  #"LUSC",
  #"KIRP",
  #"UCEC",
  #"LIHC",
  #"PRAD",
  #"HNSC",
  #"BRCA",
  "KIRC"
) 

# function: create project script inputs =======================================
makeScriptInputs <- function(project) {
  projectName <- paste0("TCGA-",project)
  tumorsFileName <- paste0("allTumors_TCGA_",project,"_fall2021.csv")
  processedFileNameKeep <- paste0(project,"tumsProcessedKeep.csv")
  processedFileNameZeroNA <- paste0(project,"tumsProcessedZeroNA.csv")
  
  scriptInputs <- list(
    projectName, 
    tumorsFileName, 
    processedFileNameKeep,
    processedFileNameZeroNA
  )
  return(scriptInputs)
}

# loop through projects to create SD ===========================================
for (project in projectList) {
  scriptInputs <- makeScriptInputs(project)
  
  # GDCquery
  queryMethylTumors <- GDCquery(project = scriptInputs[1], # input project
                                 legacy=F,
                                 data.category = "DNA Methylation",
                                 platform = c("Illumina Human Methylation 450"),
                                 sample.type = "Primary Tumor")
  
  
  GDCdownload(queryMethylTumors)
  
  tumorsPrepared <- GDCprepare(queryMethylTumors, summarizedExperiment=F)
  ###write.csv(probesClean_SD_sorted, 'keepProbes_SamplesClean.csv', row.names=F) # row.names=F prevents index addition
  write.csv(tumorsPrepared, scriptInputs[[2]]) # input tumorsFileName
  
  # read project probe data=====================================================
  
  
  tumors <- fread(scriptInputs[[2]],header=T) # tumorsFileName === redundant? no header for tumorsPrepared
  
  
  # Read in keep probes from MacIntyre and preview
  keepProbes <- read.csv("FinalMcIntyreKeepProbesAndInfo.csv")
  
  # cat("Dimensions of keepProbes data:", dim(keepProbes)) # 284482 x 24
  
  
  # probe filtering ============================================================
  
  
  # filter keep probes in tums; also remDLBCes 'chrX', 'chrY', and non-cg identifiers
  probesClean <- tumors %>% # tumorsPrepared previously 'probes'
    filter(tumors[[1]] %in% keepProbes$ID) # index column in data frame with double square bracket
  # probes[,1], previously probes[[1]]
  
  # # check number of probes in clean set
  # numProbesClean <- nrow(probesClean) # 284482 probes
  
  
  # change first column to "ProbeID" =============================================
  
  
  colnames(probesClean)[1] <- "ProbeID"
  
  
  # NA values ==================================================================
  
  
  # determine number of normal samples ==== terrible coding======================
  # numSamples <- probesClean[,11:] # subset samples only from data (beta values)
  # numSamples <- ncol(numSamples) # number of samples by counting columns in subset
  
  
  # cat("Number of samples in this data set:", numSamples) # check number of samples
  
  # remove rows with any NA values
  probesZeroNA <- na.omit(probesClean) # omits rows with NA values
  probesZeroNA_count <- nrow(probesZeroNA) # probes without any missing values; numProbesZeroNA
  # cat("Percentage of probes remaining from clean set:", probesZeroNA_count/numProbesClean*100) # 94.82744 percent of probes from 'clean' set
  
  
  # create SD column ===========================================================
  
  
  # # sorted ascending SD for keep probes (McIntyre probes)
  # probesClean_SD <- transform(probesClean, SD=apply(probesClean[,-(1:10)],1,sd)) # only apply to beta values
  # probesClean_SD_sorted <- probesClean_SD[order(probesClean_SD$SD),] # sort rows by SD (default for decreasing = FALSE)
  # # cat("Dimensions of clean probes sorted by SD:", dim(probesClean_SD_sorted))
  
  # sorted ascending SD for keep probes with zero NA 
  probesZeroNA_SD <- transform(probesZeroNA, SD=apply(probesZeroNA[,-(1:10)],1,sd)) # only apply to beta values
  probesZeroNA_SD_sorted <- probesZeroNA_SD[order(probesZeroNA_SD$SD),] # sort rows by SD (default for decreasing = FALSE)
  # cat("Dimensions of zero NA probes sorted by SD:", dim(probesZeroNA_SD_sorted))
  
  
  # write required data frames to .csv files ===================================
  
  
  # keepProbes with SD
  write.csv(probesClean_SD_sorted, scriptInputs[[3]], row.names=F) # row.names=F prevents index addition
  
  
  # probes zero NA with SD
  write.csv(probesZeroNA_SD_sorted, scriptInputs[[4]], row.names=F)
}
