##### **************** GET GSE52401DNA methyl beta values  ***************  #####

# if (!require("BiocManager", quietly = TRUE))
#   +     install.packages("BiocManager")
# BiocManager::install()
# # R 4.1.2  Bioconductor 3.14
# BiocManager::install("GEOquery")


# GEO 450k data to download==========================================================
# GSE105260: kidney n=9 normal; * no age data
# GSE79100: kidney 31 normal * Error in eval(predvars, data, env) : object 'cg01330954' not found; sum(is.na(mdat_t$cg01330954)) = 1
# GSE88883: breast n=100 * Error in eval(predvars, data, env) : object 'cg15394630' not found * missing all values for cg15394630
# GSE37754: breast n=10 (Non-tumor) * study identified DNAm patterns associated with tumor subtypes
# GSE66313: breast n=15 (Adjacent-Normal) * missing clock probe
# GSE69914: breast n=50 (incomplete meta data warning on GEO page) * no age
# GSE101961: breast n=121 normal tissue in cancer-free women ===================
# GSE61441: kidney n=46 normal * no age data; emailed
# GSE63704: lung n=43 * no age data; emailed - replied, no ages available
# GSE50874: kidney * no age, not relevant
# GSE105288: kidney but GPL10558 and GPL13534 * parsing issue
# GSE113501: Kidney * no age for tumor free samples * cg15394630 missing all
# GSE85566: lung ; asthmatic and non-asthmatic * object 'cg07856699' not found * missing cg14952449 probe all values
# GSE74214: n=18 non-diseased breast tissue; *object 'cg14952449' not found
# GSE66836: n=19 normal, n=164 lung adenocarcinoma * no age
# GSE124367: n=12 normal human breast tissue samples * no age
# GSE56044: lung 124 tumors, 12 normal * no ages
# GSE52401: lung n=244 * no ages
# GSE39279: lung n=25 normal, n=339 non-small cell lung cancer *Error in eval(predvars, data, env) : object 'cg14952449' not found
# GSE78754: breast n=11 normal adjacent, n=12 lymph node metastases * unable to find "normal adjacent" samples
# GSE101443: breast n=4 normal, 4 tumor *no age
# GSE126441: kidney n=10 normal * no age
# GSE156932: kidney n=15, benign oncocytoma vs renal cell carcinoma
# GSE70303: n=6 normal/tumor/treatments(+/-), clear cell renal cell carcinoma * no age
# GSE97466: thyroid n=50 non-neoplastic adjacent tissues * Error in eval(predvars, data, env) : object 'cg15394630' not found *** probe not in data
# GSE86961: thyroid n=41 non-neoplastic adjacent tissues * Error in eval(predvars, data, env) : object 'cg15394630' not found *** probe not in data

# GSE83842: lung adenocarcinoma n=24 * too few samples? * Error in eval(predvars, data, env) : object 'cg05533001' not found
# GSE94785: asbestos-associated DNAm changes in lung cancer n=28 fresh frozen; too few samples: non-asbestos exposed normal=14, tumor=14
# GSE67919: tumor adjacent normal tissue from BRCA patients * Error in eval(predvars, data, env) : object 'cg14952449' not found

# GSE31848: pluripotent stem cells... maybe good to test on this to see if overlap... should close to age 0? low DNAm?

###
# GSE62336: head and neck: nasopharyngeal carcinoma maybe check this?...
# GSE158075: contralateral bronchus (controls)... check this if nothing else... * only 1 specimen...
# GSE56588: liver n=10 normal, n=9 Cirrhotic * test clock vs allFeatures in detecting liver
# GSE50498: skeletal muscle
# GSE144858: blood, Alzheimer's
# GSE70977: oral rinse samples 154 cases, 72 controls * Error in eval(predvars, data, env) : object 'cg04606672' not found
# GSE94876 buccal cells of long-term smokers and moist snuff consumers 

### Hannum data: GSE40279
### GSE67751: blood, HIV +/-

library(GEOquery)
library(Biobase)
library(data.table)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 6)


gdat <- getGEO('GSE59157', 
               GSEMatrix=TRUE,
               getGPL=FALSE, # timeout when TRUE,
               destdir="C:/Users/vzeng/Documents/R/BMIF Scripts/tcga_normals/normals_zeroNA")

gdat
# $GSE113501_series_matrix.txt.gz
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 485577 features, 144 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM3107981 GSM3107982 ... GSM3108124 (144 total)
# varLabels: title geo_accession ... tnm stage:ch1 (40 total)
# varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
# pubMedIds: 30642274 
# Annotation: GPL13534 


gse <- gdat$GSE59157_series_matrix.txt.gz
mdat <- exprs(gse)
dim(mdat)
# 418811    223
mdat[1:5,1:5]
#             GSM1429528 GSM1429529 GSM1429530 GSM1429531 GSM1429532
# cg00000029  0.3218820  0.4208604  0.4009487  0.5953418  0.4294661
# cg00000108  0.9122912  0.8745204  0.8726442  0.8684039  0.8274986
# cg00000109  0.7633889  0.8409980  0.7808807  0.8246498  0.8225860
# cg00000165  0.1938463  0.1708587  0.2723594  0.3313581  0.2036975
# cg00000236  0.7247475  0.7364676  0.6582602  0.6602951  0.6673423

# check for NA
any(is.na(mdat)) # GSE40279 hannum no NA

# transpose and count missing values in cg15394630
# check if cg15394630 row exists
any(row.names(mdat) == 'cg15394630') # cg15394630 does not exist in GSE97466
mdat_t <- as.data.frame(t(mdat))
mdat_t[1:5,1:5]
mdat_t$cg01330954[1:5]
sum(is.na(mdat_t$cg01330954))
nrow(mdat_t)
# # omit rows (probes) with NA values; if transposed then samples with any NA values will be removed
# mdat_t <- na.omit(mdat_t)
# dim(mdat_t)# 29 482954
# mdat <- t(mdat_t)
# # keep rows (probes) that contain less than 10 missing values
# mdat <- mdat[rowSums(is.na(mdat)) < 10, ]
# dim(mdat) # 485512     18
# write.csv(mdat,"GSE52955_KIDNEY_exprs.csv") # need to filter common probes========

# clinical/phenotype data=======================================================
clindat <- pData(gse)
dim(clindat)
# 183  43
colnames(clindat)
# [1] "title"                   "geo_accession"           "status"                  "submission_date"        
# [5] "last_update_date"        "type"                    "channel_count"           "source_name_ch1"        
# [9] "organism_ch1"            "characteristics_ch1"     "characteristics_ch1.1"   "characteristics_ch1.2"  
# [13] "characteristics_ch1.3"   "characteristics_ch1.4"   "characteristics_ch1.5"   "molecule_ch1"           
# [17] "extract_protocol_ch1"    "label_ch1"               "label_protocol_ch1"      "taxid_ch1"              
# [21] "hyb_protocol"            "scan_protocol"           "description"             "data_processing"        
# [25] "platform_id"             "contact_name"            "contact_email"           "contact_institute"      
# [29] "contact_address"         "contact_city"            "contact_state"           "contact_zip/postal_code"
# [33] "contact_country"         "supplementary_file"      "data_row_count"          "age (y):ch1"            
# [37] "ethnicity:ch1"           "gender:ch1"              "plate:ch1"               "source:ch1"             
# [41] "tissue:ch1"  
sink("clindat.txt")
length(clindat$)
clindat[1,] # [46] is profile of sample (e.g., normal, tumor), [39] is age in months
sink()
# find unique values in column of interest
unique(clindat$`age (y):ch1`)
table(clindat$`histology:ch1`)
# 51 57 59 60 65 68 73 74 75 77 80 
# 2  2  2  2  2  4  2  2  2  2  2 
unique(clindat[[62]])
hist(unlist(clindat$`age:ch1`))
clindat_Age <- as.numeric(clindat$`age:ch1`)
hist(clindat_Age,breaks=10)

sink("clindat_file.txt")

sink()
clindat$`age (months):ch1` # remove rows with "unknown"; some ages < 20... need age transformation function?
length(clindat$`age (months):ch1`) # 95
table(clindat$`tissue:ch1`)
# 18 20 22 23 30 31 33 39 47 51 53 54 56 60 63 64 65 70 71 75 78 
# 1  1  1  1  1  1  1  1  1  1  3  1  2  4  1  2  1  2  2  1  2 

# omit rows with NA in 'age (months):ch1'=======================================
clindat <- clindat[!grepl("unknown",clindat$`asbestos exposed:ch1`),] # use ! to omit
# check after omitting "unknown" - which set is this for, one of the breast?
table(clindat$age)
# 10  11  12 144  16  20  21  23  28  29  31  33  38  39  41  42  44  47  49  54  58  62  66  71 
# 8   2   6   2   3   7   4   3   2   5   3   3   2   5   2   9   4   6   5   2   2   3   1   3 

length(clindat$`age:ch1`) # 15
clindat[1:5,c(43,54)] # age is in months, convert to years
clindatAge <- clindat[,41, drop = F]#======================
table(clindat$age)
# 13 15 17 23 36 44 48 51 58 59 68 69 70 80 
# 1  1  1  1  1  2  1  1  2  2  1  1  2  1 
dim(clindatAge) # 10  1
names(clindatAge)[1] <- "Age" # ==================
head(clindatAge)


str(clindatAge) # character type data
# convert to numeric
clindat$`age (months):ch1` <- as.numeric(as.character(clindat$`age (months):ch1`))
# rename `age (months):ch1` column to Age
names(clindat)[37] <- "Age"
# rename 'tissue:ch1' to Sample_type
names(clindat)[39] <- "Sample_type"
# convert age in months to years
clindat[39] <- clindat[39]/12
clindatAge_status <- clindat[,c(36,41)] # age:ch1 and tissue:ch1
clindatAge_status[1:5,1:2]
names(clindatAge_status)[1] <- "Age"
names(clindatAge_status)[2] <- "Tissue"
clindatAge_status[1:5,1:2]

# write.csv(clindatAge_status,"GSE94785_lung_clindatAge.csv")
# dim(mdat) 390292    100

# filter TCGA common probes=====================================================
df_commonProbes <- fread('commonProbesZeroNA.csv')
# df <- fread()

# #create list of data
# my_data <- lapply(files, fread) # apply function over list or vector

#for each data frame, subset common probes in 'commonProbesZeroNA.csv'
keepProbes <- df_commonProbes$ProbeID
commonKeep <- subset(mdat, rownames(mdat) %in% keepProbes)

# transpose and rename==========================================================
df_t <- t(commonKeep)

# both methods below do not save row names properly
# merge dfAgeID and df_t========================================================
# df_merged <- merge(clindatAge,df_t,by = 0) # merge by rows, but adds 'row.names' (see Rdocumentation for merge()) - try to remove in one line
# # rename dfSample_data rows with 'ProbeID' column
# df_mergedClean <- data.frame(df_merged, row.names = 1) # column '1' used for row names, column is removed automatically
# dim(df_mergedClean) # 15 217032
# try transform to merge and remove row.names in one 'step'
df_mergedTest <- transform(merge(clindatAge_status,df_t,by=0),row.names=Row.names,Row.names=NULL) # https://stackoverflow.com/questions/17375849/how-does-one-merge-dataframes-by-row-name-without-adding-a-row-names-column
df_mergedTest[1:5,1:5]
dim(df_mergedTest) # 144 227919
# for GSE79100, drop "Tissue" - i.e., "normal kidney"
# df_mergedTest <- df_mergedTest[,-2]
# df_mergedTest[1:5,1:5]
# dim(df_mergedTest) # 100 231023 <--- over 20K probes omitted due to missing values...
# df_mergedTest[1:5,1:5]
# dim(df_mergedTest)
# # csv of above
write.csv(df_mergedTest, 'hannum_merged.csv')
dfnew <- fread('hannum_merged.csv')
dfnew[1:5,1:5]
dim(dfnew) # 18 256744
# # subset specific sample_type (e.g., only normal)
# df_normal <- df[grepl('normal', df$sample_type),]

# Is there a generally accepted method used to impute missing DNA methylation data  (probe beta-values)? First time dealing with this situation, would appreciate any suggestions!
# #   
# I'm looking at DNA methylation (DNAm) data such as TCGA (e.g., BRCA, KIRP, KIRC, etc.). Currently trying to build use my model to predict DNAm age on test sets, but many of the data sets  are missing key probe values used in my clock/model.
# I have inspected approximately 30 GEO data series and the following TCGA cancer samples: KIRP, KIRC, LUAD, LUSC, BRCA, THCA. They are each missing key probes I am using.
# Is there a common way to impute the missing values in R?
