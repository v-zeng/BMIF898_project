# append data frames together: https://stackoverflow.com/questions/29402528/append-data-frames-together-in-a-for-loop/29419402
library(data.table)


# dfKidney <- fread('Kidney_probesAge_mergedCommon.csv')
# dfLung <- fread('Lung_probesAge_mergedCommon.csv')
dfTHCA <- fread('THCA_probesAge_mergedCommon.csv')
dfBRCA <- fread('BRCA_probesAge_mergedCommon.csv')
dfLUSC <- fread('LUSC_probesAge_mergedCommon.csv')
dfLUAD <- fread('LUAD_probesAge_mergedCommon.csv')
dfKIRC <- fread('KIRC_probesAge_mergedCommon.csv')
dfKIRP <- fread('KIRP_probesAge_mergedCommon.csv')

# other sets to merge as mergedSet2
dfUCEC <- fread('UCEC_probesAge_mergedCommon.csv')
dfPRAD <- fread('PRAD_probesAge_mergedCommon.csv')
dfHNSC <- fread('HNSC_probesAge_mergedCommon.csv')
dfLIHC <- fread('LIHC_probesAge_mergedCommon.csv')
dfCOAD <- fread('COAD_probesAge_mergedCommon.csv')

# dfUCEC[1:5,1:5]
# dfPRAD[1:5,1:5]
# dfHNSC[1:5,1:5]
# dfLIHC[1:5,1:5]
# dfCOAD[1:5,1:5]
# dim(dfCOAD) # 1 256803

nrow(dfTHCA) +# 56
nrow(dfBRCA) +# 96
nrow(dfLUSC) +# 42
nrow(dfLUAD) +# 32
nrow(dfKIRP) +# 45
nrow(dfKIRC) +# 160
nrow(dfUCEC) +# 34
nrow(dfPRAD) +# 50
nrow(dfHNSC) +# 50
nrow(dfLIHC) +# 50
nrow(dfCOAD) # 38
range(dfTHCA$Age) # 15 81
range(dfBRCA$Age) # 28 90
range(dfLUSC$Age) # 40 85
range(dfLUAD$Age) # 42 86
range(dfKIRC$Age) # 38 90
range(dfKIRP$Age) # 31 83

mergedSet <- rbind(dfUCEC,
                    dfPRAD,
                    dfHNSC,
                    dfCOAD,
                    dfLIHC,
                    dfTHCA,
                    dfLUSC,
                    dfKIRC,
                    dfKIRP,
                    dfLUAD,
                    dfBRCA
                    )

write.csv(mergedSet,'TCGA_mergedAll.csv',row.names=F)


mergedKidney <- rbind(dfKIRP,dfKIRC)
write.csv(mergedKidney,'Kidney_probesAge_mergedCommon.csv',row.names=F)
mergedLung <- rbind(dfLUSC,dfLUAD)
write.csv(mergedLung,'Lung_probesAge_mergedCommon.csv',row.names=F)

# merge 'TCGA_BR_KI_LU_TH_probesAge_mergedCommon.csv' and 'TCGA_mergedSet2.csv'
df1 <- fread("TCGA_BR_KI_LU_TH_probesAge_mergedCommon.csv")
df2 <- fread("TCGA_mergedSet2.csv")
mergedAll <- rbind(df1,df2)
write.csv(mergedAll,"TCGA_merged_all.csv",row.names=F)

# projectFiles <- list.files(pattern='_probesAge_mergedCommon.csv') # KIRC, KIRP, LUAD, LUSC, BRCA, THCA
# mergedData <- do.call(rbind,lapply(projectFiles,fread))
mergedData <- rbind(dfKidney,dfLung,dfTHCA,dfBRCA)
write.csv(mergedData,'TCGA_BR_KI_LU_TH_probesAge_mergedCommon.csv',row.names=F)

# df <- fread('Kidney_probesAge_mergedCommon.csv')
# df[1:5,1:5]
# dim(df)
# 
# df1 <- fread('Lung_probesAge_mergedCommon.csv')
# df1[1:5,1:5]
# dim(df1)

df <- fread('TCGA_mergedAll.csv')
dim(df)
df[1:5,1:5]

# sort merged file by SD=======================================================
library(data.table)
library(vroom)
df <- vroom('TCGA_mergedAll.csv')
df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
df1 <- df[,-1]
df1 <- t(df1)
df1[1:5,1:5]
# add row wise SD using apply()
dfSD <- transform(df1,SD=apply(df1,1,sd))
dfSD$SD[1:5]
dfSD_sorted <- dfSD[order(dfSD$SD),] # sort rows by SD (default for decreasing = FALSE)
# dfSD_sorted[1:5,431:432] # from mergedSet2
#           TCGA.E9.A1RI          SD
# cg02582136  0.011503152 0.001993679
# cg23664318  0.009630048 0.002059152
# cg00948102  0.010490091 0.002096094
# cg02755475  0.012293446 0.002125274
# cg21784498  0.015031815 0.002131847
# change dots to hyphens in colnames
names(dfSD_sorted) <- gsub("\\.","-",names(dfSD_sorted))
dfSD_sorted[1:5,1:5]
dfSD_sorted<-subset(dfSD_sorted,select=-(SD))
dfSD_sorted <- t(dfSD_sorted)
dfSD_sorted[1:5,1:5]
#              cg02582136  cg23664318  cg00948102  cg02755475 cg21784498
# TCGA-A4-7288 0.01180093 0.012889454 0.010632637 0.012182993 0.01321553
# TCGA-A4-7585 0.01154645 0.011025859 0.009780075 0.010916993 0.01233136
# TCGA-A4-7732 0.01326496 0.008927545 0.009578577 0.011759337 0.01124974
# TCGA-B1-A47M 0.01480051 0.011747127 0.010595445 0.012200690 0.01290663
# TCGA-BQ-5875 0.01001753 0.008148740 0.009425526 0.009246836 0.01024226

df[1:5,1:5]
dfAge <- df[,1,drop=F] # keep row names using drop=F, subset as matrix instead of vector
dfAge[1:5,1]
dim(dfSD_sorted) # 576 256801
# merge dfAgeID and df_t========================================================
df_merged <- merge(dfAge,dfSD_sorted,by = 0)
df_merged[1:5,1:5]
df_merged <-data.frame(df_merged,row.names=1)
df_merged[1:5,1:5]
write.csv(df_merged,'TCGA_mergedAll_SD.csv')

# check how many probes in dataset before omitting NA and merging with macintyre probes
df_probes <- fread("allNormals_TCGA_LUAD_fall2021.csv")
dim(df_probes)


### venn diagram================================================================
library(data.table)
library(VennDiagram)
# use preprocessed sets, use column 1 as row names and transpose
df1 <- fread("KIRC_probesAge_mergedCommon.csv")
df1[1:5,1:5]
df1 <-data.frame(df1,row.names=1)
df2 <- fread("KIRP_probesAge_mergedCommon.csv")
df3 <- fread("BRCA_probesAge_mergedCommon.csv")
df4 <- fread("LUSC_probesAge_mergedCommon.csv")
df5 <- fread("LUAD_probesAge_mergedCommon.csv")

# size of intersections between sets; use 'Reduce(intersect, list(a,b,c))' for multiple arguments
length(intersect(df1$ProbeID,df2$ProbeID))
dim(df1)[1]

# move to new plotting page
grid.newpage()

# create Venn diagram with five sets
draw.quintuple.venn(area1=dim(df1)[1], area2=dim(df2)[1], area3=dim(df3)[1],
                    area4=dim(df4)[1], area5=dim(df5)[1],
                    n12=length(intersect(df1$ProbeID,df2$ProbeID)),
                    n13=length(intersect(df1$ProbeID,df3$ProbeID)),
                    n14=length(intersect(df1$ProbeID,df4$ProbeID)),
                    n15=length(intersect(df1$ProbeID,df5$ProbeID)),
                    n23=length(intersect(df2$ProbeID,df3$ProbeID)),
                    n24=length(intersect(df2$ProbeID,df4$ProbeID)),
                    n25=length(intersect(df2$ProbeID,df5$ProbeID)),
                    n34=length(intersect(df3$ProbeID,df4$ProbeID)),
                    n35=length(intersect(df3$ProbeID,df5$ProbeID)),
                    n45=length(intersect(df4$ProbeID,df5$ProbeID)),
                    n123=length(Reduce(intersect, list(df1$ProbeID,df2$ProbeID,df3$ProbeID))),
                    n124=length(Reduce(intersect, list(df1$ProbeID,df2$ProbeID,df4$ProbeID))),
                    n125=length(Reduce(intersect, list(df1$ProbeID,df2$ProbeID,df5$ProbeID))),
                    n134=length(Reduce(intersect, list(df1$ProbeID,df3$ProbeID,df4$ProbeID))),
                    n135=length(Reduce(intersect, list(df1$ProbeID,df3$ProbeID,df5$ProbeID))),
                    n145=length(Reduce(intersect, list(df1$ProbeID,df4$ProbeID,df5$ProbeID))),
                    n234=length(Reduce(intersect, list(df2$ProbeID,df3$ProbeID,df4$ProbeID))),
                    n235=length(Reduce(intersect, list(df2$ProbeID,df3$ProbeID,df5$ProbeID))),
                    n245=length(Reduce(intersect, list(df2$ProbeID,df4$ProbeID,df5$ProbeID))),
                    n345=length(Reduce(intersect, list(df3$ProbeID,df4$ProbeID,df5$ProbeID))),
                    n1234=length(Reduce(intersect, list(df1$ProbeID,df2$ProbeID,df3$ProbeID,df4$ProbeID))),
                    n1235=length(Reduce(intersect, list(df1$ProbeID,df2$ProbeID,df3$ProbeID,df5$ProbeID))),
                    n1245=length(Reduce(intersect, list(df1$ProbeID,df2$ProbeID,df4$ProbeID,df5$ProbeID))),
                    n1345=length(Reduce(intersect, list(df1$ProbeID,df3$ProbeID,df4$ProbeID,df5$ProbeID))),
                    n2345=length(Reduce(intersect, list(df2$ProbeID,df3$ProbeID,df4$ProbeID,df5$ProbeID))),
                    n12345=length(Reduce(intersect, list(df1$ProbeID,df2$ProbeID,df3$ProbeID,df4$ProbeID,df5$ProbeID))),
                    category=c("KIRC","KIRP","BRCA","LUSC","LUAD"))
