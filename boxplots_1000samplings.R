# create histogram of age for each data set used (merged, GEO, etc.)===========
library(data.table)
histlist1 <- fread('TCGA_BR_KI_LU_TH_testMSE_lambdaMin.csv')
histlist1 <- as.numeric(histlist1)
histlist2 <- fread('TCGA_BR_KI_LU_TH_minus25KtestMSE_lambdaMin.csv')
histlist3 <- fread('TCGA_BR_KI_LU_TH_minus50KtestMSE_lambdaMin.csv')
histlist4 <- fread('TCGA_BR_KI_LU_TH_minus75KtestMSE_lambdaMin.csv')
histlist4 <- as.numeric(histlist4)
histlist5 <- fread('TCGA_BR_KI_LU_TH_minus100KtestMSE_lambdaMin.csv')
histlist6 <- fread('TCGA_BR_KI_LU_TH_minus125KtestMSE_lambdaMin.csv')
histlist7 <- fread('TCGA_BR_KI_LU_TH_minus150KtestMSE_lambdaMin.csv')
histlist7 <- as.numeric(histlist7)
histlist8 <- fread('TCGA_BR_KI_LU_TH_minus175KtestMSE_lambdaMin.csv')
histlist9 <- fread('TCGA_BR_KI_LU_TH_minus200KtestMSE_lambdaMin.csv')

# QQ plots======================================================================
dev.new(width=15,height=10,noRStudioGD = TRUE)
par(mfrow=c(1,2),pin=c(4,2))
qqnorm(histlist1, 
       ylim=c(3,9),pch=1, frame=FALSE,main="A")#, main="Normal Q-Q Plot all Features 1000 trials 10-fold CV fit test MSE")
qqline(histlist1,col="steelblue",lwd=2)
qqnorm(histlist7, 
       ylim=c(3,9),pch=1, frame=FALSE,main="B",ylab="")#, main="Normal Q-Q Plot -150K lowest variance features 1000 trials 10-fold CV fit test MSE")
qqline(histlist7,col="steelblue",lwd=2)
# 
# line = 6
# cex = 2
# las = 2
# mtext("A", side=2, line=line, cex=cex, las=las)
# # mtext("RMSE", side=1, line=2.2, cex=1, col="black")
# # mtext("B", side=2, line=line, cex=cex, las=las)
# # mtext("C", side=2, line=line, cex=cex, las=las)
# title(ylab=expression(paste("Number of lowest variance probes removed (10"^"3",")")), line=2.4)

# histograms====================================================================
# range(histlist1) # 21.60148 83.16793
hist(histlist1,
     main="Distribution of all Features 1000 trials 10-fold CV fit test MSE",
     breaks=50,
     xlim=c(20,100),
     ylim=c(0,100))
# range(histlist2) # 21.21834 92.16064
hist(histlist2, 
     main="Distribution of -25 highest variance Features 1000 trials 10-fold CV fit test MSE",
     breaks=60,
     xlim=c(20,100),
     ylim=c(0,100))
# compute unpaired two-sample t-test
# convert to RMSE first
histlist1 <- sqrt(histlist1)
histlist7 <- sqrt(histlist7)
ttest <- t.test(histlist1,histlist7,var.equal=F) # x=all features, y=-150K features
ttest

# Welch Two Sample t-test
# 
# data:  histlist1 and histlist7
# t = 4.3285, df = 1982.7, p-value = 1.576e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.5837513 1.5509504
# sample estimates:
#   mean of x mean of y 
# 34.80437  33.73702

# wilcoxon rank sum test
# reject null hypothesis, alternative hypothesis of true location shift not equal to 0
# different medians / distribution
allFeatures_RMSE <- sqrt(histlist1)
minus150K_RMSE <- sqrt(histlist7)
wilcox.test(allFeatures_RMSE,minus150K_RMSE)
# Wilcoxon rank sum test with continuity correction
# 
# data:  histlist1 and histlist7
# W = 551539, p-value = 6.575e-05
# alternative hypothesis: true location shift is not equal to 0

histlist1 <-stack(histlist1)
histlist1 <- histlist1[1]
histlist2 <-stack(histlist2)
histlist2 <- histlist2[1]
histlist3 <-stack(histlist3)
histlist3 <- histlist3[1]
histlist4 <-stack(histlist4)
histlist4 <- histlist4[1]
histlist5 <-stack(histlist5)
histlist5 <- histlist5[1]
histlist6 <-stack(histlist6)
histlist6 <- histlist6[1]
histlist7 <-stack(histlist7)
histlist7 <- histlist7[1]
histlist8 <-stack(histlist8)
histlist8 <- histlist8[1]
histlist9 <-stack(histlist9)
histlist9 <- histlist9[1]
# median(unlist(histlist1)) # 34.84726
# median(unlist(histlist2)) # 34.42873
# median(unlist(histlist3)) # 35.4503
# 
# mean(unlist(histlist1)) # 35.14379
# mean(unlist(histlist2)) # 34.66881
# mean(unlist(histlist3)) # 35.76858
# save all test MSE data into one data frame====================================
listdata <- data.frame(x1=histlist1,
                       x2=histlist2,
                       x3=histlist3,
                       x4=histlist4,
                       x5=histlist5,
                       x6=histlist6,
                       x7=histlist7,
                       x8=histlist8,
                       x9=histlist9)
# convert to RMSE
length(listdata)
listdata[9]
for (i in 1:9){
  listdata[i] = sqrt(listdata[i])
}
dim(listdata)
listdata[1:5,1:5]
median(listdata$values)
range(listdata$values.3) # 19.62670 50.71518 <--- -150K
range(listdata$values.2) # 21.40714 58.34211 <--- -100K
# calculate median values for each column
list_medians <- apply(listdata, 2, median)
list_medians
# values values.1 values.2 values.3 values.4 values.5 values.6 values.7 values.8 
# 5.876341 5.860561 5.870075 5.866863 5.855069 5.795560 5.789078 5.868769 6.047802 
# listdata$median[1:6]
# listdata[1:10,1:6]

# list means
list_means <- apply(listdata, 2, mean)
list_means
# values values.1 values.2 values.3 values.4 values.5 values.6 values.7 values.8 
# 5.880109 5.884252 5.884519 5.883096 5.869347 5.810648 5.790702 5.867459 6.043872 

# write.csv(listdata,"TCGA_BR_KI_LU_TH_1000samplingsTestMSE_merged2.csv",row.names=F)

# test RMSE boxplots=============================================================
# pdf("1000samplings_testMSE_boxplots.pdf")
# x11(width = 20, height = 20)
dev.new(width=20,height=20,noRStudioGD = TRUE)
par(mfrow=c(3,1),pin=c(5.5,2))
boxplot(listdata,
        ylim=c(4,12),
        #ylab="Number of lowest variance probes removed from subset (thousands)",
        notch = T,
        horizontal=T,
        las=1, # style of axis labels set to horizontal
        names=c("0",
                "25",
                "50",
                "75",
                "100",
                "125",
                "150",
                "175",
                "200")#,
        #xlab="RMSE"
        #main="1000 samplings 75/25 train/test MSE for varying subsets of features"
        )
# return statistics used to create boxplot for plot labels
# fiddle with y value until you have what you want
# text(x=fivenum(listdata), labels=fivenum(listdata),y=1.25) 
line = 6
cex = 2
las = 2
mtext("A", side=2, line=line, cex=cex, las=las)
# mtext("RMSE", side=1, line=2.2, cex=1, col="black")
# mtext("B", side=2, line=line, cex=cex, las=las)
# mtext("C", side=2, line=line, cex=cex, las=las)
title(ylab=expression(paste("Number of lowest variance probes removed (10"^"3",")")), line=2.4)
# mtext("Margins", side=1, line=2.5, cex=1)
# dev.off()
# test MSE histograms===========================================================

# find lower and upper range of actual values
range(histlist1)# 22.52487 64.59388
range(histlist2)# 22.20722 52.10697
range(histlist3)# 22.47476 54.59826
par(mfrow=c(3,1)) # split display screen into separate panels
hist(listdata[[1]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE all features")

hist(listdata[[2]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE minus 50K features")

hist(listdata[[3]],xlim=c(20,100),breaks=20,
     xlab="MSE",
     main="1000 samplings test MSE minus 100K features")

# df <- fread('TCGA_BR_KI_LU_TH_probesAge_mergedCommon.csv')
# df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
# library(caret)
# numFeatures_nearZeroVar <- length(nearZeroVar(df)) #
# # remove nearZeroVar columns
# dfnearZeroVar <- df[,-nearZeroVar(df)] # integer(0)

par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(20,100),
        ylab="MSE",
        notch = T,
        names=c("A",
                "B",
                "C",
                "D",
                "E"))
#main="1000 samplings 75/25 train/test MSE for varying subsets of features")

# boxplots 1000samplings reversed===============================================
library(data.table)
histlist1 <- fread('TCGA_BR_KI_LU_TH_testMSE_lambdaMin.csv')
# histlist1 <- as.numeric(histlist1)
histlist2 <- fread('TCGA_mergedSet1_minus25KtestMSE_reversed.csv')
# histlist2 <- as.numeric(histlist2)
histlist3 <- fread('TCGA_mergedSet1_minus50KtestMSE_reversed.csv')
histlist4 <- fread('TCGA_mergedSet1_minus75KtestMSE_reversed.csv')
histlist5 <- fread('TCGA_mergedSet1_minus100KtestMSE_reversed.csv')
histlist6 <- fread('TCGA_mergedSet1_minus125KtestMSE_reversed.csv')
histlist7 <- fread('TCGA_mergedSet1_minus150KtestMSE_reversed.csv')
histlist8 <- fread('TCGA_mergedSet1_minus175KtestMSE_reversed.csv')
histlist9 <- fread('TCGA_mergedSet1_minus200KtestMSE_reversed.csv')


# compute unpaired two-sample t-test
ttest <- t.test(histlist1,histlist2,var.equal=F) # x=all features, y=-150K features
ttest
# Welch Two Sample t-test
# 
# data:  histlist1 and histlist2
# t = 0.14598, df = 1997.4, p-value = 0.884
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4710187  0.5467769
# sample estimates:
#   mean of x mean of y 
# 34.80437  34.76649 

# wilcoxon rank sum test
# reject null hypothesis, alternative hypothesis of true location shift not equal to 0
# different medians / distribution
allFeatures_RMSE <- sqrt(histlist1)
minus25K_highVarProbes_RMSE <- sqrt(histlist2)
wilcox.test(allFeatures_RMSE,minus25K_highVarProbes_RMSE)
# Wilcoxon rank sum test with continuity correction
# 
# data:  allFeatures_RMSE and minus25K_highVarProbes_RMSE
# W = 503526, p-value = 0.7848
# alternative hypothesis: true location shift is not equal to 0

histlist1 <-stack(histlist1)
histlist1 <- histlist1[1]
histlist2 <-stack(histlist2)
histlist2 <- histlist2[1]
histlist3 <-stack(histlist3)
histlist3 <- histlist3[1]
histlist4 <-stack(histlist4)
histlist4 <- histlist4[1]
histlist5 <-stack(histlist5)
histlist5 <- histlist5[1]
histlist6 <-stack(histlist6)
histlist6 <- histlist6[1]
histlist7 <-stack(histlist7)
histlist7 <- histlist7[1]
histlist8 <-stack(histlist8)
histlist8 <- histlist8[1]
histlist9 <-stack(histlist9)
histlist9 <- histlist9[1]
# median(unlist(histlist1)) # 34.84726
# median(unlist(histlist2)) # 34.42873
# median(unlist(histlist3)) # 35.4503
# 
# mean(unlist(histlist1)) # 35.14379
# mean(unlist(histlist2)) # 34.66881
# mean(unlist(histlist3)) # 35.76858
# save all test MSE data into one data frame====================================
listdata <- data.frame(x1=histlist1,
                       x2=histlist2,
                       x3=histlist3,
                       x4=histlist4,
                       x5=histlist5,
                       x6=histlist6,
                       x7=histlist7,
                       x8=histlist8,
                       x9=histlist9)
# convert to RMSE
length(listdata)
listdata[9]
for (i in 1:9){
  listdata[i] = sqrt(listdata[i])
}
dim(listdata)
listdata[1:5,1:5]
median(listdata$values)
range(listdata$values.3) # 19.62670 50.71518 <--- -150K
range(listdata$values.2) # 21.40714 58.34211 <--- -100K
# calculate median values for each column
list_medians <- apply(listdata, 2, median)
list_medians
# values values.1 values.2 values.3 values.4 
# 34.53138 34.45778 34.28183 33.51342 36.57591
listdata$median[1:6]
listdata[1:10,1:6]
# list means
list_means <- apply(listdata, 2, mean)
list_means
# write.csv(listdata,"TCGA_BR_KI_LU_TH_1000samplingsTestMSE_merged2.csv",row.names=F)

# test MSE boxplots=============================================================
# pdf("1000samplings_reversed_testMSE_boxplots.pdf")
# par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(4,12),
        #ylab=expression(paste("Number of highest variance probes removed from subset (10"^"3*",")")),
        notch = T,
        horizontal=T,
        las=1, # style of axis labels set to horizontal
        names=c("0",
                "25",
                "50",
                "75",
                "100",
                "125",
                "150",
                "175",
                "200")
        # xlab="RMSE"
        #main="1000 samplings 75/25 train/test MSE for varying subsets of features"
        )
mtext("B", side=2, line=line, cex=cex, las=las)
# mtext("RMSE", side=1, line=2.2, cex=1, col="black")
title(ylab=expression(paste("Number of highest variance probes removed (10"^"3",")")), line=2.4)

# dev.off()
# test MSE histograms===========================================================

# find lower and upper range of actual values
range(histlist1)# 22.52487 64.59388
range(histlist2)# 22.20722 52.10697
range(histlist3)# 22.47476 54.59826
par(mfrow=c(2,2)) # split display screen into separate panels
hist(listdata[[1]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE all features")

hist(listdata[[2]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE minus 50K features")

hist(listdata[[3]],xlim=c(20,100),breaks=20,
     xlab="MSE",
     main="1000 samplings test MSE minus 100K features")

# df <- fread('TCGA_BR_KI_LU_TH_probesAge_mergedCommon.csv')
# df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
# library(caret)
# numFeatures_nearZeroVar <- length(nearZeroVar(df)) #
# # remove nearZeroVar columns
# dfnearZeroVar <- df[,-nearZeroVar(df)] # integer(0)

par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(20,100),
        ylab="MSE",
        notch = T,
        names=c("A",
                "B",
                "C",
                "D",
                "E"))
#main="1000 samplings 75/25 train/test MSE for varying subsets of features")

# boxplots 1000samplings 100K subsets===============================================
library(data.table)
histlist1 <- fread('TCGA_mergedSet1_0to100KtestMSE.csv') # should be from 1 to 100K***

histlist2 <- fread('TCGA_mergedSet1_25Kto125KtestMSE.csv')
# histlist2 <- as.numeric(histlist2)
histlist3 <- fread('TCGA_mergedSet1_50Kto150KtestMSE.csv')
histlist4 <- fread('TCGA_mergedSet1_75Kto175KtestMSE.csv')
histlist5 <- fread('TCGA_mergedSet1_100Kto200KtestMSE.csv')
histlist6 <- fread('TCGA_mergedSet1_125Kto225KtestMSE.csv')
# histlist6 <- as.numeric(histlist6)
histlist7 <- fread('TCGA_mergedSet1_150Kto250KtestMSE.csv')
# histlist7 <- as.numeric(histlist7)

# compute unpaired two-sample t-test
ttest <- t.test(histlist6,histlist7,var.equal=F) # x=all features, y=-150K features
ttest
# Welch Two Sample t-test
# 
# data:  histlist2 and histlist7
# t = 1.0079, df = 1979.2, p-value = 0.3136
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2344556  0.7302732
# sample estimates:
#   mean of x mean of y 
# 34.08194  33.83403

# wilcoxon rank sum test
# reject null hypothesis, alternative hypothesis of true location shift not equal to 0
# different medians / distribution
RMSE_25Kto125K <- sqrt(histlist2)
RMSE_150Kto250K <- sqrt(histlist7)
wilcox.test(RMSE_25Kto125K,RMSE_150Kto250K)
# Wilcoxon rank sum test with continuity correction
# 
# data:  RMSE_25Kto125K and RMSE_150Kto250K
# W = 508870, p-value = 0.4922
# alternative hypothesis: true location shift is not equal to 0

histlist1 <-stack(histlist1)
histlist1 <- histlist1[1]
histlist2 <-stack(histlist2)
histlist2 <- histlist2[1]
histlist3 <-stack(histlist3)
histlist3 <- histlist3[1]
histlist4 <-stack(histlist4)
histlist4 <- histlist4[1]
histlist5 <-stack(histlist5)
histlist5 <- histlist5[1]
histlist6 <-stack(histlist6)
histlist6 <- histlist6[1]
histlist7 <-stack(histlist7)
histlist7 <- histlist7[1]


# median(unlist(histlist1)) # 34.84726
# median(unlist(histlist2)) # 34.42873
# median(unlist(histlist3)) # 35.4503
# 
# mean(unlist(histlist1)) # 35.14379
# mean(unlist(histlist2)) # 34.66881
# mean(unlist(histlist3)) # 35.76858
# save all test MSE data into one data frame====================================
listdata <- data.frame(x1=histlist1,
                       x2=histlist2,
                       x3=histlist3,
                       x4=histlist4,
                       x5=histlist5,
                       x6=histlist6,
                       x7=histlist7
                       )
# convert to RMSE
length(listdata)
# listdata[8]
for (i in 1:7){
  listdata[i] = sqrt(listdata[i])
}
dim(listdata)
listdata[1:5,1:5]
median(listdata$values)
range(listdata$values.3) # 19.62670 50.71518 <--- -150K
range(listdata$values.2) # 21.40714 58.34211 <--- -100K
# calculate median values for each column
list_medians <- apply(listdata, 2, median)
list_medians
# listdata$median[1:6]
# listdata[1:10,1:6]
# 0-100K    25K-125K 50-150K  75K-175K 100K-200K 125K-225K 150K-250K 
# 7.528993 5.809412 6.421636 6.219771  6.160157  5.889370  5.796417 
# calculate average values for each column
list_means <- apply(listdata, 2, mean)
list_means
# 0-100K    25K-125K 50-150K  75K-175K 100K-200K 125K-225K 150K-250K 
# 7.536545 6.763968 6.436773 6.235975 6.163197 5.912452 5.799420 

# values from -150K lowest variance probes (median/mean)
# 5.789078 / 5.790702

# write.csv(listdata,"TCGA_BR_KI_LU_TH_1000samplingsTestMSE_merged2.csv",row.names=F)

# test MSE boxplots=============================================================
# pdf("1000samplings_100Ksubsets_testRMSE_boxplots.pdf")
# par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(4,12),
        #ylab=("Upper index position of 100,000 features in subset"),
        notch = T,
        horizontal=T,
        las=1, # style of axis labels set to horizontal
        names=c("1",
                "25",
                "50",
                "75",
                "100",
                "125",
                "150"
                )
        #xlab="RMSE"
        #main="1000 samplings 75/25 train/test RMSE for varying subsets of features"
        )
mtext("C", side=2, line=line, cex=cex, las=las)
mtext("RMSE", side=1, line=2.2, cex=1, col="black")
title(ylab=expression(paste("Starting index position of 100,000 features in subset (10"^"3",")")), line=2.4)
# dev.off()
# test MSE histograms===========================================================

# find lower and upper range of actual values
range(histlist1)# 22.52487 64.59388
range(histlist2)# 22.20722 52.10697
range(histlist3)# 22.47476 54.59826
par(mfrow=c(2,2)) # split display screen into separate panels
hist(listdata[[1]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE all features")

hist(listdata[[2]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE minus 50K features")

hist(listdata[[3]],xlim=c(20,100),breaks=20,
     xlab="MSE",
     main="1000 samplings test MSE minus 100K features")

# df <- fread('TCGA_BR_KI_LU_TH_probesAge_mergedCommon.csv')
# df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
# library(caret)
# numFeatures_nearZeroVar <- length(nearZeroVar(df)) #
# # remove nearZeroVar columns
# dfnearZeroVar <- df[,-nearZeroVar(df)] # integer(0)

par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(20,100),
        ylab="MSE",
        notch = T,
        names=c("A",
                "B",
                "C",
                "D",
                "E"))

# mergedSet3 minus low variance probes==========================================
library(data.table)
histlist1 <- fread('mergedSet3_allFeatures_testRMSE.csv')
# histlist1 <- as.numeric(histlist1)
histlist2 <- fread('mergedSet3_minus25K_testRMSE.csv')
histlist3 <- fread('mergedSet3_minus50K_testRMSE.csv')
histlist4 <- fread('mergedSet3_minus75K_testRMSE.csv')
histlist5 <- fread('mergedSet3_minus100K_testRMSE.csv')
histlist6 <- fread('mergedSet3_minus125K_testRMSE.csv')
histlist7 <- fread('mergedSet3_minus150K_testRMSE.csv')
# histlist7 <- as.numeric(histlist7)
histlist8 <- fread('mergedSet3_minus175K_testRMSE.csv')
histlist9 <- fread('mergedSet3_minus200K_testRMSE.csv')

# QQ plots======================================================================
# qqnorm(histlist1, pch=1, frame=FALSE, main="Normal Q-Q Plot all Features 1000 trials 10-fold CV fit test MSE")
# qqline(histlist1,col="steelblue",lwd=2)
# qqnorm(histlist2, pch=1, frame=FALSE, main="Normal Q-Q Plot -25K highest variance features 1000 trials 10-fold CV fit test MSE")
# qqline(histlist2,col="steelblue",lwd=2)
# histograms====================================================================
# range(histlist1) # 21.60148 83.16793
hist(histlist1,
     main="Distribution of all Features 1000 trials 10-fold CV fit test MSE",
     breaks=50,
     xlim=c(20,100),
     ylim=c(0,100))
# range(histlist2) # 21.21834 92.16064
hist(histlist2, 
     main="Distribution of -25 highest variance Features 1000 trials 10-fold CV fit test MSE",
     breaks=60,
     xlim=c(20,100),
     ylim=c(0,100))
# compute unpaired two-sample t-test
ttest <- t.test(histlist1,histlist7,var.equal=F) # x=all features, y=-150K features
ttest
# Welch Two Sample t-test
# 
# data:  histlist1 and histlist7
# t = 4.3285, df = 1982.7, p-value = 1.576e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.5837513 1.5509504
# sample estimates:
#   mean of x mean of y 
# 34.80437  33.73702

# wilcoxon rank sum test
# reject null hypothesis, alternative hypothesis of true location shift not equal to 0
# different medians / distribution
allFeatures_RMSE <- sqrt(histlist1)
minus150K_RMSE <- sqrt(histlist7)
wilcox.test(allFeatures_RMSE,minus150K_RMSE)
# Wilcoxon rank sum test with continuity correction
# 
# data:  histlist1 and histlist7
# W = 551539, p-value = 6.575e-05
# alternative hypothesis: true location shift is not equal to 0

histlist1 <-stack(histlist1)
histlist1 <- histlist1[1]
histlist2 <-stack(histlist2)
histlist2 <- histlist2[1]
histlist3 <-stack(histlist3)
histlist3 <- histlist3[1]
histlist4 <-stack(histlist4)
histlist4 <- histlist4[1]
histlist5 <-stack(histlist5)
histlist5 <- histlist5[1]
histlist6 <-stack(histlist6)
histlist6 <- histlist6[1]
histlist7 <-stack(histlist7)
histlist7 <- histlist7[1]
histlist8 <-stack(histlist8)
histlist8 <- histlist8[1]
histlist9 <-stack(histlist9)
histlist9 <- histlist9[1]
# median(unlist(histlist1)) # 34.84726
# median(unlist(histlist2)) # 34.42873
# median(unlist(histlist3)) # 35.4503
# 
# mean(unlist(histlist1)) # 35.14379
# mean(unlist(histlist2)) # 34.66881
# mean(unlist(histlist3)) # 35.76858
# save all test MSE data into one data frame====================================
listdata <- data.frame(x1=histlist1,
                       x2=histlist2,
                       x3=histlist3,
                       x4=histlist4,
                       x5=histlist5,
                       x6=histlist6,
                       x7=histlist7,
                       x8=histlist8,
                       x9=histlist9)
# convert to RMSE
# length(listdata)
# listdata[9]
for (i in 1:9){
  listdata[i] = sqrt(listdata[i])
}
dim(listdata)
listdata[1:5,1:5]
median(listdata$values)
range(listdata$values.3) # 19.62670 50.71518 <--- -150K
range(listdata$values.2) # 21.40714 58.34211 <--- -100K
range(listdata)
# calculate median values for each column
list_medians <- apply(listdata, 2, median)
list_medians
# values values.1 values.2 values.3 values.4 values.5 values.6 values.7 values.8 
# 5.876341 5.860561 5.870075 5.866863 5.855069 5.795560 5.789078 5.868769 6.047802 
# listdata$median[1:6]
# listdata[1:10,1:6]

# list means
list_means <- apply(listdata, 2, mean)
list_means
# values values.1 values.2 values.3 values.4 values.5 values.6 values.7 values.8 
# 5.880109 5.884252 5.884519 5.883096 5.869347 5.810648 5.790702 5.867459 6.043872 

# write.csv(listdata,"TCGA_BR_KI_LU_TH_1000samplingsTestMSE_merged2.csv",row.names=F)

# test RMSE boxplots=============================================================
pdf("1000samplings_mergedSet3_testRMSE_boxplots.pdf")
# par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(4,10),
        #ylab="Number of lowest variance probes removed from subset (thousands)",
        notch = T,
        horizontal=T,
        las=1, # style of axis labels set to horizontal
        names=c("0",
                "25",
                "50",
                "75",
                "100",
                "125",
                "150",
                "175",
                "200"),
        xlab="RMSE"
        #main="1000 samplings 75/25 train/test MSE for varying subsets of features"
)
title(ylab=expression(paste("Number of lowest variance probes removed (10"^"3",")")), line=2.4)
dev.off()
# test MSE histograms===========================================================

# find lower and upper range of actual values
range(histlist1)# 22.52487 64.59388
range(histlist2)# 22.20722 52.10697
range(histlist3)# 22.47476 54.59826
par(mfrow=c(2,2)) # split display screen into separate panels
hist(listdata[[1]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE all features")

hist(listdata[[2]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE minus 50K features")

hist(listdata[[3]],xlim=c(20,100),breaks=20,
     xlab="MSE",
     main="1000 samplings test MSE minus 100K features")

# df <- fread('TCGA_BR_KI_LU_TH_probesAge_mergedCommon.csv')
# df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
# library(caret)
# numFeatures_nearZeroVar <- length(nearZeroVar(df)) #
# # remove nearZeroVar columns
# dfnearZeroVar <- df[,-nearZeroVar(df)] # integer(0)

par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(20,100),
        ylab="MSE",
        notch = T,
        names=c("A",
                "B",
                "C",
                "D",
                "E"))
#main="1000 samplings 75/25 train/test MSE for varying subsets of features")

# mergedAll minus low variance probes==========================================
library(data.table)
histlist1 <- fread('mergedAll_allFeatures_testRMSE.csv')
# histlist1 <- as.numeric(histlist1)
histlist2 <- fread('mergedAll_minus25K_testRMSE.csv')
histlist3 <- fread('mergedAll_minus50K_testRMSE.csv')
histlist4 <- fread('mergedAll_minus75K_testRMSE.csv')
histlist5 <- fread('mergedAll_minus100K_testRMSE.csv')
histlist6 <- fread('mergedAll_minus125K_testRMSE.csv')
histlist7 <- fread('mergedAll_minus150K_testRMSE.csv')
# histlist7 <- as.numeric(histlist7)
histlist8 <- fread('mergedAll_minus175K_testRMSE.csv')
histlist9 <- fread('mergedAll_minus200K_testRMSE.csv')

# QQ plots======================================================================
# qqnorm(histlist1, pch=1, frame=FALSE, main="Normal Q-Q Plot all Features 1000 trials 10-fold CV fit test MSE")
# qqline(histlist1,col="steelblue",lwd=2)
# qqnorm(histlist2, pch=1, frame=FALSE, main="Normal Q-Q Plot -25K highest variance features 1000 trials 10-fold CV fit test MSE")
# qqline(histlist2,col="steelblue",lwd=2)
# histograms====================================================================
# range(histlist1) # 21.60148 83.16793
hist(histlist1,
     main="Distribution of all Features 1000 trials 10-fold CV fit test MSE",
     breaks=50,
     xlim=c(20,100),
     ylim=c(0,100))
# range(histlist2) # 21.21834 92.16064
hist(histlist2, 
     main="Distribution of -25 highest variance Features 1000 trials 10-fold CV fit test MSE",
     breaks=60,
     xlim=c(20,100),
     ylim=c(0,100))
# compute unpaired two-sample t-test
ttest <- t.test(histlist1,histlist7,var.equal=F) # x=all features, y=-150K features
ttest
# Welch Two Sample t-test
# 
# data:  histlist1 and histlist7
# t = 4.3285, df = 1982.7, p-value = 1.576e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.5837513 1.5509504
# sample estimates:
#   mean of x mean of y 
# 34.80437  33.73702

# wilcoxon rank sum test
# reject null hypothesis, alternative hypothesis of true location shift not equal to 0
# different medians / distribution
allFeatures_RMSE <- sqrt(histlist1)
minus150K_RMSE <- sqrt(histlist7)
wilcox.test(allFeatures_RMSE,minus150K_RMSE)
# Wilcoxon rank sum test with continuity correction
# 
# data:  histlist1 and histlist7
# W = 551539, p-value = 6.575e-05
# alternative hypothesis: true location shift is not equal to 0

histlist1 <-stack(histlist1)
histlist1 <- histlist1[1]
histlist2 <-stack(histlist2)
histlist2 <- histlist2[1]
histlist3 <-stack(histlist3)
histlist3 <- histlist3[1]
histlist4 <-stack(histlist4)
histlist4 <- histlist4[1]
histlist5 <-stack(histlist5)
histlist5 <- histlist5[1]
histlist6 <-stack(histlist6)
histlist6 <- histlist6[1]
histlist7 <-stack(histlist7)
histlist7 <- histlist7[1]
histlist8 <-stack(histlist8)
histlist8 <- histlist8[1]
histlist9 <-stack(histlist9)
histlist9 <- histlist9[1]
# median(unlist(histlist1)) # 34.84726
# median(unlist(histlist2)) # 34.42873
# median(unlist(histlist3)) # 35.4503
# 
# mean(unlist(histlist1)) # 35.14379
# mean(unlist(histlist2)) # 34.66881
# mean(unlist(histlist3)) # 35.76858
# save all test MSE data into one data frame====================================
listdata <- data.frame(x1=histlist1,
                       x2=histlist2,
                       x3=histlist3,
                       x4=histlist4,
                       x5=histlist5,
                       x6=histlist6,
                       x7=histlist7,
                       x8=histlist8,
                       x9=histlist9)
# convert to RMSE
# length(listdata)
# listdata[9]
for (i in 1:9){
  listdata[i] = sqrt(listdata[i])
}
dim(listdata)
listdata[1:5,1:5]
median(listdata$values)
range(listdata$values.3) # 19.62670 50.71518 <--- -150K
range(listdata$values.2) # 21.40714 58.34211 <--- -100K
range(listdata)
# calculate median values for each column
list_medians <- apply(listdata, 2, median)
list_medians
# values values.1 values.2 values.3 values.4 values.5 values.6 values.7 values.8 
# 5.876341 5.860561 5.870075 5.866863 5.855069 5.795560 5.789078 5.868769 6.047802 
# listdata$median[1:6]
# listdata[1:10,1:6]

# list means
list_means <- apply(listdata, 2, mean)
list_means
# values values.1 values.2 values.3 values.4 values.5 values.6 values.7 values.8 
# 5.880109 5.884252 5.884519 5.883096 5.869347 5.810648 5.790702 5.867459 6.043872 

# write.csv(listdata,"TCGA_BR_KI_LU_TH_1000samplingsTestMSE_merged2.csv",row.names=F)

# test RMSE boxplots=============================================================
pdf("1000samplings_mergedAll_testRMSE_boxplots.pdf")
# par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(4,10),
        #ylab="Number of lowest variance probes removed from subset (thousands)",
        notch = T,
        horizontal=T,
        las=1, # style of axis labels set to horizontal
        names=c("0",
                "25",
                "50",
                "75",
                "100",
                "125",
                "150",
                "175",
                "200"),
        xlab="RMSE"
        main="1000 trials of 10-fold CV 75/25 train/test RMSE (n=651)"
)
title(ylab=expression(paste("Number of lowest variance probes removed (10"^"3",")")), line=2.4)
dev.off()
# test MSE histograms===========================================================

# find lower and upper range of actual values
range(histlist1)# 22.52487 64.59388
range(histlist2)# 22.20722 52.10697
range(histlist3)# 22.47476 54.59826
par(mfrow=c(2,2)) # split display screen into separate panels
hist(listdata[[1]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE all features")

hist(listdata[[2]],xlim=c(20,100),breaks=30,
     xlab="MSE",
     main="1000 samplings test MSE minus 50K features")

hist(listdata[[3]],xlim=c(20,100),breaks=20,
     xlab="MSE",
     main="1000 samplings test MSE minus 100K features")

# df <- fread('TCGA_BR_KI_LU_TH_probesAge_mergedCommon.csv')
# df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
# library(caret)
# numFeatures_nearZeroVar <- length(nearZeroVar(df)) #
# # remove nearZeroVar columns
# dfnearZeroVar <- df[,-nearZeroVar(df)] # integer(0)

par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(20,100),
        ylab="MSE",
        notch = T,
        names=c("A",
                "B",
                "C",
                "D",
                "E"))

### hannum boxplots ============================================================
# create histogram of age for each data set used (merged, GEO, etc.)===========
library(data.table)
histlist1 <- fread('hannum_allFeatures_testRMSE.csv')
# histlist1 <- as.numeric(histlist1)
histlist2 <- fread('hannum_minus25K_testRMSE.csv')
histlist3 <- fread('hannum_minus50K_testRMSE.csv')
histlist4 <- fread('hannum_minus75K_testRMSE.csv')
histlist5 <- fread('hannum_minus100K_testRMSE.csv')
histlist6 <- fread('hannum_minus125K_testRMSE.csv')
histlist7 <- fread('hannum_minus150K_testRMSE.csv')
# histlist7 <- as.numeric(histlist7)
histlist8 <- fread('hannum_minus175K_testRMSE.csv')
histlist9 <- fread('hannum_minus200K_testRMSE.csv')

# QQ plots======================================================================
# qqnorm(histlist1, pch=1, frame=FALSE, main="Normal Q-Q Plot all Features 1000 trials 10-fold CV fit test MSE")
# qqline(histlist1,col="steelblue",lwd=2)
# qqnorm(histlist2, pch=1, frame=FALSE, main="Normal Q-Q Plot -25K highest variance features 1000 trials 10-fold CV fit test MSE")
# qqline(histlist2,col="steelblue",lwd=2)
# histograms====================================================================
# range(histlist1) # 21.60148 83.16793
hist(histlist1,
     main="Distribution of all Features 1000 trials 10-fold CV fit test MSE",
     breaks=50,
     xlim=c(20,100),
     ylim=c(0,100))
# range(histlist2) # 21.21834 92.16064
hist(histlist2, 
     main="Distribution of -25 highest variance Features 1000 trials 10-fold CV fit test MSE",
     breaks=60,
     xlim=c(20,100),
     ylim=c(0,100))
# compute unpaired two-sample t-test
ttest <- t.test(histlist1,histlist7,var.equal=F) # x=all features, y=-150K features
ttest
# Welch Two Sample t-test
# 
# data:  histlist1 and histlist7
# t = 4.3285, df = 1982.7, p-value = 1.576e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.5837513 1.5509504
# sample estimates:
#   mean of x mean of y 
# 34.80437  33.73702

# wilcoxon rank sum test
# reject null hypothesis, alternative hypothesis of true location shift not equal to 0
# different medians / distribution
allFeatures_RMSE <- sqrt(histlist1)
minus150K_RMSE <- sqrt(histlist7)
wilcox.test(allFeatures_RMSE,minus150K_RMSE)
# Wilcoxon rank sum test with continuity correction
# 
# data:  histlist1 and histlist7
# W = 551539, p-value = 6.575e-05
# alternative hypothesis: true location shift is not equal to 0

histlist1 <-stack(histlist1)
histlist1 <- histlist1[1]
histlist2 <-stack(histlist2)
histlist2 <- histlist2[1]
histlist3 <-stack(histlist3)
histlist3 <- histlist3[1]
histlist4 <-stack(histlist4)
histlist4 <- histlist4[1]
histlist5 <-stack(histlist5)
histlist5 <- histlist5[1]
histlist6 <-stack(histlist6)
histlist6 <- histlist6[1]
histlist7 <-stack(histlist7)
histlist7 <- histlist7[1]
histlist8 <-stack(histlist8)
histlist8 <- histlist8[1]
histlist9 <-stack(histlist9)
histlist9 <- histlist9[1]
# median(unlist(histlist1)) # 34.84726
# median(unlist(histlist2)) # 34.42873
# median(unlist(histlist3)) # 35.4503
# 
# mean(unlist(histlist1)) # 35.14379
# mean(unlist(histlist2)) # 34.66881
# mean(unlist(histlist3)) # 35.76858
# save all test MSE data into one data frame====================================
listdata <- data.frame(x1=histlist1,
                       x2=histlist2,
                       x3=histlist3,
                       x4=histlist4,
                       x5=histlist5,
                       x6=histlist6,
                       x7=histlist7,
                       x8=histlist8,
                       x9=histlist9)
# convert to RMSE
length(listdata)
listdata[9]
for (i in 1:9){
  listdata[i] = sqrt(listdata[i])
}
dim(listdata)
listdata[1:5,1:5]
median(listdata$values)
range(listdata$values.3) # 19.62670 50.71518 <--- -150K
range(listdata$values.2) # 21.40714 58.34211 <--- -100K
# calculate median values for each column
list_medians <- apply(listdata, 2, median)
list_medians
# values values.1 values.2 values.3 values.4 values.5 values.6 values.7 values.8 
# 5.876341 5.860561 5.870075 5.866863 5.855069 5.795560 5.789078 5.868769 6.047802 
# listdata$median[1:6]
# listdata[1:10,1:6]

# list means
list_means <- apply(listdata, 2, mean)
list_means
# values values.1 values.2 values.3 values.4 values.5 values.6 values.7 values.8 
# 5.880109 5.884252 5.884519 5.883096 5.869347 5.810648 5.790702 5.867459 6.043872 

# write.csv(listdata,"TCGA_BR_KI_LU_TH_1000samplingsTestMSE_merged2.csv",row.names=F)

# test RMSE boxplots=============================================================
# pdf("hannum_1000samplings_testRMSE_boxplots.pdf")
# par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(3,7),
        #ylab="Number of lowest variance probes removed from subset (thousands)",
        notch = T,
        horizontal=T,
        las=1, # style of axis labels set to horizontal
        names=c("0",
                "25",
                "50",
                "75",
                "100",
                "125",
                "150",
                "175",
                "200"),
        xlab="RMSE"
        #main="1000 samplings 75/25 train/test MSE for varying subsets of features"
)
# return statistics used to create boxplot for plot labels
# fiddle with y value until you have what you want
text(x=fivenum(listdata), labels=fivenum(listdata),y=1.25) 
title(ylab=expression(paste("Number of lowest variance probes removed (10 "^"  3",")")), line=2.4)
# dev.off()

# hannum boxplots 1000samplings 100K subsets===============================================
library(data.table)
histlist1 <- fread('hannum_0to100KtestMSE.csv')

histlist2 <- fread('hannum_25Kto125KtestMSE.csv')
# histlist2 <- as.numeric(histlist2)
histlist3 <- fread('hannum_50Kto150KtestMSE.csv')
histlist4 <- fread('hannum_75Kto175KtestMSE.csv')
histlist5 <- fread('hannum_100Kto200KtestMSE.csv')
histlist6 <- fread('hannum_125Kto225KtestMSE.csv')
# histlist6 <- as.numeric(histlist6)
histlist7 <- fread('hannum_150Kto250KtestMSE.csv')
# histlist7 <- as.numeric(histlist7)

# compute unpaired two-sample t-test
ttest <- t.test(histlist6,histlist7,var.equal=F) # x=all features, y=-150K features
ttest
# Welch Two Sample t-test
# 
# data:  histlist2 and histlist7
# t = 1.0079, df = 1979.2, p-value = 0.3136
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2344556  0.7302732
# sample estimates:
#   mean of x mean of y 
# 34.08194  33.83403

# wilcoxon rank sum test
# reject null hypothesis, alternative hypothesis of true location shift not equal to 0
# different medians / distribution
RMSE_25Kto125K <- sqrt(histlist2)
RMSE_150Kto250K <- sqrt(histlist7)
wilcox.test(RMSE_25Kto125K,RMSE_150Kto250K)
# Wilcoxon rank sum test with continuity correction
# 
# data:  RMSE_25Kto125K and RMSE_150Kto250K
# W = 508870, p-value = 0.4922
# alternative hypothesis: true location shift is not equal to 0

histlist1 <-stack(histlist1)
histlist1 <- histlist1[1]
histlist2 <-stack(histlist2)
histlist2 <- histlist2[1]
histlist3 <-stack(histlist3)
histlist3 <- histlist3[1]
histlist4 <-stack(histlist4)
histlist4 <- histlist4[1]
histlist5 <-stack(histlist5)
histlist5 <- histlist5[1]
histlist6 <-stack(histlist6)
histlist6 <- histlist6[1]
histlist7 <-stack(histlist7)
histlist7 <- histlist7[1]


# median(unlist(histlist1)) # 34.84726
# median(unlist(histlist2)) # 34.42873
# median(unlist(histlist3)) # 35.4503
# 
# mean(unlist(histlist1)) # 35.14379
# mean(unlist(histlist2)) # 34.66881
# mean(unlist(histlist3)) # 35.76858
# save all test MSE data into one data frame====================================
listdata <- data.frame(x1=histlist1,
                       x2=histlist2,
                       x3=histlist3,
                       x4=histlist4,
                       x5=histlist5,
                       x6=histlist6,
                       x7=histlist7
)
# convert to RMSE
length(listdata)
# listdata[8]
for (i in 1:7){
  listdata[i] = sqrt(listdata[i])
}
dim(listdata)
listdata[1:5,1:5]
median(listdata$values)
range(listdata$values.3) # 19.62670 50.71518 <--- -150K
range(listdata$values.2) # 21.40714 58.34211 <--- -100K
# calculate median values for each column
list_medians <- apply(listdata, 2, median)
list_medians
# listdata$median[1:6]
# listdata[1:10,1:6]
# 0-100K    25K-125K 50-150K  75K-175K 100K-200K 125K-225K 150K-250K 
# 7.528993 5.809412 6.421636 6.219771  6.160157  5.889370  5.796417 
# calculate average values for each column
list_means <- apply(listdata, 2, mean)
list_means
# 0-100K    25K-125K 50-150K  75K-175K 100K-200K 125K-225K 150K-250K 
# 7.536545 6.763968 6.436773 6.235975 6.163197 5.912452 5.799420 

# values from -150K lowest variance probes (median/mean)
# 5.789078 / 5.790702

# write.csv(listdata,"TCGA_BR_KI_LU_TH_1000samplingsTestMSE_merged2.csv",row.names=F)

# test MSE boxplots=============================================================
# pdf("hannum_100Ksubsets_testRMSE_boxplots.pdf")
par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(3,9),
        #ylab=("Upper index position of 100,000 features in subset"),
        notch = T,
        horizontal=T,
        las=1, # style of axis labels set to horizontal
        names=c("100",
                "125",
                "150",
                "175",
                "200",
                "225",
                "250"
        ),
        xlab="RMSE"
        #main="1000 samplings 75/25 train/test RMSE for varying subsets of features"
)
title(ylab=expression(paste("Upper index position of 100,000 features in subset (10 "^"  3",")")), line=2.4)
# dev.off()

# hannum boxplots 1000samplings reversed===============================================
library(data.table)
histlist1 <- fread('hannum_allFeatures_testRMSE.csv')
# histlist1 <- as.numeric(histlist1)
histlist2 <- fread('hannum_minus25KtestMSE_reversed.csv')
# histlist2 <- as.numeric(histlist2)
histlist3 <- fread('hannum_minus50KtestMSE_reversed.csv')
histlist4 <- fread('hannum_minus75KtestMSE_reversed.csv')
histlist5 <- fread('hannum_minus100KtestMSE_reversed.csv')
histlist6 <- fread('hannum_minus125KtestMSE_reversed.csv')
histlist7 <- fread('hannum_minus150KtestMSE_reversed.csv')
histlist8 <- fread('hannum_minus175KtestMSE_reversed.csv')
histlist9 <- fread('hannum_minus200KtestMSE_reversed.csv')


# compute unpaired two-sample t-test
ttest <- t.test(histlist1,histlist2,var.equal=F) # x=all features, y=-150K features
ttest
# Welch Two Sample t-test
# 
# data:  histlist1 and histlist2
# t = 0.14598, df = 1997.4, p-value = 0.884
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4710187  0.5467769
# sample estimates:
#   mean of x mean of y 
# 34.80437  34.76649 

# wilcoxon rank sum test
# reject null hypothesis, alternative hypothesis of true location shift not equal to 0
# different medians / distribution
allFeatures_RMSE <- sqrt(histlist1)
minus25K_highVarProbes_RMSE <- sqrt(histlist2)
wilcox.test(allFeatures_RMSE,minus25K_highVarProbes_RMSE)
# Wilcoxon rank sum test with continuity correction
# 
# data:  allFeatures_RMSE and minus25K_highVarProbes_RMSE
# W = 503526, p-value = 0.7848
# alternative hypothesis: true location shift is not equal to 0

histlist1 <-stack(histlist1)
histlist1 <- histlist1[1]
histlist2 <-stack(histlist2)
histlist2 <- histlist2[1]
histlist3 <-stack(histlist3)
histlist3 <- histlist3[1]
histlist4 <-stack(histlist4)
histlist4 <- histlist4[1]
histlist5 <-stack(histlist5)
histlist5 <- histlist5[1]
histlist6 <-stack(histlist6)
histlist6 <- histlist6[1]
histlist7 <-stack(histlist7)
histlist7 <- histlist7[1]
histlist8 <-stack(histlist8)
histlist8 <- histlist8[1]
histlist9 <-stack(histlist9)
histlist9 <- histlist9[1]
# median(unlist(histlist1)) # 34.84726
# median(unlist(histlist2)) # 34.42873
# median(unlist(histlist3)) # 35.4503
# 
# mean(unlist(histlist1)) # 35.14379
# mean(unlist(histlist2)) # 34.66881
# mean(unlist(histlist3)) # 35.76858
# save all test MSE data into one data frame====================================
listdata <- data.frame(x1=histlist1,
                       x2=histlist2,
                       x3=histlist3,
                       x4=histlist4,
                       x5=histlist5,
                       x6=histlist6,
                       x7=histlist7,
                       x8=histlist8,
                       x9=histlist9)
# convert to RMSE
length(listdata)
listdata[9]
for (i in 1:9){
  listdata[i] = sqrt(listdata[i])
}
dim(listdata)
listdata[1:5,1:5]
median(listdata$values)
range(listdata$values.3) # 19.62670 50.71518 <--- -150K
range(listdata$values.2) # 21.40714 58.34211 <--- -100K
# calculate median values for each column
list_medians <- apply(listdata, 2, median)
list_medians
# values values.1 values.2 values.3 values.4 
# 34.53138 34.45778 34.28183 33.51342 36.57591
listdata$median[1:6]
listdata[1:10,1:6]
# write.csv(listdata,"TCGA_BR_KI_LU_TH_1000samplingsTestMSE_merged2.csv",row.names=F)

# test MSE boxplots=============================================================
# pdf("1000samplings_reversed_testMSE_boxplots.pdf")
par(mfrow=c(1,1))
boxplot(listdata,
        ylim=c(3,10),
        #ylab=expression(paste("Number of highest variance probes removed from subset (10"^"3*",")")),
        notch = T,
        horizontal=T,
        las=1, # style of axis labels set to horizontal
        names=c("0",
                "25",
                "50",
                "75",
                "100",
                "125",
                "150",
                "175",
                "200"),
        xlab="RMSE"
        #main="1000 samplings 75/25 train/test MSE for varying subsets of features"
)
title(ylab=expression(paste("Number of highest variance probes removed (10 "^"  3",")")), line=2.4)
# dev.off()