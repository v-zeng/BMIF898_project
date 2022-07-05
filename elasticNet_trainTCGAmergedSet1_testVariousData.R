# elastic net models with low variance probes removed (50K, 100K)
# load library
# sink("BR_KI_LU_TH_elasticNet_lowVariance50K100K_removed.txt") # ============= change
library(glmnet)
library(data.table)

#age calibration function and inverse age calibration function
ageC <- function(age) {
  if (age <= 20) { OUT <- log((age+1)/(21)) }
  if (age > 20)  { OUT <- (age-20)/(21) }
  OUT }

ageC.inverse <- function(DNAmAge) {
  if (DNAmAge <= 0) { OUT <- (21)*exp(DNAmAge) - 1 }
  if (DNAmAge > 0)  { OUT <- 20 + DNAmAge*(21) }
  OUT }

# train and test with all features==============================================
# ensure repeatable results
set.seed(3)
# load data
df <- fread('TCGA_BR_KI_LU_TH_mergedSD.csv')
df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
df[1:5,1:5]
range(df$Age) # 15 90
# age calibration
df$Age <- ageC(df$Age)

x <- as.matrix(df[,-1])
y <- as.matrix(df[,1])

# divide data into training and testing sets
train_index<-sample(1:nrow(df),.75*nrow(df)) # 75% data in training
x.train<-x[train_index,]
x.test<-x[-train_index,]
y.train<-y[train_index]
y.test<-y[-train_index]


# train and test minus bottom 150K features======================================
set.seed(3)
df3 <- df[,c(1,150002:ncol(df))]
x3 <- as.matrix(df3[,-1])
y3 <- as.matrix(df3[,1])
# divide data into training and testing sets
train_index3<-sample(1:nrow(df3),.75*nrow(df3)) # 75% data in training
x.train3<-x3[train_index3,]
x.test3<-x3[-train_index3,]


# use cv.glmnet to estimate model performance and find lambda.min===============
# elastic-net regression all features===========================================
set.seed(3)
cvfitAllFeatures <- cv.glmnet(x.train, y.train, alpha=0.5, type.measure="mse", nfolds=10)
print("All Features cvfit======================")
print(cvfitAllFeatures)
lambdaAllFeatures <- cvfitAllFeatures$lambda.min
lambdaAllFeatures # 0.03005172 / re-ran this and got 0.0331767 ???
set.seed(3)
elasticFit.all<- glmnet(x, y, alpha=0.5, lambda=lambdaAllFeatures)

# extract clock CpGs for elastic net fit all features =================================
coeffs <- coef(elasticFit.all, s = "lambda.min")
# extract coefficient variable names from glmnet into data frame
clock_CpGs_allFeatures <-data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x) # https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
dim(clock_CpGs_allFeatures) # 253   2 / re-ran 231 2
clock_CpGs_allFeatures[1:5,1:2]
# name coefficient
# 1 (Intercept)  -1.9703832
# 2  cg07655450  -1.1849976
# 3  cg19724632  -0.1184105
# 4  cg23999973  -1.0340840
# 5  cg01227006  -0.4632738
write.csv(clock_CpGs_allFeatures, "clock_CpGs_allFeatures_intercept.csv",row.names=F)

# check written file for clock CpGs
clock_CpGs_allFeatures <- fread("clock_CpGs_allFeatures_intercept.csv")
clock_CpGs_allFeatures[1:5,1:2]
clock_CpGs_allFeatures <- clock_CpGs_allFeatures[-1,] # remove intercept
write.csv(clock_CpGs_allFeatures,"clock_CpGs_elasticNet_allFeatures.csv",row.names=F)

# build clock using clock CpGs==================================================
clock_CpGs_allFeatures <- fread("clock_CpGs_elasticNet_allFeatures.csv")
head(clock_CpGs_allFeatures)
clock_CpGs_allFeatures_list <- clock_CpGs_allFeatures$name
length(clock_CpGs_allFeatures_list) # 252 clock CpGs
# keep clock CpGs in test set

df <- fread('TCGA_BR_KI_LU_TH_mergedSD.csv')
df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
df_clockCpGs_allFeatures <- df[,with(df,names(df) %in% clock_CpGs_allFeatures_list)]
df_clockCpGs_allFeatures$Age <- df$Age
dim(df_clockCpGs_allFeatures)
df_clockCpGs_allFeatures[1:5,1:5]
# move age to index 1
df_clockCpGs_allFeatures<-df_clockCpGs_allFeatures[,c(which(colnames(df_clockCpGs_allFeatures)=="Age"),which(colnames(df_clockCpGs_allFeatures)!="Age"))]
df_clockCpGs_allFeatures[1:5,1:5]
# Age cg11854493 cg07060551 cg23520574 cg16053920
# TCGA-05-5420  67  0.6482925 0.17710257  0.8200992  0.8576895
# TCGA-18-3417  65  0.6209843 0.15488661  0.8774021  0.9293321
# TCGA-18-4721  74  0.6298808 0.14357998  0.9164424  0.9282929
# TCGA-18-5592  57  0.6936497 0.15999863  0.8929274  0.9470991
# TCGA-18-5595  50  0.6667959 0.09682858  0.9250100  0.9375089
write.csv(df_clockCpGs_allFeatures,"TCGA_mergedSet1_clockCpGs_allFeatures.csv")
# fit clock allFeatures=====================================================================
library(data.table)
#age calibration function and inverse age calibration function
ageC <- function(age) {
  if (age <= 20) { OUT <- log((age+1)/(21)) }
  if (age > 20)  { OUT <- (age-20)/(21) }
  OUT }

ageC.inverse <- function(DNAmAge) {
  if (DNAmAge <= 0) { OUT <- (21)*exp(DNAmAge) - 1 }
  if (DNAmAge > 0)  { OUT <- 20 + DNAmAge*(21) }
  OUT }

df_clockCpGs_allFeatures <- fread("TCGA_mergedSet1_clockCpGs_allFeatures.csv")
df_clockCpGs_allFeatures <- data.frame(df_clockCpGs_allFeatures, row.names = 1)
dim(df_clockCpGs_allFeatures) # 431 253
df_clockCpGs_allFeatures[1:5,1:5]
# histogram and density of SD positions
# need to find index of CpGs relative to McIntyre probes
lines(density(df_clockCpGs_allFeatures),col=4,lwd=2,xlim=c(0,256801)) # histogram with density line example, use for # CpGs compared in clocks not this ================================


df_clockCpGs_allFeatures$Age <- ageC(df_clockCpGs_allFeatures$Age)
clock_allFeatures <- lm(Age~., data=df_clockCpGs_allFeatures) # build clock
summary(clock_allFeatures)
# count positive and negative coefficients for all features=======================
coefficient_list_allFeatures <- summary(clock_allFeatures)$coefficients[,"Estimate"]
length(coefficient_list_allFeatures) # 253
coefficient_list_allFeatures <- coefficient_list_allFeatures[-1] # drop intercept, No. CpGs == 252
sum(coefficient_list_allFeatures < 0) # 120
sum(coefficient_list_allFeatures > 0) # 132

# Residual standard error: 0.0839 on 178 degrees of freedom
# Multiple R-squared:  0.994,	Adjusted R-squared:  0.9855 
# F-statistic: 117.2 on 252 and 178 DF,  p-value: < 2.2e-16
# test clock 
set.seed(3)
test.data <- fread('TCGA_BR_KI_LU_TH_mergedSD.csv')
test.data <- data.frame(test.data, row.names = 1)
test.data_noAge <- test.data[,-1]
test.values <- ageC.inverse(predict(clock_allFeatures,test.data_noAge))
range(test.values) # 15.10274 91.15426
range(test.data$Age) # 15 90

# pearson correlation
pCorr <- cor(test.data$Age,test.values)
pCorr # 0.9969993

# ==============================================================================
### Histogram and density of Standard Deviation positions for clock CpGs
# find indices for multiple elements (CpGs) in 256801 features
# Note: fit clocks first before running this code
dev.new(width=20,height=20,noRStudioGD = TRUE)
par(mfrow=c(1,2),pin=c(3.5,2.5))

# allFeatures
df <- fread("TCGA_BR_KI_LU_TH_mergedSD.csv")
df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
df <- df[,-1] # exclude "Age" column
dim(df) # 431 256801
CpG_names_total <- names(df)# get feature names from all 256801 probes
CpG_names_allFeatures <- names(coefficient_list_allFeatures) # get allFeatures model CpG names
length(CpG_names_allFeatures) # 252
inds_all <- which(CpG_names_total %in% c(CpG_names_allFeatures)) # find indices of allFeatures probes relative to 256801 probes
inds_all[1:5]
sum(inds_all < 150000) # 155 probes under 150K SD position mark for allFeatures clock

hist(inds_all,breaks=40,xlim=c(0,256801),freq=F,xlab="Standard Deviation Position",main="") # try adding freq=F
lines(density(inds_all),
      lwd=2,
      col="red")
text(x=50000,y=8e-06,font=2,label="A")
# mtext("RMSE", side=1, line=2.2, cex=1, col="black")
# mtext("B", side=2, line=line, cex=cex, las=las)
# mtext("C", side=2, line=line, cex=cex, las=las)
# title(xlab=expression(paste("Dummy text")), line=2.4)
# -150K
CpG_names_minus150K <- names(coefficient_list) # get allFeatures model CpG names
length(CpG_names_minus150K) # 213
inds_minus150K <- which(CpG_names_total %in% c(CpG_names_minus150K)) # find indices of allFeatures probes relative to 256801 probes
inds_minus150K[1:5]
sum(inds_minus150K>=251000) # 152


hist(inds_minus150K,breaks=20,xlim=c(0,256801),freq=F,xlab="Standard Deviation Positon",main="",ylab="") # try adding freq=F
lines(density(inds_minus150K),
      lwd=2,
      col="red")
text(x=50000,y=1.35e-05,font=2,label="B")
# see https://stackoverflow.com/questions/67166083/density-curve-on-histogram-is-flat
# dev.new(width=20,height=20,noRStudioGD = TRUE)
d <- density(inds_all)
dx <- diff(d$x)[1]
sum(d$y)*dx # 1.000962
h <- hist(inds_all)
lines(x=d$x,y=max(h$counts)*d$y/dx)
#===============================================================================

# plot chronological age vs DNAm age for GSE37754 test
cAge <- test.data$Age
plot(cAge,test.values,
     main="TCGA BR_KI_LU_TH lm() train/test 75/25 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(20,100),
     ylim=c(20,100)
)
abline(a=0,b=1)
abline(lm(test.values~cAge),col=c("blue"))

# elastic-net regression using features minus bottom 150K=======================
set.seed(3)
cvfitMinus150K <- cv.glmnet(x.train3, y.train, alpha=0.5, type.measure="mse", nfolds=10)
print("Minus 150K cvfit=============================")
print(cvfitMinus150K)
lambdaMinus150K <- cvfitMinus150K$lambda.min
lambdaMinus150K # 0.02613746; lambda value for elastic net regression  ========= reran and got 0.01576162 ???
# fit elastic net model using lambda from -150K cvfit===========================
set.seed(3)
elasticFit.minus150K <- glmnet(x3, y3, alpha=0.5, lambda=lambdaMinus150K)

# extract clock CpGs for -150K elastic net fit =================================
coeffs <- coef(elasticFit.minus150K, s = "lambda.min")
# extract coefficient variable names from glmnet into data frame
clock_CpGs <-data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x) # https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
dim(clock_CpGs) # 214   2 / 304   2
clock_CpGs[1:5,1:2]
# name  coefficient
# 1 (Intercept)  0.005026551
# 2  cg11854493 -0.451412212
# 3  cg07060551  0.653973276
# 4  cg23520574  0.021639826
# 5  cg16053920  0.104686548
write.csv(clock_CpGs, "clock_CpGs_intercept.csv",row.names=F)

# check written file for clock CpGs
clock_CpGs <- fread("clock_CpGs_intercept.csv")
clock_CpGs[1:5,1:2]
clock_CpGs <- clock_CpGs[-1,] # remove intercept
write.csv(clock_CpGs,"clock_CpGs_elasticNet_minus150K.csv",row.names=F)

# build clock using clock CpGs==================================================
clock_CpGs <- fread("clock_CpGs_elasticNet_minus150K.csv")
head(clock_CpGs)
clock_CpGs_list <- clock_CpGs$name
length(clock_CpGs_list) # 213 clock CpGs
# keep clock CpGs in test set

df <- fread('TCGA_BR_KI_LU_TH_mergedSD.csv')
df <- data.frame(df, row.names = 1) # column '1' used for row names, column is removed automatically
df_clockCpGs <- df[,with(df,names(df) %in% clock_CpGs_list)]
df_clockCpGs$Age <- df$Age
dim(df_clockCpGs)
df_clockCpGs[1:5,1:5]
# move age to index 1
df_clockCpGs<-df_clockCpGs[,c(which(colnames(df_clockCpGs)=="Age"),which(colnames(df_clockCpGs)!="Age"))]
df_clockCpGs[1:5,1:5]
# Age cg11854493 cg07060551 cg23520574 cg16053920
# TCGA-05-5420  67  0.6482925 0.17710257  0.8200992  0.8576895
# TCGA-18-3417  65  0.6209843 0.15488661  0.8774021  0.9293321
# TCGA-18-4721  74  0.6298808 0.14357998  0.9164424  0.9282929
# TCGA-18-5592  57  0.6936497 0.15999863  0.8929274  0.9470991
# TCGA-18-5595  50  0.6667959 0.09682858  0.9250100  0.9375089
write.csv(df_clockCpGs,"TCGA_mergedSet1_clockCpGs.csv")
# write.csv(df_clockCpGs,"mergedSet1_clockCpGs_age.csv")
# fit clock ====================================================================
library(data.table)
#age calibration function and inverse age calibration function
ageC <- function(age) {
        if (age <= 20) { OUT <- log((age+1)/(21)) }
        if (age > 20)  { OUT <- (age-20)/(21) }
        OUT }

ageC.inverse <- function(DNAmAge) {
        if (DNAmAge <= 0) { OUT <- (21)*exp(DNAmAge) - 1 }
        if (DNAmAge > 0)  { OUT <- 20 + DNAmAge*(21) }
        OUT }

df_clockCpGs <- fread("TCGA_mergedSet1_clockCpGs.csv")
df_clockCpGs <- data.frame(df_clockCpGs, row.names = 1)
df_clockCpGs[1:5,1:5]
dim(df_clockCpGs)
# # train/test split
# set.seed(42)
# row.num <- sample(1:nrow(df_clockCpGs), 0.75*nrow(df_clockCpGs))
# train = df_clockCpGs[row.num,]
# test = df_clockCpGs[-row.num,]
# dim(train)
# dim(test)
# transform age in train
# train$Age <- ageC(train$Age)
# train[1:5,1:5]
# fit/train clock/model
df_clockCpGs$Age <- ageC(df_clockCpGs$Age)
clock <- lm(Age~., data=df_clockCpGs) # build clock
summary(clock)
# find negative and positvely correlated CpGs
# count positive and negative coefficients for all features=======================
coefficient_list <- summary(clock)$coefficients[,"Estimate"]
length(coefficient_list) # 214
coefficient_list[1:20]
coefficient_list <- coefficient_list[-1] # drop intercept, No. CpGs == 213
sum(coefficient_list < 0) # 101
sum(coefficient_list > 0) # 112

# Residual standard error: 0.09575 on 217 degrees of freedom
# Multiple R-squared:  0.9905,	Adjusted R-squared:  0.9811 
# F-statistic: 106.1 on 213 and 217 DF,  p-value: < 2.2e-16
# test clock 
set.seed(3)
test.data <- fread('TCGA_BR_KI_LU_TH_mergedSD.csv')
test.data <- data.frame(test.data, row.names = 1)
test.data_noAge <- test.data[,-1]
test.values <- ageC.inverse(predict(clock,test.data_noAge))
range(test.values) # 13.31227 91.82250
range(test.data$Age) # 15 90
# pearson correlation
pCorr <- cor(test.data$Age,test.values)
pCorr #0.9952314
# plot chronological age vs DNAm age for GSE37754 test
cAge <- test.data$Age
plot(cAge,test.values,
     main="TCGA BR_KI_LU_TH lm() train/test 75/25 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(20,100),
     ylim=c(20,100)
)
abline(a=0,b=1)
abline(lm(test.values~cAge),col=c("blue"))
# GSE101961 test clock_allFeatures===============================================================
set.seed(3)
df_test <- fread('GSE101961_breast_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 121 256802
df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 11.11665 62.10963
range(df_test.y) # 17 76
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 4.990012
# calculate median absolute deviation of residuals
MAD <- mad((df_test.y-fitted.test)) # median absolute deviation (MAD) quantifies variation
MAD # 5.377508
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 6.165854/ 6.457379
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.8667085
# r squared
rsquared <- cor(df_test.y,fitted.test)^2 # calculate adjusted r squared instead==========================
rsquared # 0.7511836

cAge <- df_test.y

# plot predicted age (DNAm age) vs chronological age
dev.new(width=20,height=20,noRStudioGD = TRUE)
par(mfrow=c(4,2), pin=c(3.5,1.5),oma=c(2,2,0,0))
# layout(matrix(c(1,2,3,4,5,6),byrow=T,ncol=2,nrow=3))
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(0,100),
     ylim=c(0,100)
)
abline(a=0,b=1)
text(x=2,y=95,font=2,label="A)")
text(x=12,y=85,font=1,label="RMSE: 6.457379")
mtext("All-features clock", side=1, line=4, cex=1, col="black") # run this on 7th plot; "line" changes how far from plot x-label is - use this on last two bottom plots
mtext("-150K clock", side=1, line=4, cex=1, col="black") # run this on 8th plot; "line" changes how far from plot x-label is - use this on last two bottom plots
mtext("DNAm age", side=2, line=0, cex=1.5, col="black",outer=T)
mtext("Chronological age", side=1, line=0, cex=1.5, col="black",outer=T)
# GSE101961 test ===============================================================
set.seed(3)
df_test <- fread('GSE101961_breast_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 121 256802
df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 17.32979 59.07840
range(df_test.y) # 17 76
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 4.301505
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 4.552335
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 6.165854
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.8609768
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.7412811

cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(0,100),
     ylim=c(0,100)
)
abline(a=0,b=1)
text(x=2,y=95,font=2,label="B)")
text(x=12,y=85,font=1,label="RMSE: 6.165854")
# GSE79100 test data set=========================================================
df_test <- fread('GSE79100_kidney_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1) # column '1' used for row names, column is removed automatically
df_test[1:5,1:5]
df_test.x <- df_test[,-1]
df_test.y <- df_test[,1]
# predict GSE79100
fitted.GSE79100 <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg01330954' not found 

# df_test_clockCpGs <- df_test_noAge[,with(df_test_noAge, names(df_test_noAge) %in% clock_CpGs$name)]
# df_test_clockCpGs[1:5,1:5]

# GSE88883 test=================================================================
df_test_GSE88883 <- fread('GSE88883_breast_mergedClean.csv')
df_test_GSE88883<- data.frame(df_test_GSE88883, row.names = 1) 
df_test_GSE88883[1:5,1:5]
df_test.x1 <- df_test_GSE88883[,-1]
df_test.y1 <- df_test_GSE88883[,1]
fitted.GSE88883 <- ageC.inverse(predict(clock,df_test.x1)) #  Error in eval(predvars, data, env) : object 'cg15394630' not found 

# GSE37754 test=================================================================
df_test_GSE37754 <-fread('GSE37754_breast_mergedClean.csv') # study identified DNAm patterns associated with tumor subtypes
df_test_GSE37754 <- data.frame(df_test_GSE37754, row.names = 1) 
df_test_GSE37754[1:5,1:5]
df_test.x2 <- df_test_GSE37754[,-1]
df_test.y2 <- df_test_GSE37754[,1]
fitted.GSE37754 <- ageC.inverse(predict(clock,df_test.x2))
range(fitted.GSE37754) # -0.9999646 2108.0223608
range(df_test.y2) # 30 93
# plot chronological age vs DNAm age for GSE37754 test
cAge <- df_test_GSE37754$Age
plot(cAge,fitted.GSE37754,
     main="GSE37754 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,100),
     ylim=c(0,100)
     )
abline(a=0,b=1)

# GSE66313 test=================================================================
df_test_GSE66313 <-fread('GSE66313_breast_mergedClean.csv')
df_test_GSE66313 <- data.frame(df_test_GSE66313, row.names = 1) 
df_test.x3 <- df_test_GSE66313[,-1]
df_test.y3 <- df_test_GSE66313[,1]
fitted.GSE37754 <- ageC.inverse(predict(clock,df_test_GSE66313)) #  Error in eval(predvars, data, env) : object 'cg23520574' not found 
range(fitted.GSE37754)
# plot chronological age vs DNAm age for GSE37754 test
cAge <- df_test_GSE37754$Age
plot(cAge,fitted.GSE37754,
     main="GSE66313 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,100),
     ylim=c(0,100)
)
abline(a=0,b=1)
# mergedSet2 test===============================================================
df_test <-fread('TCGA_mergedSet2.csv')
df_test <- data.frame(df_test, row.names = 1) 
df_test.x <- df_test[,-1]
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 26.1682 118.1940
range(df_test.y) # 20 90

# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 14.92912 / 13.68366 (clock, clock_allFeatures)
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test)) # median absolute deviation (MAD) quantifies variation, how spread out set of data is 
MAD # 16.24691 / 17.25599 (clock, clock_allFeatures)
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 18.55652 / 16.81055
# plot chronological age vs DNAm age for GSE37754 test
cAge <- df_test.y
plot(cAge,fitted.test,
     main="TCGA_mergedSet2 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,120),
     ylim=c(0,120)
)
abline(a=0,b=1)
# PRAD test===============================================================
df_test <-fread('PRAD_probesAge_mergedCommon.csv')
dim(df_test) # 50 256803
df_test <- data.frame(df_test, row.names = 1) 
df_test.x <- df_test[,-1]
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(fitted.test) # 26.1682 100.0280 / 32.82242 90.79679 (clock, clock_allFeatures)
range(df_test.y) # 44 72

# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 9.146812 / 9.233554 (clock, clock_allFeatures)
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test)) # median absolute deviation (MAD) quantifies variation, how spread out set of data is 
MAD # 9.945353 / 14.27943 (clock, clock_allFeatures)
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 11.8528 / 11.26036
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.4452714 / 0.3694121
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.1982666 / 0.1364653
# plot chronological age vs DNAm age for GSE37754 test
cAge <- df_test.y
plot(cAge,fitted.test,
     main="clock all features PRAD n=50 RMSE: 11.26, MAE: 9.23",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,100),
     ylim=c(0,100)
)
abline(a=0,b=1)
# LIHC test===============================================================
df_test <-fread('LIHC_probesAge_mergedCommon.csv')
df_test <- data.frame(df_test, row.names = 1) 
dim(df_test) # 50
df_test.x <- df_test[,-1]
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 32.26881 85.87026 / 33.53429 107.50312 (clock, clock_allFeatures)
range(df_test.y) # 20 81

# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 8.506601 / 14.34994 (clock, clock_allFeatures)
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test)) # median absolute deviation (MAD) quantifies variation, how spread out set of data is 
MAD # 8.905077 / 12.84643 (clock, clock_allFeatures)
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 10.89658 / 17.47927
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.8372636 / 0.745903
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.7010103 / 0.5563713
# plot chronological age vs DNAm age for GSE37754 test
cAge <- df_test.y
plot(cAge,fitted.test,
     main="-150Kk clock LIHC (n=50) RMSE: 8.51 MAE: 10.9",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,110),
     ylim=c(0,110)
)
abline(a=0,b=1)
# HNSC test===============================================================
df_test <-fread('HNSC_probesAge_mergedCommon.csv')
df_test <- data.frame(df_test, row.names = 1) 
df_test.x <- df_test[,-1]
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 26.1682 100.0280 / 32.82242 90.79679 (clock, clock_allFeatures)
range(df_test.y) # 44 72

# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 22.79149 / 19.63644 (clock, clock_allFeatures)
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test)) # median absolute deviation (MAD) quantifies variation, how spread out set of data is 
MAD # 14.96063 / 16.29805 (clock, clock_allFeatures)
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 26.30662 / 22.55854
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.4722689 / 0.4631506
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.2230379 / 0.2145085
# plot chronological age vs DNAm age for GSE37754 test
cAge <- df_test.y
plot(cAge,fitted.test,
     main="allFeatures HNSC Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,120),
     ylim=c(0,120)
)
abline(a=0,b=1)
# GSE59157 test=================================================================
set.seed(3)
df_test <- fread('GSE59157_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)

df_test[1:5,1:5]
df_normal <- df_test[grepl('normal', df_test$Sample_type),]
# "normal kidney"    "nephrogenic rest" "Wilms tumour"
dim(df_normal) # 35 256803
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 2.921865 19.193066
# range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 4.513538
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 3.720575
# RMSE is SD of residuals (prediction errors); measure of how spread out residuals are/ how concentrated data is around line of best fit
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 5.729042 ; new allFeatures 6.760274
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 #
rsquared # 0.08021054


cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(0,27),
     ylim=c(0,27)
)
abline(a=0,b=1)
# allFeatures GSE59157 test=================================================================
set.seed(3)
df_test <- fread('GSE59157_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
unique(df_test$Sample_type) # "normal kidney"    "nephrogenic rest" "Wilms tumour" 
df_normal <- df_test[grepl('Wilms tumour', df_test$Sample_type),]
dim(df_normal) # 35 256803
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(df_test.y) # 0.8333333 12.0000000 
# age_in_months <- df_test.y*12
# range(age_in_months) # 10 144
range(fitted.test) # 2.921865 19.193066 ; wilms 1.178221 26.381240
# range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error; measures accuracy for continuous variables
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 5.574355
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 1.73556
# rmse measures average magnitude of the error; use together with MAE to diagnose variation in errors; greater the difference = greater variance in individual errors

RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 6.763686
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 # use adjusted r squared
rsquared # 0.06827915

cAge <- df_test.y
plot(cAge,fitted.test,
     main="allFeatures clock GSE59157 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,25),
     ylim=c(0,25)
)
abline(a=0,b=1)

# GSE74214 test=================================================================
set.seed(3)
df_test <- fread('GSE74214_breast_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
dim(df_test) #18 256802
df_test[1:5,1:5]

df_test.x <- df_normal[,-1] # not age
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 2.921865 19.193066
# range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error; measures accuracy for continuous variables
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 5.574355
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 1.73556
# rmse measures average magnitude of the error; use together with MAE to diagnose variation in errors; greater the difference = greater variance in individual errors

RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 6.763686
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 # use adjusted r squared
rsquared # 0.06827915

cAge <- df_test.y
plot(cAge,fitted.test,
     main="allFeatures clock GSE85566 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,25),
     ylim=c(0,25)
)
abline(a=0,b=1)
# GSE156932 test================================================================
set.seed(3)
df_test <- fread('GSE156932_kidney_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
unique(df_test$Tissue)
# [1] "oncocytoma"                        "normal kidney parenchyma"          "renal cell carcinoma-chromophobe" 
# [4] "clear cell renal cell carcinoma"   "hybrid oncocytic/chromophobe type" "hybrid oncocytic renal neoplasm"
df_test[1:5,1:5]
df_normal <- df_test[grepl('oncocytoma', df_test$Tissue),]
dim(df_normal) # 35 256803
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(fitted.test) # 23.76927 81.34013 / 26.46616 80.59019
range(df_test.y) # 34 79
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 6.495504/5.775443
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 8.176019/4.273967
# RMSE is SD of residuals (prediction errors); measure of how spread out residuals are/ how concentrated data is around line of best fit
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 8.028063/6.749859
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 #
rsquared # 0.7674505/0.8173732


cAge <- df_test.y
plot(cAge,fitted.test,
     main="Clock all features GSE156932 oncocytoma kidney (n=12) RMSE: 23.93 MAE: 19.51",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,130),
     ylim=c(0,130)
)
abline(a=0,b=1)

# KIRC tumors test =============================================================
library(data.table)
set.seed(3)
df_test <- fread('KIRCtumors_probesAge_merged.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 319 214

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 26.71516 138.36479
range(df_test.y) # 26 90
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 15.47573
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 15.3601
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 22.24626
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.4304003
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.1852445

cAge <- df_test.y
plot(cAge,fitted.test,
     main="KIRC tumor samples Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,100),
     ylim=c(0,100)
)
abline(a=0,b=1,col=c("red"))
abline(lm(fitted.test~df_test.y),col=c("blue")) # linear interpolation of fitted data

hist(fitted.test-df_test.y) # histogram of tumor samples' DNAm Age - Chronological Age
lines(density(fitted.test-df_test.y),col=4,lwd=2) # histogram with density line example, use for # CpGs compared in clocks not this
# KIRC tumors test =============================================================
library(data.table)
set.seed(3)
df_test <- fread('KIRCtumors_allFeaturesClock_probesAge_merged.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 319 253
df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x)) 
range(fitted.test) # -17.05909 155.08300
range(df_test.y) # 26 90
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 18.77804
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 20.55812
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 22.24626; 26.18829 (all)
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.3815521
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.145582

cAge <- df_test.y
plot(cAge,fitted.test,
     main="KIRC tumor samples Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,100),
     ylim=c(0,100)
)
abline(a=0,b=1,col=c("red"))
abline(lm(fitted.test~df_test.y),col=c("blue")) # linear interpolation of fitted data
# BRCA tumors test =============================================================
library(data.table)
set.seed(3)
df_test <- fread('BRCAtumors_probesAge_merged.csv') # BRCAtumors_probesAge_merged.csv
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 756 214

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg16053920' not found
range(fitted.test) # -20.3700 150.6424
range(df_test.y) # 26 90
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 16.55597
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 18.8321
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 22.10403
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.3887047
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.1510914

cAge <- df_test.y
plot(cAge,fitted.test,
     main="BRCA tumor samples Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(-21,151),
     ylim=c(-21,151)
)
abline(a=0,b=1,col=c("red"))
# abline(lm(fitted.test~df_test.y),col=c("blue")) # linear interpolation of fitted data
# BRCA tumors test =============================================================
library(data.table)
set.seed(3)
df_test <- fread('BRCAtumors_allFeaturesClock_probesAge_merged.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 756 253

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(fitted.test) # -37.50926 174.13549
range(df_test.y) # 26 90
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 22.75372
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 24.583
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 30.12031
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.2905858
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.08444012

cAge <- df_test.y
plot(cAge,fitted.test,
     main="allFeatures BRCA tumor samples Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(-21,151),
     ylim=c(-21,151)
)
abline(a=0,b=1,col=c("red"))
# abline(lm(fitted.test~df_test.y),col=c("blue")) # linear interpolation of fitted data
# KIRP tumors test =============================================================
library(data.table)
set.seed(3)
df_test <- fread('KIRPtumors_probesAge_merged.csv') # this has been imputed ---> not imputed, sampels with NA omitted
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 271 214

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # -10.79584 118.10849
range(df_test.y) # 28 88
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 14.23804
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 16.38601
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 18.43928
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.4184854
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.17513

cAge <- df_test.y
plot(cAge,fitted.test,
     main="KIRP tumor samples Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(-11,126),
     ylim=c(-11,126)
)
abline(a=0,b=1,col=c("red"))
# abline(lm(fitted.test~df_test.y),col=c("blue")) # linear interpolation of fitted data
# KIRP tumors test all features clock =============================================================
library(data.table)
set.seed(3)
df_test <- fread('KIRPtumors_allFeaturesClock_probesAge_merged.csv') # this has been imputed
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 271 253

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(fitted.test) # 0.6350328 125.4287608
range(df_test.y) # 28 88
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 14.23804 / 14.17823
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 16.38601 / 16.1052
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 18.43928 / 18.98291
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.4184854 / 0.3393603
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.17513 / 0.1151654

cAge <- df_test.y
plot(cAge,fitted.test,
     main="KIRP tumor samples Chronological Age vs DNAmAge allFeatures Clock",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(-11,126),
     ylim=c(-11,126)
)
abline(a=0,b=1,col=c("red"))
# abline(lm(fitted.test~df_test.y),col=c("blue")) # linear interpolation of fitted data
# LUAD tumors test =============================================================
library(data.table)
set.seed(3)
df_test <- fread('LUADtumors_probesAge_merged.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 275 262724

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg12876594' not found

range(fitted.test) # 17.32979 59.07840
range(df_test.y) # 17 76
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 4.301505
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 4.552335
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 6.165854
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.8609768
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.7412811

cAge <- df_test.y
plot(cAge,fitted.test,
     main="GSE101961 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,100),
     ylim=c(0,100)
)
abline(a=0,b=1)

# LUSC tumors test =============================================================
library(data.table)
set.seed(3)
df_test <- fread('LUSCtumors_probesAge_merged.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 275 262724

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))

# GSE97466 thyroid test =============================================================
set.seed(3)
df_test <- fread('GSE97466_thyroid_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
df_normal <- df_test[grepl('non-neoplastic', df_test$Tissue),]
dim(df_normal) # 50 252049
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg15394630' not found
range(fitted.test) # 2.921865 19.193066
# range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error; measures accuracy for continuous variables
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 5.574355
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 1.73556
# rmse measures average magnitude of the error; use together with MAE to diagnose variation in errors; greater the difference = greater variance in individual errors

RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 6.763686
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 # use adjusted r squared
rsquared # 0.06827915

cAge <- df_test.y
plot(cAge,fitted.test,
     main="allFeatures clock GSE59157 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,25),
     ylim=c(0,25)
)
abline(a=0,b=1)

# GSE97466 thyroid test =============================================================
set.seed(3)
df_test <- fread('GSE86961_thyroid_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
df_normal <- df_test[grepl('non-neoplastic', df_test$Tissue),]
dim(df_normal) # 41 252494
df_test.x <- df_normal[,-c(1,2)] # not age and tissue type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg15394630' not found
range(fitted.test) # 2.921865 19.193066
# range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error; measures accuracy for continuous variables
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 5.574355
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 1.73556
# rmse measures average magnitude of the error; use together with MAE to diagnose variation in errors; greater the difference = greater variance in individual errors

RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 6.763686
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 # use adjusted r squared
rsquared # 0.06827915

cAge <- df_test.y
plot(cAge,fitted.test,
     main="allFeatures clock GSE59157 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(0,25),
     ylim=c(0,25)
)
abline(a=0,b=1)

# GSE83842 lung test ===========================================================
set.seed(3)
# df_test <- fread('GSE83842_lung_mergedClean.csv')
df_test <- df_mergedTest
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
df_normal <- df_test[grepl('non-neoplastic', df_test$Tissue),]
dim(df_normal) # 50 252049
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg05533001' not found

# GSE67919 breast test ===========================================================
set.seed(3)
# df_test <- fread('GSE83842_lung_mergedClean.csv')
df_test <- df_mergedTest
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
df_normal <- df_test[grepl('mastectomy', df_test$Tissue),]
dim(df_normal) # 50 252049
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
df_test.x[1:5,1:5]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg14952449' not found

# GSE113501 kidney test ===========================================================
set.seed(3)
# df_test <- fread('GSE83842_lung_mergedClean.csv')
df_test <- df_mergedTest
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
df_normal <- df_test[grepl('ccRCC', df_test$Tissue),]
dim(df_normal) # 50 252049
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
df_test.x[1:5,1:5]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg15394630' not found

# GSE39279 lung test ===========================================================
set.seed(3)
# df_test <- fread('GSE83842_lung_mergedClean.csv')
df_test <- df_mergedTest
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
df_normal <- df_test[grepl('adeno', df_test$Tissue),]
dim(df_normal) # 50 252049
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
df_test.x[1:5,1:5]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg14952449' not found

# GSE88883 lung test ===========================================================
set.seed(3)
# df_test <- fread('GSE83842_lung_mergedClean.csv')
df_test <- df_mergedTest
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
df_normal <- df_test[grepl('breast', df_test$Tissue),]
dim(df_normal) # 50 252049
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
df_test.x[1:5,1:5]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg15394630' not found

# GSE85566 lung test ===========================================================
set.seed(3)
# df_test <- fread('GSE83842_lung_mergedClean.csv')
df_test <- df_mergedTest
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
df_normal <- df_test[grepl('Asthma', df_test$Tissue),]
dim(df_normal) # 50 252049
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
df_test.x[1:5,1:5]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg14952449' not found

# GSE70977 buccal (oral/pharyngeal) test ===========================================================
set.seed(3)
# df_test <- fread('GSE83842_lung_mergedClean.csv')
df_test <- df_mergedTest
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
df_normal <- df_test[grepl('Control', df_test$Tissue),]
dim(df_normal) # 50 252049
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
df_test.x[1:5,1:5]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg04606672' not found


# GSE79100 kidney normal test ===========================================================
library(data.table)
set.seed(3)
df_test <- fread('GSE79100_kidney_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 29 256587

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x)) # Error in eval(predvars, data, env) : object 'cg12876594' not found

range(fitted.test) # -14.35352  75.69424
range(df_test.y) # 17 76
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 5.395069
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 3.783796
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 9.087297
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.875151
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.7658893

cAge <- df_test.y
plot(cAge,fitted.test,
     main="GSE79100 Chronological Age vs DNAmAge",
     xlab="Chronological Age",
     ylab="DNAm Age",
     xlim=c(-15,76),
     ylim=c(-15,76)
)
abline(a=0,b=1)
# check overlapping probes from both clocks (-150K, allFeatures)================
library(data.table)
df_clockCpGs <- fread("TCGA_mergedSet1_clockCpGs.csv")
df_clockCpGs <- data.frame(df_clockCpGs, row.names = 1)
range(df_clockCpGs$Age)
hist(unlist()) # fit for Wilms vs trained scales
df_clockCpGs[1:5,1:5]
df_clockCpGs <- df_clockCpGs[,-1] # drop age
df_clockCpGs_allFeatures <- fread("TCGA_mergedSet1_clockCpGs_allFeatures.csv")
df_clockCpGs_allFeatures <- data.frame(df_clockCpGs_allFeatures, row.names = 1)
df_clockCpGs_allFeatures <- df_clockCpGs_allFeatures[,-1]
matchingCpGs <- intersect(names(df_clockCpGs),names(df_clockCpGs_allFeatures))
length(matchingCpGs) #81========= try making clock using only these or the non-overlapping?
df_mergedSet1 <- fread("TCGA_BR_KI_LU_TH_mergedSD.csv")
idx_matchingCpGs <- match(matchingCpGs,names(df_mergedSet1))

# set.seed(3)
# elasticFit.allFeatures <- glmnet(x, y, alpha=0.5, lambda=lambdaAllFeatures)
# resultsNormal_allFeatures <- predict(elasticFit.allFeatures,x.allFeatures,s=lambdaAllFeatures)
# normalMSE_allFeatures <- mean((resultsNormal_allFeatures-y.testSet)^2)
# normalMSE_allFeatures # 223.9707
# 
# 
# 
# range(resultsNormal_allFeatures) # 28.1137 105.6977
# # range(y.testSet) # 20 90
# # lowerRange <- min(y.normal_range)
# # upperRange <- max(y.normal_range)
# par(mfrow=c(2,2)) # split display screen into separate panels
# plot(resultsNormal_allFeatures,y.testSet,
#      main="All Features model predictions (lambda.min)",
#      xlab="Predicted value",
#      ylab="Actual value",
#      xlim=c(20,110),
#      ylim=c(20,110)
# )
# abline(a=0,b=1)
# 
# hist(resultsNormal_allFeatures,breaks=20,
#      xlab="Age",
#      xlim=c(20,110),
#      main="Distribution of All Features GSE79100 predictions")
# 
# hist(y.testSet,breaks = 10, xlab="Age in years",
#      xlim=c(20,110),
#      main="Distribution of ages for GSE79100 samples (n=222)")
# 
# 
# # fit and predict using Minus 150K model =======================================
# set.seed(3)
# elasticFit.minus150K <- glmnet(x3, y3, alpha=0.5, lambda=lambdaMinus150K)
# resultsNormal_minus150K <- predict(elasticFit.minus150K,x.150K,s=lambdaMinus150K,type="response")
# normalMSE_minus150K <- mean((resultsNormal_minus150K-y.testSet)^2)
# normalMSE_minus150K # 230.095
# 
# range(resultsNormal_minus150K) # 25.76657 103.49536
# 
# # range(y.testSet) # 20 90
# # lowerRange <- min(y.normal_range)
# # upperRange <- max(y.normal_range)
# par(mfrow=c(2,2)) # split display screen into separate panels
# plot(resultsNormal_minus150K,y.testSet,
#      main="Minus 150K model TCGAset2 predictions (lambda.min)",
#      xlab="Predicted value",
#      ylab="Actual value",
#      xlim=c(20,110),
#      ylim=c(20,110)
# )
# abline(a=0,b=1)
# 
# hist(resultsNormal_minus150K,breaks=20,
#      xlab="Age",
#      xlim=c(20,110),
#      main="Distribution of Minus 150K model GSE79100 predictions")
# 
# hist(y.testSet,breaks = 10, xlab="Age in years",
#      xlim=c(20,110),
#      main="Distribution of ages for GSE79100 samples (n=31)")

