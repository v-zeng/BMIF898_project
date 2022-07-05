library(data.table)
### requires fitting of allFeatures and -150K models before running code below

# GSE101961 test clock_allFeatures===============================================================
set.seed(3)
df_test <- fread('GSE101961_breast_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 121 256802
df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
age_accel <- df_test.y-fitted.test
# hist(age_accel,breaks=50)
mean(abs(age_accel)) # 4.990012
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

# GSE101961 test ===============================================================
set.seed(3)
df_test <- fread('GSE101961_breast_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 121 256802
df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
age_accel <- df_test.y-fitted.test
# hist(age_accel,breaks=50)
mean(abs(age_accel)) # 4.301505
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


# allFeatures GSE59157 test=================================================================
set.seed(3)
df_test <- fread('GSE59157_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
unique(df_test$Sample_type) # "normal kidney"    "nephrogenic rest" "Wilms tumour" 
df_normal <- df_test[grepl('normal', df_test$Sample_type),]
dim(df_normal) # 35 256803
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(df_test.y) # 0.8333333 12.0000000 
age_in_months <- df_test.y*12
range(age_in_months) # 10 144
range(fitted.test) # 3.840379 19.602462
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
     main="",
     xlab="",
     ylab="",
     xlim=c(0,35),
     ylim=c(0,35)
)
abline(a=0,b=1)
text(x=1,y=33,font=2,label="C)")
text(x=4.5,y=30,font=1,label="RMSE: 6.763686")
# GSE59157 test=================================================================
set.seed(3)
df_test <- fread('GSE59157_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
unique(df_test$Sample_type) # "normal kidney"    "nephrogenic rest" "Wilms tumour" 
df_normal <- df_test[grepl('normal', df_test$Sample_type),]
# "normal kidney"    "nephrogenic rest" "Wilms tumour"
dim(df_normal) # 35 256803
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 2.921865 19.193066
range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 4.513538
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 3.720575
# RMSE is SD of residuals (prediction errors); measure of how spread out residuals are/ how concentrated data is around line of best fit
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 5.729042
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 #
rsquared # 0.08021054


cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(0,35),
     ylim=c(0,35)
)
abline(a=0,b=1)
text(x=1,y=33,font=2,label="D)")
text(x=4.5,y=30,font=1,label="RMSE: 5.729042")

# allFeatures GSE59157 test=================================================================
set.seed(3)
df_test <- fread('GSE59157_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
unique(df_test$Sample_type) # "normal kidney"    "nephrogenic rest" "Wilms tumour" 
df_normal <- df_test[grepl('nephrogenic rest', df_test$Sample_type),]
dim(df_normal) # 35 256803
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(df_test.y) # 0.8333333 12.0000000 
age_in_months <- df_test.y*12
range(age_in_months) # 10 144
range(fitted.test) # 3.710054 18.909682
# range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error; measures accuracy for continuous variables
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 5.574355
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 1.73556
# rmse measures average magnitude of the error; use together with MAE to diagnose variation in errors; greater the difference = greater variance in individual errors

RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 6.309364
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 # use adjusted r squared
rsquared # 0.06827915

cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(0,35),
     ylim=c(0,35)
)
abline(a=0,b=1)
text(x=1,y=33,font=2,label="E)")
text(x=4.5,y=30,font=1,label="RMSE: 6.309364")
# GSE59157 test=================================================================
set.seed(3)
df_test <- fread('GSE59157_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
unique(df_test$Sample_type) # "normal kidney"    "nephrogenic rest" "Wilms tumour" 
df_normal <- df_test[grepl('nephrogenic rest', df_test$Sample_type),]
# "normal kidney"    "nephrogenic rest" "Wilms tumour"
dim(df_normal) # 35 256803
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 3.128536 19.588067
range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 4.513538
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 3.720575
# RMSE is SD of residuals (prediction errors); measure of how spread out residuals are/ how concentrated data is around line of best fit
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 7.60057
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 #
rsquared # 0.08021054


cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(0,35),
     ylim=c(0,35)
)
abline(a=0,b=1)
text(x=1,y=33,font=2,label="F)")
text(x=4.5,y=30,font=1,label="RMSE: 7.60057")
# allFeatures GSE59157 test=================================================================
set.seed(3)
df_test <- fread('GSE59157_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
unique(df_test$Sample_type) # "normal kidney"    "nephrogenic rest" "Wilms tumour" 
df_normal <- df_test[grepl('Wilms tumour', df_test$Sample_type),]
dim(df_normal) # 36 256803
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(df_test.y) # 0.8333333 12.0000000 
age_in_months <- df_test.y*12
range(age_in_months) # 10 144
range(fitted.test) # 1.178221 26.381240
# range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error; measures accuracy for continuous variables
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 5.574355
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 1.73556
# rmse measures average magnitude of the error; use together with MAE to diagnose variation in errors; greater the difference = greater variance in individual errors

RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 6.694249
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 # use adjusted r squared
rsquared # 0.06827915

cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(0,35),
     ylim=c(0,35)
)
abline(a=0,b=1)
text(x=1,y=33,font=2,label="G)")
text(x=4.5,y=30,font=1,label="RMSE: 6.694249")
mtext("All-features clock", side=1, line=4, cex=1, col="black") # run this on 7th plot; "line" changes how far from plot x-label is - use this on last two bottom plots

# GSE59157 test=================================================================
set.seed(3)
df_test <- fread('GSE59157_mergedClean.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
unique(df_test$Sample_type) # "normal kidney"    "nephrogenic rest" "Wilms tumour" 
df_normal <- df_test[grepl('Wilms tumour', df_test$Sample_type),]
# "normal kidney"    "nephrogenic rest" "Wilms tumour"
dim(df_normal) # 36 256803
df_test.x <- df_normal[,-c(1,2)] # not age and sample type
df_test.y <- df_normal[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 1.214401 25.801265
range(df_test.y) # 0.8333333 12.0000000
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 4.513538
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 3.720575
# RMSE is SD of residuals (prediction errors); measure of how spread out residuals are/ how concentrated data is around line of best fit
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 8.156166
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2 #
rsquared # 0.08021054


cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(0,35),
     ylim=c(0,35)
)
abline(a=0,b=1)
text(x=1,y=33,font=2,label="H)")
text(x=4.5,y=30,font=1,label="RMSE: 8.156166")
mtext("-150K clock", side=1, line=4, cex=1, col="black") # run this on 8th plot; "line" changes how far from plot x-label is - use this on last two bottom plots
mtext("DNAm age", side=2, line=0, cex=1.5, col="black",outer=T)
mtext("Chronological age", side=1, line=0, cex=1.5, col="black",outer=T)

### model performances on TCGA tumour data sets=================================
# KIRC tumors test =============================================================
library(data.table)
set.seed(3)
df_test <- fread('KIRCtumors_allFeaturesClock_probesAge_merged.csv') #KIRCtumors_allFeaturesClock_probesAge_merged.csv
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 319 214

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
age_accel <- fitted.test-df_test.y
dev.new(width=20,height=20,noRStudioGD = TRUE)
par(mfrow=c(3,2), pin=c(3,2),oma=c(2,2,0,0))
range(age_accel) # -61.05909  78.08300
hist(age_accel,breaks=50,main="Kirc all-features",xlim=c(-100,100))
mean(abs(age_accel)) # 16.55597
range(fitted.test) # -17.05909 155.08300
range(df_test.y) # 26 90
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 15.47573
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 15.3601
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 26.18829
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.4304003
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.1852445

cAge <- df_test.y
dev.new(width=20,height=20,noRStudioGD = TRUE)
par(mfrow=c(3,2), pin=c(4,2),oma=c(2,2,0,0))
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(-20,160),
     ylim=c(-20,160)
)
abline(a=0,b=1)
text(x=-12,y=149,font=2,label="A)")
text(x=3,y=139,font=1,label="RMSE: 26.18829")
# KIRC tumors test =============================================================
library(data.table)
set.seed(3)
df_test <- fread('KIRCtumors_probesAge_merged.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 319 253
df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
age_accel <- fitted.test-df_test.y
range(age_accel) #-29.31456  80.73565
hist(age_accel,breaks=50,main="Kirc -150K",xlim=c(-100,100))
mean(abs(age_accel)) # 16.55597
range(fitted.test) # 26.71516 138.36479
range(df_test.y) # 26 90
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 18.77804
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 20.55812
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 22.24626
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.3815521
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.145582

cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(-20,160),
     ylim=c(-20,160)
)
abline(a=0,b=1)
text(x=-12,y=149,font=2,label="B)")
text(x=3,y=139,font=1,label="RMSE: 22.24626")
# KIRP tumors test allFeatures==================================================
library(data.table)
set.seed(3)
df_test <- fread('KIRPtumors_allFeaturesClock_probesAge_merged.csv') # this has been imputed ---> not imputed, sampels with NA omitted
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 271 214

df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
age_accel <- fitted.test-df_test.y
range(age_accel) # -47.69574  83.80626
hist(age_accel,breaks=50,main="Kirp all-features",xlim=c(-100,100))
mean(abs(age_accel)) # 16.55597
range(fitted.test) # 0.6350328 125.4287608
range(df_test.y) # 28 88
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 14.23804
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 16.38601
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 18.98291
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.4184854
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.17513

cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(-11,126),
     ylim=c(-11,126)
)
abline(a=0,b=1)
text(x=-5,y=115,font=2,label="C)")
text(x=7,y=106,font=1,label="RMSE: 18.98291")
# abline(lm(fitted.test~df_test.y),col=c("blue")) # linear interpolation of fitted data
# KIRP tumors test==============================================================
library(data.table)
set.seed(3)
df_test <- fread('KIRPtumors_probesAge_merged.csv') # this has been imputed
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 271 253
df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
age_accel <- fitted.test-df_test.y
range(age_accel) # -70.79584  51.90077
hist(age_accel,breaks=50,main="KIRP -150K",xlim=c(-100,100))
mean(abs(age_accel)) # 16.55597
range(fitted.test) # -10.79584 118.10849
range(df_test.y) # 28 88
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 14.23804 / 14.17823
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 16.38601 / 16.1052
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 18.43928
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.4184854 / 0.3393603
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.17513 / 0.1151654

cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(-11,126),
     ylim=c(-11,126)
)
abline(a=0,b=1)
text(x=-5,y=115,font=2,label="D)")
text(x=7,y=106,font=1,label="RMSE: 18.43928")
# BRCA tumors test allFeatures=============================================================
library(data.table)
set.seed(3)
df_test <- fread('BRCAtumors_allFeaturesClock_probesAge_merged.csv') # BRCAtumors_probesAge_merged.csv
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 756 214
df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x)) # Error in eval(predvars, data, env) : object 'cg16053920' not found
age_accel <- fitted.test-df_test.y
range(age_accel) # -119.0753  111.6396
hist(age_accel,breaks=50,main="BRCA allFeatures",xlim=c(-120,120))
mean(abs(age_accel)) # 22.75372
range(fitted.test) # -37.50926 174.13549
range(df_test.y) # 26 90
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 16.55597
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 18.8321
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 30.12031
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.3887047
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.1510914

cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(-40,175),
     ylim=c(-40,175)
)
abline(a=0,b=1)
text(x=-30,y=165,font=2,label="E)")
text(x=-11,y=153,font=1,label="RMSE: 30.12031")
mtext("All-features clock", side=1, line=4, cex=1, col="black") # run this on 7th plot; "line" changes how far from plot x-label is - use this on last two bottom plots

# abline(lm(fitted.test~df_test.y),col=c("blue")) # linear interpolation of fitted data
# BRCA tumors test -150K=============================================================
library(data.table)
set.seed(3)
df_test <- fread('BRCAtumors_probesAge_merged.csv')
df_test <- data.frame(df_test, row.names = 1)
df_test[1:5,1:5]
dim(df_test) # 756 253
df_test.x <- df_test[,-1] # not age
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
age_accel <- fitted.test-df_test.y
range(age_accel) # -87.37000  80.25763
hist(age_accel,breaks=50,main="BRCA -150K",xlim=c(-120,120))
sum(abs(age_accel)<7) # 322
mean(abs(age_accel)) # 16.55597
range(fitted.test) # -20.3700 150.6424
range(df_test.y) # 26 90
# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 22.75372
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test))
MAD # 
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2)) # greater difference between MAE and RMSE, greater variance in individual errors
RMSE # 22.10403
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr # 0.2905858
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.08444012

cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(-40,175),
     ylim=c(-40,175)
)
abline(a=0,b=1)
text(x=-30,y=165,font=2,label="F)")
text(x=-11,y=153,font=1,label="RMSE: 22.10403")
mtext("-150K clock", side=1, line=4, cex=1, col="black") # run this on 8th plot; "line" changes how far from plot x-label is - use this on last two bottom plots
mtext("DNAm age", side=2, line=0, cex=1.5, col="black",outer=T)
mtext("Chronological age", side=1, line=0, cex=1.5, col="black",outer=T)

### model performances on TCGA PRAD and LIHC====================================
# PRAD test allFeatures===============================================================
df_test <-fread('PRAD_probesAge_mergedCommon.csv')
dim(df_test) # 50 256803
df_test <- data.frame(df_test, row.names = 1) 
df_test.x <- df_test[,-1]
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(fitted.test) # 32.82242 90.79679
range(df_test.y) # 44 72

# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 9.146812 / 9.233554 (clock, clock_allFeatures)
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test)) # median absolute deviation (MAD) quantifies variation, how spread out set of data is 
MAD # 9.945353 / 14.27943 (clock, clock_allFeatures)
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 11.26036
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.4452714 / 0.3694121
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.1982666 / 0.1364653
# plot chronological age vs DNAm age for GSE37754 test
dev.new(width=20,height=20,noRStudioGD = TRUE)
par(mfrow=c(2,2), pin=c(4,2),oma=c(2,2,0,0))
cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(25,101),
     ylim=c(25,101)
)
abline(a=0,b=1)
text(x=27,y=95,font=2,label="A)")
text(x=36,y=89,font=1,label="RMSE: 11.26036")

# PRAD test -150K===============================================================
df_test <-fread('PRAD_probesAge_mergedCommon.csv')
dim(df_test) # 50 256803
df_test <- data.frame(df_test, row.names = 1) 
df_test.x <- df_test[,-1]
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 26.1682 100.0280
range(df_test.y) # 44 72

# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 9.146812 / 9.233554 (clock, clock_allFeatures)
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test)) # median absolute deviation (MAD) quantifies variation, how spread out set of data is 
MAD # 9.945353 / 14.27943 (clock, clock_allFeatures)
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 11.8528
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.4452714 / 0.3694121
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.1982666 / 0.1364653
# plot chronological age vs DNAm age for GSE37754 test
cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(25,101),
     ylim=c(25,101)
)
abline(a=0,b=1)
text(x=27,y=95,font=2,label="B)")
text(x=35,y=89,font=1,label="RMSE: 11.8528")
# LIHC test allFeatures===============================================================
df_test <-fread('LIHC_probesAge_mergedCommon.csv')
df_test <- data.frame(df_test, row.names = 1) 
dim(df_test) # 50
df_test.x <- df_test[,-1]
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock_allFeatures,df_test.x))
range(fitted.test) # 33.53429 107.50312
range(df_test.y) # 20 81

# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 8.506601 / 14.34994 (clock, clock_allFeatures)
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test)) # median absolute deviation (MAD) quantifies variation, how spread out set of data is 
MAD # 8.905077 / 12.84643 (clock, clock_allFeatures)
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 17.47927
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.8372636 / 0.745903
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.7010103 / 0.5563713
# plot chronological age vs DNAm age for GSE37754 test
cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(30,110),
     ylim=c(30,110)
)
abline(a=0,b=1)
text(x=32,y=105,font=2,label="C)")
text(x=41,y=99,font=1,label="RMSE: 17.47927")
mtext("All-features clock", side=1, line=4, cex=1, col="black") # run this on 8th plot; "line" changes how far from plot x-label is - use this on last two bottom plots

# LIHC test -150K===============================================================
df_test <-fread('LIHC_probesAge_mergedCommon.csv')
df_test <- data.frame(df_test, row.names = 1) 
dim(df_test) # 50
df_test.x <- df_test[,-1]
df_test.y <- df_test[,1]
fitted.test <- ageC.inverse(predict(clock,df_test.x))
range(fitted.test) # 32.26881 85.87026
range(df_test.y) # 20 81

# calculate mean absolute error
MAE <- sum(abs(df_test.y - fitted.test))/length(df_test.y)
MAE # 8.506601 / 14.34994 (clock, clock_allFeatures)
# calculate median absolute deviation
MAD <- mad((df_test.y-fitted.test)) # median absolute deviation (MAD) quantifies variation, how spread out set of data is 
MAD # 8.905077 / 12.84643 (clock, clock_allFeatures)
# RMSE
RMSE <- sqrt(mean((df_test.y-fitted.test)^2))
RMSE # 10.89658
# pearson correlation
pCorr <- cor(df_test.y,fitted.test)
pCorr #0.8372636 / 0.745903
# rsquared calculation
rsquared <- cor(df_test.y,fitted.test)^2
rsquared # 0.7010103 / 0.5563713
# plot chronological age vs DNAm age for GSE37754 test
cAge <- df_test.y
plot(cAge,fitted.test,
     main="",
     xlab="",
     ylab="",
     xlim=c(30,110),
     ylim=c(30,110)
)
abline(a=0,b=1)
text(x=32,y=105,font=2,label="D)")
text(x=41,y=99,font=1,label="RMSE: 10.89658")
mtext("-150K clock", side=1, line=4, cex=1, col="black") # run this on 8th plot; "line" changes how far from plot x-label is - use this on last two bottom plots
mtext("DNAm age", side=2, line=0, cex=1.5, col="black",outer=T)
mtext("Chronological age", side=1, line=0, cex=1.5, col="black",outer=T)