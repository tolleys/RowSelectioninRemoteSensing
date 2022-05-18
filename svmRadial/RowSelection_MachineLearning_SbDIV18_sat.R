####Author: Seth Tolley
####Date:   12/03/2021
####This script is to get the heritability of the RS data
####from various row segment plot means. 
####JCEager + TCZm

####Needed Libraries
install.packages('readxl',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/", dependencies=TRUE)
library('readxl')
install.packages('dplyr',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/", dependencies=TRUE)
library("dplyr")
install.packages('caret',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/", dependencies=TRUE)
library("caret")
install.packages('kernlab',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/", dependencies=TRUE)
library("kernlab")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##################################################[SbDIV18_40cmTRIM]###################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#####SbDIV18 Yield 
SbDIV18_Yield <- read.csv("/scratch/bell/tolleys/UAV_DATA/MLPrediction/Yield/sbdiv_biomass_final_updated_Ali.csv", header=TRUE, )[,c(1,12)]
SbDIV18_Yield <- na.omit(SbDIV18_Yield)
head(SbDIV18_Yield)



SbDIV18 <- read_excel("/scratch/bell/tolleys/UAV_DATA/MLPrediction/RowSelection_AllData_sat.xlsx", sheet=1, na = "NA")
SbDIV18 = SbDIV18[,!grepl("intensity",names(SbDIV18))]; SbDIV18 = SbDIV18[,!grepl("_ht",names(SbDIV18))]; SbDIV18 = SbDIV18[,!grepl("_int",names(SbDIV18))]

SbDIV18_40cmTrim <- subset(SbDIV18, SbDIV18$Trim=="40cm")


####RS1234
SbDIV1234 <- SbDIV18_40cmTrim
RS1234 <- SbDIV1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



data1234 <- merge(SbDIV18_Yield, RS1234, by="Plot")
rownames(data1234) <- data1234$Plot
data1234 <- data1234[,-1]
data1234 <- na.omit(data1234)

# Set up Repeated k-fold Cross Validation
set.seed(24)
train_control <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions = TRUE)
svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data1234)-1))

# Fit the model 
SbDIV18_svm_RS1234 <- train(Plot_dry_wt_gramspermeter2 ~., 
                            data = data1234, 
                            method = "svmRadial", 
                            metric = "Rsquared",
                            trControl = train_control,  
                            tuneGrid = svmGrid,
                            preProcess = c("center","scale"))
#View the model
SbDIV18_svm_RS1234







####RS23
RS_data23 <- subset(SbDIV18_40cmTrim, SbDIV18_40cmTrim$RowSegment=="2"| SbDIV18_40cmTrim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")

data23 <- merge(SbDIV18_Yield, RS23, by="Plot")
rownames(data23) <- data23$Plot
data23 <- data23[,-1]
data23 <- na.omit(data23)

# Set up Repeated k-fold Cross Validation
set.seed(24)
train_control <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions = TRUE)
svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data23)-1))

# Fit the model 
SbDIV18_svm_RS23 <- train(Plot_dry_wt_gramspermeter2 ~., 
                            data = data23, 
                            method = "svmRadial", 
                            metric = "Rsquared",
                            trControl = train_control,
                            tuneGrid = svmGrid,
                            preProcess = c("center","scale"))
                           
#View the model
SbDIV18_svm_RS23






####RS14
RS_data14 <- subset(SbDIV18_40cmTrim, SbDIV18_40cmTrim$RowSegment=="1"| SbDIV18_40cmTrim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")




data14 <- merge(SbDIV18_Yield, RS14, by="Plot")
rownames(data14) <- data14$Plot
data14 <- data14[,-1]
data14 <- na.omit(data14)

# Set up Repeated k-fold Cross Validation
set.seed(24)
train_control <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions = TRUE)
svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data14)-1))


# Fit the model 
SbDIV18_svm_RS14 <- train(Plot_dry_wt_gramspermeter2 ~., 
                          data = data14, 
                          method = "svmRadial", 
                          metric = "Rsquared",
                          trControl = train_control, 
                          tuneGrid = svmGrid,
                          preProcess = c("center","scale"))
                          
#View the model
SbDIV18_svm_RS14




####RS1
RS_data1 <- subset(SbDIV18_40cmTrim, SbDIV18_40cmTrim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


data1 <- merge(SbDIV18_Yield, RS1, by="Plot")
rownames(data1) <- data1$Plot
data1 <- data1[,-1]
data1 <- na.omit(data1)

# Set up Repeated k-fold Cross Validation
set.seed(24)
train_control <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions = TRUE)
svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data1)-1))


# Fit the model 
SbDIV18_svm_RS1 <- train(Plot_dry_wt_gramspermeter2 ~., 
                          data = data1, 
                          method = "svmRadial", 
                          metric = "Rsquared",
                          trControl = train_control, 
                          tuneGrid = svmGrid,
                          preProcess = c("center","scale"))
                          
#View the model
SbDIV18_svm_RS1






####RS2
RS_data2 <- subset(SbDIV18_40cmTrim, SbDIV18_40cmTrim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


data2 <- merge(SbDIV18_Yield, RS2, by="Plot")
rownames(data2) <- data2$Plot
data2 <- data2[,-1]
data2 <- na.omit(data2)

# Set up Repeated k-fold Cross Validation
set.seed(24)
train_control <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions = TRUE)
svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data2)-1))


# Fit the model 
SbDIV18_svm_RS2 <- train(Plot_dry_wt_gramspermeter2 ~., 
                         data = data2, 
                         method = "svmRadial", 
                         metric = "Rsquared",
                         trControl = train_control, 
                         tuneGrid = svmGrid,
                         preProcess = c("center","scale"))
                         
#View the model
SbDIV18_svm_RS2




####RS3
RS_data3 <- subset(SbDIV18_40cmTrim, SbDIV18_40cmTrim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


data3 <- merge(SbDIV18_Yield, RS3, by="Plot")
rownames(data3) <- data3$Plot
data3 <- data3[,-1]
data3 <- na.omit(data3)

# Set up Repeated k-fold Cross Validation
set.seed(24)
train_control <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions = TRUE)
svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data3)-1))


# Fit the model 
SbDIV18_svm_RS3 <- train(Plot_dry_wt_gramspermeter2 ~., 
                         data = data3, 
                         method = "svmRadial", 
                         metric = "Rsquared",
                         trControl = train_control,  
                         tuneGrid = svmGrid,
                         preProcess = c("center","scale"))
#View the model
SbDIV18_svm_RS3



####RS4
RS_data4 <- subset(SbDIV18_40cmTrim, SbDIV18_40cmTrim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")

data4 <- merge(SbDIV18_Yield, RS4, by="Plot")
rownames(data4) <- data4$Plot
data4 <- data4[,-1]
data4 <- na.omit(data4)

# Set up Repeated k-fold Cross Validation
set.seed(24)
train_control <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions = TRUE)
svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data4)-1))


# Fit the model 
SbDIV18_svm_RS4 <- train(Plot_dry_wt_gramspermeter2 ~., 
                         data = data4, 
                         method = "svmRadial",
                         metric = "Rsquared",
                         trControl = train_control,  
                         tuneGrid = svmGrid,
                         preProcess = c("center","scale"))
#View the model
SbDIV18_svm_RS4



save.image(file='/scratch/bell/tolleys/UAV_DATA/MLPrediction/SbDIV18_svmRadial_workspace.RData')
