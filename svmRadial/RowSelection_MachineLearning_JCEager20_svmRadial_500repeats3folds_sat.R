####Author: Seth Tolley
####Date:   12/03/2021
####This script is to get the heritability of the RS data
####from various row segment plot means. 
####JCEager + JCEager

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
##################################################[JCEager20_40cmTRIM]###################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#####JCEager20 Yield 
#JCEager20_Yield <- read_excel("/scratch/bell/tolleys/UAV_DATA/MLPrediction//Yield/20JCEager_EndofSeason.xlsx", sheet=5)[-1,c(1,5,32)]
JCEager20_Yield <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Yield/20JCEager_EndofSeason.xlsx", sheet=5)[-1,c(1,5,32)]
JCEager20_Yield <- na.omit(JCEager20_Yield)
head(JCEager20_Yield)


#JCEager20 <- read_excel("/scratch/bell/tolleys/UAV_DATA/MLPrediction/RowSelection_AllData_sat.xlsx", sheet=4, na = "NA")
JCEager20 <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=4, na = "NA")
JCEager20 = JCEager20[,!grepl("intensity",names(JCEager20))]; JCEager20 = JCEager20[,!grepl("_ht",names(JCEager20))]; JCEager20 = JCEager20[,!grepl("_int",names(JCEager20))]

JCEager20_40cmTrim <- subset(JCEager20, JCEager20$Trim=="40cm")


####RS1234
JCEager1234 <- JCEager20_40cmTrim
RS1234 <- JCEager1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")




Data1234 <- data.frame()
# Set up Repeated k-fold Cross Validation
for(n in 1:1000){
  set.seed(24+n)
  
  data1234 <- merge(JCEager20_Yield, RS1234, by="Plot")
  rownames(data1234) <- data1234$Plot
  data1234 <- na.omit(data1234)
  data1234$Yield_kgha <- as.numeric(data1234$Yield_kgha)

  
  Temp <- subset(data1234, data1234$Type=="Temp")
  flds <- createFolds(Temp$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Temp <- cbind(flds, Temp)
  
  Trop <- subset(data1234, data1234$Type=="Trop")
  flds <- createFolds(Trop$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Trop <- cbind(flds, Trop)  
  
  #Phenotype File Prepared--------------------------------------------------------------   
  data1234 <- rbind(Temp, Trop)
  for(i in 1:3){
    training <- subset(data1234, data1234$flds!=i)
    testing <- subset(data1234, data1234$flds==i)
    
    rownames(training) <- training$Plot
    rownames(testing) <- testing$Plot
    
    training <- training[,-c(1:3)]
    testing <- testing[,-c(1:3)]
    
    # Fit the model 
    svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data1234)-1))
    JCEager20_svm_RS1234 <- train(Yield_kgha ~., 
                                  data = training, 
                                  method = "svmRadial",
                                  metric = "Rsquared",
                                  tuneGrid = svmGrid,
                                  preProcess = c("center","scale"))
    
    
    test1 <- predict(JCEager20_svm_RS1234, newdata=testing)
    
    d1234 <- cor(test1, testing$Yield_kgha)
    
    Data1234 <- rbind(Data1234,d1234)
    colnames(Data1234) <- "Correlation"
  }
  print(paste("Data1234:",n)) 
}  







####RS23
RS_data23 <- subset(JCEager20_40cmTrim, JCEager20_40cmTrim$RowSegment=="2"| JCEager20_40cmTrim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")






Data23 <- data.frame()
# Set up Repeated k-fold Cross Validation
for(n in 1:1000){
  set.seed(24+n)
  
  data23 <- merge(JCEager20_Yield, RS23, by="Plot")
  rownames(data23) <- data23$Plot
  data23 <- na.omit(data23)
  data23$Yield_kgha <- as.numeric(data23$Yield_kgha)

  
  Temp <- subset(data23, data23$Type=="Temp")
  flds <- createFolds(Temp$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Temp <- cbind(flds, Temp)
  
  Trop <- subset(data23, data23$Type=="Trop")
  flds <- createFolds(Trop$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Trop <- cbind(flds, Trop)  
  
  #Phenotype File Prepared--------------------------------------------------------------   
  data23 <- rbind(Temp, Trop)
  for(i in 1:3){
    training <- subset(data23, data23$flds!=i)
    testing <- subset(data23, data23$flds==i)
    
    rownames(training) <- training$Plot
    rownames(testing) <- testing$Plot
    
    training <- training[,-c(1:3)]
    testing <- testing[,-c(1:3)]
    
    # Fit the model 
    svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data23)-1))
    JCEager20_svm_RS23 <- train(Yield_kgha ~., 
                                  data = training, 
                                  method = "svmRadial",
                                  metric = "Rsquared",
                                  tuneGrid = svmGrid,
                                  preProcess = c("center","scale"))
    
    
    test1 <- predict(JCEager20_svm_RS23, newdata=testing)
    
    
    d23 <- cor(test1, testing$Yield_kgha)
    
    Data23 <- rbind(Data23,d23)
    colnames(Data23) <- "Correlation"
  }
  print(paste("Data23:",n))
}  






####RS14
RS_data14 <- subset(JCEager20_40cmTrim, JCEager20_40cmTrim$RowSegment=="1"| JCEager20_40cmTrim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")



data14 <- merge(JCEager20_Yield, RS14, by="Plot")
rownames(data14) <- data14$Plot
data14 <- na.omit(data14)
data14$Yield_kgha <- as.numeric(data14$Yield_kgha)
data14[1:5,1:5]


Data14 <- data.frame()

# Set up Repeated k-fold Cross Validation
for(n in 1:1000){
  set.seed(24+n)
  
  data14 <- merge(JCEager20_Yield, RS14, by="Plot")
  rownames(data14) <- data14$Plot
  data14 <- na.omit(data14)
  data14$Yield_kgha <- as.numeric(data14$Yield_kgha)

  
  Temp <- subset(data14, data14$Type=="Temp")
  flds <- createFolds(Temp$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Temp <- cbind(flds, Temp)
  
  Trop <- subset(data14, data14$Type=="Trop")
  flds <- createFolds(Trop$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Trop <- cbind(flds, Trop)  
  
  #Phenotype File Prepared--------------------------------------------------------------   
  data14 <- rbind(Temp, Trop)
  for(i in 1:3){
    training <- subset(data14, data14$flds!=i)
    testing <- subset(data14, data14$flds==i)
    
    rownames(training) <- training$Plot
    rownames(testing) <- testing$Plot
    
    training <- training[,-c(1:3)]
    testing <- testing[,-c(1:3)]
    
    # Fit the model 
    svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data14)-1))
    JCEager20_svm_RS14 <- train(Yield_kgha ~., 
                                  data = training, 
                                  method = "svmRadial",
                                  metric = "Rsquared",
                                  tuneGrid = svmGrid,
                                  preProcess = c("center","scale"))
    
    
    test1 <- predict(JCEager20_svm_RS14, newdata=testing)
    
    
    d14 <- cor(test1, testing$Yield_kgha)
    
    Data14 <- rbind(Data14,d14)
    colnames(Data14) <- "Correlation"
  }
  print(paste("Data14:",n))
}  




####RS1
RS_data1 <- subset(JCEager20_40cmTrim, JCEager20_40cmTrim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")





Data1 <- data.frame()
# Set up Repeated k-fold Cross Validation
for(n in 1:1000){
  set.seed(24+n)
  
  data1 <- merge(JCEager20_Yield, RS1, by="Plot")
  rownames(data1) <- data1$Plot
  data1 <- na.omit(data1)
  data1$Yield_kgha <- as.numeric(data1$Yield_kgha)
  
  Temp <- subset(data1, data1$Type=="Temp")
  flds <- createFolds(Temp$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Temp <- cbind(flds, Temp)
  
  Trop <- subset(data1, data1$Type=="Trop")
  flds <- createFolds(Trop$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Trop <- cbind(flds, Trop)  
  
  #Phenotype File Prepared--------------------------------------------------------------   
  data1 <- rbind(Temp, Trop)
  for(i in 1:3){
    training <- subset(data1, data1$flds!=i)
    testing <- subset(data1, data1$flds==i)
    
    rownames(training) <- training$Plot
    rownames(testing) <- testing$Plot
    
    training <- training[,-c(1:3)]
    testing <- testing[,-c(1:3)]
    
    # Fit the model 
    svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data1)-1))
    JCEager20_svm_RS1 <- train(Yield_kgha ~., 
                                  data = training, 
                                  method = "svmRadial",
                                  metric = "Rsquared",
                                  tuneGrid = svmGrid,
                                  preProcess = c("center","scale"))
    
    
    test1 <- predict(JCEager20_svm_RS1, newdata=testing)
    
    
    d1 <- cor(test1, testing$Yield_kgha)
    
    Data1 <- rbind(Data1,d1)
    colnames(Data1) <- "Correlation"
  }
  print(paste("Data1:",n))
}  









####RS2
RS_data2 <- subset(JCEager20_40cmTrim, JCEager20_40cmTrim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")





Data2 <- data.frame()

# Set up Repeated k-fold Cross Validation
for(n in 1:1000){
  set.seed(24+n)
  
  data2 <- merge(JCEager20_Yield, RS2, by="Plot")
  rownames(data2) <- data2$Plot
  data2 <- na.omit(data2)
  data2$Yield_kgha <- as.numeric(data2$Yield_kgha)

  
  Temp <- subset(data2, data2$Type=="Temp")
  flds <- createFolds(Temp$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Temp <- cbind(flds, Temp)
  
  Trop <- subset(data2, data2$Type=="Trop")
  flds <- createFolds(Trop$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Trop <- cbind(flds, Trop)  
  
  #Phenotype File Prepared--------------------------------------------------------------   
  data2 <- rbind(Temp, Trop)
  for(i in 1:3){
    training <- subset(data2, data2$flds!=i)
    testing <- subset(data2, data2$flds==i)
    
    rownames(training) <- training$Plot
    rownames(testing) <- testing$Plot
    
    training <- training[,-c(1:3)]
    testing <- testing[,-c(1:3)]
    
    # Fit the model 
    svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data2)-1))
    JCEager20_svm_RS2 <- train(Yield_kgha ~., 
                               data = training, 
                               method = "svmRadial",
                               metric = "Rsquared",
                               tuneGrid = svmGrid,
                               preProcess = c("center","scale"))
    
    
    test1 <- predict(JCEager20_svm_RS2, newdata=testing)
    
    
    d2 <- cor(test1, testing$Yield_kgha)
    
    Data2 <- rbind(Data2,d2)
    colnames(Data2) <- "Correlation"
  }
  print(paste("Data2:",n))
}  




####RS3
RS_data3 <- subset(JCEager20_40cmTrim, JCEager20_40cmTrim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


Data3 <- data.frame()

# Set up Repeated k-fold Cross Validation
for(n in 1:1000){
  set.seed(24+n)
  
  data3 <- merge(JCEager20_Yield, RS3, by="Plot")
  rownames(data3) <- data3$Plot
  data3 <- na.omit(data3)
  data3$Yield_kgha <- as.numeric(data3$Yield_kgha)
  
  Temp <- subset(data3, data3$Type=="Temp")
  flds <- createFolds(Temp$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Temp <- cbind(flds, Temp)
  
  Trop <- subset(data3, data3$Type=="Trop")
  flds <- createFolds(Trop$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Trop <- cbind(flds, Trop)  
  
  #Phenotype File Prepared--------------------------------------------------------------   
  data3 <- rbind(Temp, Trop)
  for(i in 1:3){
    training <- subset(data3, data3$flds!=i)
    testing <- subset(data3, data3$flds==i)
    
    rownames(training) <- training$Plot
    rownames(testing) <- testing$Plot
    
    training <- training[,-c(1:3)]
    testing <- testing[,-c(1:3)]
    
    # Fit the model 
    svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data3)-1))
    JCEager20_svm_RS3 <- train(Yield_kgha ~., 
                               data = training, 
                               method = "svmRadial",
                               metric = "Rsquared",
                               tuneGrid = svmGrid,
                               preProcess = c("center","scale"))
    
    
    test1 <- predict(JCEager20_svm_RS3, newdata=testing)
    
    
    d3 <- cor(test1, testing$Yield_kgha)
    
    Data3 <- rbind(Data3,d3)
    colnames(Data3) <- "Correlation"
  }
  print(paste("Data3:",n))
}  



####RS4
RS_data4 <- subset(JCEager20_40cmTrim, JCEager20_40cmTrim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


Data4 <- data.frame()

# Set up Repeated k-fold Cross Validation
for(n in 1:1000){
  set.seed(24+n)
  
  data4 <- merge(JCEager20_Yield, RS4, by="Plot")
  rownames(data4) <- data4$Plot
  data4 <- na.omit(data4)
  data4$Yield_kgha <- as.numeric(data4$Yield_kgha)
  
  Temp <- subset(data4, data4$Type=="Temp")
  flds <- createFolds(Temp$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Temp <- cbind(flds, Temp)
  
  Trop <- subset(data4, data4$Type=="Trop")
  flds <- createFolds(Trop$Yield_kgha, k=3, list=FALSE, returnTrain = FALSE)
  Trop <- cbind(flds, Trop)  
  
  #Phenotype File Prepared--------------------------------------------------------------   
  data4 <- rbind(Temp, Trop)
  for(i in 1:3){
    training <- subset(data4, data4$flds!=i)
    testing <- subset(data4, data4$flds==i)
    
    rownames(training) <- training$Plot
    rownames(testing) <- testing$Plot
    
    training <- training[,-c(1:3)]
    testing <- testing[,-c(1:3)]
    
    # Fit the model 
    svmGrid <- expand.grid(sigma= c(0.0001, 0.001, 0.01, 0.1), C=c(10, 50, 100, 150, 200, ncol(data4)-1))
    JCEager20_svm_RS4 <- train(Yield_kgha ~., 
                               data = training, 
                               method = "svmRadial",
                               metric = "Rsquared",
                               tuneGrid = svmGrid,
                               preProcess = c("center","scale"))
    
    
    test1 <- predict(JCEager20_svm_RS4, newdata=testing)
    
    
    d4 <- cor(test1, testing$Yield_kgha)
    
    Data4 <- rbind(Data4,d4)
    colnames(Data4) <- "Correlation"
  }
  print(paste("Data4:",n))
}  

#save.image(file='/scratch/bell/tolleys/UAV_DATA/MLPrediction/JCEager20_svmRadial_1000repeats3folds_workspace.RData')
save.image(file='~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/MLPrediction_Workspaces/JCEager20_svmRadial_1000repeats3folds_workspace.RData')



