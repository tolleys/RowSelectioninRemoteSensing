####Author: Seth Tolley
####Date:   12/03/2021
####This script is to get the heritability of the RS data
####from various row segment plot means. 
####JCEager + TCZm

####Needed Libraries
library('readxl')
library("dplyr")
library("lme4")
library("agricolae")
library("dplyr")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####################################################[SbDIV18_noTRIM]###################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
SbDIV18 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=1, na = "NA")

SbDIV18_NOTrim <- subset(SbDIV18, SbDIV18$Trim=="NO")

SbDIV1234 <- SbDIV18_NOTrim
RS1234 <- SbDIV1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(SbDIV18_NOTrim, SbDIV18_NOTrim$RowSegment=="2"| SbDIV18_NOTrim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(SbDIV18_NOTrim, SbDIV18_NOTrim$RowSegment=="1"| SbDIV18_NOTrim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(SbDIV18_NOTrim, SbDIV18_NOTrim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(SbDIV18_NOTrim, SbDIV18_NOTrim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(SbDIV18_NOTrim, SbDIV18_NOTrim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(SbDIV18_NOTrim, SbDIV18_NOTrim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))


SbDIV18_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2018_SbDIV/2018_SbDIV_plotinfo.csv")[-1,-7]
SbDIV18_data <- merge(SbDIV18_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
SbDIV18_data[1:30,1:15]
colnum=c(7:ncol(SbDIV18_data))

for(i in 1:1134){
  
  
  df1 <- SbDIV18_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|Entry) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/2))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4

write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_SbDIV18_Heritability_noTrim.csv",
               row.names = FALSE)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#####################################################[SbDIV18_TRIM]####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
SbDIV18 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=1, na = "NA")

SbDIV18_Trim <- subset(SbDIV18, SbDIV18$Trim=="40cm")

SbDIV1234 <- SbDIV18_Trim
RS1234 <- SbDIV1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(SbDIV18_Trim, SbDIV18_Trim$RowSegment=="2"| SbDIV18_Trim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(SbDIV18_Trim, SbDIV18_Trim$RowSegment=="1"| SbDIV18_Trim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(SbDIV18_Trim, SbDIV18_Trim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(SbDIV18_Trim, SbDIV18_Trim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(SbDIV18_Trim, SbDIV18_Trim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(SbDIV18_Trim, SbDIV18_Trim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20180604:SR700670_20180802), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))


SbDIV18_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2018_SbDIV/2018_SbDIV_plotinfo.csv")[-1,-7]
SbDIV18_data <- merge(SbDIV18_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
SbDIV18_data[1:30,1:20]
colnum=c(7:ncol(SbDIV18_data))

for(i in 1:1134){
  
  
  df1 <- SbDIV18_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|Entry) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/2))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4


write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_SbDIV18_Heritability_40cmTrim.csv",
          row.names = FALSE)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####################################################[SbDIV19_noTRIM]###################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
SbDIV19 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=2, na = "NA")

SbDIV19_NOTrim <- subset(SbDIV19, SbDIV19$Trim=="NO")

SbDIV1234 <- SbDIV19_NOTrim
RS1234 <- SbDIV1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(SbDIV19_NOTrim, SbDIV19_NOTrim$RowSegment=="2"| SbDIV19_NOTrim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(SbDIV19_NOTrim, SbDIV19_NOTrim$RowSegment=="1"| SbDIV19_NOTrim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(SbDIV19_NOTrim, SbDIV19_NOTrim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(SbDIV19_NOTrim, SbDIV19_NOTrim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(SbDIV19_NOTrim, SbDIV19_NOTrim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(SbDIV19_NOTrim, SbDIV19_NOTrim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))

SbDIV19_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2019_SbDIV/WL19_SbdivTC.csv")[-1,-7]
head(SbDIV19_PlotInfo)
SbDIV19_data <- merge(SbDIV19_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
SbDIV19_data[1:10,1:10]
colnum=c(7:ncol(SbDIV19_data))

for(i in 1:1190){
  
  
  df1 <- SbDIV19_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|X) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/2))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4


write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_SbDIV19_Heritability_noTrim.csv",
          row.names = FALSE)






### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#####################################################[SbDIV19_TRIM]####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
SbDIV19 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=2, na = "NA")


SbDIV19_Trim <- subset(SbDIV19, SbDIV19$Trim=="40cm")

SbDIV1234 <- SbDIV19_Trim
RS1234 <- SbDIV1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(SbDIV19_Trim, SbDIV19_Trim$RowSegment=="2"| SbDIV19_Trim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(SbDIV19_Trim, SbDIV19_Trim$RowSegment=="1"| SbDIV19_Trim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(SbDIV19_Trim, SbDIV19_Trim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(SbDIV19_Trim, SbDIV19_Trim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(SbDIV19_Trim, SbDIV19_Trim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(SbDIV19_Trim, SbDIV19_Trim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20190702:SR700670_20190905), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))

SbDIV19_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2019_SbDIV/WL19_SbdivTC.csv")[-1,-7]
head(SbDIV19_PlotInfo)
SbDIV19_data <- merge(SbDIV19_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
SbDIV19_data[1:10,1:10]
colnum=c(7:ncol(SbDIV19_data))

for(i in 1:1190){
  
  
  df1 <- SbDIV19_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|X) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/2))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4


write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_SbDIV19_Heritability_40cmTrim.csv",
          row.names = FALSE)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####################################################[SbDIV20_noTRIM]###################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
SbDIV20 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=3, na = "NA")


SbDIV20_NOTrim <- subset(SbDIV20, SbDIV20$Trim=="NO")

SbDIV1234 <- SbDIV20_NOTrim
RS1234 <- SbDIV1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(SbDIV20_NOTrim, SbDIV20_NOTrim$RowSegment=="2"| SbDIV20_NOTrim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(SbDIV20_NOTrim, SbDIV20_NOTrim$RowSegment=="1"| SbDIV20_NOTrim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(SbDIV20_NOTrim, SbDIV20_NOTrim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(SbDIV20_NOTrim, SbDIV20_NOTrim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(SbDIV20_NOTrim, SbDIV20_NOTrim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(SbDIV20_NOTrim, SbDIV20_NOTrim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))

SbDIV20_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2020_SbDIV/2020_SbDIV_plotinfo.csv")[-1,-c(7,8)]
head(SbDIV20_PlotInfo)
SbDIV20_data <- merge(SbDIV20_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
SbDIV20_data[1:10,1:10]
colnum=c(7:ncol(SbDIV20_data))

for(i in 1:1274){
  
  
  df1 <- SbDIV20_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|X) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/2))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4
tail(Data4)

write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_SbDIV20_Heritability_noTrim.csv",
          row.names = FALSE)






### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#####################################################[SbDIV20_TRIM]####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
SbDIV20 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=3, na = "NA")
  
SbDIV20_Trim <- subset(SbDIV20, SbDIV20$Trim=="40cm")

SbDIV1234 <- SbDIV20_Trim
RS1234 <- SbDIV1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(SbDIV20_Trim, SbDIV20_Trim$RowSegment=="2"| SbDIV20_Trim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(SbDIV20_Trim, SbDIV20_Trim$RowSegment=="1"| SbDIV20_Trim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(SbDIV20_Trim, SbDIV20_Trim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(SbDIV20_Trim, SbDIV20_Trim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(SbDIV20_Trim, SbDIV20_Trim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(SbDIV20_Trim, SbDIV20_Trim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200625:SR700670_20200813), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))

SbDIV20_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2020_SbDIV/2020_SbDIV_plotinfo.csv")[-1,-c(7,8)]
head(SbDIV20_PlotInfo)
SbDIV20_data <- merge(SbDIV20_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
SbDIV20_data[1:10,1:10]
colnum=c(7:ncol(SbDIV20_data))

for(i in 1:1274){
  
  
  df1 <- SbDIV20_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|X) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/2))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4


write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_SbDIV20_Heritability_40cmTrim.csv",
          row.names = FALSE)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####################################################[20JCEager_noTRIM]#################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
JCEager20 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=4, na = "NA")

JCEager20_NOTrim <- subset(JCEager20, JCEager20$Trim=="NO")

JCEager20_1234 <- JCEager20_NOTrim
RS1234 <- JCEager20_1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(JCEager20_NOTrim, JCEager20_NOTrim$RowSegment=="2"| JCEager20_NOTrim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(JCEager20_NOTrim, JCEager20_NOTrim$RowSegment=="1"| JCEager20_NOTrim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(JCEager20_NOTrim, JCEager20_NOTrim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(JCEager20_NOTrim, JCEager20_NOTrim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(JCEager20_NOTrim, JCEager20_NOTrim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(JCEager20_NOTrim, JCEager20_NOTrim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))

JCEager20_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2020_JCEager/2020_JCEager_plotinfo.csv")[-1,-c(7:9)]
JCEager20_data <- merge(JCEager20_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
JCEager20_data[1:10,1:10]
colnum=c(7:ncol(JCEager20_data))

for(i in 1:1687){
  
  
  df1 <- JCEager20_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|Entry) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/3))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4[1684:1687, 1:4]

write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_JCEager20_Heritability_noTrim.csv",
          row.names = FALSE)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####################################################[20JCEager_TRIM]###################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
JCEager20 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=4, na = "NA")

JCEager20_Trim <- subset(JCEager20, JCEager20$Trim=="40cm")

JCEager20_1234 <- JCEager20_Trim
RS1234 <- JCEager20_1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(JCEager20_Trim, JCEager20_Trim$RowSegment=="2"| JCEager20_Trim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(JCEager20_Trim, JCEager20_Trim$RowSegment=="1"| JCEager20_Trim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(JCEager20_Trim, JCEager20_Trim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(JCEager20_Trim, JCEager20_Trim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(JCEager20_Trim, JCEager20_Trim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(JCEager20_Trim, JCEager20_Trim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20200617:SR700670_20200826), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))

JCEager20_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2020_JCEager/2020_JCEager_plotinfo.csv")[-1,-c(7:9)]
JCEager20_data <- merge(JCEager20_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
JCEager20_data[1:10,1:10]
colnum=c(7:ncol(JCEager20_data))

for(i in 1:1687){
  
  
  df1 <- JCEager20_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|Entry) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/3))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4


write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_JCEager20_Heritability_40cmTrim.csv",
          row.names = FALSE)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####################################################[21JCEager_noTRIM]#################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
JCEager21 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=5, na = "NA")

JCEager21_NOTrim <- subset(JCEager21, JCEager21$Trim=="NO")

JCEager21_1234 <- JCEager21_NOTrim
RS1234 <- JCEager21_1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(JCEager21_NOTrim, JCEager21_NOTrim$RowSegment=="2"| JCEager21_NOTrim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(JCEager21_NOTrim, JCEager21_NOTrim$RowSegment=="1"| JCEager21_NOTrim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(JCEager21_NOTrim, JCEager21_NOTrim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(JCEager21_NOTrim, JCEager21_NOTrim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(JCEager21_NOTrim, JCEager21_NOTrim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(JCEager21_NOTrim, JCEager21_NOTrim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))

JCEager21_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2021_JCEager/2021_JCEager_plotinfo.csv")[-1,-c(7:10)]
JCEager21_data <- merge(JCEager21_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
JCEager21_data[1:10,1:10]
colnum=c(7:ncol(JCEager21_data))

for(i in 1:1778){
  
  
  df1 <- JCEager21_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|Entry) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/3))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4

write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_JCEager21_Heritability_noTrim.csv",
          row.names = FALSE)






### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####################################################[21JCEager_TRIM]####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
JCEager21 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/RowSelection_AllData_sat.xlsx", sheet=5, na = "NA")

JCEager21_Trim <- subset(JCEager21, JCEager21$Trim=="40cm")

JCEager21_1234 <- JCEager21_Trim
RS1234 <- JCEager21_1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS1234)[-1] <- paste(colnames(RS1234)[-1], "RS1234", sep="_")



RS_data23 <- subset(JCEager21_Trim, JCEager21_Trim$RowSegment=="2"| JCEager21_Trim$RowSegment=="3")
RS23 <- RS_data23 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS23)[-1] <- paste(colnames(RS23)[-1], "RS23", sep="_")


RS_data14 <- subset(JCEager21_Trim, JCEager21_Trim$RowSegment=="1"| JCEager21_Trim$RowSegment=="4")
RS14 <- RS_data14 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS14)[-1] <- paste(colnames(RS14)[-1], "RS14", sep="_")


RS_data1 <- subset(JCEager21_Trim, JCEager21_Trim$RowSegment=="1")
RS1 <- RS_data1 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS1)[-1] <- paste(colnames(RS1)[-1], "RS1", sep="_")


RS_data2 <- subset(JCEager21_Trim, JCEager21_Trim$RowSegment=="2")
RS2 <- RS_data2 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS2)[-1] <- paste(colnames(RS2)[-1], "RS2", sep="_")


RS_data3 <- subset(JCEager21_Trim, JCEager21_Trim$RowSegment=="3")
RS3 <- RS_data3 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS3)[-1] <- paste(colnames(RS3)[-1], "RS3", sep="_")


RS_data4 <- subset(JCEager21_Trim, JCEager21_Trim$RowSegment=="4")
RS4 <- RS_data4 %>%
  group_by(Plot) %>%
  summarise_at(vars(CC_20210610:SR700670_20210910), mean, na.rm = TRUE)
colnames(RS4)[-1] <- paste(colnames(RS4)[-1], "RS4", sep="_")


d1 <- merge(RS1234, RS23, by="Plot")
d2 <- merge(d1, RS14, by="Plot")
d3 <- merge(d2, RS1, by="Plot")
d4 <- merge(d3, RS2, by="Plot")
d5 <- merge(d4, RS3, by="Plot")
d6 <- merge(d5, RS4, by="Plot")

d7 <- dplyr::select(d6, -contains("int"))

JCEager21_PlotInfo <- read.csv("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Plot_Extraction/Fields/2021_JCEager/2021_JCEager_plotinfo.csv")[-1,-c(7:10)]
JCEager21_data <- merge(JCEager21_PlotInfo, d7, by='Plot')


drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
JCEager21_data[1:10,1:10]
colnum=c(7:ncol(JCEager21_data))

for(i in 1:1778){
  
  
  df1 <- JCEager21_data
  
  x=colnum[i]  #set the current [i] column as "x"
  trait=colnames(df1)[x] #sets the current column header as the trait name
  colnames(df1)[x]="y"  #renames the trait variable as "y" for the model analysis below
  
  
  #We are interested in random effects, which estimates the proportion of variation and not fixed effects. 
  #Knowing variability components allows us to calculate Heritability.
  #Random-effects terms are distinguished by vertical bars or pipes (|) 
  model <- lmer(y ~ (1|Entry) +(Rep), data=df1) #PedInfo, Year, P*Y and REP1 (nested in year) were all significant
  summary(model)
  varComp<-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
  
  
  #deleting columns in variance component dataframe (we don't need it)
  varComp<-varComp[ , !(names(varComp) %in% drops)]
  sigmaE <- varComp[2,2]
  sigmaG <- varComp[1,2]
  herit <- sigmaG/(sigmaG + (sigmaE/3))
  #set the trait name from the loop
  varComp$Trait<-trait
  #add columns to existing dataframe
  Data3 <- cbind(trait, sigmaG, sigmaE, herit)
  Data4 <- rbind(Data4, Data3)
}

Data4


write.csv(Data4, "~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Heritability/RowSelection_JCEager21_Heritability_40cmTrim.csv",
          row.names = FALSE)




