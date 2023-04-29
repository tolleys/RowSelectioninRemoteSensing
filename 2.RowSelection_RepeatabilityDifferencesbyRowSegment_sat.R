####Author: Seth Tolley
####Date:   12/03/2021
####This script is to get the LSD of the RS data
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
#################################################[Sorghum_TrimEffect]##################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
RS_Sorghum <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/RowSelection_Heritability_sat.xlsx", sheet=1, na="NA")

RS_Sorghum <- subset(RS_Sorghum, RS_Sorghum$herit>0)
head(RS_Sorghum)

Data3 <- as.data.frame(matrix(nrow=1, ncol=3))
Data4 <- data.frame(Trim = c("ANOVA", "40cm", "NO"))


for(j in 1:3){
  for(i in 1:30){
    data1 <- subset(RS_Sorghum, RS_Sorghum$GDDNum==j)
    data1.1 <- subset(data1, data1$TraitNum==i)
    Trait = paste(data1.1$Trait[[1]], data1.1$GDDNum[[1]], sep="_")
    groups = paste("groups",data1.1$Trait[[1]], data1.1$GDDNum[[1]], sep="_")
    
    
    model1 <- lm(herit ~ Trim*RowSegment*Experiment, data=data1.1)
    m1 <- anova(model1)
    print(Trait)
    print(m1)
    m2 <- LSD.test(model1, trt=c("Trim"), DFerror=anova(model1)[8,1], MSerror = anova(model1[8,3]), alpha=0.05, group=TRUE)
    
    Data3[1,1] <- m1[1,5]
    Data3[1,2] <- if(Data3[1,1]<0.001){
      "***"
    } else if (Data3[1,1]<0.01){
      "**"
    } else if (Data3[1,1]<0.05){
      "*"
    }else {
      "NS"
    }
    Data3[1,3] <- "ANOVA"
    colnames(Data3) <- c(Trait, groups, "Trim")  
    
    data2 <- m2$groups
    data2$RowSegment <- rownames(data2)
    colnames(data2) <- c(Trait, groups, "Trim")
    data2 <- rbind(Data3, data2)
    
    Data4 <- merge(Data4, data2, by="Trim")
  }
}
Data4


write.csv(Data4, "~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Heritability/RowSelection_HeritGreaterthan0Dif_Sorghum_Trimming_sat.csv")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###################################################[Maize_TrimEffect]##################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
RS_Maize <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/RowSelection_Heritability_sat.xlsx", sheet=2, na="NA")

RS_Maize <- subset(RS_Maize, RS_Maize$herit>0)


Data3 <- as.data.frame(matrix(nrow=1, ncol=3))
Data4 <- data.frame(Trim = c("ANOVA", "40cm", "NO"))

for(j in 1:3){
  for(i in 1:30){
    data1 <- subset(RS_Maize, RS_Maize$GDDNum==j)
    data1.1 <- subset(data1, data1$TraitNum==i)
    Trait = paste(data1.1$Trait[[1]], data1.1$GDDNum[[1]], sep="_")
    groups = paste("groups",data1.1$Trait[[1]], data1.1$GDDNum[[1]], sep="_")
    
    
    model1 <- lm(herit ~ Trim*RowSegment*Experiment, data=data1.1)
    m1 <- anova(model1)
    print(Trait)
    print(m1)
    m2 <- LSD.test(model1, trt=c("Trim"), DFerror=anova(model1)[8,1], MSerror = anova(model1[8,3]), alpha=0.05, group=TRUE)
    
    Data3[1,1] <- m1[1,5]
    Data3[1,2] <- if(Data3[1,1]<0.001){
      "***"
    } else if (Data3[1,1]<0.01){
      "**"
    } else if (Data3[1,1]<0.05){
      "*"
    }else {
      "NS"
    }
    Data3[1,3] <- "ANOVA"
    colnames(Data3) <- c(Trait, groups, "Trim")  
    
    data2 <- m2$groups
    data2$RowSegment <- rownames(data2)
    colnames(data2) <- c(Trait, groups, "Trim")
    data2 <- rbind(Data3, data2)
    
    Data4 <- merge(Data4, data2, by="Trim")
  }
}
Data4


write.csv(Data4, "~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Heritability/RowSelection_HeritGreaterthan0Dif_Maize_Trimming_sat.csv")













### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#######################################################[Sorghum]#######################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
RS_Sorghum <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/RowSelection_Heritability_sat.xlsx", sheet=1, na="NA")

#RS_Sorghum <- subset(RS_Sorghum, RS_Sorghum$Trim=="40cm")
RS_Sorghum <- subset(RS_Sorghum, RS_Sorghum$herit>0.0)


m1 <- lm(herit ~ RowSegment*Experiment, data=RS_Sorghum)
anova(m1)
#No RowSegment*Experiment Interaction so we combined across the years


Data3 <- as.data.frame(matrix(nrow=1, ncol=3))
Data4 <- data.frame(RowSegment = c("ANOVA", "RS23", "RS1234", "RS14", "RS1", "RS2", "RS3", "RS4"))


for(j in 1:3){
  for(i in 1:30){
  data1 <- subset(RS_Sorghum, RS_Sorghum$GDDNum==j)
  data1.1 <- subset(data1, data1$TraitNum==i)
    Trait = paste(data1.1$Trait[[1]], data1.1$GDDNum[[1]], sep="_")
    groups = paste("groups",data1.1$Trait[[1]], data1.1$GDDNum[[1]], sep="_")
  

    model1 <- lm(herit ~ Trim*RowSegment*Experiment, data=data1.1)
    m1 <- anova(model1)
    print(Trait)
    print(m1)
    m2 <- LSD.test(model1, trt=c("RowSegment"), DFerror=anova(model1)[8,1], MSerror = anova(model1[8,3]), alpha=0.05, group=TRUE)
    
    Data3[1,1] <- m1[2,5]
    Data3[1,2] <- if(Data3[1,1]<0.001){
      "***"
    } else if (Data3[1,1]<0.01){
      "**"
    } else if (Data3[1,1]<0.05){
      "*"
    }else {
      "NS"
    }
    Data3[1,3] <- "ANOVA"
    colnames(Data3) <- c(Trait, groups, "RowSegment")  
    
    data2 <- m2$groups
    data2$RowSegment <- rownames(data2)
    colnames(data2) <- c(Trait, groups, "RowSegment")
    data2 <- rbind(Data3, data2)
    
    Data4 <- merge(Data4, data2, by="RowSegment")
  }
}
Data4


write.csv(Data4, "~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Heritability/RowSelection_HeritGreaterthan0Dif_Sorghum_RowSegment_sat.csv")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#######################################################[Maize]#########################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
RS_Maize <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/RowSelection_Heritability_sat.xlsx", sheet=2, na="NA")

RS_Maize <- subset(RS_Maize, RS_Maize$herit>0.0)






Data3 <- as.data.frame(matrix(nrow=1, ncol=3))
Data4 <- data.frame(RowSegment = c("ANOVA", "RS23", "RS1234", "RS14", "RS1", "RS2", "RS3", "RS4"))


for(j in 1:3){
  for(i in 1:30){
    data1 <- subset(RS_Maize, RS_Maize$GDDNum==j)
    data1.1 <- subset(data1, data1$TraitNum==i)
    Trait = paste(data1.1$Trait[[1]], data1.1$GDDNum[[1]], sep="_")
    groups = paste("groups",data1.1$Trait[[1]], data1.1$GDDNum[[1]], sep="_")
    
    
    model1 <- lm(herit ~ Trim*RowSegment*Experiment, data=data1.1)
    m1 <- anova(model1)
    print(Trait)
    print(m1)
    m2 <- LSD.test(model1, trt=c("RowSegment"), DFerror=anova(model1)[8,1], MSerror = anova(model1[8,3]), alpha=0.05, group=TRUE)
    
    Data3[1,1] <- m1[2,5]
    Data3[1,2] <- if(Data3[1,1]<0.001){
      "***"
    } else if (Data3[1,1]<0.01){
      "**"
    } else if (Data3[1,1]<0.05){
      "*"
    }else {
      "NS"
    }
    Data3[1,3] <- "ANOVA"
    colnames(Data3) <- c(Trait, groups, "RowSegment")  
    
    data2 <- m2$groups
    data2$RowSegment <- rownames(data2)
    colnames(data2) <- c(Trait, groups, "RowSegment")
    data2 <- rbind(Data3, data2)
    
    Data4 <- merge(Data4, data2, by="RowSegment")
  }
}
Data4

write.csv(Data4, "~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Heritability/RowSelection_HeritGreaterthan0Dif_Maize_RowSegment_sat.csv")
