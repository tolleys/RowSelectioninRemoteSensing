####Author: Seth Tolley
####Date:   12/03/2021
####This script is to get the heritability of the Yieldspectral data
####from various row segment plot means. 
####JCEager + TCZm

####Needed Libraries
library('readxl')
library("dplyr")
library("lme4")
library("agricolae")
library("ggplot2")
library("RColorBrewer")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#####################################################[JCEager Loop]####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
JCEager_Yield21 <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA_Processing/Yield/21JCEager_YieldData_sat.xlsx", sheet=1)
head(JCEager_Yield21)
JCEager_Yield21$Yield_kgha <- as.numeric(JCEager_Yield21$Yield_kgha)
Yield1234 <- JCEager_Yield21
H1234 <- Yield1234 %>%
  group_by(Plot) %>%
  summarise_at(vars(Yield_kgha), mean, na.rm = TRUE)
colnames(H1234)[-1] <- paste(colnames(H1234)[-1], "RS1234", sep="_")
head(H1234)


Yield23 <- subset(JCEager_Yield21, JCEager_Yield21$RowSegment=="2"| JCEager_Yield21$RowSegment=="3")
H23 <- Yield23 %>%
  group_by(Plot) %>%
  summarise_at(vars(Yield_kgha), mean, na.rm = TRUE)
colnames(H23)[-1] <- paste(colnames(H23)[-1], "RS23", sep="_")


Yield14 <- subset(JCEager_Yield21, JCEager_Yield21$RowSegment=="1"| JCEager_Yield21$RowSegment=="4")
H14 <- Yield14 %>%
  group_by(Plot) %>%
  summarise_at(vars(Yield_kgha), mean, na.rm = TRUE)
colnames(H14)[-1] <- paste(colnames(H14)[-1], "RS14", sep="_")


Yield1 <- subset(JCEager_Yield21, JCEager_Yield21$RowSegment=="1")
H1 <- Yield1 %>%
  group_by(Plot) %>%
  summarise_at(vars(Yield_kgha), mean, na.rm = TRUE)
colnames(H1)[-1] <- paste(colnames(H1)[-1], "RS1", sep="_")


Yield2 <- subset(JCEager_Yield21, JCEager_Yield21$RowSegment=="2")
H2 <- Yield2 %>%
  group_by(Plot) %>%
  summarise_at(vars(Yield_kgha), mean, na.rm = TRUE)
colnames(H2)[-1] <- paste(colnames(H2)[-1], "RS2", sep="_")


Yield3 <- subset(JCEager_Yield21, JCEager_Yield21$RowSegment=="3")
H3 <- Yield3 %>%
  group_by(Plot) %>%
  summarise_at(vars(Yield_kgha), mean, na.rm = TRUE)
colnames(H3)[-1] <- paste(colnames(H3)[-1], "RS3", sep="_")


Yield4 <- subset(JCEager_Yield21, JCEager_Yield21$RowSegment=="4")
H4 <- Yield4 %>%
  group_by(Plot) %>%
  summarise_at(vars(Yield_kgha), mean, na.rm = TRUE)
colnames(H4)[-1] <- paste(colnames(H4)[-1], "RS4", sep="_")


d1 <- merge(H1234, H23, by="Plot")
d2 <- merge(d1, H14, by="Plot")
d3 <- merge(d2, H1, by="Plot")
d4 <- merge(d3, H2, by="Plot")
d5 <- merge(d4, H3, by="Plot")
d6 <- merge(d5, H4, by="Plot")



JCEAGER_PlotInfo <- read_excel("~/Documents/Purdue/PhD/Research/Chp3_RowSelection/DATA/2021_DATA/Reference_Data/WL21_JCEager.TCZm_Reference.xlsx", sheet=1)[-1,c(1:6)]
JCEager_data <- merge(JCEAGER_PlotInfo, d6, by='Plot')
head(JCEager_data)

drops <- c("var1","var2","sdcor") 

Data3 <- data.frame()
Data4 <- data.frame()
JCEager_data[1:10,1:10]
colnum=c(7:ncol(JCEager_data))

for(i in 1:7){
  
  
  df1 <- JCEager_data
  
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
Data4$Num <- c(1:7)
head(Data4)

ggplot(Data4, aes(x=Num, y=as.numeric(herit))) +
  geom_bar(stat="identity", colour="black", fill="steelblue4") +
  theme_bw() + 
  #scale_fill_manual(values=c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"))+
  #scale_fill_manual(values=c("#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF")) +
  #scale_fill_brewer(palette = "BuGn") +
  scale_x_discrete(name ="Row Segment", 
                   limits=c("RS1234", "RS23", "RS14", "RS1", "RS2", "RS3", "RS4")) +
  theme(plot.title=element_text(hjust=.5, face="bold", size=18)) +
  theme(axis.title = element_text(face="bold", size=14))+
  theme(axis.text = element_text(size=14)) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(face="bold", size=14)) +
  theme(legend.position = "none") +
  labs(x="Row Segment",
       y=expression(bold(Heritability~"("*sigma[g]^{2}*"/"*sigma[g]^{2}*"+"*sigma[epsilon]^{2}*")")),
       fill="Row Segments")
  
  



