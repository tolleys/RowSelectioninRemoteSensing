###Author = Seth Tolley
###Data = 10/18/2021
###Description = Envi takes a long time to extract the hyperspectral data for each plot.
###              Instead I thought I would try to write a script that will loop through my
###              ROI and through each of the wavelengths to get the information for each 
###              row segment in JCEAGER and TCZm in 2021. 



###########[Needed Libraries]#########################
library("rgdal")
library("raster")
library("reshape2")
library("ggplot2")
library("dplyr")
library('readxl')

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##################################################[JCEager 2020 Loop]##################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
dirname <- "~/Documents/Purdue/PhD/UAV_Data/DATA/JCEager/2020/VNIR"
filename <- file.path(dirname, "20200702_mosaic_f54N_JCEager.tif")
file.exists(filename)


roi <- readOGR("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Plot_Extraction/Fields/2020_JCEager/ShpFile/2020_f54N_JCEager.shp")



VNIR_Mosaic <- raster(filename, band = 100)
plot(VNIR_Mosaic)
lines(roi)


DataOutput_JCEager_noTrim <- as.data.frame(matrix(data=NA,nrow=240,ncol=139))
DataOutput_JCEager_noTrim[,1] <- roi$plot
DataOutput_JCEager_noTrim[,2] <- roi$row


DataOutput_JCEager_Trim <- as.data.frame(matrix(data=NA,nrow=240,ncol=139))
DataOutput_JCEager_Trim[,1] <- roi$plot
DataOutput_JCEager_Trim[,2] <- roi$row


ColNames <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/2021_JCEager_Hyper_ColumnNames_sat.xlsx", sheet=1)
colnames(DataOutput_JCEager_noTrim) <- colnames(ColNames)
colnames(DataOutput_JCEager_Trim) <- colnames(ColNames)


for(i in 1:137){
  
  VNIR_Mosaic <- raster(filename, band = i)
  
  for(j in 4001:4060){
    roi_j <- subset(roi, roi$plot == j)
    
    for(n in 1:4){
      #No Trim results
      roi_jn_noTrim <- subset(roi_j, roi_j$row==n)
      
      JCEager_ClippedPlots_noTrim <- crop(VNIR_Mosaic, roi_jn_noTrim)
      DataOutput_JCEager_noTrim[((j-4001)*4+n),i+2] <- mean(JCEager_ClippedPlots_noTrim@data@values, na.rm=TRUE)/10000
      
      #40cm Trim results
      roi_jn_Trim <- subset(roi_j, roi_j$row==n)
      
      if ((roi_jn_Trim@bbox[2,2]  - 0.4)- (roi_jn_Trim@bbox[2,1] +0.4) > 0){
        roi_jn_Trim@bbox[2,1] <- roi_jn_Trim@bbox[2,1] + 0.4 #Trimming 40cm from the plots
        roi_jn_Trim@bbox[2,2] <- roi_jn_Trim@bbox[2,2] - 0.4 #Trimming 40cm from the plots
      } else {
        roi_jn_Trim@bbox[2,1] <- roi_jn_Trim@bbox[2,1] #Some of the plots are not big enough to be trimmed
        roi_jn_Trim@bbox[2,2] <- roi_jn_Trim@bbox[2,2] #Used shrinking function and it never found plants.
      }
      
      JCEager_ClippedPlots_Trim <- crop(VNIR_Mosaic, roi_jn_Trim)
      DataOutput_JCEager_Trim[((j-4001)*4+n),i+2] <- mean(JCEager_ClippedPlots_noTrim@data@values, na.rm=TRUE)/10000
      
      
    }
    
  }
  print(i)
}


write.csv(DataOutput_JCEager_noTrim, '~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/Processed/JCEager/WL21_20210616JCEager_noTrim_Reflectance.csv',
          row.names = FALSE)

write.csv(DataOutput_JCEager_Trim, '~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/Processed/JCEager/WL21_20210616JCEager_40cmTrim_Reflectance.csv',
          row.names = FALSE)





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##################################################[JCEager 2021 Loop]##################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
dirname <- "~/Documents/Purdue/PhD/UAV_Data/DATA/JCEager/2021/VNIR/"
filename <- file.path(dirname, "f42_20210616_VNIR_4cm.tif")
file.exists(filename)


roi <- readOGR("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Plot_Extraction/Fields/2021_JCEager/ShpFile/2021_f42m_JCEager.shp")



VNIR_Mosaic <- raster(filename, band = 80)
plot(VNIR_Mosaic, ylim=c(4480650,4480775))
lines(roi)


DataOutput_JCEager_noTrim <- as.data.frame(matrix(data=NA,nrow=240,ncol=139))
DataOutput_JCEager_noTrim[,1] <- roi$plot
DataOutput_JCEager_noTrim[,2] <- roi$row


DataOutput_JCEager_Trim <- as.data.frame(matrix(data=NA,nrow=240,ncol=139))
DataOutput_JCEager_Trim[,1] <- roi$plot
DataOutput_JCEager_Trim[,2] <- roi$row

  
ColNames <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/2021_JCEager_Hyper_ColumnNames_sat.xlsx", sheet=1)
colnames(DataOutput_JCEager_noTrim) <- colnames(ColNames)
colnames(DataOutput_JCEager_Trim) <- colnames(ColNames)


for(i in 1:137){
    
    VNIR_Mosaic <- raster(filename, band = i)

    for(j in 4001:4060){
      roi_j <- subset(roi, roi$plot == j)
      
      for(n in 1:4){
        #No Trim results
        roi_jn_noTrim <- subset(roi_j, roi_j$row==n)
        
          JCEager_ClippedPlots_noTrim <- crop(VNIR_Mosaic, roi_jn_noTrim)
          DataOutput_JCEager_noTrim[((j-4001)*4+n),i+2] <- mean(JCEager_ClippedPlots_noTrim@data@values, na.rm=TRUE)/10000
        
          #40cm Trim results
        roi_jn_Trim <- subset(roi_j, roi_j$row==n)
        
        if ((roi_jn_Trim@bbox[2,2]  - 0.4)- (roi_jn_Trim@bbox[2,1] +0.4) > 0){
          roi_jn_Trim@bbox[2,1] <- roi_jn_Trim@bbox[2,1] + 0.4 #Trimming 40cm from the plots
          roi_jn_Trim@bbox[2,2] <- roi_jn_Trim@bbox[2,2] - 0.4 #Trimming 40cm from the plots
        } else {
          roi_jn_Trim@bbox[2,1] <- roi_jn_Trim@bbox[2,1] #Some of the plots are not big enough to be trimmed
          roi_jn_Trim@bbox[2,2] <- roi_jn_Trim@bbox[2,2] #Used shrinking function and it never found plants.
        }
      
        JCEager_ClippedPlots_Trim <- crop(VNIR_Mosaic, roi_jn_Trim)
        DataOutput_JCEager_Trim[((j-4001)*4+n),i+2] <- mean(JCEager_ClippedPlots_noTrim@data@values, na.rm=TRUE)/10000
        
      
    }

  }
      print(i)
}


write.csv(DataOutput_JCEager_noTrim, '~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/Processed/JCEager/WL21_20210616JCEager_noTrim_Reflectance.csv',
          row.names = FALSE)

write.csv(DataOutput_JCEager_Trim, '~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/Processed/JCEager/WL21_20210616JCEager_40cmTrim_Reflectance.csv',
          row.names = FALSE)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##################################################[SbDIV 2018 Loop]####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
dirname <- "~/Documents/Purdue/PhD/UAV_Data/DATA/SbDiv/2018/VNIR/"
filename <- file.path(dirname, "20180802_mosaic_f54s.tif")
file.exists(filename)




roi <- readOGR("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Plot_Extraction/Fields/2018_SbDIV/ShpFile/2018_f54s_SbDIV_NAD83.shp")

roi
VNIR_Mosaic <- raster(filename, band = 80)
plot(VNIR_Mosaic, ylim=c(4480400,4480775))
lines(roi)


DataOutput_SbDIV18 <- as.data.frame(matrix(data=NA,nrow=5040,ncol=138))
DataOutput_SbDIV18[,1] <- roi$plot
DataOutput_SbDIV18[,2] <- roi$row


ColNames <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/2018_SbDIV_Hyper_ColumnNames_sat.xlsx", sheet=1)
colnames(DataOutput_SbDIV18) <- colnames(ColNames)


for(i in 1:136){
  
  VNIR_Mosaic <- raster(filename, band = i)
  
  for(j in 4601:5860){
    roi_j <- subset(roi, roi$plot == j)
    
    for(n in 1:4){
      roi_jn <- subset(roi_j, roi_j$row==n)
      
      if ((roi_jn@bbox[2,2]  - 0.4)- (roi_jn@bbox[2,1] +0.4) > 0){
        roi_jn@bbox[2,1] <- roi_jn@bbox[2,1] + 0.4 #Trimming 40cm from the plots
        roi_jn@bbox[2,2] <- roi_jn@bbox[2,2] - 0.4 #Trimming 40cm from the plots
      } else {
        roi_jn@bbox[2,1] <- roi_jn@bbox[2,1] #Some of the plots are not big enough to be trimmed
        roi_jn@bbox[2,2] <- roi_jn@bbox[2,2] #Used shrinking function and it never found plants.
      }
      
      SbDIV18_ClippedPlots <- crop(VNIR_Mosaic, roi_jn)
      
      DataOutput_SbDIV18[((j-4601)*4+n),i+2] <- mean(SbDIV18_ClippedPlots@data@values, na.rm=TRUE)/10000
      
    }
    
  }
  print(i)
}




 write.csv(DataOutput_SbDIV18, '~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/Processed/WL18_20180802_SbDIV_Reflectance.NEWTrim.csv',
          row.names = FALSE)

 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##################################################[SbDIV 2019 Loop]####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
dirname <- "~/Documents/Purdue/PhD/UAV_Data/DATA/SbDiv/2019/VNIR/"
filename <- file.path(dirname, "20190723_mosaic_f42m.tif")
file.exists(filename)




roi <- readOGR("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Plot_Extraction/Fields/2019_SbDIV/ShpFile/2019_f42_SbDIV.shp")
               
roi

VNIR_Mosaic <- raster(filename, band = 80)
plot(VNIR_Mosaic)
lines(roi)


DataOutput_SbDIV19 <- as.data.frame(matrix(data=NA,nrow=5040,ncol=138))
DataOutput_SbDIV19[,1] <- roi$plot
DataOutput_SbDIV19[,2] <- roi$row

#Not all of the flights have the same wavelengths in 2019.
ColNames <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/2019_SbDIV_Hyper_ColumnNames_sat.xlsx", sheet=1)
colnames(DataOutput_SbDIV19) <- colnames(ColNames)


for(i in 1:136){
  
  VNIR_Mosaic <- raster(filename, band = i)
  
  for(j in 5001:6260){
    roi_j <- subset(roi, roi$plot == j)
    
    for(n in 1:4){
      roi_jn <- subset(roi_j, roi_j$row==n)
      
      #if ((roi_jn@bbox[2,2]  - 0.4)- (roi_jn@bbox[2,1] +0.4) > 0){
      #  roi_jn@bbox[2,1] <- roi_jn@bbox[2,1] + 0.4 #Trimming 40cm from the plots
      #  roi_jn@bbox[2,2] <- roi_jn@bbox[2,2] - 0.4 #Trimming 40cm from the plots
      #} else {
      #  roi_jn@bbox[2,1] <- roi_jn@bbox[2,1] #Some of the plots are not big enough to be trimmed
      #  roi_jn@bbox[2,2] <- roi_jn@bbox[2,2] #Used shrinking function and it never found plants.
      #}
      
      
      SbDIV19_ClippedPlots <- crop(VNIR_Mosaic, roi_jn)
      
      DataOutput_SbDIV19[((j-5001)*4+n),i+2] <- mean(SbDIV19_ClippedPlots@data@values, na.rm=TRUE)/10000
      #DataOutput_SbDIV19[((j-5001)*4+n),i+2] <- mean(SbDIV19_ClippedPlots@data@values, na.rm=TRUE)
      
    }
    
  }
  print(i)
}



write.csv(DataOutput_SbDIV19, '~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/Processed/WL19_20190723_SbDIV_Reflectance_NEW.csv',
          row.names = FALSE)





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##################################################[SbDIV 2020 Loop]####################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
dirname <- "~/Documents/Purdue/PhD/UAV_Data/DATA/SbDiv/2020/VNIR/"
filename <- file.path(dirname, "20200619_mosaic_f54s.tif")
file.exists(filename)




roi <- readOGR("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Plot_Extraction/Fields/2020_SbDIV/ShpFile/2020_f54s_SbDIV_NAD83.shp")

roi
VNIR_Mosaic <- raster(filename, band = 80)
plot(VNIR_Mosaic)
lines(roi)


DataOutput_SbDIV20 <- as.data.frame(matrix(data=NA,nrow=5220,ncol=138))
DataOutput_SbDIV20[,1] <- roi$plot
DataOutput_SbDIV20[,2] <- roi$row


ColNames <- read_excel("~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/2020_SbDIV_Hyper_ColumnNames_sat.xlsx", sheet=1)
colnames(DataOutput_SbDIV20) <- colnames(ColNames)


for(i in 1:136){
  
  VNIR_Mosaic <- raster(filename, band = i)
  
  for(j in 3401:4705){
    roi_j <- subset(roi, roi$plot == j)
    
    for(n in 1:4){
      roi_jn <- subset(roi_j, roi_j$row==n)
      
      #roi_jn@bbox[2,1] <- roi_jn@bbox[2,1] + 0.4 #Trimming 40cm from the plots
      #roi_jn@bbox[2,2] <- roi_jn@bbox[2,2] - 0.4 #Trimming 40cm from the plots
      
      SbDIV20_ClippedPlots <- crop(VNIR_Mosaic, roi_jn)
      
      DataOutput_SbDIV20[((j-3401)*4+n),i+2] <- mean(SbDIV20_ClippedPlots@data@values, na.rm=TRUE)/10000
      
    }
    
  }
  print(i)
}




write.csv(DataOutput_SbDIV20, '~/Documents/Purdue/PhD/UAV_Data/DATA_Processing/Hyperspectral/Processed/WL20_20200619_SbDIV_Reflectance.csv',
          row.names = FALSE)


