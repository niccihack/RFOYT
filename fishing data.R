#N.Hack
# Rockfish olive/yellowtail growth experiment
#Started 1 Feb 2017
#merging CCFRP data to my blood collection data




#Load my data
library(readr)
Fishing_data <- read_csv("C:/Users/manic/Dropbox/iRock/CCFRP blood samples/Fishing data.txt")
View(Fishing_data)
#adjust my data to match CCFRP
#because cell is coded as site.cell without a deliminiter, I have to recreate this in my data
ahackdat<-Fishing_data#adjusted Fishing data
ahackdat["SiteCell"]<-ifelse(ahackdat$Location=='PBL',paste('BL'),paste('PB'))
#adding zeros to single cell numbers
ahackdat['ncell']<-ifelse(nchar(ahackdat$Cell)==1,
                       c(paste0('0',ahackdat$Cell)),
                       ahackdat$Cell
                       )

ahackdat['sicell']<-paste(ahackdat$SiteCell,ahackdat$ncell, sep='')
#reorder to merge in correct order to match CCFRP data
ahackdat<-ahackdat[,c('Date','Location','Protection','sicell','Drift','Species','Length','TagNum','Bloodsamp')]
#merging identification columns in my data
ghackdat<-do.call(paste,as.data.frame(ahackdat[,1:7]))
#combine with unique data
ghackdat<-cbind(ghackdat,ahackdat[8:9])
#sort to line up with CCFRP data
ghackdat$ghackdat<-sort(ghackdat$ghackdat,decreasing = F)


#load CCFRP data
library(readxl)
Hack_Sample_Data <- read_excel("C:/Users/manic/Dropbox/Home computer/CalPoly/Rockfsh/Hack_Sample_Data.xlsx")
View(Hack_Sample_Data)
#change species codes to match my data
##will need to add more once I get all the data!
fccfrp<-Hack_Sample_Data#fixed CCFRP data
fccfrp$Species<-recode(fccfrp$Species,"
  'BLA'='RFBLK';
  'BLU'='RFBLU';
  'CPR'='RFCOP';
  'GPR'='RFGOP';
  'KLP'='RFKLP';
  'LCD'='LNGCD';
  'OLV'='RFOLV';
  'VER'='RFVER'
")
#change area code to match my location code
fccfrp$Area<-recode(fccfrp$Area,"
    'BL'='PBL';
    'PB'='PBN'
")
#change month to number
fccfrp$Month<-recode(fccfrp$Month,"
    'July'='07';
    'August'='08';
    'September'='09'
")
#combine to get MM/DD/YYYY format
fccfrp['Date']<-paste('0',paste(fccfrp$Month,fccfrp$Day,fccfrp$Year, sep='/'),sep='')
#reorder columns to combime in proper order
fccfrp<-fccfrp[,c('Date','Site','Area','Cell','Drift','Species','Length','TagID','StartLat','StartLon','EndLat','EndLon','Sex','Conditions')]
#merge columns inorder to combine with my data
gccfrp<-do.call(paste,as.data.frame(fccfrp[,1:7]))
#combine with unique data
gccfrp<-cbind(gccfrp,fccfrp[8:14])
#sort to line up with my data
gccfrp$gccfrp<-sort(gccfrp$gccfrp,decreasing = F)


#Combine data sets!
cfishdata<-merge(gccfrp,ghackdat,by='')#must be same length
