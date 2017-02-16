#N.Hack
# Rockfish olive/yellowtail growth experiment
#Started 1 Feb 2017
#merging CCFRP data to my blood collection data




#Load my data
library(readr)
Fishing_data <- read_csv("C:/Users/manic/Dropbox/iRock/CCFRP blood samples/Fishing data.txt")
View(Fishing_data)
#load CCFRP data
library(readxl)
Hack_Sample_Data <- read_excel("C:/Users/manic/Dropbox/Home computer/CalPoly/Rockfsh/Hack_Sample_Data.xlsx")
View(Hack_Sample_Data)


#####Formatting my data#####
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
#date column registered as date by R
ahackdat$Date<-as.Date(ahackdat$Date,format="%m/%d/%Y")
#reorder to merge in correct order to match CCFRP data
ahackdat<-ahackdat[,c('Date','Location','Protection','sicell','Drift','Species','Length','TagNum','Bloodsamp')]



#####Formatting CCFRP data#####
#change species codes to match my data
##will need to add more once I get all the data!
fccfrp<-Hack_Sample_Data#fixed CCFRP data
library(car)
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
#fix to date function
#combine to get MM/DD/YYYY format
fccfrp$Date<-paste(paste(fccfrp$Month,fccfrp$Day,fccfrp$Year, sep='/'),sep='')
fccfrp$Date<-as.Date(fccfrp$Date,format="%B/%d/%Y")
#change tagID to not be a charCter with a fucking letter 'O' at the end
fccfrp$TagID<-as.numeric(substr(fccfrp$TagID,1,nchar(fccfrp$TagID)-1))
#reorder columns to combime in proper order
fccfrp<-fccfrp[,c('Date','Area','Site','Cell','Drift','Species','Length','TagID','StartLat','StartLon','EndLat','EndLon','Sex','Conditions')]


#####Subsetting data that matches#####
#subset those with tags
tccfrp<-subset(fccfrp,!is.na(TagID))
#subset those with tag#s
thackdat<-subset(ahackdat,!is.na(TagNum))
#merge matching tag numbers
tfish<-merge(tccfrp,thackdat[c('TagNum','Bloodsamp')],by.x='TagID',by.y='TagNum',all.x=F,all.y=F)

#subset matching length fish
#remove tagged fish identified in tfish above
notagccfrp<-fccfrp[is.na(fccfrp$TagID),]
#merge columns in CCFRP data in order to combine with my data
gccfrp<-as.data.frame(do.call(paste,as.data.frame((notagccfrp[,1:7]))))
colnames(gccfrp)<-'fishID'
#combine with unique data
gccfrp<-cbind(gccfrp,notagccfrp[8:14])
#remove tagged fish identified in tfish above in my data also
notaghackdat<-ahackdat[is.na(ahackdat$TagNum),]
#merging identification columns in my data
ghackdat<-do.call(paste,as.data.frame(notaghackdat[,1:7]))
#combine with unique data
ghackdat<-cbind(ghackdat,notaghackdat[9])
#merge matching fish!
mfish<-merge(gccfrp,ghackdat,by.x='fishID',by.y='ghackdat',all.x=F,all.y=F)

#final subset that does not match
#my data excluding tagged fish and matching fish (ghackdat already has tagged fish removed)
whatfish<-ghackdat[!(ghackdat$ghackdat %in% mfish$fishID),]
#count records with only unique fish
onefish<-as.data.frame(table(whatfish$ghackdat))
#find records that occur more than once
morefish<-onefish[onefish$Freq > 1,]




# how to deal with duplicates


