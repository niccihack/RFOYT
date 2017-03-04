#N.Hack
# Rockfish olive/yellowtail growth experiment
#Started 28 Feb 2017
#statistics



#Load data
library(readr)
RFOYTgrowth <- read_csv("~/RFOYT/RFOYT/RFOYTgrowth.csv")
#make new dataframe to play with
growth<-as.data.frame(RFOYTgrowth)


#Calculate condition factor


#Calculate somatic growth rates
#reorganize data into long form
library(tidyr)
b = tidyr::gather(growth, 'data.type','value', 5:16)
b$mass.sl = gsub(".*\\-",'',b$data.type)#make separate data type column - mass or SL
b$day = as.integer(substr( gsub("\\-.*","",b$data.type), 4,5))#make separate day column
# generate all unique combinations of days 
c = combn(unique(b$day), 2)
dat = data.frame()#make empty dataframe for loop

for (d in (1: ncol(c))){
  e = c[,d]
  #pull out the mass for those two days: 
  f = dplyr::filter(b, day == max(e) & mass.sl == 'mass')$value
  g = dplyr::filter(b, day == min(e) & mass.sl == 'mass')$value
  sgr = log ( (f-g)/(max(e)-min(e))*100 )#somatic growth rate formula
  treat = dplyr::filter(b, day == max(e) & mass.sl == 'mass')$treatment#pull out treatment
  IGF1 = dplyr::filter(b, day == max(e) & mass.sl == 'mass')$IGF1#and IGF1
  PITtag = dplyr::filter(b, day == max(e) & mass.sl == 'mass')$PITtag#and tag number
  Comb = paste('Day',e[1], e[2], sep = '-')#make title
  dat = rbind(dat, data.frame(PITtag, treat, IGF1, sgr, Comb))#combine all into dataframe
}


 

#####Pairwise ANOVAs#####


#####Pearson's correlations#####



####MANIVA#####