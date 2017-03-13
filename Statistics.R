#N.Hack
# Rockfish olive/yellowtail growth experiment
#Started 28 Feb 2017
#statistics



#Load data
library(readr)
RFOYTgrowth <- read_csv("~/RFOYT/RFOYT/RFOYTgrowth.csv")
#make new dataframe to play with
growth<-as.data.frame(RFOYTgrowth)


#reorganize data into long form
library(tidyr)
b = tidyr::gather(growth, 'data.type','value', 5:16)
b$mass.sl = gsub(".*\\-",'',b$data.type)#make separate data type column - mass or SL
b$day = as.integer(substr( gsub("\\-.*","",b$data.type), 4,5))#make separate day column

#calculate condition factor
dat.cf<-data.frame()
for (d in (1:length(b))){
  m = dplyr::filter(b, mass.sl == 'mass')$value
  l = dplyr::filter(b, mass.sl == 'SL')$value
  cf = (m/(l^3))*100 #condition factor formula
  dat.cf = rbind(dat.cf, data.frame(b$PITtag, b$day, cf))#combine all into dataframe
}
# generate all unique combinations of days 
c = combn(unique(b$day), 2)


#Calculate somatic growth rates for mass
dat.mass = data.frame()#make empty dataframe for loop
for (d in (1: ncol(c))){
  e = c[,d]
  #pull out the mass for those two days: 
  f = dplyr::filter(b, day == max(e) & mass.sl == 'mass')$value
  g = dplyr::filter(b, day == min(e) & mass.sl == 'mass')$value
  sgr.mass = log ( (f-g)/(max(e)-min(e))*100 )#somatic growth rate formula
  Comb = paste('Day',e[1], e[2], sep = '-')#make title
  dat.mass = rbind(dat.mass, data.frame(b$PITtag, b$treatment, b$IGF1, sgr.mass, Comb))#combine all into dataframe
}


#calculate SGR for SL
dat.sl = data.frame()#make empty dataframe for loop
for (d in (1: ncol(c))){
  e = c[,d]
  #pull out the sl for those two days: 
  f = dplyr::filter(b, day == max(e) & mass.sl == 'SL')$value
  g = dplyr::filter(b, day == min(e) & mass.sl == 'SL')$value
  sgr.SL = log ( (f-g)/(max(e)-min(e))*100 )#somatic growth rate formula
  Comb = paste('Day',e[1], e[2], sep = '-')#make title
  dat.sl = rbind(dat.sl, data.frame(b$PITtag, b$treatment, b$IGF1, sgr.SL, Comb))#combine all into dataframe
}

 

#####Pairwise ANOVAs#####

##mass
mass<-subset(b, mass.sl=='mass')#subset mass
library(ggplot2)
mass.plot=ggplot(mass, aes(treatment,value,colour=treatment))+ geom_boxplot()#look at data
#check analysis of variance
myaov<-function(h){
  aov.mass<-aov(lm(value~treatment,data=mass,day==h,na.action=na.omit))
}
#apply to each day
aov.stats <- sapply(unique(mass$day), myaov, simplify=FALSE)
#check normality and residuals if necessary
#plot(aov.mass[])
#then run anova
myanova<-function(i){
  anova.mass<-anova((lm(value~treatment,data=mass,day==i,na.action=na.omit)))
}
anova.stats<- sapply(unique(mass$day), myanova, simplify=FALSE)
names(anova.stats)<-paste('Day -',unique(mass$day))#rename outputs to include days
capture.output(anova.stats,file='anova.stats.csv')#export file

#SL



#####Pearson's correlations#####



####MANIVA#####