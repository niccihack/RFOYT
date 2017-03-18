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

#calculate Condition factor
dat.cf<-data.frame()
m = dplyr::filter(b, mass.sl == 'mass')$value
l = dplyr::filter(b, mass.sl == 'SL')$value
cf = (m/(l^3))*100 #condition factor formula
dat.cf = rbind(dat.cf, data.frame(b$PITtag, b$day, b$treatment, cf))#combine all into dataframe


# generate all unique combinations of days 
comb = combn(unique(b$day), 2)

#Calculate SGR for mass
#subset out mass
fish.mass <- dplyr::filter(b, mass.sl == 'mass')
tag <- dplyr::filter(fish.mass, day == 0) %>% dplyr::select(.,PITtag)
treatment <- dplyr::filter(fish.mass, day == 0) %>% dplyr::select(.,treatment)
IGF1 <- dplyr::filter(fish.mass, day == 0) %>% dplyr::select(.,IGF1)
dat.mass<-data.frame()
for (d in 1:ncol(comb)){
  fday = comb[2,d]
  iday = comb[1,d]
  f = dplyr::filter(fish.mass, day == fday)$value
  i = dplyr::filter(fish.mass, day == iday)$value
  sgr.mass = log (((f-i)/(fday - iday))*100 )#somatic growth rate formula
  Comb = paste('Day',iday, fday, sep = '-')#make title
  #rbind 
  dat.mass = rbind(dat.mass, data.frame(tag, treatment, IGF1, Comb, sgr.mass))#combine all into dataframe  
}


#calculate SGR for SL
fish.sl <- dplyr::filter(b, mass.sl == 'SL')
dat.sl<-data.frame()
for (d in 1:ncol(comb)){
  fday = comb[2,d]
  iday = comb[1,d]
  f = dplyr::filter(fish.sl, day == fday)$value
  i = dplyr::filter(fish.sl, day == iday)$value
  sgr.sl = log (((f-i)/(fday - iday))*100 )#somatic growth rate formula
  Comb = paste('Day',iday, fday, sep = '-')#make title
  #rbind 
  dat.sl = rbind(dat.sl, data.frame(tag, treatment, IGF1, Comb, sgr.sl))#combine all into dataframe  
}



#calculate condition factor
dat.sl<-data.frame()
for (d in 1:ncol(comb)){
  fday = comb[2,d]
  iday = comb[1,d]
  f = dplyr::filter(fish.sl, day == fday)$value
  i = dplyr::filter(fish.sl, day == iday)$value
  sgr.sl = log (((f-i)/(fday - iday))*100 )#somatic growth rate formula
  Comb = paste('Day',iday, fday, sep = '-')#make title
  #rbind 
  dat.sl = rbind(dat.sl, data.frame(tag, treatment, IGF1, Comb, sgr.sl))#combine all into dataframe  
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
plot(aov.mass[])
#then run anova
myanova<-function(i){
  anova.mass<-anova((lm(value~treatment,data=mass,day==i,na.action=na.omit)))
}
anova.stats<- sapply(unique(mass$day), myanova, simplify=FALSE)
names(anova.stats)<-paste('Day -',unique(mass$day))#rename outputs to include days
capture.output(anova.stats,file='anova.stats.mass.csv')#export file

#SL
sl<-subset(b, mass.sl=='SL')#subset SL
library(ggplot2)
sl.plot=ggplot(sl, aes(treatment,value,colour=treatment))+ geom_boxplot()#look at data
#check analysis of variance
myaov<-function(h){
  aov.sl<-aov(lm(value~treatment,data=sl,day==h,na.action=na.omit))
}
#apply to each day
aov.stats <- sapply(unique(sl$day), myaov, simplify=FALSE)
#check normality and residuals if necessary
#plot(aov.sl[])
#then run anova
myanova<-function(i){
  anova.sl<-anova((lm(value~treatment,data=sl,day==i,na.action=na.omit)))
}
anova.stats<- sapply(unique(sl$day), myanova, simplify=FALSE)
names(anova.stats)<-paste('Day -',unique(sl$day))#rename outputs to include days
capture.output(anova.stats,file='anova.stats.sl.csv')#export file

#CF
cf.plot=ggplot(dat.cf, aes(b.treatment,cf,colour=b.treatment))+ geom_boxplot()#look at data  
#check analysis of variance
myaov<-function(h){
  aov.cf<-aov(lm(cf~b.treatment,data=dat.cf,b.day==h,na.action=na.omit))
}
#apply to each day
aov.stats <- sapply(unique(dat.cf$b.day), myaov, simplify=FALSE)
#check normality and residuals if necessary
#plot(aov.sl[])
#then run anova
myanova<-function(i){
  anova.cf<-anova((lm(cf~b.treatment,data=dat.cf,b.day==i,na.action=na.omit)))
}
anova.stats<- sapply(unique(dat.cf$b.day), myanova, simplify=FALSE)
names(anova.stats)<-paste('Day -',unique(dat.cf$b.day))#rename outputs to include days
capture.output(anova.stats,file='anova.stats.cf.csv')#export file

#####Pearson's correlations#####
#must remove Inf and NAs since pearson's cannot deal with them
cor.sgr<-dplyr::filter(dat.sl, sgr.SL!=-Inf & !is.na(b.IGF1))


for (k in unique(cor.sgr$Comb)){
  tempx<-cor.sgr[cor.sgr$Comb==k,]
  t
}

cor.sgr %>%
  dplyr::select(b.IGF1,sgr.SL) %>%
  dplyr::filter()

cor(cor.sgr[c(3,4)])#need to do for each time interval
cor.test(~ b.IGF1 + sgr.SL, data=cor.sgr[c(3,4)])

cor.sgr %>% 
  dplyr::group_by(Comb, b.treatment) %>% 
  
mycor<- function(j){
  cor.test(~ b.IGF1 + sgr.SL, data=j)
}





####MANOVA#####