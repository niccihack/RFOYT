#N.Hack
# Rockfish olive/yellowtail growth experiment
#Started 28 Feb 2017
#statistics



#Load data
library(readr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(car)
library(ggplot2)
library(lmerTest)

RFOYTgrowth <- read.csv("~/RFOYT/RFOYT/RFOYTgrowth.csv")
#make new dataframe to play with
growth<-as.data.frame(RFOYTgrowth)
#remove rows b[b$PITtag != 63]


#reorganize data into long form
b = tidyr::gather(growth, 'data.type','value', 5:16)
b$mass.sl = gsub(".*\\.",'',b$data.type)#make separate data type column - mass or SL
b$day = as.integer(substr( gsub("\\-.*","",b$data.type), 4,5))#make separate day column


# generate all unique combinations of days 
comb = combn(unique(b$day), 2)

#Calculate SGR for mass
#subset out mass
fish.mass <- dplyr::filter(b, mass.sl == 'mass')
#pull out tag, treatment, and igf1 values for one day (no repeats)
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

#calculate Condition factor
dat.cf<-data.frame()
day <- fish.mass$day
tank<- fish.mass$tank
m = dplyr::filter(b, mass.sl == 'mass')$value
l = dplyr::filter(b, mass.sl == 'SL')$value
cf = (m/(l^3))*100 #condition factor formula
dat.cf = rbind(dat.cf, data.frame(tag, day, treatment, tank, IGF1, cf))#combine all into dataframe

#####Pairwise ANOVAs#####

##mass
mass<-subset(b, mass.sl=='mass')#subset mass
mass$treatment <- as.factor(mass$treatment)
mass.plot=ggplot(mass, aes(treatment,value,colour=treatment))+ geom_boxplot()+ylab('mass')+ theme_classic()#remove background
ggsave(filename = 'mass.plot.pdf')
#then run anova
myanova<-function(i){
  anova.mass<-anova((lm(value~treatment,data=mass,day==i,na.action=na.omit)))
}
anova.stats<- sapply(unique(mass$day), myanova, simplify=FALSE)
names(anova.stats)<-paste('Day -',unique(mass$day))#rename outputs to include days
capture.output(anova.stats,file='anova.stats.mass.txt')#export file

#SL
sl<-subset(b, mass.sl=='SL')#subset SL
sl.plot=ggplot(sl, aes(treatment,value,colour=treatment))+ geom_boxplot()+ylab('SL') + theme_classic()
ggsave(filename = 'sl.plot.pdf')
#then run anova
myanova<-function(i){
  anova.sl<-anova((lm(value~treatment,data=sl,day==i,na.action=na.omit)))
}
anova.stats<- sapply(unique(sl$day), myanova, simplify=FALSE)
names(anova.stats)<-paste('Day -',unique(sl$day))#rename outputs to include days
capture.output(anova.stats,file='anova.stats.sl.csv')#export file

#CF
cf.plot=ggplot(dat.cf, aes(treatment,cf,colour=treatment))+ geom_boxplot() + theme_classic()
ggsave(filename = 'cf.plot.pdf')
#then run anova
myanova<-function(i){
  anova.cf<-anova((lm(cf~treatment,data=dat.cf,day==i,na.action=na.omit)))
}
anova.stats<- sapply(unique(dat.cf$day), myanova, simplify=FALSE)
names(anova.stats)<-paste('Day -',unique(dat.cf$day))#rename outputs to include days
capture.output(anova.stats,file='anova.stats.cf.csv')#export file

#####Pearson's correlations#####

#SL
#must remove Inf and NAs since pearson's cannot deal with them
cor.sgr<-dplyr::filter(dat.sl, sgr.sl!=-Inf & !is.na(IGF1))
cor.sl <- list()
for (k in unique(cor.sgr$Comb)){
  sl = dplyr::filter(cor.sgr, Comb==k)$sgr.sl#extract sgr of sl
  igf1 = dplyr::filter(cor.sgr, Comb==k)$IGF1#extract IGF concentrations
  correlation = cor.test(sl,igf1,method='pearson')#calculate pearson's correlations
  correlation$data.name = as.character(k)#rename to include day interval
  cor.sl = cbind(cor.sl,correlation)#combine into one database
}
capture.output(cor.sl,file='cor.sl.csv')#export file


#Mass
#must remove Inf and NAs since pearson's cannot deal with them
cor.sgr<-dplyr::filter(dat.mass, sgr.mass!=-Inf & !is.na(IGF1))
cor.mass <- list()
for (k in unique(cor.sgr$Comb)){
  mass = dplyr::filter(cor.sgr, Comb==k)$sgr.mass#extract sgr of mass
  igf1 = dplyr::filter(cor.sgr, Comb==k)$IGF1#extract IGF concentrations
  correlation = cor.test(mass,igf1,method='pearson')#calculate pearson's correlations
  correlation$data.name = as.character(k)#rename to include day interval
  cor.mass = cbind(cor.mass,correlation)#combine into one database
}
capture.output(cor.mass,file='cor.mass.csv')#export file

#CF
cor.cf <- list()
for (l in unique(dat.cf$day)){#only have to run look per day
  igf1 = dplyr::filter(dat.cf, day == l)$IGF1
  cf = dplyr::filter(dat.cf, day == l)$cf
  correlation = cor.test(igf1,cf,method='pearson')
  correlation$data.name = as.character(l)
  cor.cf = cbind(cor.cf,correlation)
}
capture.output(cor.cf, file='cor.cf.csv')


####Repeated measures ANOVA#####

#mass
b.mass = dplyr::filter(b, mass.sl == 'mass')#subset mass
ggplot(b.mass, aes(x = day, y =value, group = PITtag, color = treatment))+
  geom_line()+ylab('mass (g)')+
  scale_x_continuous(expand = c(0,0),limits = c(0,98),breaks = c(0,24,48,75,91,98))+ 
  theme_classic()#plot masses over time for each fish
mean.mass = group_by(b.mass,day,treatment) %>% summarise(mean_mass=mean(value),sem_mass = sd(value)/sqrt(length(value)),n=length(value))
mass.rmANOVA.plot = ggplot(mean.mass, aes(x = day, y =mean_mass, color = treatment))+
  geom_point()+ylab('mass (g)')+
  geom_line()+
  geom_errorbar(aes(ymin=mean_mass - sem_mass, ymax = mean_mass + sem_mass), width = 1.5) +
  scale_x_continuous(expand = c(0,0),limits = c(-10,108),breaks = c(0,24,48,75,91,98))+ 
  theme_classic()#plot masses over time for each fish
ggsave(filename = 'mass.rmANOVA.pdf')
b.mass$day = as.factor(b.mass$day)#change days to factors
b.mass$tank =as.factor(b.mass$tank)#change tank number to factor
b.mass$PITtag = as.factor(b.mass$PITtag)#change tag number to factors
mass.lmer2 =  lmer(value~day*treatment + (1|tank/PITtag), data = b.mass, REML = F) #full model
mass.lmer = lmer(value~day*treatment+(1|PITtag), data = b.mass) #test tank effect
anova(mass.lmer2,mass.lmer)#tank not significant so excluded
mass.lm = lm(value~day*treatment, data = b.mass)#simplest model
anova(mass.lmer,mass.lm)#PITtag significant so included
summary(mass.lmer)
anova.mass.lmer = anova(mass.lmer)
capture.output(anova.mass.lmer,file='anova.mass.lmer.csv')#export file


#SL
b.length = dplyr::filter(b, mass.sl == 'SL')#subset length
ggplot(b.length, aes(x = day, y =value, group = PITtag, color = treatment))+
  geom_line()+ylab('length (cm)')+
  scale_x_continuous(expand = c(0,0),limits = c(0,98),breaks = c(0,24,48,75,91,98))+ 
  theme_classic()#plot masses over time for each fish
mean.length = group_by(b.length,day,treatment) %>% summarise(mean_length=mean(value),sem_length = sd(value)/sqrt(length(value)),n=length(value))
sl.rmANOVA.plot = ggplot(mean.length, aes(x = day, y = mean_length, color = treatment))+
  geom_line()+ylab('length (cm)')+
  geom_point()+
  geom_errorbar(aes(ymin=mean_length - sem_length, ymax = mean_length + sem_length), width = 1.5) +
  scale_x_continuous(expand = c(0,0),limits = c(-10,108),breaks = c(0,24,48,75,91,98))+ 
  theme_classic()#plot masses over time for each fish
ggsave(filename = 'sl.rmANOVA.pdf')
b.length$day = as.factor(b.length$day)#change days to factors
b.length$tank =as.factor(b.length$tank)#change tank number to factor
b.length$PITtag = as.factor(b.length$PITtag)#change tag number to factors
length.lmer = lmer(value~day*treatment + (1|tank/PITtag), data = b.length, REML = F) #full model
length.lmer2 = lmer(value~day*treatment + (1|PITtag), data = b.length, REML = F)#model without tank affects
length.lm = lm(value~day*treatment, data=b.length)#simplest model
anova(length.lmer,length.lmer2)#compare full model and wihtout tanks to see if tank is significant
anova(length.lmer2,length.lm)#compare model wihtouttanks to simple model to see PITtag affect
summary(length.lmer2) #no tank effect so excuded but PITtag significant
anova.length.lmer = anova(length.lmer2)
capture.output(anova.length.lmer,file='anova.length.lmer.csv')

#CF
ggplot(dat.cf, aes(x = day, y =cf, group = PITtag, color = treatment))+
  geom_line()+ylab('CF')+
  scale_x_continuous(expand = c(0,0),limits = c(0,98),breaks = c(0,24,48,75,91,98))+ 
  theme_classic()#plot masses over time for each fish
mean.cf =  group_by(dat.cf,day,treatment) %>% summarise(average_cf=mean(cf),sem_cf = sd(cf)/sqrt(length(cf)),n=length(cf))
cf.rmANOVA.plot = 
  ggplot(mean.cf, aes(x = day, y = average_cf, color = treatment))+
  geom_line()+ylab('length (cm)')+
  geom_point()+
  geom_errorbar(aes(ymin=average_cf - sem_cf, ymax = average_cf + sem_cf), width = 1.5) +
  scale_x_continuous(expand = c(0,0),limits = c(-10,108),breaks = c(0,24,48,75,91,98))+ 
  theme_classic()#plot masses over time for each fish
ggsave(filename = 'cf.rmANOVA.pdf')
dat.cf$day = as.factor(dat.cf$day)#change days to factors
dat.cf$tank = as.factor(dat.cf$tank)#change tank number to factor
dat.cf$PITtag = as.factor(dat.cf$PITtag)#change tag number to factors
#lmer doesn't work for some reason so use lme instead to get p-values
library(nlme)
cf.lme = lme(cf~day*treatment, random=~1|tank/PITtag, data = dat.cf) #test tank effect
cf.lme2 = lme(cf~day*treatment, random=~1|PITtag, data = dat.cf)#model with no tank effect
anova(cf.lme, cf.lme2)#compare models to see tank effect
anova.cf.lme = anova(cf.lme)
capture.output(anova.cf.lme,file = 'anova.cf.lme.csv')


#looking at normality
ggplot(b.mass, aes(log(value), color = treatment))+geom_histogram()
shapiro.test((b.length$value))
qqnorm(log(b.mass$value))
qqline(log(b.mass$value), col='red')
library(nortest)
ad.test(log(b.length$value))
qqnorm(log(b.length$value))
