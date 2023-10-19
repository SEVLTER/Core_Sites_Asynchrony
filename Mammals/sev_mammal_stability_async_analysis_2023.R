####### Goal: Analysis of stability of mammals
####### with asynchrony of plants at Sevilleta between ecosystems
### J Rudgers October 2023

# give R a blank slate
rm(list=ls(all=TRUE)) 

# load libraries
library(nlme)
library(tidyverse)
library(plyr)
library(codyn)
library(emmeans)  
library(reshape2)  
library(car)
library(visreg)

# Set working directory
setwd("C:/Users/jrcassvc/Desktop/SEV Data/Asynchrony/")

# Read in data
species_sync<-read.csv("mammal_species_async_stability.csv",stringsAsFactors = T)
total_sync<-read.csv("mammal_total_async_stability.csv",stringsAsFactors = T)

#######  Q: Are mammals more stable with greater plant asynchrony? #######
hist(total_sync$abund_stability09)
summary(total_sync$abund_stability09)
hist(total_sync$abund_stability09.2)
summary(total_sync$abund_stability09.2) # has greater range... use this metric

### LOREAU METRIC
# test synchrony alone
# n.s.
m.tot.Loreau<-lm(abund_stability09.2 ~ synchrony,data=total_sync)
Anova(m.tot.Loreau,type=2)
# does synchrony effect vary between ecosystems - low sample size!
# no.
m.tot.Loreau.ecosystem<-lm(abund_stability09.2 ~ synchrony*ecosystem,data=total_sync)
Anova(m.tot.Loreau.ecosystem,type=2)
# does synchrony effect vary between seasons - low sample size!
# no.
m.tot.Loreau.season<-lm(abund_stability09.2 ~ synchrony*season,data=total_sync)
Anova(m.tot.Loreau.season,type=2)

# Loreau make graph of all samples
g.tot.Loreau<-ggplot(data=total_sync, aes(x=synchrony, y=abund_stability09.2, color=ecosystem))+
  geom_smooth(formula=y~x,method="lm", color="black")+
  geom_point(aes(color=ecosystem, shape=season),size=3)+
  scale_color_manual(values=c("grey30","olivedrab4"))+
  scale_shape_manual(values = c(15, 16)) + 
  xlab("Synchrony (Loreau metric)")+
  ylab(bquote('Stability in mammals (1/'*~CV^ 2*')'))+
  theme_minimal()+
  theme(legend.position = "right")
g.tot.Loreau       
ggsave("loreau_tot_mammals.jpg",dpi=500,width=8,height=6)

#### GROSS METRIC
# test synchrony alone
# n.s.
m.tot.Gross<-lm(abund_stability09.2 ~ synchrony.Gross,data=total_sync)
Anova(m.tot.Gross,type=2)

# does synchrony effect vary between ecosystems - low sample size!
# yes.
m.tot.Gross.ecosystem<-lm(abund_stability09.2 ~ synchrony.Gross*ecosystem,data=total_sync)
Anova(m.tot.Gross.ecosystem,type=2)
emtrends(m.tot.Gross.ecosystem, var="synchrony.Gross", ~ecosystem)

# Gross metric - by ecosystem - make graph of all samples
g.tot.Gross<-ggplot(data=total_sync, aes(x=synchrony.Gross, y=abund_stability09.2, fill=ecosystem, color=ecosystem))+
 facet_grid(cols=vars(ecosystem), scales="free_x")+
  geom_point(aes(color=ecosystem, shape=season),size=3)+
  geom_smooth(aes(color=ecosystem), formula=y~x,method="lm")+
  scale_color_manual(values=c("grey30","olivedrab4"))+
  scale_fill_manual(values=c("grey30","olivedrab4"))+
  scale_shape_manual(values = c(15, 16)) + 
  xlab("Synchrony (Gross metric)")+
  ylab(bquote('Stability in mammals (1/'*~CV^ 2*')'))+
  theme_minimal()+
  theme(legend.position = "right")
  #ylim(1,2.5)
g.tot.Gross   
ggsave("gross_tot_mammals.jpg",dpi=500,width=8,height=4)

# does synchrony effect vary between seasons - low sample size!
# no.
m.tot.Gross.season<-lm(abund_stability09.2 ~ synchrony.Gross*season,data=total_sync)
Anova(m.tot.Gross.season,type=2)

#######  Q: Are mammals more stable with greater plant biomass stability? #######
hist(total_sync$biomass_stability)
summary(total_sync$biomass_stability)
hist(total_sync$biomass_stability2)
summary(total_sync$biomass_stability2) # has greater range... use this metric

# test stability alone
# n.s.
m.tot.stability<-lm(abund_stability09.2 ~ biomass_stability2,data=total_sync)
Anova(m.tot.stability,type=2)
# does stability effect vary between ecosystems - low sample size!
# no.
m.tot.stability.ecosystem<-lm(abund_stability09.2 ~ biomass_stability2*ecosystem,data=total_sync)
Anova(m.tot.stability.ecosystem,type=2)
# does stability effect vary between seasons - low sample size!
# no.
m.tot.stability.season<-lm(abund_stability09.2 ~ biomass_stability2*season,data=total_sync)
Anova(m.tot.stability.season,type=2)

# make graph of all samples
g.tot.stability<-ggplot(data=total_sync, aes(x=biomass_stability2, y=abund_stability09.2, color=ecosystem))+
  geom_smooth(formula=y~x,method="lm", color="black")+
  geom_point(aes(color=ecosystem, shape=season),size=3)+
  scale_color_manual(values=c("grey30","olivedrab4"))+
  scale_shape_manual(values = c(15, 16)) + 
  xlab(bquote('Plant biomass stability (1/'*~CV^ 2*')'))+
  ylab(bquote('Stability in mammals (1/'*~CV^ 2*')'))+
  theme_minimal()+
  theme(legend.position = "right")
g.tot.stability       
ggsave("plant_biomass_stability_tot_mammals.jpg",dpi=500,width=8,height=6)


#######  Q: Are certain species more stable with greater plant asynchrony? #######
# PGFV
PGFV<-subset(species_sync, species=="PGFV")
summary(PGFV)
###  
# test synchrony alone
# n.s.
m.PGFV.Gross<-lm(species_stability09.2 ~ synchrony.Gross,data=PGFV)
Anova(m.PGFV.Gross,type=2)

# does synchrony effect vary between ecosystems - low sample size!
# yes. Gross
m.PGFV.Gross.ecosystem<-lm(species_stability09.2 ~ synchrony.Gross*ecosystem,data=PGFV)
Anova(m.PGFV.Gross.ecosystem,type=2)
emtrends(m.PGFV.Gross.ecosystem,var="synchrony.Gross",~ecosystem)
# no. Loreau
m.PGFV.Loreau.ecosystem<-lm(species_stability09.2 ~ synchrony*ecosystem,data=PGFV)
Anova(m.PGFV.Loreau.ecosystem,type=2)

# make graph PGFV
g.PGFV.Gross<-ggplot(data=PGFV, aes(x=synchrony.Gross, y=species_stability09.2, fill=ecosystem, color=ecosystem))+
  facet_grid(cols=vars(ecosystem), scales="free_x")+
  geom_point(aes(color=ecosystem, shape=season),size=3)+
  geom_smooth(aes(color=ecosystem), formula=y~x,method="lm")+
  scale_color_manual(values=c("grey30","olivedrab4"))+
  scale_fill_manual(values=c("grey30","olivedrab4"))+
  scale_shape_manual(values = c(15, 16)) + 
  xlab("Synchrony (Gross metric)")+
  ylab(expression(paste("Stability in ",italic("Perognathus flavus")," abundance")))+
  theme_minimal()+
  theme(legend.position = "right")
#ylim(1,2.5)
g.PGFV.Gross   
ggsave("gross_PGFV_ecosystem.jpg",dpi=500,width=7,height=4)

#######  Q: Are certain species more abundant with greater plant asynchrony? #######
# PGFV
PGFV<-subset(species_sync, species=="PGFV")
summary(PGFV)
###  
# test synchrony alone
# yes.
m.PGFV.Gross.abund<-lm(species_mu09 ~ synchrony.Gross,data=PGFV)
Anova(m.PGFV.Gross.abund,type=2)
m.PGFV.Loreau.abund<-lm(species_mu09 ~ synchrony,data=PGFV)
Anova(m.PGFV.Loreau.abund,type=2)

# does synchrony effect vary between ecosystems - low sample size!
# no.
m.PGFV.Gross.ecosystem.abund<-lm(species_mu09 ~ synchrony.Gross*ecosystem,data=PGFV)
Anova(m.PGFV.Gross.ecosystem.abund,type=2)
# yes? Loreau
m.PGFV.Loreau.ecosystem.abund<-lm(species_mu09 ~ synchrony*ecosystem,data=PGFV)
Anova(m.PGFV.Loreau.ecosystem.abund,type=2)

# make graph PGFV abundance
g.PGFV.Gross.abund<-ggplot(data=PGFV, aes(x=synchrony.Gross, y=species_mu09, color=ecosystem))+
 # facet_grid(cols=vars(ecosystem), scales="free_x")+
  geom_point(aes(color=ecosystem, shape=season),size=3)+
  geom_smooth(formula=y~x,method="lm",color="black")+
  scale_color_manual(values=c("grey30","olivedrab4"))+
 # scale_fill_manual(values=c("grey30","olivedrab4"))+
  scale_shape_manual(values = c(15, 16)) + 
  xlab("Synchrony (Gross metric)")+
  ylab(expression(paste(italic("Perognathus flavus")," abundance per km2")))+
  theme_minimal()+
  theme(legend.position = "right")
#ylim(1,2.5)
g.PGFV.Gross.abund  
ggsave("gross_PGFV_abund_ecosystem.jpg",dpi=500,width=7,height=4)



# DIOR
DIOR<-subset(species_sync, species=="DIOR")
summary(DIOR)
###  
# test synchrony alone
# yes?
m.DIOR.Gross.abund<-lm(species_mu09 ~ synchrony.Gross,data=DIOR)
Anova(m.DIOR.Gross.abund,type=2)
# no.
m.DIOR.Loreau.abund<-lm(species_mu09 ~ synchrony,data=DIOR)
Anova(m.DIOR.Loreau.abund,type=2)

# does synchrony effect vary between ecosystems - low sample size!
# no.
m.DIOR.Gross.ecosystem.abund<-lm(species_mu09 ~ synchrony.Gross*ecosystem,data=DIOR)
Anova(m.DIOR.Gross.ecosystem.abund,type=2)
# no. Loreau
m.DIOR.Loreau.ecosystem.abund<-lm(species_mu09 ~ synchrony*ecosystem,data=DIOR)
Anova(m.DIOR.Loreau.ecosystem.abund,type=2)

# make graph DIOR abundance
g.DIOR.Gross.abund<-ggplot(data=DIOR, aes(x=synchrony.Gross, y=species_mu09, color=ecosystem))+
  # facet_grid(cols=vars(ecosystem), scales="free_x")+
  geom_point(aes(color=ecosystem, shape=season),size=3)+
  geom_smooth(formula=y~x,method="lm",color="black")+
  scale_color_manual(values=c("grey30","olivedrab4"))+
  # scale_fill_manual(values=c("grey30","olivedrab4"))+
  scale_shape_manual(values = c(15, 16)) + 
  xlab("Synchrony (Gross metric)")+
  ylab(expression(paste(italic("Dipodomys merriami")," abundance per km2")))+
  theme_minimal()+
  theme(legend.position = "right")
#ylim(1,2.5)
g.DIOR.Gross.abund  
ggsave("gross_DIOR_abund_ecosystem.jpg",dpi=500,width=7,height=4)