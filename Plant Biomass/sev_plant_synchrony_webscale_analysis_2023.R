####### Goal: Analysis of synchrony and stability of plant communities at Sevilleta between sites/ecosystems/treatments
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
sync<-read.csv("asynchrony_stability_allwebs.csv",stringsAsFactors = T)


### does the stability of productivity (total biomass)
### decline with greater synchrony in the temporal dynamics of plant species? --------------

#filter data to controls
# only incl quads with >10 y
sync_controls<-sync
summary(sync_controls$year_count)

# look at response variable
hist(sync_controls$biomass_stability) # has better distribution than biomass_stability2
hist(sync_controls$biomass_stability2)
summary(sync_controls$biomass_stability)
summary(sync_controls$web_ID)
summary(sync_controls$year_count)

# SIMPLE LINEAR MODEL - Loreau metric
m.loreau<-lm(biomass_stability ~ synchrony*season*ecosystem, data=sync_controls)
Anova(m.loreau,type=2)

# quick visuals
visreg(m.loreau, "synchrony", by="season")
visreg(m.loreau, "synchrony", by="ecosystem")

# estimate and compare slopes
emtrends(m.loreau, var="synchrony", ~ecosystem|season)
emtrends(m.loreau, var="synchrony", ~ecosystem)
#pairs(emtrends(m.loreau, var="synchrony",  ~ecosystem|season))

# assumption 1 - normality of residuals
hist(resid(m.loreau))
qqnorm(resid(m.loreau))

# assumption 2 - homogeneity of variances
plot(m.loreau)

# graphic LOREAU
sync2 <- sync_controls %>% dplyr::mutate(ecosystem = factor(ecosystem, levels=c("Pinon woodland", "Plains grassland", "Desert grassland", "Desert shrubland"))) 
summary(sync2$ecosystem)
g.loreau.all<-ggplot(data=sync2, aes(x=synchrony, y=biomass_stability, fill=ecosystem))+
  facet_grid(rows=vars(ecosystem))+
  geom_point(aes(color=ecosystem))+
  geom_smooth(aes(group=ecosystem),formula=y~x,color="black",method="lm")+
  scale_color_manual(values=c("darkslategrey", "steelblue3", "grey30","olivedrab4"))+
  scale_fill_manual(values=c("darkslategrey", "steelblue3", "grey30","olivedrab4"))+
  xlab("Synchrony (Loreau metric)")+
  ylab(bquote('Stability in aboveground production (g  '*~m^-2~y^-1*')'))+
  theme_minimal()+
  theme(legend.position = "none")
g.loreau.all          
ggsave("loreau_web_all.jpg",dpi=500,width=5,height=8)


g.loreau.all<-ggplot(data=sync2, aes(x=synchrony, y=biomass_stability, fill=ecosystem))+
  #facet_grid(rows=vars(ecosystem))+
  geom_point(aes(color=ecosystem))+
  geom_smooth(aes(color=ecosystem),formula=y~x,method="lm")+
  scale_color_manual(values=c("darkslategrey", "steelblue3", "grey30","olivedrab4"))+
  scale_fill_manual(values=c("darkslategrey", "steelblue3", "grey30","olivedrab4"))+
  xlab("Synchrony (Loreau metric)")+
  ylab(bquote('Stability in aboveground production (g  '*~m^-2~y^-1*')'))+
  theme_minimal()+
  theme(legend.position = "none")
g.loreau.all          
ggsave("loreau_web_all_no_facet.jpg",dpi=500,width=8,height=8)


############# Gross metric
# SIMPLE LINEAR MODEL - gross metric
m.gross<-lm(biomass_stability ~ synchrony.Gross*season*ecosystem, data=sync_controls)
Anova(m.gross,type=2)

# quick visuals
visreg(m.gross, "synchrony.Gross", by="season")
visreg(m.gross, "synchrony.Gross", by="ecosystem")

# estimate and compare slopes
emtrends(m.gross, var="synchrony.Gross", ~ecosystem|season)
emtrends(m.gross, var="synchrony.Gross", ~ecosystem)
pairs(emtrends(m.gross, var="synchrony.Gross",  ~ecosystem|season))

# assumption 1 - normality of residuals
hist(resid(m.gross))
qqnorm(resid(m.gross))

# assumption 2 - homogeneity of variances
plot(m.gross)

# graphic LOREAU
sync2 <- sync_controls %>% dplyr::mutate(ecosystem = factor(ecosystem, levels=c("Pinon woodland", "Plains grassland", "Desert grassland", "Desert shrubland"))) 
summary(sync2$ecosystem)
g.loreau.all<-ggplot(data=sync2, aes(x=synchrony.Gross, y=biomass_stability, fill=ecosystem))+
  facet_grid(rows=vars(ecosystem))+
  geom_point(aes(color=ecosystem))+
  geom_smooth(aes(group=ecosystem),formula=y~x,color="black",method="lm")+
  scale_color_manual(values=c("darkslategrey", "steelblue3", "grey30","olivedrab4"))+
  scale_fill_manual(values=c("darkslategrey", "steelblue3", "grey30","olivedrab4"))+
  xlab("Synchrony (Gross metric)")+
  ylab(bquote('Stability in aboveground production (g  '*~m^-2~y^-1*')'))+
  theme_minimal()+
  theme(legend.position = "none")
g.loreau.all          
ggsave("gross_web_all.jpg",dpi=500,width=5,height=8)
