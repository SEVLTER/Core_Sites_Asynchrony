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
sync<-read.csv("asynchrony_stability_allquads.csv",stringsAsFactors = T)


### does the stability of productivity (total biomass)
### decline with greater synchrony in the temporal dynamics of plant species? --------------

#filter data to controls
# only incl quads with >10 y
sync_controls<-subset(sync, treatment=="C")
summary(sync_controls$year_count)

# look at response variable
hist(sync_controls$biomass_stability) # has better distribution than biomass_stability2
hist(sync_controls$biomass_stability2)
summary(sync_controls$biomass_stability)

# SIMPLE LINEAR MODEL - Loreau metric
m.loreau<-lm(biomass_stability ~ synchrony*season*ecosystem, data=sync_controls)
Anova(m.loreau,type=2)

# quick visuals
visreg(m.loreau, "synchrony", by="season")
visreg(m.loreau, "synchrony", by="ecosystem")

# estimate and compare slopes
emtrends(m.loreau, var="synchrony", ~ecosystem|season)
emtrends(m.loreau, var="synchrony", ~ecosystem)
pairs(emtrends(m.loreau, var="synchrony",  ~ecosystem|season))

# assumption 1 - normality of residuals
hist(resid(m.loreau))
qqnorm(resid(m.loreau))

# assumption 2 - homogeneity of variances
plot(m.loreau)

# MIXED EFFECTS MODEL - Loreau metric
#how are ecosystems represented by sites?
summary(sync_controls$ecosystem:sync_controls$site)
#  Plains grassland | core = 60, EDGE = 20
#  Plains-Desert grassland ecotone | warming = 20, NutNet = 10, tower_east= 20, mixed_grass = 60, fertilizer = 80

# mixed model
m.mixed.loreau <- lme(biomass_stability ~ synchrony*season*ecosystem, random = ~1|site, data = sync_controls)
Anova(m.mixed.loreau, type=2)

# quick visuals
visreg(m.mixed.loreau, "synchrony", by="season")
visreg(m.mixed.loreau, "synchrony", by="ecosystem")

# estimate and compare slopes
emtrends(m.mixed.loreau, var="synchrony", ~ecosystem)
emtrends(m.mixed.loreau, var="synchrony", ~ecosystem|season)
pairs(emtrends(m.mixed.loreau, var="synchrony",  ~ecosystem|season))
pairs(emtrends(m.mixed.loreau, var="synchrony",  ~season|ecosystem))

# assumption 1 - normality of residuals
hist(resid(m.mixed.loreau))
qqnorm(resid(m.mixed.loreau))

# assumption 2 - homogeneity of variances
plot(m.mixed.loreau)

# graphic
g.loreau<-ggplot(data=sync_controls, aes(x=synchrony, y=biomass_stability, group=ecosystem))+
          facet_grid(rows=vars(ecosystem),cols=vars(season))+
            geom_point()+
            geom_smooth(aes(group=ecosystem),formula=y~x,color="black",method="lm")+
            #scale_color_manual(values=c("olivedrab4","grey30","khaki3","steelblue3"))+
            xlab("Synchrony (Loreau metric)")+
            ylab(bquote('Stability in aboveground production (g  '*~m^-2~y^-1*')'))+
            theme_minimal()+
           theme(legend.position = "none")
g.loreau          
ggsave("controls_loreau_raw.jpg",dpi=500,width=4,height=10)

#only fall
sync_controls_fall<-subset(sync_controls, season=="fall")
summary(sync_controls_fall$ecosystem)
sync_controls_fall2 <- sync_controls_fall %>% dplyr::mutate(ecosystem = factor(ecosystem, levels=c("Pinon woodland", "Plains grassland", "Plains-Desert grassland ecotone", "Desert grassland", "Desert grass-shrub ecotone", "Desert shrubland"))) 

      
g.loreau.fall<-ggplot(data=sync_controls_fall2, aes(x=synchrony, y=biomass_stability, fill=ecosystem))+
  facet_grid(rows=vars(ecosystem))+
  geom_point(aes(color=ecosystem))+
  geom_smooth(aes(group=ecosystem),formula=y~x,color="black",method="lm")+
  scale_color_manual(values=c("darkslategrey", "steelblue3", "lightblue4","grey70","darkgreen", "olivedrab4"))+
  scale_fill_manual(values=c("darkslategrey", "steelblue3", "lightblue4","grey70","darkgreen", "olivedrab4"))+
  xlab("Synchrony (Loreau metric)")+
  ylab(bquote('Stability in aboveground production (g  '*~m^-2~y^-1*')'))+
  theme_minimal()+
  theme(legend.position = "none")
g.loreau.fall          
ggsave("controls_loreau_raw_fall.jpg",dpi=500,width=4,height=12)


sync_controls2 <- sync_controls %>% dplyr::mutate(ecosystem = factor(ecosystem, levels=c("Pinon woodland", "Plains grassland", "Plains-Desert grassland ecotone", "Desert grassland", "Desert grass-shrub ecotone", "Desert shrubland"))) 
summary(sync_controls2$ecosystem)
g.loreau.all<-ggplot(data=sync_controls2, aes(x=synchrony, y=biomass_stability, fill=ecosystem))+
  facet_grid(rows=vars(ecosystem))+
  geom_point(aes(color=ecosystem))+
  geom_smooth(aes(group=ecosystem),formula=y~x,color="black",method="lm")+
  scale_color_manual(values=c("darkslategrey", "steelblue3", "lightblue4","grey70","darkgreen", "olivedrab4"))+
  scale_fill_manual(values=c("darkslategrey", "steelblue3", "lightblue4","grey70","darkgreen", "olivedrab4"))+
  xlab("Synchrony (Loreau metric)")+
  ylab(bquote('Stability in aboveground production (g  '*~m^-2~y^-1*')'))+
  theme_minimal()+
  theme(legend.position = "none")
g.loreau.all          
ggsave("controls_loreau_raw_all.jpg",dpi=500,width=4,height=12)

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

# MIXED EFFECTS MODEL - gross metric
#how are ecosystems represented by sites?
summary(sync_controls$ecosystem:sync_controls$site)
#  Plains grassland | core = 60, EDGE = 20
#  Plains-Desert grassland ecotone | warming = 20, NutNet = 10, tower_east= 20, mixed_grass = 60, fertilizer = 80

# mixed model
m.mixed.gross <- lme(biomass_stability ~ synchrony.Gross*season*ecosystem, random = ~1|site, data = sync_controls)
Anova(m.mixed.gross, type=2)
emtrends(m.mixed.gross, var="synchrony.Gross", ~ecosystem)

g.gross.all<-ggplot(data=sync_controls2, aes(x=synchrony.Gross, y=biomass_stability, fill=ecosystem))+
  facet_grid(rows=vars(ecosystem))+
  geom_point(aes(color=ecosystem))+
  geom_smooth(aes(group=ecosystem),formula=y~x,color="black",method="lm")+
  scale_color_manual(values=c("darkslategrey", "steelblue3", "lightblue4","grey70","darkgreen", "olivedrab4"))+
  scale_fill_manual(values=c("darkslategrey", "steelblue3", "lightblue4","grey70","darkgreen", "olivedrab4"))+
  xlab("Synchrony (Gross metric)")+
  ylab(bquote('Stability in aboveground production (g  '*~m^-2~y^-1*')'))+
  theme_minimal()+
  theme(legend.position = "none")
g.gross.all          
ggsave("controls_gross_raw_all.jpg",dpi=500,width=4,height=12)
