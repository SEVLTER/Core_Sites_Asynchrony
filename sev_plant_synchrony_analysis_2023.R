
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
sync_controls<-subset(sync, treatment=="C")

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
            ylab(bquote('Stability in aboveground production (g  '*~m^-3~y^-1*')'))+
            theme_minimal()+
           theme(legend.position = "none")
g.loreau          
ggsave("controls_loreau_raw.jpg",dpi=500,width=4,height=10)
