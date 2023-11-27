####### Goal: Calculate synchrony of communities at Sevilleta between sites/ecosystems/treatments
### J Rudgers August 2023
rm(list=ls(all=TRUE)) 

# load libraries
library(nlme)
library(tidyverse)
library(plyr)
library(codyn)
library(emmeans)  
library(reshape2)  
library(vegan)
library(car)
library(visreg)

# data steps
# Set working directory
setwd("C:/Users/jrcassvc/Desktop/SEV Data/Asynchrony/")

################ Set up data steps -----------------------
# import data from EDI 
# Or use this function to import
mass<-read.csv("C:/Users/jrcassvc/Desktop/SEV Data/NPP/Sevilleta_allbiomass_05Feb2023.csv",stringsAsFactors = T) #
str(mass)

# create unique identifier for quadrat - need to adjust this line based on the focal site or experiment
# for core sites need site, web, plot, quad to ID the quadrat - this includes everything, clunky, but comprehensive
mass$quad_ID<-as.factor(paste(mass$site,mass$web,mass$transect,mass$treatment,mass$block,mass$plot,mass$subplot,mass$quad,sep="_"))

# cast the data to create columns for each plant species so zeros are added for dates where species was absent  
# for all quads
mass_species<-dcast(mass, year +site +season +season.precip +GDD +SPEI.comp +SiteCluster +MetStation +web +transect +treatment +block +plot +subplot +quad +quad_ID ~ kartez, sum, value.var="biomass.BM", fill=0)

# get rid of column called EMPTY - this was a placeholder for quadrats that had bare ground, not a real species
mass_species <- mass_species %>% select(-EMPTY)
# combine two components of estimated LATR2 biomass
mass_species$LATR2 <-mass_species$LATR2+mass_species$STEM
mass_species <- mass_species %>% select(-STEM)
summary(mass_species)

# calculate total biomass
mass_species$totmass<-rowSums(mass_species[,17:338])

# add richness, diversity, evenness
mass_species$richness<-specnumber(mass_species[,17:338])
mass_species$shannonH<-diversity(mass_species[,17:338])
mass_species$evenness<-mass_species$shannonH/log(mass_species$richness,1)

# create a continuous time variable
mass_species$seas<-as.numeric(as.character(recode_factor(mass_species$season,fall="0.5",spring="0")))
summary(mass_species$seas)
mass_species$year_seas<-as.factor(mass_species$year+mass_species$seas)                            
summary(mass_species$year_seas)
mass_species$time<-as.numeric(as.factor(mass_species$year_seas))
summary(mass_species$time)
mass_species$year.f<-as.factor(mass_species$year)

########## ANPP Biomass Stability Data Steps ----------------------------
totmass<-mass_species %>% select(season,seas,site,SiteCluster,MetStation,web,transect,treatment,block,plot,subplot,quad,quad_ID,year.f,totmass)
totmass_cast<-dcast(totmass, season+seas+site+SiteCluster+MetStation+web+transect+treatment+block+plot+subplot+quad+quad_ID ~ year.f, sum, value.var="totmass",fill=-9999)

# fix problem with fill in dcast
totmass_cast[totmass_cast== "-9999"] <- "NA"
totmass_cast <- totmass_cast %>% mutate_if(is.character, as.numeric) #warning is not important
summary(totmass_cast)

# calculate totmass sd, mean, CV
totmass_cast$sd<-apply(totmass_cast[,14:37], 1, sd, na.rm=TRUE)
totmass_cast$mu<-apply(totmass_cast[,14:37], 1, mean, na.rm=TRUE)
totmass_cast$CV<-totmass_cast$sd/totmass_cast$mu
totmass_cast$biomass_stability<-1/(totmass_cast$CV)
totmass_cast$biomass_stability2<-1/((totmass_cast$CV)^2)

# possible number of years (to 2023) = 24
# so you can filter dataset later to only quads that have a minimum number of years recorded
totmass_cast$year_count<-rowSums(!is.na(totmass_cast[,14:37]))
summary(totmass_cast)
# 69 quadrats of 3568 (2%) with zero mean biomass, these have NAs for CV and stability metrics

##########  Asynchrony Data Steps -----------------------------
# melt to make a row for each species in each quad, differs from original datafile because now includes zeros
mass_species_melt<-melt(mass_species,id=c("season","seas","site", "SiteCluster","MetStation","year","year.f","time","year_seas","season.precip","GDD","SPEI.comp",
                                          "web","transect","treatment","block","plot","subplot","quad","quad_ID"),variable.name="kartez",value.name="biomass.BM")
summary(mass_species_melt$biomass.BM)
summary(mass_species_melt)
levels(mass_species_melt$kartez)
# will need to filter out totmass, richness, shannonH, evenness if we just want the species matrix

# remove unknowns, they will not be consistent in their ID through time and likely appear only once
mass_species_melt<-filter(mass_species_melt, kartez!="UKN1" & kartez!="UNK" & kartez!="UNK 1"  & kartez!="UNK_516_MV_C" & kartez!="UNK1" & kartez!="UNK11" & kartez!="UNK12"
                        &  kartez!="UNK13" & kartez!="UNK2" & kartez!="UNK3" & kartez!="UNK4" & kartez!="UNKASTER" & kartez!="UNKFORB1" 
                        &  kartez!="UNKN" & kartez!="UNKN1" & kartez!="UNKN2" & kartez!="UNKNOWN" & kartez!="UNKNOWN 2" & kartez!="UNKNOWN BRASS" & kartez!="UNKOWN1" & kartez!="UNKSHRUB")

# now recast data to make each row a species X quadrat X season combination, with columns for all years
# this step may be slow
mass_species_cast<-dcast(mass_species_melt, season+seas+site+SiteCluster+MetStation+web+transect+treatment+block+plot+subplot+quad+quad_ID+kartez ~ year.f, sum, value.var = "biomass.BM",fill=-9999)
# fix problem with fill in dcast
mass_species_cast[mass_species_cast== "-9999"] <- "NA"
mass_species_cast <- mass_species_cast %>% mutate_if(is.character, as.numeric) #warning is not important
summary(mass_species_cast)

#remove all species that were never observed in a quad
mass_species_cast$sum<-rowSums(mass_species_cast[,15:38],na.rm=T)
summary(mass_species_cast$sum)
mass_species_stability <-mass_species_cast %>%  filter(sum>0)

#count up number of years of observations for each species/quad/season (row)
mass_species_stability$years_count<-rowSums(!is.na(mass_species_stability[,15:38]))
summary(mass_species_stability$years_count)

# calculate mean long-term species biomass separately by season, also sd in species mass and CV
mass_species_stability$species_sd<-apply(mass_species_stability[,15:38], 1, sd, na.rm=TRUE)
mass_species_stability$species_mu<-apply(mass_species_stability[,15:38], 1, mean, na.rm=TRUE)
mass_species_stability$species_CV<-mass_species_stability$species_sd/mass_species_stability$species_mu

#calculate stability as CV and 1/CV2
mass_species_stability$species_stability<-1/(mass_species_stability$species_CV)
mass_species_stability$species_stability2<-1/((mass_species_stability$species_CV)^2)

# calculate total richness over all years so that synchrony analysis can exclude quads with only ever 1 species
str(mass_species_stability)
levels(mass_species_stability$kartez)

richness_stability <- data.frame(mass_species_stability  %>% filter(kartez=="richness") 
                                 %>%  select(quad_ID, season, mu_richness=species_mu, quad_richness=sum/years_count, CV_quad_richness=species_CV))       
summary(richness_stability)
levels(richness_stability$quad_ID)

msrs <- merge(mass_species_stability, richness_stability, by.x=c("quad_ID","season"), by.y=c("quad_ID","season"), all.x=T)
summary(msrs)
write.csv(msrs, "plants_mass_richness_stability.csv")

# prepare data for asynchrony metrics
# remove quads where only one species was present over the full time series
# remove quads with <10 years of observation
levels(msrs$kartez)
species_async <- as.data.frame( msrs  %>% filter(quad_richness>1 & years_count>10 & kartez!="totmass" & kartez!="richness" & kartez!="shannonH" & kartez!="evenness") )
summary(species_async$kartez)
species_async <- species_async[,1:38]
# there were only 6 quads that had only 1 species through time

# pivot to long format
SEV_plant <- pivot_longer(species_async, cols = '1999':'2022', names_to = "year", values_to = "abundance")
summary(SEV_plant$site)
summary(SEV_plant$kartez)


#create a variable for ecosystem type
SEV_plant$ecosystem<-recode_factor(SEV_plant$site, 
                                  core_black="Desert grassland", 
                                  core_blue="Plains grassland",
                                  core_creosote="Desert shrubland",
                                  core_PJ="Pinon woodland",
                                  crust_creosote="Desert grass-shrub ecotone",
                                  crust_grass="Plains-Desert ecotone",
                                  crust_PJ="Pinon-juniper ecotone",
                                  EDGE_black ="Desert grassland",     
                                  EDGE_blue="Plains grassland",
                                  fertilizer="Plains-Desert grassland ecotone",
                                  grassland_burn="Plains-Desert grassland ecotone",
                                  grid_black="Desert grassland",  
                                  grid_blue ="Plains grassland",
                                  grid_creosote ="Desert shrubland",      
                                  iso_web="Desert grass-shrub ecotone",
                                  JSAV="Juniper savanna",
                                  meanvar_black="Desert grassland",
                                  meanvar_blue ="Plains grassland",  
                                  meanvar_creosote="Desert shrubland",    
                                  meanvar_jsav="Juniper savanna",       
                                  meanvar_pj="Pinon woodland",      
                                  mixed_grass="Plains-Desert grassland ecotone",    
                                  mixed_shrub="Desert grass-shrub ecotone",             
                                  MRME="Desert grassland",
                                  NutNet="Plains-Desert grassland ecotone",         
                                  PJControl="Pinon woodland",         
                                  PJGirdle ="Pinon-juniper ecotone",      
                                  tower_east="Plains-Desert grassland ecotone",       
                                  tower_west="Desert grassland",          
                                  warming="Plains-Desert grassland ecotone")  

quad_key <- SEV_plant%>%
            dplyr::select(site, SiteCluster, ecosystem, treatment, block, plot, subplot, web, transect, quad_ID)%>%
            unique()


SEV_spring <- subset(SEV_plant, season == "spring")
SEV_spring$year<-as.numeric(SEV_spring$year)
SEV_fall <- subset(SEV_plant, season == "fall")
SEV_fall$year<-as.numeric(SEV_fall$year)

### Synchrony Metric Calculations #####################################################################
#SPRING
spring_synch <- synchrony(
  df = SEV_spring,
  time.var = "year",
  species.var = "kartez",
  abundance.var = "abundance",
  metric = "Loreau", #"Gross" 
  replicate.var = "quad_ID"
)%>%
  left_join(quad_key, by = "quad_ID")

mod <- lme(synchrony~ecosystem, random = ~1|site, data = spring_synch)
anova(mod)
a1 <- aov(synchrony~site, data = spring_synch)
TukeyHSD(x = a1, 'site', conf.level=0.95)
emmeans(mod, list(pairwise ~ ecosystem), adjust = "tukey")

spring_synch_summary <- spring_synch%>%
  ddply(.(site),function(x)data.frame(
    synchrony.mean = mean(x$synchrony),
    synchrony.confint = qt(0.975, df=length(x$synchrony)-1)*sd(x$synchrony, na.rm = TRUE)/sqrt(length(x$synchrony)-1),
    synchrony.se = sd(x$synchrony, na.rm = TRUE)/sqrt(length(x$synchrony))
  ))

ggplot(spring_synch_summary, aes(site, synchrony.mean))+
geom_pointrange(aes(ymin = synchrony.mean-synchrony.confint, ymax = synchrony.mean+synchrony.confint), size = 1.5)+
theme_bw()

# SPRING merge quadrat scale async metric with stability (by quadrat)
stabil_spring<-totmass_cast %>% filter(season=="spring" & year_count>10)  %>%  select(quad_ID, sd, mu, CV, biomass_stability, biomass_stability2, year_count) 
async_stabil_spring<-merge(spring_synch,stabil_spring, by.x="quad_ID", by.y="quad_ID" )
write.csv(async_stabil_spring, "async_stabil_spring_allquads.csv")
# spring: 659 quadrats with >1 species and at least 10 years of data
async_stabil_spring_controls<-async_stabil_spring %>% filter(treatment=="C")
# N=368 quadrats were controls


### SPRING FLORA: does stability of biomass decline with greater synchrony in the temporal dynamics of plant species?
hist(async_stabil_spring$biomass_stability) # has better distribution than biomass_stability2
m.spring<-lm(biomass_stability ~ synchrony*ecosystem, data=async_stabil_spring_controls)
Anova(m.spring,type=2)
visreg(m.spring, "synchrony", by="ecosystem")
emtrends(m.spring, var="synchrony", ~ecosystem)
pairs(emtrends(m.spring, var="synchrony",  ~ecosystem))


#FALL
fall_synch <- synchrony(
  df = SEV_fall,
  time.var = "year",
  species.var = "kartez",
  abundance.var = "abundance",
  metric = "Loreau", #"Gross" 
  replicate.var = "quad_ID"
)%>%
  left_join(quad_key, by = "quad_ID")

mod <- lme(synchrony~ecosystem, random = ~1|site, data = fall_synch)
anova(mod)

fall_synch_summary <- fall_synch%>%
  ddply(.(site),function(x)data.frame(
    synchrony.mean = mean(x$synchrony),
    synchrony.confint = qt(0.975, df=length(x$synchrony)-1)*sd(x$synchrony, na.rm = TRUE)/sqrt(length(x$synchrony)-1),
    synchrony.se = sd(x$synchrony, na.rm = TRUE)/sqrt(length(x$synchrony))
  ))

ggplot(fall_synch_summary, aes(site, synchrony.mean))+
geom_pointrange(aes(ymin = synchrony.mean-synchrony.confint, ymax = synchrony.mean+synchrony.confint), size = 1.5)+
theme_bw()


# fall merge quadrat scale async metric with stability (by quadrat)
stabil_fall<-totmass_cast %>% filter(season=="fall" & year_count>10)  %>%  select(quad_ID, sd, mu, CV, biomass_stability, biomass_stability2, year_count) 
async_stabil_fall<-merge(fall_synch,stabil_fall, by.x="quad_ID", by.y="quad_ID" )
write.csv(async_stabil_fall, "async_stabil_fall_allquads.csv")


# fall: 739 quadrats with > 1 species and at least 10 years of data
async_stabil_fall_controls<-async_stabil_fall %>% filter(treatment=="C")
# N=370 quadrats were controls



#### Rank Abundance Change #####################################################################
#FALL
fall_RAC <- SEV_fall%>%
  subset(abundance != 0)%>%
  RAC_change(
  time.var = "year",
  species.var = "kartez",
  abundance.var = "abundance",
  replicate.var = "quad_ID",
  reference.time = NULL
)%>%
  left_join(quad_key, by = "quad_ID")

#richness change
mod <- lme(richness_change~ecosystem, random = ~1 |site/quad_ID, data = fall_RAC)
anova(mod)
emmeans(mod, list(pairwise ~ ecosystem), adjust = "tukey")

#rank change
mod <- lme(rank_change~ecosystem, random = ~1 | site/quad_ID, data = fall_RAC)
anova(mod)
emmeans(mod, list(pairwise ~ ecosystem), adjust = "tukey")


#SPRING
spring_RAC <- SEV_spring%>%
  subset(abundance != 0)%>%
  RAC_change(
    time.var = "year",
    species.var = "kartez",
    abundance.var = "abundance",
    replicate.var = "quad_ID",
    reference.time = NULL
  )%>%
  left_join(quad_key, by = "quad_ID")

#richness change
mod <- lme(richness_change~ecosystem, random = ~1|site/quad_ID, data = spring_RAC)
anova(mod)
emmeans(mod, list(pairwise ~ ecosystem), adjust = "tukey")

#rank change
mod <- lme(rank_change~ecosystem, random = ~1|site/quad_ID, data = spring_RAC)
anova(mod)
emmeans(mod, list(pairwise ~ ecosystem), adjust = "tukey")


##### subset for TOWERS - do with >5 years of data ####################################
towers_async <- as.data.frame( msrs  %>% filter(site=="tower_east"|site=="tower_west")
                                 %>% filter(quad_richness>1 & years_count>5 & kartez!="totmass" 
                                                & kartez!="richness" & kartez!="shannonH" & kartez!="evenness"))

summary(towers_async$site)

towers_async <- towers_async[,1:38]

# pivot to long format
SEV_towers <- pivot_longer(towers_async, cols = '1999':'2022', names_to = "year", values_to = "abundance")
summary(SEV_towers$site)



#create a variable for ecosystem type
SEV_towers$ecosystem<-recode_factor(SEV_towers$site, 
                                   tower_east="Plains-Desert grassland ecotone",       
                                   tower_west="Desert grassland")  

quad_key_towers <- SEV_towers%>%
  dplyr::select(site, SiteCluster, ecosystem, treatment, transect, quad_ID)%>%
  unique()


SEV_spring_towers <- subset(SEV_towers, season == "spring")
SEV_spring_towers$year<-as.numeric(SEV_spring_towers$year)
SEV_fall_towers <- subset(SEV_towers, season == "fall")
SEV_fall_towers$year<-as.numeric(SEV_fall_towers$year)

### TOWERS   Synchrony Metric Calculations #####################################################################
#SPRING
spring_synch_towers <- synchrony(
  df = SEV_spring_towers,
  time.var = "year",
  species.var = "kartez",
  abundance.var = "abundance",
  metric = "Loreau", #"Gross" 
  replicate.var = "quad_ID"
)%>%
  left_join(quad_key_towers, by = "quad_ID")

# SPRING merge quadrat scale async metric with stability (by quadrat)
stabil_spring_towers<-totmass_cast %>% filter(season=="spring" & year_count>5)  %>%  
  filter(site=="tower_east"|site=="tower_west")  %>%  select(quad_ID, sd, mu, CV, biomass_stability, biomass_stability2, year_count) 
async_stabil_spring_towers<-merge(spring_synch_towers,stabil_spring_towers, by.x="quad_ID", by.y="quad_ID" )
# spring towers: 40 quadrats with >1 species and at least 5 years of data
write.csv(async_stabil_spring_towers, "async_stabil_spring_towers.csv")


#FALL
fall_synch_towers <- synchrony(
  df = SEV_fall_towers,
  time.var = "year",
  species.var = "kartez",
  abundance.var = "abundance",
  metric = "Loreau", #"Gross" 
  replicate.var = "quad_ID"
)%>%
  left_join(quad_key_towers, by = "quad_ID")


# fall merge quadrat scale async metric with stability (by quadrat)
stabil_fall_towers<-totmass_cast %>% filter(season=="fall" & year_count>5)  %>%  
  filter(site=="tower_east"|site=="tower_west")  %>%  select(quad_ID, sd, mu, CV, biomass_stability, biomass_stability2, year_count) 
async_stabil_fall_towers<-merge(fall_synch_towers,stabil_fall_towers, by.x="quad_ID", by.y="quad_ID" )
# fall towers: 40 quadrats with >1 species and at least 5 years of data
write.csv(async_stabil_fall_towers, "async_stabil_fall_towers.csv")




