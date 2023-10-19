###### DATA WRANGLING SCRIPT FOR SEV RODENTS ######
###### Date: 19 November 2021, revised 18 October 2023 ######
rm(list=ls(all=TRUE)) #give R a blank slate

library(tidyverse)
library(reshape2)
library(vegan)
#set working directory --> adjust to your computer
setwd("C:/Users/jrcassvc/Desktop/SEV Data/Asynchrony/")

######  SEV MAMMAL DATA WRANGLING -----------------------------------
#read in mammal data: 
SEV<-read.csv("sev008_rodentpopns_20231002.csv",stringsAsFactors = T)
summary(SEV)
summary(SEV$species)

# convert variables to factor
SEV <- SEV %>% 
  mutate(season = as.factor(season),web = as.factor(web))
# fix issue with season = 0 - the next line finds those
summary(as.factor(SEV$year):SEV$season)

#exclude recaps 
SEV<-subset(SEV,recap=="N")

### FURTHER SUBSET OPTIONS - some are currently disabled ###
#exclude trap = NA (these are cases when trap number was not written down)
#SEV<-subset(SEV,!(is.na(trap)))

#cull to summer/fall samples only 
summary(SEV$season)
SEV<-subset(SEV,season!="0")
summary(SEV)

#remove squirrels 
SEV<-subset(SEV,species!="SPSP")
SEV<-subset(SEV,species!="AMIN")
SEV<-subset(SEV,species!="AMLE")

#remove chipmunks 
SEV<-subset(SEV,species!="EUQU")
SEV<-subset(SEV,species!="EUDO")

#remove rabbits and birds and lizard
SEV<-subset(SEV,species!="BUNNY"&species!="SYAU")
SEV<-subset(SEV,species!="BIRD")
SEV<-subset(SEV,species!="LIZARD")

#remove species with no or low abundance
#SEV<-subset(SEV,species!="NEXX"&species!="NOCR"
#            &species!="ONXX"&species!="OTVA"
#            &species!="MYGA"&species!="MUMU"&species!="MIME")

# not sure if other species should be removed?
# may want to remove unidentified to species (e.g., DIXX)
summary(SEV$species)

#create web/grid code
SEV$web_ID<-as.factor(paste(SEV$location,SEV$web,sep="_"))

#create a season column that lumps summer data with fall for early years
summary(SEV$season)
SEV$season.fs<-recode_factor(SEV$season,'2'="fall",'1'="spring",'3'="fall")
summary(SEV)

#create a count variable so that mammals can be summarized
SEV$count<-as.numeric(1)

#reshape so that species are columns
#for now these will be total tallies, we will get per trap night later
SEVmatrix1<-dcast(SEV,location+year+season.fs+season+web+web_ID~species,value.var="count",function(x) sum(x),fill=0)

#look at new dataframe (939 X 35)
summary(SEVmatrix1)

#normalize by number of trap nights
#create a separate vector of the number of trapping nights for each web
#SEV$YS<-as.factor(paste(SEV$year,SEV$season,sep="_"))
#library(tidyverse)
TN<-SEV %>% 
  dplyr::group_by(year,season,web_ID,web) %>% 
  dplyr::summarise(max_nights=max(night,na.rm=TRUE))

#merge TN into SEVmatrix
SEVmatrix2<-merge(TN,SEVmatrix1,by=c("year","season","web_ID","web"))
summary(SEVmatrix2)

#subset to species matrix only
SEVspp<-SEVmatrix2[,8:36]
head(SEVspp)

#and divide by trap nights and web area to get indiv/km2
SEV_TN<-SEVspp %>%
  mutate_all(list(~(. /(SEVmatrix2$max_nights*0.031416))))
                  
summary(SEV_TN)

#get rid of trap night count 
SEVmatrix3<-select(SEVmatrix2,-max_nights)

#put back together adjusted counts with variable names
mammals_async<-cbind(SEVmatrix3[,1:6],SEV_TN)

#look at data (939 X 35)
summary(mammals_async)

#sum total mammal abundance
mammals_async$total<-rowSums(mammals_async[,7:35])

# add richness, diversity, evenness
mammals_async$richness<-specnumber(mammals_async[,7:35])
mammals_async$shannonH<-diversity(mammals_async[,7:35])
mammals_async$evenness<-mammals_async$shannonH/log(mammals_async$richness,1)

#remove mammal webs that have no plant data or inconsistency in sampling over time
mammals_async <- mammals_async %>% filter(web_ID=="core_black_1"|web_ID=="core_black_4"|web_ID=="core_black_5"|
                                            web_ID=="core_creosote_1"|web_ID=="core_creosote_3"|web_ID=="core_creosote_5")
#rename ecosystems
summary(mammals_async$location)
mammals_async$ecosystem <- recode_factor(mammals_async$location,core_black="Desert grassland", 
                               core_creosote="Desert shrubland")

# rename seasons
mammals_async$season <- mammals_async$season.fs
summary(mammals_async$season)
mammals_async <- mammals_async %>% dplyr::select(-season.fs)

######## TOTAL MAMMAL ABUNDANCE STABILITY ----------------------------------------------------
totabund<-mammals_async %>% select(year,season,web_ID,web,ecosystem,total)
totabund_cast<-dcast(totabund, season+web_ID+web+ecosystem ~ year, sum, value.var="total",fill=-9999)

# fix problem with fill in dcast
totabund_cast[totabund_cast== "-9999"] <- "NA"
totabund_cast <- totabund_cast %>% mutate_if(is.character, as.numeric) #warning is not important
summary(totabund_cast)

# calculate totabund sd, mean, CV
#totabund_cast$sd89<-apply(totabund_cast[,5:39], 1, sd, na.rm=TRUE)
#totabund_cast$mu89<-apply(totabund_cast[,14:39], 1, mean, na.rm=TRUE)
#totabund_cast$CV89<-totabund_cast$sd89/totabund_cast$mu89
#totabund_cast$biomass_stability89<-1/(totabund_cast$CV89)
#totabund_cast$biomass_stability89.2<-1/((totabund_cast$CV89)^2)

# get metrics from 2009 - 2022 to match with years of consecutive trapping AND plant data ends 2022
totabund_cast$sd09<-apply(totabund_cast[,25:38], 1, sd, na.rm=TRUE)
totabund_cast$mu09<-apply(totabund_cast[,25:38], 1, mean, na.rm=TRUE)
totabund_cast$CV09<-totabund_cast$sd09/totabund_cast$mu09
totabund_cast$abund_stability09<-1/(totabund_cast$CV09)
totabund_cast$abund_stability09.2<-1/((totabund_cast$CV09)^2)

# possible number of years (to 2022) = 14
# so you can filter dataset later to only quads that have a minimum number of years recorded
totabund_cast$year_count89<-rowSums(!is.na(totabund_cast[,5:38]))
totabund_cast$year_count09<-rowSums(!is.na(totabund_cast[,25:38]))
summary(totabund_cast)


###### MAMMAL SPECIES STABILITY ----------------------------------------------------
# melt to make a row for each species in each quad, differs from original datafile because now includes zeros
mammals_melt<-melt(mammals_async,id=c("year","season","web_ID","web","location","ecosystem"),variable.name="species",value.name="indiv_per_km2")
summary(mammals_melt$indiv_per_km2)
summary(mammals_melt)

# now recast data to make each row a species X web X season combination, with columns for all years
# this step may be slow
mammals_cast<-dcast(mammals_melt, season+web_ID+web+location+ecosystem+species ~ year, sum, value.var = "indiv_per_km2",fill=-9999)
# fix problem with fill in dcast
mammals_cast[mammals_cast== "-9999"] <- "NA"
mammals_cast <- mammals_cast %>% mutate_if(is.character, as.numeric) #warning is not important
summary(mammals_cast)

#remove all species that were never observed in a web
mammals_cast$sum1989<-rowSums(mammals_cast[,7:41],na.rm=T)
summary(mammals_cast$sum1989)
mammals_stability <-mammals_cast %>%  filter(sum1989>0)

#count up number of years of observations for each species/web/season (row)
mammals_stability$years_count89<-rowSums(!is.na(mammals_stability[,7:41]))
mammals_stability$years_count09<-rowSums(!is.na(mammals_stability[,27:40]))
summary(mammals_stability)

# calculate mean long-term species biomass separately by season, also sd in species mass and CV
#mammals_stability$species_sd89<-apply(mammals_stability[,7:42], 1, sd, na.rm=TRUE)
#mammals_stability$species_mu89<-apply(mammals_stability[,7:42], 1, mean, na.rm=TRUE)
#mammals_stability$species_CV89<-mammals_stability$species_sd89/mammals_stability$species_mu89

mammals_stability$species_sd09<-apply(mammals_stability[,27:40], 1, sd, na.rm=TRUE)
mammals_stability$species_mu09<-apply(mammals_stability[,27:40], 1, mean, na.rm=TRUE)
mammals_stability$species_CV09<-mammals_stability$species_sd09/mammals_stability$species_mu09

#calculate stability as CV and 1/CV2
#mammals_stability$species_stability89<-1/(mammals_stability$species_CV89)
#mammals_stability$species_stability89.2<-1/((mammals_stability$species_CV89)^2)
mammals_stability$species_stability09<-1/(mammals_stability$species_CV09)
mammals_stability$species_stability09.2<-1/((mammals_stability$species_CV09)^2)


#########  PAIR WITH PLANT DATA ----------------------------------------------
# import plant  asynchrony data at web scale
# created in sev_plant_synchrony_webscale_mammals_09-22.R
sync<-read.csv("asynchrony_stability_allwebs_mammals.csv",stringsAsFactors = T)
summary(sync$web_ID)
summary(mammals_stability$web_ID)

# remove plant webs that have no mammal data matched by year
sync <- sync %>% filter(web_ID=="core_black_1_NA"|web_ID=="core_black_4_NA"|web_ID=="core_black_5_NA"|
                     web_ID=="core_creosote_1_NA"|web_ID=="core_creosote_3_NA"|web_ID=="core_creosote_5_NA")

summary(sync$web_ID)

# make web_ID match format from plant data
totabund_cast$web_ID<-as.factor(paste(totabund_cast$web_ID,"NA",sep="_"))
mammals_stability$web_ID<-as.factor(paste(mammals_stability$web_ID,"NA",sep="_"))
summary(totabund_cast$web_ID) 
summary(mammals_stability$web_ID)

# MERGE
mammal_species_async <- merge(mammals_stability, sync, by.x=c("web_ID","season","ecosystem","web"), by.y=c("web_ID","season","ecosystem","web"), all.x=T)
summary(mammal_species_async)
write.csv(mammal_species_async, "mammal_species_async_stability.csv")

mammal_total_async <- merge(totabund_cast, sync, by.x=c("web_ID","season","ecosystem","web"), by.y=c("web_ID","season","ecosystem","web"), all.x=T)
summary(mammal_total_async)
write.csv(mammal_total_async, "mammal_total_async_stability.csv")
