# SSHI ANALYSIS SOCKEYE STAGE III
# A.K. Teffer
#### Sockeye salmon productivity versus infection profiles

#### Load packages and set directory
setwd("~/Documents.nosync/DFO PDF/Data/SSHI-sockeye")

library(lme4)
library(rstanarm) # https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html
library(ggplot2)
library(plotrix)
library(tidyverse)
library(gridExtra)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(shinystan)
library(data.table)
library(base)


#### Read in data, clean, standardize - SW metric averaged across all stocks each year per agent  
inf_agt_resid_data_gl <- read.csv("data/global_ONNE_productivity_infection_analysis.csv")
head(inf_agt_resid_data_gl)

# Data cleaning
inf_agt_resid_data_gl$mean_load_all <- inf_agt_resid_data_gl$mean_load #add col for mean_load_all
inf_agt_resid_data_gl$mean_load_all[is.na(inf_agt_resid_data_gl$mean_load_all)] <- 0 #replace NA with 0

# Standardize and incorporate into df
inf_std <- plyr::ddply(inf_agt_resid_data_gl, c("agent"),function(x) {
  scaled_prev <- scale(x$prev)
  scaled_load <- scale(x$mean_load_all)
  xx <- data.frame(scaled_prev, scaled_load)
})
inf_agt_resid_data_gl$prev_std <- inf_std[,2]
inf_agt_resid_data_gl$load_std <- inf_std[,3]

# Add Stock column to update later
inf_agt_resid_data_gl$Stock <- inf_agt_resid_data_gl$Stock_Analysis
head(inf_agt_resid_data_gl)

# Create "agent" object
agents <- unique(inf_agt_resid_data_gl$agent)

# Plot total detections of each agent by year - note variable prevalence across agents and years
samplesperagent.sw<-inf_agt_resid_data_gl %>% 
  group_by(agent, Year) %>%
  summarise(N. = sum(N.))
ggplot(data=samplesperagent.sw, aes(x=reorder(agent, N.), y=N., fill=factor(Year)))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("Stock")+
  ylab("sample totals")

# Plot raw data by: 
## Prevalence

ggplot(inf_agt_resid_data_gl,aes(prev, resid_value, color=Stock, shape=factor(Year)))+
  geom_smooth(aes(prev, resid_value, group=Stock), method = "lm", se=F, size=.2)+
  geom_point()+
  scale_shape_manual(values=1:nlevels(factor(inf_agt_resid_data_gl$Year))) +
  facet_wrap(~ agent,nrow=5, scales = "free")+
  xlab("prevalence")+
  ylab("residual")+
  theme_bw()

## Load
ggplot(inf_agt_resid_data_gl,aes(log10(mean_load), resid_value, color=Stock, shape=factor(Year)))+
  geom_smooth(aes(log10(mean_load), resid_value, group=Stock), method = "lm", se=F, size=.2)+
  geom_point()+
  scale_shape_manual(values=1:nlevels(factor(inf_agt_resid_data_gl$Year))) +
  facet_wrap(~ agent,nrow=5, scales="free")+
  xlab("log10 load")+
  ylab("residual")+
  theme_bw()


# Investigate variability of agents by Latitude
#CODE IN PROCESS#

## Plot raw agent data by latitude
### In spsu, calculate count/sampled per year
ic_mul.spsu<-spsu[spsu]
spsu.year <-spsu %>% group_by(Year) #create object to be summarized by year

##ic_mul
all.ic_mul.sw =
  data.frame(
    spsu.year %>% 
      summarise(
        mean(Latitude[Latitude!=0], na.rm=TRUE),
        min(Latitude[Latitude!=0], na.rm=TRUE),
        max(Latitude[Latitude!=0], na.rm=TRUE),
        length(which(ic_mul!="NA")), #samples
        length(which(ic_mul>0)), #positive detections
        length(which(ic_mul>0))/length(which(!is.na(ic_mul))),  #calculates prevalence
        mean(ic_mul[ic_mul!=0], na.rm=TRUE),
        (length(which(ic_mul>0)) / length(which(!is.na(ic_mul)))) * mean(ic_mul[ic_mul!=0], na.rm=TRUE)
      )
  )
names(all.ic_mul.sw) <- c("Year", "Latitude", "minLat", "maxLat", "N", "N+", "prev", "mean_load", "prevload") #rename columns
all.ic_mul.sw$brood_year<-all.ic_mul.sw$Year-2

ggplot(all.ic_mul.sw,aes(prev, Latitude, color=factor(brood_year))) +
  geom_point() +
  geom_linerange(aes(ymin = minLat, ymax = maxLat)) +
  coord_flip()
ggplot(all.ic_mul.sw,aes(Latitude, brood_year)) +
  geom_point()

ggplot(spsu,aes(Latitude, log10(ic_mul), shape=Zone, color=factor(Year))) +
  geom_point() +
  geom_smooth(aes(Latitude, log10(ic_mul)), method = "lm", se=F, size=.2) 


