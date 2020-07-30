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


#Bring in pathogen data without LOD and clean
all <- read.csv("data/ONNE metadata no LOD_3.27.2019.csv",header=TRUE)
all2 <- droplevels(all[!(all$Stock_Region=="") ,])#remove any fish without stock assignment for now
#take out rare stock regions and no stock ID
temp2<-droplevels(all2[-which(all2$Stock_Area=="Northern") ,]) #Northern out
dim(temp2)
temp3<-droplevels(temp2[-which(temp2$Stock_Area=="QCI") ,]) #QCI out
dim(temp3)
major<-droplevels(temp3[-which(temp3$Stock_Area=="TransBoundary") ,]) #Transboundary out
dim(major)
#only major stock regions sampled in SW
sw.major<-droplevels(major[(major$SWFW=="SW") ,])#SW only major
str(sw.major)
fw.major<-droplevels(major[(major$SWFW=="FW") ,])#FW only major
str(fw.data)
#Reduce temporal period to spring-summer, remove samples from WCVI, 2018 and high latitudes (>51.5 lat) 
spsu1<-droplevels(sw.major[-which(sw.major$SEASON1=="Overwinter") ,]) #remove winter 
dim(spsu1)
spsu2<-droplevels(spsu1[-which(spsu1$SEASON1=="Fall") ,]) # remove fall 
dim(spsu2) ##spsu is only spring and summer
spsu3<-droplevels(spsu2[-which(spsu2$Zone=="WCVI") ,]) #remove WCVI 
dim(spsu3)
spsu4<-droplevels(spsu3[-which(spsu3$Year=="2018") ,]) #remove 2018 
dim(spsu4)
spsu<-droplevels(spsu4[-which(spsu4$Latitude > 51.5) ,]) # remove high latitude samples
dim(spsu)
#Change stock names
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Early Stuart"] <- "E.Stuart"
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Late Shuswap"] <- "L.Shuswap"
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Harrison-Widgeon"] <- "Harrison"
sw.data <- spsu #rename
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Early Stuart"] <- "E.Stuart"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Late Shuswap"] <- "L.Shuswap"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Harrison-Widgeon"] <- "Harrison"
fw.data <- fw.major #rename

sst <- read.csv("data/SST anomalies ONNE by Stock.csv")
names(sst) <- c("Year", "Stock_Analysis", "sst_anom")

# Create "agent" object
agents <- unique(inf_agt_resid_data_gl$agent)
years <- unique(inf_agt_resid_data_gl$brood_year)

#Bring in truncated resdiduals
trnc_resid<-read.csv("data/survival_indices_truncated.csv", head=TRUE)
str(trnc_resid)
trnc_resid$Year <- trnc_resid$brood_year+2
names(trnc_resid) <- c("orderID", "Stock_Analysis", "brood_year", "metric", "resid_value", "Year") #rename columns
##create object with just SRR_resid metric to align with infection data
trnc_resid_srr <- trnc_resid[trnc_resid$metric=="SR_resid",]

# Plot sampled fish per stock by year
samplesperstock.sw<-sw.data %>% 
  group_by(Stock, Year) %>%
  count(Year)

jpeg(filename='figs/Fig_Total fish sampled by stock per year_SW.jpg', 
     width=480, height=800, quality=75)
ggplot(data=samplesperstock.sw, aes(x=reorder(Stock, n), y=n, fill=factor(Year)))+
  geom_bar(stat="identity") +
  labs(fill="Sampling year") +
  coord_flip()+
  xlab("Stock")+
  ylab("Total fish sampled")
dev.off()

jpeg(filename='figs/Fig_Total fish sampled by stock per year_SW_yearY.jpg', 
     width=800, height=500, quality=75)
ggplot(data=samplesperstock.sw, aes(x=factor(Year), y=n, fill=Stock))+
  geom_bar(stat="identity") +
  labs(fill="Stock") +
  coord_flip()+
  xlab("Year")+
  ylab("Total fish sampled")
dev.off()

## Freshwater totals by stock and year
samplesperstock.fw<-fw.data %>% 
  group_by(Stock, Year) %>%
  count(Year)

jpeg(filename='figs/Fig_Total fish sampled by stock per year_FW_yearY.jpg', 
     width=800, height=500, quality=75)
ggplot(data=samplesperstock.fw, aes(x=factor(Year), y=n, fill=Stock))+
  geom_bar(stat="identity") +
  labs(fill="Stock") +
  coord_flip()+
  xlab("Year")+
  ylab("Total fish sampled")
dev.off()

# Investigate variability of agents by Latitude
#CODE IN PROCESS#

## Plot raw agent data by latitude
### In sw.data, calculate count/sampled per year

sw.data.year <-sw.data %>% group_by(Year) #create object to be summarized by year

##ic_mul
all.ic_mul.sw =
  data.frame(
    sw.data.year %>% 
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
  geom_linerange(aes(ymin = minLat, ymax = maxLat))

ggplot(sw.data,aes(Latitude, log10(ic_mul), shape=Zone, color=factor(Year))) +
  geom_point() +
  geom_smooth(aes(Latitude, log10(ic_mul)), method = "lm", se=F, size=.2) 


## ic_mul FW influence
fw.data.year <-fw.data %>% group_by(Year) #create object to be summarized by year

all.ic_mul.fw =
  data.frame(
    fw.data.year %>% 
      summarise(
        length(which(ic_mul!="NA")), #samples
        length(which(ic_mul>0)), #positive detections
        length(which(ic_mul>0))/length(which(!is.na(ic_mul))),  #calculates prevalence
        mean(ic_mul[ic_mul!=0], na.rm=TRUE),
        (length(which(ic_mul>0)) / length(which(!is.na(ic_mul)))) * mean(ic_mul[ic_mul!=0], na.rm=TRUE)
      )
  )
names(all.ic_mul.fw) <- c("Year", "N", "N+", "prev", "mean_load", "prevload") #rename columns
all.ic_mul.fw$brood_year <- all.ic_mul.fw$Year-2
ic_mul.resid.fw <- merge(trnc_resid_srr, all.ic_mul.fw, by = c("brood_year", "Year"))

jpeg(filename='figs/Fig_ic_mul FW prev corr w SR resid.jpg', 
     width=480, height=500, quality=75)
ggplot(ic_mul.resid.fw, aes(prev, resid_value, color=Stock_Analysis)) +
  geom_point(aes(color=Stock_Analysis)) +
  geom_smooth(aes(prev, resid_value), method = "lm", se=F, size=.2) +
  labs(y = "Stock-recruitment residuals",x = "Freshwater prevalence per year", 
       title=expression(paste(italic("I. multifiliis"))), color="Stock")
dev.off()

# Stats for ic_mul
mod_ic_mul_fw <- lmer(resid_value ~ prev + (1 | Stock_Analysis), ic_mul.resid.fw)
summary(mod_ic_mul_fw)
mod_ic_mul_fw_null <- lmer(resid_value ~ (1 | Stock_Analysis), ic_mul.resid.fw)
summary(mod_ic_mul_fw_null)

## te_mar
all.te_mar.sw =
  data.frame(
    sw.data.year %>% 
      summarise(
        mean(Latitude[Latitude!=0], na.rm=TRUE),
        min(Latitude[Latitude!=0], na.rm=TRUE),
        max(Latitude[Latitude!=0], na.rm=TRUE),
        length(which(te_mar!="NA")), #samples
        length(which(te_mar>0)), #positive detections
        length(which(te_mar>0))/length(which(!is.na(te_mar))),  #calculates prevalence
        mean(te_mar[te_mar!=0], na.rm=TRUE),
        (length(which(te_mar>0)) / length(which(!is.na(te_mar)))) * mean(te_mar[te_mar!=0], na.rm=TRUE)
      )
  )
names(all.te_mar.sw) <- c("Year", "Latitude", "minLat", "maxLat", "N", "N+", "prev", "mean_load", "prevload") #rename columns
all.te_mar.sw$brood_year<-all.te_mar.sw$Year-2

ggplot(all.te_mar.sw,aes(prev, Latitude, color=factor(brood_year))) +
  geom_point() +
  geom_linerange(aes(ymin = minLat, ymax = maxLat))
ggplot(sw.data,aes(Latitude, log10(te_mar), shape=Zone, color=factor(Year))) +
  geom_point() +
  geom_smooth(aes(Latitude, log10(te_mar)), method = "lm", se=F, size=.2) 

ggplot(sw.data, aes(Latitude, log10(te_mar))) +
  geom_point() 
#See maps in GitHub - te_mar detections centeres in northern SOG, DI, JS
#Possible migration conditions, density (transmission), or exposure to farmed salmon (if te_mar an issue)
#generally low prevalence agent


##pa_ther
sw.data.year.st <-sw.data %>% group_by(Year, Stock_Analysis) #create object to be summarized by year

all.pa_ther.sw =
  data.frame(
    sw.data.year %>% 
      summarise(
        mean(Latitude[Latitude!=0], na.rm=TRUE),
        min(Latitude[Latitude!=0], na.rm=TRUE),
        max(Latitude[Latitude!=0], na.rm=TRUE),
        length(which(pa_ther!="NA")), #samples
        length(which(pa_ther>0)), #positive detections
        length(which(pa_ther>0))/length(which(!is.na(pa_ther))),  #calculates prevalence
        mean(pa_ther[pa_ther!=0], na.rm=TRUE),
        (length(which(pa_ther>0)) / length(which(!is.na(pa_ther)))) * mean(pa_ther[pa_ther!=0], na.rm=TRUE)
      )
  )
names(all.pa_ther.sw) <- c("Year", "Stock_Analysis", "Latitude", "minLat", "maxLat", "N", "N+", "prev", "mean_load", "prevload") #rename columns
all.pa_ther.sw$brood_year<-all.pa_ther.sw$Year-2

ggplot(all.pa_ther.sw,aes(prev, Latitude, color=factor(brood_year))) +
  geom_point() +
  geom_linerange(aes(ymin = minLat, ymax = maxLat))

ggplot(sw.data,aes(Latitude, log10(pa_ther), shape=Zone, color=factor(Year))) +
  geom_point() +
  geom_smooth(aes(Latitude, log10(pa_ther)), method = "lm", se=F, size=.2) 

# examine temperature correlation
pa_ther.sst <- merge(all.pa_ther.sw, sst, by = c("Stock_Analysis", "Year"))
ggplot(pa_ther.sst,aes(sst_anom, log10(mean_load), color=factor(brood_year))) +
  geom_point() 
ggplot(pa_ther.sst,aes(sst_anom, log10(mean_load), color=Stock_Analysis)) +
  geom_point() 

ggplot(pa_ther.sst,aes(sst_anom, prev, color=factor(brood_year))) +
  geom_point() 

jpeg(filename='figs/Fig_pa_ther prev corr w SST by stock.jpg', 
     width=480, height=600, quality=75)
ggplot(pa_ther.sst,aes(sst_anom, prev, color=Stock_Analysis)) +
  geom_point() +
  geom_smooth(aes(sst_anom, prev), method = "lm", se=F, size=.2) 
dev.off()

