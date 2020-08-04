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
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Late Stuart"] <- "L.Stuart"
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Harrison-Widgeon"] <- "Harrison"
sw.data <- spsu #rename
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Early Stuart"] <- "E.Stuart"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Late Shuswap"] <- "L.Shuswap"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Late Stuart"] <- "L.Stuart"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Harrison-Widgeon"] <- "Harrison"
fw.data <- fw.major #rename

#### Brood table data
#Run code from sst_anomalies.R and exploratory_stage_1_analysis.Rmd for sst and brood table data 
#Add column with stock name by ID#


#### SST data
raw.clim <- read.csv("data/sst_yr_1_stock_anomalies.csv")
early.sst <- clim.wgt.avg(brood.table = brood_table,
                          env.data = raw.clim,
                          env.covar = "sst_anomaly",
                          type = "first_year",
                          out.covar = "early_sst") 
head(early.sst)
stock.ids2 <- brood_table[,1:2]
stock.ids <- stock.ids2[!duplicated(stock.ids2),]
early.sst <- merge(early.sst, stock.ids, by="Stock.ID")
names(early.sst) <- c("Stock.ID", "brood_year", "sst_anom", "Stock_Analysis")

# Create "agent" and "years" objects
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
mod_ic_mul_fw <- lmer(resid_value ~ prev + (1 | Stock_Analysis), ic_mul.resid.fw, REML=F)
summary(mod_ic_mul_fw)
mod_ic_mul_fw_null <- lmer(resid_value ~ (1 | Stock_Analysis), ic_mul.resid.fw, REML=F)
summary(mod_ic_mul_fw_null)
anova(mod_ic_mul_fw, mod_ic_mul_fw_null)

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

# stats te_mar
names(sw.data)
te_mar.data <- na.omit(sw.data[,c(5,22,79)])
te_mar.data$log.te_mar<- log10(te_mar.data$te_mar)
te_mar.data$log.te_mar[te_mar.data$log.te_mar == "-Inf"] <- 0

mod.te_mar.lat <- lmer(log.te_mar ~ Latitude + (1 | Year), te_mar.data, REML=F)
summary(mod.te_mar.lat)
mod.te_mar.lat.null <- lmer(log.te_mar ~ (1 | Year), te_mar.data, REML=F)
summary(mod.te_mar.lat.null)
anova(mod.te_mar.lat, mod.te_mar.lat.null)

#See maps in GitHub - te_mar detections centeres in northern SOG, DI, JS
#Possible migration conditions, density (transmission), or exposure to farmed salmon (if te_mar an issue)
#generally low prevalence agent

###############################################################################################
# LATITUDE AND TEMPERATURE CORRELATIONS

## pa_ther
sw.data.year.st <-sw.data %>% group_by(Year, Stock_Analysis) #create object to be summarized by year
all.pa_ther.sw =
  data.frame(
    sw.data.year.st %>% 
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
pa_ther.sst <- merge(all.pa_ther.sw, early.sst, by = c("Stock_Analysis", "brood_year"), all.x=TRUE)

## Latitude plots
ggplot(all.pa_ther.sw,aes(prev, Latitude, color=factor(brood_year))) +
  geom_point() +
  geom_linerange(aes(ymin = minLat, ymax = maxLat))
ggplot(sw.data,aes(Latitude, log10(pa_ther), color=factor(Year))) +
  geom_point() +
  geom_smooth(aes(Latitude, log10(pa_ther)), method = "lm", se=F, size=.2) 

## Temperature plots
ggplot(pa_ther.sst) +
  geom_point(aes(sst_anom, log10(mean_load), color=factor(Stock_Analysis))) +
  geom_smooth(data=pa_ther.sst, aes(sst_anom, log10(mean_load)), color="black", method = "lm", se=T, size=.2) +
  labs(x="SST anomaly", y="Log mean load", title=expression(italic("P. theridion")), color="Stock")

jpeg(filename='figs/Fig_pa_ther load corr w SST by year.jpg', 
     width=400, height=400, quality=75)
ggplot(pa_ther.sst) +
  geom_point(aes(sst_anom, log10(mean_load), color=factor(Year))) +
  geom_smooth(data=pa_ther.sst, aes(sst_anom, log10(mean_load)), color="black", method = "lm", se=T, size=.2) +
  labs(x="SST anomaly", y="Log mean load", title=expression(italic("P. theridion")), color="Sampling year")
dev.off()

jpeg(filename='figs/Fig_pa_ther prev corr w SST by year.jpg', 
     width=400, height=400, quality=75)
ggplot(pa_ther.sst) +
  geom_point(aes(sst_anom, prev, color=factor(Year))) +
  geom_smooth(data=pa_ther.sst, aes(sst_anom, prev), color="black", method = "lm", se=T, size=.2) +
  labs(x="SST anomaly", y="Prevalence", title=expression(italic("P. theridion")), color="Sampling year")
dev.off()


## te_mar
sw.data.year.st <-sw.data %>% group_by(Year, Stock_Analysis) #create object to be summarized by year
all.te_mar.sw =
  data.frame(
    sw.data.year.st %>% 
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
names(all.te_mar.sw) <- c("Year", "Stock_Analysis", "Latitude", "minLat", "maxLat", "N", "N+", "prev", "mean_load", "prevload") #rename columns

all.te_mar.sw$brood_year<-all.te_mar.sw$Year-2
te_mar.sst <- merge(all.te_mar.sw, early.sst, by = c("Stock_Analysis", "brood_year"), all.x=TRUE)

## Latitude plots
ggplot(all.te_mar.sw,aes(prev, Latitude, color=factor(Year))) +
  geom_point() +
  geom_linerange(aes(ymin = minLat, ymax = maxLat))+
  labs(x="Prevalence", y="Latitude", title=expression(italic("Tenacibaculum maritinum")), color="Sampling year")

ggplot(sw.data) +
  geom_point(aes(Latitude, log10(te_mar), color=factor(Year))) +
  geom_smooth(aes(Latitude, log10(te_mar)), method = "lm", se=F, size=.2, color="black") +
  labs(x="Log mean load", y="Latitude", title=expression(italic("Tenacibaculum maritinum")), color="Sampling year")

## Temperature plots
jpeg(filename='figs/Fig_te_mar load corr w SST by year.jpg', 
     width=400, height=400, quality=75)
ggplot(te_mar.sst) +
  geom_point(aes(sst_anom, log10(mean_load), color=factor(Year))) +
  geom_smooth(data=te_mar.sst, aes(sst_anom, log10(mean_load)), color="black", method = "lm", se=T, size=.2) +
  labs(x="SST anomaly", y="Log mean load", title=expression(italic("Tenacibaculum maritinum")), color="Sampling year")
dev.off()

jpeg(filename='figs/Fig_te_mar prev corr w SST by year.jpg', 
     width=400, height=400, quality=75)
ggplot(te_mar.sst) +
  geom_point(aes(sst_anom, prev, color=factor(Year))) +
  geom_smooth(data=te_mar.sst, aes(sst_anom, prev), color="black", method = "lm", se=T, size=.2) +
  labs(x="SST anomaly", y="Prevalence", title=expression(italic("Tenacibaculum maritinum")), color="Sampling year")
dev.off()


## ic_mul
sw.data.year.st <-sw.data %>% group_by(Year, Stock_Analysis) #create object to be summarized by year
sw.data.year <-sw.data %>% group_by(Year) #create object to be summarized by year
all.ic_mul.sw =
  data.frame(
    sw.data.year.st %>% 
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
names(all.ic_mul.sw) <- c("Year", "Stock_Analysis", "Latitude", "minLat", "maxLat", "N", "N+", "prev", "mean_load", "prevload") #rename columns
all.ic_mul.sw$brood_year<-all.ic_mul.sw$Year-2
ic_mul.sst <- merge(all.ic_mul.sw, early.sst, by = c("Stock_Analysis", "brood_year"), all.x=TRUE)

### Not stock specific
all.ic_mul.sw.global =
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
names(all.ic_mul.sw.global) <- c("Year", "Latitude", "minLat", "maxLat", "N", "N+", "prev", "mean_load", "prevload") #rename columns
all.ic_mul.sw.global$brood_year<-all.ic_mul.sw.global$Year-2
ic_mul.sst.global <- merge(all.ic_mul.sw.global, early.sst, by = "brood_year", all.x=TRUE)
ic_mul.sst.global = do.call(rbind,  lapply(split(early.sst, early.sst$Stock_Analysis), function(x)
  merge(x, all.ic_mul.sw.global, by = "brood_year", all = F)))

## Latitude plots
ggplot(all.ic_mul.sw,aes(prev, Latitude, color=factor(Year))) +
  geom_point() +
  geom_linerange(aes(ymin = minLat, ymax = maxLat))+
  labs(x="Prevalence", y="Latitude", title=expression(italic("Ichthyophthirius multifiliis")), color="Sampling year")
ggplot(sw.data) +
  geom_point(aes(Latitude, log10(ic_mul), color=factor(Year))) +
  geom_smooth(aes(Latitude, log10(ic_mul)), method = "lm", se=T, size=.2, color="black") +
  labs(y="Log mean load", x="Latitude", title=expression(italic("Ichthyophthirius multifiliis")), color="Sampling year")

ggplot(all.ic_mul.sw,aes(prev, Latitude, color=factor(Stock_Analysis))) +
  geom_point() +
  geom_linerange(aes(ymin = minLat, ymax = maxLat))+
  labs(x="Prevalence", y="Latitude", title=expression(italic("Ichthyophthirius multifiliis")), color="Stock")
ggplot(sw.data) +
  geom_point(aes(Latitude, log10(ic_mul), color=factor(Stock_Analysis))) +
  geom_smooth(aes(Latitude, log10(ic_mul)), method = "lm", se=F, size=.2, color="black") +
  labs(y="Log mean load", x="Latitude", title=expression(italic("Ichthyophthirius multifiliis")), color="Stock")

# Not stock specific
ggplot(all.ic_mul.sw.global,aes(prev, Latitude, color=factor(Year))) +
  geom_point() +
  geom_linerange(aes(ymin = minLat, ymax = maxLat))+
  labs(x="Prevalence", y="Latitude", title=expression(italic("Ichthyophthirius multifiliis")), color="Sampling year")


## Temperature plots
jpeg(filename='figs/Fig_ic_mul load corr w SST by year.jpg', 
     width=400, height=400, quality=75)
ggplot(ic_mul.sst) +
  geom_point(aes(sst_anom, log10(mean_load), color=factor(Year))) +
  geom_smooth(data=ic_mul.sst, aes(sst_anom, log10(mean_load)), color="black", method = "lm", se=T, size=.2) +
  labs(x="SST anomaly", y="Log mean load", title=expression(italic("Ichthyophthirius multifiliis")), color="Sampling year")
dev.off()

jpeg(filename='figs/Fig_ic_mul prev corr w SST by year.jpg', 
     width=400, height=400, quality=75)
ggplot(ic_mul.sst) +
  geom_point(aes(sst_anom, prev, color=factor(Year))) +
  geom_smooth(data=ic_mul.sst, aes(sst_anom, prev), color="black", method = "lm", se=T, size=.2) +
  labs(x="SST anomaly", y="Prevalence", title=expression(italic("Ichthyophthirius multifiliis")), color="Sampling year")
dev.off()

ggplot(ic_mul.sst.global) +
  geom_point(aes(sst_anom, prev, color=factor(Year))) +
  geom_smooth(data=ic_mul.sst.global, aes(sst_anom, prev), color="black", method = "lm", se=T, size=.2) +
  labs(x="SST anomaly", y="Prevalence", title=expression(italic("Ichthyophthirius multifiliis")), color="Sampling year")
dev.off()

ggplot(ic_mul.sst.global) +
  geom_point(aes(sst_anom, log10(mean_load), color=factor(Year))) +
  geom_smooth(data=ic_mul.sst.global, aes(sst_anom, log10(mean_load)), color="black", method = "lm", se=T, size=.2) +
  labs(x="SST anomaly", y="Log mean load", title=expression(italic("Ichthyophthirius multifiliis")), color="Sampling year")
dev.off()

