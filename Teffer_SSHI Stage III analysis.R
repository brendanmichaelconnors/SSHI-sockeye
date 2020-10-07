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

# Clean FW data
fw3<-droplevels(fw.major[-which(fw.major$Stock_Analysis=="VI") ,]) #remove VI 
dim(fw3)
fw4<-droplevels(fw3[-which(fw3$Stock_Analysis=="Central Coast") ,]) #remove CC
dim(fw4)
fw.data <- fw4 #rename

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
  group_by(Stock_Analysis, Year) %>%
  count(Year)

jpeg(filename='figs/Fig_Total fish sampled by stock per year_FW_yearY.jpg', 
     width=800, height=500, quality=75)
ggplot(data=samplesperstock.fw, aes(x=factor(Year), y=n, fill=Stock_Analysis))+
  geom_bar(stat="identity") +
  labs(fill="Stock") +
  coord_flip()+
  xlab("Year")+
  ylab("Total fish sampled")
dev.off()


## Freshwater totals by stock
samplesperstock.fw.st<-fw.data %>% 
  group_by(Stock_Analysis) %>%
  count(Stock_Analysis)
# Stocks sampled in FW
jpeg(filename='figs/Fig_Total fish sampled by Stock_Analysis_FW.jpg', 
     width=800, height=500, quality=75)
ggplot(data=samplesperstock.fw.st, aes(x=reorder(Stock_Analysis, n), y=n))+
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Stock") +
  ylab("Total fish sampled")
dev.off()

ggplot(data=fw.data, aes(x=Stock_Analysis, y=Year))+
  geom_point(stat="identity") +
  coord_flip() +
  xlab("Stock") +
  ylab("FW Longitude")

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
fw.data.year <-fw.data %>% group_by(Year, Stock_Analysis) #create object to be summarized by year

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
names(all.ic_mul.fw) <- c("Year","Stock_Analysis","N","N+","prev","mean_load","prevload") #rename columns
all.ic_mul.fw$brood_year <- all.ic_mul.fw$Year-2
ic_mul.resid.fw <- merge(trnc_resid_srr, all.ic_mul.fw, by = c("brood_year", "Stock_Analysis"))

jpeg(filename='figs/Fig_ic_mul FW prev corr w SR resid.jpg', 
     width=480, height=500, quality=75)
ggplot(ic_mul.resid.fw, aes(prev, resid_value, color=Stock_Analysis)) +
  geom_point(aes(color=Stock_Analysis)) +
  geom_smooth(aes(prev, resid_value), method = "lm", se=F, size=.2) +
  labs(y = "Stock-recruitment residuals",x = "Freshwater prevalence per year", 
       title=expression(paste(italic("I. multifiliis"))), color="Stock")
dev.off()

# Stats for ic_mul
## STAN model for ic_mul FW
# PREVALENCE
fw.ic_mul.FW <- stan_lmer(resid_value ~ 0 + prev + (prev|Stock_Analysis) +(1|brood_year), 
          data = ic_mul.resid.fw,
          adapt_delta=0.95,
          REML = F)
summary(fw.ic_mul.FW)

ind_coef.mean <- summary(fw.ic_mul.FW, 
        pars = c("prev"),
        probs = c(0.025,0.25,0.5,0.75, 0.975),
        digits = 2)
fw.mod.post.ic_mul.mean <- as.matrix(ind_coef.mean[1, c(4:8)])
fw.mod.post.ic_mul.mean <- data.frame(t(fw.mod.post.ic_mul.mean))

ind_coef <- summary(fw.ic_mul.FW, 
        regex_pars = c("b\\[\\prev Stock_Analysis\\:"),
        probs = c(0.025,0.25,0.5,0.75, 0.975),
        digits = 2)
fw.mod.post.ic_mul <- data.frame(ind_coef[c(1:15),c(4:8)])
fw.ic_mul.slp <- rbind(fw.mod.post.ic_mul.mean, fw.mod.post.ic_mul)
fw.ic_mul.slp <- rownames_to_column(fw.ic_mul.slp)
ggplot(fw.ic_mul.slp) +
  geom_hline(yintercept = 0, linetype = "dashed", col="blue")+
  geom_linerange(aes(x = reorder(rowname, -X50.), ymax = X75., ymin = X25.), size=1.5, col="black") +
  geom_linerange(aes(x = rowname, ymax = X97.5., ymin = X2.5.), col="black") +
  geom_point(aes(x = rowname, y = X50.), size = 3) +
  #ylim(-0.8,0.8)+
  coord_flip()

temp2 <- data.frame(fw.ic_mul.FW)
temp3 <- temp2[,grepl("prev",names(temp2))] #include only columns with "prev" in name
temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
temp5 <- temp4[2001:4000,] #remove warm up iterations
temp6 <- data.frame(temp5[,1] + temp5[,2:16]) #add the stock-specific draws to global column
temp7 <- cbind(temp5[,1],temp6) #bind with global column
para_name2 <- colnames(temp5) #create an object with column names
colnames(temp7) <- colnames(temp5) #assign names to columns
temp8 <- as.matrix(apply(temp7, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975)))
temp9 <- t(temp8)
temp10<- as.matrix(temp9)
colnames(temp10) <- c("2.5","25","50","75","97.5") #assign names to columns
write.csv(temp10, file="data/Stock specific slopes_prev ic_mul_FW_by stock prev.csv")
stspslp.ic_mulfw <- read.csv("data/Stock specific slopes_prev ic_mul_FW_by stock prev.csv")
stspslp.ic_mulfw$stock <- substr(stspslp.ic_mulfw$X, 23, 33)
stspslp.ic_mulfw$stock <- substr(stspslp.ic_mulfw$stock, 1, nchar(stspslp.ic_mulfw$stock)-1)
stspslp.ic_mulfw$stock <- sub("^$", "Global", stspslp.ic_mulfw$stock)
colnames(stspslp.ic_mulfw) <- c("X","X2.5","X25","X50","X75","X97.5","stock")

## Plot
jpeg(filename='figs/Fig_SSHI ONNE_stock sp slope_ic_mul_FW.jpg', 
     width=480, height=500, quality=75)
ggplot(stspslp.ic_mulfw) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stspslp.ic_mulfw[stspslp.ic_mulfw$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stspslp.ic_mulfw[stspslp.ic_mulfw$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stspslp.ic_mulfw[stspslp.ic_mulfw$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size") +
  coord_flip()
dev.off()

# Proportion <0
model<-as.matrix(fw.ic_mul.FW)
model2<-model[2001:4000,]
param.prop0.ic_mul.fw <- (colSums(model2 < 0))/2000 #proportion <0
minfw <- colSums(model2 > 0) #total above zero


# INTENSITY - will not run without divergent transitions
fw.ic_mul.FW.load <- stan_lmer(resid_value ~ 0 + mean_load + (mean_load|Stock_Analysis) +(1|brood_year), 
                          data = ic_mul.resid.fw,
                          adapt_delta=0.99,
                          REML = F)
summary(fw.ic_mul.FW.load)

ind_coef.mean.load <- summary(fw.ic_mul.FW.load, 
                         pars = c("mean_load"),
                         probs = c(0.025,0.25,0.5,0.75, 0.975),
                         digits = 2)
fw.mod.post.ic_mul.mean.load <- as.matrix(ind_coef.mean.load[1, c(4:8)])
fw.mod.post.ic_mul.mean.load <- data.frame(t(fw.mod.post.ic_mul.mean.load))

ind_coef <- summary(fw.ic_mul.FW.load, 
                    regex_pars = c("b\\[\\mean_load Stock_Analysis\\:"),
                    probs = c(0.025,0.25,0.5,0.75, 0.975),
                    digits = 2)
fw.mod.post.ic_mul.load <- data.frame(ind_coef[c(1:18),c(4:8)])
fw.ic_mul.slp.load <- rbind(fw.mod.post.ic_mul.mean, fw.mod.post.ic_mul)
fw.ic_mul.slp.load <- rownames_to_column(fw.ic_mul.slp.load)
ggplot(fw.ic_mul.slp.load) +
  geom_hline(yintercept = 0, linetype = "dashed", col="blue")+
  geom_linerange(aes(x = reorder(rowname, -X50.), ymax = X75., ymin = X25.), size=1.5, col="black") +
  geom_linerange(aes(x = rowname, ymax = X97.5., ymin = X2.5.), col="black") +
  geom_point(aes(x = rowname, y = X50.), size = 3) +
  #ylim(-0.8,0.8)+
  coord_flip()

temp2 <- data.frame(fw.ic_mul.FW.load)
temp3 <- temp2[,grepl("mean_load",names(temp2))] #include only columns with "mean_load" in name
temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
temp5 <- temp4[2001:4000,] #remove warm up iterations
temp6 <- data.frame(temp5[,1] + temp5[,2:19]) #add the stock-specific draws to global column
temp7 <- cbind(temp5[,1],temp6) #bind with global column
para_name2 <- colnames(temp5) #create an object with column names
colnames(temp7) <- colnames(temp5) #assign names to columns
temp8 <- as.matrix(apply(temp7, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975)))
temp9 <- t(temp8)
temp10<- as.matrix(temp9)
colnames(temp10) <- c("2.5","25","50","75","97.5") #assign names to columns
write.csv(temp10, file="data/Stock specific slopes_load ic_mul_FW.csv")
stspslp.ic_mulfw.load <- read.csv("data/Stock specific slopes_load ic_mul_FW.csv")
stspslp.ic_mulfw.load$stock <- substr(stspslp.ic_mulfw.load$X, 28, 38)
stspslp.ic_mulfw.load$stock <- substr(stspslp.ic_mulfw.load$stock, 1, nchar(stspslp.ic_mulfw.load$stock)-1)
stspslp.ic_mulfw.load$stock <- sub("^$", "Global", stspslp.ic_mulfw.load$stock)
colnames(stspslp.ic_mulfw.load) <- c("X","X2.5","X25","X50","X75","X97.5","stock")

## Plot
jpeg(filename='figs/Fig_SSHI ONNE_stock sp slope_ic_mul_FW_LOAD.jpg', 
     width=480, height=500, quality=75)
ggplot(stspslp.ic_mulfw.load) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stspslp.ic_mulfw.load[stspslp.ic_mulfw.load$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stspslp.ic_mulfw.load[stspslp.ic_mulfw.load$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stspslp.ic_mulfw.load[stspslp.ic_mulfw.load$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size") +
  coord_flip()
dev.off()

# Proportion <0
model<-as.matrix(fw.ic_mul.FW.load)
model2<-model[2001:4000,]
param.prop0.ic_mul.fw <- (colSums(model2 < 0))/2000 #proportion <0
minfw <- colSums(model2 > 0) #total above zero






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
ggplot(all.pa_ther.sw,aes(Latitude, prev)) +
  geom_point() +
  geom_smooth(aes(Latitude, prev), method = "loess", se=F, size=.2) 
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

pa_ther.mod.sst.load <- lmer(log10(mean_load) ~ sst_anom + sst_anom|Year, data = pa_ther.sst, REML=T)
summary(pa_ther.mod.sst.load)
pa_ther.mod.sst <- lmer(prev ~ sst_anom + sst_anom|Year, data = pa_ther.sst, REML=F)
summary(pa_ther.mod.sst)

library(mgcv)
install.packages(tidymv)
gam.pa_ther.sst <- gam(prev ~ sst_anom, data = pa_ther.sst)
summary(gam.pa_ther.sst)
plot_smooths(gam.pa_ther.sst)
  

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


### Examining loads relative to distance to farms and Fraser River mouth
ONNE_with_dists <- readRDS("~/Documents.nosync/DFO PDF/Data/SSHI-sockeye/ONNE_with_dists.rds")
names(ONNE_with_dists)

### te_mar
ggplot(ONNE_with_dists)+
  geom_point(aes(nearest_AC_sway, log10(te_mar)))+
  geom_smooth(data=ONNE_with_dists, aes(nearest_AC_sway, log10(te_mar)), method="lm")
ggplot(ONNE_with_dists)+
  geom_point(aes(fm_dist, log10(te_mar)))+
  geom_smooth(data=ONNE_with_dists, aes(fm_dist, log10(te_mar)), method="lm")


### pa_ther
ggplot(ONNE_with_dists)+
  geom_point(aes(nearest_AC_sway, log10(pa_ther)))+
  geom_smooth(data=ONNE_with_dists, aes(nearest_AC_sway, log10(pa_ther)), method="lm")
ggplot(ONNE_with_dists)+
  geom_point(aes(fm_dist, log10(pa_ther)))+
  geom_smooth(data=ONNE_with_dists, aes(fm_dist, log10(pa_ther)), method="lm")

### ic_mul
ggplot(ONNE_with_dists)+
  geom_point(aes(nearest_AC_sway, log10(ic_mul)))+
  geom_smooth(data=ONNE_with_dists, aes(nearest_AC_sway, log10(ic_mul)), method="lm")
ggplot(ONNE_with_dists)+
  geom_point(aes(fm_dist, log10(ic_mul)))+
  geom_smooth(data=ONNE_with_dists, aes(fm_dist, log10(ic_mul)), method="lm")


##Define breaks
ONNE_with_dists$te_mar.ct <- as.numeric(ifelse(ONNE_with_dists$te_mar > 0, 1, 0))
ONNE_with_dists$brks <- as.numeric(cut_number(ONNE_with_dists$fm_dist,30))
ONNE_with_dists <- data.frame(ONNE_with_dists)
ONNE_with_dists.brks <-ONNE_with_dists %>% group_by(brks, te_mar.ct) #create object to be summarized by brks

ggplot(data=ONNE_with_dists) +
  geom_bar(aes(x=brks), col="black", position=position_stack()) +
  geom_bar(stat="identity", aes(x=brks, y=te_mar.ct, fill=factor(Year)), col="black", position=position_stack()) 

te_mar.dist1 <- aggregate(ONNE_with_dists$te_mar.ct,
                         by=list(dist.bin=ONNE_with_dists$brks), 
                         FUN=sum, na.rm=TRUE, na.action=NULL)
te_mar.dist2 <- aggregate(ONNE_with_dists$te_mar.ct,
                          by=list(dist.bin=ONNE_with_dists$brks), 
                          count, na.rm=TRUE, na.action=NULL)

library(dplyr)
temp <-
  ONNE_with_dists %>%
  group_by(brks, te_mar.ct) %>%
  count(brks)

temp =
  data.frame(
    ONNE_with_dists.brks %>% 
      summarise(
        length(te_mar.ct>0, na.rm=TRUE),
        length(brks, na.rm=TRUE)
  )
)
names(temp) <- c("brks", "te_mar.pos", "n") #rename columns



ggplot(data=ONNE_with_dists, aes(x = fm_dist, fill=brks)) +
  geom_dotplot(stackgroups = TRUE, binwidth = 4000, method = "histodot") +
  scale_fill_manual(values=colorRampPalette(c("white", "red"))( length(ONNE_with_dists$brks)))
head(ONNE_with_dists)


### Weight-length
