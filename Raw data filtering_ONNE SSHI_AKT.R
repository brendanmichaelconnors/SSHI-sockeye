## ONNE metadata filtering code for Andrew Bateman
## 200811 - Amy K Teffer

## Set your working directory to where you will source the files from (change the file path below)
## I keep my csv data files in a "data" folder within the project folder, so the code can call from there to load it
setwd("~/Documents.nosync/DFO PDF/Data/SSHI-sockeye")

#Bring in pathogen data without LOD and remove samples with no stock ID (may be updated data available!)
all <- read.csv("data/ONNE metadata no LOD_3.27.2019.csv",header=TRUE)
all2 <- droplevels(all[!(all$Stock_Region=="") ,])#remove any fish without stock assignment for now

#take out rare stock regions (northern, QCI, Transboundary)
temp2<-droplevels(all2[-which(all2$Stock_Area=="Northern") ,]) #Northern out
dim(temp2)
temp3<-droplevels(temp2[-which(temp2$Stock_Area=="QCI") ,]) #QCI out
dim(temp3)
major<-droplevels(temp3[-which(temp3$Stock_Area=="TransBoundary") ,]) #Transboundary out
dim(major)

## Divide up SW and FW collected samples
sw.major<-droplevels(major[(major$SWFW=="SW") ,])#SW only
head(sw.major)
fw.major<-droplevels(major[(major$SWFW=="FW") ,])#FW only 
head(fw.data)

## Temporally limit: Reduce sampling period to spring-summer only
spsu1<-droplevels(sw.major[-which(sw.major$SEASON1=="Overwinter") ,]) #remove winter 
dim(spsu1)
spsu2<-droplevels(spsu1[-which(spsu1$SEASON1=="Fall") ,]) # remove fall 
dim(spsu2) 
spsu3<-droplevels(spsu2[-which(spsu2$Year=="2018") ,]) #remove 2018 - no SR data yet
dim(spsu3)

## Geographically limit: Remove samples from WCVI, 2018 and high latitudes (>51.5 lat) 
spsu4<-droplevels(spsu3[-which(spsu3$Zone=="WCVI") ,]) #remove WCVI 
dim(spsu4)
spsu<-droplevels(spsu4[-which(spsu4$Latitude > 51.5) ,]) # remove high latitude samples
dim(spsu)

## Change stock names to align with SR data in both SW and FW data
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

## Examine output
head(sw.data)
head(fw.data)
