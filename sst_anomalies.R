########################################################################################
# sst_anomalies.R
# download and process SST data and calculate early marine anomaly 
# largely compliments of M. Malick (NOAA) see: https://michaelmalick.com/
#########################################################################################

## download and process SST data UNCOMMENT IF YOU NEED TO DOWNLOAD RAW SST
#ersst::sst_download(years = 1950:2018,
#                    months = 1:12,
#                    save.dir = "./data/sst_raw/",
#                    version = 5)
#
#sst.raw.full <- ersst::sst_load(years = 1950:2018,
#                                months = 1:12,
#                                read.dir = "./data/sst_raw/",
#                                version = 5)
#
#sst.raw.np <- ersst::sst_subset_space(sst.raw.full,
#                                      lat.min = 36,
#                                      lat.max = 80,
#                                      lon.min = 170,
#                                      lon.max = 250)
#
#sst.raw.df <- ersst::sst_dataframe(sst.raw.np)
#
#write.csv(sst.raw.df, "./data/sst_raw.csv", row.names = FALSE)

sst.raw <- read.csv("./data/sst_raw.csv")
head(sst.raw)
tail(sst.raw)
sapply(sst.raw, class)
summary(sst.raw)

#------------------------------------------------------------------------------#
#   Calculate SST anomalies and average across specified period and region 
#------------------------------------------------------------------------------#

## Calculate SST anomalies
sst.anom <- sst.anomaly(sst.raw, ref.years = 1950:2010)
head(sst.anom)
tail(sst.anom)
summary(sst.anom)
sapply(sst.anom, class)

## Convert longitude to match salmon data
sst.anom$lon2 <- ifelse(sst.anom$lon > 180, sst.anom$lon - 360, sst.anom$lon)

write.csv(sst.anom, "data/sst_raw_anomalies.csv", row.names=F)

## Read in population summary table (includes many "extra" stocks)
summary_table <- read.csv("data/master_stock_info.csv", header=T)
head(summary_table)


## Calculate average SST anomaly within area where stock spends first few months of marine life 
sst_yr_1_stock_anomalies <- sst.averager(summary_table, sst.anom, distance = 400)
colnames(sst_yr_1_stock_anomalies) <- c("Year","sst_raw","sst_anomaly","Stock.ID")

write.csv(sst_yr_1_stock_anomalies, "data/sst_yr_1_stock_anomalies.csv", row.names=F)
