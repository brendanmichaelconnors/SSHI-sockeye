########################################################################################
# functions.R
# most functions compliments of M. Malick (NOAA) see: https://michaelmalick.com/
# ########################################################################################


## consec_years --------------------------------------------
consec_years <- function(x) {
    ## Find longest consective years
    ##
    ## This function takes as input a vector of years with NA values and finds
    ## the longest set of consecutive years without an NA. The function returns
    ## the vector of consecutive years.
    ##
    ## x = vector of years

    n <- ifelse(is.na(x), 0, x)
    brk <- c(0, which(diff(n) != 1), length(n))
    vec <- lapply(seq(length(brk) - 1),
                  function(i) n[(brk[i] + 1):brk[i+1]])
    consec <- vec[[which.max(lapply(vec, length))]]
    ret <- n[n %in% consec]

    return(ret)
}

if(FALSE) {

    n <- c(NA, 1950:1955)
    consec_years(n)

    n <- c(1949, NA, 1951:1955)
    consec_years(n)

    n <- c(NA, 1951:1955)
    consec_years(n)

    n <- c(1951:1955, NA)
    consec_years(n)

    n <- c(NA, 1950:1955, NA, 1957:1965, NA, 1967)
    consec_years(n)

    n <- c(1950:1955, 1957:1965, NA, 1967:1970)
    consec_years(n)

}



## load_rdata ----------------------------------------------
load_rdata <- function(path = "./output/", verbose = TRUE) {
    ## Load saved .RData output
    ##
    ## This function loads all *.RData files found in the input directory.
    ##
    ## path = path to directory to read output
    ## verbose = logical, print status

    if(!dir.exists(path))
        stop("Input path doesn't exist")

    fls <- list.files(path = path)

    for(i in fls) {

        if(verbose)
            cat("Reading file", which(fls == i), "of",
                length(fls), "...", "\n")

        load(paste(path, i, sep = ""), envir = .GlobalEnv)
    }

    if(verbose) {
        if(length(fls) == 0)
            cat("No files to read", "\n")
        else
            cat("Done!", "\n")
    }
}



## clim.wgt.avg --------------------------------------------
clim.wgt.avg <- function(brood.table,
                         env.data,
                         env.covar,
                         type,
                         out.covar = "index") {
    ## Calculate climate index weighted by age-class
    ##
    ## Create data frame with estimates of climate weighted by the proportion of
    ## recruits for each age-class. Indices can be created for either the first
    ## ocean year or second ocean year.
    ##
    ## brood.table = sockeye salmon brood table (need columns $BY $Stock.ID $Rx.x)
    ## env.data = data.frame of pink data (need column ($Year)
    ## env.covar = column name in `pink.data` for calc covar with
    ## type = c("first_year", "second_year")
    ## out.covar = name of column for new covar in output data.frame


    if("year" %in% names(env.data))
        names(env.data)[names(env.data) == "year"] <- "Year"

    env.mat <- matrix(NA, length(unique(brood.table$BY)),
                      length(unique(brood.table$Stock.ID)),
                      dimnames=list(unique(brood.table$BY),unique(brood.table$Stock.ID)))
    env.mat <- env.mat[order(row.names(env.mat)), ]

    for (i in unique(brood.table$Stock.ID)){
        brood <- subset(brood.table, Stock.ID == i)
        if(type == "first_year")
            climate <- subset(env.data, Stock.ID == i)
        for (j in unique(brood$BY)){

            if(type == "first_year") {
                env.mat[as.character(j),as.character(i)] <-
                    brood$ocean_0[brood$BY == j] * climate[climate$Year == j+1, env.covar] +
                    brood$ocean_1[brood$BY == j] * climate[climate$Year == j+2, env.covar] +
                    brood$ocean_2[brood$BY == j] * climate[climate$Year == j+3, env.covar] 
            }
            if(type == "second_year") {
                env.mat[as.character(j),as.character(i)] <-
                    brood$ocean_0[brood$BY == j] * env.data[env.data$Year == j+2, env.covar] +
                    brood$ocean_1[brood$BY == j] * env.data[env.data$Year == j+3, env.covar] +
                    brood$ocean_2[brood$BY == j] * env.data[env.data$Year == j+4, env.covar] +
                    brood$ocean_3[brood$BY == j] * env.data[env.data$Year == j+5, env.covar] +
                    brood$ocean_4[brood$BY == j] * env.data[env.data$Year == j+6, env.covar]
            }
        } # end j loop
    } # end i loop
    long.df <- reshape2::melt(env.mat,
                              measure.vars=c((min(brood.table$BY):max(brood.table$BY)),
                                             unique(brood.table$Stock.ID)))
    colnames(long.df) <- c("BY","Stock.ID",out.covar)
    long.df <- long.df[complete.cases(long.df), ]
    return(long.df)
}




## pink.wgt.avg --------------------------------------------
pink.wgt.avg <- function(brood.table,
                         pink.data,
                         pink.covar,
                         type,
                         out.covar = "pink_index") {
    ## Calculate pink abundance index weighted by age-class
    ##
    ## Create data frame with estimates of competitors weighted by
    ## the proportion of recruits for each age-class. Competitor index is either
    ## (1) cumulative competitor abundance while at sea, (2) just second year
    ## competitor abundance, or (3) a geographically varying second year
    ## competitor index (by+4 for WC and GOA stocks and by+3 for BS stocks)
    ##
    ## brood.table = sockeye salmon brood table (need columns $BY $Stock.ID $Rx.x)
    ## pink.data = data.frame of pink data (need column ($Year)
    ## pink.covar = column name in `pink.data` for calc covar with
    ## type = c("cumulative", "second_year", "geographic")
    ## out.covar = name of column for new covar in output data.frame

    np.pink <- matrix(NA,length(unique(brood.table$BY)),
                          length(unique(brood.table$Stock.ID)),
                          dimnames=list(unique(brood.table$BY), unique(brood.table$Stock.ID)))
    np.pink <- np.pink[order(row.names(np.pink)), ]


    for (i in unique(brood.table$Stock.ID)){
        brood <- subset(brood.table, Stock.ID == i)
        for (j in unique(brood$BY)){
            if(type == "second_year") {
                np.pink[as.character(j),as.character(i)] <- pink.data[pink.data$Year == j+4, pink.covar]            }
        } # end j loop
    } # end i loop

    long.df <- reshape2::melt(np.pink,
                              measure.vars=c((min(brood.table $BY):max(brood.table $BY)),
                                             unique(brood.table$Stock.ID)))
    colnames(long.df) <- c("BY","Stock.ID",out.covar)
    return(long.df)
}

if(FALSE) {
    bt.raw <- read.csv("data/master_brood_table.csv", header=T)
    bt.complete <- bt.raw[complete.cases(bt.raw),]
    raw.comp <- read.csv(file="data-downloaded/pink_abundance_2017_12_08.csv",
                         header=TRUE)

    cu <- pink.wgt.avg(bt.complete, raw.comp, "Total",
                       "cumulative", "cum_np_pinks")
    np <- pink.wgt.avg(bt.complete, raw.comp, "Total",
                       "second_year", "np_pinks")
    go <- pink.wgt.avg(bt.complete, raw.comp, "Total",
                       "geographic", "geo_np_pinks")

    # all.equal(cu, long.cum.np.pink)
    # all.equal(np, long.np.pink)
    # all.equal(go, long.geo.np.pink)
}


## fill.time.series ----------------------------------------
fill.time.series <- function(data) {

    ## Fill salmon data time series so that all brood years are consecutive,
    ## filling missing years with NA.
    ##
    ## This function takes as input, a brood table with columns `Stock.ID`,
    ## `Stock`, and `BY` and fills in any non-consecutive BY within a stocks
    ## times series with NA values for all other columns present in `data`.
    ##
    ## When filtering the data by the `use` column, data values in the middle of
    ## the time series for a particular salmon stocks also get removed if `use`
    ## is set to 0. This function adds back in those data points, setting them
    ## to NA.

    id <- unique(data$Stock.ID)
    lst <- vector("list", length(id))
    for(i in seq_along(id)) {
        sub <- data[data$Stock.ID == id[i], ]
        BY <- min(sub$BY):max(sub$BY)
        Stock.ID <- unique(sub$Stock.ID)
        Stock <- unique(sub$Stock)
        df <- data.frame(Stock.ID = Stock.ID, Stock = Stock, BY = BY,
                         stringsAsFactors = FALSE)
        lst[[i]] <- merge(df, sub, by = c("Stock.ID", "Stock", "BY"), all.x = TRUE)
    }
    df <- do.call("rbind", c(lst, make.row.names = FALSE))

    ## Don't want NA in these columns
    out <- plyr::ddply(df, .(Stock.ID), transform,
                       Region = unique(na.omit(Region)),
                       Sub.Region = unique(na.omit(Sub.Region)),
                       Ocean.Region = unique(na.omit(Ocean.Region)),
                       Lat = unique(na.omit(Lat)),
                       Lon = unique(na.omit(Lon)))
    return(out)
}



## get.npgo ------------------------------------------------
get.npgo <- function(years) {

    ## This function takes as input a range of years and downloads and processes
    ## the NPGO index. The output of the function is a dataframe in 'long'
    ## format with a column for year, month, and the NPGO index.
    ##
    ## years = vector of years

    if(min(years) < 1950)
        stop("Earliest NPGO year is 1950")

    npgo    <- read.table("http://www.o3d.org/npgo/npgo.php", sep = "\t",
                          strip.white = TRUE)
    n.npgo  <- length(npgo[ , 1])
    npgo    <- npgo[4:n.npgo, ]
    n.npgo  <- length(npgo)
    rm.tail <- n.npgo - 3
    npgo    <- npgo[1:rm.tail]
    npgo    <- as.character(npgo)
    npgo    <- strsplit(npgo, "  ")
    npgo    <- do.call("rbind", npgo)
    npgo    <- data.frame(year = as.numeric(npgo[ , 1]),
                          month = as.numeric(npgo[ , 2]),
                          npgo = as.numeric(npgo[ , 3]))
    npgo    <- npgo[npgo$year >= min(years) & npgo$year <= max(years), ]

    return(npgo)
}

if(FALSE) {

    get.npgo(1950:1950)
    get.npgo(1950:2013)

    npgo <- get.npgo(1950:2016)
    head(npgo)
    tail(npgo)
    sapply(npgo, class)
    summary(npgo)

}



## sst.averager --------------------------------------------
sst.averager <- function(info, sst, distance = 400) {

    ## This function takes as input a data.frame of sst data output from the
    ## sst.anomaly() function and computes sst averages for each stock only
    ## including sst grid cells that are within a specified distance from the
    ## ocean entry location of a particular stock.
    ##
    ## Different months are used to calculate the SST average based on where the
    ## salmon stock enters the ocean:
    ##  WA, BC, SEAK: April-July
    ##  GOA: May-August
    ##  BB, AYK: Jun-September
    ##
    ## The function outputs a data frame with columns:
    ##  $year = one year for each year in input SST data and stock.id
    ##  $stock.id = id number for a salmon stock
    ##  $sst = averaged raw SST values
    ##  $sst.anom = averaged SST anomalies
    ##
    ## Function arguments:
    ##   info = stock.info data.frame w/ stock.id number, lon, and lat, should
    ##          have one row per stock
    ##   sst = sst data output from sst.anomaly()
    ##   distance = distance in km from ocean entry location of stock a grid
    ##              cell can be to be included in the averaging. This distance
    ##              is measured to the center of the SST grid cell

    info$stock.id <- info$Stock.ID
    stock.id <- info$stock.id
    cells    <- unique(subset(sst, select = c(id, lat, lon2)))
    cells    <- cells[order(cells$id), ]
    n.cells  <- length(cells[ , 1])
    row.names(cells) <- NULL

    sst.out <- vector("list", length(stock.id))
    for(i in seq_along(stock.id)) {

        info.i     <- info[i , ]
        stock.id.i <- info.i$stock.id
        lat.i      <- info.i$lat
        lon.i      <- info.i$lon # Brendan removed "* -1" from this line of code; longitude was already converted to negative degrees east

        dist <- rep(NA, n.cells)
        for(j in 1:n.cells)
            dist[j] <- haversine(lat.i, lon.i, cells$lat[j], cells$lon2[j])

        cells.sub <- cells[which(dist <= distance), ]
        sst.sub   <- sst[sst$id %in% cells.sub$id, ]

        if(stock.id.i <= 137)
            months <- 4:7  ## WA, BC, SEAK

        if(stock.id.i > 137 & stock.id.i <= 152)
            months <- 5:8  ## GOA

        if(stock.id.i > 152)
            months <- 6:9  ## BB and AYK

        sst.sub.mnths <- sst.sub[sst.sub$month %in% months, ]

        sst.avg <- ddply(sst.sub.mnths, .(year), summarize,
                               sst = mean(sst, na.rm = TRUE),
                               sst.anom = mean(sst.anom, na.rm = TRUE))

        sst.avg$stock.id <- stock.id.i

        sst.out[[i]] <- sst.avg
    }
    sst.out <- rbind.fill(sst.out)
    return(sst.out)
}



## haversine -----------------------------------------------
haversine <- function(lat1, lon1, lat2, lon2) {
    ## This function computes the great circle distance between two points given
    ## their latitiude and longitude (in decimal degrees) using the haversine
    ## formula. The output is the distance between the two points in km.
    ##
    ## lat1 = latitude of first point
    ## lon1 = longitude of first point
    ## lat2 = latitude of second point
    ## lon2 = longitude of second point

    # Convert degrees to radians
    lat1 <- lat1 * pi / 180
    lon1 <- lon1 * pi / 180
    lat2 <- lat2 * pi / 180
    lon2 <- lon2 * pi / 180

    R <- 6371 # earth mean radius [km]
    delta.lon <- (lon2 - lon1)
    delta.lat <- (lat2 - lat1)
    a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.lon/2)^2
    d <- 2 * R * asin(min(1, sqrt(a)))

    return(d) # distance in km
}

if(FALSE) {

    haversine(48, -147, 49, -154)

}



## enviro.avg.months ---------------------------------------
enviro.avg.months <- function(data,
                              first.month,
                              last.month,
                              avg.var,
                              month.var = "month",
                              year.var = "year",
                              grid.id.var = NULL,
                              lat.var = "lat",
                              lon.var = "lon",
                              type = avg.var) {
    ## Average an environmental variable over the specified months within a year
    ##
    ## Input data should be in a "long" format with a column for year, month,
    ## and the environmental variable to average. `first.month` and
    ## `last.month` give the first and last months of a continuous range to
    ## average the environmental variable over.
    ##
    ## If `first.month` is less than `last.month` the environmental variable is
    ## averaged over the months first.month:last.month within each year. For
    ## example, if `first.month` = 3 and `last.month` = 4, the environmental
    ## variable will be averaged over Mar and Apr for each year.
    ##
    ## If `first.month` equals `last.month`, that month is returned with no
    ## averaging.
    ##
    ## If `first.month` is greater than `last.month` the environmental variable
    ## is averaged for year t starting in `first.month` of year t - 1 and ending
    ## in `last.month` of year t. The output year corresponds to the year
    ## January occurs within the average. For example if `first.month` = 12 and
    ## `last.month` = 3, then the average for the environmental variable will
    ## occur over Dec, Jan, Feb, March and the year is specified by the year
    ## for Jan occurs in.
    ##
    ## The function outputs a data.frame with a `year` column and an `index`
    ## column. When `first.month` is greater than `last.month`, the output
    ## data.frame will have one less year than the input data.frame with no
    ## value for the minimum year in the input data frame.
    ##
    ## If `grid.id.var` is non-null, the averaging is done on a per-grid-cell
    ## basis within a year.
    ##
    ## data = a data.frame
    ## first.month = numeric giving the month to start annual average
    ## last.month = numeric giving the month to stop annual average
    ## avg.var = string, column name in `data` of the variable to average
    ## month.var = string, column name in `data` of the month variable
    ## year.var = string, column name in `data` of the year variable
    ## grid.id.var = string, column name in `data` of the grid cell id column
    ## lon.var = string, column name in `data` of for longitude, only used if
    ##           grid.id.var is non-null
    ## lat.var = string, column name in `data` of for latitude, only used if
    ##           grid.id.var is non-null
    ## type = string, value to set a `type` column in the output data.frame.
    ##        Useful if you want to rbind multiple indices together

    if(!is.data.frame(data))
        stop("Input data is not a data.frame", call. = FALSE)

    if(!first.month %in% 1:12 | !last.month %in% 1:12)
        stop("Months not between 1 and 12", call. = FALSE)

    if(!is.numeric(data[ , month.var]) | !is.numeric(data[ , year.var]))
       stop("Month variable must be numeric", call. = FALSE)

    if(first.month < last.month | first.month == last.month) {
        months <- first.month:last.month
        df <- data[data[ , month.var] %in% months, ]
    }

    if(first.month > last.month) {

        ## Remove months prior to `first.month` in first year and
        ## months after `last.month` in last year
        min.yr <- min(data[ , year.var])
        max.yr <- max(data[ , year.var])
        min.rm <- which(data[ , year.var] == min.yr &
                        data[ , month.var] < first.month)
        max.rm <- which(data[ , year.var] == max.yr &
                        data[ , month.var] > last.month)
        sub <- data[-c(min.rm, max.rm), ]

        ## Remove months not being averaged over
        months <- c(first.month:12, 1:last.month)
        sub2 <- sub[sub[ , month.var] %in% months, ]

        ## Create new year index to average over
        sp <- split(sub2, sub2[ , year.var])
        lst <- lapply(sp, function(x) {
                   x$yr.avg <- ifelse(x[ , month.var] %in% first.month:12,
                                      x[ , year.var] + 1, x[ , year.var])
                   return(x)
               })
        df <- do.call("rbind", c(lst, make.row.names = FALSE))
        df[ , year.var] <- df$yr.avg

    }

    ## Calculate averages
    if(is.null(grid.id.var)) {
        sp.avg <- split(df, df[ , year.var])
        lst.avg <- lapply(sp.avg, function(x) {
                    data.frame(year = unique(x[ , year.var]),
                               index = mean(x[ , avg.var]))
                })
        enviro <- do.call("rbind", c(lst.avg, make.row.names = FALSE))
    } else {
        sp.avg <- split(df, list(df[ , grid.id.var], df[ , year.var]))
        lst.avg <- lapply(sp.avg, function(x) {
                    data.frame(year = unique(x[ , year.var]),
                               id = unique(x[ , grid.id.var]),
                               lon = unique(x[ , lon.var]),
                               lat = unique(x[ , lat.var]),
                               index = mean(x[ , avg.var]))
                })
        enviro <- do.call("rbind", c(lst.avg, make.row.names = FALSE))
    }

    enviro$type <- type

    char.months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
                     "sep", "oct", "nov", "dec")
    enviro$months <- paste(char.months[first.month],
                           char.months[last.month], sep = "-")
    enviro$months <- as.factor(enviro$months)

    return(enviro)
}

if(FALSE) {

    set.seed(101)
    yr <- c(rep(1950, 24), rep(1951, 24))
    id <- c(rep(1, 12), rep(2, 12), rep(1, 12), rep(2, 12))
    mt <- rep(1:12, 4)
    lt <- c(rep(48, 12), rep(50, 12), rep(48, 12), rep(50, 12))
    ln <- c(rep(120, 12), rep(122, 12), rep(120, 12), rep(122, 12))
    vl <- rnorm(48)
    df <- data.frame(year = yr, month = mt, lat = lt, lon = ln, id = id, value = vl)

    enviro.avg.months(df, 1, 12, "value", grid.id.var = NULL)
    enviro.avg.months(df, 1, 12, "value", grid.id.var = "id")

    enviro.avg.months(df, 1, 3, "value", grid.id.var = NULL)
    enviro.avg.months(df, 1, 3, "value", grid.id.var = "id")

    enviro.avg.months(df, 12, 3, "value", grid.id.var = NULL)
    enviro.avg.months(df, 12, 3, "value", grid.id.var = "id")

    enviro.avg.months(df, 1, 3, "value")
    enviro.avg.months(df, 3, 3, "value")
    enviro.avg.months(df, 12, 3, "value")
    enviro.avg.months(df, 6, 5, "value", group = "test")
}



## get.pdo -------------------------------------------------
get.pdo <- function(years = 1900:2016) {

    ## This function takes as input a range of years and downloads and processes
    ## the PDO index. The output of the function is a dataframe in 'long' format
    ## with a column for year, month, and the PDO index.
    ##
    ## years = vector of years

    if(min(years) < 1900)
        stop("Earliest PDO year is 1900")

    ## Read data
    pdo.r <- readLines("http://jisao.washington.edu/pdo/PDO.latest.txt")

    ## Trim header and footer
    start.line <- grep("^YEAR", pdo.r)
    end.line <- grep(paste0("^", max(years)), pdo.r)
    pdo.s <- pdo.r[start.line:end.line]

    ## Split strings
    pdo.c <- gsub("[/*]", "", pdo.s)
    pdo.c <- gsub("   ", "  ", pdo.c)
    pdo.c <- gsub("    ", "  ", pdo.c)
    pdo.c <- gsub("      ", "  ", pdo.c)
    pdo.c <- strsplit(pdo.c, "  ")
    pdo.c <- do.call("rbind", pdo.c)

    ## Convert to data.frame
    df.names <- trimws(tolower(pdo.c[1, ]), "both")
    pdo.d <- pdo.c[2:nrow(pdo.c), ]
    pdo.d <- apply(pdo.d, 2, function(x) as.numeric(x))
    pdo.d <- as.data.frame(pdo.d)
    names(pdo.d) <- df.names

    ## Reshape into "long" format
    months <- df.names[2:length(df.names)]
    pdo <- reshape(pdo.d, direction = "long",
                   varying = months,
                   v.names = "pdo",
                   times = months,
                   timevar = "month")
    pdo <- pdo[order(pdo$year), ]
    pdo <- pdo[pdo$year >= min(years) & pdo$year <= max(years), ]
    row.names(pdo) <- NULL
    pdo$id <- NULL

    ## Convert months to numerics
    pdo$month <- ifelse(pdo$month == "jan", 1, pdo$month)
    pdo$month <- ifelse(pdo$month == "feb", 2, pdo$month)
    pdo$month <- ifelse(pdo$month == "mar", 3, pdo$month)
    pdo$month <- ifelse(pdo$month == "apr", 4, pdo$month)
    pdo$month <- ifelse(pdo$month == "may", 5, pdo$month)
    pdo$month <- ifelse(pdo$month == "jun", 6, pdo$month)
    pdo$month <- ifelse(pdo$month == "jul", 7, pdo$month)
    pdo$month <- ifelse(pdo$month == "aug", 8, pdo$month)
    pdo$month <- ifelse(pdo$month == "sep", 9, pdo$month)
    pdo$month <- ifelse(pdo$month == "oct", 10, pdo$month)
    pdo$month <- ifelse(pdo$month == "nov", 11, pdo$month)
    pdo$month <- ifelse(pdo$month == "dec", 12, pdo$month)
    pdo$month <- as.numeric(pdo$month)

    return(pdo)
}

if(FALSE) {

    get.pdo(1950:1950)
    get.pdo(1950:2013)

    pdo <- get.pdo(1900:2016)
    head(pdo)
    tail(pdo)
    sapply(pdo, class)
    summary(pdo)

}


## sst.anomaly ---------------------------------------------
sst.anomaly <- function(data, ref.years) {
    ## Calculate monthly, per grid cell, anomalies of SST
    ##
    ## Anomalies are calculated as the difference between a grid cell specific
    ## SST value for a given year/month and the long-term monthly mean (defined
    ## by ref.years) for that grid cell. This follows the methods outlined in
    ## Mueter et al. 2002, CJFAS (https://doi.org/10.1139/f02-020).
    ##
    ## data = data.frame of SST data
    ## ref.years = reference years in which to calculate the long-term mean,
    ##             should be a continuous sequence, e.g., 1950:2016

    ## make sure ref.years is a continuous sequence
    if(length(ref.years) != length(min(ref.years):max(ref.years)))
        stop("years vector is not a sequence ascending by 1")

    ## make sure ref.years are available in input data
    if(sum(ref.years %in% data$year) != length(ref.years))
        stop("ref.years not contained in input data")

    ## subset reference years from data
    ref.sst <- data[data$year >= min(ref.years) & data$year <= max(ref.years), ]

    ## calculate monthly long-term mean for each grid cell
    ## NA's are removed in the calculation of long-term mean
    mnth.avg <- aggregate(sst ~ month + id, data = ref.sst,
                          function(x) mean(x, na.rm = TRUE),
                          na.action = na.pass)
    names(mnth.avg)[names(mnth.avg) == 'sst'] <- 'long.avg'
    sst.merge <- merge(data, mnth.avg)
    sst.merge <- sst.merge[order(sst.merge$year,
                                 sst.merge$month,
                                 sst.merge$lat,
                                 sst.merge$lon), ]

    ## calculate sst anomaly
    sst.merge$anom <- sst.merge$sst - sst.merge$long.avg
    row.names(sst.merge) <- NULL
    sst <- data.frame(year = sst.merge$year,
                      month = sst.merge$month,
                      lon = sst.merge$lon,
                      lat = sst.merge$lat,
                      id = sst.merge$id,
                      sst = sst.merge$sst,
                      sst.anom = sst.merge$anom)
    return(sst)
}



## ggplot: theme_sleek -------------------------------------
# see: https://github.com/seananderson/ggsidekick
theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = rel(1.1)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}



## lattice: theme.mjm --------------------------------------
theme.mjm <- function(fontsize = 11, ...) {
    ## see latticeExtra::custom.theme
    ## trellis.par.get()

    col.symbol <- "#00A0EDFF"

    col.superpose     <- c("#00A0EDFF", "#DE732DFF", "#00AB06FF", "#BE6CF5FF",
                           "#A09500FF", "#F352B1FF", "#00B3A3FF")
    col.superpose.pol <- c("#00B9FFFF", "#FA8F59FF", "#0DC649FF", "#D58CFFFF",
                           "#BBAF00FF", "#FF75CAFF", "#00CCBDFF")

    col.regions <- c("#004B9FFF", "#004AA2FF", "#0049A4FF", "#0047A6FF",
                     "#0046A8FF", "#0045A9FF", "#0043ABFF", "#0042ADFF",
                     "#0041AEFF", "#003FB0FF", "#003EB1FF", "#003CB2FF",
                     "#003AB4FF", "#0039B5FF", "#0B37B6FF", "#2936B7FF",
                     "#3934B7FF", "#4633B8FF", "#5031B9FF", "#592FBAFF",
                     "#622EBAFF", "#692DBBFF", "#702BBBFF", "#772ABBFF",
                     "#7D29BCFF", "#8328BCFF", "#8827BCFF", "#8D26BCFF",
                     "#9325BCFF", "#9725BCFF", "#9C25BCFF", "#A125BBFF",
                     "#A525BBFF", "#A925BBFF", "#AD26BAFF", "#B127BAFF",
                     "#B528B9FF", "#B929B9FF", "#BC2BB8FF", "#C02DB7FF",
                     "#C32EB6FF", "#C730B5FF", "#CA32B4FF", "#CD35B3FF",
                     "#D037B2FF", "#D339B1FF", "#D63CB0FF", "#D93EAFFF",
                     "#DB41ADFF", "#DE44ACFF", "#E146ABFF", "#E349A9FF",
                     "#E54CA8FF", "#E84FA6FF", "#EA52A4FF", "#EC54A3FF",
                     "#EE57A1FF", "#F05A9FFF", "#F25D9DFF", "#F4609BFF",
                     "#F66399FF", "#F86697FF", "#FA6994FF", "#FB6D92FF",
                     "#FD7090FF", "#FE738DFF", "#FF768BFF", "#FF7988FF",
                     "#FF7C86FF", "#FF7F83FF", "#FF8280FF", "#FF867EFF",
                     "#FF897BFF", "#FF8C78FF", "#FF8F75FF", "#FF9271FF",
                     "#FF966EFF", "#FF996BFF", "#FF9C68FF", "#FF9F64FF",
                     "#FFA261FF", "#FFA65DFF", "#FFA959FF", "#FFAC56FF",
                     "#FFAF52FF", "#FFB34EFF", "#FFB64AFF", "#FFB946FF",
                     "#FFBC42FF", "#FFC03EFF", "#FFC339FF", "#FFC635FF",
                     "#FFC931FF", "#FFCD2DFF", "#FFD029FF", "#FFD325FF",
                     "#FFD722FF", "#FFDA1FFF", "#FFDD1CFF", "#FFE01BFF")
    lwd <- 0.8

    theme <- list(
        fontsize          = list(text = fontsize),
        par.main.text     = list(font = 1),
        strip.background  = list(col  = c("grey93", "grey93")),
        strip.shingle     = list(col  = c("grey65", "grey65")),
        strip.border      = list(col  = "grey50", lwd = rep(lwd, 7)),
        axis.components   = list(right  = list(tck = 0.5),
                                 top    = list(tck = 0.5),
                                 left   = list(tck = 0.5),
                                 bottom = list(tck = 0.5)),
        axis.line         = list(col = "grey50", lwd = lwd),
        plot.symbol       = list(pch = 16, col = col.symbol),
        plot.line         = list(pch = 16, col = "grey10"),
        add.line          = list(col = "grey10"),
        superpose.symbol  = list(col = col.superpose, pch = 16),
        superpose.line    = list(col = col.superpose),
        regions           = list(col = col.regions),
        dot.symbol        = list(col = col.symbol),
        plot.polygon      = list(col = "grey50", border = "white"),
        superpose.polygon = list(col = col.superpose.pol, border = "white"),
        box.rectangle     = list(col = "grey40"),
        box.dot           = list(col = "grey10"),
        box.umbrella      = list(col = "grey40", lty = 1),
        box.3d            = list(col = "grey70", lwd = lwd))
    modifyList(modifyList(standard.theme("pdf"), theme), simpleTheme(...))
}
