##### How do the stock-specific metric approaches (independent and full models) compare with the global metric independent models?

#stock-specific indep
coefs_stan_psi <- read.csv("data/prev_coefs_stan_bystock_>4 sample sizes.csv")
colnames(coefs_stan_psi)<-c("agent", "lower","25","mid","75","upper")
coefs_stan_psi$mod<-"stock-specific indep"
coefs_stan_psi$metric<-"prev"

coefs_stan_lsi <- read.csv("data/load_coefs_stan_bystock_>4 sample sizes.csv")
colnames(coefs_stan_lsi)<-c("agent", "lower","25","mid","75","upper")
coefs_stan_lsi$mod<-"stock-specific indep"
coefs_stan_lsi$metric<-"load"

#global metric in SOG, independent models by agent
coefs_stan_pgi <- read.csv("data/prev_coefs_stan_global.csv")
colnames(coefs_stan_pgi)<-c("agent", "lower","25","mid","75","upper")
coefs_stan_pgi$mod<-"global indep"
coefs_stan_pgi$metric<-"prev"
coefs_stan_pgi<-droplevels(coefs_stan_pgi[!(coefs_stan_pgi$agent=="ku_thy"),]) #remove rows with extra agents
coefs_stan_pgi<-droplevels(coefs_stan_pgi[!(coefs_stan_pgi$agent=="smallUK"),]) #remove rows with extra agents
coefs_stan_pgi<-droplevels(coefs_stan_pgi[!(coefs_stan_pgi$agent=="prv"),]) #remove rows with extra agents
coefs_stan_pgi <- coefs_stan_pgi[order(-coefs_stan_pgi[,4]),]

coefs_stan_lgi <- read.csv("data/int_coefs_stan_global.csv")
colnames(coefs_stan_lgi)<-c("agent", "lower","25","mid","75","upper")
coefs_stan_lgi$mod<-"global indep"
coefs_stan_lgi$metric<-"load"
coefs_stan_lgi<-droplevels(coefs_stan_lgi[!(coefs_stan_lgi$agent=="ku_thy"),]) #remove rows with extra agents
coefs_stan_lgi<-droplevels(coefs_stan_lgi[!(coefs_stan_lgi$agent=="smallUK"),]) #remove rows with extra agents
coefs_stan_lgi<-droplevels(coefs_stan_lgi[!(coefs_stan_lgi$agent=="prv"),]) #remove rows with extra agents
coefs_stan_lgi <- coefs_stan_lgi[order(-coefs_stan_lgi[,4]),]

coefs_stan_psi<-coefs_stan_psi[order(match(coefs_stan_psi[,1],coefs_stan_pgi[,1])),]
coefs_stan_psf<-coefs_stan_psf[order(match(coefs_stan_psf[,1],coefs_stan_pgi[,1])),]
coefs_stan_lsi<-coefs_stan_lsi[order(match(coefs_stan_lsi[,1],coefs_stan_lgi[,1])),]
coefs_stan_lsf<-coefs_stan_lsf[order(match(coefs_stan_lsf[,1],coefs_stan_lgi[,1])),]




#### Plot prevalence coefficients
jpeg(filename='figs/Fig_SSHI ONNE productivity model type comparison.jpg', 
     width=480, height=500, quality=75)
agents <- unique(coefs_stan_lsi$agent)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))

#add PSI
plotCI(x = coefs_stan_psi[,4], col="green",
       y = seq(1,length(agents))+.25,
       li = (coefs_stan_psi[,2]),
       ui = (coefs_stan_psi[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-3,3),
       pch = 16,
       scol = "green")
plotCI(x = coefs_stan_psi[,4], col="green",
       y = seq(1,length(agents))+.25,
       li = (coefs_stan_psi[,3]),
       ui = (coefs_stan_psi[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "green")

#add PGI
plotCI(x = coefs_stan_pgi[,4], col="blue",
       y = seq(1,length(agents)),
       li = (coefs_stan_pgi[,2]),
       ui = (coefs_stan_pgi[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-3,3),
       add = TRUE,
       pch = 16,
       scol = "blue")
plotCI(x = coefs_stan_pgi[,4], col="blue",
       y = seq(1,length(agents)),
       li = (coefs_stan_pgi[,3]),
       ui = (coefs_stan_pgi[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "blue")

text(rep(-2.75,length(agents)), 
     seq(1,length(agents)), 
     labels = (coefs_stan_psi[,1]), 
     pos = 4,
     font = 2,
     cex=0.95)
legend(
  x=0.5, # x coordinate of the top left of the legend
  y=21.5,
  legend=c("Stock-specific prevalence", "Global prevalence"), 
  pch=21,
  pt.bg=c("green","blue"))
axis(1, at = c(-2, -1, -0.5, 0, 0.5, 1, 2))
abline(v = 0, lty = 2)
box(col="grey") 
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Prevalence",3,line=0.25)
dev.off()



# INTENSITY
#### Plot load coefficients
jpeg(filename='figs/Fig_SSHI ONNE productivity model type comparison_load.jpg', 
     width=480, height=500, quality=75)
agents <- unique(coefs_stan_lsi$agent)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))

#add PSI
plotCI(x = coefs_stan_lsi[,4], col="green",
       y = seq(1,length(agents))+.25,
       li = (coefs_stan_lsi[,2]),
       ui = (coefs_stan_lsi[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-3,4),
       pch = 16,
       scol = "green")
plotCI(x = coefs_stan_lsi[,4], col="green",
       y = seq(1,length(agents))+.25,
       li = (coefs_stan_lsi[,3]),
       ui = (coefs_stan_lsi[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "green")

#add PGI
plotCI(x = coefs_stan_lgi[,4], col="blue",
       y = seq(1,length(agents)),
       li = (coefs_stan_lgi[,2]),
       ui = (coefs_stan_lgi[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-3,4),
       add = TRUE,
       pch = 16,
       scol = "blue")
plotCI(x = coefs_stan_lgi[,4], col="blue",
       y = seq(1,length(agents)),
       li = (coefs_stan_lgi[,3]),
       ui = (coefs_stan_lgi[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "blue")

text(rep(-2.75,length(agents)), 
     seq(1,length(agents)), 
     labels = (coefs_stan_lsi[,1]), 
     pos = 4,
     font = 2,
     cex=0.95)
legend(
  x=1.5, # x coordinate of the top left of the legend
  y=21,
  legend=c("Stock-specific intensity", "Global intensity"), 
  pch=21,
  pt.bg=c("green","blue"))
axis(1, at = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3))
abline(v = 0, lty = 2)
box(col="grey") 
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Intensity",3,line=0.25)
dev.off()
