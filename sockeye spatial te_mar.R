# A. Bass
## spatial evaluation for agents in sw

# load packages
library(postpack)
library(StatonMisc)
library(lme4)
library(tidyr)
library(dplyr)
library(raster)


dat<-sw.data

dat$yr_mo<-format(as.POSIXct(strptime(dat$Date,format="%d-%b-%y")),format="%Y-%m")

# we need a spatial dataset
dat2<-dat%>%
  filter(!is.na(Latitude),
         !is.na(Longitude))

coordinates<-dat2[,c(23,22)]
metadata<-dat2[,-c(23,22)]
trans<- CRS("+proj=longlat +datum=WGS84")
skpoi<-SpatialPointsDataFrame(coords=coordinates, data=metadata, proj4string=trans) # define GCS

r<-readRDS("data/aqua_sway_rasterstack.rds") ##need this
projection(r)<-trans

# match ac_dist to points
dat2$nearest_AC_sway<-NA
for(i in unique(dat2$yr_mo)){
  r_layer<-match(paste(substr(i,1,4),sep=".",substr(i,6,7)),substr(names(r),11,17))
  r_points<-which(dat2$yr_mo==i)
  if(length(r_points)>0){
    ac_dist<-raster::extract(r[[r_layer]],matrix(coordinates(skpoi)[r_points,],ncol=2),method="bilinear")
    dat2$nearest_AC_sway[r_points]<-ac_dist
  }
}

# going back in to deal with points that were not covered using buffering
need_aqua<-which(is.na(dat2$nearest_AC_sway) & dat2$Longitude>-130) # only worried about those in the SDM area
for(i in unique(dat2$yr_mo[need_aqua])){
  r_layer<-match(paste(substr(i,1,4),sep=".",substr(i,6,7)),substr(names(r),11,17))
  r_points<-which(dat2$yr_mo==i & is.na(dat2$nearest_AC_sway))
  if(length(r_points)>0){
    ac_dist<-raster::extract(r[[r_layer]],matrix(coordinates(skpoi)[r_points,],ncol=2), buffer=5000,fun=mean,na.rm=TRUE) #using buffer up to 5 km
    dat2$nearest_AC_sway[r_points]<-ac_dist
  }
}

windows()
  ggplot(dat2,aes(y=Latitude, x=Longitude, col=nearest_AC_sway)) +
  geom_point() + xlim(-129,-123) + ylim(47.5,52)

plot(dat2$nearest_AC_sway,log(dat2$te_mar))
fit<-lm(log(dat2$te_mar+1)~dat2$nearest_AC_sway)
abline(fit)

dat2$te_mar_pres<-ifelse(dat2$te_mar>0,1,0)

### here is the PRV code from Gideon's paper
dat2$ac_km<-dat2$nearest_AC_sway/1000 # convert to km
dat2$jaz<-as.factor(dat2$Stock_Group) # going to use this as the random stock effect

dat3<-dat2 %>%
  filter(!is.na(ac_km),
         !is.na(te_mar),
         !is.na(jaz)) %>%
  arrange(jaz)

# compile data into a list for JAGS
pred_x = seq(0,200,10)
jags_data = list(
  n_obs = nrow(dat3), 
  x = dat3$ac_km, 
  y = dat3$te_mar_pres, 
  grp = as.numeric(dat3$jaz), 
  n_grp = length(unique(dat3$jaz)),
  pred_x = pred_x,
  n_pred = length(pred_x)
)

##### STEP 2: SPECIFY JAGS MODEL CODE #####
jags_model = function() {
  # HYPER-PRIORS (PRIORS ON THE HYPERDIST)
  B0 ~ dnorm(0, 1e-3)
  B1 ~ dnorm(0, 1e-3)
  sig_B0 ~ dunif(0, 10)
  sig_B1 ~ dunif(0, 10)
  tau_B0 <- 1/sig_B0^2
  tau_B1 <- 1/sig_B1^2
  
  # JAZ / CU EFFECTS
  for (j in 1:n_grp) {
    b0[j] ~ dnorm(B0, tau_B0)
    b1[j] ~ dnorm(B1, tau_B1)
  }
  
  # LIKELIHOOD
  for(i in 1:n_obs){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- b0[grp[i]] + b1[grp[i]]*x[i]
  }
  
  # DERIVED QUANTITIES
  for(j in 1:n_grp){
    for(i in 1:n_pred){
      logit(pred_p[i,j])<-b0[j] + b1[j]*pred_x[i]  
    }
  }  
  
}

# write model to a text file
jags_file = "model.txt"
write_model(jags_model, jags_file)

##### STEP 3: SPECIFY INITIAL VALUES #####
# fit model with R to get reasonable initial values
fit = with(jags_data, lme4::glmer(y ~ x + (1|grp) + (0 + x|grp),family="binomial"))
coefs = lme4::fixef(fit)
coef_names = names(coefs)
int = stringr::str_detect(coef_names, "Intercept")
slope = stringr::str_detect(coef_names, "x")
b0_fit = unname(coefs[int]) 
b1_fit = unname(coefs[slope]) 

jags_inits = function(nc) {
  inits = list()
  for (c in 1:nc) {
    inits[[c]] = list(
      B0 = rnorm(1, b0_fit, 1e-3),
      B1 = rnorm(1, b1_fit, 1e-3),
      sig_B0 = runif(1, 0, 1),
      sig_B1 = runif(1, 0, 1),
      b0 = rnorm(jags_data$n_grp, b0_fit, abs(b0_fit) * 0.2),
      b1 = rnorm(jags_data$n_grp, b1_fit, abs(b1_fit) * 0.2)
    )
  }
  return(inits)
}

##### STEP 4: SET NODES TO MONITOR #####
jags_params = c("B0","B1","b0", "b1", "pred_p")

##### STEP 5: SPECIFY MCMC DIMENSIONS #####
jags_dims = c(
  ni = 10000,  # number of post-burn-in samples per chain
  nb = 30000,  # number of burn-in samples
  nt = 5,     # thinning rate
  nc = 2      # number of chains
)

##### STEP 6: RUN THE MODEL WITH JAGS #####
# should take less than 1 minute
post = jagsUI::jags.basic(
  data = jags_data,
  model.file = jags_file,
  inits = jags_inits(jags_dims["nc"]),
  parameters.to.save = jags_params,
  n.adapt = 1000,
  n.iter = sum(jags_dims[c("ni", "nb")]),
  n.thin = jags_dims["nt"],
  n.burnin = jags_dims["nb"],
  n.chains = jags_dims["nc"],
  parallel = F
)

##### STEP 7: CONVERGENCE DIAGNOSTICS #####
# view convergence diagnostic summaries for nodes with priors
t(post_summ(post, "b", ess = T, Rhat = T)[c("Rhat", "ess"),])

# view diagnostic plots
diag_plots(post, "b", ext_device = T)

b0_est = postpack::post_summ(post, "B0|b0",p_summ=c(0.5,0.025,0.25,0.75,0.95))
b1_est = postpack::post_summ(post, "B1|b1",p_summ=c(0.5,0.025,0.25,0.75,0.95))

##### STEP 8: MAKE INFERENCE #####
# view boxplots of posterior slope and intercept for each stock
ext_device(h = 8, w = 4)
par(mfrow = c(2,1), mar = c(2,2,2,2))
boxplot(post_subset(post, "b0", matrix = T), outline = F, col = "grey")
boxplot(post_subset(post, "b1", matrix = T), outline = F, col = "grey")

# view fitted curves for each year
pred_p = post_summ(post, "pred_p")
mean_pred_p = array_format(pred_p["mean",])

ext_device(h = 4, w = 4)
par(mar = c(4,4,1,1), xaxs = "i", yaxs = "i")
plot(1,1,type = "n", xlim = range(pred_x), ylim = c(0,1), xlab = "distance from aquaculture", ylab = "Probability positive")
sapply(1:jags_data$n_grp, function(t) lines(mean_pred_p[,t] ~ pred_x, col = t))

mean<-exp(b0_est[1,1] + b1_est[1,1]*pred_x)/(1 + exp(b0_est[1,1] + b1_est[1,1]*pred_x))
llci<-exp(b0_est[4,1] + b1_est[4,1]*pred_x)/(1 + exp(b0_est[4,1] + b1_est[4,1]*pred_x))
lci<-exp(b0_est[5,1] + b1_est[5,1]*pred_x)/(1 + exp(b0_est[5,1] + b1_est[5,1]*pred_x))
uci<-exp(b0_est[6,1] + b1_est[6,1]*pred_x)/(1 + exp(b0_est[6,1] + b1_est[6,1]*pred_x))
uuci<-exp(b0_est[7,1] + b1_est[7,1]*pred_x)/(1 + exp(b0_est[7,1] + b1_est[7,1]*pred_x))


# do it Andrew's way
core<-post_subset(post,p="^B",matrix=T)
coreplus<-post_subset(post,p="^b",matrix=T)

preds<-data.frame(matrix(NA,ncol=length(pred_x),nrow=nrow(core)))
preds.RE<-preds

for(i in 1:nrow(core)){
  preds[i,]=plogis(core[i,1]+pred_x*core[i,2])
  stock=sample(1:9,1)
  stock_int=stock
  stock_slp=stock+9
  preds.RE[i,]=plogis(coreplus[i,stock_int] + pred_x*coreplus[i,stock_slp])
}

CIs_ss = apply(preds,2,function(x){quantile(x, probs=c(0.025,0.05,0.5,0.95,0.975))})    
CIs.RE_ss = apply(preds.RE,2,function(x){quantile(x, probs=c(0.025,0.05,0.5,0.95,0.975))}) 

ss_res<-cbind(pred_x,CIs_ss[3,],CIs_ss[1,],CIs_ss[5,],CIs.RE_ss[1,],CIs.RE_ss[5,])
colnames(ss_res)<-c("pred_x","mean","lci","uci","llci","uuci")

#saveRDS(ss_res,"te_mar_res.rds")

breaks<-seq(0,150000,15000)

h_ch<-hist(dat3$nearest_AC_sway[dat3$te_mar>0],breaks=breaks,plot=F)
h2_ch<-hist(dat3$nearest_AC_sway[!is.na(dat3$te_mar)],breaks=breaks,plot=F)
perc_temar<-h_ch$counts/h2_ch$counts

LCI<-rep(NA,10)
UCI<-rep(NA,10)
SE<-rep(NA,10)
for(i in 1:10){
  SE[i]=sqrt(perc_temar[i]*(1-perc_temar[i])/h2_ch$counts[i])
  LCI[i]=perc_temar[i]-qnorm(.975)*SE[i]
  UCI[i]=perc_temar[i]+qnorm(.975)*SE[i]  
}
ss<-data.frame(cbind(perc_temar,SE,h2_ch$counts))
ss$alph<-1-SE
ss$centers<-(breaks[1:10]+7500)/1000

scaleFUN <- function(x) sprintf("%.1f", x)

bays_ss<-as.data.frame(ss_res[1:16,])
#bays_ss<-as.data.frame(readRDS("C:/Users/artie/OneDrive/PBS_postdoc/GM PRV paper/paper PRV/distPRV_ss.rds"))

p1<-ggplot(data = ss, aes(x=centers,y=perc_temar)) + 
  geom_ribbon(inherit.aes=FALSE,data=bays_ss,aes(x=pred_x,ymin=llci,ymax=uuci),fill="#4682B475") +
  geom_ribbon(inherit.aes=FALSE,data=bays_ss,aes(x=pred_x,ymin=lci,ymax=uci),fill="#4682B465",alpha=0.5) +
  geom_bar(aes(alpha=alph),stat="identity") +
  geom_text(aes(x=centers,y=0.15,label=V3))+
  geom_line(data=bays_ss, aes(x=pred_x,y=mean),
            color="#4682B4ff",size=1.5) +
  scale_alpha(range=c(0.1,1)) + coord_cartesian(ylim=c(0,0.17)) +
  theme_classic() + #labs(title="Spring and summer",tag="B")+
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        legend.position = "none") +
  ylab(expression('Proportion '*italic(Tenacibaculum~maritimum)*' positive')) +
  xlab("Distance from nearest active salmon netpen")


