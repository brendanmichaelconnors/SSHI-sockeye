# SSHI ANALYSIS SOCKEYE STAGE II
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

# Create objects for analysis
agents <- unique(inf_agt_resid_data_gl$agent)
stocks <- unique(inf_agt_resid_data_gl$Stock)
years <- unique(inf_agt_resid_data_gl$Year)

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


# STAN Approach for Multi-level Modeling

# PREVALENCE - INDEPENDENT MODELS by AGENT

## SW metric averaged across all stocks per year - Independent models

### Create files for each agent
for(i in unique(inf_agt_resid_data_gl$agent)) {
  nam <- paste("df", i, sep = ".")
  assign(nam, inf_agt_resid_data_gl[inf_agt_resid_data_gl$agent==i,])
}

### Loop for STAN independent models
for(i in agents){
  data <- subset(inf_agt_resid_data_gl, agent==i)
  nam <- paste("mod", i, sep = ".")
  assign(nam, stan_lmer(resid_value ~ 0 + prev_std + (prev_std|Stock) +(1|Year), 
                        data = data,
                        adapt_delta=0.95,
                        REML = F))
}

# uninformed priors - for example:
prior_summary(mod.arena2)

## Derive coefficient estimates and save in .csv file
coefs_stan <- matrix(NA,
                     nrow = length(agents),
                     ncol = 5,
                     dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      pars = c("prev_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan, file="data/prev_coefs_stan_global_indep mod.csv")

# Load estimates from file (if not running full model) and assign rownames
coefs_stan <- read.csv("data/prev_coefs_stan_global_indep mod.csv")
rownames(coefs_stan) <- coefs_stan[,1]
coefs_stan <- coefs_stan[,-1]  

# Plot effect size per agent
coefs_order <- coefs_stan[order(coefs_stan[,3]),]
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))

plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,1]),
       ui = (coefs_order[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-1.75,1),
       pch = 16,
       scol = "grey")

plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,2]),
       ui = (coefs_order[,4]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "grey")

text(rep(-1.75,length(agents)), 
     seq(1,length(agents)), 
     labels = rownames(coefs_order), 
     pos = 4,
     font = 2,
     cex=0.95)

axis(1, at = c(-1, -0.5, 0, 0.5, 1))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Prevalence",3,line=0.25)


## Derive posterior estimates by stock

### Intercepts
coefs_stan_stk_int <- matrix(NA,
                             nrow = length(stocks),
                             ncol = 6,
                             dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\(\\Intercept) Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_int", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_int[i] <- cbind(ind_coef[c(1:18),c(4:8)], paste(i))))
}

### Slopes
coefs_stan_stk_slp <- matrix(NA,
                             nrow = length(stocks),
                             ncol = 6,
                             dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\prev_std Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_slp", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_slp[i] <- cbind(ind_coef[c(1:18),c(4:8)], paste(i))))
}

### Year intercepts
coefs_stan_stk_year <- matrix(NA,
                              nrow = length(years),
                              ncol = 6,
                              dimnames = list(years,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = "Year",
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_year", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_year[i] <- cbind(ind_coef[c(1:8),c(4:8)], paste(i))))
}

### Sigmas
coefs_stan_stk_sig <- matrix(NA,
                             nrow = length(stocks),
                             ncol = 6,
                             dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("sigma","Sigma"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_sig", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_sig[i] <- cbind(ind_coef[c(1:5),c(4:8)], paste(i))))
}

### Rhat and Neff
sims <-as.matrix(mod.arena2) #extract parameter names from any agent model 
dim(sims)
para_name <- c(colnames(sims), "mean_PPD", "log-posterior")
para_name

#### Rhat
coefs_stan_stk_rhat <- matrix(NA,
                              nrow = length(para_name),
                              ncol = 2,
                              dimnames = list(para_name,c("Rhat","agent")))
for(i in agents){
  model<-get(paste("mod",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  nam <- paste("coefs_stan_stk_rhat", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_rhat[i] <- cbind(ind_coef[c(1:52),1], paste(i), paste("Rhat"))))
}

#### Neff
coefs_stan_stk_neff <- matrix(NA,
                              nrow = length(para_name),
                              ncol = 2,
                              dimnames = list(para_name,c("Neff","agent")))
for(i in agents){
  model<-get(paste("mod",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  nam <- paste("coefs_stan_stk_neff", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_neff[i] <- cbind(ind_coef[c(1:52),1], paste(i), paste("Neff"))))
}

#### Rbind all posteriors and save as .csv file
post_int_slpintsig <- rbind(coefs_stan_stk_int.arena2,
                            coefs_stan_stk_int.c_b_cys,
                            coefs_stan_stk_int.ce_sha,
                            coefs_stan_stk_int.de_sal,
                            coefs_stan_stk_int.fa_mar,
                            coefs_stan_stk_int.fl_psy,
                            coefs_stan_stk_int.ic_hof,
                            coefs_stan_stk_int.ic_mul,
                            coefs_stan_stk_int.ku_thy,
                            coefs_stan_stk_int.lo_sal,
                            coefs_stan_stk_int.my_arc,
                            coefs_stan_stk_int.pa_kab,
                            coefs_stan_stk_int.pa_min,
                            coefs_stan_stk_int.pa_pse,
                            coefs_stan_stk_int.pa_ther,
                            coefs_stan_stk_int.prv,
                            coefs_stan_stk_int.pspv,
                            coefs_stan_stk_int.rlo,
                            coefs_stan_stk_int.sch,
                            coefs_stan_stk_int.smallUK,
                            coefs_stan_stk_int.sp_des,
                            coefs_stan_stk_int.te_bry,
                            coefs_stan_stk_int.te_mar,
                            coefs_stan_stk_int.ven,
                            coefs_stan_stk_slp.arena2,
                            coefs_stan_stk_slp.c_b_cys,
                            coefs_stan_stk_slp.ce_sha,
                            coefs_stan_stk_slp.de_sal,
                            coefs_stan_stk_slp.fa_mar,
                            coefs_stan_stk_slp.fl_psy,
                            coefs_stan_stk_slp.ic_hof,
                            coefs_stan_stk_slp.ic_mul,
                            coefs_stan_stk_slp.ku_thy,
                            coefs_stan_stk_slp.lo_sal,
                            coefs_stan_stk_slp.my_arc,
                            coefs_stan_stk_slp.pa_kab,
                            coefs_stan_stk_slp.pa_min,
                            coefs_stan_stk_slp.pa_pse,
                            coefs_stan_stk_slp.pa_ther,
                            coefs_stan_stk_slp.prv,
                            coefs_stan_stk_slp.pspv,
                            coefs_stan_stk_slp.rlo,
                            coefs_stan_stk_slp.sch,
                            coefs_stan_stk_slp.smallUK,
                            coefs_stan_stk_slp.sp_des,
                            coefs_stan_stk_slp.te_bry,
                            coefs_stan_stk_slp.te_mar,
                            coefs_stan_stk_slp.ven,
                            coefs_stan_stk_sig.arena2,
                            coefs_stan_stk_sig.c_b_cys,
                            coefs_stan_stk_sig.ce_sha,
                            coefs_stan_stk_sig.de_sal,
                            coefs_stan_stk_sig.fa_mar,
                            coefs_stan_stk_sig.fl_psy,
                            coefs_stan_stk_sig.ic_hof,
                            coefs_stan_stk_sig.ic_mul,
                            coefs_stan_stk_sig.ku_thy,
                            coefs_stan_stk_sig.lo_sal,
                            coefs_stan_stk_sig.my_arc,
                            coefs_stan_stk_sig.pa_kab,
                            coefs_stan_stk_sig.pa_min,
                            coefs_stan_stk_sig.pa_pse,
                            coefs_stan_stk_sig.pa_ther,
                            coefs_stan_stk_sig.prv,
                            coefs_stan_stk_sig.pspv,
                            coefs_stan_stk_sig.rlo,
                            coefs_stan_stk_sig.sch,
                            coefs_stan_stk_sig.smallUK,
                            coefs_stan_stk_sig.sp_des,
                            coefs_stan_stk_sig.te_bry,
                            coefs_stan_stk_sig.te_mar,
                            coefs_stan_stk_sig.ven,
                            coefs_stan_stk_year.arena2,
                            coefs_stan_stk_year.c_b_cys,
                            coefs_stan_stk_year.ce_sha,
                            coefs_stan_stk_year.de_sal,
                            coefs_stan_stk_year.fa_mar,
                            coefs_stan_stk_year.fl_psy,
                            coefs_stan_stk_year.ic_hof,
                            coefs_stan_stk_year.ic_mul,
                            coefs_stan_stk_year.ku_thy,
                            coefs_stan_stk_year.lo_sal,
                            coefs_stan_stk_year.my_arc,
                            coefs_stan_stk_year.pa_kab,
                            coefs_stan_stk_year.pa_min,
                            coefs_stan_stk_year.pa_pse,
                            coefs_stan_stk_year.pa_ther,
                            coefs_stan_stk_year.prv,
                            coefs_stan_stk_year.pspv,
                            coefs_stan_stk_year.rlo,
                            coefs_stan_stk_year.sch,
                            coefs_stan_stk_year.smallUK,
                            coefs_stan_stk_year.sp_des,
                            coefs_stan_stk_year.te_bry,
                            coefs_stan_stk_year.te_mar,
                            coefs_stan_stk_year.ven)

write.csv(post_int_slpintsig, file="data/Posterior distributions_Int Slp Sig_global_indep mod_prev.csv")

#### Rbind all convergence parameters and save as .csv file
post_rhatneff_prev <- rbind(coefs_stan_stk_rhat.arena2,
                            coefs_stan_stk_rhat.c_b_cys,
                            coefs_stan_stk_rhat.ce_sha,
                            coefs_stan_stk_rhat.de_sal,
                            coefs_stan_stk_rhat.fa_mar,
                            coefs_stan_stk_rhat.fl_psy,
                            coefs_stan_stk_rhat.ic_hof,
                            coefs_stan_stk_rhat.ic_mul,
                            coefs_stan_stk_rhat.ku_thy,
                            coefs_stan_stk_rhat.lo_sal,
                            coefs_stan_stk_rhat.my_arc,
                            coefs_stan_stk_rhat.pa_kab,
                            coefs_stan_stk_rhat.pa_min,
                            coefs_stan_stk_rhat.pa_pse,
                            coefs_stan_stk_rhat.pa_ther,
                            coefs_stan_stk_rhat.prv,
                            coefs_stan_stk_rhat.pspv,
                            coefs_stan_stk_rhat.rlo,
                            coefs_stan_stk_rhat.sch,
                            coefs_stan_stk_rhat.smallUK,
                            coefs_stan_stk_rhat.sp_des,
                            coefs_stan_stk_rhat.te_bry,
                            coefs_stan_stk_rhat.te_mar,
                            coefs_stan_stk_rhat.ven,
                            coefs_stan_stk_neff.arena2,
                            coefs_stan_stk_neff.c_b_cys,
                            coefs_stan_stk_neff.ce_sha,
                            coefs_stan_stk_neff.de_sal,
                            coefs_stan_stk_neff.fa_mar,
                            coefs_stan_stk_neff.fl_psy,
                            coefs_stan_stk_neff.ic_hof,
                            coefs_stan_stk_neff.ic_mul,
                            coefs_stan_stk_neff.ku_thy,
                            coefs_stan_stk_neff.lo_sal,
                            coefs_stan_stk_neff.my_arc,
                            coefs_stan_stk_neff.pa_kab,
                            coefs_stan_stk_neff.pa_min,
                            coefs_stan_stk_neff.pa_pse,
                            coefs_stan_stk_neff.pa_ther,
                            coefs_stan_stk_neff.prv,
                            coefs_stan_stk_neff.pspv,
                            coefs_stan_stk_neff.rlo,
                            coefs_stan_stk_neff.sch,
                            coefs_stan_stk_neff.smallUK,
                            coefs_stan_stk_neff.sp_des,
                            coefs_stan_stk_neff.te_bry,
                            coefs_stan_stk_neff.te_mar,
                            coefs_stan_stk_neff.ven)

write.csv(post_rhatneff_prev, file="data/Posterior distributions_Rhat and Neff_global_indep mod_prev.csv")

### Plot posteriors per agent model from files
post_all <- read.csv("data/Posterior distributions_Int Slp Sig_global_indep mod_prev.csv")
post_all <- droplevels(post_all[!post_all$X.1 == "smallUK",]) #remove smallUK from analysis
post_agents <- read.csv("data/prev_coefs_stan_global_indep mod.csv")
post_agents <- droplevels(post_agents[!post_agents$X == "smallUK",])

## Plot Posteriors for all agents
#tiff('Fig_SSHI ONNE Pathogen Productivity_Agent slopes_Prevalence.tiff', 
#     units="in", width=5, height=6, res=300)
ggplot(post_agents) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(X, -mid), ymax = X75, ymin = X25), size=1.5, col="darkblue") +
  geom_linerange(aes(x = X, ymax = upper, ymin = lower), col="darkblue") +
  geom_point(aes(x = X, y = mid), size = 3, col="darkblue")+
  labs(x ="Infectious agents", y = "Effect size", title="Prevalence")+
  scale_x_discrete(labels=c("ic_mul" = "I. multifiliis", 
                            "te_mar" = "T. maritinum",
                            "pa_ther" = "P. theridion",
                            "fl_psy" = "F. psychrophilum",
                            "sch" = "Ca. S. salmonis",
                            "te_bry" = "T. bryosalmonae",
                            "pa_kab" = "P. kabatai",
                            "c_b_cys" = "Ca. B. cysticola",
                            "pa_min" = "P. minibicornis",
                            "arena2" = "Arenavirus",
                            "fa_mar" = "F. margolisi",
                            "my_arc" = "M. arcticus",
                            "ven" = "VENV",
                            "ic_hof" = "I. hoferi",
                            "lo_sal" = "L. salmonae",
                            "rlo" = "RLO",
                            "sp_des" = "S. destruens",
                            "ku_thy" = "K. thyrsites",
                            "prv" = "PRV",
                            "pspv" = "PSPV",
                            "ce_sha" = "C. shasta",
                            "pa_pse" = "P. pseudobranchicola",
                            "de_sal" = "D. salmonis"))+
  theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(hjust = 0.5))+
  coord_flip()
#dev.off()

## Plots per agent
#### Extract output from agent model
post_ic_mul <- post_all[post_all$X.1=="ic_mul",]

## Extract Posterior slopes by Stock
post_ic_mul_stockslp <- post_ic_mul[grep("prev_std Stock", post_ic_mul$X) ,]

#### Plot
#tiff('Fig_SSHI ONNE Pathogen Productivity_ic_mul prev by stock.tiff', 
#     units="in", width=5, height=6, res=300)
ggplot(post_ic_mul_stockslp) +
  geom_hline(yintercept = 0, linetype = "dashed", col="blue")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="black") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="black") +
  geom_point(aes(x = X, y = X50.), size = 3) +
  ylim(-0.8,0.8)+
  coord_flip()
#dev.off()

## Extract Posterior intercepts for Stocks - ic_mul example
post_ic_mul_stockint <- post_ic_mul[c(1:18) ,]
## Plot
ggplot(post_ic_mul_stockint) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()

## Extract sigma values
post_ic_mul_stocksig <- post_ic_mul[grep("igma", post_ic_mul$X) ,]
## Plot
ggplot(post_ic_mul_stocksig) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = X, ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()

## Extract Posterior intercepts for years
post_ic_mul_year <- post_ic_mul[grep("b\\[\\(\\Intercept) Year", post_ic_mul$X) ,]
## Plot
ggplot(post_ic_mul_year) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()


## Calculate stock-specific Posterior estimates by adding 
## stock-specific draws to "global" (averaged across stocks) and then averaging
param <- colnames(data.frame(mod.arena2)) #create object of parameters using any model
temp2 <- data.frame(mod.arena2)
temp3 <- temp2[,grepl("prev",names(temp2))] #include only columns with "prev" in name (posterior slopes/intercepts)
temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove sigma-related columns 
temp <- temp4[2001:4000,] #remove warm up iterations
temp.add <- data.frame(temp[,1] + temp[,2:19]) #add the stock-specific draws to global column
df <- cbind(temp4[,1], temp.add) #bind with global column
para_name2 <- colnames(temp4) #create an object with column names
colnames(df) <- colnames(temp4) #assign names to columns

## Create a matrix and loop
stk.spec.slope <- matrix(NA,
                         nrow = 19,
                         ncol = 6,
                         dimnames = list(para_name2,c("2.5","25","50","75","97.5","agent")))

for(i in agents){
  model<-get(paste("mod",i, sep="."))
  temp2 <- data.frame(model)
  temp3 <- temp2[,grepl("prev",names(temp2))] #include only columns with "prev" in name
  temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
  temp5 <- temp4[2001:4000,] #remove warm up iterations
  temp6 <- data.frame(temp5[,1] + temp5[,2:19]) #add the stock-specific draws to global column
  temp7 <- cbind(temp5[,1],temp6) #bind with global column
  para_name2 <- colnames(temp5) #create an object with column names
  colnames(temp7) <- colnames(temp5) #assign names to columns
  nam <- paste("stk.spec.slope", i, sep = ".")
  temp8 <- as.matrix(apply(temp7, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975)))
  temp9 <- t(temp8)
  temp10<- as.matrix(cbind(temp9, paste(i)))
  colnames(temp10) <- c("2.5","25","50","75","97.5","agent") #assign names to columns
  as.matrix(assign(nam, stk.spec.slope[i] <- temp10))
}

#### Rbind all posteriors and save as .csv file
stk.spec.slope.all <- rbind(stk.spec.slope.arena2,
                            stk.spec.slope.c_b_cys,
                            stk.spec.slope.ce_sha,
                            stk.spec.slope.de_sal,
                            stk.spec.slope.fa_mar,
                            stk.spec.slope.fl_psy,
                            stk.spec.slope.ic_hof,
                            stk.spec.slope.ic_mul,
                            stk.spec.slope.ku_thy,
                            stk.spec.slope.lo_sal,
                            stk.spec.slope.my_arc,
                            stk.spec.slope.pa_kab,
                            stk.spec.slope.pa_min,
                            stk.spec.slope.pa_pse,
                            stk.spec.slope.pa_ther,
                            stk.spec.slope.prv,
                            stk.spec.slope.pspv,
                            stk.spec.slope.rlo,
                            stk.spec.slope.sch,
                            stk.spec.slope.sp_des,
                            stk.spec.slope.te_bry,
                            stk.spec.slope.te_mar,
                            stk.spec.slope.ven)
write.csv(stk.spec.slope.all, file="data/Stock specific slopes_prev.csv")

##### READ IN DATA FROM FILE
stspslp <- read.csv("data/Stock specific slopes_prev.csv")

## Extract stock-specific slopes - ic_mul model
stk.spec.ic_mul <-stspslp[stspslp$agent=="ic_mul",]

ggplot(stk.spec.ic_mul) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(X, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.ic_mul[stk.spec.ic_mul$X=="Total",], aes(x = X, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.ic_mul[stk.spec.ic_mul$X=="Total",], aes(x = X, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = X, y = X50), size = 2) +
  geom_point(data=stk.spec.ic_mul[stk.spec.ic_mul$X=="Total",], aes(x = X, y = X50), size = 3) +
  labs(x="Stocks", y="Effect size") +
  coord_flip()

## Plot stock-specific slopes - te_mar
stk.spec.te_mar <-stspslp[stspslp$agent=="te_mar",]
ggplot(stk.spec.te_mar) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.te_mar[stk.spec.te_mar$stock=="Total",], aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.te_mar[stk.spec.te_mar$stock=="Total",], aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.te_mar[stk.spec.te_mar$stock=="Total",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stocks", y="Effect size") +
  coord_flip()

## Plot stock-specific slopes - pa_ther
stk.spec.pa_ther <-stspslp[stspslp$agent=="pa_ther",]
ggplot(stk.spec.pa_ther) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.pa_ther[stk.spec.pa_ther$stock=="Total",], aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.pa_ther[stk.spec.pa_ther$stock=="Total",], aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.pa_ther[stk.spec.pa_ther$stock=="Total",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stocks", y="Effect size") +
  coord_flip()


#################################################################################
#################################################################################
# LOAD - INDEPENDENT MODELS by AGENT

## Global metric - Independent models

### Loop for STAN independent models
agents <- unique(inf_agt_resid_data_gl$agent)
for(i in agents){
  data <- subset(inf_agt_resid_data_gl, agent==i)
  nam <- paste("mod.load", i, sep = ".")
  assign(nam, stan_lmer(resid_value ~ 0 + load_std + (load_std|Stock) +(1|Year), 
                        data = data,
                        adapt_delta=0.95,
                        REML = F))
}

## Loop to derive coefficient estimates
coefs_stan <- matrix(NA,
                     nrow = length(agents),
                     ncol = 5,
                     dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- summary(model, 
                      pars = c("load_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan, file="load_coefs_stan_global_indep mod.csv")

#load estimates from file
agents <- unique(inf_agt_resid_data_gl$agent)
coefs_stan <- read.csv("load_coefs_stan_global_indep mod.csv")
rownames(coefs_stan) <- coefs_stan[,1]
coefs_stan <- coefs_stan[,-1] # drop first column with agent names

#Plot effect size of agents
coefs_order <- coefs_stan[order(coefs_stan[,3]),]
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))

plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,1]),
       ui = (coefs_order[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-1.75,1),
       pch = 16,
       scol = "grey")

plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,2]),
       ui = (coefs_order[,4]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "grey")

text(rep(-1.75,length(agents)), 
     seq(1,length(agents)), 
     labels = rownames(coefs_order), 
     pos = 4,
     font = 2,
     cex=0.95)

axis(1, at = c(-1, -0.5, 0, 0.5, 1))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("effect size",1,line=2.2, cex=1.1)
mtext("Intensity",3,line=0.25)


## Loop to derive posteriors by stock
stocks<-unique(inf_agt_resid_data_gl$Stock)
### Intercepts
coefs_stan_stk_int_load <- matrix(NA,
                                  nrow = length(stocks),
                                  ncol = 6,
                                  dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\(\\Intercept) Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_int_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_int_load[i] <- cbind(ind_coef[c(1:18),c(4:8)], paste(i))))
}

### Slopes
coefs_stan_stk_slp_load <- matrix(NA,
                                  nrow = length(stocks),
                                  ncol = 6,
                                  dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\load_std Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_slp_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_slp_load[i] <- cbind(ind_coef[c(1:18),c(4:8)], paste(i))))
}

### Year intercepts
years<-unique(inf_agt_resid_data_gl$Year)
coefs_stan_stk_year_load <- matrix(NA,
                                   nrow = length(years),
                                   ncol = 6,
                                   dimnames = list(years,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = "b\\[\\(\\Intercept) Year",
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_year_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_year_load[i] <- cbind(ind_coef[c(1:8),c(4:8)], paste(i))))
}

### Sigmas
coefs_stan_stk_sig_load <- matrix(NA,
                                  nrow = length(stocks),
                                  ncol = 6,
                                  dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("sigma","Sigma"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_sig_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_sig_load[i] <- cbind(ind_coef[c(1:5),c(4:8)], paste(i))))
}

### Rhat and Neff
#extract parameter names
sims<-as.matrix(mod.load.arena2) #use any model to get parameter names, they should be the same
dim(sims)
para_name <- c(colnames(sims), "mean_PPD", "log-posterior")
para_name

coefs_stan_stk_rhat_load <- matrix(NA,
                                   nrow = length(para_name),
                                   ncol = 2,
                                   dimnames = list(para_name,c("Rhat","agent")))
for(i in agents){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  nam <- paste("coefs_stan_stk_rhat_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_rhat_load[i] <- cbind(ind_coef[c(1:52),1], paste(i), paste("Rhat"))))
}

coefs_stan_stk_neff_load <- matrix(NA,
                                   nrow = length(para_name),
                                   ncol = 2,
                                   dimnames = list(para_name,c("Neff","agent")))
for(i in agents){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  nam <- paste("coefs_stan_stk_neff_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_neff_load[i] <- cbind(ind_coef[c(1:52),1], paste(i), paste("Neff"))))
}

post_int_slpintsig_load <- rbind(coefs_stan_stk_int_load.arena2,
                                 coefs_stan_stk_int_load.c_b_cys,
                                 coefs_stan_stk_int_load.ce_sha,
                                 coefs_stan_stk_int_load.de_sal,
                                 coefs_stan_stk_int_load.fa_mar,
                                 coefs_stan_stk_int_load.fl_psy,
                                 coefs_stan_stk_int_load.ic_hof,
                                 coefs_stan_stk_int_load.ic_mul,
                                 coefs_stan_stk_int_load.ku_thy,
                                 coefs_stan_stk_int_load.lo_sal,
                                 coefs_stan_stk_int_load.my_arc,
                                 coefs_stan_stk_int_load.pa_kab,
                                 coefs_stan_stk_int_load.pa_min,
                                 coefs_stan_stk_int_load.pa_pse,
                                 coefs_stan_stk_int_load.pa_ther,
                                 coefs_stan_stk_int_load.prv,
                                 coefs_stan_stk_int_load.pspv,
                                 coefs_stan_stk_int_load.rlo,
                                 coefs_stan_stk_int_load.sch,
                                 coefs_stan_stk_int_load.smallUK,
                                 coefs_stan_stk_int_load.sp_des,
                                 coefs_stan_stk_int_load.te_bry,
                                 coefs_stan_stk_int_load.te_mar,
                                 coefs_stan_stk_int_load.ven,
                                 coefs_stan_stk_slp_load.arena2,
                                 coefs_stan_stk_slp_load.c_b_cys,
                                 coefs_stan_stk_slp_load.ce_sha,
                                 coefs_stan_stk_slp_load.de_sal,
                                 coefs_stan_stk_slp_load.fa_mar,
                                 coefs_stan_stk_slp_load.fl_psy,
                                 coefs_stan_stk_slp_load.ic_hof,
                                 coefs_stan_stk_slp_load.ic_mul,
                                 coefs_stan_stk_slp_load.ku_thy,
                                 coefs_stan_stk_slp_load.lo_sal,
                                 coefs_stan_stk_slp_load.my_arc,
                                 coefs_stan_stk_slp_load.pa_kab,
                                 coefs_stan_stk_slp_load.pa_min,
                                 coefs_stan_stk_slp_load.pa_pse,
                                 coefs_stan_stk_slp_load.pa_ther,
                                 coefs_stan_stk_slp_load.prv,
                                 coefs_stan_stk_slp_load.pspv,
                                 coefs_stan_stk_slp_load.rlo,
                                 coefs_stan_stk_slp_load.sch,
                                 coefs_stan_stk_slp_load.smallUK,
                                 coefs_stan_stk_slp_load.sp_des,
                                 coefs_stan_stk_slp_load.te_bry,
                                 coefs_stan_stk_slp_load.te_mar,
                                 coefs_stan_stk_slp_load.ven,
                                 coefs_stan_stk_sig_load.arena2,
                                 coefs_stan_stk_sig_load.c_b_cys,
                                 coefs_stan_stk_sig_load.ce_sha,
                                 coefs_stan_stk_sig_load.de_sal,
                                 coefs_stan_stk_sig_load.fa_mar,
                                 coefs_stan_stk_sig_load.fl_psy,
                                 coefs_stan_stk_sig_load.ic_hof,
                                 coefs_stan_stk_sig_load.ic_mul,
                                 coefs_stan_stk_sig_load.ku_thy,
                                 coefs_stan_stk_sig_load.lo_sal,
                                 coefs_stan_stk_sig_load.my_arc,
                                 coefs_stan_stk_sig_load.pa_kab,
                                 coefs_stan_stk_sig_load.pa_min,
                                 coefs_stan_stk_sig_load.pa_pse,
                                 coefs_stan_stk_sig_load.pa_ther,
                                 coefs_stan_stk_sig_load.prv,
                                 coefs_stan_stk_sig_load.pspv,
                                 coefs_stan_stk_sig_load.rlo,
                                 coefs_stan_stk_sig_load.sch,
                                 coefs_stan_stk_sig_load.smallUK,
                                 coefs_stan_stk_sig_load.sp_des,
                                 coefs_stan_stk_sig_load.te_bry,
                                 coefs_stan_stk_sig_load.te_mar,
                                 coefs_stan_stk_sig_load.ven,
                                 coefs_stan_stk_year_load.arena2,
                                 coefs_stan_stk_year_load.c_b_cys,
                                 coefs_stan_stk_year_load.ce_sha,
                                 coefs_stan_stk_year_load.de_sal,
                                 coefs_stan_stk_year_load.fa_mar,
                                 coefs_stan_stk_year_load.fl_psy,
                                 coefs_stan_stk_year_load.ic_hof,
                                 coefs_stan_stk_year_load.ic_mul,
                                 coefs_stan_stk_year_load.ku_thy,
                                 coefs_stan_stk_year_load.lo_sal,
                                 coefs_stan_stk_year_load.my_arc,
                                 coefs_stan_stk_year_load.pa_kab,
                                 coefs_stan_stk_year_load.pa_min,
                                 coefs_stan_stk_year_load.pa_pse,
                                 coefs_stan_stk_year_load.pa_ther,
                                 coefs_stan_stk_year_load.prv,
                                 coefs_stan_stk_year_load.pspv,
                                 coefs_stan_stk_year_load.rlo,
                                 coefs_stan_stk_year_load.sch,
                                 coefs_stan_stk_year_load.sp_des,
                                 coefs_stan_stk_year_load.te_bry,
                                 coefs_stan_stk_year_load.te_mar,
                                 coefs_stan_stk_year_load.ven)
write.csv(post_int_slpintsig_load, file="Posterior distributions_Int Slp Sig_global_indep mod_load.csv")

#### Rbind all convergence parameters and write as csv
post_rhatneff_load <- rbind(coefs_stan_stk_rhat_load.arena2,
                            coefs_stan_stk_rhat_load.c_b_cys,
                            coefs_stan_stk_rhat_load.ce_sha,
                            coefs_stan_stk_rhat_load.de_sal,
                            coefs_stan_stk_rhat_load.fa_mar,
                            coefs_stan_stk_rhat_load.fl_psy,
                            coefs_stan_stk_rhat_load.ic_hof,
                            coefs_stan_stk_rhat_load.ic_mul,
                            coefs_stan_stk_rhat_load.ku_thy,
                            coefs_stan_stk_rhat_load.lo_sal,
                            coefs_stan_stk_rhat_load.my_arc,
                            coefs_stan_stk_rhat_load.pa_kab,
                            coefs_stan_stk_rhat_load.pa_min,
                            coefs_stan_stk_rhat_load.pa_pse,
                            coefs_stan_stk_rhat_load.pa_ther,
                            coefs_stan_stk_rhat_load.prv,
                            coefs_stan_stk_rhat_load.pspv,
                            coefs_stan_stk_rhat_load.rlo,
                            coefs_stan_stk_rhat_load.sch,
                            coefs_stan_stk_rhat_load.smallUK,
                            coefs_stan_stk_rhat_load.sp_des,
                            coefs_stan_stk_rhat_load.te_bry,
                            coefs_stan_stk_rhat_load.te_mar,
                            coefs_stan_stk_rhat_load.ven,
                            coefs_stan_stk_neff_load.arena2,
                            coefs_stan_stk_neff_load.c_b_cys,
                            coefs_stan_stk_neff_load.ce_sha,
                            coefs_stan_stk_neff_load.de_sal,
                            coefs_stan_stk_neff_load.fa_mar,
                            coefs_stan_stk_neff_load.fl_psy,
                            coefs_stan_stk_neff_load.ic_hof,
                            coefs_stan_stk_neff_load.ic_mul,
                            coefs_stan_stk_neff_load.ku_thy,
                            coefs_stan_stk_neff_load.lo_sal,
                            coefs_stan_stk_neff_load.my_arc,
                            coefs_stan_stk_neff_load.pa_kab,
                            coefs_stan_stk_neff_load.pa_min,
                            coefs_stan_stk_neff_load.pa_pse,
                            coefs_stan_stk_neff_load.pa_ther,
                            coefs_stan_stk_neff_load.prv,
                            coefs_stan_stk_neff_load.pspv,
                            coefs_stan_stk_neff_load.rlo,
                            coefs_stan_stk_neff_load.sch,
                            coefs_stan_stk_neff_load.sp_des,
                            coefs_stan_stk_neff_load.te_bry,
                            coefs_stan_stk_neff_load.te_mar,
                            coefs_stan_stk_neff_load.ven)
write.csv(post_rhatneff_load, file="Posterior distributions_Rhat and Neff_global_indep mod_load.csv")

#Plot posteriors
posterior <- as.matrix(mod.load.te_mar)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_intervals(posterior,
               prob = 0.8) + plot_title
summary(mod.load.arena2)[,"Rhat"]



## Plot from file
post_all_load <- read.csv("Posterior distributions_Int Slp Sig_global_indep mod_load.csv")
post_all_load <- droplevels(post_all_load[!post_all_load$X.1 == "smallUK",])
post_agents_load <- read.csv("load_coefs_stan_global_indep mod.csv")
post_agents_load <- droplevels(post_agents_load[!post_agents_load$X == "smallUK",])

## Plot Posterior for all agents
ggplot(post_agents_load) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(X, -mid), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = upper, ymin = lower), col="gray") +
  geom_point(aes(x = X, y = mid), size = 2)+
  coord_flip()

tiff('Fig_SSHI ONNE Pathogen Productivity_Agent slopes_Load.tiff', 
     units="in", width=5, height=6, res=300)
ggplot(post_agents_load) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(X, -mid), ymax = X75, ymin = X25), size=1.5, col="darkblue") +
  geom_linerange(aes(x = X, ymax = upper, ymin = lower), col="darkblue") +
  geom_point(aes(x = X, y = mid), size = 3, col="darkblue")+
  labs(x ="Infectious agents", y = "Global effect size", title="Intensity")+
  scale_x_discrete(labels=c("ic_mul" = "I. multifiliis", 
                            "te_mar" = "T. maritinum",
                            "pa_ther" = "P. theridion",
                            "fl_psy" = "F. psychrophilum",
                            "sch" = "Ca. S. salmonis",
                            "te_bry" = "T. bryosalmonae",
                            "pa_kab" = "P. kabatai",
                            "c_b_cys" = "Ca. B. cysticola",
                            "pa_min" = "P. minibicornis",
                            "arena2" = "Arenavirus",
                            "fa_mar" = "F. margolisi",
                            "my_arc" = "M. arcticus",
                            "ven" = "VENV",
                            "ic_hof" = "I. hoferi",
                            "lo_sal" = "L. salmonae",
                            "rlo" = "RLO",
                            "sp_des" = "S. destruens",
                            "ku_thy" = "K. thyrsites",
                            "prv" = "PRV",
                            "pspv" = "PSPV",
                            "ce_sha" = "C. shasta",
                            "pa_pse" = "P. pseudobranchicola",
                            "de_sal" = "D. salmonis"))+
  theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(hjust = 0.5))+
  coord_flip()
dev.off()

# Plots for smallUK
post_smallUK_load <- post_all_load[post_all_load$X.1=="smallUK",]
## Plot Posterior slopes for Stocks
post_smallUK_load_stockslp <- post_smallUK_load[grep("load_std Stock", post_smallUK_load$X) ,]
ggplot(post_smallUK_load_stockslp) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()
## Plot Posterior intercepts for Stocks
post_smallUK_load_stockint <- post_smallUK_load[c(1:18) ,]
ggplot(post_smallUK_load_stockint) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()
## Plot sigma values
post_smallUK_load_stocksig <- post_smallUK_load[grep("igma", post_smallUK_load$X) ,]
ggplot(post_smallUK_load_stocksig) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = X, ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()
## Plot Posterior intercepts for year
post_smallUK_load_year <- post_smallUK_load[grep("b\\[\\(\\Intercept) Year", post_smallUK_load$X) ,]
ggplot(post_smallUK_load_year) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()

# Plots for te_mar
post_te_mar_load <- post_all_load[post_all_load$X.1=="te_mar",]
## Plot Posterior slopes for Stocks
post_te_mar_load_stockslp <- post_te_mar_load[grep("load_std Stock", post_te_mar_load$X) ,]
tiff('Fig_SSHI ONNE Pathogen Productivity_te_mar load by stock.tiff', 
     units="in", width=5, height=6, res=300)
ggplot(post_te_mar_load_stockslp) +
  geom_hline(yintercept = 0, linetype = "dashed", col="blue")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="black") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="black") +
  geom_point(aes(x = X, y = X50.), size = 3) +
  ylim(-0.8,0.8)+
  coord_flip()
dev.off()

## Plot Posterior intercepts for Stocks
post_te_mar_load_stockint <- post_te_mar_load[c(1:18) ,]
ggplot(post_te_mar_load_stockint) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()
## Plot sigma values
post_te_mar_load_stocksig <- post_te_mar_load[grep("igma", post_te_mar_load$X) ,]
ggplot(post_te_mar_load_stocksig) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = X, ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()
## Plot Posterior intercepts for year
post_te_mar_load_year <- post_te_mar_load[grep("b\\[\\(\\Intercept) Year", post_te_mar_load$X) ,]
ggplot(post_te_mar_load_year) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()




#################################################################################
## Loop to derive proportions of posterior draws <0 per model - prevalence
param<-colnames(data.frame(mod.arena2)) #create object of parameters in model

param.prop0 <- matrix(NA,
                      ncol = length(param),
                      nrow = length(agents),
                      dimnames = list(agents,param))

for (i in agents){
  model<-as.matrix(get(paste("mod.",i, sep="")))
  model2<-model[2001:4000,]
  param.prop0[i,] <- (colSums(model2 < 0))/2000
}

write.csv(param.prop0, file="Percent of posterior draws less than 0_prev.csv")


## Loop to derive proportions of posterior draws <0 per model - load
param.load<-colnames(data.frame(mod.load.arena2)) #create object of parameters in model

param.prop0.load <- matrix(NA,
                           ncol = length(param.load),
                           nrow = length(agents),
                           dimnames = list(agents,param.load))

for (i in agents){
  model<-as.matrix(get(paste("mod.load.",i, sep="")))
  model2<-model[2001:4000,]
  param.prop0.load[i,] <- (colSums(model2 < 0))/2000
}

write.csv(param.prop0.load, file="Percent of posterior draws less than 0_load.csv")

#################################################################################

#################################################################################
## Loop to calculate stock-specific agent estimates - LOAD
param <- colnames(data.frame(mod.load.arena2)) #create object of parameters in model
temp2 <- data.frame(mod.load.arena2)
temp3 <- temp2[,grepl("load",names(temp2))] #include only columns with "prev" in name
temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
temp <- temp4[2001:4000,] #remove warm up iterations
temp.add <- data.frame(temp[,1] + temp[,2:19]) #add the stock-specific draws to global column
df <- cbind(temp4[,1], temp.add) #bind with global column
para_name2 <- colnames(temp4) #create an object with column names
colnames(df) <- colnames(temp4) #assign names to columns

## Create a matrix and loop

stk.spec.slope.load <- matrix(NA,
                              nrow = 19,
                              ncol = 6,
                              dimnames = list(para_name2,c("2.5","25","50","75","97.5","agent")))

for(i in agents){
  model<-get(paste("mod",i, sep="."))
  temp2 <- data.frame(model)
  temp3 <- temp2[,grepl("prev",names(temp2))] #include only columns with "prev" in name
  temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
  temp5 <- temp4[2001:4000,] #remove warm up iterations
  temp6 <- data.frame(temp5[,1] + temp5[,2:19]) #add the stock-specific draws to global column
  temp7 <- cbind(temp5[,1],temp6) #bind with global column
  para_name2 <- colnames(temp5) #create an object with column names
  colnames(temp7) <- colnames(temp5) #assign names to columns
  nam <- paste("stk.spec.slope.load", i, sep = ".")
  temp8 <- as.matrix(apply(temp7, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975)))
  temp9 <- t(temp8)
  temp10<- as.matrix(cbind(temp9, paste(i)))
  colnames(temp10) <- c("2.5","25","50","75","97.5","agent") #assign names to columns
  as.matrix(assign(nam, stk.spec.slope.load[i] <- temp10))
}

#### Rbind all posteriors and write as csv
stk.spec.slope.all.load <- rbind(stk.spec.slope.load.arena2,
                                 stk.spec.slope.load.c_b_cys,
                                 stk.spec.slope.load.ce_sha,
                                 stk.spec.slope.load.de_sal,
                                 stk.spec.slope.load.fa_mar,
                                 stk.spec.slope.load.fl_psy,
                                 stk.spec.slope.load.ic_hof,
                                 stk.spec.slope.load.ic_mul,
                                 stk.spec.slope.load.ku_thy,
                                 stk.spec.slope.load.lo_sal,
                                 stk.spec.slope.load.my_arc,
                                 stk.spec.slope.load.pa_kab,
                                 stk.spec.slope.load.pa_min,
                                 stk.spec.slope.load.pa_pse,
                                 stk.spec.slope.load.pa_ther,
                                 stk.spec.slope.load.prv,
                                 stk.spec.slope.load.pspv,
                                 stk.spec.slope.load.rlo,
                                 stk.spec.slope.load.sch,
                                 stk.spec.slope.load.sp_des,
                                 stk.spec.slope.load.te_bry,
                                 stk.spec.slope.load.te_mar,
                                 stk.spec.slope.load.ven)

write.csv(stk.spec.slope.all.load, file="Stock specific slopes_load.csv")

##### READ IN DATA FROM FILE
stspslp.load<-read.csv("Stock specific slopes_load.csv")

## Plot stock-specific slopes - te_mar
stk.spec.te_mar.load <-stspslp.load[stspslp.load$agent=="te_mar",]
ggplot(stk.spec.te_mar.load) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.te_mar.load[stk.spec.te_mar.load$stock=="Total",], aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.te_mar.load[stk.spec.te_mar.load$stock=="Total",], aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.te_mar.load[stk.spec.te_mar.load$stock=="Total",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stocks", y="Effect size") +
  coord_flip()
