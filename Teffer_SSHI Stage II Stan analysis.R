#SSHI ANALYSIS SOCKEYE STAGE II
#Teffer
#### Load packages
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


#### Read in data and clean
#Global metric (averaged across all stocks in SW)  
inf_agt_resid_data_gl <- read.csv("global_ONNE_productivity_infection_analysis.csv")
inf_agt_resid_data_gl$mean_load_all <- inf_agt_resid_data_gl$mean_load #add col for mean_load_all
inf_agt_resid_data_gl$mean_load_all[is.na(inf_agt_resid_data_gl$mean_load_all)] <- 0 #replace NA with 0

inf_std <- plyr::ddply(inf_agt_resid_data_gl, c("agent"),function(x) {
  scaled_prev <- scale(x$prev)
  scaled_load <- scale(x$mean_load_all)
  xx <- data.frame(scaled_prev, scaled_load)
})

inf_agt_resid_data_gl$prev_std <- inf_std[,2]
inf_agt_resid_data_gl$load_std <- inf_std[,3]
agents <- unique(inf_agt_resid_data_gl$agent)
inf_agt_resid_data_gl$Stock <- inf_agt_resid_data_gl$Stock_Analysis
head(inf_agt_resid_data_gl)


#plot raw data - Prevalence
ggplot(inf_agt_resid_data_gl,aes(prev, resid_value, color=Stock, shape=factor(Year)))+
  geom_smooth(aes(prev, resid_value, group=Stock), method = "lm", se=F, size=.2)+
  geom_point()+
  scale_shape_manual(values=1:nlevels(factor(inf_agt_resid_data_gl$Year))) +
  facet_wrap(~ agent,ncol=4, scales = "free")+
  xlab("prevalence")+
  ylab("residual")+
  theme_bw()

#plot raw data - Load
ggplot(inf_agt_resid_data_gl,aes(log10(mean_load), resid_value, color=Stock, shape=factor(Year)))+
  geom_smooth(aes(log10(mean_load), resid_value, group=Stock), method = "lm", se=F, size=.2)+
  geom_point()+
  scale_shape_manual(values=1:nlevels(factor(inf_agt_resid_data_gl$Year))) +
  facet_wrap(~ agent,ncol=4, scales="free")+
  xlab("log10 load")+
  ylab("residual")+
  theme_bw()


### STAN Approach for multi-level modeling
#### Global SW infection metric as full model (includes all stocks and agents)
#### Prevalence

model_ind_stan_all_gl <- stan_lmer(resid_value ~ 0 + agent:prev_std + (agent:prev_std|Stock) + (1|Year), 
                                   data = inf_agt_resid_data_gl,
                                   REML = F)

summary(model_ind_stan_all_gl)

#extract parameter names
sims<-as.matrix(model_ind_stan_all_gl)
dim(sims)
para_name <- colnames(sims)
para_name

#### Extract coefficients

agents <- unique(inf_agt_resid_data_gl$agent)
coefs_stan_gl <- matrix(NA,
                        nrow = length(agents),
                        ncol = 5,
                        dimnames = list(agents,c("lower","25","mid","75","upper")))
ind_coef_all_gl <- summary(model_ind_stan_all_gl, 
                           pars = c("agentarena2:prev_std",
                                    "agentc_b_cys:prev_std",
                                    "agentce_sha:prev_std",
                                    "agentde_sal:prev_std",
                                    "agentfa_mar:prev_std",
                                    "agentfl_psy:prev_std",
                                    "agentic_hof:prev_std",
                                    "agentic_mul:prev_std",
                                    "agentku_thy:prev_std",
                                    "agentlo_sal:prev_std",
                                    "agentmy_arc:prev_std",
                                    "agentpa_kab:prev_std",
                                    "agentpa_min:prev_std",
                                    "agentpa_pse:prev_std",
                                    "agentpa_ther:prev_std",
                                    "agentprv:prev_std",
                                    "agentpspv:prev_std",
                                    "agentrlo:prev_std",
                                    "agentsch:prev_std",
                                    "agentsmallUK:prev_std",
                                    "agentsp_des:prev_std",
                                    "agentte_bry:prev_std",
                                    "agentte_mar:prev_std",
                                    "agentven:prev_std"),
                           probs = c(0.025,0.25,0.5,0.75, 0.975),
                           digits = 2)
coefs_stan_gl <- ind_coef_all_gl[,c(4:8)]
rownames(coefs_stan_gl) <- agents
write.csv(coefs_stan_gl, "prev_coefs_stan_full model_global_sw.csv")


#### Plot effect sizes

#load parameter estimates
coefs_stan <- read.csv("prev_coefs_stan_full model_global_sw.csv")
rownames(coefs_stan) <- coefs_stan[,1]
coefs_stan <- coefs_stan[,-1] # drop first column with agent names
colnames(coefs_stan)<-c("lower","25","mid","75","upper")
coefs_stan

#Plot
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
mtext("Prevalence",3,line=0.25)