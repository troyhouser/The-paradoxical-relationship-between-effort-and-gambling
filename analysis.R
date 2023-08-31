setwd("~/Dropbox (University of Oregon)/effort_data/")
source("load_prep_data.R")
data = load_data()

############################################################
############################################################
############################################################

### RUN MODELS

############################################################
############################################################
############################################################
source("models.R")

## lists for model results
rh0=rexpo=p0=hboid=h0=qh=gh=expo=
  rgh=rhboid=rp0=rqh=list()
## vector for number of trials per subject
ntrials = c()
bics = matrix(0,length(S),12) # matrix for bic scores
new_df = rep(NA,6) # dataframe for trial-unique discounts

for(s in 1:length(S)){
  ix = data[data$Participant.Private.ID==S[s],]
  ix = ix[!is.na(ix$effortChoice),]
  ix$N = nrow(ix)
  ntrials[s] = nrow(ix)
  versions[s] = ix$version[1]
  effort_selfreport[s] = ix$effortSelfReport[1]
  happy_selfreport[s] = ix$happySelfReport[1]
  if(versions[s]==3) ix$X1 = -ix$X1
  for(n in 1:nrow(ix)){
    ix$s2_u[n] = max(1e-5,var(c(ix$X1[1:n],ix$X2[1:n]),na.rm=T))
  }
  dat = ix
  rh0[[s]] = optim(par = rep(0.1,2), fn = RH0, lower = c(1e-6,1e-6), upper = c(50,100), 
                  method = "L-BFGS-B")
  h0[[s]] = optim(par = rep(0.1,2), fn = H0, lower = c(1e-6,1e-6), upper = c(100,50), 
                  method = "L-BFGS-B")
  qh[[s]] = optim(par = rep(0.1,3), fn = QH, lower = c(0,0,1e-6), upper = c(1,1,100), 
                  method = "L-BFGS-B")
  p0[[s]] = optim(par = rep(0,3), fn = P0, lower = c(0,0,0), upper = c(100,50,50), 
                   method = "L-BFGS-B")
  hboid[[s]] = optim(par = rep(0,3), fn = HBoid, lower = c(0,0,0), upper = c(100,50,50), 
                  method = "L-BFGS-B")
  gh[[s]] = optim(par = rep(0.1,3), fn = GH, lower = c(1e-6,1e-6,1e-6), upper = c(50,100,50), 
                 method = "L-BFGS-B")
  expo[[s]] = optim(par = rep(0.1,2), fn = EXPO, lower = c(1e-6,1e-6), upper = c(100,50), 
                  method = "L-BFGS-B")
  rexpo[[s]] = optim(par = rep(0.1,2), fn = REXPO, lower = c(1e-6,1e-6), upper = c(50,100), 
                  method = "L-BFGS-B")
  rgh[[s]] = optim(par = rep(0.1,3), fn = RGH, lower = c(1e-6,1e-6,1e-6), upper = c(50,100,50), 
                  method = "L-BFGS-B")
  rhboid[[s]] = optim(par = rep(0.1,3), fn = RHBoid, lower = rep(1e-6,3), upper = c(100,50,1), 
                     method = "L-BFGS-B")
  rp0[[s]] = optim(par = rep(0.1,3), fn = RP0, lower = rep(1e-6), upper = c(50,100,1), 
                  method = "L-BFGS-B")
  
  bics[s,1] = 2*rh0[[s]]$value+log(ntrials[s])*2
  bics[s,2] = 2*rexpo[[s]]$value+log(ntrials[s])*2
  bics[s,3] = 2*h0[[s]]$value+log(ntrials[s])*2
  bics[s,4] = 2*qh[[s]]$value+log(ntrials[s])*3
  bics[s,5] = 2*rand()+log(ntrials[s])
  bics[s,6] = 2*p0[[s]]$value+log(ntrials[s])*3
  bics[s,7] = 2*hboid[[s]]$value+log(ntrials[s])*3
  bics[s,8] = 2*gh[[s]]$value+log(ntrials[s])*3
  bics[s,9] = 2*expo[[s]]$value+log(ntrials[s])*2
  bics[s,10] = 2*rhboid[[s]]$value+log(ntrials[s])*3
  bics[s,11] = 2*rgh[[s]]$value+log(ntrials[s])*3
  bics[s,12] = 2*rp0[[s]]$value+log(ntrials[s])*3
  
  parsREXPO = c(rexpo[[s]]$par[1],rexpo[[s]]$par[2])
  rh0_k = get_k_for_RH0(parsREXPO[1],parsREXPO[2])
  new_df = rbind(new_df,cbind(ix$Participant.Private.ID,ix$risk,ix$LL,
                              ix$version,ix$effortSelfReport,
                              rh0_k))
  
}

# mean & sd of each model
colMeans(bics);apply(bics,2,sd)

## test for sig difference between #1 & #2 and #1 & #3 models
t.test(bics[,2],bics[,9],paired=T);lsr::cohensD(bics[,2],bics[,9])
t.test(bics[,2],bics[,3],paired=T);lsr::cohensD(bics[,2],bics[,3])

# get pxp scores
library(qpcR)
bic_weights = akaike.weights(bics)
bic_weights = matrix(bic_weights$weights,nrow(bics),12,byrow=F)
pxp = bmsR::VB_bms(log(bic_weights))
############################################################
############################################################

### plot model fits

############################################################
############################################################

bic_df = data.frame(BIC = c(bics),
                    model = rep(c("BH","BExp","H","QH",
                                  "Rand","P","HD","GH",
                                  "Exp","BHD","BGH","BP"),each=nrow(bics)),
                    sub = rep(S,12),version = rep(factor(versions),12))
bic_df = bic_df[!bic_df$model=="QH",]
colors = c("#4133FF","#3348FF","#3389FF","#33ACFF","#33D2FF",
           "white",
           "#FFD933","#FFCB33","#FFC133","#FFAF33","#FF9B33")
bic_df$model = factor(bic_df$model,levels=c("H","GH","Exp","HD","P",
                                            "Rand",
                                            "BH","BGH","BExp","BHD","BP"))
p1 = ggplot(bic_df,aes(x=model,y=BIC,fill=model))+
  #geom_line(aes(group=sub),alpha=0.1)+
  geom_bar(stat="summary",position="dodge",col="black")+
  geom_jitter(width=0.2,shape=21)+
  scale_y_continuous(limits = c(0,56),expand= c(0,0))+
  scale_fill_manual(values=colors)+
  geom_errorbar(stat="summary",position=position_dodge(0.9),width=0.2,col="black",linewidth=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=15))+
  ylab("BIC")+xlab("model")+
  theme(legend.position = "None")+
  theme(text=element_text(size=20))
p1
ggsave(plot=p1,filename="model_fits.png",units="px",width=2000,height=1500,
       path="~/Dropbox (University of Oregon)/effort_data/plots/")

pxp_df = data.frame(PXP = pxp$pxp,
                    model = c("BH","BExp","H","QH",
                              "Rand","P","HD","GH",
                              "Exp","BHD","BGH","BP"))
pxp_df = pxp_df[!pxp_df$model=="QH",]
pxp_df$model = factor(pxp_df$model,levels=c("H","GH","Exp","HD","P",
                                            "Rand",
                                            "BH","BGH","BExp","BHD","BP"))
p2 = ggplot(pxp_df,aes(x=PXP,y=model,fill=model))+
  geom_bar(stat="identity",col="black")+
  scale_x_continuous(limits = c(0,1.009),expand= c(0,0))+
  scale_fill_manual(values=colors)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=15))+
  ylab("model")+xlab("PXP")+
  theme(legend.position = "None")+
  theme(text=element_text(size=20))
p2  
ggsave(plot=p2,filename="pxp.png",units="px",width=1200,height=1200,
       path="~/Dropbox (University of Oregon)/effort_data/plots/")

############################################################
############################################################

#### parameter recovery

############################################################
############################################################

rexpo_beta = unlist(lapply(rexpo, function(x) x$par[1]))
rexpo_alpha = unlist(lapply(rexpo, function(x) x$par[2]))
sims = 500
sim_data = function(beta1,alpha1){
  props = matrix(0,length(S),sims)
  for(s in 1:length(S)){
    beta = beta1[s]; alpha = alpha1[s]
    ix = data[data$Participant.Private.ID==S[s],]
    ix = ix[!is.na(ix$effortChoice),]
    ix$N = nrow(ix)
    ntrials[s] = nrow(ix)
    versions[s] = ix$version[1]
    if(versions[s]==3) ix$X1 = -ix$X1
    for(n in 1:nrow(ix)){
      ix$s2_u[n] = max(1e-5,var(c(ix$X1[1:n],ix$X2[1:n]),na.rm=T))
    }
    dat = ix
    r = cbind(dat$X2,dat$X1)
    r[r==0]=1e-6
    t = cbind(dat$T2,dat$T1)
    s2_e = dat$s2_u/(abs(r)*beta)
    k = s2_e/dat$s2_u
    V = r*exp(-k*t)
    P = pnorm(alpha*(V[,1]-V[,2]))
    
    for(k in 1:sims){
      choices = c()
      for(i in 1:length(P)){
        choices[i] = sample(0:1,1,prob=c(1-P[i],P[i]))
      }
      props[s,k] = mean(choices,na.rm=T)
    }
  }
  return(props)
}
actualLL = aggregate(LL~Participant.Private.ID,data,mean)
simLL = sim_data(rexpo_beta,rexpo_alpha)
avgSimLL = rowMeans(simLL)
versions[versions==2]=1
parm_recov = data.frame(actualLL = actualLL$LL,
                        simLL = avgSimLL,
                        sub = actualLL$Participant.Private.ID,
                        framing = factor(versions))

p0 = ggplot(parm_recov,aes(x=simLL,y=actualLL,fill=framing,col=framing))+
  geom_smooth(method="lm")+
  geom_point(shape=21,size=2,col="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=15))+
  ylab("actual proportion of larger effort choices")+xlab("predicted proportion of larger effort choices")+
  theme(legend.position = c(0.75,0.2))+
  scale_fill_manual(values=c("purple","#ED5BB6"),labels=c("gains","losses"))+
  scale_color_manual(values=c("purple","#ED5BB6"),labels=c("gains","losses"))+
  theme(text=element_text(size=20))
p0
ggsave(plot=p0,filename="parm_recov.png",units="px",width=1400,height=1400,
       path="~/Dropbox (University of Oregon)/effort_data/plots/")

cor.test(parm_recov$actualLL,parm_recov$simLL)
cor.test(parm_recov$actualLL[parm_recov$framing==1],parm_recov$simLL[parm_recov$framing==1])
cor.test(parm_recov$actualLL[parm_recov$framing==3],parm_recov$simLL[parm_recov$framing==3])

#####################
############################################################
############################################################

###### stats and plots for simpson's paradox

############################################################
############################################################
############################################################

###########################################################################
#mixed model
mm1a = lmerTest::lmer(risk~LL*factor(version)+(1+LL|Participant.Private.ID),data,REML=F)
mm1aCs = summary(mm1a)#negative within subject effect
# bootstrap the CIs
mm1boot <- bootMer(mm1a, FUN_bootMer,
                   nsim = 1000, type = "parametric", .progress =
                     "txt", PBargs = list(style = 3))
mixedCI=confint(mm1boot)

#simple model
lm1 = lm(risk~HE*versions,df)#positive group level effect
lm1Cs=summary(lm1)

#plot effects
plot_models = data.frame(coefs = c(mm1aCs$coefficients[,1],lm1Cs$coefficients[,1]),
                         LB = c(mixedCI[,1],lm1Cs$coefficients[,1]-lm1Cs$coefficients[,2]),
                         UB = c(mixedCI[,2],lm1Cs$coefficients[,1]+lm1Cs$coefficients[,2]),
                         predictor = rep(c("intercept","high-effort",
                                       "framing","effort x framing"),2),
                         model = rep(c("mixed","simple"),each=4))
plot_models = plot_models[!plot_models$predictor=="intercept",]
p100=ggplot(plot_models,aes(x=coefs,y=predictor,fill=model))+
  geom_bar(stat="identity",position="dodge",col="black")+
  geom_errorbar(aes(x=coefs,xmin=LB,xmax=UB,y=predictor),
                position=position_dodge(0.9),width=0.15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("predictor")+xlab("coefficient")+
  theme(legend.position = c(0.75,0.5))+
  scale_fill_manual(values=c("gray","white"))+
  theme(text=element_text(size=20))

ggsave(plot=p100,filename="paradox.png",units="px",width=2000,height=1400,
       path="~/Dropbox (University of Oregon)/effort_data/plots/")

### get aggregated variables of interest
LE = aggregate(LL~Participant.Private.ID,data,mean)## proportion of high-effort choices
risky = aggregate(risk~Participant.Private.ID,data,mean)## proportion of gambles
sr = aggregate(effortSelfReport~Participant.Private.ID,data,mean)## mean self-reports
new_df = new_df[-1,]
new_df = as.data.frame(new_df)
colnames(new_df) = c("sub","risk","LE","version","sr","REXPO_k")
new_df$REXPO_k = log(new_df$REXPO_k)## convert discount factors to logarithm
kagg = aggregate(REXPO_k~version+sub,new_df,mean)## mean log(k)
colnames(kagg)[3] = "k"## change column name to k
for(i in 1:length(S)){ ## get framing group per subject
  d = data[data$Participant.Private.ID==S[i],]
  versions[i] = d$version[1]
}
df = data.frame(HE = LE$LL,risk = risky$risk,
                beta = rexpo_beta, alpha = rexpo_alpha,
                sr = sr$effortSelfReport,
                framing = factor(versions),
                k = kagg$k,
                sub = actualLL$Participant.Private.ID)

##############################################################
##############################################################
##############################################################
##############################################################

## reanalyze simpson's paradox with discount factors

##############################################################
##############################################################
##############################################################

library(lme4)
mm3 = lmerTest::lmer(risk~REXPO_k*factor(version)+(1|sub),new_df)
mm3Cs=summary(mm3)
mm3boot <- bootMer(mm3, FUN_bootMer,
                   nsim = 1000, type = "parametric", .progress =
                     "txt", PBargs = list(style = 3))
mixedCI2=confint(mm3boot)

lm2 = lm(risk~k*framing,df)
lm2Cs=summary(lm2)

plot_models2 = data.frame(coefs = c(mm3Cs$coefficients[,1],lm2Cs$coefficients[,1]),
                         LB = c(mixedCI2[,1],lm2Cs$coefficients[,1]-lm2Cs$coefficients[,2]),
                         UB = c(mixedCI2[,2],lm2Cs$coefficients[,1]+lm2Cs$coefficients[,2]),
                         predictor = rep(c("intercept","high-effort",
                                           "framing","effort x framing"),2),
                         model = rep(c("mixed","simple"),each=4))
plot_models2 = plot_models2[!plot_models2$predictor=="intercept",]
p101=ggplot(plot_models2,aes(x=coefs,y=predictor,fill=model))+
  geom_bar(stat="identity",position="dodge",col="black")+
  geom_errorbar(aes(x=coefs,xmin=LB,xmax=UB,y=predictor),
                position=position_dodge(0.9),width=0.15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("predictor")+xlab("coefficient")+
  theme(legend.position = c(0.8,0.85))+
  scale_fill_manual(values=c("gray","white"))+
  theme(text=element_text(size=20))
p101
ggsave(plot=p101,filename="paradox_wK.png",units="px",width=2000,height=1400,
       path="~/Dropbox (University of Oregon)/effort_data/plots/")

##############################################################
##############################################################
##############################################################

##### exploratory reward sensitivity analyses

##############################################################
##############################################################
##############################################################
cor.test(df$beta,df$alpha)
cor.test(df$beta,df$risk)
cor.test(df$beta[df$framing==1],df$risk[df$framing==1])
cor.test(df$beta[df$framing==3],df$risk[df$framing==3])
cor.test(df$beta,df$sr)
cor.test(df$beta[df$framing==1],df$sr[df$framing==1])
cor.test(df$beta[df$framing==3],df$sr[df$framing==3])

p3=ggplot(df,aes(x=beta,y=risk,col=framing,fill=framing))+
  geom_smooth(method="lm")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("proportion of gambles")+xlab(expression(paste("reward sensitivity (",beta,")")))+
  theme(legend.position = c(0.2,0.25))+
  scale_fill_manual(values=c("gray","white"))+
  theme(text=element_text(size=20))+
  scale_fill_manual(values=c("purple","#ED5BB6"),labels=c("gains","losses"))+
  scale_color_manual(values=c("purple","#ED5BB6"),labels=c("gains","losses"))
ggsave(plot=p3,filename="beta_risk.png",units="px",width=1400,height=1400,
       path="~/Dropbox (University of Oregon)/effort_data/plots/")
p4=ggplot(df,aes(x=beta,y=sr,col=framing,fill=framing))+
  geom_smooth(method="lm")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("self-reported effort")+xlab(expression(paste("reward sensitivity (",beta,")")))+
  theme(legend.position = c(0.2,0.3))+
  theme(text=element_text(size=20))+
  scale_fill_manual(values=c("purple","#ED5BB6"),labels=c("gains","losses"))+
  scale_color_manual(values=c("purple","#ED5BB6"),labels=c("gains","losses"))
p4
ggsave(plot=p4,filename="beta_sr.png",units="px",width=1400,height=1400,
       path="~/Dropbox (University of Oregon)/effort_data/plots/")
