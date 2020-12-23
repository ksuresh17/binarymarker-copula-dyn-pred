###############################################
#Simulation for binary biomarker (Scenario 1a)
###############################################
#Produces predicted probabilties (for Bias and Variance)
#Generate data under inpsection time marker measurement 
#Markov model with a single baseline covariate
#Structure data as a longitudinal data set

rm(list=ls())

#Set working directory to where "BinaryPredictionFunctions.R" is located 
#setwd()
source("BinaryPredictionFunctions.R")

#Load in required packages 
library(SmoothHazard,quietly=TRUE)
library(lava,quietly=TRUE)
library(prodlim,quietly=TRUE)
library(dynpred,quietly=TRUE)
library(reshape2,quietly=TRUE)
library(pbivnorm,quietly=TRUE)
library(mstate,quietly=TRUE)
library(splines,quietly=TRUE)
library(timeROC,quietly = TRUE)
library(pec,quietly=TRUE)


#Specify inputs 
n<-1000                  #number of patients to simulate
Xval1<-1                #baseline covaraite value 
thor<-3                    #prediction window for creating landmrk data set 
w_predict<-thor            #prediction window to compute prediction probabilities
LMs <- seq(0,7,by=1)  #landmark times 
LMx<-LMs
insp.rate<-0.5

#Data Generation 
#Parameters to change
msmt<-"IT"
ds<-"Long"
scale.illtime=15      #Weibull scale for latent illness time (0->1 transition)
shape.illtime=1.15    #Weibull shape for latent illness time (0->1 transition)
scale.lifetime=12.5   #Weibull scale for latent life time (0->2 transition)
shape.lifetime=1.15   #Weibull shape for latent life time (0->2 transition)
scale.waittime=10     #Weibull scale for latent wait time (1->2 transition)
shape.waittime=1.15   #Weibull shape for latent wait time (1->2 transition)

scale.censtime=80
shape.censtime=1
cens_horiz<-15

#Baseline covariate X1: Stronger exposure effect on death in diseased subjects
beta_X1_01=1
beta_X1_02=0.5
beta_X1_12=2

seed<-17 #Change for each simulation 

set.seed(seed)
set <- sort(sample(unique(1:n), 500))

#Simulate baseline covariate X1 with 40% prevalence
X1_sim<-rbinom(n,size=1,prob=.5)

#Create dataset
data<-data.frame(latent.illtime=rep(NA,n))

data$latent.illtime<-rweibull(n=n,shape=shape.illtime,scale=scale.illtime)/((exp(beta_X1_01*X1_sim))^(1/shape.illtime))
data$latent.lifetime<-rweibull(n=n,shape=shape.lifetime,scale=scale.lifetime)/((exp(beta_X1_02*X1_sim))^(1/shape.lifetime))
data$latent.waittime<-(-log(runif(n))*scale.waittime^shape.waittime/exp(beta_X1_12*X1_sim)+data$latent.illtime^(shape.waittime))^(1/shape.waittime)

#administrative censoring at 20 years 
data$censtime<-runif(n=n,0,cens_horiz)
data$censtime<-ifelse(data$censtime>=cens_horiz,cens_horiz,data$censtime)

#Simulate inspection times 
n.inspections=30
illness.known.at.death<-0

if (n.inspections>0)
  for (k in 1:n.inspections){
    data[paste("inspection", k, sep="")] <- NA
  }
ipos <- grep("inspection[0-9]+",names(data))
for(i in 1:nrow(data))
{
  data[i,ipos]<-cumsum(rexp(n=n.inspections,insp.rate))
}

# construct illtime and true illness status
data$illtime <- data$latent.illtime
data$illstatus <- 1*(data$illtime<=data$latent.lifetime)
data$illtime[data$illtime>data$latent.lifetime] <- 0

# construct lifetime
data$lifetime <- data$latent.lifetime
data$lifetime[data$illstatus==1] <- data$latent.waittime[data$illstatus==1]

#see simulateIDM in SmoothHazard package
iframe<-data[,ipos]
interval<-data.frame("L"=NA,"R"=NA,"seen.ill"=NA)
for(i in 1:n)
{
  itimes<-iframe[i,]
  ## remove inspection times that are larger than the individual lifetime
  itimes <- itimes[itimes<data$lifetime[i]]
  ## and those larger than the right censoring time
  itimes <- itimes[itimes<data$censtime[i]]
  ## if all inspection times are censored set a single one at 0
  if (length(itimes)==0) itimes <- 0
  ## mark the last inspection time 
  last.inspection <- itimes[length(itimes)]
  ## find the interval where illness happens
  if (data$illstatus[i]==1){
    ## subject was ill
    if (data$illtime[i] > last.inspection){
      ## no illness observed at last inspection
      if (data$censtime[i]<data$lifetime[i]){
        ## right censored: no illness observed
        ret<- c(last.inspection,data$censtime[i],0)
      }else{
        ## user option decides if illness is recorded at death
        ret<-  c(last.inspection,data$lifetime[i],illness.known.at.death)
      }
    }else{ ## illtime is smaller or equal to last inspection time
      if (length(itimes)==1){
        ret<-c(0,itimes,1)
      } else{
        hit <- prodlim::sindex(eval.times=data$illtime[i],
                               jump.times=itimes,
                               strict=TRUE)
        ret<- c(c(0,itimes)[c(1+hit,2+hit)],1)
      }
    }
  } else {
    ## subject was never ill
    if (data$censtime[i]<data$lifetime[i]){
      ## right censored: no illness observed until last inspection
      ret<- c(last.inspection,data$censtime[i],0)
    } else{
      ## no illness observed until death
      if(illness.known.at.death==1){
        ret<-  c(data$lifetime[i],data$lifetime[i],0)
      }else{
        ## no illness observed until last inspection provide [L=last insp.; R=lifetime]
        ret<-  c(last.inspection,data$lifetime[i],0)
      }
    }
  }
  interval[i,]<-ret
}

data <- cbind(data,interval)

data$seen.exit <- 1*(data$lifetime<data$censtime)
data$observed.lifetime <- pmin(data$lifetime,data$censtime)
data$observed.illtime <- pmin(data$R,data$censtime)

data$X1<-X1_sim

data$id<-1:nrow(data)
data$id<-factor(data$id)
data$survtime<-data$observed.lifetime
data$survstatus<-data$seen.exit

long_data <- reshape(data, 
                     varying = paste("inspection", 1:n.inspections, sep=""), 
                     v.names = "inspection.time",
                     timevar = "inspection", 
                     times = c(1:n.inspections), 
                     direction = "long")
long_data<-subset(long_data,long_data$inspection.time<=long_data$survtime)
long_data$illness<-ifelse(long_data$latent.illtime>=long_data$inspection.time,0,1)

long_data_temp<-data[,which(!names(data)%in%paste("inspection", 1:n.inspections, sep=""))]
long_data_temp$inspection.time=0
long_data_temp$inspection=0
long_data_temp$illness=0

long_data<-rbind(long_data,long_data_temp)

long_data$wsurvtime<-pmin(long_data$inspection.time+thor,long_data$survtime)
long_data$survstatus<-ifelse(long_data$lifetime<=long_data$wsurvtime,1,0)

long_data<-long_data[order(long_data$id,long_data$inspection.time),]

long_data$LM<-long_data$inspection.time

#Select training set 
train_data <-subset(long_data,!long_data$id%in%set) 
test_data <-subset(long_data,long_data$id%in%set) 

train_data.id<-train_data[!duplicated(train_data$id)&train_data$inspection.time==0,]
test_data.id<-test_data[!duplicated(test_data$id),]

LMdata<-train_data
LMdata<-LMdata[order(LMdata$id,LMdata$inspection.time),]

tt<-sort(unique(LMdata$wsurvtime[LMdata$survstatus==1]))
LMdata$LM<-LMdata$inspection.time
LMdata$Tstart<-LMdata$inspection.time

LMdata2<-survSplit(Surv(Tstart,wsurvtime,survstatus)~.,data=LMdata,cut=tt,end="wsurvtime",start="Tstart",event="survstatus")

#Define covariates used in landmark models 
LMdata$illnessX<-LMdata$illness*(LMdata$wsurvtime-LMdata$inspection.time)
LMdata$illnessX2<-LMdata$illness*(LMdata$wsurvtime-LMdata$inspection.time)^2
LMdata$illness_tau<-LMdata$illness*(LMdata$inspection.time)
LMdata$illness_tau2<-LMdata$illness*(LMdata$inspection.time)^2
LMdata$X1_int<-LMdata$illness*LMdata$X1
LMdata$V<-LMdata$illness*LMdata$observed.illtime

LMdata2$illnessX<-LMdata2$illness*(LMdata2$wsurvtime-LMdata2$inspection.time)
LMdata2$illnessX2<-LMdata2$illness*(LMdata2$wsurvtime-LMdata2$inspection.time)^2
LMdata2$illness_tau<-LMdata2$illness*(LMdata2$inspection.time)
LMdata2$illness_tau2<-LMdata2$illness*(LMdata2$inspection.time)^2
LMdata2$X1_int<-LMdata2$illness*LMdata2$X1
LMdata2$V<-LMdata2$illness*LMdata2$observed.illtime

g1<-function(t) t
g2<-function(t) t^2 
LMdata2$LM1<-g1(LMdata2$LM)
LMdata2$LM2<-g2(LMdata2$LM)

#######################################################
#Landmark Models 
#######################################################
#LM2: extended super model 
LMsupercox2 <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ illness + illness_tau + illness_tau2 + LM1+LM2+ X1 +cluster(id), data=LMdata2, method="breslow")
bh_supercox2<-basehaz(LMsupercox2,centered=FALSE)

#LM3: non-proportional hazards extended super model
LMsupercox3<-coxph(Surv(Tstart,wsurvtime,survstatus)~illness+illnessX+illnessX2+LM1+LM2+X1+cluster(id),data=LMdata2,method="breslow")
bh_supercox3<-basehaz(LMsupercox3,centered=FALSE)

#LM4: non-proportional hazards extended super model with functions of tau and s
LMsupercox4<-coxph(Surv(Tstart,wsurvtime,survstatus)~illness+illnessX+illnessX2+illness_tau+illness_tau2+LM1+LM2+X1+cluster(id),data=LMdata2,method="breslow")
bh_supercox4<-basehaz(LMsupercox4,centered=FALSE)

#LM2Int: LM2 + interactions
LMsupercox2_int <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ illness + illness_tau + illness_tau2 + LM1+LM2+ X1+ +X1_int+cluster(id), data=LMdata2, method="breslow")
bh_supercox2_int<-basehaz(LMsupercox2_int,centered=FALSE)

#LM3Int: LM3 + interactions
LMsupercox3_int<-coxph(Surv(Tstart,wsurvtime,survstatus)~illness+illnessX+illnessX2+LM1+LM2+X1 +X1_int+cluster(id),data=LMdata2,method="breslow")
bh_supercox3_int<-basehaz(LMsupercox3_int,centered=FALSE)

#LM4Int: LM4 + interactions
LMsupercox4_int<-coxph(Surv(Tstart,wsurvtime,survstatus)~illness+illnessX+illnessX2+illness_tau+illness_tau2+LM1+LM2+X1+X1_int+cluster(id),data=LMdata2,method="breslow")
bh_supercox4_int<-basehaz(LMsupercox4_int,centered=FALSE)

#######################################################
#Joint Models 
#######################################################
#MM: Markov model 
fit_MM <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1,
              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1,
              formula12=Hist(time=observed.lifetime,event=seen.exit)~X1,data=train_data.id,
              conf.int=FALSE,method="Weib") 

#MMCox: semi-parametric Markov model (Cox proportional hazards) 
#Create dataset
data_fit<-data.frame(matrix(nrow=nrow(train_data.id),ncol=6))
names(data_fit)<-c("times1","delta","times2","time","status","X1")
data_fit$times1<-ifelse(train_data.id$seen.ill==1,train_data.id$observed.illtime,train_data.id$observed.lifetime)
data_fit$delta<-train_data.id$seen.ill
data_fit$time<-train_data.id$observed.lifetime
data_fit$times2<-ifelse(train_data.id$seen.ill==1,train_data.id$observed.lifetime-train_data.id$observed.illtime,0)
data_fit$status<-train_data.id$seen.exit
data_fit$X1<-train_data.id$X1

tmat<-trans.illdeath() 
data_fit_long<-msprep(data=data_fit,
                      trans=tmat,
                      time=c(NA,"times1","time"),
                      status=c(NA,"delta","status"),
                      keep=c("X1","times2"))
data_fit_long$X1<-as.factor(data_fit_long$X1)
data_fit_long<-expand.covs(data_fit_long,covs=c("X1"))

fit_MMCox<-coxph(Surv(Tstart,Tstop,status)~X11.1+X11.2+X11.3+strata(trans),data=data_fit_long,method="breslow")

##########################################################################
#Copula functions
##########################################################################
##################################
#Marginal model for T: h(t)=h0(t)exp(X)
#"survtime": minimum of censtime and lifetime 
#"survstatus": indicator of death
#"X1": covariate
##################################
#Cox PH for T
cox_T<-coxph(Surv(survtime,seen.exit)~X1,data=train_data.id)

weib_T<-survreg(Surv(survtime,seen.exit)~X1,data=train_data.id,dist="weibull")

shape.weib<-1/weib_T$scale
scale.weib<-exp(weib_T$coef[1])
coef.weib<- -weib_T$coef[-1]/weib_T$scale

F_T_func_cond<-function(dat)
{
  xdata1<-as.numeric(dat[1])
  tt<-as.numeric(dat[2])
  tt_LM<-as.numeric(dat[3])
  
  bh<-list(basehaz(cox_T,FALSE))
 
  bet<-list(cox_T$coef)
 
  names(bet)<-c("cox_T")

  sfi<-bh
  sfi<-lapply(sfi,function(x) data.frame(x,haz0=diff(c(0,x$hazard))))
  
  sfi[[1]]["haz"]<-sfi[[1]]$haz0*as.numeric(exp(bet[["cox_T"]][1]*xdata1))
  
  sfi<-lapply(sfi,function(x) data.frame(x,Haz=cumsum(x$haz)))
  tmp<-lapply(sfi,function(x) evalstep(x$time,x$Haz,c(tt,tt_LM),subst=0))
  Fw<-lapply(tmp,function(x) ((1-exp(-x[1]))-(1-exp(-x[2])))/exp(-x[2]))  
  
  Fw_weib_tt<-exp(-(tt/scale.weib)^shape.weib*exp(coef.weib*xdata1))
  Fw_weib_tt_LM<-exp(-(tt_LM/scale.weib)^shape.weib*exp(coef.weib*xdata1))
  Fw_weib<-((1-Fw_weib_tt)-(1-Fw_weib_tt_LM))/Fw_weib_tt_LM
  
  ret<-c(unlist(Fw),Fw_weib)
  
  return(ret)
}

FT_cond_vec<-apply(train_data[c("X1","survtime","inspection.time")],1,F_T_func_cond)
FT_cond_vec<-t(FT_cond_vec)
FT_names<-c("cox","weib")
colnames(FT_cond_vec)<-FT_names

##################################
#Model for Z*: N(mu=alpha0+alpha1*tau+alpha2*tau^2+alpha3*X,1)
#Model Z* using a probit 
#Z=I(Z*<0)
##################################
F_Zstar<-function(zstar,mean_Zstar,log.prob=FALSE)
{
  pnorm(zstar,mean=mean_Zstar,sd=1,log.p=log.prob)
}

#Models for Z 
const<-illness~1
noLM<-illness~X1
nocovs<-illness~LM
simp<-illness~LM+X1
int<-illness~LM*X1

mods<-list(const,noLM,nocovs,simp,int)
mods_names<-c("const","noLM","nocovs","simp","int")
names(mods)<-mods_names

for(mean_string in mods_names)
{
  mod_mean_Zstar<-mods[[mean_string]]
  mod_Zstar<-glm(mod_mean_Zstar,data=train_data,family=binomial(link="probit"))
  assign(paste0("mod_Zstar_",mean_string),mod_Zstar)
  assign(paste0("mod_Zstar_mean",mean_string),mod_Zstar$coef)
}

##################################
#Joint distribution
##################################
FZstarT<-function(FT_cond,zstar,mean_Zstar,rho)
{
  a<-qnorm(FT_cond,mean=0,sd=1)
  b<-qnorm(F_Zstar(zstar,mean_Zstar),mean=0,sd=1)
  ret<-ifelse(a==-Inf,0,
              ifelse(a==Inf,F_Zstar(zstar,mean_Zstar),
                     ifelse(b==Inf,FT_cond,pbivnorm(x=a,y=b,rho=rho))))
  return(ret)
}

##################################
#Maximize joint likelihood
##################################
#g: parameter estimates
#coef_dat: model matrix for rho parameter 
#FTcond: P(T<t|T>tau) 
loglik_ZT<-function(g,coef_dat,FTcond,meanZ,dat)
{
  #Reparameterize rho using Fisher's z-transformation since rho (-1,1)
  eta<-coef_dat%*%g
  rho<-1-2/(1+exp(2*eta))
  
  q1<-qnorm(FTcond,mean=0,sd=1)
  q2<-qnorm(F_Zstar(0,meanZ,log.prob=TRUE),mean=0,sd=1,log.p=TRUE)
  
  #L1: P(T=t,Z=0) 
  L1<-pnorm((q2-rho*q1)/sqrt(1-rho^2))
  
  #L2: P(T>t,Z=0)
  L2<-F_Zstar(0,meanZ)-FZstarT(FTcond,0,meanZ,rho)
  L2<-ifelse(L2<0,0,L2)
  
  #L3: P(T=t,Z=1)
  L3<-pnorm(-(q2-rho*q1)/sqrt(1-rho^2))
  
  #L4: P(T>t,Z=1)
  L4<-(1-F_Zstar(0,meanZ))-FTcond+FZstarT(FTcond,0,meanZ,rho)
  L4<-ifelse(L4<0,0,L4)
  
  L1_cont<-log(L1[dat$illness==0&dat$seen.exit==1])
  L2_cont<-log(L2[dat$illness==0&dat$seen.exit==0])
  L3_cont<-log(L3[dat$illness==1&dat$seen.exit==1])
  L4_cont<-log(L4[dat$illness==1&dat$seen.exit==0])
  
  if(any(is.na(-(sum(L1_cont)+sum(L2_cont)+sum(L3_cont)+sum(L4_cont)))))
  {
    print(g)
    stop()
  }
  -(sum(L1_cont,na.rm=TRUE)+sum(L2_cont,na.rm=TRUE)+sum(L3_cont,na.rm=TRUE)+sum(L4_cont,na.rm=TRUE))
  ret<--(sum(L1_cont)+sum(L2_cont)+sum(L3_cont)+sum(L4_cont))
  return(ret)
}

for(rho_string in mods_names)
{
  for(mean_string in mods_names)
  {
    for(FT_string in FT_names)
    {
      print(c(rho_string,mean_string,FT_string))
      modZ<-get(paste0("mod_Zstar_",mean_string))
      meanZ<-predict(modZ,newdata=train_data)
      
      coef_dat<-model.matrix(mods[[rho_string]],data=train_data)
      
      start_val<-rep(0,ncol(coef_dat))
      temp<-optim(start_val,loglik_ZT,
                  coef_dat=coef_dat,
                  FTcond=FT_cond_vec[,FT_string],
                  meanZ=meanZ,
                  dat=train_data,
                  method="Nelder-Mead") #,control=list(maxit=5000,trace=FALSE))
      assign(paste0("MLE_rho",rho_string,"_mean",mean_string,"_",FT_string),temp$par)
    }
  }
}


########################################################
#Prediction function 
########################################################
models<-c(models,cop_labs)
nmodels<-length(models)

BS_full<-NULL 

for(t0 in LMx)
{
  # t0<-1
  print(t0)
  
  pred1<-summary(survfit(Surv(survtime,survstatus)~1,data=train_data.id),time=t0)$surv
  pred2<-summary(survfit(Surv(survtime,survstatus)~1,data=train_data.id),time=t0+w_predict)$surv
  pred.full<-(pred1-pred2)/pred1
  
  #Subset of data with those who are still alive at time t0
  sub_dat<-subset(test_data,test_data$survtime>t0)
  sub_dat.id<-subset(test_data.id,test_data.id$survtime>t0)
  
  #select those with event time>LM and their data up to LM 
  sub<-subset(sub_dat,sub_dat$inspection.time<=t0) 
  sub<-sub[order(sub$id,sub$inspection.time),]
  #select last entry for LM models 
  sub_LM<-do.call("rbind", as.list(by(sub,sub$id,tail,n=1)))
  sub_LM$LM<-t0
  
  data.temp<-rep(pred.full,nrow(sub_dat.id))
  
  zdata=sub_LM$illness
  xdata=sub_LM$X1 #as.numeric(levels(sub_LM$group))[sub_LM$group]
  data.temp<-cbind(data.temp,BSpredict_Truth(zdata=zdata,t0=t0,xdata1=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LM2(bh=bh_supercox2,bet=LMsupercox2$coef,zdata=zdata,t0=t0,xdata1=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LM3(bh=bh_supercox3,bet=LMsupercox3$coef,zdata=zdata,t0=t0,xdata1=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LM4(bh=bh_supercox4,bet=LMsupercox4$coef,zdata=zdata,t0=t0,xdata1=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LM2Int(bh=bh_supercox2_int,bet=LMsupercox2_int$coef,zdata=zdata,t0=t0,xdata1=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LM3Int(bh=bh_supercox3_int,bet=LMsupercox3_int$coef,zdata=zdata,t0=t0,xdata1=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LM4Int(bh=bh_supercox4_int,bet=LMsupercox4_int$coef,zdata=zdata,t0=t0,xdata1=xdata,w_predict))
  
  data.temp<-cbind(data.temp,BSpredict_MM(fit=fit_MM,zdata=zdata,t0=t0,xdata1=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_MMCox(fitCox=fit_MMCox,zdata=zdata,t0=t0,xdata1=xdata,w_predict))

    for(mean_string in mods_names) 
  {
    for(FT_string in c("cox","weib"))
      {
        for(rho_string in mods_names) 
        {
  print(c(mean_string,FT_string,rho_string))
  data.temp<-cbind(data.temp,BSpredict_Copula(zdata=zdata,t0=t0,
                                              xdata=xdata,
                                              w_predict=w_predict,mean_string=mean_string,
                                              FT_string=FT_string,rho_string=rho_string))
  
        }
      }
  }

  df.temp<-data.frame(sub_LM,data.temp)
  
  write.table(df.temp,paste0(dir0,"/Results_mid/df",tail_str),row.names=FALSE,append=TRUE,col.names=FALSE)
  
  BS_full<-rbind(BS_full,df.temp)
}


names(BS_full)<-c(names(sub_LM),models)

#########################
#RMSE
#########################
RMSE_all<-sqrt(colMeans((BS_full-BS_full$Truth)^2))[models]

RMSE_X<-NULL
RMSE_LM<-NULL
RMSE_LM_X<-NULL
for(i in LMx)
{
  RMSE_sub<-subset(BS_full,BS_full$LM==i)
  RMSE_LM<-rbind(RMSE_LM,c(LM=i,sqrt(colMeans((RMSE_sub-RMSE_sub$Truth)^2))[models]))
  for(j in c(0,1))
  {
    RMSE_sub<-subset(BS_full,BS_full$X1==j&BS_full$LM==i)
    RMSE_LM_X<-rbind(RMSE_LM_X,c(LM=i,X1=j,sqrt(colMeans((RMSE_sub-RMSE_sub$Truth)^2))[models]))
    
    RMSE_sub<-subset(BS_full,BS_full$X1==j)
    RMSE_X<-rbind(RMSE_X,c(X1=j,sqrt(colMeans((RMSE_sub-RMSE_sub$Truth)^2))[models]))
  }
}

#########################
#AUC and BS 
#########################
df_BS<-data.frame(matrix(nrow=length(LMx),ncol=nmodels))
df_AUC<-data.frame(matrix(nrow=length(LMx),ncol=nmodels))
for(j in 1:nmodels)
{
  print(j)
  if(all(is.na(BS_full[models[j]]))==TRUE)
  {
    df_BS[,j]<-NA
    df_AUC[,j]<-NA
  } else {
    pred.error<-PE(models[j],LMx,thor,BS_full)
    df_BS[,j]<-pred.error[[1]]
    df_AUC[,j]<-pred.error[[2]]
  }
}

names(df_BS)<-models
names(df_AUC)<-models


