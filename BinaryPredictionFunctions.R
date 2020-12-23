#Prediction functions 

#######################################################
#Copula Model 
#######################################################
BSpredict_Copula<-function(zdata,t0,xdata1,w_predict,mean_string,FT_string,rho_string)
{
  Fw<-NULL
  tti<-t0
  
  for(i in 1:length(zdata))
  {
    covdat<-data.frame(LM=tti,X1=xdata1[i],illness=zdata[i])
    
    #modZ
    modZ<-get(paste0("mod_Zstar_",mean_string))
    meanZ<-predict(modZ,newdata=covdat)
    
    rho_param<-get(paste0("MLE_rho",rho_string,"_mean",mean_string,"_",FT_string))
    
    coef_dat<-model.matrix(mods[[rho_string]],data=covdat)
    
    eta_Ystar_LM<-coef_dat%*%rho_param
    rho_LM<-1-2/(1+exp(2*eta_Ystar_LM))
    
    F_T_temp<-F_T_func_cond(c(xdata1[i],tti+w_predict,tti))
    names(F_T_temp)<-FT_names
    FT_tau_w_cond<-F_T_temp[FT_string]
    
    num_0<-FZstarT(zstar=0,FT_cond=FT_tau_w_cond,mean_Zstar=meanZ,rho=rho_LM)
    denom_0<-F_Zstar(zstar=0,meanZ)
    
    num_1<-FT_tau_w_cond-num_0
    denom_1<-1-denom_0 #(1-FT_tau)-d1_0+n2_0
    
    if(zdata[i]==0)
    {
      ret<-num_0/denom_0
      ret<-ifelse(ret<0,0,ret)
      Fw<-c(Fw,ret)
    }
    else
    {
      ret<-num_1/denom_1
      ret<-ifelse(ret<0,0,ret)
      Fw<-c(Fw,ret)
    }
  }
  return(Fw)
}

######################################################
#Compute AUC and BS 
######################################################
models<-c("Null","Truth","LM2","LM3","LM4","LM2Int","LM3Int","LM4Int","MM","MMCox")

cop_labs<-NULL
for(mean_string in mods_names) 
{
  for(FT_string in c("cox","weib"))
  {
    for(rho_string in mods_names) 
    {
      cop_labs<-c(cop_labs,paste0(mean_string,"_",FT_string,"_",rho_string))
    }
  }
}
#######################################################
#Landmark Models 
#######################################################
#LM1 
BSpredict_LM1<-function(bh,bet,zdata,t0,xdata1,w_predict)
{
  tti <- t0
  sfi<-subset(bh,bh$strata==paste0("LM=",tti))
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw <- NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(tti) + bet[3]*zdata[i]*(tti)^2+ 
                                bet[4]*xdata1[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LM1Int 
BSpredict_LM1Int<-function(bh,bet,zdata,t0,xdata1,w_predict)
{
  tti<-t0
  sfi<-subset(bh,bh$strata==paste0("LM=",tti))
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw <- NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(tti) + bet[3]*zdata[i]*(tti)^2+ 
                                bet[4]*xdata1[i]+
                                bet[5]*xdata1[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}


#LSM1 
BSpredict_LSM1<-function(bh,bet,zdata,t0,xdata1,w_predict,v_predict)
{
  tti <- t0
  sfi<-subset(bh,bh$strata==paste0("LM=",tti))
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw <- NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(tti) + bet[3]*zdata[i]*(tti)^2+ 
                                bet[4]*xdata1[i]+bet[5]*v_predict[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LM2 
BSpredict_LM2<-function(bh,bet,zdata,t0,xdata1,w_predict)
{
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw<-NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*g1(tti) + bet[3]*zdata[i]*g2(tti) +
                                bet[4]*g1(tti)+bet[5]*g2(tti)+bet[6]*xdata1[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LM2Int
BSpredict_LM2Int<-function(bh,bet,zdata,t0,xdata1,w_predict)
{
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw<-NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*g1(tti) + bet[3]*zdata[i]*g2(tti) +
                                bet[4]*g1(tti)+bet[5]*g2(tti)+bet[6]*xdata1[i]+
                                bet[7]*xdata1[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}


#LSM2
BSpredict_LSM2<-function(bh,bet,zdata,t0,xdata1,w_predict,v_predict)
{
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw<-NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*g1(tti) + bet[3]*zdata[i]*g2(tti) +
                                bet[4]*g1(tti)+bet[5]*g2(tti)+bet[6]*xdata1[i]+bet[7]*v_predict[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LM3
BSpredict_LM3<-function(bh,bet,zdata,t0,xdata1,w_predict)
{
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw<-NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(sfi$time-tti) + bet[3]*zdata[i]*(sfi$time-tti)^2
                              + bet[4]*g1(tti) + bet[5]*g2(tti)+bet[6]*xdata1[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LM3Int 
BSpredict_LM3Int<-function(bh,bet,zdata,t0,xdata1,w_predict)
{
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw<-NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(sfi$time-tti) + bet[3]*zdata[i]*(sfi$time-tti)^2
                              + bet[4]*g1(tti) + bet[5]*g2(tti)+bet[6]*xdata1[i]
                              + bet[7]*xdata1[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}


#LSM3 
BSpredict_LSM3<-function(bh,bet,zdata,t0,xdata1,w_predict,v_predict)
{
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw<-NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(sfi$time-tti) + bet[3]*zdata[i]*(sfi$time-tti)^2
                              + bet[4]*g1(tti) + bet[5]*g2(tti)+bet[6]*xdata1[i]+bet[7]*v_predict[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}


#LM4
BSpredict_LM4<-function(bh,bet,zdata,t0,xdata1,w_predict)
{
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw<-NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(sfi$time-tti) + bet[3]*zdata[i]*(sfi$time-tti)^2
                              + bet[4]*zdata[i]*tti + bet[5]*zdata[i]*tti^2
                              + bet[6]*g1(tti) + bet[7]*g2(tti)+bet[8]*xdata1[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LMInt4
BSpredict_LM4Int<-function(bh,bet,zdata,t0,xdata1,w_predict)
{
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw<-NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(sfi$time-tti) + bet[3]*zdata[i]*(sfi$time-tti)^2
                              + bet[4]*zdata[i]*tti + bet[5]*zdata[i]*tti^2
                              + bet[6]*g1(tti) + bet[7]*g2(tti)+bet[8]*xdata1[i]
                              + bet[9]*xdata1[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LSM4
BSpredict_LSM4<-function(bh,bet,zdata,t0,xdata1,w_predict,v_predict)
{
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  Fw<-NULL
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(sfi$time-tti) + bet[3]*zdata[i]*(sfi$time-tti)^2
                              + bet[4]*zdata[i]*tti + bet[5]*zdata[i]*tti^2
                              + bet[6]*g1(tti) + bet[7]*g2(tti)+bet[8]*xdata1[i]+bet[9]*v_predict[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#######################################################
#Joint Models 
#######################################################
#MM
BSpredict_MM<-function(fit,zdata,t0,xdata1,w_predict)
{
  params<-fit$modelPar
  betas<-fit$coef
  
  scale01<-1/params[2] #Weibull scale for latent illness time
  shape01<-params[1]    #Weibull shape for latent illness time
  scale02<-1/params[4]  #Weibull scale for latent life time
  shape02<-params[3]   #Weibull shape for latent life time
  scale12<-1/params[6]  #Weibull scale for latent wait time
  shape12<-params[5]   #Weibull shape for latent wait time
  beta_X1_01<-betas[1]
  beta_X1_02<-betas[2]
  beta_X1_12<-betas[3]
  
  lambda12<-function(x,xdata1) {shape12*(1/scale12)^(shape12)*(x)^(shape12-1)*exp(beta_X1_12*xdata1)}
  lambda01<-function(x,xdata1) {shape01*(1/scale01)^(shape01)*(x)^(shape01-1)*exp(beta_X1_01*xdata1)}
  lambda02<-function(x,xdata1) {shape02*(1/scale02)^(shape02)*(x)^(shape02-1)*exp(beta_X1_02*xdata1)}
  
  H12<-function(x,xdata1) {-pweibull(x, scale=scale12, shape=shape12, lower = FALSE, log = TRUE)*exp(beta_X1_12*xdata1)}
  H01<-function(x,xdata1) {-pweibull(x, scale=scale01, shape=shape01, lower = FALSE, log = TRUE)*exp(beta_X1_01*xdata1)}
  H02<-function(x,xdata1) {-pweibull(x, scale=scale02, shape=shape02, lower = FALSE, log = TRUE)*exp(beta_X1_02*xdata1)}
  
  tti<-t0
  Fw<-NULL
  
  resdist0<-function(v,xdata1) {exp(-((H02(v,xdata1)-H02(tti,xdata1))+(H01(v,xdata1)-H01(tti,xdata1))))*lambda01(v,xdata1)*
      exp(-(H12(tti+w_predict,xdata1)-H12(v,xdata1)))}
  
  Fw_return<-function(xdata1,zdata)
  {
    Fw1<-1-exp(-(H12(tti+w_predict,xdata1)-H12(tti,xdata1)))
    
    temp<-seq(tti,tti+w_predict,by=0.001)
    temp2<-sum(resdist0(temp,xdata1))*0.001
    
    Fw0<-1-(exp(-((H02(tti+w_predict,xdata1)-H02(tti,xdata1))+(H01(tti+w_predict,xdata1)-H01(tti,xdata1))))+temp2)
    
    return(zdata*Fw1+(1-zdata)*Fw0)
  }
  
  for(i in 1:length(zdata))
  {
    Fw<-c(Fw,Fw_return(xdata1[i],zdata[i]))
  }
  return(Fw)
}

#MSM
BSpredict_MSM<-function(betas,zdata,t0,xdata1,w_predict,v_predict)
{
  scale01<-betas[1]   #Weibull scale for latent illness time
  shape01<-betas[2]    #Weibull shape for latent illness time
  scale02<-betas[3]  #Weibull scale for latent life time
  shape02<-betas[4]   #Weibull shape for latent life time
  scale12<-betas[5]  #Weibull scale for latent wait time
  shape12<-betas[6]   #Weibull shape for latent wait time
  beta_X1_01<-betas[7]
  beta_X1_02<-betas[8]
  beta_X1_12<-betas[9]
  beta_v<-betas[10]
  
  lambda12<-function(x,xdata1,v_predict) {shape12*(1/scale12)^(shape12)*(x)^(shape12-1)*exp(beta_X1_12*xdata1+beta_v*v_predict)}
  lambda01<-function(x,xdata1) {shape01*(1/scale01)^(shape01)*(x)^(shape01-1)*exp(beta_X1_01*xdata1)}
  lambda02<-function(x,xdata1) {shape02*(1/scale02)^(shape02)*(x)^(shape02-1)*exp(beta_X1_02*xdata1)}
  
  H12<-function(x,xdata1,v_predict) {-pweibull(x, scale=scale12, shape=shape12, lower = FALSE, log = TRUE)*exp(beta_X1_12*xdata1+beta_v*v_predict)}
  H01<-function(x,xdata1) {-pweibull(x, scale=scale01, shape=shape01, lower = FALSE, log = TRUE)*exp(beta_X1_01*xdata1)}
  H02<-function(x,xdata1) {-pweibull(x, scale=scale02, shape=shape02, lower = FALSE, log = TRUE)*exp(beta_X1_02*xdata1)}
  
  tti<-t0
  Fw<-NULL
  
  resdist0<-function(v,xdata1) {exp(-((H02(v,xdata1)-H02(tti,xdata1))+(H01(v,xdata1)-H01(tti,xdata1))))*lambda01(v,xdata1)*
      exp(-(H12(tti+w_predict,xdata1,v)-H12(v,xdata1,v)))}
  
  Fw_return<-function(xdata1,zdata,v_predict)
  {
    Fw1<-1-exp(-(H12(tti+w_predict,xdata1,v_predict)-H12(tti,xdata1,v_predict)))
    
    temp<-seq(tti,tti+w_predict,by=0.001)
    temp2<-sum(resdist0(temp,xdata1))*0.001
    
    Fw0<-1-(exp(-((H02(tti+w_predict,xdata1)-H02(tti,xdata1))+(H01(tti+w_predict,xdata1)-H01(tti,xdata1))))+temp2)
    
    return(zdata*Fw1+(1-zdata)*Fw0)
  }
  
  for(i in 1:length(zdata))
  {
    Fw<-c(Fw,Fw_return(xdata1[i],zdata[i],v_predict[i]))
  }
  return(Fw)
}

#MMCox
BSpredict_MMCox<-function(fitCox,zdata,t0,xdata1,w_predict)
{
  
  predict_df_MMCox<-data.frame(X11.1=rep(0,3),X11.2=rep(0,3),X11.3=rep(0,3),
                               trans=c(1:3),strata=c(1:3))
  
  fitCox_sum<-msfit(object=fitCox,predict_df_MMCox,trans=tmat)
  
  beta_Cox<-fitCox$coefficients
  
  sfi_01<-subset(fitCox_sum$Haz,trans==1)
  sfi_01$haz<-diff(c(0,sfi_01$Haz))
  sfi_02<-subset(fitCox_sum$Haz,trans==2)
  sfi_12<-subset(fitCox_sum$Haz,trans==3)
  Fw<-NULL
  tti<-t0
  
  Fw_MMCox<-function(xdata1,zdata)
  {
    sfi_01$Haz_cov<-sfi_01$Haz*exp(beta_Cox[1]*xdata1)
    sfi_01$haz_cov<-sfi_01$haz*exp(beta_Cox[1]*xdata1)
    
    sfi_02$Haz_cov<-sfi_02$Haz*exp(beta_Cox[2]*xdata1)
    
    sfi_12$Haz_cov<-sfi_12$Haz*exp(beta_Cox[3]*xdata1)
    
    tmp1<-evalstep(sfi_01$time,sfi_01$Haz_cov,c(tti,tti+w_predict),subst=0)
    tmp2<-evalstep(sfi_02$time,sfi_02$Haz_cov,c(tti,tti+w_predict),subst=0)
    Fw0_1<-exp(-((tmp1[2]-tmp1[1])+(tmp2[2]-tmp2[1])))
    
    tmp1<-evalstep(sfi_12$time,sfi_12$Haz_cov,c(tti,tti+w_predict),subst=0)
    Fw1<-1-exp(-(tmp1[2]-tmp1[1]))
    
    times<-subset(sfi_01$time,sfi_01$time>=tti&sfi_01$time<=tti+w_predict)
    temp_k<-NULL
    for(k in times)
    {
      tmp1<-evalstep(sfi_01$time,sfi_01$Haz_cov,c(tti,k),subst=0)
      tmp2<-evalstep(sfi_02$time,sfi_02$Haz_cov,c(tti,k),subst=0)
      int1<-exp(-((tmp1[2]-tmp1[1])+(tmp2[2]-tmp2[1])))
      int2<-sfi_01[sfi_01$time==k,]$haz_cov
      tmp3<-evalstep(sfi_12$time,sfi_12$Haz_cov,c(k,tti+w_predict),subst=0)
      int3<-exp(-(tmp3[2]-tmp3[1]))
      temp_k<-c(temp_k,int1*int2*int3)
    }
    Fw0_2<-sum(temp_k)
    Fw0<-1-(Fw0_1+Fw0_2)
    return(zdata*Fw1+(1-zdata)*Fw0)
  }
  
  FZ1_X1<-Fw_MMCox(xdata1=1,zdata=1)
  FZ0_X1<-Fw_MMCox(xdata1=1,zdata=0)
  FZ1_X0<-Fw_MMCox(xdata1=0,zdata=1)
  FZ0_X0<-Fw_MMCox(xdata1=0,zdata=0)
  
  for(i in 1:length(zdata))
  {
    Fw<-c(Fw,zdata[i]*xdata1[i]*FZ1_X1+
            (1-zdata[i])*xdata1[i]*FZ0_X1+
            zdata[i]*(1-xdata1[i])*FZ1_X0+
            (1-zdata[i])*(1-xdata1[i])*FZ0_X0)
  }
  return(Fw)
}

#MSMCox
BSpredict_MSMCox<-function(fitCox,zdata,t0,xdata1,w_predict,v_predict)
{
  predict_df_MSMCox<-data.frame(X11.1=rep(0,3),X11.2=rep(0,3),X11.3=rep(0,3),
                                trans=c(1:3),V=rep(0,3),strata=c(1:3))
  
  fitCox_sum<-msfit(object=fitCox,predict_df_MSMCox,trans=tmat)
  
  beta_Cox<-fitCox$coefficients
  
  sfi_01<-subset(fitCox_sum$Haz,trans==1)
  sfi_01$haz<-diff(c(0,sfi_01$Haz))
  sfi_02<-subset(fitCox_sum$Haz,trans==2)
  sfi_12<-subset(fitCox_sum$Haz,trans==3)
  
  Fw<-NULL
  tti<-t0
  for(i in 1:length(zdata))
  {
    Fw0_1<-NULL
    Fw0_2<-NULL
    Fw1<-NULL
    
    sfi_01$Haz_cov<-sfi_01$Haz*exp(beta_Cox[1]*xdata1[i])
    sfi_02$Haz_cov<-sfi_02$Haz*exp(beta_Cox[2]*xdata1[i])
    sfi_12$Haz_cov<-sfi_12$Haz*exp(beta_Cox[3]*xdata1[i]+beta_Cox[4]*v_predict[i])
    sfi_01$haz_cov<-sfi_01$haz*exp(beta_Cox[1]*xdata1[i])
    
    if(zdata[i]==1){
      tmp1<-evalstep(sfi_12$time,sfi_12$Haz_cov,c(tti,tti+w_predict),subst=0)
      Fw1<-1-exp(-(tmp1[2]-tmp1[1]))
      Fw0<-0
    } else {
      tmp1<-evalstep(sfi_01$time,sfi_01$Haz_cov,c(tti,tti+w_predict),subst=0)
      tmp2<-evalstep(sfi_02$time,sfi_02$Haz_cov,c(tti,tti+w_predict),subst=0)
      Fw0_1<-exp(-((tmp1[2]-tmp1[1])+(tmp2[2]-tmp2[1])))
      
      times<-subset(sfi_01$time,sfi_01$time>=tti&sfi_01$time<=tti+w_predict)
      temp_k<-NULL
      for(k in times)
      {
        tmp1<-evalstep(sfi_01$time,sfi_01$Haz_cov,c(tti,k),subst=0)
        tmp2<-evalstep(sfi_02$time,sfi_02$Haz_cov,c(tti,k),subst=0)
        int1<-exp(-((tmp1[2]-tmp1[1])+(tmp2[2]-tmp2[1])))
        int2<-sfi_01[sfi_01$time==k,]$haz_cov
        sfi_12$Haz_cov<-sfi_12$Haz*exp(beta_Cox[3]*xdata1[i]+beta_Cox[4]*k)
        tmp3<-evalstep(sfi_12$time,sfi_12$Haz_cov,c(k,tti+w_predict),subst=0)
        int3<-exp(-(tmp3[2]-tmp3[1]))
        temp_k<-c(temp_k,int1*int2*int3)
      }
      Fw0_2<-sum(temp_k)
      Fw0<-1-(Fw0_1+Fw0_2)
      Fw1<-0
    }
    Fw<-c(Fw,zdata[i]*Fw1+(1-zdata[i])*Fw0)
  }
  return(Fw)
}

#SMM
BSpredict_SMM<-function(betas,zdata,t0,xdata1,w_predict,v_predict)
{
  scale01<-betas[1]   #Weibull scale for latent illness time
  shape01<-betas[2]    #Weibull shape for latent illness time
  scale02<-betas[3]  #Weibull scale for latent life time
  shape02<-betas[4]   #Weibull shape for latent life time
  scale12<-betas[5]  #Weibull scale for latent wait time
  shape12<-betas[6]   #Weibull shape for latent wait time
  beta_X1_01<-betas[7]
  beta_X1_02<-betas[8]
  beta_X1_12<-betas[9]
  
  lambda12<-function(x,xdata1,v) {shape12*(1/scale12)^(shape12)*(x-v)^(shape12-1)*exp(beta_X1_12*xdata1)}
  lambda01<-function(x,xdata1) {shape01*(1/scale01)^(shape01)*(x)^(shape01-1)*exp(beta_X1_01*xdata1)}
  lambda02<-function(x,xdata1) {shape02*(1/scale02)^(shape02)*(x)^(shape02-1)*exp(beta_X1_02*xdata1)}
  
  H12<-function(x,xdata1,v) {-pweibull(x-v, scale=scale12, shape=shape12, lower = FALSE, log = TRUE)*exp(beta_X1_12*xdata1)}
  H01<-function(x,xdata1) {-pweibull(x, scale=scale01, shape=shape01, lower = FALSE, log = TRUE)*exp(beta_X1_01*xdata1)}
  H02<-function(x,xdata1) {-pweibull(x, scale=scale02, shape=shape02, lower = FALSE, log = TRUE)*exp(beta_X1_02*xdata1)}
  
  Fw<-NULL
  tti<-t0
  
  resdist0<-function(v,xdata1) {exp(-((H02(v,xdata1)-H02(tti,xdata1))+(H01(v,xdata1)-H01(tti,xdata1))))*lambda01(v,xdata1)*
      exp(-(H12(tti+w_predict,xdata1,v)-H12(v,xdata1,v)))}
  
  Fw_return<-function(xdata1,zdata,v_predict)
  {
    Fw1<-1-exp(-(H12(tti+w_predict,xdata1,v_predict)-H12(tti,xdata1,v_predict)))
    
    temp<-seq(tti,tti+w_predict,by=0.001)
    temp2<-sum(resdist0(temp,xdata1))*0.001
    
    Fw0<-1-(exp(-((H02(tti+w_predict,xdata1)-H02(tti,xdata1))+(H01(tti+w_predict,xdata1)-H01(tti,xdata1))))+temp2)
    
    return(zdata*Fw1+(1-zdata)*Fw0)
  }
  
  for(i in 1:length(zdata))
  {
    Fw<-c(Fw,Fw_return(xdata1[i],zdata[i],v_predict[i]))
  }
  return(Fw)
}

#Truth
BSpredict_Truth<-function(zdata,t0,xdata1,w_predict)
{
  scale01<-scale.illtime #Weibull scale for latent illness time
  shape01<-shape.illtime    #Weibull shape for latent illness time
  scale02<-scale.lifetime  #Weibull scale for latent life time
  shape02<-shape.lifetime   #Weibull shape for latent life time
  scale12<-scale.waittime  #Weibull scale for latent wait time
  shape12<-shape.waittime   #Weibull shape for latent wait time
  
  tti<-t0
  Fw<-NULL
  lambda12<-function(x,xdata1) {shape12*(1/scale12)^(shape12)*(x)^(shape12-1)*exp(beta_X1_12*xdata1)}
  lambda01<-function(x,xdata1) {shape01*(1/scale01)^(shape01)*(x)^(shape01-1)*exp(beta_X1_01*xdata1)}
  lambda02<-function(x,xdata1) {shape02*(1/scale02)^(shape02)*(x)^(shape02-1)*exp(beta_X1_02*xdata1)}
  
  
  # H12<-function(x) {(x/scale12)^shape12}
  H12<-function(x,xdata1) {-pweibull(x, scale=scale12, shape=shape12, lower = FALSE, log = TRUE)*exp(beta_X1_12*xdata1)}
  H01<-function(x,xdata1) {-pweibull(x, scale=scale01, shape=shape01, lower = FALSE, log = TRUE)*exp(beta_X1_01*xdata1)}
  H02<-function(x,xdata1) {-pweibull(x, scale=scale02, shape=shape02, lower = FALSE, log = TRUE)*exp(beta_X1_02*xdata1)}
  
  # NO clock reset
  resdist0<-function(v,xdata1) {exp(-((H02(v,xdata1)-H02(tti,xdata1))+(H01(v,xdata1)-H01(tti,xdata1))))*lambda01(v,xdata1)*
      exp(-(H12(tti+w_predict,xdata1)-H12(v,xdata1)))}
  
  Fw_return<-function(xdata1,zdata)
  {
    Fw1<-1-exp(-(H12(tti+w_predict,xdata1)-H12(tti,xdata1)))
    
    temp<-seq(tti,tti+w_predict,by=0.001)
    temp2<-sum(resdist0(temp,xdata1))*0.001
    
    Fw0<-1-(exp(-((H02(tti+w_predict,xdata1)-H02(tti,xdata1))+(H01(tti+w_predict,xdata1)-H01(tti,xdata1))))+temp2)
    
    return(zdata*Fw1+(1-zdata)*Fw0)
  }
  
  for(i in 1:length(zdata))
  {
    Fw<-c(Fw,Fw_return(xdata1[i],zdata[i]))
  }
  return(Fw)
}


#Truth (SemiMarkov)
BSpredict_Truth_SM<-function(zdata,t0,xdata1,w_predict,v_predict)
{
  scale01_MSM_SM<-scale.illtime #Weibull scale for latent illness time
  shape01_MSM_SM<-shape.illtime    #Weibull shape for latent illness time
  scale02_MSM_SM<-scale.lifetime  #Weibull scale for latent life time
  shape02_MSM_SM<-shape.lifetime   #Weibull shape for latent life time
  scale12_MSM_SM<-scale.waittime  #Weibull scale for latent wait time
  shape12_MSM_SM<-shape.waittime   #Weibull shape for latent wait time
  
  tti<-t0
  Fw<-NULL
  lambda12_MSM_SM<-function(x,xdata1,v) {shape12_MSM_SM*(1/scale12_MSM_SM)^(shape12_MSM_SM)*(x-v)^(shape12_MSM_SM-1)*exp(beta_X1_12*xdata1)}
  lambda01_MSM_SM<-function(x,xdata1) {shape01_MSM_SM*(1/scale01_MSM_SM)^(shape01_MSM_SM)*(x)^(shape01_MSM_SM-1)*exp(beta_X1_01*xdata1)}
  lambda02_MSM_SM<-function(x,xdata1) {shape02_MSM_SM*(1/scale02_MSM_SM)^(shape02_MSM_SM)*(x)^(shape02_MSM_SM-1)*exp(beta_X1_02*xdata1)}
  
  H12_MSM_SM<-function(x,xdata1,v) {-pweibull(x-v, scale=scale12_MSM_SM, shape=shape12_MSM_SM, lower = FALSE, log = TRUE)*exp(beta_X1_12*xdata1)}
  H01_MSM_SM<-function(x,xdata1) {-pweibull(x, scale=scale01_MSM_SM, shape=shape01_MSM_SM, lower = FALSE, log = TRUE)*exp(beta_X1_01*xdata1)}
  H02_MSM_SM<-function(x,xdata1) {-pweibull(x, scale=scale02_MSM_SM, shape=shape02_MSM_SM, lower = FALSE, log = TRUE)*exp(beta_X1_02*xdata1)}
  
  # NO clock reset
  resdist0<-function(v,xdata1) {exp(-((H02_MSM_SM(v,xdata1)-H02_MSM_SM(tti,xdata1))+(H01_MSM_SM(v,xdata1)-H01_MSM_SM(tti,xdata1))))*lambda01_MSM_SM(v,xdata1)*
      exp(-(H12_MSM_SM(tti+w_predict,xdata1,v)-H12_MSM_SM(v,xdata1,v)))}
  
  Fw_return<-function(xdata1,zdata,v_predict)
  {
    Fw1<-1-exp(-(H12_MSM_SM(tti+w_predict,xdata1,v_predict)-H12_MSM_SM(tti,xdata1,v_predict)))
    
    temp<-seq(tti,tti+w_predict,by=0.001)
    temp2<-sum(resdist0(temp,xdata1))*0.001
    
    Fw0<-1-(exp(-((H02_MSM_SM(tti+w_predict,xdata1)-H02_MSM_SM(tti,xdata1))+(H01_MSM_SM(tti+w_predict,xdata1)-H01_MSM_SM(tti,xdata1))))+temp2)
    
    return(zdata*Fw1+(1-zdata)*Fw0)
  }
  
  for(i in 1:length(zdata))
  {
    Fw<-c(Fw,Fw_return(xdata1[i],zdata[i],v_predict[i]))
  }
  return(Fw)
}

#######################################################
#Code obtianed from Blanche et al. (2015) for computation of dynamic BS and AUC
#######################################################
#Brier Score function from Blanche et al. (2015)
BS <- function(timepoints,times,status,pred,cause=1){ 
  n <- length(times)
  n_times <- length(timepoints)
  timepoints <- timepoints[order(timepoints)]
  times_names <- paste("t=",timepoints,sep="")
  # output initialisation 
  BS <- rep(NA,n_times)
  CumInci <- rep(NA,n_times)
  surv <- rep(NA,n_times)
  Stats <- matrix(NA,nrow=n_times,ncol=4)
  hit1_all <- matrix(NA,nrow=n,ncol=n_times)
  hit2_all <- matrix(NA,nrow=n,ncol=n_times)
  epsilon_i <- matrix(NA,nrow=n,ncol=n_times)
  #adds name to outputs
  names(BS) <- times_names
  names(CumInci) <- times_names
  names(surv) <- times_names
  colnames(Stats) <- c("Cases","survivor at t","Other events at t","Censored at t")
  rownames(Stats) <- times_names
  colnames(epsilon_i) <- times_names
  colnames(hit1_all) <-  times_names
  colnames(hit2_all)  <- times_names 
  # we need to order to use the ipcw() function of the pec package
  #browser()
  order_T <- order(times)
  times <-  times[order_T]
  delta  <-  status[order_T]
  pred <-  pred[order_T,,drop=FALSE]
  #compute KM weights
  weights <- ipcw(Surv(failure_time,status)~1,
                  data=data.frame(failure_time=times,status=as.numeric(delta!=0)),
                  method="marginal",times=timepoints,subjectTimes=times,subjectTimesLag=1)
  Mat_data <- cbind(times,delta,as.numeric(delta==0))
  colnames(Mat_data) <- c("T","delta","indic_Cens")
  # computate weights of cases
  Weights_cases_all <- 1/(weights$IPCW.subjectTimes*n)
  # loop on all time points
  for(t in 1:n_times){
    Cases <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]==cause)
    Controls_1 <- (Mat_data[,"T"]> timepoints[t] )
    Controls_2 <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]!=cause & Mat_data[,"delta"]!=0)  
    # compute weights
    Weights_controls_1 <- rep(1/(weights$IPCW.times[t]*n),times=n)
    Weights_cases <- Weights_cases_all
    Weights_controls_2 <- Weights_cases_all
    Weights_cases[!Cases] <- 0
    Weights_controls_1[!Controls_1] <- 0
    Weights_controls_2[!Controls_2] <- 0   
    #compute outputs
    CumInci[t] <- c(sum(Weights_cases))
    surv[t] <- c(sum(Weights_controls_1))
    Stats[t,] <- c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2)) 
    hit1_all[,t] <- (Weights_controls_1*((pred[,t])^2))*n
    hit2_all[,t] <- (Weights_cases*((1-pred[,t])^2) + Weights_controls_2*((pred[,t])^2))*n
    BS[t] <- (sum(hit1_all[,t]) +sum(hit2_all[,t]))/n
  } 
  
  out <- list(BS=BS,res=(hit1_all+hit2_all),
              CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats,timepoints=timepoints
  )
  class(out) <- "ipcwEBS"
  out 
}

PE<-function(model,LMx,w_predict,BS_dat)
{    
  AUCst.M1 <- rep(NA,length(LMx))
  BrierS.s.M1 <- rep(NA,length(LMx))
  for (s in LMx){
    # print(s)
    # Create landmark data set
    d.s<-BS_dat[,c("observed.lifetime","seen.exit")]
    d.s$Pred.s.M1<-BS_dat[,model]
    d.s<-subset(d.s,BS_dat$LM==s)  #d.s$observed.lifetime>s&
    d.s$time.s<-d.s$observed.lifetime-s
    # AUC and BS for prediction based on M1
    # estimate ROC curve and AUC
    ROC.s.M1<-timeROC(T=d.s$time.s,
                      delta=d.s$seen.exit,
                      marker=d.s$Pred.s.M1,
                      cause=1,weighting="marginal",
                      times=c(w_predict),
                      iid=TRUE)
    # estimate expected Brier score
    BS.s.M1 <- BS(timepoints=c(w_predict),
                  times=d.s$time.s,
                  status=d.s$seen.exit,
                  pred=as.matrix(d.s$Pred.s.M1),
                  cause=1)
    # save useful results (estimates, s.e. and iid decompositions)
    BrierS.s.M1[which(s==LMx)] <- BS.s.M1$BS # BS estimate
    AUCst.M1[which(s==LMx)] <- ROC.s.M1$AUC[2] # AUC estimate
  }
  return(list(BrierS.s.M1,AUCst.M1))
}

