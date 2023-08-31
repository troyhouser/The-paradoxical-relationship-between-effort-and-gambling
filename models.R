############################################################
############################################################
### functions
############################################################
############################################################

# function for taking log prob ensuring no Inf return
safelog = function(x){
  x[x==0]=1e-100
  y=log(x)
  return(y)
}

## rational inattention model 2 parameters
RH0 = function(pars){
  beta = pars[1]; alpha = pars[2]
  r = cbind(dat$X2,dat$X1)
  r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/(abs(r)*beta)
  k = s2_e/dat$s2_u
  V = abs(r)/(1+(k*t))
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
REXPO = function(pars){
  beta = pars[1]; alpha = pars[2]
  r = cbind(dat$X2,dat$X1)
  r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/(abs(r)*beta)
  k = s2_e/dat$s2_u
  V = r*exp(-k*t)
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}

## hyperbolic model
H0 = function(pars){
  alpha = pars[1]; k = pars[2]
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6;r=abs(r)
  t = cbind(dat$T2,dat$T1)
  V = r/(1+(k*t))
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}

# generalized hyperbolic model
GH = function(pars){
  alpha = pars[1]; k = pars[2]; v = pars[3]
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6;r=abs(r)
  t = cbind(dat$T2,dat$T1)
  V = r/(1+(k*t)^(v/k))
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
## rational inattention model 2 parameters
RGH = function(pars){
  beta = pars[1]; alpha = pars[2]; v = pars[3]
  r = cbind(dat$X2,dat$X1);r=abs(r)
  r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/(abs(r)*beta)
  k = s2_e/dat$s2_u
  V = r/(1+(k*t)^(v/k))
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
## exponential model
EXPO = function(pars){
  alpha = pars[1]; k = pars[2]
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6;r=abs(r)
  t = cbind(dat$T2,dat$T1)
  V = r*exp(-k*t)
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}

#hyperboloid model
HBoid = function(pars){
  alpha = pars[1]; k = pars[2]; S = pars[3]
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6;r=abs(r)
  t = cbind(dat$T2,dat$T1)
  V = r/(1+(k*(t^S)))
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
RHBoid = function(pars){
  beta = pars[1]; alpha = pars[2]; S = pars[3]
  r = cbind(dat$X2,dat$X1)
  r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/(abs(r)*beta)
  k = s2_e/dat$s2_u
  V = abs(r)/(1+(k*(t^S)))
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
# power model 
P0 = function(pars){
  alpha = pars[1]; k = pars[2]; S = pars[3]
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6;r=abs(r)
  t = cbind(dat$T2,dat$T1)
  V = r - (k*(t^S))
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
RP0 = function(pars){
  beta = pars[1]; alpha = pars[2]; S = pars[3]
  r = cbind(dat$X2,dat$X1)
  r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/(abs(r)*beta)
  k = s2_e/dat$s2_u
  V = r - (k*(t^S))
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
## quasi hyperbolic model
QH = function(pars){
  d = pars[1]; b = pars[2]; alpha = pars[3]
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6;r=abs(r)
  t = cbind(dat$T2,dat$T1)
  V = r*b*(d^t)
  P = pnorm(alpha*(V[,1]-V[,2]))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}

rand = function(){
  lik = safelog(0.5)*dat$LL+safelog(0.5)*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}

get_k_for_RH0 = function(beta,alpha){
  r = cbind(dat$X2,dat$X1)
  r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/(abs(r)*beta)
  k = s2_e/dat$s2_u
  return(k[,1])
}

FUN_bootMer <- function(fit){
  return(fixef(fit)) 
}