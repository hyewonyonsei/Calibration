library(survival)
library(survey)
library(dplyr)
library(MASS)

# Data generation
# Creating a full cohort dataset: c1, c2 - 2 covariates used in the model, 
# c3 - auxiliary cov, associated with c2
#) Values arbitrarily assigned by HW
num.pop = 1000; num.subc = 2
Sigma = matrix(c(1,0.5,0.5,1), nrow=2)
u2 = 0; u3 = 0
beta = c(0.5, 1, 1.2)
lamzero = c(1, 1.5)
gam = c(rep(1,num.subc))

u <- c(); z <- c()
u <- runif(num.pop)
z1 <- rbinom(num.pop,1,0.5) # Binomial
z_cts <- mvrnorm(n=num.pop, c(u2,u3), Sigma) # Multivariate normal
z2 <- z_cts[,1]
z3 <- z_cts[,2] #) z2 & z3 -> correlated data
z <- cbind(z1, z2, z3)
p = ncol(z)

#) Time generation to adopt the MATLAB code
theta = 0.25 

#) Extension of Clayton and Cuzick (1985)
suma <- matrix(rep(0,num.pop),nrow=num.pop)
tt1 = (-log(1-u)*(exp(-z%*%beta)*(1/lamzero[1])))^(1/gam[1])
tt = tt1
if (num.subc>1) {
  for (i in 2:num.subc) {
    ui <- runif(num.pop)
    z1i <- rbinom(num.pop,1,0.5) # bernoulli
    z_ctsi <- mvrnorm(n=num.pop, c(u2,u3), Sigma) # multivariate normal
    zi <- cbind(z1i, z_ctsi)
    suma <- suma + exp((tt[,i-1]^gam[i])*((exp(z[,((i-2)*p+1):((i-1)*p)]%*%beta))*(lamzero[1]/theta)))
    b1 <- log((i-1)-suma+(suma-(i-2))*((1-ui)^(-(theta+i-1)^(-1))))
    tttem = ((theta/lamzero[i])*b1*(exp(-zi%*%beta)))^(1/gam[i])
    tt <- cbind(tt,tttem)
    u <- cbind(u, ui)
    z <- cbind(z, zi)
  }
}

#) Censoring
# cmax = 2 
censor <- matrix(nrow=num.pop, ncol=num.subc)
delta <- matrix(nrow=num.pop, ncol=num.subc)
time_obs <- matrix(nrow=num.pop, ncol=num.subc)
for (i in 1:num.subc) {
  censor[,i] = runif(num.pop,0,1) #) Modest censoring assumption by Cai and Shen (2000)
  delta[,i] = (tt[,i]<=censor[,i])*1
  time_obs[,i] = tt[,i]*delta[,i] + censor[,i]*(1-delta[,i])
}
id = 1:num.pop
colnames(z) <- paste0("c",c(rep(1:num.subc,each=p)),".",c(rep(1:p)))
data.full <- data.frame(id = id, time1=time_obs[,1], time2=time_obs[,2], 
                        delta1 = delta[,1], delta2 = delta[,2], 
                        c1.1 = z[,1], c1.2 = z[,2], c1.3 = z[,3],
                        c2.1 = z[,4], c2.2 = z[,5], c2.3 = z[,6])
head(data.full)

data.full$case = c(rep(0,num.pop))            # case indicator
#data.full$case[data.full$time < cmax] = 1 cbind(data.full$case,data.full$delta)
data.full$case = as.numeric(data.full$delta1==1 | data.full$delta2==1)

num.case = length(unique(data.full$id[data.full$case==1]))
alpha.long = num.subc/(num.pop-num.case)
data.full$subcind = c(rep(0,num.pop))         # subcohort indicator             
data.full$weight = c(rep(1/alpha.long,num.pop))    # weight
data.full$in.cch = c(rep(0,num.pop))          # case-cohort indicator
data.full$prob = c(rep(alpha.long,num.pop))        # sampling probabilites 
data.full$strt = c(rep(2,num.pop))    # strt = 1: cases, = 2: controls 
data.full$strt[data.full$case==1] = 1

#) Stratification + Twophase
# Creating indicators and conduct sampling ------------------------------------;
idx.case <- data.full$id[data.full$case==1]                      # idx.case: ids for cases only 
idx.cont2 <- data.full$id[data.full$case==0 & data.full$strt==2] # idx.cont2: ids for censored subjects in the cohort
# sampling subcohort controls
idx.scont2 <- sample(idx.cont2,size=num.subc)
idx.case2 <- sample(idx.case,size=ceiling(num.case*0.8))
# idx.sample <- sort(cbind(t(idx.case),t(idx.scont2)))           # idx.sample: ids for cch sample
idx.sample<-sort(cbind(t(idx.case2),t(idx.scont2))) 

data.full$contind = c(rep(0,num.pop))
data.full$contind[idx.scont2] <- 1                           # contind: indicator for censored subjecrts in the subcohort  
# Creating cch sampling indicator (T or F)
data.full$in.cch <- F                             
data.full$in.cch[idx.sample]<-T 
# cch sample size
num.in.cch = length(idx.sample)
dstrat <- twophase(id=list(~1,~1), strata=list(NULL,~strt), subset = ~in.cch, data=data.full)
data.sampled <- model.frame(dstrat)
head(data.sampled)

#) To a long form
idl <- c(rep(1:num.pop, each=num.subc))
timel <- as.vector(t(time_obs))
deltal <- as.vector(t(delta))
ztemp <- matrix(nrow=num.pop*num.subc, ncol=p)
for(i in 1:p) {
  colnames(ztemp) <- paste0("z",1:p)
  ztemp[,i] <- as.vector(t(cbind(z[,c(i,i+p)])))
}
data.long <- data.frame(id = idl, time=timel, delta = deltal,
                        c1 = ztemp[,1], c2 = ztemp[,2], c3 = ztemp[,3])

seq1 <- seq(from=1, to=nrow(data.long), by=2)
seq2 <- seq(from=2, to=nrow(data.long), by=2)
data.long$case <- c()
for (i in seq1) {
  if (data.long$delta[i]==1 || data.long$delta[i+1]==1) {
    data.long$case[i]=1
    data.long$case[i+1]=1
  } else {
    data.long$case[i]=0
    data.long$case[i+1]=0
  }
}
num.case = length(unique(data.long$id[data.long$case==1]))
alpha.long = num.subc/(num.pop-num.case)
data.long$subcind = c(rep(0,num.pop*num.subc))         # subcohort indicator             
data.long$weight = c(rep(1/alpha.long,num.pop*num.subc))    # weight
data.long$in.cch = c(rep(0,num.pop*num.subc))          # case-cohort indicator
data.long$prob = c(rep(alpha.long,num.pop*num.subc))        # sampling probabilites 
data.long$strt = c(rep(2,num.pop*num.subc))    # strt = 1: cases, = 2: controls 
data.long$strt[data.long$case==1] = 1

data.long = data.frame(id = idl, time=timel, delta = deltal,
                       c1 = ztemp[,1], c2 = ztemp[,2], c3 = ztemp[,3],
                       subcind=data.long$subcind, weight = data.long$weight, in.cch = data.long$in.cch,
                       prob = data.long$prob, case = data.long$case, strt = data.long$strt)

#) Subcohort
sub.full <- subset(data.full, data.full$id %in% data.sampled$id==1)
sub.full.num <- nrow(sub.full)
sub.idl <- c(rep(sub.full$id, each=num.subc))
sub.timel <- as.vector(t(cbind(sub.full$time1, sub.full$time2)))
sub.deltal <- as.vector(t(cbind(sub.full$delta1, sub.full$delta2)))
sub.ztemp <- matrix(nrow=sub.full.num*num.subc, ncol=p)
for(i in 1:p) {
  colnames(sub.ztemp) <- paste0("z",1:p)
  sub.ztemp[,i] <- as.vector(t(sub.full[,c((i+5),(i+p+5))]))
}
sub.long <- data.frame(id = sub.idl, time = sub.timel, delta = sub.deltal,
                        c1 = sub.ztemp[,1], c2 = sub.ztemp[,2], c3 = sub.ztemp[,3])

seq <- seq(from=1, to=nrow(sub.long), by=2)
sub.long$case <- c()
for (i in seq) {
  if (sub.long$delta[i]==1 | sub.long$delta[i+1]==1) {
    sub.long$case[i]=1
    sub.long$case[i+1]=1
  } else {
    sub.long$case[i]=0
    sub.long$case[i+1]=0
  }
}
sub.num.case = length(unique(sub.long$id[sub.long$case==1]))
alpha.sub = num.subc/(num.pop-sub.num.case)
sub.long$subcind = c(rep(0,sub.full.num*num.subc))         # subcohort indicator             
sub.long$weight = c(rep(1/alpha.long,sub.full.num*num.subc))    # weight
sub.long$in.cch = c(rep(0,sub.full.num*num.subc))          # case-cohort indicator
sub.long$prob = c(rep(alpha.long,sub.full.num*num.subc))        # sampling probabilites 
sub.long$strt = c(rep(2,sub.full.num*num.subc))    # strt = 1: cases, = 2: controls 
sub.long$strt[sub.long$case==1] = 1

sub.long <- data.frame(id = sub.idl, time = sub.timel, delta = sub.deltal,
                       c1 = sub.ztemp[,1], c2 = sub.ztemp[,2], c3 = sub.ztemp[,3],
                       subcind=sub.long$subcind, weight = sub.long$weight, in.cch = sub.long$in.cch,
                       prob = sub.long$prob, case = sub.long$case, strt = sub.long$strt)

#) Imputation
sub.w <- rep(weights(dstrat), each=2)
Hmodel <- glm(c2 ~ c3 + c1 + time, weight=sub.w, data=sub.long)
data.long$estH <- as.numeric(predict(Hmodel, type="response", newdata = data.long, se=F)) # using the predicted values only
data.long$estH[data.long$in.cch==T] <- data.long$c2[data.long$in.cch==T]  # in case using observed values+predicted values
# surv <- Surv(data.long$time,data.long$delta)
# temp <- array(data=NA)
# temp <- as.matrix(cbind(data.long$c1, data.long$estH, data.long$c3))
imputed <- list()
imputed <- coxph(Surv(time, delta) ~ c1 + estH, data=data.long) # creating auxiliary variables by fitting a full cohort data with imputed values
# summary(imputed)

#) Auxiliary variable
#################Calibrating######################
resid <- residuals(imputed, "dfbeta")
invD = qr.solve(imputed$var) #) To Check!
db = resid%*%invD #+1            # 1 was added to avoid a computation problem
colnames(db) <- paste("db",1:ncol(db),sep="")
data.long.db <- cbind(data.long,db) #head(data.long.db)

#) Back to the short form
db.short <- cbind(matrix(data.long.db$db1, ncol=num.subc), matrix(data.long.db$db2, ncol=num.subc))
colnames(db.short) <- paste0("db",c(rep(1:2,each=2)),".",1:2)
db.sum <- cbind(db.short[,1]+db.short[,2],db.short[,3]+db.short[,4])
colnames(db.sum) <- paste0("db",1:2)
data.db <- cbind(data.full, db.sum)
head(data.db)
str(data.db)

#) Calibration
dstrt <- twophase(id=list(~1,~1), strata=list(NULL,~strt), subset=~in.cch, data=data.db)
dcal <- calibrate(dstrt, formula=make.formula(colnames(db.sum)), pop=c(`(Intercept)`=num.pop, colSums(db.sum)), calfun="linear", eps=0.0001)
sample.f <- model.frame(dstrt)
calw <- rep(weights(dcal), each=num.subc)
sub.f <- subset(data.long, data.long$id %in% sample.f$id==1)
data.long.f <- cbind(sub.f, calw) #str(data.long.f)

#) Fitting to Cox regression for multiple failure time data
which(calw<0)
model_cal <- coxph(Surv(time, delta) ~ c1+c2+c3, weights=calw, data=data.long.f)

#############db for variance estimation##############  
data.fit = data.full[data.full$in.cch==T,] 
v.surv = Surv(data.fit$time,data.fit$delta)
v.temp=as.matrix(cbind(data.fit$c1, data.fit$c2))  
v.fit = ahaz(v.surv,v.temp,weights = data.fit$weight)
#v.resid = array(data=NA)
v.resid <- residuals(v.fit)
#v.invD = array(data=NA)
v.invD <- qr.solve(v.fit$D)
#v.db = array(data=NA)
v.db = v.resid%*%v.invD 
colnames(v.db)<-paste("db",1:ncol(db),sep="")
#data.v.db = array(data=NA)
data.v.db = cbind(data.fit,v.db) #head(data.db)
v.num.cont = length(data.v.db$delta[data.v.db$delta==0])
#v.est = 0
#vv.db = cbind(data.v.db$db1[data.v.db$contind==1],data.v.db$db2[data.v.db$contind==1])
vv.db = cbind(data.v.db$db1,data.v.db$db2) #colSums(vv.db)
v.est = t(vv.db)%*%vv.db
v.est2 = v.invD%*%v.fit$B%*%v.invD + (1-alpha)*t(vv.db[data.v.db$contind==1,])%*%vv.db[data.v.db$contind==1,]

####################################
#w = data.fit$weight
#dstrt = list()
#dstrt<-twophase(id=list(~1,~1),strata=list(NULL,~strt),subset=~in.cch,fpc=list(NULL,~fpc),data=data.db)
dstrt<-twophase(id=list(~1,~1),strata=list(NULL,~strt),subset=~in.cch,data=data.db)
#dcal = list()
#dcal<-calibrate(dstrt,formula=make.formula(colnames(db)),pop=c(`(Intercept)`=1000,colSums(db)),calfun="linear",eps=0.0001)
dcal<-calibrate(dstrt,formula=make.formula(colnames(db)),pop=c(`(Intercept)`=num.pop,colSums(db)),calfun="raking",eps=0.0001)

# Using SRS for generated data - Once status=1 (117 patients), selected + 20% of remainder (80*0.2=16 patients) = 113 are selected
set.seed(2019)
stid.full <- unique(data.full$id[data.full$delta==1])
nid.full <- unique(data.full$id[!(data.full$id %in% stid.full)])
subid.full <- unique(c(sample(nid.full, length(nid.full)*0.2),stid.full))
data.sub <- data.full[data.full$id %in% subid.full,]
# str(data.sub); head(data.sub)
n.full <- length(unique(data.full$id)); subn.full <- length(unique(subid.full))