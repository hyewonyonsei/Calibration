# @kkanghye 
# hyewonyonsei/Calibration
# 248 lines (219 sloc)  10.5 KB

library(survival)
library(survey)
library(dplyr)
library(MASS)

calibration <- function(num.pop, num.subc, k, Sigma, u2, u3, beta, lamzero, gam, theta) {
  # Data generation
  # Creating a full cohort dataset: c1, c2 - 2 covariates used in the model, 
  # c3 - auxiliary cov, associated with c2
  ##### Values arbitrarily assigned by HW #####
  u <- c(); z <- c()
  u <- runif(num.pop)
  z1 <- rbinom(num.pop,1,0.5) # Binomial
  z_cts <- mvrnorm(n=num.pop, c(u2,u3), Sigma) # Multivariate normal
  z2 <- z_cts[,1]
  z3 <- z_cts[,2] #) z2 & z3 -> correlated data
  z <- cbind(z1, z2, z3)
  p = ncol(z)
  
  ##### Time generation to adopt the MATLAB code - Extension of Clayton and Cuzick (1985) #####
  suma <- matrix(rep(0,num.pop),nrow=num.pop)
  # tt <- c()
  tt1 = (-log(1-u)*(exp(-z%*%beta)*(1/lamzero[1])))^(1/gam[1])
  tt = tt1
  if (k>1) {
    for (i in 2:k) {
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
  
  ##### Censoring #####
  cmax = 0.025
  censor <- matrix(nrow=num.pop, ncol=k)
  delta <- matrix(nrow=num.pop, ncol=k)
  time_obs <- matrix(nrow=num.pop, ncol=k)
  for (i in 1:k) {
    censor[,i] = cmax*runif(num.pop,0,1) #) Modest censoring assumption by Cai and Shen (2000)
    delta[,i] = (tt[,i]<=censor[,i])*1
    time_obs[,i] = tt[,i]*delta[,i] + censor[,i]*(1-delta[,i])
  }
  id = 1:num.pop
  colnames(z) <- paste0("c",c(rep(1:k,each=p)),".",c(rep(1:p)))
  data.full <- data.frame(id = id, time1=time_obs[,1], time2=time_obs[,2], 
                          delta1 = delta[,1], delta2 = delta[,2], 
                          c1.1 = z[,1], c1.2 = z[,2], c1.3 = z[,3],
                          c2.1 = z[,4], c2.2 = z[,5], c2.3 = z[,6])
  
  data.full$case = c(rep(0,num.pop))            # case indicator
  # data.full$case[data.full$time < cmax] = 1 cbind(data.full$case,data.full$delta)
  data.full$case = as.numeric(data.full$delta1==1 | data.full$delta2==1)
  num.case = length(unique(data.full$id[data.full$case==1]))
  num.subc = (num.pop-num.case)*0.1  #subcohort 10% -> alpha=10
  
  alpha = num.subc/(num.pop-num.case)
  data.full$subcind = c(rep(0,num.pop))         # subcohort indicator             
  data.full$weight = c(rep(1/alpha,num.pop))    # weight
  data.full$weight[data.full$case==1] = 1
  data.full$in.cch = c(rep(0,num.pop))          # case-cohort indicator
  data.full$prob = c(rep(alpha,num.pop))        # sampling probabilites 
  data.full$strt = c(rep(2,num.pop))    # strt = 1: cases, = 2: controls 
  data.full$strt[data.full$case==1] = 1
  
  ##### Stratification + Twophase #####
  # Creating indicators and conduct sampling ------------------------------------;
  idx.case <- data.full$id[data.full$case==1]                      # idx.case: ids for cases only 
  idx.cont <- data.full$id[data.full$case==0 & data.full$strt==2] # idx.cont2: ids for censored subjects in the cohort
  # sampling subcohort controls
  idx.scont2<-sample(idx.cont, size=num.subc)
  #idx.case2<-sample(idx.case,size=ceiling(num.case*0.8))
  idx.sample<-sort(c(t(idx.case),t(idx.scont2)))           # idx.sample: ids for cch sample
  #idx.sample<-sort(cbind(t(idx.case2),t(idx.scont2))) 
  
  # data.full$contind = c(rep(0,num.pop))
  # data.full$contind[idx.scont2] <- 1                           # contind: indicator for censored subjecrts in the subcohort  
  
  # Creating cch sampling indicator (T or F)
  data.full$in.cch <- F
  data.full$in.cch[idx.sample] <- T
  num.in.cch = length(idx.sample) # cch sample size
  
  dstrat <- twophase(id=list(~1,~1), strata=list(NULL,~strt), subset = ~in.cch, data=data.full)
  # head(data.sampled); str(data.sampled)
  # data.full$subcind <- as.numeric(data.full$id %in% idx.sample)
  
  ##### To a long form #####
  idl <- c(rep(1:num.pop, each=k))
  timel <- as.vector(t(time_obs))
  deltal <- as.vector(t(delta))
  ztemp <- matrix(nrow=num.pop*k, ncol=p)
  for(i in 1:p) {
    colnames(ztemp) <- paste0("z",1:p)
    ztemp[,i] <- as.vector(t(cbind(z[,c(i,i+p)])))
  }
  data.long <- data.frame(id = idl, time=timel, delta = deltal,
                          c1 = ztemp[,1], c2 = ztemp[,2], c3 = ztemp[,3])
  
  sequ <- seq(from=1, to=nrow(data.long), by=2)
  data.long$case <- c()
  for (i in sequ) {
    if (data.long$delta[i]==1 | data.long$delta[i+1]==1) {
      data.long$case[i]=1
      data.long$case[i+1]=1
    } else {
      data.long$case[i]=0
      data.long$case[i+1]=0
    }
  }
  # data.long$subcind <- as.vector(t(data.full$subcind))
  data.long$weight <- rep(data.full$weight, each=k)
  data.long$in.cch <- rep(data.full$in.cch, each=k)
  data.long$prob <- rep(data.full$prob, each=k)
  data.long$strt <- rep(data.full$strt, each=k)
  
  data.long = data.frame(id = idl, time=timel, delta = deltal,
                         c1 = ztemp[,1], c2 = ztemp[,2], c3 = ztemp[,3],
                         subcind=data.long$subcind, weight = data.long$weight, in.cch = data.long$in.cch,
                         prob = data.long$prob, case = data.long$case, strt = data.long$strt)
  
  ##### Case-cohort + long form #####
  cch.full <- subset(data.full, data.full$id %in% idx.sample)
  cch.full.num <- nrow(cch.full)
  cch.long <- subset(data.long, data.long$id %in% idx.sample)
  cch.long.num <- nrow(cch.long)
  
  ##### Imputation #####
  Hmodel <- glm(c2 ~ c3, weights = weight, data=cch.long)
  # summary(Hmodel)
  data.long$estH <- as.numeric(predict(Hmodel, type="response", newdata = data.long, se=F)) # using the predicted values only
  # length(which(data.long$in.cch==T))
  data.long$estH[data.long$in.cch==T] <- data.long$c2[data.long$in.cch==T]  # in case using observed values+predicted values
  # length(which(data.long$estH==data.long$c2))
  
  imputed <- list()
  imputed <- coxph(Surv(time, delta) ~ c1 + estH + c3 + cluster(id), data=data.long) # creating auxiliary variables by fitting a full cohort data with imputed values
  # summary(imputed)
  
  ##### Auxiliary variable #####
  # Calibrating
  resid <- residuals(imputed, "dfbeta")
  # invD = qr.solve(imputed$var) 
  # db = resid%*%invD #+1            # 1 was added to avoid a computation problem
  colnames(db) <- paste("db",1:ncol(db),sep="")
  # head(db)
  data.long.db <- cbind(data.long,db) #head(data.long.db)
  
  ##### Back to the short form #####
  db.short <- cbind(matrix(data.long.db$db1,ncol=k,byrow=TRUE),matrix(data.long.db$db2,ncol=k,byrow=TRUE),
                    matrix(data.long.db$db3,ncol=k,byrow=TRUE))
  colnames(db.short) <- paste0("db",c(rep(1:ncol(db),each=k)),".",1:k)
  
  db.sum <- cbind(db.short[,1]+db.short[,2],db.short[,3]+db.short[,4],db.short[,5]+db.short[,6])
  colnames(db.sum) <- paste0("db",1:3)
  sum.db <- cbind(data.full, db.sum)
  # colSums(db.sum)
  # head(data.db); str(data.db)
  
  dataxx <- subset(sum.db, sum.db$case==0)
  ss <- colSums(dataxx[,c(17,18,19)])
  xxnum <- nrow(dataxx)
  dataxx$db1 <- dataxx$db1-rep(ss[1]/xxnum,xxnum)
  dataxx$db2 <- dataxx$db2-rep(ss[2]/xxnum,xxnum)
  dataxx$db3 <- dataxx$db3-rep(ss[3]/xxnum,xxnum)
  
  ##### Calibration #####
  dstrt <- twophase(id=list(~1,~1), strata=list(NULL,~strt), subset=~in.cch, data=sum.db)
  dcal <- calibrate(dstrt, formula=make.formula(colnames(db.sum)), pop=c(`(Intercept)`=0, colSums(db.sum)), calfun="raking", eps=0.0001)
  data.dstrt <- cbind(model.frame(dstrt), calw=weights(dcal))
  # head(weights(dcal)); head(sample.f)
  # calw <- rep(weights(dcal), each=k)
  # cch.final <- cbind(cch.long, calw)
  cch.final <- merge(cch.long, data.dstrt[,c(1,20)])
  
  designx <- twophase(id=list(~1,~1), strata=list(NULL,~strt), subset=~in.cch, data=dataxx)
  calx <- calibrate(designx, formula=make.formula(colnames(db.sum)),
                    pop=c(`(Intercept)`=0, colSums(db.sum)), calfun="raking",
                    eps=0.00001)
  datax <- cbind(model.frame(designx), calwx=weights(calx))
  datax2 <- cbind(subset(sum.db, sum.db$case==1), calwx=rep(1, nrow(subset(cch.full, cch.full$case==1))))
  data.cx <- rbind(datax, datax2)
  data.mx <- as.data.frame(merge(cch.long, data.cx[,c(1,20)]))
  # which(data.mx$calwx<0)
 
  ##### Fitting to Cox regression for multiple failure time data #####
  base_model <- coxph(Surv(time, delta) ~ c1+c2+c3, weights=weight, data=data.long.f)
  # summary(base_model)
 
  model_cal <- coxph(Surv(time, delta) ~ c1+c2+c3+cluster(id), weights=calw, data=cch.final)
  model_cal2 <- coxph(Surv(time, delta) ~ c1+c2+c3+cluster(id), weights=calwx, data=data.mx)
  # summary(model_cal)
  
  compare <- c(beta, base_model$coefficients, model_cal$coefficients)
  return(compare)
  # return(compare2)
}

# set.seed(20190419)
# num.pop = 1000; num.subc = 2; k=2
# Sigma = matrix(c(1,0.5,0.5,1), nrow=2)
# u2 = 0; u3 = 0
# beta = c(0.5, 1, 1.2)
# lamzero = c(1, 1.5)
# gam = c(rep(1,k))
# theta = 0.25 

result <- function(data) {
  beta.mean <- matrix(apply(data,2,mean), ncol=3)
  beta.var <- matrix(apply(data,2,var), ncol=3)
  colnames(beta.mean) <- c("True", "Traditional", "Calibrated")
  rownames(beta.mean) <- c("beta1", "beta2", "beta3")
  colnames(beta.var) <- c("True", "Traditional", "Calibrated")
  rownames(beta.var) <- c("beta1", "beta2", "beta3")
  return(list(mean=beta.mean,variance=beta.var))
}

Sigma = matrix(c(1,0.7,0.7,1), nrow=2); beta = c(0.5, 1, 1.2)
lamzero = c(1, 1.5); gam = c(rep(1,2)); theta = 0.25 

set.seed(2019)
beta.compare <- matrix(rep(0,1000*9), ncol=9)
time.elapse <- system.time(
for (i in 1:1000) {
  beta.compare[i,] <- calibration(1000,100,2,Sigma,0,0,beta,lamzero,gam,0.25)
}
)

result(beta.compare)
