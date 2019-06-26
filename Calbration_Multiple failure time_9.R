library(survival)
library(survey)
library(dplyr)
library(MASS)

set.seed(20190523)
datagen <- function(num.pop, num.subc, k, Sigma, u2, u3, beta, lamzero, gam, theta, cmax) {
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
  p = 3
  
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
  censor <- matrix(nrow=num.pop, ncol=k)
  delta <- matrix(nrow=num.pop, ncol=k)
  time_obs <- matrix(nrow=num.pop, ncol=k)
  for (i in 1:k) {
    censor[,i] = cmax[i]*runif(num.pop,0,0.1) #) Modest censoring assumption by Cai and Shen (2000)
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
  data.long$weight <- rep(data.full$weight, each=k)
  data.long$in.cch <- rep(data.full$in.cch, each=k)
  data.long$prob <- rep(data.full$prob, each=k)
  data.long$strt <- rep(data.full$strt, each=k)
  
  # idl <- c(rep(data.full$weight, each=k))
  data.long = data.frame(id = idl, time=timel, delta = deltal,
                         c1 = ztemp[,1], c2 = ztemp[,2], c3 = ztemp[,3],
                         weight = data.long$weight, in.cch = data.long$in.cch,
                         prob = data.long$prob, case = data.long$case, strt = data.long$strt)
  
  ##### Subcohort + long form #####0
  cch.full <- subset(data.full, data.full$id %in% idx.sample)
  cch.full.num <- nrow(cch.full)
  cch.long <- subset(data.long, data.long$id %in% idx.sample)
  cch.long.num <- nrow(cch.long)
  return(list(data.full, data.long, cch.full, cch.long)) 
}

calibration <- function(data, op1, op2, op3, op4) {
  data.full <- data[[1]]
  data.long <- data[[2]]
  cch.full <- data[[3]]
  cch.long <- data[[4]]
  p=3; k=2
  ccc <- data[[4]]
  nrow(ccc[ccc$strt==2,])
  ##### Imputation #####
  if (op1==1) {glmdata=cch.long} else if (op1==2) {glmdata=cch.long[cch.long$strt==2,]}
  if (op2==1) {Hmodel <- glm(c2 ~ c3, weight=weight, data=glmdata)} 
  if (op2==2) {Hmodel <- glm(c2 ~ c3+time, weight=weight, data=glmdata)}
  # summary(Hmodel)
  data.long$estH <- as.numeric(predict(Hmodel, type="response", newdata = data.long, se=F)) # using the predicted values only
  # length(which(data.long$in.cch==T))
  if (op3==1) {data.long$estH[data.long$in.cch==T] <- data.long$c2[data.long$in.cch==T]}  # in case using observed values+predicted values}
  # length(which(data.long$estH==data.long$c2))
  
  imputed <- list()
  imputed <- coxph(Surv(time, delta) ~ c1 + estH + c3 + cluster(id), data=data.long) # creating auxiliary variables by fitting a full cohort data with imputed values
  # summary(imputed)
  
  ##### Auxiliary variable #####
  # Calibrating
  resid <- residuals(imputed, "dfbeta")
  #invD = qr.solve(imputed$var) #) To Check!
  #db = resid%*%invD #+1            # 1 was added to avoid a computation problem
  db <- resid + 1
  colnames(db) <- paste("db",1:ncol(db),sep="")
  data.long.db <- cbind(data.long,db) #head(data.long.db)
  
  ##### Back to the short form #####
  db.short <- cbind(matrix(data.long.db$db1,ncol=k,byrow=TRUE),matrix(data.long.db$db2,ncol=k,byrow=TRUE),
                    matrix(data.long.db$db3,ncol=k,byrow=TRUE))
  colnames(db.short) <- paste0("db",c(rep(1:ncol(db),each=k)),".",1:k)
  data.short.db <- cbind(data.full, db.short)
  
  db.sum <- cbind(data.short.db$db1.1+data.short.db$db1.2,
                  data.short.db$db2.1+data.short.db$db2.2,
                  data.short.db$db3.1+data.short.db$db3.2)
  colnames(db.sum) <- paste0("dbs",1:p)
  sum.db <- cbind(data.short.db, db.sum)
  # colSums(db.sum)
  # head(data.db); str(data.db)
  
  ##### Calibration #####
  if (op4==1) {
    dstrt <- twophase(id=list(~1,~1), strata=list(NULL,~strt), subset=~in.cch, data=sum.db)
    dcal <- calibrate(dstrt, formula=make.formula(colnames(db.sum)), 
                      pop=c(`(Intercept)`=nrow(data.full), colSums(db.sum)), 
                      calfun="raking", eps=0.0001)
    calw <- rep(weights(dcal), each=k)
  } else if (op4==2) {
    # when using only controls ----------------------;
    # local calibration: excluding strata sampled completely;
    datac = data.full[data.full$strt==2,]
    # #imputed.ctl = ahaz(Surv(data.ctl$time,data.ctl$delta), as.matrix(cbind(data.ctl$c1, data.ctl$estH))) # creating auxiliary variables by fitting a full cohort data with imputed values
    # #resid.ctl <- residuals(imputed.ctl)
    # #invD.ctl <- qr.solve(imputed.ctl$D)
    # #db.ctl = resid.ctl%*%invD.ctl +1            # 1 was added to avoid a computation problem
    dbc = db.sum[data.full$strt==2,]           # 1 was added to avoid a computation problem
    colnames(dbc) <- paste("dbc",1:ncol(dbc),sep="")
    data.dbc <- cbind(datac, dbc)
    
    dstrtc <- twophase(id=list(~1,~1),subset=~in.cch, data=data.dbc)
    # dcalc <-calibrate(dstrtc,formula=make.formula(colnames(db)),pop=c(`(Intercept)`=1000,colSums(db)),calfun="linear",eps=0.0001)
    # dcalc <-calibrate(dstrtc,formula=make.formula(colnames(db.ctl)),pop=c(`(Intercept)`=num.pop-num.case,colSums(db.ctl)),calfun="raking",eps=0.0001)
    ss <- nrow(data.full[data.full$strt==2,])
    dcalc <- calibrate(dstrtc,formula=make.formula(colnames(dbc)),pop=c(`(Intercept)`=ss,colSums(dbc)),
                       calfun="raking",eps=0.0001)
    cch.long$calweight = 1
    cch.long[cch.long$strt==2,]$calweight=rep(weights(dcalc), each=k)
  }
  
  base_model <- coxph(Surv(time, delta) ~ c1+c2+c3+cluster(id), weights=weight, data=cch.long)
  # summary(base_model)
  
  if (op4==1) {model_cal <- coxph(Surv(time, delta) ~ c1+c2+c3+cluster(id), weights=calw, data=cch.long)}
  if (op4==2) {model_cal <- coxph(Surv(time, delta) ~ c1+c2+c3+cluster(id), weights=calweight, data=cch.long)}
  # summary(model_cal)
  
  compare <- c(beta, base_model$coefficients, model_cal$coefficients)
  # colnames(compare) <- c("beta0", "Traditional", "Calibrated") 
  return(compare)
}

result <- function(data) {
  beta.mean <- matrix(apply(data,2,mean), ncol=3)
  beta.var <- matrix(apply(data,2,var), ncol=3)
  colnames(beta.mean) <- c("True", "Traditional", "Calibrated")
  rownames(beta.mean) <- c("beta1", "beta2", "beta3")
  colnames(beta.var) <- c("True", "Traditional", "Calibrated")
  rownames(beta.var) <- c("beta1", "beta2", "beta3")
  return(list(mean=beta.mean,variance=beta.var))
}

# rm(list="result1.1"); rm(list=ls())

# For all simulation
# num.pop = 1000; num.subc = 100; k=2
# Sigma = matrix(c(1,0.8,0.8,1), nrow=2); u2 = 0; u3 = 0
# lamzero = c(1, 1.5); gam = c(rep(1, 2)); theta = 0.25
# beta = c(0.5, 1, 1.2)
Sigma = matrix(c(1,0.8,0.8,1), nrow=2); lamzero = c(1, 1.5); gam = c(rep(1, 2))
beta = c(0.5, 1, 1.2)
nsim = 2000

# Option 1 op1=1 if glm is conducted for all subjects; =2 if glm is conducted only for control
# Option 2 op2=1 if only c3 is used in glm; =2 if time is added to glm
# Option 3 op3=1 if observed c2 + predicted c2 is used; =2 if predicted c2 is used
# Option 4 op4=1 if all data are used for calibration; =2 if only control is used for calibration

# Simulation 1: glm, only c3, predicted c2 + two calibration options (option 4) 
sim1.1 <- matrix(rep(0,nsim*9), ncol=9)
sim1.2 <- matrix(rep(0,nsim*9), ncol=9)
time.elapse <- system.time(
  for (i in 1:nsim) {
    data <- datagen(1000,100,2,Sigma,0,0,beta,lamzero,gam,0.25,0.35)
    sim1.1[i,] <- calibration(data,1,1,2,1)
    sim1.2[i,] <- calibration(data,1,1,2,2)
  }
)
cat("sim1.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim1.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim1.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim1.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)

# Compare the results with various combinations of paremeter choices
# See the result for both op4=1 and 2 when changing a parameter
# Simulation 2,3: censoring 80, subcohort 100 + rho 0.5 / 0.9
sim2.1 <- matrix(rep(0,nsim*9), ncol=9); sim2.2 <- matrix(rep(0,nsim*9), ncol=9)
sim3.1 <- matrix(rep(0,nsim*9), ncol=9); sim3.2 <- matrix(rep(0,nsim*9), ncol=9)
Sigma2 = matrix(c(1,0.5,0.5,1), nrow=2); Sigma3 = matrix(c(1,0.9,0.9,1), nrow=2)
time.elapse <- system.time(
  for (i in 1:nsim) {
    data2 <- datagen(1000,100,2,Sigma2,0,0,beta,lamzero,gam,0.25,0.35)
    data3 <- datagen(1000,100,2,Sigma3,0,0,beta,lamzero,gam,0.25,0.35)
    sim2.1[i,] <- calibration(data2,1,1,2,1)
    sim2.2[i,] <- calibration(data3,1,1,2,1)
    sim3.1[i,] <- calibration(data2,1,1,2,2)
    sim3.2[i,] <- calibration(data3,1,1,2,2)
  }
)
cat("sim2.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim2.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim2.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim3.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim3.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim2.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim3.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim3.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)

# Simulation 4,5,6: censoring 80 + sub cohort 100 -> 200 + rho 0.8 / 0.5 / 0.9
sim4.1 <- matrix(rep(0,nsim*9), ncol=9); sim4.2 <- matrix(rep(0,nsim*9), ncol=9)
sim5.1 <- matrix(rep(0,nsim*9), ncol=9); sim5.2 <- matrix(rep(0,nsim*9), ncol=9)
sim6.1 <- matrix(rep(0,nsim*9), ncol=9); sim6.2 <- matrix(rep(0,nsim*9), ncol=9)
time.elapse <- system.time(
  for (i in 1:nsim) {
    data4 <- datagen(1000,200,2,Sigma,0,0,beta,lamzero,gam,0.25,0.35)
    data5 <- datagen(1000,200,2,Sigma2,0,0,beta,lamzero,gam,0.25,0.35)
    data6 <- datagen(1000,200,2,Sigma3,0,0,beta,lamzero,gam,0.25,0.35)
    sim4.1[i,] <- calibration(data4,1,1,2,1)
    sim4.2[i,] <- calibration(data4,1,1,2,2)
    sim5.1[i,] <- calibration(data5,1,1,2,1)
    sim5.2[i,] <- calibration(data5,1,1,2,2)
    sim6.1[i,] <- calibration(data6,1,1,2,1)
    sim6.2[i,] <- calibration(data6,1,1,2,2)
  }
)
cat("sim4.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim4.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim4.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim4.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim5.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim5.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim5.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim5.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim6.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim6.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim6.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim6.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)

# Simulation 7,8,9: censoring 80 (cmax=0.35) -> 90 (cmax=0.12) + sub cohort 100, rho 0.8 / 0.5 / 0.9
sim7.1 <- matrix(rep(0,nsim*9), ncol=9); sim7.2 <- matrix(rep(0,nsim*9), ncol=9)
sim8.1 <- matrix(rep(0,nsim*9), ncol=9); sim8.2 <- matrix(rep(0,nsim*9), ncol=9)
sim9.1 <- matrix(rep(0,nsim*9), ncol=9); sim9.2 <- matrix(rep(0,nsim*9), ncol=9)
time.elapse <- system.time(
  for (i in 1:nsim) {
    data7 <- datagen(1000,200,2,Sigma,0,0,beta,lamzero,gam,0.25,0.35)
    data8 <- datagen(1000,200,2,Sigma2,0,0,beta,lamzero,gam,0.25,0.35)
    data9 <- datagen(1000,200,2,Sigma3,0,0,beta,lamzero,gam,0.25,0.35)
    sim7.1[i,] <- calibration(data7,1,1,2,1)
    sim7.2[i,] <- calibration(data7,1,1,2,2)
    sim8.1[i,] <- calibration(data8,1,1,2,1)
    sim8.2[i,] <- calibration(data8,1,1,2,2)
    sim9.1[i,] <- calibration(data9,1,1,2,1)
    sim9.2[i,] <- calibration(data9,1,1,2,2)
  }
)
cat("sim7.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim7.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim7.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim7.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim8.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim8.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim8.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim8.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim9.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim9.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim9.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim9.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)

# Simulation 10,11,12: censoring 80 (cmax=0.35) -> 90 (cmax=0.12) + sub cohort 100 -> 200, rho 0.8 / 0.5 / 0.9
sim10.1 <- matrix(rep(0,nsim*9), ncol=9); sim10.2 <- matrix(rep(0,nsim*9), ncol=9)
sim11.1 <- matrix(rep(0,nsim*9), ncol=9); sim11.2 <- matrix(rep(0,nsim*9), ncol=9)
sim12.1 <- matrix(rep(0,nsim*9), ncol=9); sim12.2 <- matrix(rep(0,nsim*9), ncol=9)
time.elapse <- system.time(
  for (i in 1:nsim) {
    data10 <- datagen(1000,200,2,Sigma,0,0,beta,lamzero,gam,0.25,0.35)
    data11 <- datagen(1000,200,2,Sigma2,0,0,beta,lamzero,gam,0.25,0.35)
    data12 <- datagen(1000,200,2,Sigma3,0,0,beta,lamzero,gam,0.25,0.35)
    sim10.1[i,] <- calibration(data10,1,1,2,1)
    sim10.2[i,] <- calibration(data10,1,1,2,2)
    sim11.1[i,] <- calibration(data11,1,1,2,1)
    sim11.2[i,] <- calibration(data11,1,1,2,2)
    sim12.1[i,] <- calibration(data12,1,1,2,1)
    sim12.2[i,] <- calibration(data12,1,1,2,2)
  }
)
cat("sim10.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim10.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim10.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim10.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim11.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim11.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim11.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim11.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim12.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim12.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim12.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim12.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)

# Simulation 20: Simulation 1 + theta 1.5 (paper)
sim20.1 <- matrix(rep(0,nsim*9), ncol=9)
sim20.2 <- matrix(rep(0,nsim*9), ncol=9)
time.elapse <- system.time(
  for (i in 1:nsim) {
    data20 <- datagen(1000,100,2,Sigma,0,0,beta,lamzero,gam,1.5,0.35)
    sim20.1[i,] <- calibration(data20,1,1,2,1)
    sim20.2[i,] <- calibration(data20,1,1,2,1)
  }
)
cat("sim20.1","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim20.1),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)
cat("sim20.2","\n",file="C:/Users/absol/Desktop/sim.txt",append=TRUE);capture.output(result(sim20.2),file="C:/Users/absol/Desktop/sim.txt",append=TRUE)