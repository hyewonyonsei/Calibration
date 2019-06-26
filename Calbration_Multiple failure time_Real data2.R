library(survival)
library(survey)
library(dplyr)
library(MASS)

set.seed(20190613)
# Case = 629, Subject of subcohort = 983
bhs <- read.csv("bhs.csv", header=TRUE) # str(bhs)
nna <- function(data) {
  return(sum(is.na(data)))
}
apply(bhs,2,nna)
bhs2 <- bhs[-c(which(is.na(bhs$DIABRX)),which(is.na(bhs$RXHYPER))),] # Case = 627, Subject of subcohort = 978
bhs2$id <- 1:nrow(bhs2)
# length(which(is.na(bhs2$FERRITIN)==FALSE & bhs2$SUBCOH==1))
# length(which(bhs2$STRCENS==1 & bhs2$CHDCENS==1 & is.na(bhs2$FERRITIN)==FALSE))
# length(which(bhs2$STRCENS==1 & is.na(bhs2$FERRITIN)==FALSE))
# length(which(bhs2$CHDCENS==1 & is.na(bhs2$FERRITIN)==FALSE))
# coln <- cbind(1:ncol(bhs2),colnames(bhs2))
# bhs2.sub <- bhs2[which(bhs2$SUBCOH==1),]
# matrix(format(apply(bhs2.sub,2,mean), scientific=FALSE))

bhs.temp <- bhs2[,c(24,8,9,10,11,16,17,20,1:7,12:14)]
a <- length(which(is.na(bhs.temp$FERRITIN) & bhs.temp$SUBCOH==1))/length(which(bhs.temp$SUBCOH==1))
b <- length(which(is.na(bhs.temp$FERRITIN) & bhs.temp$SUBCOH==0))
minus <- sample(which(is.na(bhs.temp$FERRITIN) & bhs.temp$SUBCOH==0),round(a*b))
bhs.full <- bhs.temp[-minus,]

n.full = nrow(bhs.full)
bhs.full$case <- c(rep(0,n.full)) # case indicator
bhs.full$case <- as.numeric(bhs.full$CHDCENS==1 | bhs.full$STRCENS==1)
num.case = length(unique(bhs.full$id[bhs.full$case==1]))
# num.subc = 608*0.1

##### Stratification + Twophase #####
# Creating indicators and conduct sampling ------------------------------------;
idx.case <- bhs.full$id[bhs.full$case==1]                      # idx.case: ids for cases only 
# idx.cont <- bhs.full$id[bhs.full$case==0 & bhs.full$strt==2] # idx.cont2: ids for censored subjects in the cohort
idx.cont <- bhs.full$id[is.na(bhs.full$FERRITIN)==FALSE]
idx.sample <- sort(c(t(idx.case),t(idx.cont))) # idx.sample: ids for cch sample

# Creating cch sampling indicator (T or F)
num.in.cch = length(idx.sample) # cch sample size
bhs.full$in.cch <- F
bhs.full$in.cch[bhs.full$id %in% idx.sample] <- T

num.subc = length(idx.cont)
alpha = num.subc/(n.full-num.case)
bhs.full$weight <- c(rep(1/alpha,n.full))    # weight
bhs.full$weight[bhs.full$case==1] = 1
bhs.full$prob = c(rep(alpha,n.full))         # sampling probabilites
bhs.full$strt = c(rep(2,nrow(bhs.full)))     # strt = 1: cases, = 2: controls 
bhs.full$strt[bhs.full$case==1] = 1

dstrat <- twophase(id=list(~1,~1), strata=list(NULL,~strt), subset = ~in.cch, data=bhs.full)

##### GLM ##### 
# full <- lm(FERRITIN~﻿AGE+SMOKE+SBP+CHOL+TRIGLYCE+HAEMOGLO+BMI+DIABRX+RXHYPER+SEX, data=full)
# null <- lm(FERRITIN~1, data=full)
# full <- lm(LOGFERR~﻿AGE+SMOKE+SBP+CHOL+TRIGLYCE+HAEMOGLO+BMI+DIABRX+RXHYPER+SEX, data=bhs.full)
# null <- lm(LOGFERR~1, data=bhs.full)
# step(null, direction="both", scope=list(upper=full))
# Hmodelf <- glm(FERRITIN~SMOKE+SBP+CHOL+TRIGLYCE+HAEMOGLO+BMI+DIABRX+SEX, data=bhs.full)
glm.model <- glm(LOGFERR~﻿AGE+SBP+CHOL+HAEMOGLO+BMI+DIABRX+SEX, data=bhs.full)

##### Imputation ##### 
# using the predicted values only
bhs.full$estH <- as.numeric(predict(glm.model, type="response", newdata = bhs.full, se=F))

##### To a long form #####
idl <- rep(bhs.full$id, each=2)
timel <- as.vector(t(cbind(bhs.full$CHDTIME, bhs.full$STRTIME)))
deltal <- as.vector(t(cbind(bhs.full$CHDCENS,bhs.full$STRCENS)))
z <- bhs.full[,7:18]
ztemp <- matrix(nrow=n.full*2,ncol=ncol(z))
for (i in 1:ncol(z)) {
  ztemp[,i] <- rep(z[,i],each=2)
}
colnames(ztemp) <- c("FERRITIN","LOGFERR","AGE","SMOKE","SBP","CHOL","TRIGLYCE","HAEMOGLO","BMI","DIABRX","RXHYPER","SEX")
bhs.long <- data.frame(id = idl, time=timel, delta = deltal, ztemp)
bhs.long$case <- as.vector(t(cbind(bhs.full$CHDCENS, bhs.full$STRCENS)))
bhs.long$weight <- rep(bhs.full$weight, each=2)
bhs.long$in.cch <- rep(bhs.full$in.cch, each=2)
bhs.long$prob <- rep(bhs.full$prob, each=2)
bhs.long$strt <- rep(bhs.full$strt, each=2)
bhs.long$SUBCOH <- rep(bhs.full$SUBCOH, each=2)
bhs.long$estH <- rep(bhs.full$estH, each=2)

# Calibration
# Creating auxiliary variables by fitting a full cohort data with imputed values
# GLM: ﻿AGE+SBP+CHOL+HAEMOGLO+BMI+DIABRX+SEX
# imputed <- coxph(Surv(time, delta)~estH+AGE+SMOKE+SBP+CHOL+TRIGLYCE+HAEMOGLO+BMI+DIABRX+RXHYPER+SEX+cluster(id), data=bhs.long)
# imputed <- coxph(Surv(time, delta)~estH+TRIGLYCE+RXHYPER+cluster(id), data=bhs.long)
# imputed <- coxph(Surv(time, delta)~estH+AGE+TRIGLYCE+HAEMOGLO+DIABRX+RXHYPER+cluster(id), data=bhs.long)
imputed <- coxph(Surv(time, delta)~estH+AGE+TRIGLYCE+HAEMOGLO+RXHYPER+cluster(id), data=bhs.long)
# summary(imputed)

##### Auxiliary variable #####
# Calibration
resid <- residuals(imputed, "dfbeta")
db <- resid
colnames(db) <- paste("db",1:ncol(db),sep="")
bhs.long <- cbind(bhs.long, db)

##### Back to the short form #####
db.short <- matrix(db[,1],ncol=2,byrow=TRUE)
for (i in 2:ncol(db)) {
  db.short <- cbind(db.short, matrix(db[,i],ncol=2,byrow=TRUE))
}
db.sum <- matrix(ncol=ncol(db),nrow=nrow(bhs.full))
for (i in 1:ncol(db)) {
  db.sum[,i] <- db.short[,(i*2-1)]+db.short[,(i*2)]+1  # 1 is added for the calculation purpose
}
# apply(db.sum,2,sum)
colnames(db.sum) <- paste0("db",c(1:ncol(db)))
bhs.full <- cbind(bhs.full,db.sum)

##### Calibration #####
dstrt <- twophase(id=list(~1,~1), strata=list(NULL,~strt), subset=~in.cch, data=bhs.full)
dcal <- calibrate(dstrt, formula=make.formula(colnames(db.sum)),pop=c(`(Intercept)`=nrow(bhs.full), colSums(db.sum)),
                  calfun="raking", eps=0.0001)
calw <- rep(weights(dcal), each=2)

# when using only controls ----------------------;
# local calibration: excluding strata sampled completely;
bhs.control <- bhs.full[bhs.full$strt==2,]
dbc <- db.sum[bhs.full$strt==2,]+1
# apply(dbc,2,sum)
colnames(dbc) <- paste("dbc",1:ncol(dbc),sep="")
bhs.control <- cbind(bhs.control,dbc)
n.control <- nrow(bhs.control)
dstrtc <- twophase(id=list(~1,~1),subset=~in.cch, data=bhs.control)
dcalc <- calibrate(dstrtc,formula=make.formula(colnames(dbc)),
                   pop=c(`(Intercept)`=n.control,colSums(dbc)),
                   calfun="raking",eps=0.0001)

##### Subcohort + long form #####
cch.full <- subset(bhs.full, bhs.full$id %in% idx.sample)
cch.full.num <- nrow(cch.full)
cch.long <- subset(bhs.long, bhs.long$id %in% idx.sample)
cch.long.num <- nrow(cch.long)
cch.long$calweight = 1
cch.long[cch.long$strt==2,]$calweight <- rep(weights(dcalc), each=2)

##### Final Result #####
# GLM: ﻿AGE+SBP+CHOL+HAEMOGLO+BMI+DIABRX+SEX
# base_model <- coxph(Surv(time, delta)~estH+AGE+SMOKE+SBP+CHOL+TRIGLYCE+HAEMOGLO+BMI+DIABRX+RXHYPER+SEX+cluster(id),
#                     weights=weight, data=cch.long)
base_model <- coxph(Surv(time, delta)~estH+AGE+TRIGLYCE+HAEMOGLO+RXHYPER+cluster(id),
                    weights=weight, data=cch.long)
model_cal1 <- coxph(Surv(time, delta)~estH+AGE+TRIGLYCE+HAEMOGLO+RXHYPER+cluster(id),
                    weights=calw, data=cch.long)
model_cal2 <- coxph(Surv(time, delta)~estH+AGE+TRIGLYCE+HAEMOGLO+RXHYPER+cluster(id),
                    weights=calweight, data=cch.long)
# summary(base_model); summary(model_cal1)

base <- summary(base_model); cal1 <- summary(model_cal1); cal2 <- summary(model_cal2)
beta <- cbind(base$coefficients[,1], cal1$coefficients[,1], cal2$coefficients[,1])
colnames(beta) <- c("Traditional", "Cal_All", "Cal_Control")
se <- cbind(base$coefficients[,3], cal1$coefficients[,3], cal2$coefficients[,3])
colnames(se) <- c("Traditional", "Cal_All", "Cal_Control")
se_compare <- cbind(1, base$coefficients[,3]/cal1$coefficients[,3], base$coefficients[,3]/cal2$coefficients[,3])
se_compare