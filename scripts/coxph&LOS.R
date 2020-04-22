###########################################################
##
## caesarean ssi data excess length of stay analysis
##
## N Green
## 11-02-13
##
###########################################################
#
# LOSop: operation -> departure time
# days2inf : operation -> onset of infection time
# LOSpreop: admission -> operation
# detection: 1- inpatient
#
#
#

library(foreign)  # read.dta()
library(survival)
library(mvna)
library(etm)
library(MASS)

setwd("H:\\caesarean SSI LOS")

source(".\\R code\\bootstrapping_SE.R") # standard error function for LOS

## original stata data format file
#survData <- read.dta(".\\data\\C section 2009 4110 Jan 13.dta")
#write.table(survData,".\\data\\data.txt")

survData <- read.table(".\\raw data\\data.txt")

##############
# preprocess #
##############

## look at duplicate id entries
## 831 1447 1483
x <- survData[survData$id==survData[831,"id"],]
fix(x); View(x)

## number of infections detected in hospital
det <- survData$days2inf<survData$LOSop
sum(det, na.rm=TRUE)  # 23

# who are they?
x <- survData[det & !is.na(survData$days2inf),]
fix(x); View(x)
# there are 5 additional cases all with code 5-reported by midwife or outpatient clinic 

## number of inpatient recorded infections
det <- survData$detection==1
sum(det, na.rm=TRUE) # 20
# who are they?
x <- survData[det & !is.na(survData$days2inf),]
fix(x); View(x)

## 3 of these have departure at the same day as onset of SSI
## which are not picked up above


########################
# time dependent array #
########################

## detected infection in hosptial patients only
survData.inf <- survData[survData$detection==1 & !is.na(survData$days2inf),]
## everyone else
survData.mix <- survData[survData$detection!=1 | is.na(survData$detection),]

## number of detected infected in hospital patients
ninf <- nrow(survData.inf)

## for equal times for subsequent events add error
## e.g. when admission date = discharge date
epsilon <- 0.1

originalid <- 0L  				# original array single-line (wide) format patient row id

## new multiple-line (long) time-dependent format array
## split an infected patients record across lines
## separate before and after infection
tvdata.inf <- data.frame(id=rep(1:ninf,each=2),
                         tstart=NA,tstop=NA,inf=NA, 
                         survData.inf[rep(1:ninf,each=2),
                                      c("age", "antimicpx", "asascore", "bmicat", "bloodloss", "complicatedcsd",
                                        "csdtype", "diab", "LOSpreop", "LOSar", "durationoperation",
                                        "durationofactivelabour", "woundclass")])

## number of new rows for infected patients
numRows.new <- 2*ninf

oddRows <- seq(1,numRows.new,by=2)
tvdata.inf[oddRows,"tstart"] <- 0    # admission start time always zero
tvdata.inf[,"inf"] <- rep(0:1,length.out=numRows.new)  	# not infected yet

## loop through each patient because need to adjust for same day events
if (numRows.new > 0){ # only for infected cases
  for (i in oddRows){		# fill two rows at a time
     
    originalid <- originalid+1
    
    # 1) admission -> infection
	  tvdata.inf[i,"tstop"] <- ifelse(survData.inf$days2inf[originalid]==0, epsilon, survData.inf$days2inf[originalid])	# so start<stop due to same day events

        
    # 2) infection -> death/discharge
    tvdata.inf[i+1,"tstart"] <- tvdata.inf[i,"tstop"]     # start next state at the end time of first state
    
    # stop time is the event time death, discharge or administrative censoring
    # to ensure start<stop add epsilon due to same day events
    tvdata.inf[i+1,"tstop"]  <- ifelse(survData.inf$LOSop[originalid]==survData.inf$days2inf[originalid],
                                       survData.inf$LOSop[originalid]+epsilon,
                                       survData.inf$LOSop[originalid])
  }
}



## construct an equivalent matrix for the uninfected patients
## include arbitrary unique IDs (otherwise bootstrapper ignores patient)
## split event into death and discharge and can still recover censored data

tvdata.mix <- data.frame(id=originalid+10+seq_along(survData.mix$ssi),
                         tstart=0,
                         tstop=ifelse(survData.mix$LOSop>0, survData.mix$LOSop, epsilon),
                         inf=0,
                         survData.mix[,c("age", "antimicpx", "asascore", "bmicat", "bloodloss", "complicatedcsd",
                                         "csdtype", "diab", "LOSpreop", "LOSar", "durationoperation",
                                         "durationofactivelabour", "woundclass")])

## combine both arrays
tvdata <- rbind(tvdata.inf,tvdata.mix)

#save(tvdata, file=".\\processed data\\tvdata.RData")
#load(".\\processed data\\tvdata.RData"))


##################
# Cox regression #
##################


## drop rows containing any NAs
tvdata.naomit <- na.omit(tvdata)

## minimum model
## only time of infection as covariate
summary(coxph(Surv(tstart, tstop, rep.int(1,length(tstop)))~cluster(id)+inf, data=tvdata))

## `saturated' model
## all identified covariates included
PHmodel <- coxph(Surv(tstart, tstop, rep.int(1,length(tstop)))~
                   cluster(id)+inf+antimicpx+asascore+bmicat+bloodloss+complicatedcsd+csdtype+diab+LOSpreop+durationoperation+durationofactivelabour+woundclass, data=tvdata.naomit)
PHmodel

## Akaike information criteria (AIC) deviances
## http://en.wikipedia.org/wiki/Akaike_information_criterion
outcx <- stepAIC(PHmodel, direction="both")
outcx <- stepAIC(PHmodel, scope = list(lower = ~cluster(id)+inf), direction="both")

## selected subset
PHmodel <- coxph(Surv(tstart, tstop, rep.int(1,length(tstop)))~cluster(id)+inf+LOSpreop+asascore+durationofactivelabour+bloodloss+csdtype, data=tvdata)
PHmodel

plot(survfit(Surv(tstart, tstop, rep.int(1,length(tstop)))~inf, data=tvdata), main="", lty=c(1,2), col=c(1,2), xlab="time (days)", ylab="S_hat")
legend("topright", legend=c("Control","inf"), lty = c(1,2), col=c(1,2), bty = "n", cex=0.7)

## specific survival curve comparisons
#PHmodelx1 <- survfit(PHmodel, list(inf=1,antimicpx=,asascore=,bmicat=,bloodloss=,complicatedcsd=,csdtype=,diab=,LOSpreop=,durationoperation=,durationofactivelabour=,woundclass=))
#PHmodelx2 <- survfit(PHmodel, list(inf=0,antimicpx=,asascore=,bmicat=,bloodloss=,complicatedcsd=,csdtype=,diab=,LOSpreop=,durationoperation=,durationofactivelabour=,woundclass=))
#plot(PHmodelx1$time, PHmodelx1$surv, type="l", ylim=c(0,1), xlab="Time to Infection", ylab="Surv Prob")
#lines(PHmodelx2$time, PHmodelx2$surv, lty=2)



#########################
# excess length of stay #
#########################

tra <- matrix(FALSE, 3, 3, dimnames = list(as.character(0:2), as.character(0:2)))  # admission-infection-death/discharge
tra[1, 2:3] <- TRUE
tra[2, 3] <- TRUE


## split dataset by length of time between admission and operation
## operation on the same day as admission or longer
## can't remember why I do this??
tvdataSplit <- split(tvdata, list(tvdata$LOSpreop>0))
names(tvdataSplit)<-c("LOSpreop0", "LOSpreopOver0")
head(tvdataSplit[["LOSpreop0"]])



## choose which (stratified) data to use
##--------------------------------------

## whole sample ##

tvdata.LOS <- tvdata

## create etm format array
numCases.inf <- sum(tvdata.LOS$inf==1)
numCases.mix <- nrow(tvdata.LOS)-(2*numCases.inf)

id <- tvdata.LOS$id
entry <- tvdata.LOS$tstart  	
exit  <- tvdata.LOS$tstop
from <- c(rep(0:1,numCases.inf), rep(0,numCases.mix))		
to <- c(rep(1:2,numCases.inf), rep(2,numCases.mix))
to2 <- to			# ifelse(tvdata$disch==1,3,to)

## multistate model array
msm.data <- data.frame(id,entry,exit,from,to,to2)

mvna.data <- mvna(msm.data, c("0","1","2"), tra, cens.name="cens")
etm.data <- etm(msm.data, c("0","1","2"), tra, "cens", s=0)
clos.data <- clos(etm.data)
clos.se <- sqrt(var(boot.clos(data=msm.data, state.names=c("0","1","2"), tra=tra, cens.name="cens", 0, nboot = 20)))

plot(mvna.data)
summary(etm.data)
xyplot(etm.data)
plot(clos.data)
clos.data   # exxcess length of stay
clos.se

x <- boot.clos(data=msm.data, state.names=c("0","1","2"), tra=tra, cens.name="cens", 0, nboot = 50)
qqnorm(x); qqline(x)


## same day delivery as admission ##

tvdata.LOS <- tvdataSplit$LOSpreop0

## create etm format array
numCases.inf <- sum(tvdata.LOS$inf==1)
numCases.mix <- nrow(tvdata.LOS)-2*numCases.inf

id <- tvdata.LOS$id
entry <- tvdata.LOS$tstart    
exit  <- tvdata.LOS$tstop
from <- c(rep(0:1,numCases.inf), rep(0,numCases.mix))		
to <- c(rep(1:2,numCases.inf), rep(2,numCases.mix))
to2 <- to			# ifelse(tvdata$disch==1,3,to)

## multistate model array
msm.data <- data.frame(id,entry,exit,from,to,to2)

mvna.data <- mvna(msm.data, c("0","1","2"), tra, cens.name="cens")
etm.data <- etm(msm.data, c("0","1","2"), tra, "cens", s=0)
clos.data <- clos(etm.data)
clos.se <- sqrt(var(boot.clos(data=msm.data, state.names=c("0","1","2"), tra=tra, cens.name="cens", 0, nboot = 20)))

plot(mvna.data)
summary(etm.data)
xyplot(etm.data)
plot(clos.data)
clos.data   # exxcess length of stay
clos.se



## longer than same day delivery after admission ##

tvdata.LOS <- tvdataSplit$LOSpreopOver0

## create etm format array
numCases.inf <- sum(tvdata.LOS$inf==1)
numCases.mix <- nrow(tvdata.LOS)-2*numCases.inf

id <- tvdata.LOS$id
entry <- tvdata.LOS$tstart    
exit  <- tvdata.LOS$tstop
from <- c(rep(0:1,numCases.inf), rep(0,numCases.mix))		
to <- c(rep(1:2,numCases.inf), rep(2,numCases.mix))
to2 <- to			# ifelse(tvdata$disch==1,3,to)

## multistate model array
msm.data <- data.frame(id,entry,exit,from,to,to2)

mvna.data <- mvna(msm.data, c("0","1","2"), tra, cens.name="cens")
etm.data <- etm(msm.data, c("0","1","2"), tra, "cens", s=0)
clos.data <- clos(etm.data)
clos.se <- sqrt(var(boot.clos(data=msm.data, state.names=c("0","1","2"), tra=tra, cens.name="cens", 0, nboot = 20)))

plot(mvna.data)
summary(etm.data)
xyplot(etm.data)
plot(clos.data)
clos.data   # exxcess length of stay
clos.se
