##
## Matching caesarean script
##
## N.Green
## 17/04/13
##


## suppress warning message
options(warn=-1)

setwd("H:\\caesarean SSI LOS")

survData <- read.table(".\\raw data\\data.txt")

## read-in matching look-up table
## and convert to list type
matches <- as.data.frame( t(read.table(".\\processed data\\Match LOS Mar 13.txt", fill=TRUE)))
colnames(matches) <- matches[1,]
matches <- matches[-1,]
matches <- as.list(matches)
matches <- lapply(matches, function(x) x[!is.na(x)] )
matches

## dummy test data
#matches$a <- c("b","c","d")
#survData <- data.frame(serialnumber=c("a","b","c","d"), LOSar=c(1,2,3,4))


######################
## helper functions ##
######################

meandiff <- function(matches, survData){
  ## difference between each infected case and MEAN of matched samples
  ## so doesn't account for difference in within matches sample sizes
  ## LOSar: time from admission to departure
  ## returns: vector
  
  ## drop empty vectors in list
  matches <- matches[sapply(matches, function(x) length(x) > 0)]
  
  diff <- c(NULL, NULL)
  for (serialno in names(matches)){
    diff <- rbind(diff, c(serialno, survData$LOSar[survData$serialnumber==serialno]-
                                    mean(survData$LOSar[survData$serialnumber%in%matches[[serialno]]]))
    )
  }

diff
}


listdiff <- function(matches, survData){
  ## difference between each infected case and EACH of matched samples
  ## LOSar: time from admission to departure
  ## returns: list
  
  ## drop empty vectors in list
  matches <- matches[sapply(matches, function(x) length(x) > 0)]
  
  diff <- list()
  for (serialno in names(matches)){
    diff[[serialno]] <- survData$LOSar[survData$serialnumber==serialno] -
                        survData$LOSar[survData$serialnumber%in%matches[[serialno]]]
    }
  
  ## add some 'jiggle'
  diff <- lapply(diff, function(x) x+rnorm(1,0,0.01))
  
diff
}


matches.jk <- function(matches, survData){
  ## standard error
  ## by jackknifing matching table (leave-one-out)
  ## returns: numeric
  
  out <- NULL
  n <- 1
  
  ## drop empty vectors in list
  matches <- matches[sapply(matches, function(x) length(x) > 0)]
  
  for (i in 1:length(matches)){
    for (j in 1:length(matches[[i]])){
      matches.temp <- matches
      matches.temp[[i]] <- matches.temp[[i]][-j]
      ## aggregated mean
#      out[n] <- mean(as.numeric(meandiff(matches.temp, survData)[,2]), na.rm=TRUE)
      ## weighted mean
      out[n] <- mean(unlist(listdiff(matches.temp, survData)), na.rm=TRUE)
      
      n <- n + 1
    }
  }

out
}



###################
## data analysis ##
###################

z <- 1.96
z.vec <- c(-z, z)


## using original matching table
##------------------------------

matches.jk.orig <- suppressWarnings(matches.jk(matches, survData))

## each within group difference
indiv.mean.orig <- meandiff(matches, survData)

## overall expected excess LOS and standard error
overall.mean.orig <- mean(as.numeric(indiv.mean.orig[,2]))
overall.se.orig <- sqrt(var(matches.jk.orig))
overall.CI.orig <- overall.mean.orig + z.vec*overall.se.orig

overall.mean.orig
overall.se.orig

## weighted mean
wmean.orig <- mean(unlist(listdiff(matches, survData)))
wmean.orig

overall.CI.orig <- wmean.orig + z.vec*overall.se.orig
overall.CI.orig


## include matching by overall LOS prior to infection
##---------------------------------------------------
## LOSar >= LOSpreop+days2inf
## i.e. (time from admission to departure) >= (time from admission to operation) + (time from operation to infection)

## update matching look-up table
## by removing non-matches
matches.ar <- matches
for (serialno in names(matches)){
	matches.ar[[serialno]] <- matches.ar[[serialno]][survData$LOSar[survData$serialnumber%in%matches.ar[[serialno]]] >=
					 	                                      (survData$LOSpreop[survData$serialnumber==serialno]+
						                                       survData$days2inf[survData$serialnumber==serialno]) ]
}
matches.ar

matches.jk.ar <- suppressWarnings(matches.jk(matches.ar, survData))

## each within group difference
indiv.mean.ar <- meandiff(matches.ar, survData)

## overall expected excess LOS and standard error including LOSar
overall.mean.ar <- mean(as.numeric(indiv.mean.ar[,2]), na.rm=TRUE)
overall.se.ar <- sqrt(var(matches.jk.ar))
#overall.CI.ar <- overall.mean.ar + z.vec*overall.se.ar

overall.mean.ar
overall.se.ar

## weighted mean
wmean.ar <- mean(unlist(listdiff(matches.ar, survData)))
wmean.ar

overall.CI.ar <- wmean.ar + z.vec*overall.se.ar
overall.CI.ar

## include matching by postop length of stay prior to infection
##-------------------------------------------------------------
## LOSop >= days2inf
## i.e. (time from operation to discharge) >= (time from operation to infection)
## this is a slightly stronger condition than above

## update matching look-up table
matches.op <- matches
for (serialno in names(matches)){
  matches.op[[serialno]] <- matches.op[[serialno]][survData$LOSop[survData$serialnumber%in%matches.op[[serialno]]] >=
                                                     survData$days2inf[survData$serialnumber==serialno] ]
}
matches.op

matches.jk.op <- suppressWarnings(matches.jk(matches.op, survData))

## each within group difference
indiv.mean.op <- meandiff(matches.op, survData)

## overall expected excess LOS and standard error including LOSop
overall.mean.op <- mean(as.numeric(indiv.mean.op[,2]), na.rm=TRUE)
overall.se.op <- sqrt(var(matches.jk.op))
#overall.CI.op <- overall.mean.op + z.vec*overall.se.op

overall.mean.op
overall.se.op


## weighted mean
wmean.op <- mean(unlist(listdiff(matches.op, survData)))
wmean.op
overall.CI.op <- wmean.op + z.vec*overall.se.op
overall.CI.op

###########
## plots ##
###########

hist(matches.jk.orig, 20)
lines(density(matches.jk.orig, bw="SJ", adjust=4))
qqnorm(matches.jk.orig); qqline(matches.jk.orig)

hist(matches.jk.ar, 20)
lines(density(matches.jk.ar, bw="SJ", adjust=4))
qqnorm(matches.jk.ar); qqline(matches.jk.ar)

hist(matches.jk.op, 20)
lines(density(matches.jk.op, bw="SJ", adjust=4))
qqnorm(matches.jk.op); qqline(matches.jk.op)

#qqplot(out, rt(1000,df=1))
#qqline(out, rt(1000,df=1))
#hist(rt(173,df=1), 200)

