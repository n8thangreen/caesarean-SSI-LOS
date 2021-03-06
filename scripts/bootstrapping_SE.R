##################################
#
# Compute bootstrapped SE
# for multistate model
# excess length of stay (LOS)
# package:etm example
#
##################################

## TODO ##
## function does what it should but gives incorrect answers
## or at least different to e.phi
boot.clos.alt <- function(msm.data, state.names, tra, cens.name, s = 0, nboot=10){
  #
  #
  #
  res <- double(nboot)
  etm.data <- etm::etm(data=msm.data, state.names, tra, cens.name, s, cova = FALSE)
  
  for (i in seq_len(nboot)){
    windex <- sample(x=clos(etm.data)$w.time, replace=TRUE, prob=clos(etm.data)$weights)
    ntimes <-  sapply(windex, function(x) which(clos(etm.data)$time==x))
    res[i] <- mean(clos(etm.data)$phi.case[ntimes]-clos(etm.data)$phi.control[ntimes])
  }
}
## END FUNCTION ##


boot.clos <- function(data, state.names, tra, cens.name, s = 0, nboot) {
  #
  # generate a bootstrap sample of expected LOS
  # then used to compute bootstrapped SE
  # data: msm data in etm format
  # nboot: number of bootstrap samples. Other arguments are as in etm()
  #
  # call: se <- sqrt(var(boot.clos(msm.data, c("0","1","2"), tra, NULL, 0, nboot = 10)))
  
    res <- double(nboot)
    for (i in seq_len(nboot)) {
      ####################
      # bootstrap sample #
      ####################
      ## the same size as original
      #index <- sample(unique(data$id), replace = TRUE)  
      
      ## or a defined fraction of original sample to speed-up computation (stratified by inf/non-inf)
      #index <- sample(unique(data$id), replace = TRUE, size=length(unique(data$id))/1000)  
      index.inf <- sample(x=unique(data$id[data$to==1]), replace = TRUE)  #, size=length(unique(data$id[data$to==1]))/100)
      index.mix <- sample(x=unique(data$id[data$from==0&data$to==2]), replace = TRUE) #, size=length(unique(data$id[data$from==0&data$to==2]))/10) 
      index <- c(index.inf, index.mix)
      
      ## rows corresponding to each bootstrapped patient
      ## and new ids
      linds <- sapply(index, function(x) which(x==data$id))
      indrep <- sapply(linds, function(y) length(y))

      inds <- unlist(linds)
      new.id <- rep(seq_along(index), indrep)
     

      ## create combined array of new ids and associated patient data
      dboot <- cbind(data[inds, ], new.id)
      ## replace ids
      dboot$id <- dboot$new.id
      ## remove old ids
      dboot$new.id <- NULL
      
      tr.prob <- etm::etm(dboot, state.names, tra, cens.name, s, cova = FALSE)
      res[i] <- etm::clos(tr.prob)$e.phi
    }
res
}

## problem resampling becasue the non-HCAI BSI sample is so much bigger than the HCAI BSI sample
## so that probability of sample an infected patient is very small
## for some of the organims groups where the size is in single figures calculating the se is not sensible
## sample a subset of the whole sample also reduces the proportion of infective further
## sampling from each group separately would have the same problem unless the sampling probabilities are different
## in this case, an adjustment would be needed


boot.clos2 <- function(data, state.names, tra, cens.name, s = 0, nboot) {
  # data.table version of boot_clos()
  # faster runtime
  # data: class(data.table)
  #
  res <- double(nboot)
  
  ## set id as the searchable key column
  setkey(data, "id")
  
  for (i in seq_len(nboot)) {
    
    index <- sample(unique(data[, id]), replace=TRUE)
    
    ## join data.tables for each patient together
    dboot <- data[J(index)]
    tr.prob <- etm(dboot, state.names, tra, cens.name, s, cova = FALSE)
    res[i] <- etm::clos(tr.prob)$e.phi
  }
  res
}
