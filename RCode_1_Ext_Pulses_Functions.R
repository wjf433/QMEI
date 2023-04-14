# -------------------------------------------------------------------------- #
#
# These R code functions were modified from the paper
# Wang, S., & Zhong, L. (2018). Estimating the number of pulses in a mass extinction. 
# Paleobiology, 44(2), 199-218. doi:10.1017/pab.2016.30, 
#
# Parallelisation was done by Prof. Steve Wang and only work on macOS or Unix operating systems
# -------------------------------------------------------------------------- #
#
# Version 22/03/2023
#
# -------------------------------------------------------------------------- #
# This code requires the functions within the knnflex and parallel packages

#If you need to install the "knnflex" package use the following two lines of code.
install.packages("remotes")
remotes::install_github("Dasonk/knnflex")

# The "parallel" package should be pre-installed with R.
library("knnflex")
library("parallel")

# ------------------------- FUNCTION DEFINITIONS --------------------------- #

knn.probability <- function(train, test, y, dist.matrix, k=1, ties.meth="min") {

# number of predictions to make
n<-length(test)

# sort the indexes for the training and test sets
if (is.unsorted(train)) train<-sort(train)
if (is.unsorted(test)) test<-sort(test)

# only need the rows for the test data and columns
#for the training data
d<-dist.matrix[test,train]

# ensure y is a factor
y<-as.factor(y)

# only need the responses for the training data
if (length(y)>length(train)) y<-y[train]

# calculate closest neighbors and
# return aggregate response for the k closest neighbors
if (n==1) {
  d<-rank(d, ties.method = ties.meth)
  x<-classprob(y[d <= k])
  x<-data.frame(x)
  names(x)<-test
  row.names(x)<-levels(y)
  return(x)
  }
else {
  d<-t(apply(d,1,function(x) rank(x,ties.method=ties.meth)))
  x<-apply(d,1,function(x) classprob(y[x<=k]))
  row.names(x)<-levels(y)
  return(x)
  }
}


knn.dist <- function(x, dist.meth="euclidean", p=2) {
# create a distance matrix using all values in the data
d<-as.matrix(dist(x,dist.meth,p))
# fix for some small high persision errors
round(d,digits=15)
}


classprob <- function(x){
  x <- as.factor(x)
  n <- nlevels(x)
  votes <- rep(0, n)
  for (i in 1:length(x)) votes[as.integer(x[i])] <- votes[as.integer(x[i])]+1
  votes/length(x) }


rangechart <- function(data, theta=NA, orientup=1, addbase=0, ylab="")  {
  # draw range chart, with true extinctions given by theta (use NA if not wanted)
  # orientup = 1 if larger values are at the top, 0 if at the bottom
  # addbase = amount to add to y-axis labels
  ord <- order( data[1,])
  data <- data[,ord]           # sort data by highest occurrence
  theta <- theta[ord]
  numtaxa <- dim(data)[2]
  maxocc <- dim(data)[1]
  ymax <- max(data[1,])
  if(orientup)     ylims <- c(-.3, ymax*1.1) + addbase
  if(!orientup)    ylims <- c(1.1*ymax, -.3) + addbase
  plot( NULL, pch="", xlab="", ylab=ylab, xlim=c(0,numtaxa+1), ylim=ylims, 
        bty="L", xaxt="n", cex.axis=.8)
  for(taxon in 1:numtaxa)  {
    # draw taxon lines
    lines( c(taxon,taxon), c(0,data[1,taxon])+addbase, lwd=.9, col="darkgray" )
    # plot fossil horizons
    points(rep(taxon,maxocc), data[,taxon]+addbase, pch=16, cex=.5)
    # plot extinction horizons
    points(taxon, theta[taxon]+addbase, pch="-", cex=2.7, col="red")
  }
}


getpartitions <- function(n, m) {
# enumerate all vectors that have m elements with sum n.
# n must >= m
partitionTable <- list()
  for (i in 1:m) {
    for (j in i:(n-(m-i))) {
      tag = paste("m",i,"n",j)
      if (!is.null(partitionTable[[tag]]))
        next      
      if (i==1) {
        vec <- array(j, c(1,1))
      } else {
        vec <- array(0, c(0,i))
        for (k in (i-1):(j-1)) {
          tag = paste("m",i-1,"n",k);
          subvec <- partitionTable[[tag]]
          subvec <- cbind(c(j-k), subvec)
          vec <- rbind(vec, subvec)
        }
      }
      tag = paste("m",i,"n",j)
      partitionTable[[tag]] <- vec      
    }
  }
  tag = paste("m",m,"n",n)     
  return(partitionTable[[tag]])
}


simnumocc <- function(numtaxa, mu, minocc=3, maxocc)
# Simulate number of fossil occurrences per taxon
# numtaxa: number of taxa
# mu: Poisson parameter; controls mean number of occurrences per taxon
# maxocc: max number of occurrences allowed per taxa
{
  # Generate number of occurrences for each taxa;
  #   check to make sure none have <minocc or >maxocc horizons;
  #   if so, adjust them to minocc or maxocc, respectively
  numocc <- rpois(numtaxa, mu)
  numocc[numocc<minocc] <- minocc
  numocc[numocc>maxocc] <- maxocc
  return(numocc)
}


simdata <- function(numtaxa=10, numocc, thetavec, lambdavec)
# This is a lower-level function called by createdataset
# Returns a matrix of fossil occurrences according to given arguments
# numtaxa: number of taxa
# numocc: vector of number of occurrences per taxon
# thetavec: vector of true extinction horizons/times
# lambdavec: vector of lambdas (recovery potential parameter)
{
  rrefbeta <- function(n,lambda)  {
    if(lambda<=0)  { 
      return(rbeta(n, 1, 1-lambda))
    }  else  return(rbeta(n, 1+lambda, 1))
  }

  # Simulate the occurrences for each taxa
  maxocc <- max(numocc)
  data <- matrix(NA, maxocc, numtaxa)
  for(i in 1:numtaxa)
    data[,i] <- c( sort( rrefbeta(numocc[i],lambdavec[i])*thetavec[i], decr=T),
                   rep(NA, maxocc-numocc[i]) )
  return(round(data,1))
}


createdataset <- function(numtaxa, samplesizes=NULL, mu, numpulses=1, maxtheta=100, 
                          maxocc, lambdamean=0, lambdaSD=0, OUTPUT=0)  {
# Calls the lower-level function simdata
# Must specify either samplesizes or mu to control sample size (#occurrences/taxon)
# maxocc = maximum possible number of occurrences (can have fewer occurrences)
# mu = average number of occurrences per taxon (Poisson mean)
# lambdamean = mean of parameter for reflected Beta recovery potential (scalar)
# lambdaSD = SD of parameter for reflected Beta recovery potential (scalar)
# warning: does not check if argument values are consistent


    # read in list of partitions from global variables
    listOfPartitions <- partitionlist[[numpulses]]
    numpartitions <- dim(listOfPartitions)[1] 

    # set pulse locations depending on number of pulses
    # values are random but constrained to fall with some spacing from each other
    if(numpulses==1)  pulselocations <- maxtheta
    if(numpulses==2)  pulselocations <- c(runif(1,.2,.8), 1) * maxtheta
    if(numpulses>=3)  
      pulselocations <- c( seq(0,1,by=1/numpulses)[-c(1,numpulses+1)] + 
                           runif(numpulses-1,-.25/numpulses,.25/numpulses), 1 ) * maxtheta

    # choose which taxa go extinct at which pulses
    # note: the difference between pulselocations and thetas is that 
    #       the former has numpulses elements, whereas the latter has numtaxa elements 
    currrow <- sample(1:numpartitions, 1)
    currpartition <- listOfPartitions[currrow,]         # choose a random partition
    thetas <- rep(pulselocations, times=currpartition)   # fill in the thetas by taxon

    # assign sample sizes if provided
    if(!is.null(samplesizes))  n <- samplesizes
    # else get random sample sizes
    if(is.null(samplesizes))   
      n <- simnumocc(numtaxa=numtaxa, mu=mu, maxocc=maxocc) 

    # choose random lambdas 
    lambdas <- rnorm(numtaxa, lambdamean, lambdaSD)

    if(OUTPUT)  {  
      cat("\nCreate dataset\n")
      cat("  numpulses\t", numpulses,"\n")
      cat("  theta \t", round(thetas,1),"\n")
      cat("  n     \t", n,"\n")
    }
    data <- simdata(numtaxa=numtaxa, numocc=n, thetavec=thetas, lambdavec=lambdas)
    data <- data[1:max(n),]                  # omit excess rows of NAs
    return(data)
}


calcIC <- function(data, maxnumpulses=4, WRITE=0)  {
# calculates the max likelihood, AIC, & BIC over all scenarios for each number of pulses 

  # prepare dataset
  numtaxa <- dim(data)[2]                         # number of taxa 
  n <- rep(NA,numtaxa)                            # vector of sample sizes for each taxon
  ord <- order(data[1,])
  data <- data[,ord]                              # sort taxa by highest occurrence
  y <- data[1,]                                   # vector of highest occurrences
  n <- apply( 1-apply(data,2,is.na), 2, sum)      # get sample sizes
  theta <- rep(y[numtaxa],numtaxa)                # initialize all pulses to highest occ. 
  best <- matrix(NA,maxnumpulses,numtaxa+5)       # rows = pulses, columns = best fit info
  if(WRITE)  filename <- "outIC.txt"              # name of output file
  colnames(best) <- c(paste("th",1:numtaxa, sep=""),"loglik","aic","aicw","bic","bicw")


  # run 1-pulse case separately
  pulses <- 1
  comb <- 1
  loglik <- sum( log(n*(y/theta)^(n-1)) )
  # save results
  best[pulses,1:numtaxa] <- theta
  best[pulses, "loglik"] <- loglik
  if(WRITE)  write(c(pulses, comb, theta, loglik), file=filename,sep="\t", ncol=1+numtaxa)


  # loop over multiple-pulse cases
  for(pulses in 2:maxnumpulses)  {
    currmaxlik <- -Inf                                    # reset max likelihood
    currpartitions <- partitionlist[[pulses]]             # read in from global variable
    numpartitions <- dim(currpartitions)[1]

    # loop over all possible partitions for this number of pulses
    for(i in 1:numpartitions)  {
      # select possible pulse locations, which must coincide w/ highest occ. of some taxon
      pulselocs <- y[cumsum(currpartitions[i,])]
      # assign taxa to most likely pulse locations
      theta <- rep(pulselocs, times=currpartitions[i,])

      # calculate likelihood and check if this scenario is better
      loglik <- sum( log(n*(y/theta)^(n-1)) )
      if(loglik > currmaxlik)  {             # if this is the new max, then save results
        currmaxlik <- loglik
        best[pulses,1:numtaxa] <- theta
        best[pulses,"loglik"] <- loglik
      }

      # write results and reset pulse locations for next iteration
      if(WRITE)  write(c(pulses,comb,theta,loglik),sep="\t",file=filename,ncol=1+numtaxa,append=T)
      theta <- rep(y[numtaxa],numtaxa)
    }
  }

  # calculate AIC and BIC
  for(pulses in 1:maxnumpulses)  {
    loglik <- best[pulses, "loglik"]
    best[pulses, "aic"] <- -2*loglik + pulses*2 + (2*pulses*(pulses+1))/(sum(n)-pulses-1) 
    best[pulses, "bic"] <- -2*loglik + pulses*log(sum(n))
  }
  # calculate AIC weights
  temp <- best[,"aic"] - min(best[,"aic"])
  temp <- exp(-1/2*temp)
  best[,"aicw"] <- temp / sum(temp)
  # calculate BIC weights
  temp <- best[,"bic"] - min(best[,"bic"])
  temp <- exp(-1/2*temp)
  best[,"bicw"] <- temp / sum(temp)  

  rm(temp);    return(best)
}


howmanypulses <- function(data, numsims=300, maxnumpulses=4, kval=20)  {

  # get sample sizes and number of taxa
  n <- apply( 1-apply(data,2,is.na), 2, sum)
  numtaxa <- length(n)


  # simulate training data based on the given sample sizes

  # initialize
  theta <- rep(NA, numtaxa)                         # vector of extinction positions
  maxtheta <- 100
  pulselocations <- rep(NA,numtaxa)    
  weightsAIC <- weightsBIC <- matrix(NA, maxnumpulses*numsims,maxnumpulses)
  index <- 1
  train <- vector("list", numsims)                  # list to store training datasets

  # generate datasets with the same parameters as the passed-in dataset
  cat("\nSimulating training datasets")
  for(numpulses in 1:maxnumpulses)  {               # loop over number of pulses
    cat("\n  pulses =", numpulses, ": ")

    # create training datasets
    for(rep in 1:numsims) 
      train[[rep]] <- createdataset(numtaxa=numtaxa, samplesizes=n, numpulses=numpulses, 
                                    maxtheta=maxtheta, maxocc=max(n)) 
                                    
    # calculate likelihood, AIC, and BIC (parallelized)
    bestlist <- mclapply(train, calcIC, maxnumpulses=maxnumpulses, mc.cores=detectCores())
    
    # store AIC and BIC weights
    for(rep in 1:numsims)  {
      best <- bestlist[[rep]]
      weights <- t( best[,c("aicw","bicw")] )
      weightsAIC[index,] <- weights["aicw",]
      weightsBIC[index,] <- weights["bicw",]
      index <- index + 1
    } 

  }  # end loop over pulses

  # results
  truth <- rep(1:maxnumpulses, each=numsims)
  trainset <- cbind(log(weightsAIC), log(weightsBIC))

  
  # get vector of AIC/BIC for the test point 
  best <- calcIC(data, maxnumpulses=maxnumpulses)
  weights <- t( best[,c("aicw","bicw")] )
  testaic <- t(as.matrix(log(weights["aicw",]),1,numtaxa))
  testbic <- t(as.matrix(log(weights["bicw",]),1,numtaxa))
  testpoint <- cbind(testaic,testbic)  

  # find nearest neighbors of testpoint in the training set
  x <- rbind(trainset, testpoint)
  kdist <- knn.dist(x)
  probs <- knn.probability(1:(maxnumpulses*numsims), maxnumpulses*numsims+1, 
                           y=as.factor(truth), kdist, k=kval)
  # return results
  cat("\n\n")
  temp <- cbind(best,probs)
  colnames(temp)[numtaxa+6] <- "conf"
  return(temp)
}


plotpulses <- function(data, pulseout, maxplots=dim(pulseout)[1], orientup=1, addbase=0, units="")  {
  # pulseout is the output from howmanypulses, data is the raw dataset
  # maxnumplots = how many columns of plots to show
  # orientup = 1 if high values should be at the top of the plot
  numtaxa <- dim(data)[2]
  ord <- order( data[1,])
  data <- data[,ord]           # sort data by highest occurrence

  # Range charts for each number of pulses
  par(omi=c(2,1.5,1,0)*.3, mai=c(0,2,2.3,1.5)*.2)
  layout(matrix(1:(maxplots*2), 2,maxplots, byrow=T))
  for(pulsenum in 1:maxplots)  {
    rangechart(data, pulseout[pulsenum,1:numtaxa], orientup=orientup, addbase=addbase)
    title(main=ifelse(pulsenum==1, "1 pulse", paste(pulsenum,"pulses")))
    if(pulsenum==1)  mtext(units, side=2, line=3, cex=.7)
  }

  # Confidence levels barplot
  for(pulsenum in 1:maxplots)  {
    barplot(as.matrix(pulseout[pulsenum,"conf"]), axes=F, ylim=c(0,1), xlim=c(0,1), 
            border=NA, width=.3, space=1.2, names.arg="", col="darkgray")
    axis(side=1, lwd.ticks=0, labels=F)
    if(pulsenum==1)  mtext("confidence", side=2, line=3, cex=.7)
    mtext(as.matrix(pulseout[pulsenum,"conf"]), side=1, line=.5, cex=.6)
  }
}


pulsesCI <- function(probs, conf=.9)  {
# takes a vector of probabilities/confidences and returns a CI for the number of pulses
# vector of probs assumed to sum to 1 (not checked)
# the CI is discrete and may include non-consecutive values
# if there are ties, the smaller number of pulses will be chosen
# confidence level may be > conf, due to discreteness
# returns CI for number of pulses and the achieved confidence level
  ord <- order(probs, decreasing=T)
  probsord <- probs[ord]
  cumprobs <- cumsum(probsord)
  i <- 1
  while( (cumprobs[i]<conf) && !isTRUE(all.equal(cumprobs[i],conf)) )
    i <- i + 1
  last <- i                              # index of last value in CI
  first <- 1                             # index of first value in CI 
  CI <- ord[first:last]
  return(CI)
}



# ----------------------- END FUNCTION DEFINITIONS -------------------------- #