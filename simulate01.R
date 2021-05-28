# Simulation code for the beta-mix model.
# In all cases, we have to provide N and P, plus the specific parameters
# for each scenario. The simulation scenario is determined by the simtype
# variable.

genData <- function(mod, P, N, ...) {
  cat("Generating data...\n")
  simdat <- switch (mod,
                 lmSimple = simulateLinearModelSimple(P, N),
                 lmCor = simulateLinearModel1(P, N,...),
                 ar1 = simulateAR1(P, N,corcoef), 
                 blockar1 = simulateBlockAR1(P, N, ...),
                 band = simulateBand(P, N, ...),
                 cls = simulateCluster(P, N, ...),
                 rbk = simulateRandomBlocks(P, N, ...),
                 hub = simulateHub(P, N, ...),
                 cycle = simulateCycle(P, N, ...)
  )
  cat("Done.\n")
  simdat
}



# Linear model with K-1 predictors which are correlated with X1, 
# plus some uncorrelated ones
simulateLinearModel1 <- function(P, N, cortype="ar1", corcoef=0.5, sig=0.1,
                                 coefs=c(1.3,6,4,2.7), cols=c(1,30,100),K=15) {
  X <- matrix(runif(N*P),nrow=N, ncol=P)
  sigma <- sqrt((1-corcoef^2)/(12*corcoef^2)) # for standard uniform X1
  TM <- Matrix(F, P+1, P+1)
  for (jj in 2:K) {
    if (cortype == "ar1") {
      X[,jj] <- X[,jj-1] + rnorm(N,0,sigma)
    }
    if (cortype == "hub") {
      X[,jj] <- X[,1] + rnorm(N,0,sigma)
      #      TM[2,3:K] <- TM[3:K,2] <- TRUE
    }
    TM[1:(K+1),1:(K+1)] <- TRUE
  }
  TM[1,cols+1] <- TM[cols+1,1] <- TRUE
  diag(TM) <- FALSE
  Y <- coefs[1] + X[,cols]%*%matrix(coefs[-1],ncol=1) + rnorm(N,0,sig)
  list(dataMat=cbind(Y,X), trueStructure=TM)
}


simulateLinearModelSimple <- function(P, N, sig=0.1, coefs=c(1.3,6,4,3), 
                                      cols=c(1,30,100)) {
  X <- matrix(runif(N*P),nrow=N, ncol=P)
  Y <- coefs[1] + X[,cols]%*%matrix(coefs[-1],ncol=1) + rnorm(N,0,sig)
  TM <- Matrix(F, P+1, P+1)
  TM[1,cols+1] <- TM[cols+1,1] <- TRUE
  diag(TM) <- FALSE
  list(dataMat=cbind(Y,X), trueStructure=TM)
}


simulateHub <- function(P, N, hubsize=20, corcoef=0.9) {
  hubs <- P/hubsize
  Sigma <- diag(1,hubsize)
  Sigma[1,] <- Sigma[,1] <- corcoef
  Sigma <- kronecker(diag(1,hubs),Sigma)
  TM <- Matrix(Sigma != 0)
  diag(TM) <- FALSE
  Sigma <- make.positive.definite(Sigma)
  list(dataMat=mvrnorm(N, rep(0,P), Sigma), trueStructure=TM)
}


simulateAR1 <- function(P, N, corcoef=0.9) {
  Sigma <- toeplitz(1:P)
  TM <- Matrix(Sigma == 2)
  for (i in 2:P) {
    Sigma[which(Sigma == i)] <- corcoef^(i-1)
  }
  list(dataMat=mvrnorm(N, rep(0,P), Sigma), trueStructure=TM)
}


simulateBlockAR1 <- function(P, N, K, corcoef=0.9) {
  MTM0 = simulateAR1(P-K,N,0)
  MTM = simulateAR1(K,N, corcoef)
  TM <- Matrix(0,P,P)
  TM[1:K,1:K] <- MTM$trueStructure
  list(dataMat=cbind(MTM$dataMat,MTM0$dataMat), trueStructure=TM)
}


simulateBand <- function(P, N, corcoef=0.3, BW=30) {
  Sigma <- toeplitz(1:P)
  for (i in 2:BW) {
    Sigma[which(Sigma == i)] <- corcoef
  }
  for (i in (BW+1):P) {
    Sigma[which(Sigma == i)] <- 0
  }
  TM <- Sigma != 0
  diag(TM) <- FALSE
  Sigma <- make.positive.definite(Sigma)
  list(dataMat=mvrnorm(N, rep(0,P), Sigma), trueStructure=TM)
}


simulateCluster <- function(P, N, clustsize=20, corcoef=0.3) {
  clusts <- P/clustsize
  Sigma <- matrix(corcoef,clustsize,clustsize)
  diag(Sigma) <- 1
  Sigma <- kronecker(diag(1,clusts),Sigma)
  TM <- Matrix(Sigma != 0)
  diag(TM) <- FALSE
  Sigma <- make.positive.definite(Sigma)
  list(dataMat=mvrnorm(N, rep(0,P), Sigma), trueStructure=TM)
}


simulateCycle <- function(P, N, cyclesize=20, corcoef=0.3) {
  cycles <- P/cyclesize
  Sigma <- diag(1, cyclesize)
  for (i in 1:(cyclesize-1))
    Sigma[i,i+1] <- Sigma[i+1,i] <- corcoef
  Sigma[cyclesize,1] <- Sigma[1,cyclesize] <- corcoef
  Sigma <- kronecker(diag(1,cycles),Sigma)
  TM <- Matrix(Sigma != 0)
  diag(TM) <- FALSE
  Sigma <- make.positive.definite(Sigma)
  list(dataMat=mvrnorm(N, rep(0,P), Sigma), trueStructure=TM)
}


simulateRandomBlocks <-function(P,N, g=20,  corcoef=0.3) {
  Sigma <- matrix(0,P,P)
  brks <- c(0,sort(sample(2:(P-1),g-1, replace = FALSE)),P)
  for (i in 1:g) {
    sz <- brks[i+1] - brks[i] #+ 1
    rng <- (brks[i]+1):brks[i+1]
    Sigma[rng,rng] <- matrix(corcoef,sz,sz)
  }
  diag(Sigma) <- 1
  TM <- Matrix(Sigma != 0)
  diag(TM) <- FALSE
  Sigma <- make.positive.definite(Sigma)
  list(dataMat=mvrnorm(N, rep(0,P), Sigma), trueStructure=TM)
}

