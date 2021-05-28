library("betaMix")
source("simulate01.R")
library("igraph")
library("lqmm") # for make.positive.definite
library("MASS")

reps <- 1:30
unlink("results.txt")

plotResults <- function(res, ylimtop=5) {
  # par(mfrow=c(1,2))
  ccc <- seq(0.001,0.999,length=1000)
  hist(res$z_j,freq=F,breaks=300, border="grey",main="",ylim=c(0,ylimtop),xlim=c(0,1),
       xlab=expression(sin^{2} ~ (theta)))
  lines(ccc,(res$p0)*dbeta(ccc,(N-1)/2,0.5),col=3, lwd=4)
  lines(ccc,(1-res$p0)*dbeta(ccc/res$nonNullMax,res$ahat,res$bhat),lwd=1,col=2)
  lines(ccc, (1-res$p0)*dbeta(ccc/res$nonNullMax,res$ahat,res$bhat)+
          (res$p0)*dbeta(ccc,res$etahat,0.5), col=4, lwd=2)
  cccsig <- ccc[which(ccc < res$ppthr)]
  rect(0,0, res$ppthr, 6, col='#FF7F5020', border = "orange")
  #  cat(res$ahat, res$bhat, res$etahat,(N-1)/2,"\n")
}

##listOfSeeds <- sample(1000000:9999999,20000)
##save(listOfSeeds, file="listOfSeeds.RData")

load("listOfSeeds.RData")
seedCt <- 0


# 1. some linear models where one variable is designated as the "response"
# and we're interested in finding which predictors are associated with it
# (variable selection in the large P small n setting)
#
for (N in c(200, 500)) {
  for (P in c(500, 1000)) {
    for (repno in reps){
      seedCt <- seedCt+1; set.seed(listOfSeeds[seedCt])
      gendat <- genData('lmSimple', P, N)
      res <- betaMix(gendat$dataMat, delta = 1/choose(P,2), ppr = 0.01, ind = T)
      if(repno == 1)
        plotResults(res)
      zz <- sin(res$angleMat)^2
      zz0 <- Matrix(zz < res$ppthr)
      diag(zz0) <- FALSE
      # Matrix::image(zz0[1:250,1:250])
      cat("Total True",sum(gendat$trueStructure[1,]),
          "TP", sum(zz0[1,] & gendat$trueStructure[1,]),
          "FP", sum(zz0[1,] & !gendat$trueStructure[1,]),"\n")
      cat(paste('lmSimple', P, N, "", repno, sum(gendat$trueStructure[1,]),
                 sum(zz0[1,] & gendat$trueStructure[1,]),
                 sum(zz0[1,] & !gendat$trueStructure[1,]),
                listOfSeeds[seedCt], sep = "\t"),"\n", 
          file="results.txt", append=TRUE)
    }
  }
}


# 2. Still a linear model with one designated response variable, but the predictors
# are correlated (simulateLinearModel1 has more options to control simulation settings)
#
for (N in c(200, 500)) {
  for (P in c(500, 1000)) {
    for (cortype in c("ar1","hub")) {
      for (corcoef in c(0.3, 0.9)) {
        for (repno in reps){
          seedCt <- seedCt+1; set.seed(listOfSeeds[seedCt])
          gendat <- genData('lmCor', P, N, cortype, corcoef)
          res <- betaMix(gendat$dataMat, delta = 1/choose(P,2), ppr = 0.01, ind = T)
          if(repno == 1)
            plotResults(res)
          zz <- sin(res$angleMat)^2
          zz0 <- Matrix(zz < res$ppthr)
          diag(zz0) <- FALSE
          # Matrix::image(zz0[1:250,1:250])
          cat("Total True",sum(gendat$trueStructure[1,]),
              "TP", sum(zz0[1,] & gendat$trueStructure[1,]),"FP", sum(zz0[1,] & !gendat$trueStructure[1,]),"\n")
          cat(paste('lmCor', P, N, paste0(cortype,corcoef), repno, sum(gendat$trueStructure[1,]),
                    sum(zz0[1,] & gendat$trueStructure[1,]),
                    sum(zz0[1,] & !gendat$trueStructure[1,]),
                    listOfSeeds[seedCt], sep = "\t"),"\n", 
              file="results.txt", append=TRUE)
        }
      }
    }
  }
}


# 3. P repeated measurements for N subjects, with an AR(1) correlation structure
#
for (N in c(200, 500)) {
  for (P in c(500, 1000)) {
    for (corcoef in c(0.3, 0.9)) {
      for (repno in reps){
        seedCt <- seedCt+1; set.seed(listOfSeeds[seedCt])
        gendat <- genData('ar1', P, N, corcoef)
        res <- betaMix(gendat$dataMat, delta = 1/choose(P,2), ppr = 0.01, ind = T)
        if(repno == 1)
          plotResults(res)
        zz <- sin(res$angleMat)^2
        zz0 <- Matrix(zz < res$ppthr)
        diag(zz0) <- FALSE
        # Matrix::image(zz0[1:250,1:250])
        cat("Total True",sum(gendat$trueStructure)/2,
            "TP", sum(zz0 & (toeplitz(1:P)==2) & gendat$trueStructure)/2,
            "FP", sum(zz0 & (toeplitz(1:P)==2) & !gendat$trueStructure)/2,"\n")
        cat(paste('ar1', P, N, corcoef, repno, sum(gendat$trueStructure)/2,
                  sum(zz0 & (toeplitz(1:P)==2) & gendat$trueStructure)/2,
                  sum(zz0 & (toeplitz(1:P)==2) & !gendat$trueStructure)/2,
                  listOfSeeds[seedCt], sep = "\t"),"\n", 
            file="results.txt", append=TRUE)
      }
    }
  }
}


# 4. Similar (P repeated measurements for N subjects) but only the first K have 
# an AR(1) correlation structure, and the rest are iid
#
for (N in c(200, 500)) {
  for (P in c(500, 1000)) {
    for (corcoef in c(0.3, 0.9)) {
      for (repno in reps){
        seedCt <- seedCt+1; set.seed(listOfSeeds[seedCt])
        gendat <- genData('blockar1', P=P, N=N, K=20, corcoef=corcoef)
        res <- betaMix(gendat$dataMat, delta = 1/choose(P,2), ppr = 0.01, ind = T)
        if(repno == 1)
          plotResults(res)
        zz <- sin(res$angleMat)^2
        zz0 <- Matrix(zz < res$ppthr)
        diag(zz0) <- FALSE
        # Matrix::image(zz0[1:250,1:250])
        cat("Total True",sum(gendat$trueStructure)/2,
            "TP", sum(zz0 & (toeplitz(1:P)==2) & gendat$trueStructure)/2,
            "FP", sum(zz0 & (toeplitz(1:P)==2) & !gendat$trueStructure)/2,"\n")
        cat(paste('blockar1', P, N, paste0("K20",corcoef), repno, sum(gendat$trueStructure)/2,
                  sum(zz0 & (toeplitz(1:P)==2) & gendat$trueStructure)/2,
                  sum(zz0 & (toeplitz(1:P)==2) & !gendat$trueStructure)/2,
                  listOfSeeds[seedCt], sep = "\t"),"\n", 
            file="results.txt", append=TRUE)
      }
    }
  }
}


# 5. A band correlation structure
#
for (N in c(200, 500)) {
  for (P in c(500, 1000)) {
    for (corcoef in c(0.3, 0.9)) {
      for (BW in c(30,150)){
        for (repno in reps){
          seedCt <- seedCt+1; set.seed(listOfSeeds[seedCt])
          gendat <- genData('band', P=P, N=N, corcoef=corcoef, BW=BW)
          res <- betaMix(gendat$dataMat, delta = 1/choose(P,2), ppr = 0.01, ind = T)
          if(repno == 1)
            plotResults(res)
          zz <- sin(res$angleMat)^2
          zz0 <- Matrix(zz < res$ppthr)
          diag(zz0) <- FALSE
          # Matrix::image(zz0[1:250,1:250])
          cat("Total True",sum(gendat$trueStructure)/2,
              "TP", sum(zz0 & gendat$trueStructure)/2,
              "FP", sum(zz0 & !gendat$trueStructure)/2,"\n")
          cat(paste('band', P, N, paste0(BW,corcoef), repno, sum(gendat$trueStructure)/2,
                    sum(zz0 & gendat$trueStructure)/2,
                    sum(zz0 & !gendat$trueStructure)/2,
                    listOfSeeds[seedCt], sep = "\t"),"\n", 
              file="results.txt", append=TRUE)
        }
      }
    }
  }
}


# 6. Cluster configurations
#
for (N in c(200, 500)) {
  for (P in c(500, 1000)) {
    for (corcoef in c(0.3, 0.9)) {
      for (clustsize in c(P/20, P/10, P/2)) {
        for (repno in reps){
          seedCt <- seedCt+1; set.seed(listOfSeeds[seedCt])
          gendat <- genData('cls', P=P, N=N, corcoef=corcoef,clustsize=clustsize)
          res <- betaMix(gendat$dataMat, delta = 1/choose(P,2), ppr = 0.01, ind = T)
          if(repno == 1)
            plotResults(res)
          zz <- sin(res$angleMat)^2
          zz0 <- Matrix(zz < res$ppthr)
          diag(zz0) <- FALSE
          # Matrix::image(zz0[1:250,1:250])
          cat("Total True",sum(gendat$trueStructure)/2,
              "TP", sum(zz0  & gendat$trueStructure)/2,
              "FP", sum(zz0 & !gendat$trueStructure)/2,"\n")
          cat(paste('cls', P, N, paste0(clustsize,corcoef), repno, sum(gendat$trueStructure)/2,
                    sum(zz0 & gendat$trueStructure)/2,
                    sum(zz0 & !gendat$trueStructure)/2,
                    listOfSeeds[seedCt], sep = "\t"),"\n", 
              file="results.txt", append=TRUE)
        }
      }
    }
  }
}


# 7. similar to 6, but the cluster sizes are not fixed
#
for (N in c(200, 500)) {
  for (P in c(500, 1000)) {
    for (corcoef in c(0.3, 0.9)) {
      for (repno in reps){
        seedCt <- seedCt+1; set.seed(listOfSeeds[seedCt])
        gendat <- genData('rbk', P=P, N=N, g=40, corcoef=corcoef)
        res <- betaMix(gendat$dataMat, delta = 1/choose(P,2), ppr = 0.01, ind = T)
        if(repno == 1)
          plotResults(res)
        zz <- sin(res$angleMat)^2
        zz0 <- Matrix(zz < res$ppthr)
        diag(zz0) <- FALSE
        # Matrix::image(zz0[1:250,1:250])
        cat("Total True",sum(gendat$trueStructure)/2,
            "TP", sum(zz0  & gendat$trueStructure)/2,
            "FP", sum(zz0 & !gendat$trueStructure)/2,"\n")
        cat(paste('rbk', P, N, paste0("g40",corcoef), repno, sum(gendat$trueStructure)/2,
                  sum(zz0 & gendat$trueStructure)/2,
                  sum(zz0 & !gendat$trueStructure)/2,
                  listOfSeeds[seedCt], sep = "\t"),"\n", 
            file="results.txt", append=TRUE)
      }
    }
  }
}


# 8. Hub configurations (unlike the cluster where each node in the cluster is
# connected to all other nodes in the cluster, in the hub, all are connected to
# one central node)
#
for (N in c(200, 500)) {
  for (P in c(500, 1000)) {
    for (corcoef in c(0.3, 0.9)) {
      for (hubsize in c(P/50, P/20)) {
        for (repno in reps){
          seedCt <- seedCt+1; set.seed(listOfSeeds[seedCt])
          gendat <- genData('hub', P=P, N=N, corcoef=corcoef,hubsize=hubsize)
          res <- betaMix(gendat$dataMat, delta = 1/choose(P,2), ppr = 0.01, ind = T)
          if(repno == 1)
            plotResults(res)
          zz <- sin(res$angleMat)^2
          zz0 <- Matrix(zz < res$ppthr)
          diag(zz0) <- FALSE
          # Matrix::image(zz0[1:250,1:250])
          cat("Total True",sum(gendat$trueStructure)/2,
              "TP", sum(zz0  & gendat$trueStructure)/2,
              "FP", sum(zz0 & !gendat$trueStructure)/2,"\n")
          cat(paste('hub', P, N, paste0(hubsize,corcoef), repno, sum(gendat$trueStructure)/2,
                    sum(zz0 & gendat$trueStructure)/2,
                    sum(zz0 & !gendat$trueStructure)/2,
                    listOfSeeds[seedCt], sep = "\t"),"\n", 
              file="results.txt", append=TRUE)
        }
      }
    }
  }
}


# 9. Cycle configurations - each cluster forms a cycle, e.g. 1--2--3--...--25--1
#
for (N in c(200, 500)) {
  for (P in c(500, 1000)) {
    for (corcoef in c(0.3, 0.9)) {
      for (cyclesize in c(P/20, P/10, P/2)) {
        for (repno in reps){
          seedCt <- seedCt+1; set.seed(listOfSeeds[seedCt])
          gendat <- genData('cycle', P=P, N=N, corcoef=corcoef,cyclesize=cyclesize)
          res <- betaMix(gendat$dataMat, delta = 1/choose(P,2), ppr = 0.01, ind = T)
          if(repno == 1)
            plotResults(res)
          zz <- sin(res$angleMat)^2
          zz0 <- Matrix(zz < res$ppthr)
          diag(zz0) <- FALSE
          # Matrix::image(zz0[1:250,1:250])
          cat("Total True",sum(gendat$trueStructure)/2,
              "TP", sum(zz0  & gendat$trueStructure)/2,
              "FP", sum(zz0 & !gendat$trueStructure)/2,"\n")
          cat(paste('cycle', P, N, paste0(cyclesize,corcoef), repno, sum(gendat$trueStructure)/2,
                    sum(zz0 & gendat$trueStructure)/2,
                    sum(zz0 & !gendat$trueStructure)/2,
                    listOfSeeds[seedCt], sep = "\t"),"\n", 
              file="results.txt", append=TRUE)
        }
      }
    }
  }
}

