rm(list=ls())
library('betaMix')
library("sp")
library(Matrix)
# colors https://www.canva.com/learn/100-color-combinations/
cols <- rep("#DA70D633",26)
cols[c(2,4,15,23,24)] <- "#1995AD33"
#cols[c(25)] <- "#C99E1033"
cols[c(6,10,11,17,20)] <- "#9B4F0F33" 
cols[c(7,12,14,16,26)] <- "#3F681C33"
cols[c(3,5,9,12,13,25)] <- "#FFBB0033"


# cols <- rep("#DA70D633",21)
# cols[c(8)] <- "#C99E1033"
# cols[-c(8)] <- "#0095AD33"
# 
# cols[c(2,3,9,21)] <- "#C99E1033"
# cols[c(8,19)] <- "#9B4F0F33"
# cols[c(5,15)] <- "#3F681C33"
# cols[-c(1,2,3,5,8,9,10,11,14,15,16,19,21)] <- "#0095AD33"
# 

reduceMat <- function() {
  CountsReduced <- Matrix(0,nrow=nrow(routes), ncol=nrow(birds))
  for (i in 1:max(gc$clustNo)) {
    cls <- which(gc$clustNo == i)
    ctr <- which(gc$clustNo == i & gc$iscenter == 1)
    if(length(cls) > 1) {
      CountsReduced[ctr,] <- colSums(Counts[cls,])
    } else {
      CountsReduced[ctr,] <- Counts[cls,]
    }
  }
  
  exc <- which(rowSums(CountsReduced) == 0)
  routesReduced <- routes[-exc,]
  CountsReduced <- CountsReduced[-exc,]
  excBirds <- which(colSums(CountsReduced) == 0)
  CountsReduced <- CountsReduced[,-excBirds]
  logCounts <- log10(1+CountsReduced)
  logCounts <- logCounts+ matrix(rnorm(prod(dim(logCounts)),0,0.1),ncol=ncol(logCounts),
                                 nrow=nrow(logCounts))
  rownames(logCounts) <- routes$uid[-exc]
  colnames(logCounts) <- birds[-excBirds,3]
  logCounts
}

polyarea <- function(x,y,slices=60, eps=0.01,
                     colpoly=rgb(red=1, green = 0, blue = 0, alpha = 0.2)) {
  ctrs <- c(mean(x),mean(y))
  lengths <- (x-ctrs[1])^2+(y-ctrs[2])^2
  angles <- atan2((x-ctrs[1]),(y-ctrs[2])) + pi
  anglesord <- order(angles)
  angles <- angles[anglesord]
  lengths <- lengths[anglesord]
  x <- x[anglesord]
  y <- y[anglesord]
  #plot(x,y)
  #points(ctrs[1],ctrs[2],col=2)
  # in every slice, keep just the farthest point
  keepPt <- c()
  for (i in 1:slices) {
    inSlice <- intersect(which(angles >= (i-1)*2*pi/slices),which(angles < i*2*pi/slices))
    if (length(inSlice) > 0) {
      keepPt <- c(keepPt, inSlice[which.max(lengths[inSlice])])
    }
  }
  #points(x[keepPt],y[keepPt],col=3)
  x <- x[keepPt]
  x <- c(x,x[1])
  y <- y[keepPt]
  y <- c(y,y[1])
  angles <- angles[keepPt]
  lengths <- lengths[keepPt]
  pts <- 1:3
  #points(x[1],y[1],cex=2,col=2)
  keepPt <- c()
  for (i in 1:(length(x)-2)) {
    A <- abs(det(matrix(c(x[pts],y[pts],rep(1,3)),ncol=3,nrow=3)))
    B <- abs(det(matrix(c(x[pts[1]],ctrs[1],x[pts[3]],y[pts[1]],ctrs[2],y[pts[3]],
                          rep(1,3)),ncol=3,nrow=3)))
    #cat(i,length(x),x[pts],y[pts],A,B,"\n")
    if (A/B < eps) {
      pts <- c(pts[1],pts[3], pts[3]+1)
    } else {
      pts <- (pts[1]+1):(pts[1]+3)
      keepPt <- c(keepPt, pts[1])
    }
  }
  polygon(x[keepPt],y[keepPt],col=colpoly,border=substr(colpoly,1,7), lwd=1)
}


year <- 2015
load("birds2015.RData")
M <- as.matrix(t(logCounts))
N <- nrow(logCounts); P <- ncol(logCounts)
res <- betaMix(M, delta = 1e-5, ppr=0.001,subsamplesize=30000)
plotFittedBetaMix(M, res)
A = Matrix(sin(res$angleMat)^2 < res$ppthr)
diag(A) = FALSE

gcBirds <- graphComponents(A,minCtr = 3, type = -1)
gcBirds$labels <- colnames(A)
summarizeClusters(gcBirds)[,1:9]

CL <- collapsedGraph(A,gcBirds)

routesincls <- gcBirds$labels[which(gcBirds$clustNo == 0)]
routesincls <- rownames(gc)[which(routes$uid  %in% routesincls)]
incls <- which(routes$uid %in% routesincls)

#pdf("birds2015.pdf",width=6, height = 4.2)
plot(uscanada$x[inccoord],uscanada$y[inccoord],cex=0.1,col="grey85",main="",
     xlab="Latitude",ylab="Longitude", axes=F)
axis(1); axis(2)
# to show all the sampling routes:
#points(routes[,7], routes[,6],cex=0.4,main=clsnum, col="grey66", pch=19)
ctrs <- matrix(0,nrow=nrow(CL), ncol=2)
set.seed(119934)
for (clsnum in 1:max(gcBirds$clustNo)) {
  routesincls <- gcBirds$labels[which(gcBirds$clustNo == clsnum)]
  routesincls <- rownames(gc)[which(routes$uid  %in% routesincls)]
  incls <- which(routes$uid %in% routesincls)
  rnum  <- which(rownames(CL) == sprintf("CLS%d",clsnum))
  ctrs[rnum,] <- c(median(routes[incls,7]),median(routes[incls,6]))
  #points(routes[incls,7],routes[incls,6],cex=0.5,pch=19,col=substr(cols[clsnum],1,7))
  x0 <- routes[incls,7]
  y0 <- routes[incls,6]
  while(length(x0) < 2000) {
    x0 <- c(x0,x0+rnorm(length(x0),0,0.2))
    y0 <- c(y0,y0+rnorm(length(y0),0,0.2))
  }
  if (cols[clsnum] == "#FFBB0033") {
    text(ctrs[rnum,1],ctrs[rnum,2],clsnum,cex=.8,col="black",font=2)
  } else {
    text(ctrs[rnum,1],ctrs[rnum,2],clsnum,cex=.8,col=substr(cols[clsnum],1,7),font=2)
  }
  polyarea(x0,y0,colpoly =cols[clsnum], slices=40, eps = 0.00000001)
}
text(-110,65,year,cex = 1.2)
# to show points not in any cluster, uncomment:
#isolated <- which(gcBirds$clustNo==0)
#points(routes[isolated,7], routes[isolated,6],cex=0.4,main=clsnum, col="red", pch=19)
#dev.off()
