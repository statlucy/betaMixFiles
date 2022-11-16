#
# the analysis of the radar data with betaMix

library(betaMix)
library(igraph)

# read the data, and remove variables with zero variance
# the last column is the true classification
dat <- read.csv("ionosphere.data",header=F) # 351 by 35
grp <- dat$V35
print(table(grp))

M <- matrix(unlist(dat[,-35]),nrow=351,ncol=34)
exc <- which(apply(M,2,sd) == 0) # column2 is all zeros
exc = c(1,2)
dat <- dat[,-exc]
M <- matrix(unlist(dat[,-33]),nrow=351,ncol=32)
# we're going to treat the observation as "variables" and
# classify them by their similarity to other "variables"
# based on their vector of 33 measurements
M <- t(M)

# create a training dataset and a test data set
set.seed(228831)
ssize <- 120
sset <- sort(union(sample(which(grp=="b"),ssize/2),
                   sample(which(grp=="g"),ssize/2)))
compsset <- setdiff(1:ncol(M), sset)

# run the betaMix algorithm on all the data in order to get the 
# links between predictors. 
res <- betaMix(M, maxalpha = 1e-5, ppr=0.001,ind=T, subsamplesize = 30000)
plotFittedBetaMix(res, yLim = 10)
A <-  getAdjMat(res)
B <- as.matrix(A)*cos(res$angleMat)
B[lower.tri(B)] <- 0
gc1 <- graphComponents(A,type=1,minCtr = 3)
plotBitmapCC(A,gc1,orderByCluster = T)
summarizeClusters(gc1)
for (i in 0:max(gc1$clustNo)) {
  print(table(grp[which(gc1$clustNo ==i)]))
}

# find for each observation from the test set how many of
# its neighbors among the training set are in each class,
# and assign it the majority class.
cls1 <- rep("g", length=length(compsset))
for (i in compsset) {
  nbrs <- intersect(sset, which(A[i, ] == 1))
  if (length(nbrs) > 0) { # can get better results if using > 3
    cls1[i] <- names(which.max(table(grp[nbrs])))
  } else {
    cls1[i] <- "b"
  }
}
print(table(cls1[compsset], grp[compsset]))

cls1 <- rep("", length=ncol(A))
for (i in compsset) {
  nbrs <- intersect(sset, which(A[i, ] == 1))
  if (length(nbrs) > 0) {
    cls1[i] <- names(which.max(table(grp[nbrs])))
  } else {
    cls1[i] <- "b"
  }
}
print(table(cls1[compsset], grp[compsset]))
colseq <- rep("slateblue",length(compsset))
colseq[which(cls1[compsset] == "b")] <- "red"
# show just the test set
plot(graph.adjacency(A[compsset,compsset], mode="undirected"), 
     vertex.label=grp[compsset],
     vertex.label.cex=0.9, vertex.size=0, vertex.frame.color="white",
     vertex.label.color=colseq, edge.color='grey80',  asp=1)


colseq <- rep("grey66",ncol(M))
colseq[which(cls1 == "g")] <- "springgreen4"
colseq[which(cls1 == "b")] <- "red"
grp1 <- grp
grp1[compsset] <- toupper(grp[compsset])
# show all the data
plot(graph.adjacency(A, mode="undirected"), 
     vertex.label=grp1,
     vertex.label.cex=0.9, vertex.size=0, vertex.frame.color="white",
     vertex.label.color=colseq, edge.color='grey90',  asp=1)



# an unsupervised approach: add the class as a variable, and find the
# network:
Mall <- rbind(M, as.numeric(as.factor(grp)))
resall <- betaMix(Mall, maxalpha = 1e-5, ppr=0.001,ind=T, subsamplesize = 30000)
plotFittedBetaMix(resall, yLim = 10)
A <-  getAdjMat(res)
B <- as.matrix(A)*cos(resall$angleMat)
B[lower.tri(B)] <- 0
gc2 <- graphComponents(A,type=1,minCtr = 3)
summarizeClusters(gc2)
for (i in 1:max(gc2$clustNo)) {
  plotCluster(A,i,gc2,nodecol = Mall[33,]+2)
  print(table(grp[which(gc2$clustNo ==i)]))
}

colseq=rep("slateblue",ncol(M))
colseq[which(grp == "b")] <- "red"
plot(graph.adjacency(A, mode="undirected"), 
     vertex.label=grp,
     vertex.label.cex=0.9, vertex.size=0, vertex.frame.color="white",
     vertex.label.color=colseq, edge.color='grey80',  asp=1)

