#
# The analysis of the riboflavin data using the betaMix method
#

library(betaMix)
library(igraph)

fname = "riboflavin1.csv"
M <- read.csv(fname,sep="\t")
N <- nrow(M); P <- ncol(M)
res <- betaMix(as.matrix(M),delta = 1e-6,ppr = 1e-3, subsamplesize=50000, ind = TRUE)

#  pdf("riboflavinFit.pdf",width = 6, height = 6)
plotFittedBetaMix(res)
#  dev.off()
#  pdf("riboflavinLinks.pdf",width = 8, height = 6)
zz0 = getAdjMat(res)

b12nbrs <- c(1,which(zz0[1,]))
Msub <- zz0[b12nbrs,b12nbrs]
plot(graph.adjacency(Msub, mode="undirected"),  vertex.label= rownames(Msub),
     vertex.label.cex=0.7, vertex.size=0.1, vertex.label.color='blue',
     edge.color='gray90',  asp=.5)
#  dev.off()
