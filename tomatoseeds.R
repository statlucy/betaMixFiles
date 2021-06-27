#https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000561&StudyType=MS&ResultType=1
#
# The metabolomics data analysis in the betaMix paper
#
library("betaMix")
library("igraph")
library("xtable")

filenames <- c("ST000561_AN000862.txt", "ST000561_AN000863.txt")
group <- c("Dry Seeds", "6 hours Imbibed Seeds")

# use edgefinder and betamix
bm <- bmclustinfo <- Mat <- list()
for (grp in 1:2) {
  cat("\n\n=================\n",grp,"\n\n")
  mtbraw <- readLines(filenames[grp])
  hdr <- which(substr(mtbraw,1,7) == "Samples")
  endrow <- which(mtbraw == "MS_METABOLITE_DATA_END")
  mtbnames <- rep("",endrow-(hdr+2))
  subjects <- unlist(strsplit(mtbraw[hdr+1],"\t"))[-1]
  Mall <- matrix(0,nrow=(endrow-(hdr+2)), ncol=length(subjects))
  for (i in (hdr+2):(endrow-1)) {
    dat = unlist(strsplit(mtbraw[i],"\t"))
    mtbnames[i-hdr-1] <- dat[1]
    Mall[i-hdr-1,] <- as.numeric(dat[2:(1+length(subjects))])
  }
  
  exc <- which(substr(mtbnames,1,3) == "RI_")
  mtbnames <- mtbnames[-exc]
  Mall <- t(Mall[-exc,])
  
  M <- log2(Mall)
  Mat[[grp]] <- M
  N <- nrow(M); P <- ncol(M)
  res <- betaMix(M,ind=TRUE, delta=1e-4,ppr=0.01)
  plotFittedBetaMix(res)
  bm[[grp]] <- res
  zz0 = getAdjMat(res)
  #image(zz0)
  bmclustinfo[[grp]] <- graphComponents(zz0,minCtr = 2,type=1)
  bmclustinfo[[grp]]$labels <- mtbnames
  for (cls in 1:max(bmclustinfo[[grp]]$clustNo)) {
    plotCluster(zz0,cls,bmclustinfo[[grp]], labels=TRUE)
  }
  #print(summarizeClusters(bmclustinfo[[grp]]))
  plotBitmapCC(zz0,bmclustinfo[[grp]],orderByCluster = TRUE)
  plot(graph.adjacency(zz0, mode="undirected"),  vertex.label= mtbnames,
       vertex.label.cex=0.7, vertex.size=0.1, vertex.label.color='blue',
       edge.color='green',  asp=.5)
}

print(summarizeClusters(bmclustinfo[[1]])[,1:9])
print(summarizeClusters(bmclustinfo[[2]])[,1:9])

zz1 = Matrix(sin(bm[[1]]$angleMat)^2 < bm[[1]]$ppthr)
zz2 = Matrix(sin(bm[[2]]$angleMat)^2 < bm[[2]]$ppthr)
diag(zz1) <- diag(zz2) <- FALSE

mtbnames2 = mtbnames
mtbnames2[which(substr(mtbnames2,1,3) == "RI_")] = ""

par(mfrow=c(1,2))
plot(graph.adjacency(zz1, mode="undirected"),  vertex.label= mtbnames2,
     vertex.label.cex=0.8, vertex.size=0.1, vertex.label.color='blue',
     edge.color='green',  asp=1)
plot(graph.adjacency(zz2, mode="undirected"),  vertex.label= mtbnames2,
     vertex.label.cex=0.8, vertex.size=0.1, vertex.label.color='brown',
     edge.color='gray',  asp=1)
par(mfrow=c(1,1))

ttests <- unlist(apply(Mat[[1]]-Mat[[2]],1,t.test))
FoldChange <- as.numeric(ttests[which(names(ttests) == "estimate.mean of x")])
sigchanges <- intersect(which(p.adjust(as.numeric(ttests[which(names(ttests) == "p.value")]),
                                       method = "holm") < 0.01),
                        which(abs(FoldChange) > 1))
dfFC <- data.frame(mtbnames[sigchanges],FoldChange[sigchanges])
colnames(dfFC) <- c("Metabolite","FoldChange")
dfFC <- dfFC[order(dfFC$FoldChange),]
xtable(dfFC)


for (cls in 1:max(bmclustinfo[[2]]$clustNo)) {
  ids <- which(bmclustinfo[[2]]$clustNo == cls)
  tmp2 <- zz2[ids, ids]
  tmp1 <- zz1[ids, ids]
  thetas <- seq(0,2*pi,length=(length(ids)+1))
  posx <- cos(thetas)
  posy <- sin(thetas)
  plot(posx, posy, cex = 0, col=0,  axes = F, xlab = "", ylab = "",
       xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))
  text(posx, posy, bmclustinfo[[2]]$labels[ids], pos = ifelse(posy>0,3,1), cex = 0.8, srt=30)
  for (i in 1:ncol(tmp2)) {
    nbrs2 <- setdiff(which(tmp2[i, ] == 1), 1:i)
    if (length(nbrs2) > 0) {
      for (j in nbrs2) {
        lines(c(posx[i], posx[j]), c(posy[i], posy[j]), col = "grey80",
              lwd = 1)
      }
      nbrs1 <- setdiff(which(tmp1[i, ] == 1), 1:i)
      if (length(nbrs1) > 0) {
        for (j in nbrs1) {
          lines(c(posx[i], posx[j]), c(posy[i], posy[j]), col = "orange",
                lwd = 2, lty=2)
        }
      }
    }
  }
  points(posx, posy, cex = ifelse(ids %in% sigchanges,1.5,0.75), 
         pch = ifelse(FoldChange[ids]>0, 24,25), 
         col = ifelse(ids %in% sigchanges,"red","blue"))
  
}

#pdf("tomatoClusters.pdf", width=6, height = 6)
par(mfrow=c(2,2), mar=c(0.5,0.5,1,0))
for (cls in 1:4) {
  ids <- which(bmclustinfo[[2]]$clustNo == cls)
  tmp2 <- zz2[ids, ids]
  tmp1 <- zz1[ids, ids]
  thetas <- seq(0,2*pi,length=(length(ids)+1))
  posx <- cos(thetas)
  posy <- sin(thetas)
  plot(posx, posy, cex = 0, col=0,  axes = F, xlab = "", ylab = "",
       xlim=c(-1.3,1.3), ylim=c(-1.3,1.3))
  for (i in 1:ncol(tmp2)) {
    nbrs2 <- setdiff(which(tmp2[i, ] == 1), 1:i)
    if (length(nbrs2) > 0) {
      for (j in nbrs2) {
        lines(c(posx[i], posx[j]), c(posy[i], posy[j]), col = "grey80",
              lwd = 1)
      }
      nbrs1 <- setdiff(which(tmp1[i, ] == 1), 1:i)
      if (length(nbrs1) > 0) {
        for (j in nbrs1) {
          lines(c(posx[i], posx[j]), c(posy[i], posy[j]), col = "orange",
                lwd = 2, lty=2)
        }
      }
    }
  }
  text(posx, posy, bmclustinfo[[2]]$labels[ids], pos = ifelse(posy>0,3,1), cex = 0.8, srt=30)
  points(posx, posy, cex = ifelse(ids %in% sigchanges,1.5,0.75), 
         pch = ifelse(FoldChange[ids]>0, 24,25), 
         col = ifelse(ids %in% sigchanges,"red","blue"))
  
}
#dev.off()
