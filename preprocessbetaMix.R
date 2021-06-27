# Preparing the riboflavin data to be loaded to sqlite and demonstrating the two possible methods (in-memory, and SQL-based)

switch (codeChunk,
        '1' = {
          library("betaMix")
          dat <- read.csv("riboflavin1.csv", sep="\t")
          #dim(dat) 
        }, 
        '2' = {
          pdf("tmp/histB2.pdf",width = 4, height = 4)
          hist(dat[,1],col="grey66",border="white", xlab="log vitamin B2 production rate", main="") 
          dev.off()
        },
        '3' = {
          #### In-memory method
          resMem <- betaMix(dat,ind=TRUE, subsamplesize = 100000, delta = 1/choose(ncol(dat),2),ppr=0.001, msg=FALSE)
          pdf("tmp/fittedInMem.pdf",width = 4, height = 4)
          plotFittedBetaMix(resMem,yLim=30)
          dev.off()
          shortSummary(resMem)
          adjMat1 <- getAdjMat(resMem)
          #length(which(adjMat1[1,]  > 0))
        },
        '4' = {
          #### SQL-based method
          outdir <- "tmp" # where the batch-load file and the SQL database will be saved
          P <- ncol(dat)
          n <- nrow(dat)
          zeroVar <- 1e-5 # The threshold which will be used to decide if the variance of a vector is effectively 0.
          dbname <- "riboflavin" # A database name
          description <- "The riboflavin data, B{\"u}hlmann, Kalisch, and Meier 2014"
          fin <- sprintf("%s/%s_fmt.txt",outdir,dbname) # the re-formatted data file which will be used by calcCorr.
          if(file.exists(fin))
            file.remove(fin) # make sure we start a new file
          MaxNameLen <- max(nchar(colnames(dat))) + 1
          cat(MaxNameLen,"\n")
          dat2 <- t(dat)
          # Create a fixed-width file:
          nchr <- rep(0,nrow(dat2))
          for (i in 1:nrow(dat2)) {
            numval <- paste0(formatC(dat2[i,],width=16,digits=10), collapse = "")
            cat(paste0(sprintf("%-13s",rownames(dat2)[i]),numval,collapse = ""),"\n",
                file = fin, append = T, sep = "")
            nchr[i] <- nchar(paste0(sprintf("%-13s",rownames(dat2)[i]),numval))
          }
          cat(unique(nchr)+1,"\n")
          # all lines must have the same length! If not, stop here and check.
          if (length(unique(nchr)) == 1)
            recSize <- unique(nchr)+1
          
          cat(fin,P,n,zeroVar,format(Sys.time(), "%Y%M%d%H%M%S"),
              description,sep="|", file=sprintf("%s/%s_meta.txt",outdir,dbname))
          cat(paste(1:P,rownames(dat2), sep ="|"),sep="\n", file=sprintf("%s/%s_names.txt", outdir,dbname))
          
          fout <- sprintf("%s/%s_cors.txt",outdir,dbname)
          calcCorr(c(fin,fout), MaxNameLen, recSize, P,  n)
        },
        '5' = {
          cat(sprintf("CREATE TABLE metadata ( inputfile TEXT NOT NULL,  p INTEGER, n INTEGER, zerovar DOUBLE, date TEXT, description TEXT);
.import  %s/%s_meta.txt metadata
CREATE TABLE predictors (num INTEGER, name TEXT NOT NULL);
.import  %s/%s_names.txt predictors
CREATE TABLE correlations (node1 INTEGER, node2 INTEGER,corr DOUBLE,zij DOUBLE);
.import %s/%s_cors.txt correlations\n",outdir,dbname,outdir,dbname,outdir,dbname), file=sprintf("%s/sqlbatch",outdir))
        },
        '6' = {
          dbfile <- sprintf("%s/%s",outdir,dbname)
          resSQL <- betaMix(dbname=dbfile, ind=TRUE, subsamplesize = 100000, delta = 1/choose(ncol(dat),2), ppr=0.001, msg=FALSE)
          adjMat2 <- getAdjMat(resSQL,dbname=dbfile)
        },
        '7' = {
          b2nbrs <- c(1, which(adjMat2[1,]>0))
          adjMat3 <- getAdjMat(resSQL, dbname=dbfile, nodes = b2nbrs)
          Msub <- adjMat3[b2nbrs,b2nbrs]
          library(igraph)
          pdf("tmp/B2nbrs.pdf",width = 6, height = 6)
          plot(graph.adjacency(Msub, mode="undirected"),  vertex.label= rownames(Msub),  vertex.label.cex=0.7, vertex.size=0.1, vertex.label.color='blue', edge.color='grey90',  asp=1)
          dev.off()
        },
        cat("default\n")
)
