library(igraph)
library("betaMix")
load("./congress.RData")

# use betaMix to find network of congress members by voting patterns.
# each year is handled separately, since members of congress are up for vote every two years

year <- 2008
datYear <- dat[which(dat$year == year),]
datYear$vote_legislator.text = droplevels(datYear$vote_legislator.text)
# table(datYear$vote_legislator.text)
# most members have the same number of roll calls but some don't, and are excluded if not
# participated in too many votes

datYear$vote_result <- droplevels(datYear$vote_result)
print(table(datYear$vote_result))
# can be null, Aye, Nay, No, Not Voting, Present, Yea (or )

cat("total roll calls", length(unique(datYear$roll_number)),"\n")
rolls <- sort(unique(datYear$roll_number))
issues <- rep("",length(rolls))
members <- sort(unique(datYear$vote_legislator.text))[-1] # remove "" as name

# build the voting matrix, members in the rows, rolls in the columns
M <- matrix(0,nrow=length(members),ncol=length(rolls))
rownames(M) <- members
party <- rep("",length(members))
state <- rep("",length(members))

for (i in 1:length(members)) {
  party[i] <- as.character(datYear$vote_legislator.party[which(datYear$vote_legislator.text == members[i])][1])
  state[i] <- as.character(datYear$vote_legislator.state[which(datYear$vote_legislator.text == members[i])][1])
}
print(table(party))
print(table(state))
majority <- ifelse(sum(party == "D") > sum(party == "R"),"D","R")

for (roll in 1:length(rolls)) {
  rollrows <- which(datYear$roll_number == roll)
  issues[roll] <- as.character(datYear$issue[rollrows[1]])
  rollnumeric <- rep(0,length(rollrows))
  votingmems <- which(members %in% datYear$vote_legislator.text[rollrows])
  if (length(votingmems) < 100)
    next
  rollnumeric[which(tolower(datYear$vote_result[rollrows]) %in% c("aye","yea","yes"))] <- 0.5
  rollnumeric[which(tolower(datYear$vote_result[rollrows]) %in% c("nay","no"))] <- -0.5
  M[votingmems,roll] <- rollnumeric
}

# remove unanimous votes
excCol <- which(apply(M,2,sd) == 0)
if (length(excCol) > 0) {
  M <- M[,-excCol]
}

# remove membrs who have voted less than 10% of the times
excRow <- which(apply(M,1,quantile, 0.7) == 0)
if (length(excRow) > 0) {
  M <- M[-excRow,]
  members <- members[-excRow]
  party <- party[-excRow]
  cat("Excluding:\n")
  for (j in 1:length(excRow)) {
    cat("\t",as.character(members[excRow[j]]),party[excRow[j]],state[excRow[j]],"\n")
  }
} 

issues <- issues[-excCol]
# some issues correspond to multiple votes. Keep just one vote for each:
M <- M[,-which(duplicated(issues))]
issues <- issues[-which(duplicated(issues))]
# and remove the following issues:
exc <- which(issues %in% c("ADJOURN","JOURNAL","MOTION"))
if (length(exc) > 0) {
  M <- M[,-exc]
  issues <- issues[-exc]
}

cat("total issues", length(unique(issues)),"\n")

M <- t(M) # put the congress members in the columns and the issues in rows
# run the algorithm:
res <- betaMix(M, delta = 1e-4, ppr=0.001)
plotFittedBetaMix(M,res, yLim = 16)
A <-  Matrix(sin(res$angleMat)^2 < res$ppthr)
diag(A) <- FALSE
gc1 <- graphComponents(A)
summarizeClusters(gc1) # two large blocks

for (cnum in 1:max(gc1$clustNo)) {
  clscomp <- table(party[which(gc1$clustNo == cnum)])
  print(clscomp)
  gctmp <- gc1
  nodecols <- ifelse(party == "R","red","blue")
}

Alpha <- rep(0.6, ncol(M))
vrtxsize <- rep(0.3,ncol(M))
vrtxsize[which(gc1$clustNo == 0)] <- 0.9
Alpha[which(gc1$clustNo == 0)] <- 1
misfits <- union(intersect(which(gc1$clustNo == 1), which(party == "R")),
                 intersect(which(gc1$clustNo == 2), which(party == "D")))
vrtxsize[misfits] <- 0.9
Alpha[misfits] <- 1

vrtxlabel <- party
vrtxlabel[which(gc1$clustNo == 0)] <- paste0(members[which(gc1$clustNo == 0)]," (",party[which(gc1$clustNo == 0)],")")
colseq <- rgb(0, 0.2, 1, Alpha)
colseq[which(party == "R")] <- rgb(1, 0, 0, Alpha[which(party == "R")])

reord <- c(setdiff(1:ncol(M), misfits), misfits)

plot(graph.adjacency(A[reord,reord], mode="undirected"), 
     vertex.label=vrtxlabel[reord], 
     vertex.label.cex=vrtxsize[reord], vertex.size=0, vertex.frame.color="white",
     vertex.label.color=colseq[reord], edge.color='grey80',  asp=1)
