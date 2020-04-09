#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))

parser <- ArgumentParser()

parser$add_argument("-a", "--alpha", type = "double", default = 0.5,
                    help = "shape parameter of gamma distribution for simulating PRC2 movement")
parser$add_argument("-b", "--beta", type = "double", default = 0.001,
                    help = "rate parameter of gamma distribution for simulating PRC2 movement; smaller beta allows longer-range movement of PRC2")
parser$add_argument("--regionLength", type = "integer", default = 5000,
                    help = "length of simulated genomic region")
parser$add_argument("--PRC2Peak1", type = "integer", default = c(1490,1510), nargs = '+',
                    help = "start/end coordinates indicating where a PRC2 binding peak exists")
parser$add_argument("--PRC2Peak2", type = "integer", default = c(3190,3210), nargs = '+',
                    help = "start/end coordinates where a PRC2 binding peak exists")
parser$add_argument("--pme1", type = "double", default = 0.9, 
                    help = "normal conversion rate from me0 to me1")
parser$add_argument("--pme2", type = "double", default = 0.25, 
                    help = "normal conversion rate from me1 to me2")
parser$add_argument("--pme2adj", type = "double", default = 0.5, 
                    help = "increased conversion rate from me1 to me2 by allosteric activation")
parser$add_argument("--pme3", type = "double", default = 0.01, 
                    help = "normal conversion rate from me2 to me3")
parser$add_argument("--pme3adj", type = "double", default = 0.15, 
                    help = "increased conversion rate from me2 to me3 by allosteric activation")
parser$add_argument("--enhancedTime", type = "integer", default = 2, 
                    help = "time for allosterically activated PRC2 to maintain its enhanced catalytic ability; set to 0 to remove this effect")
parser$add_argument("-N", "--Npop", type = "integer", default = 100, 
                    help = "number of cells")
parser$add_argument("--Nenzyme1", type = "integer", default = 8, 
                    help = "number of PRC2 molecules initially bound at peak1")
parser$add_argument("--Nenzyme2", type = "integer", default = 8, 
                    help = "number of PRC2 molecules initially bound at peak2")
parser$add_argument("--Nstep", type = "integer", default = 5000, 
                    help = "total simulation time")
parser$add_argument("--period", type = "integer", default = 1000,
                    help = "cell division interval")
parser$add_argument("--TSSTES", type = "integer", default = c(2100,2300,3500,3750), nargs = '+', 
                    help = "start/end coordinates for actively transcribed regions")
parser$add_argument("--K36Me3", type = "integer", default = c(2100,2300,3500,3750), nargs = '+', 
                    help = "start/end coordinates for regions with K36me3 marks; generally overlapping with actively transcribed regions")
parser$add_argument("--K36Me2", type = "integer", default = c(2301,2900), nargs = '+',
                    help = "start/end coordinates for regions with K36me2 marks")
parser$add_argument("--K36me3harm01", type = "double", default = 5, 
                    help = "K36me3 decreases me0 to me1 conversion rate by K36me3harm01-fold")
parser$add_argument("--K36me3harm12", type = "double", default = 7.5, 
                    help = "K36me3 decreases me1 to me2 conversion rate by K36me3harm12-fold")
parser$add_argument("--K36me3harm23", type = "double", default = 15, 
                    help = "K36me3 decreases me2 to me3 conversion rate by K36me3harm23-fold")
parser$add_argument("--K36me2harm01", type = "double", default = 2, 
                    help = "K36me2 decreases me0 to me1 conversion rate by K36me2harm01-fold")
parser$add_argument("--K36me2harm12", type = "double", default = 3, 
                    help = "K36me2 decreases me1 to me2 conversion rate by K36me2harm12-fold")
parser$add_argument("--K36me2harm23", type = "double", default = 5, 
                    help = "K36me2 decreases me2 to me3 conversion rate by K36me2harm23-fold")
parser$add_argument("-E", "--Expression", type = "logical", default = FALSE, 
                    help = "whether to include removal of methylated histones in actively transcribed regions upon active gene expression")
parser$add_argument("--periodExpression", type = "integer", default = 100, 
                    help = "gene expression interval")
parser$add_argument("--mutRate", type = "double", default = 0.05, 
                    help = "rate of K27M mutants; set to 0 to simulate a K27M-free scenario")
parser$add_argument("--K27Mharm01", type = "double", default = 2, 
                    help = "K27M decreases me0 to me1 conversion rate by K27Mharm01-fold")
parser$add_argument("--K27Mharm12", type = "double", default = 5, 
                    help = "K27M decreases me1 to me2 conversion rate by K27Mharm12-fold")
parser$add_argument("--K27Mharm23", type = "double", default = 10, 
                    help = "K27M decreases me2 to me3 conversion rate by K27Mharm23-fold")
parser$add_argument("--K27MharmPermanent", type = "logical", default = TRUE, 
                    help = "whether K27M permanently damages the catalytic capacity of PRC2 and prevents future allosteric activation")
parser$add_argument("--sequestrationTime", type = "integer", default = 5, 
                    help = "time to detain PRC2 on K27M; set to 0 to remove this effect")
parser$add_argument("--normalBackRate", type = "double", default = 0, 
                    help = "probability that a PRC2 automatically dissociates from the genome and gets recruited back to the nucleation site")
parser$add_argument("--K27MBackRate", type = "double", default = 0, 
                    help = "increased probability that a PRC2 returns to the nucleation site when the molecule encounters K27M; set to 0 to remove this effect")
parser$add_argument("--targetedRegion", type = "integer", default = c(1001,4000), nargs = '+', 
                    help = "start/end coordinates for the region used for generating summary statistics; cropping out regions close to the boundaries")
parser$add_argument("--shortRangeProb", type = "double", default = 0.95, 
                    help = "probability that a PRC2 takes small steps; set to 0 to remove short-range movement")
parser$add_argument("--shortRangeLimit", type = "integer", default = 8, 
                    help = "maximum distance a PRC2 molecule can travel when it takes small steps")
parser$add_argument("--shortRangeCentrality", type = "double", default = 0.5, 
                    help = "centrality of the short-range movement; probability of a PRC2 remains at its current locus is shortRangeCentrality to the power of shortRangeLimit")
parser$add_argument("--K27MintroduceTime", type = "integer", default = 0, 
                    help = "time point when K27M is introduced; set to 0 to introduce the mutants from the beginning of simulation")
parser$add_argument("--equiPeriod", type = "integer", default = 1,
                    help = "number of last time point(s) to be considered as equilibrium stage; proportion of marks will be reported by averaging across these time points. E.g. set to 1 to report for only the last time point; set to 100 to report average proportion across the last 100 time points")
parser$add_argument("-o", "--output", type = "character", default = "ModificationPattern",
                    help = "prefix for output file names")


args <- parser$parse_args()

args$alpha -> alpha 
args$beta -> beta
args$regionLength -> regionLength
args$PRC2Peak1[1]:args$PRC2Peak1[2] -> PRC2Peak1
args$PRC2Peak2[1]:args$PRC2Peak2[2] -> PRC2Peak2
args$pme1 -> pme1
args$pme2 -> pme2
args$pme2adj -> pme2adj
args$pme3 -> pme3
args$pme3adj -> pme3adj 
args$enhancedTime -> enhancedTime
args$Npop -> Npop
args$Nenzyme1 -> Nenzyme1
args$Nenzyme2 -> Nenzyme2
args$Nstep -> Nstep
args$period -> period
unlist(lapply(seq(2,length(args$TSSTES),by=2),function(x) args$TSSTES[x-1]:args$TSSTES[x])) -> TSSTES
unlist(lapply(seq(2,length(args$K36Me3),by=2),function(x) args$K36Me3[x-1]:args$K36Me3[x])) -> K36Me3
unlist(lapply(seq(2,length(args$K36Me2),by=2),function(x) args$K36Me2[x-1]:args$K36Me2[x])) -> K36Me2
args$K36me3harm01 -> K36me3harm01
args$K36me3harm12 -> K36me3harm12
args$K36me3harm23 -> K36me3harm23
args$K36me2harm01 -> K36me2harm01
args$K36me2harm12 -> K36me2harm12
args$K36me2harm23 -> K36me2harm23
args$Expression -> Expression
args$periodExpression -> periodExpression
args$mutRate -> mutRate
args$K27Mharm01 -> K27Mharm01
args$K27Mharm12 -> K27Mharm12 
args$K27Mharm23 -> K27Mharm23
args$K27MharmPermanent -> K27MharmPermanent
args$sequestrationTime -> sequestrationTime
args$normalBackRate -> normalBackRate
args$K27MBackRate -> K27MBackRate
unlist(lapply(seq(2,length(args$targetedRegion),by=2),function(x) args$targetedRegion[x-1]:args$targetedRegion[x])) -> targetedRegion
args$shortRangeProb -> shortRangeProb
args$shortRangeLimit -> shortRangeLimit
args$shortRangeCentrality -> shortRangeCentrality
args$K27MintroduceTime -> K27MintroduceTime
args$equiPeriod -> equiPeriod
args$output -> filename

mutantRegion <- lapply(1:Npop, function(x) sample(seq(1:regionLength),round(mutRate * regionLength))) ### histones containing H3K27M are randomly distributed

pme1weak <- pme1/K36me3harm01 ## me0 -> me1
pme1weak2 <- pme1/K36me2harm01 ## me0 -> me1

posList <- list()
PMe2List <- list()
PMe3List <- list()
simulateMoveK27M <- function(POS,NSTEP,PMe2,PMe3,leftBound,rightBound,CELL) {
  ori <- POS
  oriPMe2 <- PMe2
  oriPMe3 <- PMe3
  POSout <- rep(0, NSTEP + 1)
  PMe2out <- rep(0, NSTEP + 1)
  PMe3out <- rep(0, NSTEP + 1)
  POSout[1] <- ori
  PMe2out[1] <- oriPMe2
  PMe3out[1] <- oriPMe3
  for (t in 1:NSTEP) {
    if (rbinom(1,1,shortRangeProb)) {
      POSout[t+1] <- POSout[t]+round(rbinom(1,shortRangeLimit,1 - shortRangeCentrality))*ifelse(rbinom(1,1,0.5),1,-1)
    }
    else {
      POSout[t+1] <- POSout[t]+round(rgamma(1,alpha,beta))*ifelse(rbinom(1,1,0.5),1,-1)
    }
    PMe2out[t+1] <- PMe2out[t]
    PMe3out[t+1] <- PMe3out[t]
    if (POSout[t+1]<leftBound|POSout[t+1]>rightBound) { # If a PRC2 molecule gets out of the region, a new molecule will bind to the PRC2 binding site and will have a 50% chance of having full enzymatic capacity
      PMe2out[t+1] <- sample(c(PMe2out[t],oriPMe2),1)
      PMe3out[t+1] <- sample(c(PMe3out[t],oriPMe3),1)
    }
    if (POSout[t+1] %in% mutantRegion[[CELL]] & PMe2out[t+1]==oriPMe2 & t >= K27MintroduceTime) {
      PMe2out[t+1] <- oriPMe2/K27Mharm12                                
      PMe3out[t+1] <- oriPMe3/K27Mharm23
    }
    if (POSout[t+1] %in% mutantRegion[[CELL]] & t >= K27MintroduceTime) {
      POSout[t+1] <- ifelse(rbinom(1,1,K27MBackRate),ori,POSout[t+1])
    }
    POSout[t+1] <- ifelse(rbinom(1,1,normalBackRate),ori,POSout[t+1])
    POSout[t+1] <- ifelse(POSout[t+1]<leftBound,ori,ifelse(POSout[t+1]>rightBound,ori,POSout[t+1]))
  }
  return(c(POSout,PMe2out,PMe3out))
}
countnum <- function(ls, num){ # function to count number in vector equal to "num"
  length(ls[ls==num])
}
divisionRound <- sample(1:period,Npop,replace = TRUE) # cell division does not occur at the same time for the population
divisionSeries <- c(divisionRound)
for (T in 1:(Nstep/period-1)) {
  divisionSeries <- rbind(divisionSeries, divisionRound+T*period) # each cell divides every 2,000 time units
}

if (Expression) {
  geneExpressionRound <- sample(1:periodExpression,Npop,replace = TRUE) # gene expression does not occur at the same time for the population
  expressionSeries <- c(geneExpressionRound)
  for (T in 1:(Nstep/periodExpression-1)) {
    expressionSeries <- rbind(expressionSeries, geneExpressionRound+T*periodExpression) # each cell expresses gene every 100 time units
  }
}

for (cell in 1:Npop) {
  cat("simulating PRC2 activity for cell:",cell,collapse="\n")
  if (length(PRC2Peak1)==1) {
    peak1 <- rep(PRC2Peak1,Nenzyme1)
  }
  else {
    peak1 <- sample(PRC2Peak1,Nenzyme1,replace = T)
  }
  if (length(PRC2Peak2)==1) {
    peak2 <- rep(PRC2Peak2,Nenzyme2)
  }
  else {
    peak2 <- sample(PRC2Peak2,Nenzyme2,replace = T)
  }
  pos <- matrix(c(peak1,peak2),1,(Nenzyme1 + Nenzyme2))
  temporaryRes <- apply(pos, 2, function(x) simulateMoveK27M(x,Nstep,pme2,pme3,1,regionLength,cell))
  posList[[cell]] <- temporaryRes[1:(Nstep+1),]
  PMe2List[[cell]] <- temporaryRes[(Nstep+2):(2*Nstep+2),]
  PMe3List[[cell]] <- temporaryRes[(2*Nstep+3):(3*Nstep+3),]
}

checkpoint <- list()
checkpointDistribution <- list()

for (round in 1:Npop) {
  cat("simulating H3K27me dynamics for cell:",round,collapse="\n")
  pos <- posList[[round]]
  PME2 <- PMe2List[[round]]
  PME3 <- PMe3List[[round]]
  posBase <- pos
  genomicRegion <- rep(0,regionLength+2) # add one "pseudo"-position at each end which will be deleted eventually; allow detection for adjacent me3 on both edges
  genomicRegionHarbor <- c()
  markDistribution <- c()
  for (T in 1:(Nstep+1)) {
    if (T %in% divisionSeries[,round]) {
      pos[T,] <- ifelse(rbinom(Nenzyme1 + Nenzyme2,1,0.5),pos[T,],pos[1,])
      genomicRegion[-unique(pos[T,])+1] <- ifelse(rbinom(length(genomicRegion[-unique(pos[T,])+1]),1,0.5),genomicRegion[-unique(pos[T,])+1],0)
    }
    if (Expression) {
      if (T %in% expressionSeries[,round]) {
        genomicRegion[TSSTES + 1] <- ifelse(rbinom(length(genomicRegion[TSSTES+1]),1,0.5),0,genomicRegion[TSSTES + 1])
        genomicRegion[TSSTES + 1] <- ifelse(rbinom(length(genomicRegion[TSSTES+1]),1,0.5),0,genomicRegion[TSSTES + 1])
      }
    }
    
    for (reg in 1:length(pos[T,])) {
      if (genomicRegion[posBase[T,reg]]==3 | genomicRegion[posBase[T,reg]+2]==3) {
        if (K27MharmPermanent) {
          if (PME2[T,reg] == pme2 & PME3[T,reg] == pme3) {
            PME2[T:min((T+enhancedTime),Nstep+1),reg] <- pme2adj
            PME3[T:min((T+enhancedTime),Nstep+1),reg] <- pme3adj
          }
        }
        else {
          if (PME2[T,reg] == pme2 & PME3[T,reg] == pme3) {
            PME2[T:min((T+enhancedTime),Nstep+1),reg] <- pme2adj
            PME3[T:min((T+enhancedTime),Nstep+1),reg] <- pme3adj
          }
          else {
            if (posBase[T,reg] %in% mutantRegion[[round]] & T >= K27MintroduceTime) {
              PME2[T:min((T+enhancedTime),Nstep+1),reg] <- pme2adj / K27Mharm12
              PME3[T:min((T+enhancedTime),Nstep+1),reg] <- pme3adj / K27Mharm23
            }
          }
        }
      }
      if (posBase[T,reg] %in% mutantRegion[[round]] & T >= K27MintroduceTime) {
        pos[T:min((T+sequestrationTime),Nstep+1),reg] <- pos[T,reg]
      }
      if (pos[T,reg] %in% K36Me3) {
        if (genomicRegion[pos[T,reg]+1]==0) {
          genomicRegion[pos[T,reg]+1]<-genomicRegion[pos[T,reg]+1]+rbinom(1,1,pme1weak / ifelse(pos[T,reg] %in% mutantRegion[[round]] & T >= K27MintroduceTime,K27Mharm01,1))
        }
        if (genomicRegion[pos[T,reg]+1]==1) {
          genomicRegion[pos[T,reg]+1]<-genomicRegion[pos[T,reg]+1]+rbinom(1,1,PME2[T,reg]/K36me3harm12)
        }
        if (genomicRegion[pos[T,reg]+1]==2) {
          genomicRegion[pos[T,reg]+1]<-genomicRegion[pos[T,reg]+1]+rbinom(1,1,PME3[T,reg]/K36me3harm23)
        }
      } 
      if (pos[T,reg] %in% K36Me2) {
        if (genomicRegion[pos[T,reg]+1]==0) {
          genomicRegion[pos[T,reg]+1]<-genomicRegion[pos[T,reg]+1]+rbinom(1,1,pme1weak2 / ifelse(pos[T,reg] %in% mutantRegion[[round]] & T >= K27MintroduceTime,K27Mharm01,1))
        }
        if (genomicRegion[pos[T,reg]+1]==1) {
          genomicRegion[pos[T,reg]+1]<-genomicRegion[pos[T,reg]+1]+rbinom(1,1,PME2[T,reg]/K36me2harm12)
        }
        if (genomicRegion[pos[T,reg]+1]==2) {
          genomicRegion[pos[T,reg]+1]<-genomicRegion[pos[T,reg]+1]+rbinom(1,1,PME3[T,reg]/K36me2harm23)
        }
      } 
      if (!pos[T,reg] %in% K36Me3 & !pos[T,reg] %in% K36Me2) {
        if (genomicRegion[pos[T,reg]+1]==0) {
          genomicRegion[pos[T,reg]+1]<-genomicRegion[pos[T,reg]+1]+rbinom(1,1,pme1 / ifelse(pos[T,reg] %in% mutantRegion[[round]] & T >= K27MintroduceTime,K27Mharm01,1))
        }
        if (genomicRegion[pos[T,reg]+1]==1) {
          genomicRegion[pos[T,reg]+1]<-genomicRegion[pos[T,reg]+1]+rbinom(1,1,PME2[T,reg])
        }
        if (genomicRegion[pos[T,reg]+1]==2) {
          genomicRegion[pos[T,reg]+1]<-genomicRegion[pos[T,reg]+1]+rbinom(1,1,PME3[T,reg])
        }
      }
    }
    if (T >= K27MintroduceTime) {
      genomicRegion[mutantRegion[[round]] + 1] <- 0
    }
    wholeGenomicRegion <- genomicRegion[-c(1,regionLength+2)]
    wholeGenomicRegion <- wholeGenomicRegion[targetedRegion]
    markDistribution <- rbind.data.frame(markDistribution,
                                         c(countnum(wholeGenomicRegion,0),
                                           countnum(wholeGenomicRegion,1),
                                           countnum(wholeGenomicRegion,2),
                                           countnum(wholeGenomicRegion,3)))
    if ((T-1)%%100==0) {
      genomicRegionHarbor <- rbind(genomicRegionHarbor, genomicRegion)
    }
    if ((T-1)%%1000==0 & T > 1) {
      cat("processed",T-1,"iterations",collapse="\n")
    }
  }
  genomicRegionHarbor <- genomicRegionHarbor[,-c(1,regionLength+2)]
  genomicRegionHarbor <- genomicRegionHarbor[,targetedRegion]
  checkpoint[[round]] <- genomicRegionHarbor
  checkpointDistribution[[round]] <- markDistribution
  posList[[round]] <- pos
}

pdf(paste0(filename,".PRC2.occupancy.pdf"))
hist(unlist(posList),breaks = regionLength,xlab = "",main = "",ylab = "Cumulative occupancy of PRC2")
dev.off()

plotDistFun <- function(check=(Nstep/100+1)) {
  sumAllCell <- matrix(nrow = length(targetedRegion), ncol = Npop)
  for (cell in 1:Npop) {
    sumAllCell[,cell] <- t(checkpoint[[cell]][check,]) # remove the 1st and last position which are always 0
  }
  me1 <- apply(sumAllCell, 1, function(x) countnum(x,1)) # store num of me1 at each position in a vector
  me2 <- apply(sumAllCell, 1, function(x) countnum(x,2)) # store num of me2 at each position in a vector
  me3 <- apply(sumAllCell, 1, function(x) countnum(x,3)) # store num of me3 at each position in a vector
  par(mfrow=c(3,1))
  barplot(me3,ylim = c(0,Npop),xlab = "H3K27me3",main = paste0("T = ",check*100-99))
  barplot(me2,ylim = c(0,Npop),xlab = "H3K27me2")
  barplot(me1,ylim = c(0,Npop),xlab = "H3K27me1")
}
pdf(paste0(filename,".final.pattern.pdf"))
plotDistFun((Nstep/100+1)) # visualize methylation pattern at specified time point [taking values from 1 to (total time/100 + 1)]. Default is the end-stage.
dev.off()

#par(mfrow=c(1,1))
## plot change in total number of marks
newDist <- checkpointDistribution[[1]] + checkpointDistribution[[2]]
for (j in 3:Npop) {
  newDist <- newDist + checkpointDistribution[[j]]
}
#newDist <- checkpointDistribution[[1]] # if there is only one cell, use this sentence instead
colnames(newDist) <- c("Me0","Me1","Me2","Me3")
newDist$Time <- seq(1,Nstep+1)
newDat <- data.frame(`Number of methylated histones`=newDist$Me0,Time=newDist$Time,Mark="Unmethylated H3K27")
newDat <- rbind(newDat, data.frame(`Number of methylated histones`=newDist$Me1,Time=newDist$Time,Mark="H3K27me1"))
newDat <- rbind(newDat, data.frame(`Number of methylated histones`=newDist$Me2,Time=newDist$Time,Mark="H3K27me2"))
newDat <- rbind(newDat, data.frame(`Number of methylated histones`=newDist$Me3,Time=newDist$Time,Mark="H3K27me3"))
newDat$Number.of.methylated.histones <- newDat$Number.of.methylated.histones / Npop / length(targetedRegion)

pdf(paste0(filename,".cumulative.amount.pdf"))
ggplot(newDat,aes(x=Time,y=Number.of.methylated.histones,colour=Mark)) + geom_line() + theme_bw() + ylab("Fraction")
dev.off()
tmp1 <- newDist[1:Nstep,1:4]
tmp2 <- newDist[2:(Nstep+1),1:4]
firstDerive <- tmp2-tmp1 # get first derivative
firstDerive$Time <- seq(1,Nstep)
newFD <- data.frame(`Change of methylated histones`=firstDerive$Me0,Time=firstDerive$Time,Mark="Unmethylated H3K27")
newFD <- rbind(newFD, data.frame(`Change of methylated histones`=firstDerive$Me1,Time=firstDerive$Time,Mark="H3K27me1"))
newFD <- rbind(newFD, data.frame(`Change of methylated histones`=firstDerive$Me2,Time=firstDerive$Time,Mark="H3K27me2"))
newFD <- rbind(newFD, data.frame(`Change of methylated histones`=firstDerive$Me3,Time=firstDerive$Time,Mark="H3K27me3"))

pdf(paste0(filename,".deposition.rate.pdf"))
p <- ggplot(newFD,aes(x=Time,y=Change.of.methylated.histones,colour=Mark)) + geom_smooth() + theme_bw() + ylab("Change in number of histones")
suppressMessages(print(p))
dev.off()

cat("\n")
cat("########",collapse = "\n")
cat("proportion of marks at equilibrium:", collapse="\n")
cat("H3K27me0","H3K27me1","H3K27me2","H3K27me3",collapse="\n",sep = "\t")
cat(tapply(newDat$Number.of.methylated.histones[newDat$Time > max(newDat$Time) - equiPeriod],newDat$Mark[newDat$Time > max(newDat$Time) - equiPeriod],mean),collapse = "\n",sep = "\t")
