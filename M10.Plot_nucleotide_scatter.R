#http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
library("ggplot2")
library("reshape2")
library("scales")
library("gtable")
library("plyr")
args <- commandArgs(TRUE)
#setwd("C:/_LocalData/Programming/R/ggplot/GFP")
#Read in data
filename<-args[1]
Rep1 <- read.delim(filename,header=T)

if (sum(Rep1$nc) < 1){
  print("File  is empty, skipping. ")
}
if (sum(Rep1$nc) > 1){
  #Delete consensus NT
  Rep1 <- within(Rep1, X.T[WT_NT=="T"] <- NA)
  Rep1 <- within(Rep1, X.A[WT_NT=="A"] <- NA)
  Rep1 <- within(Rep1, X.G[WT_NT=="G"] <- NA)
  Rep1 <- within(Rep1, X.C[WT_NT=="C"] <- NA)
  #Repair HR core
  if (grepl("stop", filename) | grepl("Stop", filename)){
    Rep1$X.T[48] <- NA
    Rep1$X.T[49] <- NA
    Rep1$X.A[50] <- NA
    Rep1$X.T[51] <- NA
    Rep1$X.T[51] <- NA
    Rep1$X.A[52] <- NA
    Rep1$X.A[53] <- NA
    Rep1$X.G[53] <- NA
  }
  if (grepl("sil", filename) | grepl("Sil", filename) | grepl("INTERPOP", filename)){
    Rep1$X.G[53] <- NA
  }
  
  #Extract Frequencies
  Freqs1<-as.data.frame(as.vector(cbind(Rep1[,"X.A"],Rep1[,"X.C"],Rep1[,"X.T"],Rep1[,"X.G"])))
  
  #Add positional information and merge data
  Data1 <- cbind(rep(seq(281,378),4),
                 Freqs1,c(rep("A",98),rep("C",98),rep("T",98),rep("G",98)),
                 rep(c(rep("L-Flank",35),rep("L-Var",12),rep("Core",6),rep("R-Var",12),rep("R-Flank",33)),4),
                 rep(1,392))
  
  colnames(Data1)<-c("Pos","Freq","NT","Region","Rep")
  CumulData <- as.data.frame(rbind(Data1))
  for (i in 1:ncol(CumulData)) CumulData[,1] <- as.numeric(as.character(CumulData[,1]))
  for (i in 1:ncol(CumulData)) CumulData[,2] <- as.numeric(as.character(CumulData[,2]))
  CumulData$Rep <- as.integer(CumulData$Rep)
  
  # Plot scatter plot
  dataAverageA<-cbind(CumulData[1:95,2]+CumulData[96:190,2]+CumulData[191:285,2])/3
  dataAverageA[is.na(dataAverageA)] <- 0
  dataAverageC<-cbind(CumulData[286:380,2]+CumulData[381:475,2]+CumulData[476:570,2])/3
  dataAverageC[is.na(dataAverageC)] <- 0
  dataAverageG<-cbind(CumulData[571:665,2]+CumulData[666:760,2]+CumulData[761:855,2])/3
  dataAverageG[is.na(dataAverageG)] <- 0
  dataAverageT<-cbind(CumulData[856:950,2]+CumulData[951:1045,2]+CumulData[1046:1140,2])/3
  dataAverageT[is.na(dataAverageT)] <- 0
  dataAverageCum<-data.frame(Pos=rep(seq(from = 281, to = 375), 4),
                             Region=c(rep('L-Flank',35),rep('L-Var',12),rep('Core',6),rep('R-Var',12),rep('R-Flank',30)),
                             NT=c(rep('A',95),rep('C',95),rep('G',95),rep('T',95)),
                             Freq=cbind(dataAverageA+dataAverageC+dataAverageG+dataAverageT)/3)
  dataCum <- data.frame(Pos=rep(seq(from = 281, to = 375), 4),
                        Region=c(rep('L-Flank',35),rep('L-Var',12),rep('Core',6),rep('R-Var',12),rep('R-Flank',30)),
                        NT=c(rep('A',95),rep('C',95),rep('G',95),rep('T',95)), 
                        Freq=as.vector(cbind(dataAverageA,dataAverageC,dataAverageG,dataAverageT)))
  dataCum <- dataCum[!(apply(dataCum, 1, function(y) any(y == 0))),]
  dataCum[dataCum==0] <- NA
  CumulData[CumulData==100]<-0
  
  outputfile<-paste(filename,"_stackedbar_plot.pdf", sep ="")
  p <- ggplot(data=CumulData, aes(x=Pos, y=Freq, group=NT)) + 
    labs(title= outputfile ) +
    ylab("%") +
    theme(plot.title = element_text(hjust=0.5)) +
    geom_point(aes(color=NT)) + scale_color_brewer(palette="Spectral") +  
    geom_hline(aes(yintercept=2), colour="#999999") + 
    geom_vline(aes(xintercept=315.5), colour="#999999", linetype="dashed") + 
    geom_vline(aes(xintercept=327.5), colour="#999999", linetype="dashed") + 
    geom_vline(aes(xintercept=333.5), colour="#999999", linetype="dashed") + 
    geom_vline(aes(xintercept=345.5), colour="#999999", linetype="dashed")
  ggsave(outputfile)
  graphics.off()  
  
}