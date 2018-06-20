library("ggplot2")
library("reshape2")
library("scales")
library("gtable")
library("plyr")

setwd("~/Desktop")
setwd("Z:/Martijn/GFP/2018January/post_analysis")
Inputfile <- na.omit(read.csv(file="2018January_allmutations.csv",header=T,as.is = TRUE))
Inputfile <- Inputfile[,c(1:4,6,7,5)]
Input.mut <- Inputfile[rep(row.names(Inputfile),Inputfile$Count),1:6]



subzet <- subset(Input.mut,ref != alt) #remove silent mutations, 40,190/55,724 are missense
subzet <- subzet[!(subzet$Sample %in% c('GFP_Aug_g8_VarSil12_Neg','GFP_Aug_g8_VarSil12_Int','GFP_Aug_g8_VarSil12_Pos')),] #remove g8
subzet$Sample <- factor(subzet$Sample, levels = c("Stop_Oligo","GFP_March_VarStop12_Neg","GFP_Nov_VarStop12_Neg_br1","GFP_Nov_VarStop12_Neg_br2",
                                                  "GFP_March_VarSil12_Neg","GFP_Nov_VarSil12_Neg_br1","GFP_Nov_VarSil12_Neg_br2",
                                                  "GFP_GFPhetINTERPOP_Mixed","GFP_March_VarSil12_Int","GFP_Nov_VarSil12_Int_br1","GFP_Nov_VarSil12_Int_br2",
                                                  "GFP_March_VarSil12_Pos", "GFP_Nov_VarSil12_Pos_br1","GFP_Nov_VarSil12_Pos_br2")) #reorder
unique(unlist(subzet$Sample, use.names = F))
SampleSizes<-count(subzet,'Sample')

outputfile<-"GFP_mutations_violin.pdf"
p <- ggplot(subzet, aes(Sample,ddG))
p <- p + geom_violin(aes(fill = Sample))
#p <- p + geom_jitter(height = 0, width = 0.1, size = 1,shape=4,alpha=0.5)
p <- p + stat_summary(fun.y=mean, geom="point", shape=23, fill="red", size=2)
p <- p + labs(x="Sample", y="ddG")
#for (i in 1:length(SampleSizes[,1])){
#  p <- p + annotate("text", x=i, y=9.5, size = 3, label = paste('N= ',SampleSizes[i,2],sep=''))
#}
p <- p + theme(axis.text.x = element_text(angle = 90, size=10, hjust=0.95,vjust=0))
p <- p + theme(legend.position="none")
p <- p + geom_hline(aes(yintercept=0.8), colour="cyan", linetype="dashed")
p <- p + geom_hline(aes(yintercept=1.6), colour="magenta", linetype="dashed")
p
ggsave(outputfile)

###Summed
Subzetsum<-subzet
for (i in c("GFP_March_VarStop12_Neg","GFP_Nov_VarStop12_Neg_br1","GFP_Nov_VarStop12_Neg_br2")){ 
  Subzetsum<-data.frame(lapply(Subzetsum,function(x){
    gsub(i,'Stop_Neg',x)
  }))
}
for (i in c("GFP_March_VarSil12_Neg","GFP_Nov_VarSil12_Neg_br1","GFP_Nov_VarSil12_Neg_br2")){ 
  Subzetsum<-data.frame(lapply(Subzetsum,function(x){
    gsub(i,'Silent_Neg',x)
  }))
}
for (i in c("GFP_GFPhetINTERPOP_Mixed","GFP_March_VarSil12_Int","GFP_Nov_VarSil12_Int_br1","GFP_Nov_VarSil12_Int_br2")){
  Subzetsum<-data.frame(lapply(Subzetsum,function(x){
    gsub(i,'Silent_Int',x)
  }))
}
for (i in c("GFP_March_VarSil12_Pos", "GFP_Nov_VarSil12_Pos_br1","GFP_Nov_VarSil12_Pos_br2")){ 
  Subzetsum<-data.frame(lapply(Subzetsum,function(x){
    gsub(i,'Silent_Pos',x)
  }))
}

ddGlist<-read.csv("GFP_ddG_values.csv")
ddGlist<-as.data.frame(cbind(as.character(ddGlist$Substitution),as.character(ddGlist$Substitution),
                              as.character(ddGlist$Substitution),as.character(ddGlist$Substitution),
                              as.character(ddGlist$Substitution),ddGlist$Calculated...G))
colnames(ddGlist)<-c("Mutation","ref","pos","alt","Sample","ddG")
ddGlist$Sample<-"Possible_AAsubs"
Subzetsum<-as.data.frame(rbind(Subzetsum,ddGlist))

Subzetsum$Sample <- factor(Subzetsum$Sample, levels = c('Possible_AAsubs','Stop_Oligo','Stop_Neg','Silent_Neg','Silent_Int','Silent_Pos'))
Subzetsum$Sample[is.na(Subzetsum$Sample)]<-"Stop_Oligo"
Subzetsum$ddG <- as.numeric(as.character(Subzetsum$ddG))
SampleSizes<-count(Subzetsum,'Sample') #from plyr package

outputfile<-"GFP_mutations_violin_sum.pdf"
p <- ggplot(Subzetsum, aes(Sample,ddG))
p <- p + geom_violin(aes(fill = Sample))
#p <- p + geom_jitter(height = 0, width = 0.1, size = 1,shape=4,alpha=0.5)
p <- p + stat_summary(fun.y=mean, geom="point", shape=23, fill="red", size=2)
p <- p + labs(x="Sample", y="ddG")
for (i in 1:length(SampleSizes[,1])){
  p <- p + annotate("text", x=i, y=9.5, size = 3, label = paste('N= ',SampleSizes[i,2],sep=''))
}
p <- p + theme(axis.text.x = element_text(angle = 90, size=10, hjust=0.95,vjust=0))
p <- p + theme(legend.position="none")
p <- p + geom_hline(aes(yintercept=0.8), colour="cyan", linetype="dashed")
p <- p + geom_hline(aes(yintercept=1.6), colour="magenta", linetype="dashed")
p
ggsave(outputfile)


#######
# install.packages('venneuler')
#library(venneuler)

SilInt <- as.character(Subzetsum[Subzetsum$Sample == "Silent_Int",]$Mutation)
SilPos <- as.character(Subzetsum[Subzetsum$Sample == "Silent_Pos",]$Mutation)
SilNeg <- as.character(Subzetsum[Subzetsum$Sample == "Silent_Neg",]$Mutation)
StopNeg <- as.character(Subzetsum[Subzetsum$Sample == "Stop_Neg",]$Mutation)
Oligo <- as.character(Subzetsum[Subzetsum$Sample == "Stop_Oligo",]$Mutation)
Lists <- list(SilInt,SilNeg,SilPos,StopNeg,Oligo)  #put the word vectors into a list to supply lapply
items <- sort(unique(unlist(Lists)))   #put in alphabetical order
MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=5)  #make a matrix of 0's
colnames(MAT) <- c("SilInt","SilNeg","SilPos","StopNeg","Oligo")
rownames(MAT) <- items
MAT<-MAT[- grep("\\*", rownames(MAT)),]

lapply(seq_along(Lists), function(i) {   #fill the matrix
  MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
})
MAT   #look at the results

library(venneuler)
v <- venneuler(MAT)
plot(v)

Lists <- list(StopNeg,Oligo)  #put the word vectors into a list to supply lapply
items <- sort(unique(unlist(Lists)))   #put in alphabetical order
MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=2)  #make a matrix of 0's
colnames(MAT) <- c("StopNeg","Oligo")
rownames(MAT) <- items
lapply(seq_along(Lists), function(i) {   #fill the matrix
  MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
})
MAT   #look at the results
MAT<-MAT[- grep("\\*", rownames(MAT)),]
library(venneuler)
v <- venneuler(MAT)
plot(v)
