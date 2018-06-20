#read in file
args <- commandArgs(TRUE)
filename <- args[1]
mutationfile<-read.delim(filename, stringsAsFactors = FALSE)

#for filename conversion
filencount<-20
if (grepl("stop", filename) | grepl("Stop", filename)){
  filencount<-21
}
#isolate read with single mutation
mutlist <- as.data.frame(mutationfile$AA_MUT_HR, mutationfile$ID, stringsAsFactors = FALSE)
mutlist <- cbind(rownames(mutlist), mutlist)
rownames(mutlist) <- NULL
colnames(mutlist) <- c("ID","mutation")

#separate fields
mutlist <- mutlist[nchar(mutlist$mutation) <= 5, ] 
mutlist2 <- stringr::str_split_fixed(mutlist$mutation,',',6)
#melt into single row and remove emtpy cells
avector <- c(mutlist2[,1],mutlist2[,2],mutlist2[,3],mutlist2[,4],mutlist2[,5],mutlist2[,6])
avector <- avector[avector != ""]

#count unique mutations
mutcount <- as.data.frame(tapply(avector, avector, length))
mutcount <- as.data.frame(mutcount[rownames(mutcount) != "AA_MUT_HR",])
colnames(mutcount) <- c("frequency")
mutcount$mutation<-rownames(mutcount)
mutcount$mutation[mutcount$mutation == '82'] <- 'H82.'

#format rows
mutcount$ref <- substring(gsub("[[:digit:]]","",mutcount$mutation),1,1)
mutcount$pos <- as.numeric(as.character(gsub("[^0-9\\]", "", mutcount$mutation)))
mutcount$alt <- substring(gsub("[[:digit:]]","",mutcount$mutation),2,2)
mutcount<- as.data.frame(cbind(mutcount$mutation,mutcount$ref,mutcount$pos,mutcount$alt,mutcount$frequency))
colnames(mutcount)<- c("mutation","ref","pos","alt","frequency")

#load and deal with ddG
ddGfile<-read.csv('GFP_ddG_values.csv', stringsAsFactors = FALSE)
mutcount$ddG[1] <- NA
for (i in 1:nrow(mutcount)){
  for (j in c(1:160)){
    if (mutcount$mutation[i] == ddGfile$Substitution[j]){
      mutcount$ddG[i] <- ddGfile$Calculated...G[j]
    }
  }
}
mutcount$sample<-substr(filename,1,filencount)

#isolate only variable positions
mutcount$pos <- as.numeric(as.character(mutcount$pos))
varcount <- mutcount[mutcount$pos > 105 & mutcount$pos < 116,]
varcount$sample<-substr(filename,1,21)
#write file
write.csv(mutcount, paste(substr(filename,1,filencount),".mutfreq.all.CSV",sep=""))
write.csv(varcount, paste(substr(filename,1,filencount),".mutfreq.var.CSV",sep=""))