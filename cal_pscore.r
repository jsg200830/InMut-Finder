
args <- commandArgs(trailingOnly = TRUE)
filename=args[1]
data <- read.table(filename, header = FALSE)
cnts=data[,3]
mean(cnts)
scr=ppois(cnts, lambda=mean(cnts),  lower=FALSE)
data[,5]=scr
data=data[order(data[,3], decreasing = TRUE),]
names(data)=c("location","gap","counts","gene","P-value")
filename=paste(filename,"_score.csv", sep = "")
write.csv(data, file = filename,row.names=FALSE)
