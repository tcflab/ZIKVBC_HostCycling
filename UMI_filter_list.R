#### Load packages ####
list.of.packages <- c("dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(dplyr)

#### Enter the sample name and output directory for the sample
args<-commandArgs(trailingOnly = TRUE)
outputdir<-args[1]
samplename<-args[2]
inputRNA<-args[4]
tsvfile1<-paste("UMI_groups_",samplename,sep="")
tsvfile2<-paste(tsvfile1,".tsv",sep="")
tsvpath<-paste(outputdir,tsvfile2,sep="/")

#### Read tsv file ####
tsv<-read.delim(tsvpath, header=TRUE, sep="")

#### Subset to UMIs observed at least 3 times ####
tsv<-tsv[c(7,8)]

tsv_3filter<-subset(tsv,tsv$final_umi_count>2)

#### Filter rows to list only one UMI per group ####
oneumi<-distinct(tsv_3filter)

#### Export list of UMIs as a .lst file ####
list<-oneumi[c(1)]
setwd(outputdir)
write.table(list, file="passumifilter.lst", row.names = FALSE, col.names = FALSE, quote = FALSE)

#### Calculate UMI conversion rate ####
umi_consensus_count<-nrow(oneumi)

conversion_rate<-(as.numeric(umi_consensus_count)/(10^(as.numeric(inputRNA))))*100

conv_rate<-c(conversion_rate)
conv_rate<-t(conv_rate)
colnames(conv_rate)<-c("converation_rate")
conv_rate<-as.data.frame(conv_rate)
write.csv(conversion_rate, file="conversionrate.csv", row.names=FALSE)
