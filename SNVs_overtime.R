## load required packages
list.of.packages <- c("ggplot2","dplyr","splitstackshape","RcppRoll","magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(splitstackshape))
suppressPackageStartupMessages(library(RcppRoll))
suppressPackageStartupMessages(library(magrittr))


######  Mosquito Lineage By Passage     ######################################################

p1_L1<-read.table(file="/PATH/TO/rep1and2_refalignMP1_1.vcf", sep="\t", header=FALSE)
p1_L1<-cSplit(p1_L1,"V8",sep=";", type.convert=FALSE)
p1_L1<-cSplit(p1_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L1$SNVID<-paste(p1_L1$V2,p1_L1$V4,sep=";")
p1_L1$SNVID<-paste(p1_L1$SNVID,p1_L1$V5,sep=">")
p1_L1<-p1_L1[,c(12,11)]
colnames(p1_L1)<-c("SNVID","freq")
p1_L1<-subset(p1_L1,p1_L1$freq>=0.003)

p3_L1<-read.table(file="/PATH/TO/rep1and2_refalignMP3_1.vcf", sep="\t", header=FALSE)
p3_L1<-cSplit(p3_L1,"V8",sep=";", type.convert=FALSE)
p3_L1<-cSplit(p3_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L1$SNVID<-paste(p3_L1$V2,p3_L1$V4,sep=";")
p3_L1$SNVID<-paste(p3_L1$SNVID,p3_L1$V5,sep=">")
p3_L1<-p3_L1[,c(12,11)]
colnames(p3_L1)<-c("SNVID","freq")
p3_L1<-subset(p3_L1,p3_L1$freq>=0.003)

p13_L1<-merge(p1_L1, p3_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L1)<-c("SNVID","freq1","freq3")

p4_L1<-read.table(file="/PATH/TO/rep1and2_refalignMP4_1.vcf", sep="\t", header=FALSE)
p4_L1<-cSplit(p4_L1,"V8",sep=";", type.convert=FALSE)
p4_L1<-cSplit(p4_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L1$SNVID<-paste(p4_L1$V2,p4_L1$V4,sep=";")
p4_L1$SNVID<-paste(p4_L1$SNVID,p4_L1$V5,sep=">")
p4_L1<-p4_L1[,c(12,11)]
colnames(p4_L1)<-c("SNVID","freq")
p4_L1<-subset(p4_L1,p4_L1$freq>=0.003)

p14_L1<-merge(p13_L1, p4_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L1)<-c("SNVID","freq1","freq3","freq4")

p10_L1<-read.table(file="/PATH/TO/rep1and2_refalignMP10_1.vcf", sep="\t", header=FALSE)
p10_L1<-cSplit(p10_L1,"V8",sep=";", type.convert=FALSE)
p10_L1<-cSplit(p10_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L1$SNVID<-paste(p10_L1$V2,p10_L1$V4,sep=";")
p10_L1$SNVID<-paste(p10_L1$SNVID,p10_L1$V5,sep=">")
p10_L1<-p10_L1[,c(12,11)]
colnames(p10_L1)<-c("SNVID","freq")
p10_L1<-subset(p10_L1,p10_L1$freq>=0.003)

p110_L1<-merge(p14_L1, p10_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L1)<-c("SNVID","freq1","freq2","freq3","freq4")

p11_L1<-read.table(file="/PATH/TO/rep1and2_refalignMP10_A.vcf", sep="\t", header=FALSE)
p11_L1<-cSplit(p11_L1,"V8",sep=";", type.convert=FALSE)
p11_L1<-cSplit(p11_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L1$SNVID<-paste(p11_L1$V2,p11_L1$V4,sep=";")
p11_L1$SNVID<-paste(p11_L1$SNVID,p11_L1$V5,sep=">")
p11_L1<-p11_L1[,c(12,11)]
colnames(p11_L1)<-c("SNVID","freq")
p11_L1<-subset(p11_L1,p11_L1$freq>=0.003)

p111_L1<-merge(p110_L1, p11_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L1)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p111_L1[is.na(p111_L1)]<-0


p1_L2<-read.table(file="/PATH/TO/rep1and2_refalignMP1_2.vcf", sep="\t", header=FALSE)
p1_L2<-cSplit(p1_L2,"V8",sep=";", type.convert=FALSE)
p1_L2<-cSplit(p1_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L2$SNVID<-paste(p1_L2$V2,p1_L2$V4,sep=";")
p1_L2$SNVID<-paste(p1_L2$SNVID,p1_L2$V5,sep=">")
p1_L2<-p1_L2[,c(12,11)]
colnames(p1_L2)<-c("SNVID","freq")
p1_L2<-subset(p1_L2,p1_L2$freq>=0.003)

p3_L2<-read.table(file="/PATH/TO/rep1and2_refalignMP3_2.vcf", sep="\t", header=FALSE)
p3_L2<-cSplit(p3_L2,"V8",sep=";", type.convert=FALSE)
p3_L2<-cSplit(p3_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L2$SNVID<-paste(p3_L2$V2,p3_L2$V4,sep=";")
p3_L2$SNVID<-paste(p3_L2$SNVID,p3_L2$V5,sep=">")
p3_L2<-p3_L2[,c(12,11)]
colnames(p3_L2)<-c("SNVID","freq")
p3_L2<-subset(p3_L2,p3_L2$freq>=0.003)

p13_L2<-merge(p1_L2, p3_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L2)<-c("SNVID","freq1","freq3")

p4_L2<-read.table(file="/PATH/TO/rep1and2_refalignMP4_2.vcf", sep="\t", header=FALSE)
p4_L2<-cSplit(p4_L2,"V8",sep=";", type.convert=FALSE)
p4_L2<-cSplit(p4_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L2$SNVID<-paste(p4_L2$V2,p4_L2$V4,sep=";")
p4_L2$SNVID<-paste(p4_L2$SNVID,p4_L2$V5,sep=">")
p4_L2<-p4_L2[,c(12,11)]
colnames(p4_L2)<-c("SNVID","freq")
p4_L2<-subset(p4_L2,p4_L2$freq>=0.003)

p14_L2<-merge(p13_L2, p4_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L2)<-c("SNVID","freq1","freq3","freq4")

p10_L2<-read.table(file="/PATH/TO/rep1and2_refalignMP10_2.vcf", sep="\t", header=FALSE)
p10_L2<-cSplit(p10_L2,"V8",sep=";", type.convert=FALSE)
p10_L2<-cSplit(p10_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L2$SNVID<-paste(p10_L2$V2,p10_L2$V4,sep=";")
p10_L2$SNVID<-paste(p10_L2$SNVID,p10_L2$V5,sep=">")
p10_L2<-p10_L2[,c(12,11)]
colnames(p10_L2)<-c("SNVID","freq")
p10_L2<-subset(p10_L2,p10_L2$freq>=0.003)

p110_L2<-merge(p14_L2, p10_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L2)<-c("SNVID","freq1","freq2","freq3","freq4")

p11_L2<-read.table(file="/PATH/TO/rep1and2_refalignMP10_B.vcf", sep="\t", header=FALSE)
p11_L2<-cSplit(p11_L2,"V8",sep=";", type.convert=FALSE)
p11_L2<-cSplit(p11_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L2$SNVID<-paste(p11_L2$V2,p11_L2$V4,sep=";")
p11_L2$SNVID<-paste(p11_L2$SNVID,p11_L2$V5,sep=">")
p11_L2<-p11_L2[,c(12,11)]
colnames(p11_L2)<-c("SNVID","freq")
p11_L2<-subset(p11_L2,p11_L2$freq>=0.003)

p111_L2<-merge(p110_L2, p11_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L2)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p111_L2[is.na(p111_L2)]<-0


p1_L3<-read.table(file="/PATH/TO/rep1and2_refalignMP1_3.vcf", sep="\t", header=FALSE)
p1_L3<-cSplit(p1_L3,"V8",sep=";", type.convert=FALSE)
p1_L3<-cSplit(p1_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L3$SNVID<-paste(p1_L3$V2,p1_L3$V4,sep=";")
p1_L3$SNVID<-paste(p1_L3$SNVID,p1_L3$V5,sep=">")
p1_L3<-p1_L3[,c(12,11)]
colnames(p1_L3)<-c("SNVID","freq")
p1_L3<-subset(p1_L3,p1_L3$freq>=0.003)

p3_L3<-read.table(file="/PATH/TO/rep1and2_refalignMP3_3.vcf", sep="\t", header=FALSE)
p3_L3<-cSplit(p3_L3,"V8",sep=";", type.convert=FALSE)
p3_L3<-cSplit(p3_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L3$SNVID<-paste(p3_L3$V2,p3_L3$V4,sep=";")
p3_L3$SNVID<-paste(p3_L3$SNVID,p3_L3$V5,sep=">")
p3_L3<-p3_L3[,c(12,11)]
colnames(p3_L3)<-c("SNVID","freq")
p3_L3<-subset(p3_L3,p3_L3$freq>=0.003)

p13_L3<-merge(p1_L3, p3_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L3)<-c("SNVID","freq1","freq3")

p4_L3<-read.table(file="/PATH/TO/rep1and2_refalignMP4_3.vcf", sep="\t", header=FALSE)
p4_L3<-cSplit(p4_L3,"V8",sep=";", type.convert=FALSE)
p4_L3<-cSplit(p4_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L3$SNVID<-paste(p4_L3$V2,p4_L3$V4,sep=";")
p4_L3$SNVID<-paste(p4_L3$SNVID,p4_L3$V5,sep=">")
p4_L3<-p4_L3[,c(12,11)]
colnames(p4_L3)<-c("SNVID","freq")
p4_L3<-subset(p4_L3,p4_L3$freq>=0.003)

p14_L3<-merge(p13_L3, p4_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L3)<-c("SNVID","freq1","freq3","freq4")

p10_L3<-read.table(file="/PATH/TO/rep1and2_refalignMP10_3.vcf", sep="\t", header=FALSE)
p10_L3<-cSplit(p10_L3,"V8",sep=";", type.convert=FALSE)
p10_L3<-cSplit(p10_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L3$SNVID<-paste(p10_L3$V2,p10_L3$V4,sep=";")
p10_L3$SNVID<-paste(p10_L3$SNVID,p10_L3$V5,sep=">")
p10_L3<-p10_L3[,c(12,11)]
colnames(p10_L3)<-c("SNVID","freq")
p10_L3<-subset(p10_L3,p10_L3$freq>=0.003)

p110_L3<-merge(p14_L3, p10_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L3)<-c("SNVID","freq1","freq3","freq4","freq10")

p11_L3<-read.table(file="/PATH/TO/rep1and2_refalignMP10_C.vcf", sep="\t", header=FALSE)
p11_L3<-cSplit(p11_L3,"V8",sep=";", type.convert=FALSE)
p11_L3<-cSplit(p11_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L3$SNVID<-paste(p11_L3$V2,p11_L3$V4,sep=";")
p11_L3$SNVID<-paste(p11_L3$SNVID,p11_L3$V5,sep=">")
p11_L3<-p11_L3[,c(12,11)]
colnames(p11_L3)<-c("SNVID","freq")
p11_L3<-subset(p11_L3,p11_L3$freq>=0.003)

p111_L3<-merge(p110_L3, p11_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L3)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p111_L3[is.na(p111_L3)]<-0

p111_L3<-merge(p110_L3, p11_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L3)<-c("SNVID","freq1","freq3","freq4","freq10","freq11")

p111_L3[is.na(p111_L3)]<-0

p1_L4<-read.table(file="/PATH/TO/rep1and2_refalignMP1_4.vcf", sep="\t", header=FALSE)
p1_L4<-cSplit(p1_L4,"V8",sep=";", type.convert=FALSE)
p1_L4<-cSplit(p1_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L4$SNVID<-paste(p1_L4$V2,p1_L4$V4,sep=";")
p1_L4$SNVID<-paste(p1_L4$SNVID,p1_L4$V5,sep=">")
p1_L4<-p1_L4[,c(12,11)]
colnames(p1_L4)<-c("SNVID","freq")
p1_L4<-subset(p1_L4,p1_L4$freq>=0.003)

p3_L4<-read.table(file="/PATH/TO/rep1and2_refalignMP3_4.vcf", sep="\t", header=FALSE)
p3_L4<-cSplit(p3_L4,"V8",sep=";", type.convert=FALSE)
p3_L4<-cSplit(p3_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L4$SNVID<-paste(p3_L4$V2,p3_L4$V4,sep=";")
p3_L4$SNVID<-paste(p3_L4$SNVID,p3_L4$V5,sep=">")
p3_L4<-p3_L4[,c(12,11)]
colnames(p3_L4)<-c("SNVID","freq")
p3_L4<-subset(p3_L4,p3_L4$freq>=0.003)

p13_L4<-merge(p1_L4, p3_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L4)<-c("SNVID","freq1","freq3")

p4_L4<-read.table(file="/PATH/TO/rep1and2_refalignMP4_4.vcf", sep="\t", header=FALSE)
p4_L4<-cSplit(p4_L4,"V8",sep=";", type.convert=FALSE)
p4_L4<-cSplit(p4_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L4$SNVID<-paste(p4_L4$V2,p4_L4$V4,sep=";")
p4_L4$SNVID<-paste(p4_L4$SNVID,p4_L4$V5,sep=">")
p4_L4<-p4_L4[,c(12,11)]
colnames(p4_L4)<-c("SNVID","freq")
p4_L4<-subset(p4_L4,p4_L4$freq>=0.003)

p14_L4<-merge(p13_L4, p4_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L4)<-c("SNVID","freq1","freq3","freq4")

p10_L4<-read.table(file="/PATH/TO/rep1and2_refalignMP10_4.vcf", sep="\t", header=FALSE)
p10_L4<-cSplit(p10_L4,"V8",sep=";", type.convert=FALSE)
p10_L4<-cSplit(p10_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L4$SNVID<-paste(p10_L4$V2,p10_L4$V4,sep=";")
p10_L4$SNVID<-paste(p10_L4$SNVID,p10_L4$V5,sep=">")
p10_L4<-p10_L4[,c(12,11)]
colnames(p10_L4)<-c("SNVID","freq")
p10_L4<-subset(p10_L4,p10_L4$freq>=0.003)

p110_L4<-merge(p14_L4, p10_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L4)<-c("SNVID","freq1","freq2","freq3","freq4")

p11_L4<-read.table(file="/PATH/TO/rep1and2_refalignMP10_D.vcf", sep="\t", header=FALSE)
p11_L4<-cSplit(p11_L4,"V8",sep=";", type.convert=FALSE)
p11_L4<-cSplit(p11_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L4$SNVID<-paste(p11_L4$V2,p11_L4$V4,sep=";")
p11_L4$SNVID<-paste(p11_L4$SNVID,p11_L4$V5,sep=">")
p11_L4<-p11_L4[,c(12,11)]
colnames(p11_L4)<-c("SNVID","freq")
p11_L4<-subset(p11_L4,p11_L4$freq>=0.003)

p111_L4<-merge(p110_L4, p11_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L4)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p111_L4[is.na(p111_L4)]<-0


p1_L5<-read.table(file="/PATH/TO/rep1and2_refalignMP1_5.vcf", sep="\t", header=FALSE)
p1_L5<-cSplit(p1_L5,"V8",sep=";", type.convert=FALSE)
p1_L5<-cSplit(p1_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L5$SNVID<-paste(p1_L5$V2,p1_L5$V4,sep=";")
p1_L5$SNVID<-paste(p1_L5$SNVID,p1_L5$V5,sep=">")
p1_L5<-p1_L5[,c(12,11)]
colnames(p1_L5)<-c("SNVID","freq")
p1_L5<-subset(p1_L5,p1_L5$freq>=0.003)

p3_L5<-read.table(file="/PATH/TO/rep1and2_refalignMP3_5.vcf", sep="\t", header=FALSE)
p3_L5<-cSplit(p3_L5,"V8",sep=";", type.convert=FALSE)
p3_L5<-cSplit(p3_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L5$SNVID<-paste(p3_L5$V2,p3_L5$V4,sep=";")
p3_L5$SNVID<-paste(p3_L5$SNVID,p3_L5$V5,sep=">")
p3_L5<-p3_L5[,c(12,11)]
colnames(p3_L5)<-c("SNVID","freq")
p3_L5<-subset(p3_L5,p3_L5$freq>=0.003)

p13_L5<-merge(p1_L5, p3_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L5)<-c("SNVID","freq1","freq3")

p4_L5<-read.table(file="/PATH/TO/rep1and2_refalignMP4_5.vcf", sep="\t", header=FALSE)
p4_L5<-cSplit(p4_L5,"V8",sep=";", type.convert=FALSE)
p4_L5<-cSplit(p4_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L5$SNVID<-paste(p4_L5$V2,p4_L5$V4,sep=";")
p4_L5$SNVID<-paste(p4_L5$SNVID,p4_L5$V5,sep=">")
p4_L5<-p4_L5[,c(12,11)]
colnames(p4_L5)<-c("SNVID","freq")
p4_L5<-subset(p4_L5,p4_L5$freq>=0.003)

p14_L5<-merge(p13_L5, p4_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L5)<-c("SNVID","freq1","freq3","freq4")

p10_L5<-read.table(file="/PATH/TO/rep1and2_refalignMP10_5.vcf", sep="\t", header=FALSE)
p10_L5<-cSplit(p10_L5,"V8",sep=";", type.convert=FALSE)
p10_L5<-cSplit(p10_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L5$SNVID<-paste(p10_L5$V2,p10_L5$V4,sep=";")
p10_L5$SNVID<-paste(p10_L5$SNVID,p10_L5$V5,sep=">")
p10_L5<-p10_L5[,c(12,11)]
colnames(p10_L5)<-c("SNVID","freq")
p10_L5<-subset(p10_L5,p10_L5$freq>=0.003)

p110_L5<-merge(p14_L5, p10_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L5)<-c("SNVID","freq1","freq2","freq3","freq4")

p11_L5<-read.table(file="/PATH/TO/rep1and2_refalignMP10_E.vcf", sep="\t", header=FALSE)
p11_L5<-cSplit(p11_L5,"V8",sep=";", type.convert=FALSE)
p11_L5<-cSplit(p11_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L5$SNVID<-paste(p11_L5$V2,p11_L5$V4,sep=";")
p11_L5$SNVID<-paste(p11_L5$SNVID,p11_L5$V5,sep=">")
p11_L5<-p11_L5[,c(12,11)]
colnames(p11_L5)<-c("SNVID","freq")
p11_L5<-subset(p11_L5,p11_L5$freq>=0.003)

p111_L5<-merge(p110_L5, p11_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L5)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p111_L5[is.na(p111_L5)]<-0


L1andL2<-merge(p111_L1, p111_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1andL2)<-c("SNVID","p1_L1","p3_L1","p4_L1","p10_L1","p11_L1","p1_L2","p3_L2","p4_L2","p10_L2","p11_L2")
L1thruL3<-merge(L1andL2,p111_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1thruL3)<-c("SNVID","p1_L1","p3_L1","p4_L1","p10_L1","p11_L1","p1_L2","p3_L2","p4_L2","p10_L2","p11_L2","p1_L3","p3_L3","p4_L3","p10_L3","p11_L3")
L1thruL4<-merge(L1thruL3,p111_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1thruL4)<-c("SNVID","p1_L1","p3_L1","p4_L1","p10_L1","p11_L1","p1_L2","p3_L2","p4_L2","p10_L2","p11_L2","p1_L3","p3_L3","p4_L3","p10_L3","p11_L3","p1_L4","p3_L4","p4_L4","p10_L4","p11_L4")
Mosquito_L1thruL5<-merge(L1thruL4,p111_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(Mosquito_L1thruL5)<-c("SNVID","p1_L1","p3_L1","p4_L1","p10_L1","p11_L1","p1_L2","p3_L2","p4_L2","p10_L2","p11_L2","p1_L3","p3_L3","p4_L3","p10_L3","p11_L3","p1_L4","p3_L4","p4_L4","p10_L4","p11_L4","p1_L5","p3_L5","p4_L5","p10_L5","p11_L5")

Mosquito_L1thruL5[is.na(Mosquito_L1thruL5)]<-0

Lineage1_p1<-Mosquito_L1thruL5[,c(1,2)]
Lineage1_p1$passage<-1
colnames(Lineage1_p1)<-c("SNVID","freq","passage")
Lineage1_p3<-Mosquito_L1thruL5[,c(1,3)]
Lineage1_p3$passage<-3
colnames(Lineage1_p3)<-c("SNVID","freq","passage")
Lineage1_p4<-Mosquito_L1thruL5[,c(1,4)]
Lineage1_p4$passage<-4
colnames(Lineage1_p4)<-c("SNVID","freq","passage")
Lineage1_p10<-Mosquito_L1thruL5[,c(1,5)]
Lineage1_p10$passage<-10
colnames(Lineage1_p10)<-c("SNVID","freq","passage")
Lineage1_p11<-Mosquito_L1thruL5[,c(1,6)]
Lineage1_p11$passage<-11
colnames(Lineage1_p11)<-c("SNVID","freq","passage")

Lineage1<-rbind(Lineage1_p1,Lineage1_p3)
Lineage1<-rbind(Lineage1,Lineage1_p4)
Lineage1<-rbind(Lineage1,Lineage1_p10)
Mosquito_Lineage1<-rbind(Lineage1,Lineage1_p11)

Mosquito_Lineage1$freq<-as.numeric(Mosquito_Lineage1$freq)
Mosquito_Lineage1$passage<-as.numeric(Mosquito_Lineage1$passage)

Lineage2_p1<-Mosquito_L1thruL5[,c(1,7)]
Lineage2_p1$passage<-1
colnames(Lineage2_p1)<-c("SNVID","freq","passage")
Lineage2_p3<-Mosquito_L1thruL5[,c(1,8)]
Lineage2_p3$passage<-3
colnames(Lineage2_p3)<-c("SNVID","freq","passage")
Lineage2_p4<-Mosquito_L1thruL5[,c(1,9)]
Lineage2_p4$passage<-4
colnames(Lineage2_p4)<-c("SNVID","freq","passage")
Lineage2_p10<-Mosquito_L1thruL5[,c(1,10)]
Lineage2_p10$passage<-10
colnames(Lineage2_p10)<-c("SNVID","freq","passage")
Lineage2_p11<-Mosquito_L1thruL5[,c(1,11)]
Lineage2_p11$passage<-11
colnames(Lineage2_p11)<-c("SNVID","freq","passage")

Lineage2<-rbind(Lineage2_p1,Lineage2_p3)
Lineage2<-rbind(Lineage2,Lineage2_p4)
Lineage2<-rbind(Lineage2,Lineage2_p10)
Mosquito_Lineage2<-rbind(Lineage2,Lineage2_p11)

Mosquito_Lineage2$freq<-as.numeric(Mosquito_Lineage2$freq)
Mosquito_Lineage2$passage<-as.numeric(Mosquito_Lineage2$passage)

Lineage3_p1<-Mosquito_L1thruL5[,c(1,12)]
Lineage3_p1$passage<-1
colnames(Lineage3_p1)<-c("SNVID","freq","passage")
Lineage3_p3<-Mosquito_L1thruL5[,c(1,13)]
Lineage3_p3$passage<-3
colnames(Lineage3_p3)<-c("SNVID","freq","passage")
Lineage3_p4<-Mosquito_L1thruL5[,c(1,14)]
Lineage3_p4$passage<-4
colnames(Lineage3_p4)<-c("SNVID","freq","passage")
Lineage3_p10<-Mosquito_L1thruL5[,c(1,15)]
Lineage3_p10$passage<-10
colnames(Lineage3_p10)<-c("SNVID","freq","passage")
Lineage3_p11<-Mosquito_L1thruL5[,c(1,16)]
Lineage3_p11$passage<-11
colnames(Lineage3_p11)<-c("SNVID","freq","passage")

Lineage3<-rbind(Lineage3_p1,Lineage3_p3)
Lineage3<-rbind(Lineage3,Lineage3_p4)
Lineage3<-rbind(Lineage3,Lineage3_p10)
Mosquito_Lineage3<-rbind(Lineage3,Lineage3_p11)

Mosquito_Lineage3$freq<-as.numeric(Mosquito_Lineage3$freq)
Mosquito_Lineage3$passage<-as.numeric(Mosquito_Lineage3$passage)

Lineage4_p1<-Mosquito_L1thruL5[,c(1,17)]
Lineage4_p1$passage<-1
colnames(Lineage4_p1)<-c("SNVID","freq","passage")
Lineage4_p3<-Mosquito_L1thruL5[,c(1,18)]
Lineage4_p3$passage<-3
colnames(Lineage4_p3)<-c("SNVID","freq","passage")
Lineage4_p4<-Mosquito_L1thruL5[,c(1,19)]
Lineage4_p4$passage<-4
colnames(Lineage4_p4)<-c("SNVID","freq","passage")
Lineage4_p10<-Mosquito_L1thruL5[,c(1,20)]
Lineage4_p10$passage<-10
colnames(Lineage4_p10)<-c("SNVID","freq","passage")
Lineage4_p11<-Mosquito_L1thruL5[,c(1,21)]
Lineage4_p11$passage<-11
colnames(Lineage4_p11)<-c("SNVID","freq","passage")

Lineage4<-rbind(Lineage4_p1,Lineage4_p3)
Lineage4<-rbind(Lineage4,Lineage4_p4)
Lineage4<-rbind(Lineage4,Lineage4_p10)
Mosquito_Lineage4<-rbind(Lineage4,Lineage4_p11)

Mosquito_Lineage4$freq<-as.numeric(Mosquito_Lineage4$freq)
Mosquito_Lineage4$passage<-as.numeric(Mosquito_Lineage4$passage)

Lineage5_p1<-Mosquito_L1thruL5[,c(1,22)]
Lineage5_p1$passage<-1
colnames(Lineage5_p1)<-c("SNVID","freq","passage")
Lineage5_p3<-Mosquito_L1thruL5[,c(1,23)]
Lineage5_p3$passage<-3
colnames(Lineage5_p3)<-c("SNVID","freq","passage")
Lineage5_p4<-Mosquito_L1thruL5[,c(1,24)]
Lineage5_p4$passage<-4
colnames(Lineage5_p4)<-c("SNVID","freq","passage")
Lineage5_p10<-Mosquito_L1thruL5[,c(1,25)]
Lineage5_p10$passage<-10
colnames(Lineage5_p10)<-c("SNVID","freq","passage")
Lineage5_p11<-Mosquito_L1thruL5[,c(1,26)]
Lineage5_p11$passage<-11
colnames(Lineage5_p11)<-c("SNVID","freq","passage")

Lineage5<-rbind(Lineage5_p1,Lineage5_p3)
Lineage5<-rbind(Lineage5,Lineage5_p4)
Lineage5<-rbind(Lineage5,Lineage5_p10)
Mosquito_Lineage5<-rbind(Lineage5,Lineage5_p11)

Mosquito_Lineage5$freq<-as.numeric(Mosquito_Lineage5$freq)
Mosquito_Lineage5$passage<-as.numeric(Mosquito_Lineage5$passage)



######  Mouse Lineage By Passage     ######################################################
p1_L1<-read.table(file="/PATH/TO/rep1and2_refalignSCP1_1.vcf", sep="\t", header=FALSE)
p1_L1<-cSplit(p1_L1,"V8",sep=";", type.convert=FALSE)
p1_L1<-cSplit(p1_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L1$SNVID<-paste(p1_L1$V2,p1_L1$V4,sep=";")
p1_L1$SNVID<-paste(p1_L1$SNVID,p1_L1$V5,sep=">")
p1_L1<-p1_L1[,c(12,11)]
colnames(p1_L1)<-c("SNVID","freq")
p1_L1<-subset(p1_L1,p1_L1$freq>=0.003)

p3_L1<-read.table(file="/PATH/TO/rep1and2_refalignSCP3_1.vcf", sep="\t", header=FALSE)
p3_L1<-cSplit(p3_L1,"V8",sep=";", type.convert=FALSE)
p3_L1<-cSplit(p3_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L1$SNVID<-paste(p3_L1$V2,p3_L1$V4,sep=";")
p3_L1$SNVID<-paste(p3_L1$SNVID,p3_L1$V5,sep=">")
p3_L1<-p3_L1[,c(12,11)]
colnames(p3_L1)<-c("SNVID","freq")
p3_L1<-subset(p3_L1,p3_L1$freq>=0.003)

p13_L1<-merge(p1_L1, p3_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L1)<-c("SNVID","freq1","freq3")

p4_L1<-read.table(file="/PATH/TO/rep1and2_refalignSCP4_1.vcf", sep="\t", header=FALSE)
p4_L1<-cSplit(p4_L1,"V8",sep=";", type.convert=FALSE)
p4_L1<-cSplit(p4_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L1$SNVID<-paste(p4_L1$V2,p4_L1$V4,sep=";")
p4_L1$SNVID<-paste(p4_L1$SNVID,p4_L1$V5,sep=">")
p4_L1<-p4_L1[,c(12,11)]
colnames(p4_L1)<-c("SNVID","freq")
p4_L1<-subset(p4_L1,p4_L1$freq>=0.003)

p14_L1<-merge(p13_L1, p4_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L1)<-c("SNVID","freq1","freq3","freq4")

p10_L1<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_1.vcf", sep="\t", header=FALSE)
p10_L1<-cSplit(p10_L1,"V8",sep=";", type.convert=FALSE)
p10_L1<-cSplit(p10_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L1$SNVID<-paste(p10_L1$V2,p10_L1$V4,sep=";")
p10_L1$SNVID<-paste(p10_L1$SNVID,p10_L1$V5,sep=">")
p10_L1<-p10_L1[,c(12,11)]
colnames(p10_L1)<-c("SNVID","freq")
p10_L1<-subset(p10_L1,p10_L1$freq>=0.003)

p110_L1<-merge(p14_L1, p10_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L1)<-c("SNVID","freq1","freq2","freq3","freq4")

p11_L1<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_A.vcf", sep="\t", header=FALSE)
p11_L1<-cSplit(p11_L1,"V8",sep=";", type.convert=FALSE)
p11_L1<-cSplit(p11_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L1$SNVID<-paste(p11_L1$V2,p11_L1$V4,sep=";")
p11_L1$SNVID<-paste(p11_L1$SNVID,p11_L1$V5,sep=">")
p11_L1<-p11_L1[,c(12,11)]
colnames(p11_L1)<-c("SNVID","freq")
p11_L1<-subset(p11_L1,p11_L1$freq>=0.003)

p111_L1<-merge(p110_L1, p11_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L1)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p111_L1[is.na(p111_L1)]<-0


p1_L2<-read.table(file="/PATH/TO/rep1and2_refalignSCP1_2.vcf", sep="\t", header=FALSE)
p1_L2<-cSplit(p1_L2,"V8",sep=";", type.convert=FALSE)
p1_L2<-cSplit(p1_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L2$SNVID<-paste(p1_L2$V2,p1_L2$V4,sep=";")
p1_L2$SNVID<-paste(p1_L2$SNVID,p1_L2$V5,sep=">")
p1_L2<-p1_L2[,c(12,11)]
colnames(p1_L2)<-c("SNVID","freq")
p1_L2<-subset(p1_L2,p1_L2$freq>=0.003)

p3_L2<-read.table(file="/PATH/TO/rep1and2_refalignSCP3_2.vcf", sep="\t", header=FALSE)
p3_L2<-cSplit(p3_L2,"V8",sep=";", type.convert=FALSE)
p3_L2<-cSplit(p3_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L2$SNVID<-paste(p3_L2$V2,p3_L2$V4,sep=";")
p3_L2$SNVID<-paste(p3_L2$SNVID,p3_L2$V5,sep=">")
p3_L2<-p3_L2[,c(12,11)]
colnames(p3_L2)<-c("SNVID","freq")
p3_L2<-subset(p3_L2,p3_L2$freq>=0.003)

p13_L2<-merge(p1_L2, p3_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L2)<-c("SNVID","freq1","freq3")

p4_L2<-read.table(file="/PATH/TO/rep1and2_refalignSCP4_2.vcf", sep="\t", header=FALSE)
p4_L2<-cSplit(p4_L2,"V8",sep=";", type.convert=FALSE)
p4_L2<-cSplit(p4_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L2$SNVID<-paste(p4_L2$V2,p4_L2$V4,sep=";")
p4_L2$SNVID<-paste(p4_L2$SNVID,p4_L2$V5,sep=">")
p4_L2<-p4_L2[,c(12,11)]
colnames(p4_L2)<-c("SNVID","freq")
p4_L2<-subset(p4_L2,p4_L2$freq>=0.003)

p14_L2<-merge(p13_L2, p4_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L2)<-c("SNVID","freq1","freq3","freq4")

p10_L2<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_2.vcf", sep="\t", header=FALSE)
p10_L2<-cSplit(p10_L2,"V8",sep=";", type.convert=FALSE)
p10_L2<-cSplit(p10_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L2$SNVID<-paste(p10_L2$V2,p10_L2$V4,sep=";")
p10_L2$SNVID<-paste(p10_L2$SNVID,p10_L2$V5,sep=">")
p10_L2<-p10_L2[,c(12,11)]
colnames(p10_L2)<-c("SNVID","freq")
p10_L2<-subset(p10_L2,p10_L2$freq>=0.003)

p110_L2<-merge(p14_L2, p10_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L2)<-c("SNVID","freq1","freq2","freq3","freq4")

p11_L2<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_B.vcf", sep="\t", header=FALSE)
p11_L2<-cSplit(p11_L2,"V8",sep=";", type.convert=FALSE)
p11_L2<-cSplit(p11_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L2$SNVID<-paste(p11_L2$V2,p11_L2$V4,sep=";")
p11_L2$SNVID<-paste(p11_L2$SNVID,p11_L2$V5,sep=">")
p11_L2<-p11_L2[,c(12,11)]
colnames(p11_L2)<-c("SNVID","freq")
p11_L2<-subset(p11_L2,p11_L2$freq>=0.003)

p111_L2<-merge(p110_L2, p11_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L2)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p111_L2[is.na(p111_L2)]<-0


p1_L3<-read.table(file="/PATH/TO/rep1and2_refalignSCP1_3.vcf", sep="\t", header=FALSE)
p1_L3<-cSplit(p1_L3,"V8",sep=";", type.convert=FALSE)
p1_L3<-cSplit(p1_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L3$SNVID<-paste(p1_L3$V2,p1_L3$V4,sep=";")
p1_L3$SNVID<-paste(p1_L3$SNVID,p1_L3$V5,sep=">")
p1_L3<-p1_L3[,c(12,11)]
colnames(p1_L3)<-c("SNVID","freq")
p1_L3<-subset(p1_L3,p1_L3$freq>=0.003)

p3_L3<-read.table(file="/PATH/TO/rep1and2_refalignSCP3_3.vcf", sep="\t", header=FALSE)
p3_L3<-cSplit(p3_L3,"V8",sep=";", type.convert=FALSE)
p3_L3<-cSplit(p3_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L3$SNVID<-paste(p3_L3$V2,p3_L3$V4,sep=";")
p3_L3$SNVID<-paste(p3_L3$SNVID,p3_L3$V5,sep=">")
p3_L3<-p3_L3[,c(12,11)]
colnames(p3_L3)<-c("SNVID","freq")
p3_L3<-subset(p3_L3,p3_L3$freq>=0.003)

p13_L3<-merge(p1_L3, p3_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L3)<-c("SNVID","freq1","freq3")

p4_L3<-read.table(file="/PATH/TO/rep1and2_refalignSCP4_3.vcf", sep="\t", header=FALSE)
p4_L3<-cSplit(p4_L3,"V8",sep=";", type.convert=FALSE)
p4_L3<-cSplit(p4_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L3$SNVID<-paste(p4_L3$V2,p4_L3$V4,sep=";")
p4_L3$SNVID<-paste(p4_L3$SNVID,p4_L3$V5,sep=">")
p4_L3<-p4_L3[,c(12,11)]
colnames(p4_L3)<-c("SNVID","freq")
p4_L3<-subset(p4_L3,p4_L3$freq>=0.003)

p14_L3<-merge(p13_L3, p4_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L3)<-c("SNVID","freq1","freq3","freq4")

p10_L3<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_3.vcf", sep="\t", header=FALSE)
p10_L3<-cSplit(p10_L3,"V8",sep=";", type.convert=FALSE)
p10_L3<-cSplit(p10_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L3$SNVID<-paste(p10_L3$V2,p10_L3$V4,sep=";")
p10_L3$SNVID<-paste(p10_L3$SNVID,p10_L3$V5,sep=">")
p10_L3<-p10_L3[,c(12,11)]
colnames(p10_L3)<-c("SNVID","freq")
p10_L3<-subset(p10_L3,p10_L3$freq>=0.003)

p110_L3<-merge(p14_L3, p10_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L3)<-c("SNVID","freq1","freq3","freq4","freq10")

p11_L3<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_C.vcf", sep="\t", header=FALSE)
p11_L3<-cSplit(p11_L3,"V8",sep=";", type.convert=FALSE)
p11_L3<-cSplit(p11_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L3$SNVID<-paste(p11_L3$V2,p11_L3$V4,sep=";")
p11_L3$SNVID<-paste(p11_L3$SNVID,p11_L3$V5,sep=">")
p11_L3<-p11_L3[,c(12,11)]
colnames(p11_L3)<-c("SNVID","freq")
p11_L3<-subset(p11_L3,p11_L3$freq>=0.003)

p111_L3<-merge(p110_L3, p11_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L3)<-c("SNVID","freq1","freq3","freq4","freq10","freq11")

p111_L3[is.na(p111_L3)]<-0

p1_L4<-read.table(file="/PATH/TO/rep1and2_refalignSCP1_4.vcf", sep="\t", header=FALSE)
p1_L4<-cSplit(p1_L4,"V8",sep=";", type.convert=FALSE)
p1_L4<-cSplit(p1_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L4$SNVID<-paste(p1_L4$V2,p1_L4$V4,sep=";")
p1_L4$SNVID<-paste(p1_L4$SNVID,p1_L4$V5,sep=">")
p1_L4<-p1_L4[,c(12,11)]
colnames(p1_L4)<-c("SNVID","freq")
p1_L4<-subset(p1_L4,p1_L4$freq>=0.003)

p3_L4<-read.table(file="/PATH/TO/rep1and2_refalignSCP3_4.vcf", sep="\t", header=FALSE)
p3_L4<-cSplit(p3_L4,"V8",sep=";", type.convert=FALSE)
p3_L4<-cSplit(p3_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L4$SNVID<-paste(p3_L4$V2,p3_L4$V4,sep=";")
p3_L4$SNVID<-paste(p3_L4$SNVID,p3_L4$V5,sep=">")
p3_L4<-p3_L4[,c(12,11)]
colnames(p3_L4)<-c("SNVID","freq")
p3_L4<-subset(p3_L4,p3_L4$freq>=0.003)

p13_L4<-merge(p1_L4, p3_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L4)<-c("SNVID","freq1","freq3")

p4_L4<-read.table(file="/PATH/TO/rep1and2_refalignSCP4_4.vcf", sep="\t", header=FALSE)
p4_L4<-cSplit(p4_L4,"V8",sep=";", type.convert=FALSE)
p4_L4<-cSplit(p4_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L4$SNVID<-paste(p4_L4$V2,p4_L4$V4,sep=";")
p4_L4$SNVID<-paste(p4_L4$SNVID,p4_L4$V5,sep=">")
p4_L4<-p4_L4[,c(12,11)]
colnames(p4_L4)<-c("SNVID","freq")
p4_L4<-subset(p4_L4,p4_L4$freq>=0.003)

p14_L4<-merge(p13_L4, p4_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L4)<-c("SNVID","freq1","freq3","freq4")

p10_L4<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_4.vcf", sep="\t", header=FALSE)
p10_L4<-cSplit(p10_L4,"V8",sep=";", type.convert=FALSE)
p10_L4<-cSplit(p10_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L4$SNVID<-paste(p10_L4$V2,p10_L4$V4,sep=";")
p10_L4$SNVID<-paste(p10_L4$SNVID,p10_L4$V5,sep=">")
p10_L4<-p10_L4[,c(12,11)]
colnames(p10_L4)<-c("SNVID","freq")
p10_L4<-subset(p10_L4,p10_L4$freq>=0.003)

p110_L4<-merge(p14_L4, p10_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L4)<-c("SNVID","freq1","freq2","freq3","freq4")

p11_L4<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_D.vcf", sep="\t", header=FALSE)
p11_L4<-cSplit(p11_L4,"V8",sep=";", type.convert=FALSE)
p11_L4<-cSplit(p11_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L4$SNVID<-paste(p11_L4$V2,p11_L4$V4,sep=";")
p11_L4$SNVID<-paste(p11_L4$SNVID,p11_L4$V5,sep=">")
p11_L4<-p11_L4[,c(12,11)]
colnames(p11_L4)<-c("SNVID","freq")
p11_L4<-subset(p11_L4,p11_L4$freq>=0.003)

p111_L4<-merge(p110_L4, p11_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L4)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p111_L4[is.na(p111_L4)]<-0


p1_L5<-read.table(file="/PATH/TO/rep1and2_refalignSCP1_5.vcf", sep="\t", header=FALSE)
p1_L5<-cSplit(p1_L5,"V8",sep=";", type.convert=FALSE)
p1_L5<-cSplit(p1_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L5$SNVID<-paste(p1_L5$V2,p1_L5$V4,sep=";")
p1_L5$SNVID<-paste(p1_L5$SNVID,p1_L5$V5,sep=">")
p1_L5<-p1_L5[,c(12,11)]
colnames(p1_L5)<-c("SNVID","freq")
p1_L5<-subset(p1_L5,p1_L5$freq>=0.003)

p3_L5<-read.table(file="/PATH/TO/rep1and2_refalignSCP3_5.vcf", sep="\t", header=FALSE)
p3_L5<-cSplit(p3_L5,"V8",sep=";", type.convert=FALSE)
p3_L5<-cSplit(p3_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L5$SNVID<-paste(p3_L5$V2,p3_L5$V4,sep=";")
p3_L5$SNVID<-paste(p3_L5$SNVID,p3_L5$V5,sep=">")
p3_L5<-p3_L5[,c(12,11)]
colnames(p3_L5)<-c("SNVID","freq")
p3_L5<-subset(p3_L5,p3_L5$freq>=0.003)

p13_L5<-merge(p1_L5, p3_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L5)<-c("SNVID","freq1","freq3")

p4_L5<-read.table(file="/PATH/TO/rep1and2_refalignSCP4_5.vcf", sep="\t", header=FALSE)
p4_L5<-cSplit(p4_L5,"V8",sep=";", type.convert=FALSE)
p4_L5<-cSplit(p4_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L5$SNVID<-paste(p4_L5$V2,p4_L5$V4,sep=";")
p4_L5$SNVID<-paste(p4_L5$SNVID,p4_L5$V5,sep=">")
p4_L5<-p4_L5[,c(12,11)]
colnames(p4_L5)<-c("SNVID","freq")
p4_L5<-subset(p4_L5,p4_L5$freq>=0.003)

p14_L5<-merge(p13_L5, p4_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L5)<-c("SNVID","freq1","freq3","freq4")

p10_L5<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_5.vcf", sep="\t", header=FALSE)
p10_L5<-cSplit(p10_L5,"V8",sep=";", type.convert=FALSE)
p10_L5<-cSplit(p10_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L5$SNVID<-paste(p10_L5$V2,p10_L5$V4,sep=";")
p10_L5$SNVID<-paste(p10_L5$SNVID,p10_L5$V5,sep=">")
p10_L5<-p10_L5[,c(12,11)]
colnames(p10_L5)<-c("SNVID","freq")
p10_L5<-subset(p10_L5,p10_L5$freq>=0.003)

p110_L5<-merge(p14_L5, p10_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L5)<-c("SNVID","freq1","freq2","freq3","freq4")

p11_L5<-read.table(file="/PATH/TO/rep1and2_refalignSCP10_E.vcf", sep="\t", header=FALSE)
p11_L5<-cSplit(p11_L5,"V8",sep=";", type.convert=FALSE)
p11_L5<-cSplit(p11_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L5$SNVID<-paste(p11_L5$V2,p11_L5$V4,sep=";")
p11_L5$SNVID<-paste(p11_L5$SNVID,p11_L5$V5,sep=">")
p11_L5<-p11_L5[,c(12,11)]
colnames(p11_L5)<-c("SNVID","freq")
p11_L5<-subset(p11_L5,p11_L5$freq>=0.003)

p111_L5<-merge(p110_L5, p11_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L5)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p111_L5[is.na(p111_L5)]<-0


L1andL2<-merge(p111_L1, p111_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1andL2)<-c("SNVID","p1_L1","p3_L1","p4_L1","p10_L1","p11_L1","p1_L2","p3_L2","p4_L2","p10_L2","p11_L2")
L1thruL3<-merge(L1andL2,p111_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1thruL3)<-c("SNVID","p1_L1","p3_L1","p4_L1","p10_L1","p11_L1","p1_L2","p3_L2","p4_L2","p10_L2","p11_L2","p1_L3","p3_L3","p4_L3","p10_L3","p11_L3")
L1thruL4<-merge(L1thruL3,p111_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1thruL4)<-c("SNVID","p1_L1","p3_L1","p4_L1","p10_L1","p11_L1","p1_L2","p3_L2","p4_L2","p10_L2","p11_L2","p1_L3","p3_L3","p4_L3","p10_L3","p11_L3","p1_L4","p3_L4","p4_L4","p10_L4","p11_L4")
Mouse_L1thruL5<-merge(L1thruL4,p111_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(Mouse_L1thruL5)<-c("SNVID","p1_L1","p3_L1","p4_L1","p10_L1","p11_L1","p1_L2","p3_L2","p4_L2","p10_L2","p11_L2","p1_L3","p3_L3","p4_L3","p10_L3","p11_L3","p1_L4","p3_L4","p4_L4","p10_L4","p11_L4","p1_L5","p3_L5","p4_L5","p10_L5","p11_L5")

Mouse_L1thruL5[is.na(Mouse_L1thruL5)]<-0

Lineage1_p1<-Mouse_L1thruL5[,c(1,2)]
Lineage1_p1$passage<-1
colnames(Lineage1_p1)<-c("SNVID","freq","passage")
Lineage1_p3<-Mouse_L1thruL5[,c(1,3)]
Lineage1_p3$passage<-3
colnames(Lineage1_p3)<-c("SNVID","freq","passage")
Lineage1_p4<-Mouse_L1thruL5[,c(1,4)]
Lineage1_p4$passage<-4
colnames(Lineage1_p4)<-c("SNVID","freq","passage")
Lineage1_p10<-Mouse_L1thruL5[,c(1,5)]
Lineage1_p10$passage<-10
colnames(Lineage1_p10)<-c("SNVID","freq","passage")
Lineage1_p11<-Mouse_L1thruL5[,c(1,6)]
Lineage1_p11$passage<-11
colnames(Lineage1_p11)<-c("SNVID","freq","passage")

Lineage1<-rbind(Lineage1_p1,Lineage1_p3)
Lineage1<-rbind(Lineage1,Lineage1_p4)
Lineage1<-rbind(Lineage1,Lineage1_p10)
Mouse_Lineage1<-rbind(Lineage1,Lineage1_p11)

Mouse_Lineage1$freq<-as.numeric(Mouse_Lineage1$freq)
Mouse_Lineage1$passage<-as.numeric(Mouse_Lineage1$passage)

Lineage2_p1<-Mouse_L1thruL5[,c(1,7)]
Lineage2_p1$passage<-1
colnames(Lineage2_p1)<-c("SNVID","freq","passage")
Lineage2_p3<-Mouse_L1thruL5[,c(1,8)]
Lineage2_p3$passage<-3
colnames(Lineage2_p3)<-c("SNVID","freq","passage")
Lineage2_p4<-Mouse_L1thruL5[,c(1,9)]
Lineage2_p4$passage<-4
colnames(Lineage2_p4)<-c("SNVID","freq","passage")
Lineage2_p10<-Mouse_L1thruL5[,c(1,10)]
Lineage2_p10$passage<-10
colnames(Lineage2_p10)<-c("SNVID","freq","passage")
Lineage2_p11<-Mouse_L1thruL5[,c(1,11)]
Lineage2_p11$passage<-11
colnames(Lineage2_p11)<-c("SNVID","freq","passage")

Lineage2<-rbind(Lineage2_p1,Lineage2_p3)
Lineage2<-rbind(Lineage2,Lineage2_p4)
Lineage2<-rbind(Lineage2,Lineage2_p10)
Mouse_Lineage2<-rbind(Lineage2,Lineage2_p11)

Mouse_Lineage2$freq<-as.numeric(Mouse_Lineage2$freq)
Mouse_Lineage2$passage<-as.numeric(Mouse_Lineage2$passage)

Lineage3_p1<-Mouse_L1thruL5[,c(1,12)]
Lineage3_p1$passage<-1
colnames(Lineage3_p1)<-c("SNVID","freq","passage")
Lineage3_p3<-Mouse_L1thruL5[,c(1,13)]
Lineage3_p3$passage<-3
colnames(Lineage3_p3)<-c("SNVID","freq","passage")
Lineage3_p4<-Mouse_L1thruL5[,c(1,14)]
Lineage3_p4$passage<-4
colnames(Lineage3_p4)<-c("SNVID","freq","passage")
Lineage3_p10<-Mouse_L1thruL5[,c(1,15)]
Lineage3_p10$passage<-10
colnames(Lineage3_p10)<-c("SNVID","freq","passage")
Lineage3_p11<-Mouse_L1thruL5[,c(1,16)]
Lineage3_p11$passage<-11
colnames(Lineage3_p11)<-c("SNVID","freq","passage")

Lineage3<-rbind(Lineage3_p1,Lineage3_p3)
Lineage3<-rbind(Lineage3,Lineage3_p4)
Lineage3<-rbind(Lineage3,Lineage3_p10)
Mouse_Lineage3<-rbind(Lineage3,Lineage3_p11)

Mouse_Lineage3$freq<-as.numeric(Mouse_Lineage3$freq)
Mouse_Lineage3$passage<-as.numeric(Mouse_Lineage3$passage)

Lineage4_p1<-Mouse_L1thruL5[,c(1,17)]
Lineage4_p1$passage<-1
colnames(Lineage4_p1)<-c("SNVID","freq","passage")
Lineage4_p3<-Mouse_L1thruL5[,c(1,18)]
Lineage4_p3$passage<-3
colnames(Lineage4_p3)<-c("SNVID","freq","passage")
Lineage4_p4<-Mouse_L1thruL5[,c(1,19)]
Lineage4_p4$passage<-4
colnames(Lineage4_p4)<-c("SNVID","freq","passage")
Lineage4_p10<-Mouse_L1thruL5[,c(1,20)]
Lineage4_p10$passage<-10
colnames(Lineage4_p10)<-c("SNVID","freq","passage")
Lineage4_p11<-Mouse_L1thruL5[,c(1,21)]
Lineage4_p11$passage<-11
colnames(Lineage4_p11)<-c("SNVID","freq","passage")

Lineage4<-rbind(Lineage4_p1,Lineage4_p3)
Lineage4<-rbind(Lineage4,Lineage4_p4)
Lineage4<-rbind(Lineage4,Lineage4_p10)
Mouse_Lineage4<-rbind(Lineage4,Lineage4_p11)

Mouse_Lineage4$freq<-as.numeric(Mouse_Lineage4$freq)
Mouse_Lineage4$passage<-as.numeric(Mouse_Lineage4$passage)

Lineage5_p1<-Mouse_L1thruL5[,c(1,22)]
Lineage5_p1$passage<-1
colnames(Lineage5_p1)<-c("SNVID","freq","passage")
Lineage5_p3<-Mouse_L1thruL5[,c(1,23)]
Lineage5_p3$passage<-3
colnames(Lineage5_p3)<-c("SNVID","freq","passage")
Lineage5_p4<-Mouse_L1thruL5[,c(1,24)]
Lineage5_p4$passage<-4
colnames(Lineage5_p4)<-c("SNVID","freq","passage")
Lineage5_p10<-Mouse_L1thruL5[,c(1,25)]
Lineage5_p10$passage<-10
colnames(Lineage5_p10)<-c("SNVID","freq","passage")
Lineage5_p11<-Mouse_L1thruL5[,c(1,26)]
Lineage5_p11$passage<-11
colnames(Lineage5_p11)<-c("SNVID","freq","passage")

Lineage5<-rbind(Lineage5_p1,Lineage5_p3)
Lineage5<-rbind(Lineage5,Lineage5_p4)
Lineage5<-rbind(Lineage5,Lineage5_p10)
Mouse_Lineage5<-rbind(Lineage5,Lineage5_p11)

Mouse_Lineage5$freq<-as.numeric(Mouse_Lineage5$freq)
Mouse_Lineage5$passage<-as.numeric(Mouse_Lineage5$passage)


######  Alternate Lineage By Passage     ######################################################

p1_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP1_1.vcf", sep="\t", header=FALSE)
p1_L1<-cSplit(p1_L1,"V8",sep=";", type.convert=FALSE)
p1_L1<-cSplit(p1_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L1$SNVID<-paste(p1_L1$V2,p1_L1$V4,sep=";")
p1_L1$SNVID<-paste(p1_L1$SNVID,p1_L1$V5,sep=">")
p1_L1<-p1_L1[,c(12,11)]
colnames(p1_L1)<-c("SNVID","freq")
p1_L1<-subset(p1_L1,p1_L1$freq>=0.003)

p2_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP2_1_3B.vcf", sep="\t", header=FALSE)
p2_L1<-cSplit(p2_L1,"V8",sep=";", type.convert=FALSE)
p2_L1<-cSplit(p2_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2_L1$SNVID<-paste(p2_L1$V2,p2_L1$V4,sep=";")
p2_L1$SNVID<-paste(p2_L1$SNVID,p2_L1$V5,sep=">")
p2_L1<-p2_L1[,c(12,11)]
colnames(p2_L1)<-c("SNVID","freq")
p2_L1<-subset(p2_L1,p2_L1$freq>=0.003)

p12_L1<-merge(p1_L1, p2_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12_L1)<-c("SNVID","freq1","freq2B")

p2S_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP2_1_3S.vcf", sep="\t", header=FALSE)
p2S_L1<-cSplit(p2S_L1,"V8",sep=";", type.convert=FALSE)
p2S_L1<-cSplit(p2S_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2S_L1$SNVID<-paste(p2S_L1$V2,p2S_L1$V4,sep=";")
p2S_L1$SNVID<-paste(p2S_L1$SNVID,p2S_L1$V5,sep=">")
p2S_L1<-p2S_L1[,c(12,11)]
colnames(p2S_L1)<-c("SNVID","freq")
p2S_L1<-subset(p2S_L1,p2S_L1$freq>=0.003)

p12S_L1<-merge(p12_L1, p2S_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12S_L1)<-c("SNVID","freq1","freq2B","freq2S")

p3_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP3_1.vcf", sep="\t", header=FALSE)
p3_L1<-cSplit(p3_L1,"V8",sep=";", type.convert=FALSE)
p3_L1<-cSplit(p3_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L1$SNVID<-paste(p3_L1$V2,p3_L1$V4,sep=";")
p3_L1$SNVID<-paste(p3_L1$SNVID,p3_L1$V5,sep=">")
p3_L1<-p3_L1[,c(12,11)]
colnames(p3_L1)<-c("SNVID","freq")
p3_L1<-subset(p3_L1,p3_L1$freq>=0.003)

p13_L1<-merge(p12S_L1, p3_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L1)<-c("SNVID","freq1","freq2B","freq2S","freq3")

p4_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP4_1_B.vcf", sep="\t", header=FALSE)
p4_L1<-cSplit(p4_L1,"V8",sep=";", type.convert=FALSE)
p4_L1<-cSplit(p4_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L1$SNVID<-paste(p4_L1$V2,p4_L1$V4,sep=";")
p4_L1$SNVID<-paste(p4_L1$SNVID,p4_L1$V5,sep=">")
p4_L1<-p4_L1[,c(12,11)]
colnames(p4_L1)<-c("SNVID","freq")
p4_L1<-subset(p4_L1,p4_L1$freq>=0.003)

p14_L1<-merge(p13_L1, p4_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L1)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4")

p5_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP5_1.vcf", sep="\t", header=FALSE)
p5_L1<-cSplit(p5_L1,"V8",sep=";", type.convert=FALSE)
p5_L1<-cSplit(p5_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p5_L1$SNVID<-paste(p5_L1$V2,p5_L1$V4,sep=";")
p5_L1$SNVID<-paste(p5_L1$SNVID,p5_L1$V5,sep=">")
p5_L1<-p5_L1[,c(12,11)]
colnames(p5_L1)<-c("SNVID","freq")
p5_L1<-subset(p5_L1,p5_L1$freq>=0.003)

p15_L1<-merge(p14_L1, p5_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p15_L1)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5")

p6_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP6_1_B.vcf", sep="\t", header=FALSE)
p6_L1<-cSplit(p6_L1,"V8",sep=";", type.convert=FALSE)
p6_L1<-cSplit(p6_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p6_L1$SNVID<-paste(p6_L1$V2,p6_L1$V4,sep=";")
p6_L1$SNVID<-paste(p6_L1$SNVID,p6_L1$V5,sep=">")
p6_L1<-p6_L1[,c(12,11)]
colnames(p6_L1)<-c("SNVID","freq")
p6_L1<-subset(p6_L1,p6_L1$freq>=0.003)

p16_L1<-merge(p15_L1, p6_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p16_L1)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6")

p7_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP7_1.vcf", sep="\t", header=FALSE)
p7_L1<-cSplit(p7_L1,"V8",sep=";", type.convert=FALSE)
p7_L1<-cSplit(p7_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p7_L1$SNVID<-paste(p7_L1$V2,p7_L1$V4,sep=";")
p7_L1$SNVID<-paste(p7_L1$SNVID,p7_L1$V5,sep=">")
p7_L1<-p7_L1[,c(12,11)]
colnames(p7_L1)<-c("SNVID","freq")
p7_L1<-subset(p7_L1,p7_L1$freq>=0.003)

p17_L1<-merge(p16_L1, p7_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p17_L1)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7")

p8_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP8_1_B.vcf", sep="\t", header=FALSE)
p8_L1<-cSplit(p8_L1,"V8",sep=";", type.convert=FALSE)
p8_L1<-cSplit(p8_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p8_L1$SNVID<-paste(p8_L1$V2,p8_L1$V4,sep=";")
p8_L1$SNVID<-paste(p8_L1$SNVID,p8_L1$V5,sep=">")
p8_L1<-p8_L1[,c(12,11)]
colnames(p8_L1)<-c("SNVID","freq")
p8_L1<-subset(p8_L1,p8_L1$freq>=0.003)

p18_L1<-merge(p17_L1, p8_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p18_L1)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8")

p9_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP9_1.vcf", sep="\t", header=FALSE)
p9_L1<-cSplit(p9_L1,"V8",sep=";", type.convert=FALSE)
p9_L1<-cSplit(p9_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p9_L1$SNVID<-paste(p9_L1$V2,p9_L1$V4,sep=";")
p9_L1$SNVID<-paste(p9_L1$SNVID,p9_L1$V5,sep=">")
p9_L1<-p9_L1[,c(12,11)]
colnames(p9_L1)<-c("SNVID","freq")
p9_L1<-subset(p9_L1,p9_L1$freq>=0.003)

p19_L1<-merge(p18_L1, p9_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p19_L1)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8","freq9")

p10_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP10_1_B.vcf", sep="\t", header=FALSE)
p10_L1<-cSplit(p10_L1,"V8",sep=";", type.convert=FALSE)
p10_L1<-cSplit(p10_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L1$SNVID<-paste(p10_L1$V2,p10_L1$V4,sep=";")
p10_L1$SNVID<-paste(p10_L1$SNVID,p10_L1$V5,sep=">")
p10_L1<-p10_L1[,c(12,11)]
colnames(p10_L1)<-c("SNVID","freq")
p10_L1<-subset(p10_L1,p10_L1$freq>=0.003)

p110_L1<-merge(p19_L1, p10_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L1)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10")

p11_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP10_A.vcf", sep="\t", header=FALSE)
p11_L1<-cSplit(p11_L1,"V8",sep=";", type.convert=FALSE)
p11_L1<-cSplit(p11_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L1$SNVID<-paste(p11_L1$V2,p11_L1$V4,sep=";")
p11_L1$SNVID<-paste(p11_L1$SNVID,p11_L1$V5,sep=">")
p11_L1<-p11_L1[,c(12,11)]
colnames(p11_L1)<-c("SNVID","freq")
p11_L1<-subset(p11_L1,p11_L1$freq>=0.003)

p111_L1<-merge(p110_L1, p11_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L1)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10","freqCC")

p111_L1[is.na(p111_L1)]<-0


p1_L2<-read.table(file="/PATH/TO/rep1and2_refalignAP1_2.vcf", sep="\t", header=FALSE)
p1_L2<-cSplit(p1_L2,"V8",sep=";", type.convert=FALSE)
p1_L2<-cSplit(p1_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L2$SNVID<-paste(p1_L2$V2,p1_L2$V4,sep=";")
p1_L2$SNVID<-paste(p1_L2$SNVID,p1_L2$V5,sep=">")
p1_L2<-p1_L2[,c(12,11)]
colnames(p1_L2)<-c("SNVID","freq")
p1_L2<-subset(p1_L2,p1_L2$freq>=0.003)

p2_L2<-read.table(file="/PATH/TO/rep1and2_refalignAP2_2_3B.vcf", sep="\t", header=FALSE)
p2_L2<-cSplit(p2_L2,"V8",sep=";", type.convert=FALSE)
p2_L2<-cSplit(p2_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2_L2$SNVID<-paste(p2_L2$V2,p2_L2$V4,sep=";")
p2_L2$SNVID<-paste(p2_L2$SNVID,p2_L2$V5,sep=">")
p2_L2<-p2_L2[,c(12,11)]
colnames(p2_L2)<-c("SNVID","freq")
p2_L2<-subset(p2_L2,p2_L2$freq>=0.003)

p12_L2<-merge(p1_L2, p2_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12_L2)<-c("SNVID","freq1","freq2")

p2S_L2<-read.table(file="/PATH/TO/rep1and2_refalignAP2_2_3B.vcf", sep="\t", header=FALSE)
p2S_L2<-cSplit(p2S_L2,"V8",sep=";", type.convert=FALSE)
p2S_L2<-cSplit(p2S_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2S_L2$SNVID<-paste(p2S_L2$V2,p2S_L2$V4,sep=";")
p2S_L2$SNVID<-paste(p2S_L2$SNVID,p2S_L2$V5,sep=">")
p2S_L2<-p2S_L2[,c(12,11)]
colnames(p2S_L2)<-c("SNVID","freq")
p2S_L2<-subset(p2S_L2,p2S_L2$freq>=0.003)

p12S_L2<-merge(p12_L2, p2S_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12S_L2)<-c("SNVID","freq1","freq2B","freq2S")

p3_L2<-read.table(file="/PATH/TO/rep1and2_refalignAP3_2.vcf", sep="\t", header=FALSE)
p3_L2<-cSplit(p3_L2,"V8",sep=";", type.convert=FALSE)
p3_L2<-cSplit(p3_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L2$SNVID<-paste(p3_L2$V2,p3_L2$V4,sep=";")
p3_L2$SNVID<-paste(p3_L2$SNVID,p3_L2$V5,sep=">")
p3_L2<-p3_L2[,c(12,11)]
colnames(p3_L2)<-c("SNVID","freq")
p3_L2<-subset(p3_L2,p3_L2$freq>=0.003)

p13_L2<-merge(p12S_L2, p3_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L2)<-c("SNVID","freq1","freq2B","freq2S","freq3")

p13_L2[is.na(p13_L2)]<-0


p1_L3<-read.table(file="/PATH/TO/rep1and2_refalignAP1_3.vcf", sep="\t", header=FALSE)
p1_L3<-cSplit(p1_L3,"V8",sep=";", type.convert=FALSE)
p1_L3<-cSplit(p1_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L3$SNVID<-paste(p1_L3$V2,p1_L3$V4,sep=";")
p1_L3$SNVID<-paste(p1_L3$SNVID,p1_L3$V5,sep=">")
p1_L3<-p1_L3[,c(12,11)]
colnames(p1_L3)<-c("SNVID","freq")
p1_L3<-subset(p1_L3,p1_L3$freq>=0.003)

p2_L3<-read.table(file="/PATH/TO/rep1and2_refalignAP3_3.vcf", sep="\t", header=FALSE)
p2_L3<-cSplit(p2_L3,"V8",sep=";", type.convert=FALSE)
p2_L3<-cSplit(p2_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2_L3$SNVID<-paste(p2_L3$V2,p2_L3$V4,sep=";")
p2_L3$SNVID<-paste(p2_L3$SNVID,p2_L3$V5,sep=">")
p2_L3<-p2_L3[,c(12,11)]
colnames(p2_L3)<-c("SNVID","freq")
p2_L3<-subset(p2_L3,p2_L3$freq>=0.003)

p12_L3<-merge(p1_L3, p2_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12_L3)<-c("SNVID","freq1","freq2")

p2S_L3<-read.table(file="/PATH/TO/rep1and2_refalignAP3_3.vcf", sep="\t", header=FALSE)
p2S_L3<-cSplit(p2S_L3,"V8",sep=";", type.convert=FALSE)
p2S_L3<-cSplit(p2S_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2S_L3$SNVID<-paste(p2S_L3$V2,p2S_L3$V4,sep=";")
p2S_L3$SNVID<-paste(p2S_L3$SNVID,p2S_L3$V5,sep=">")
p2S_L3<-p2S_L3[,c(12,11)]
colnames(p2S_L3)<-c("SNVID","freq")
p2S_L3<-subset(p2S_L3,p2S_L3$freq>=0.003)

p12S_L3<-merge(p12_L3, p2S_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12S_L3)<-c("SNVID","freq1","freq2B","freq2S")

p3_L3<-read.table(file="/PATH/TO/rep1and2_refalignAP3_3.vcf", sep="\t", header=FALSE)
p3_L3<-cSplit(p3_L3,"V8",sep=";", type.convert=FALSE)
p3_L3<-cSplit(p3_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L3$SNVID<-paste(p3_L3$V2,p3_L3$V4,sep=";")
p3_L3$SNVID<-paste(p3_L3$SNVID,p3_L3$V5,sep=">")
p3_L3<-p3_L3[,c(12,11)]
colnames(p3_L3)<-c("SNVID","freq")
p3_L3<-subset(p3_L3,p3_L3$freq>=0.003)

p13_L3<-merge(p12S_L3, p3_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L3)<-c("SNVID","freq1","freq2B","freq2S","freq3")

p13_L3[is.na(p13_L3)]<-0

p1_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP1_4.vcf", sep="\t", header=FALSE)
p1_L4<-cSplit(p1_L4,"V8",sep=";", type.convert=FALSE)
p1_L4<-cSplit(p1_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L4$SNVID<-paste(p1_L4$V2,p1_L4$V4,sep=";")
p1_L4$SNVID<-paste(p1_L4$SNVID,p1_L4$V5,sep=">")
p1_L4<-p1_L4[,c(12,11)]
colnames(p1_L4)<-c("SNVID","freq")
p1_L4<-subset(p1_L4,p1_L4$freq>=0.003)

p2_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP2_4_6B.vcf", sep="\t", header=FALSE)
p2_L4<-cSplit(p2_L4,"V8",sep=";", type.convert=FALSE)
p2_L4<-cSplit(p2_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2_L4$SNVID<-paste(p2_L4$V2,p2_L4$V4,sep=";")
p2_L4$SNVID<-paste(p2_L4$SNVID,p2_L4$V5,sep=">")
p2_L4<-p2_L4[,c(12,11)]
colnames(p2_L4)<-c("SNVID","freq")
p2_L4<-subset(p2_L4,p2_L4$freq>=0.003)

p12_L4<-merge(p1_L4, p2_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12_L4)<-c("SNVID","freq1","freq2")

p2S_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP2_4_6B.vcf", sep="\t", header=FALSE)
p2S_L4<-cSplit(p2S_L4,"V8",sep=";", type.convert=FALSE)
p2S_L4<-cSplit(p2S_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2S_L4$SNVID<-paste(p2S_L4$V2,p2S_L4$V4,sep=";")
p2S_L4$SNVID<-paste(p2S_L4$SNVID,p2S_L4$V5,sep=">")
p2S_L4<-p2S_L4[,c(12,11)]
colnames(p2S_L4)<-c("SNVID","freq")
p2S_L4<-subset(p2S_L4,p2S_L4$freq>=0.003)

p12S_L4<-merge(p12_L4, p2S_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12S_L4)<-c("SNVID","freq1","freq2B","freq2S")

p3_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP3_4.vcf", sep="\t", header=FALSE)
p3_L4<-cSplit(p3_L4,"V8",sep=";", type.convert=FALSE)
p3_L4<-cSplit(p3_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L4$SNVID<-paste(p3_L4$V2,p3_L4$V4,sep=";")
p3_L4$SNVID<-paste(p3_L4$SNVID,p3_L4$V5,sep=">")
p3_L4<-p3_L4[,c(12,11)]
colnames(p3_L4)<-c("SNVID","freq")
p3_L4<-subset(p3_L4,p3_L4$freq>=0.003)

p13_L4<-merge(p12S_L4, p3_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L4)<-c("SNVID","freq1","freq2B","freq2S","freq3")

p4_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP4_4_B.vcf", sep="\t", header=FALSE)
p4_L4<-cSplit(p4_L4,"V8",sep=";", type.convert=FALSE)
p4_L4<-cSplit(p4_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L4$SNVID<-paste(p4_L4$V2,p4_L4$V4,sep=";")
p4_L4$SNVID<-paste(p4_L4$SNVID,p4_L4$V5,sep=">")
p4_L4<-p4_L4[,c(12,11)]
colnames(p4_L4)<-c("SNVID","freq")
p4_L4<-subset(p4_L4,p4_L4$freq>=0.003)

p14_L4<-merge(p13_L4, p4_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L4)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4")

p5_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP5_4.vcf", sep="\t", header=FALSE)
p5_L4<-cSplit(p5_L4,"V8",sep=";", type.convert=FALSE)
p5_L4<-cSplit(p5_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p5_L4$SNVID<-paste(p5_L4$V2,p5_L4$V4,sep=";")
p5_L4$SNVID<-paste(p5_L4$SNVID,p5_L4$V5,sep=">")
p5_L4<-p5_L4[,c(12,11)]
colnames(p5_L4)<-c("SNVID","freq")
p5_L4<-subset(p5_L4,p5_L4$freq>=0.003)

p15_L4<-merge(p14_L4, p5_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p15_L4)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5")

p6_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP6_4_B.vcf", sep="\t", header=FALSE)
p6_L4<-cSplit(p6_L4,"V8",sep=";", type.convert=FALSE)
p6_L4<-cSplit(p6_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p6_L4$SNVID<-paste(p6_L4$V2,p6_L4$V4,sep=";")
p6_L4$SNVID<-paste(p6_L4$SNVID,p6_L4$V5,sep=">")
p6_L4<-p6_L4[,c(12,11)]
colnames(p6_L4)<-c("SNVID","freq")
p6_L4<-subset(p6_L4,p6_L4$freq>=0.003)

p16_L4<-merge(p15_L4, p6_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p16_L4)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6")

p7_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP7_4.vcf", sep="\t", header=FALSE)
p7_L4<-cSplit(p7_L4,"V8",sep=";", type.convert=FALSE)
p7_L4<-cSplit(p7_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p7_L4$SNVID<-paste(p7_L4$V2,p7_L4$V4,sep=";")
p7_L4$SNVID<-paste(p7_L4$SNVID,p7_L4$V5,sep=">")
p7_L4<-p7_L4[,c(12,11)]
colnames(p7_L4)<-c("SNVID","freq")
p7_L4<-subset(p7_L4,p7_L4$freq>=0.003)

p17_L4<-merge(p16_L4, p7_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p17_L4)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7")

p8_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP8_4_B.vcf", sep="\t", header=FALSE)
p8_L4<-cSplit(p8_L4,"V8",sep=";", type.convert=FALSE)
p8_L4<-cSplit(p8_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p8_L4$SNVID<-paste(p8_L4$V2,p8_L4$V4,sep=";")
p8_L4$SNVID<-paste(p8_L4$SNVID,p8_L4$V5,sep=">")
p8_L4<-p8_L4[,c(12,11)]
colnames(p8_L4)<-c("SNVID","freq")
p8_L4<-subset(p8_L4,p8_L4$freq>=0.003)

p18_L4<-merge(p17_L4, p8_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p18_L4)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8")

p9_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP9_4.vcf", sep="\t", header=FALSE)
p9_L4<-cSplit(p9_L4,"V8",sep=";", type.convert=FALSE)
p9_L4<-cSplit(p9_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p9_L4$SNVID<-paste(p9_L4$V2,p9_L4$V4,sep=";")
p9_L4$SNVID<-paste(p9_L4$SNVID,p9_L4$V5,sep=">")
p9_L4<-p9_L4[,c(12,11)]
colnames(p9_L4)<-c("SNVID","freq")
p9_L4<-subset(p9_L4,p9_L4$freq>=0.003)

p19_L4<-merge(p18_L4, p9_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p19_L4)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8","freq9")

p10_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP10_4_B.vcf", sep="\t", header=FALSE)
p10_L4<-cSplit(p10_L4,"V8",sep=";", type.convert=FALSE)
p10_L4<-cSplit(p10_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L4$SNVID<-paste(p10_L4$V2,p10_L4$V4,sep=";")
p10_L4$SNVID<-paste(p10_L4$SNVID,p10_L4$V5,sep=">")
p10_L4<-p10_L4[,c(12,11)]
colnames(p10_L4)<-c("SNVID","freq")
p10_L4<-subset(p10_L4,p10_L4$freq>=0.003)

p110_L4<-merge(p19_L4, p10_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L4)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10")

p11_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP10_D.vcf", sep="\t", header=FALSE)
p11_L4<-cSplit(p11_L4,"V8",sep=";", type.convert=FALSE)
p11_L4<-cSplit(p11_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L4$SNVID<-paste(p11_L4$V2,p11_L4$V4,sep=";")
p11_L4$SNVID<-paste(p11_L4$SNVID,p11_L4$V5,sep=">")
p11_L4<-p11_L4[,c(12,11)]
colnames(p11_L4)<-c("SNVID","freq")
p11_L4<-subset(p11_L4,p11_L4$freq>=0.003)

p111_L4<-merge(p110_L4, p11_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L4)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10","freqCC")

p111_L4[is.na(p111_L4)]<-0


p1_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP1_5.vcf", sep="\t", header=FALSE)
p1_L5<-cSplit(p1_L5,"V8",sep=";", type.convert=FALSE)
p1_L5<-cSplit(p1_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L5$SNVID<-paste(p1_L5$V2,p1_L5$V4,sep=";")
p1_L5$SNVID<-paste(p1_L5$SNVID,p1_L5$V5,sep=">")
p1_L5<-p1_L5[,c(12,11)]
colnames(p1_L5)<-c("SNVID","freq")
p1_L5<-subset(p1_L5,p1_L5$freq>=0.003)

p2_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP3_5.vcf", sep="\t", header=FALSE)
p2_L5<-cSplit(p2_L5,"V8",sep=";", type.convert=FALSE)
p2_L5<-cSplit(p2_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2_L5$SNVID<-paste(p2_L5$V2,p2_L5$V4,sep=";")
p2_L5$SNVID<-paste(p2_L5$SNVID,p2_L5$V5,sep=">")
p2_L5<-p2_L5[,c(12,11)]
colnames(p2_L5)<-c("SNVID","freq")
p2_L5<-subset(p2_L5,p2_L5$freq>=0.003)

p12_L5<-merge(p1_L5, p2_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12_L5)<-c("SNVID","freq1","freq2")

p2S_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP3_5.vcf", sep="\t", header=FALSE)
p2S_L5<-cSplit(p2S_L5,"V8",sep=";", type.convert=FALSE)
p2S_L5<-cSplit(p2S_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2S_L5$SNVID<-paste(p2S_L5$V2,p2S_L5$V4,sep=";")
p2S_L5$SNVID<-paste(p2S_L5$SNVID,p2S_L5$V5,sep=">")
p2S_L5<-p2S_L5[,c(12,11)]
colnames(p2S_L5)<-c("SNVID","freq")
p2S_L5<-subset(p2S_L5,p2S_L5$freq>=0.003)

p12S_L5<-merge(p12_L5, p2S_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12S_L5)<-c("SNVID","freq1","freq2B","freq2S")

p3_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP3_5.vcf", sep="\t", header=FALSE)
p3_L5<-cSplit(p3_L5,"V8",sep=";", type.convert=FALSE)
p3_L5<-cSplit(p3_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L5$SNVID<-paste(p3_L5$V2,p3_L5$V4,sep=";")
p3_L5$SNVID<-paste(p3_L5$SNVID,p3_L5$V5,sep=">")
p3_L5<-p3_L5[,c(12,11)]
colnames(p3_L5)<-c("SNVID","freq")
p3_L5<-subset(p3_L5,p3_L5$freq>=0.003)

p13_L5<-merge(p12S_L5, p3_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L5)<-c("SNVID","freq1","freq2B","freq2S","freq3")

p4_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP4_5_B.vcf", sep="\t", header=FALSE)
p4_L5<-cSplit(p4_L5,"V8",sep=";", type.convert=FALSE)
p4_L5<-cSplit(p4_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L5$SNVID<-paste(p4_L5$V2,p4_L5$V4,sep=";")
p4_L5$SNVID<-paste(p4_L5$SNVID,p4_L5$V5,sep=">")
p4_L5<-p4_L5[,c(12,11)]
colnames(p4_L5)<-c("SNVID","freq")
p4_L5<-subset(p4_L5,p4_L5$freq>=0.003)

p14_L5<-merge(p13_L5, p4_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L5)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4")

p5_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP5_5.vcf", sep="\t", header=FALSE)
p5_L5<-cSplit(p5_L5,"V8",sep=";", type.convert=FALSE)
p5_L5<-cSplit(p5_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p5_L5$SNVID<-paste(p5_L5$V2,p5_L5$V4,sep=";")
p5_L5$SNVID<-paste(p5_L5$SNVID,p5_L5$V5,sep=">")
p5_L5<-p5_L5[,c(12,11)]
colnames(p5_L5)<-c("SNVID","freq")
p5_L5<-subset(p5_L5,p5_L5$freq>=0.003)

p15_L5<-merge(p14_L5, p5_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p15_L5)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5")

p6_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP6_5_B.vcf", sep="\t", header=FALSE)
p6_L5<-cSplit(p6_L5,"V8",sep=";", type.convert=FALSE)
p6_L5<-cSplit(p6_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p6_L5$SNVID<-paste(p6_L5$V2,p6_L5$V4,sep=";")
p6_L5$SNVID<-paste(p6_L5$SNVID,p6_L5$V5,sep=">")
p6_L5<-p6_L5[,c(12,11)]
colnames(p6_L5)<-c("SNVID","freq")
p6_L5<-subset(p6_L5,p6_L5$freq>=0.003)

p16_L5<-merge(p15_L5, p6_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p16_L5)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6")

p7_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP7_5.vcf", sep="\t", header=FALSE)
p7_L5<-cSplit(p7_L5,"V8",sep=";", type.convert=FALSE)
p7_L5<-cSplit(p7_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p7_L5$SNVID<-paste(p7_L5$V2,p7_L5$V4,sep=";")
p7_L5$SNVID<-paste(p7_L5$SNVID,p7_L5$V5,sep=">")
p7_L5<-p7_L5[,c(12,11)]
colnames(p7_L5)<-c("SNVID","freq")
p7_L5<-subset(p7_L5,p7_L5$freq>=0.003)

p17_L5<-merge(p16_L5, p7_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p17_L5)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7")

p8_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP8_5_B.vcf", sep="\t", header=FALSE)
p8_L5<-cSplit(p8_L5,"V8",sep=";", type.convert=FALSE)
p8_L5<-cSplit(p8_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p8_L5$SNVID<-paste(p8_L5$V2,p8_L5$V4,sep=";")
p8_L5$SNVID<-paste(p8_L5$SNVID,p8_L5$V5,sep=">")
p8_L5<-p8_L5[,c(12,11)]
colnames(p8_L5)<-c("SNVID","freq")
p8_L5<-subset(p8_L5,p8_L5$freq>=0.003)

p18_L5<-merge(p17_L5, p8_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p18_L5)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8")

p9_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP9_5.vcf", sep="\t", header=FALSE)
p9_L5<-cSplit(p9_L5,"V8",sep=";", type.convert=FALSE)
p9_L5<-cSplit(p9_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p9_L5$SNVID<-paste(p9_L5$V2,p9_L5$V4,sep=";")
p9_L5$SNVID<-paste(p9_L5$SNVID,p9_L5$V5,sep=">")
p9_L5<-p9_L5[,c(12,11)]
colnames(p9_L5)<-c("SNVID","freq")
p9_L5<-subset(p9_L5,p9_L5$freq>=0.003)

p19_L5<-merge(p18_L5, p9_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p19_L5)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8","freq9")

p10_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP10_5_B.vcf", sep="\t", header=FALSE)
p10_L5<-cSplit(p10_L5,"V8",sep=";", type.convert=FALSE)
p10_L5<-cSplit(p10_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L5$SNVID<-paste(p10_L5$V2,p10_L5$V4,sep=";")
p10_L5$SNVID<-paste(p10_L5$SNVID,p10_L5$V5,sep=">")
p10_L5<-p10_L5[,c(12,11)]
colnames(p10_L5)<-c("SNVID","freq")
p10_L5<-subset(p10_L5,p10_L5$freq>=0.003)

p110_L5<-merge(p19_L5, p10_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L5)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10")

p11_L5<-read.table(file="/Volumes/HD1/R21_ZIKVConspecificPassage/Alternate_Passage/WGS/AP10_E_rep1_v3.1/reference_aligned/lofreq_reference_AP10_E_1.vcf", sep="\t", header=FALSE)
p11_L5<-cSplit(p11_L5,"V8",sep=";", type.convert=FALSE)
p11_L5<-cSplit(p11_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L5$SNVID<-paste(p11_L5$V2,p11_L5$V4,sep=";")
p11_L5$SNVID<-paste(p11_L5$SNVID,p11_L5$V5,sep=">")
p11_L5<-p11_L5[,c(14,13)]
colnames(p11_L5)<-c("SNVID","freq")
p11_L5<-subset(p11_L5,p11_L5$freq>=0.003)

p111_L5<-merge(p110_L5, p11_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L5)<-c("SNVID","freq1","freq2B","freq2S","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10","freqCC")

p111_L5[is.na(p111_L5)]<-0


L1andL2<-merge(p111_L1, p13_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1andL2)<-c("SNVID","p1_L1","p2B_L1","p2S_L1","p3_L1","p4_L1","p5_L1","p6_L1","p7_L1","p8_L1","p9_L1","p10_L1","p11_L1","p1_L2","p2B_L2","p2S_L2","p3_L2")
L1thruL3<-merge(L1andL2,p13_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1thruL3)<-c("SNVID","p1_L1","p2B_L1","p2S_L1","p3_L1","p4_L1","p5_L1","p6_L1","p7_L1","p8_L1","p9_L1","p10_L1","p11_L1","p1_L2","p2B_L2","p2S_L2","p3_L2","p1_L3","p2B_L3","p2S_L3","p3_L3")
L1thruL4<-merge(L1thruL3,p111_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1thruL4)<-c("SNVID","p1_L1","p2B_L1","p2S_L1","p3_L1","p4_L1","p5_L1","p6_L1","p7_L1","p8_L1","p9_L1","p10_L1","p11_L1","p1_L2","p2B_L2","p2S_L2","p3_L2","p1_L3","p2B_L3","p2S_L3","p3_L3","p1_L4","p2B_L4","p2S_L4","p3_L4","p4_L4","p5_L4","p6_L4","p7_L4","p8_L4","p9_L4","p10_L4","p11_L4")
Alternate_L1thruL5<-merge(L1thruL4,p111_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(Alternate_L1thruL5)<-c("SNVID","p1_L1","p2B_L1","p2S_L1","p3_L1","p4_L1","p5_L1","p6_L1","p7_L1","p8_L1","p9_L1","p10_L1","p11_L1","p1_L2","p2B_L2","p2S_L2","p3_L2","p1_L3","p2B_L3","p2S_L3","p3_L3","p1_L4","p2B_L4","p2S_L4","p3_L4","p4_L4","p5_L4","p6_L4","p7_L4","p8_L4","p9_L4","p10_L4","p11_L4","p1_L5","p2B_L5","p2S_L5","p3_L5","p4_L5","p5_L5","p6_L5","p7_L5","p8_L5","p9_L5","p10_L5","p11_L5")

Alternate_L1thruL5[is.na(Alternate_L1thruL5)]<-0

Lineage1_p1<-Alternate_L1thruL5[,c(1,2)]
Lineage1_p1$passage<-1
colnames(Lineage1_p1)<-c("SNVID","freq","passage")
Lineage1_p2B<-Alternate_L1thruL5[,c(1,3)]
Lineage1_p2B$passage<-2
colnames(Lineage1_p2B)<-c("SNVID","freq","passage")
Lineage1_p3<-Alternate_L1thruL5[,c(1,5)]
Lineage1_p3$passage<-3
colnames(Lineage1_p3)<-c("SNVID","freq","passage")
Lineage1_p4<-Alternate_L1thruL5[,c(1,6)]
Lineage1_p4$passage<-4
colnames(Lineage1_p4)<-c("SNVID","freq","passage")
Lineage1_p5<-Alternate_L1thruL5[,c(1,7)]
Lineage1_p5$passage<-5
colnames(Lineage1_p5)<-c("SNVID","freq","passage")
Lineage1_p6<-Alternate_L1thruL5[,c(1,8)]
Lineage1_p6$passage<-6
colnames(Lineage1_p6)<-c("SNVID","freq","passage")
Lineage1_p7<-Alternate_L1thruL5[,c(1,9)]
Lineage1_p7$passage<-7
colnames(Lineage1_p7)<-c("SNVID","freq","passage")
Lineage1_p8<-Alternate_L1thruL5[,c(1,10)]
Lineage1_p8$passage<-8
colnames(Lineage1_p8)<-c("SNVID","freq","passage")
Lineage1_p9<-Alternate_L1thruL5[,c(1,11)]
Lineage1_p9$passage<-9
colnames(Lineage1_p9)<-c("SNVID","freq","passage")
Lineage1_p10<-Alternate_L1thruL5[,c(1,12)]
Lineage1_p10$passage<-10
colnames(Lineage1_p10)<-c("SNVID","freq","passage")
Lineage1_p11<-Alternate_L1thruL5[,c(1,13)]
Lineage1_p11$passage<-11
colnames(Lineage1_p11)<-c("SNVID","freq","passage")

Alternate_Lineage1<-rbind(Lineage1_p1,Lineage1_p2B)
Alternate_Lineage1<-rbind(Alternate_Lineage1,Lineage1_p3)
Alternate_Lineage1<-rbind(Alternate_Lineage1,Lineage1_p4)
Alternate_Lineage1<-rbind(Alternate_Lineage1,Lineage1_p5)
Alternate_Lineage1<-rbind(Alternate_Lineage1,Lineage1_p6)
Alternate_Lineage1<-rbind(Alternate_Lineage1,Lineage1_p7)
Alternate_Lineage1<-rbind(Alternate_Lineage1,Lineage1_p8)
Alternate_Lineage1<-rbind(Alternate_Lineage1,Lineage1_p9)
Alternate_Lineage1<-rbind(Alternate_Lineage1,Lineage1_p10)
Alternate_Lineage1<-rbind(Alternate_Lineage1,Lineage1_p11)

Alternate_Lineage1$freq<-as.numeric(Alternate_Lineage1$freq)
Alternate_Lineage1$passage<-as.numeric(Alternate_Lineage1$passage)

Lineage2_p1<-Alternate_L1thruL5[,c(1,14)]
Lineage2_p1$passage<-1
colnames(Lineage2_p1)<-c("SNVID","freq","passage")
Lineage2_p2B<-Alternate_L1thruL5[,c(1,15)]
Lineage2_p2B$passage<-2
colnames(Lineage2_p2B)<-c("SNVID","freq","passage")
Lineage2_p3<-Alternate_L1thruL5[,c(1,17)]
Lineage2_p3$passage<-3
colnames(Lineage2_p3)<-c("SNVID","freq","passage")

Alternate_Lineage2<-rbind(Lineage2_p1,Lineage2_p2B)
Alternate_Lineage2<-rbind(Alternate_Lineage2,Lineage2_p3)

Alternate_Lineage2$freq<-as.numeric(Alternate_Lineage2$freq)
Alternate_Lineage2$passage<-as.numeric(Alternate_Lineage2$passage)

Lineage3_p1<-Alternate_L1thruL5[,c(1,18)]
Lineage3_p1$passage<-1
colnames(Lineage3_p1)<-c("SNVID","freq","passage")
Lineage3_p2B<-Alternate_L1thruL5[,c(1,19)]
Lineage3_p2B$passage<-2
colnames(Lineage3_p2B)<-c("SNVID","freq","passage")
Lineage3_p3<-Alternate_L1thruL5[,c(1,21)]
Lineage3_p3$passage<-3
colnames(Lineage3_p3)<-c("SNVID","freq","passage")

Alternate_Lineage3<-rbind(Lineage3_p1,Lineage3_p2B)
Alternate_Lineage3<-rbind(Alternate_Lineage3,Lineage3_p3)

Alternate_Lineage3$freq<-as.numeric(Alternate_Lineage3$freq)
Alternate_Lineage3$passage<-as.numeric(Alternate_Lineage3$passage)

Lineage4_p1<-Alternate_L1thruL5[,c(1,22)]
Lineage4_p1$passage<-1
colnames(Lineage4_p1)<-c("SNVID","freq","passage")
Lineage4_p2B<-Alternate_L1thruL5[,c(1,23)]
Lineage4_p2B$passage<-2
colnames(Lineage4_p2B)<-c("SNVID","freq","passage")
Lineage4_p3<-Alternate_L1thruL5[,c(1,25)]
Lineage4_p3$passage<-3
colnames(Lineage4_p3)<-c("SNVID","freq","passage")
Lineage4_p4<-Alternate_L1thruL5[,c(1,26)]
Lineage4_p4$passage<-4
colnames(Lineage4_p4)<-c("SNVID","freq","passage")
Lineage4_p5<-Alternate_L1thruL5[,c(1,27)]
Lineage4_p5$passage<-5
colnames(Lineage4_p5)<-c("SNVID","freq","passage")
Lineage4_p6<-Alternate_L1thruL5[,c(1,28)]
Lineage4_p6$passage<-6
colnames(Lineage4_p6)<-c("SNVID","freq","passage")
Lineage4_p7<-Alternate_L1thruL5[,c(1,29)]
Lineage4_p7$passage<-7
colnames(Lineage4_p7)<-c("SNVID","freq","passage")
Lineage4_p8<-Alternate_L1thruL5[,c(1,30)]
Lineage4_p8$passage<-8
colnames(Lineage4_p8)<-c("SNVID","freq","passage")
Lineage4_p9<-Alternate_L1thruL5[,c(1,31)]
Lineage4_p9$passage<-9
colnames(Lineage4_p9)<-c("SNVID","freq","passage")
Lineage4_p10<-Alternate_L1thruL5[,c(1,32)]
Lineage4_p10$passage<-10
colnames(Lineage4_p10)<-c("SNVID","freq","passage")
Lineage4_p11<-Alternate_L1thruL5[,c(1,33)]
Lineage4_p11$passage<-11
colnames(Lineage4_p11)<-c("SNVID","freq","passage")

Lineage4<-rbind(Lineage4_p1,Lineage4_p2B)
Lineage4<-rbind(Lineage4,Lineage4_p3)
Lineage4<-rbind(Lineage4,Lineage4_p4)
Lineage4<-rbind(Lineage4,Lineage4_p5)
Lineage4<-rbind(Lineage4,Lineage4_p6)
Lineage4<-rbind(Lineage4,Lineage4_p7)
Lineage4<-rbind(Lineage4,Lineage4_p8)
Lineage4<-rbind(Lineage4,Lineage4_p9)
Lineage4<-rbind(Lineage4,Lineage4_p10)
Alternate_Lineage4<-rbind(Lineage4,Lineage4_p11)

Alternate_Lineage4$freq<-as.numeric(Alternate_Lineage4$freq)
Alternate_Lineage4$passage<-as.numeric(Alternate_Lineage4$passage)

Lineage5_p1<-Alternate_L1thruL5[,c(1,34)]
Lineage5_p1$passage<-1
colnames(Lineage5_p1)<-c("SNVID","freq","passage")
Lineage5_p2B<-Alternate_L1thruL5[,c(1,35)]
Lineage5_p2B$passage<-2
colnames(Lineage5_p2B)<-c("SNVID","freq","passage")
Lineage5_p3<-Alternate_L1thruL5[,c(1,37)]
Lineage5_p3$passage<-3
colnames(Lineage5_p3)<-c("SNVID","freq","passage")
Lineage5_p4<-Alternate_L1thruL5[,c(1,38)]
Lineage5_p4$passage<-4
colnames(Lineage5_p4)<-c("SNVID","freq","passage")
Lineage5_p5<-Alternate_L1thruL5[,c(1,39)]
Lineage5_p5$passage<-5
colnames(Lineage5_p5)<-c("SNVID","freq","passage")
Lineage5_p6<-Alternate_L1thruL5[,c(1,40)]
Lineage5_p6$passage<-6
colnames(Lineage5_p6)<-c("SNVID","freq","passage")
Lineage5_p7<-Alternate_L1thruL5[,c(1,41)]
Lineage5_p7$passage<-7
colnames(Lineage5_p7)<-c("SNVID","freq","passage")
Lineage5_p8<-Alternate_L1thruL5[,c(1,42)]
Lineage5_p8$passage<-8
colnames(Lineage5_p8)<-c("SNVID","freq","passage")
Lineage5_p9<-Alternate_L1thruL5[,c(1,43)]
Lineage5_p9$passage<-9
colnames(Lineage5_p9)<-c("SNVID","freq","passage")
Lineage5_p10<-Alternate_L1thruL5[,c(1,44)]
Lineage5_p10$passage<-10
colnames(Lineage5_p10)<-c("SNVID","freq","passage")
Lineage5_p11<-Alternate_L1thruL5[,c(1,45)]
Lineage5_p11$passage<-11
colnames(Lineage5_p11)<-c("SNVID","freq","passage")

Lineage5<-rbind(Lineage5_p1,Lineage5_p2B)
Lineage5<-rbind(Lineage5,Lineage5_p3)
Lineage5<-rbind(Lineage5,Lineage5_p4)
Lineage5<-rbind(Lineage5,Lineage5_p5)
Lineage5<-rbind(Lineage5,Lineage5_p6)
Lineage5<-rbind(Lineage5,Lineage5_p7)
Lineage5<-rbind(Lineage5,Lineage5_p8)
Lineage5<-rbind(Lineage5,Lineage5_p9)
Lineage5<-rbind(Lineage5,Lineage5_p10)
Alternate_Lineage5<-rbind(Lineage5,Lineage5_p11)

Alternate_Lineage5$freq<-as.numeric(Alternate_Lineage5$freq)
Alternate_Lineage5$passage<-as.numeric(Alternate_Lineage5$passage)


#### Plot SNVs over time as line plots ####
### Mosquito - As line plot

x<-(rep(c("red","green","blue","yellow","purple","grey","orange","brown","turquoise","pink","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","cyan","darkgrey"),113))

ggplot(data=Mosquito_Lineage1, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linA_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Mosquito_Lineage2, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linB_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Mosquito_Lineage3, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linC_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Mosquito_Lineage4, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linD_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Mosquito_Lineage5, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linE_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


### Mouse - As line plot
x<-(rep(c("red","green","blue","yellow","purple","grey","orange","brown","turquoise","pink","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","cyan","darkgrey"),113))

ggplot(data=Mouse_Lineage1, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linA_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Mouse_Lineage2, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linB_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Mouse_Lineage3, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linC_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Mouse_Lineage4, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linD_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Mouse_Lineage5, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0), limits = c(1,10)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linE_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


### Alternate - As line plot
Alternate_Lineage3<-subset(Alternate_Lineage3,Alternate_Lineage3$passage!=2)
Alternate_Lineage5<-subset(Alternate_Lineage5,Alternate_Lineage5$passage!=2)

x<-(rep(c("red","green","blue","purple","yellow","grey","orange","brown","turquoise","pink","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","cyan","darkgrey"),113))

ggplot(data=Alternate_Lineage1, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1, 10), breaks=seq(1,10,by=1), labels = c("1","2B","3","4B","5","6B","7","8B","9","10B"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/AP_linA_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Alternate_Lineage2, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1, 10), breaks=seq(1,10,by=1), labels = c("1","2B","3","4B","5","6B","7","8B","9","10B"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/AP_linB_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Alternate_Lineage3, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1, 10), breaks=seq(1,10,by=1), labels = c("1","2B","3","4B","5","6B","7","8B","9","10B"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/AP_linC_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Alternate_Lineage4, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1, 10), breaks=seq(1,10,by=1), labels = c("1","2B","3","4B","5","6B","7","8B","9","10B"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/AP_linD_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

ggplot(data=Alternate_Lineage5, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1,10), breaks=seq(1,10,by=1), labels = c("1","2B","3","4B","5","6B","7","8B","9","10B"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/AP_linE_SNVs_overtime_line.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


#### Plot P10 to CC ####
## Mosquitoes ##
ggplot(data=Mosquito_Lineage1, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linA_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Mosquito_Lineage2, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linB_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Mosquito_Lineage3, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linC_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Mosquito_Lineage4, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linD_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Mosquito_Lineage5, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/MP_linE_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)


## Mice ##
ggplot(data=Mouse_Lineage1, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linA_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Mouse_Lineage2, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linB_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Mouse_Lineage3, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linC_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Mouse_Lineage4, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linD_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Mouse_Lineage5, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(10,11,by=1), labels = c("10","CC"), expand = c(0, 0), limits = c(10,11)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/SCP_linE_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)


## Alternate Passages ##
ggplot(data=Alternate_Lineage1, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(10,11), breaks=seq(10,11,by=1), labels = c("10B","CC"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/AP_linA_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Alternate_Lineage4, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(10,11), breaks=seq(10,11,by=1), labels = c("10B","CC"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/AP_linD_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

ggplot(data=Alternate_Lineage5, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=SNVID), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(10,11), breaks=seq(10,11,by=1), labels = c("10B","CC"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("/PATH/TO/AP_linE_SNVs_P10toCC.png",plot=last_plot(),device = png(), width=2, height=6, dpi = 300)

