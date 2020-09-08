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
colnames(p12_L1)<-c("SNVID","freq1","freq2")

p3_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP3_1.vcf", sep="\t", header=FALSE)
p3_L1<-cSplit(p3_L1,"V8",sep=";", type.convert=FALSE)
p3_L1<-cSplit(p3_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L1$SNVID<-paste(p3_L1$V2,p3_L1$V4,sep=";")
p3_L1$SNVID<-paste(p3_L1$SNVID,p3_L1$V5,sep=">")
p3_L1<-p3_L1[,c(12,11)]
colnames(p3_L1)<-c("SNVID","freq")
p3_L1<-subset(p3_L1,p3_L1$freq>=0.003)

p13_L1<-merge(p12_L1, p3_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L1)<-c("SNVID","freq1","freq2","freq3")

p4_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP4_1_B.vcf", sep="\t", header=FALSE)
p4_L1<-cSplit(p4_L1,"V8",sep=";", type.convert=FALSE)
p4_L1<-cSplit(p4_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L1$SNVID<-paste(p4_L1$V2,p4_L1$V4,sep=";")
p4_L1$SNVID<-paste(p4_L1$SNVID,p4_L1$V5,sep=">")
p4_L1<-p4_L1[,c(12,11)]
colnames(p4_L1)<-c("SNVID","freq")
p4_L1<-subset(p4_L1,p4_L1$freq>=0.003)

p14_L1<-merge(p13_L1, p4_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L1)<-c("SNVID","freq1","freq2","freq3","freq4")

p5_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP5_1.vcf", sep="\t", header=FALSE)
p5_L1<-cSplit(p5_L1,"V8",sep=";", type.convert=FALSE)
p5_L1<-cSplit(p5_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p5_L1$SNVID<-paste(p5_L1$V2,p5_L1$V4,sep=";")
p5_L1$SNVID<-paste(p5_L1$SNVID,p5_L1$V5,sep=">")
p5_L1<-p5_L1[,c(12,11)]
colnames(p5_L1)<-c("SNVID","freq")
p5_L1<-subset(p5_L1,p5_L1$freq>=0.003)

p15_L1<-merge(p14_L1, p5_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p15_L1)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p6_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP6_1_B.vcf", sep="\t", header=FALSE)
p6_L1<-cSplit(p6_L1,"V8",sep=";", type.convert=FALSE)
p6_L1<-cSplit(p6_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p6_L1$SNVID<-paste(p6_L1$V2,p6_L1$V4,sep=";")
p6_L1$SNVID<-paste(p6_L1$SNVID,p6_L1$V5,sep=">")
p6_L1<-p6_L1[,c(12,11)]
colnames(p6_L1)<-c("SNVID","freq")
p6_L1<-subset(p6_L1,p6_L1$freq>=0.003)

p16_L1<-merge(p15_L1, p6_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p16_L1)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6")

p7_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP7_1.vcf", sep="\t", header=FALSE)
p7_L1<-cSplit(p7_L1,"V8",sep=";", type.convert=FALSE)
p7_L1<-cSplit(p7_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p7_L1$SNVID<-paste(p7_L1$V2,p7_L1$V4,sep=";")
p7_L1$SNVID<-paste(p7_L1$SNVID,p7_L1$V5,sep=">")
p7_L1<-p7_L1[,c(12,11)]
colnames(p7_L1)<-c("SNVID","freq")
p7_L1<-subset(p7_L1,p7_L1$freq>=0.003)

p17_L1<-merge(p16_L1, p7_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p17_L1)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7")

p8_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP8_1_B.vcf", sep="\t", header=FALSE)
p8_L1<-cSplit(p8_L1,"V8",sep=";", type.convert=FALSE)
p8_L1<-cSplit(p8_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p8_L1$SNVID<-paste(p8_L1$V2,p8_L1$V4,sep=";")
p8_L1$SNVID<-paste(p8_L1$SNVID,p8_L1$V5,sep=">")
p8_L1<-p8_L1[,c(12,11)]
colnames(p8_L1)<-c("SNVID","freq")
p8_L1<-subset(p8_L1,p8_L1$freq>=0.003)

p18_L1<-merge(p17_L1, p8_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p18_L1)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8")

p9_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP9_1.vcf", sep="\t", header=FALSE)
p9_L1<-cSplit(p9_L1,"V8",sep=";", type.convert=FALSE)
p9_L1<-cSplit(p9_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p9_L1$SNVID<-paste(p9_L1$V2,p9_L1$V4,sep=";")
p9_L1$SNVID<-paste(p9_L1$SNVID,p9_L1$V5,sep=">")
p9_L1<-p9_L1[,c(12,11)]
colnames(p9_L1)<-c("SNVID","freq")
p9_L1<-subset(p9_L1,p9_L1$freq>=0.003)

p19_L1<-merge(p18_L1, p9_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p19_L1)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8","freq9")

p10_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP10_1_B.vcf", sep="\t", header=FALSE)
p10_L1<-cSplit(p10_L1,"V8",sep=";", type.convert=FALSE)
p10_L1<-cSplit(p10_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L1$SNVID<-paste(p10_L1$V2,p10_L1$V4,sep=";")
p10_L1$SNVID<-paste(p10_L1$SNVID,p10_L1$V5,sep=">")
p10_L1<-p10_L1[,c(12,11)]
colnames(p10_L1)<-c("SNVID","freq")
p10_L1<-subset(p10_L1,p10_L1$freq>=0.003)

p110_L1<-merge(p19_L1, p10_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L1)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10")

p11_L1<-read.table(file="/PATH/TO/rep1and2_refalignAP10_A.vcf", sep="\t", header=FALSE)
p11_L1<-cSplit(p11_L1,"V8",sep=";", type.convert=FALSE)
p11_L1<-cSplit(p11_L1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L1$SNVID<-paste(p11_L1$V2,p11_L1$V4,sep=";")
p11_L1$SNVID<-paste(p11_L1$SNVID,p11_L1$V5,sep=">")
p11_L1<-p11_L1[,c(12,11)]
colnames(p11_L1)<-c("SNVID","freq")
p11_L1<-subset(p11_L1,p11_L1$freq>=0.003)

p111_L1<-merge(p110_L1, p11_L1, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L1)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10","freqCC")

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

p3_L2<-read.table(file="/PATH/TO/rep1and2_refalignAP3_2.vcf", sep="\t", header=FALSE)
p3_L2<-cSplit(p3_L2,"V8",sep=";", type.convert=FALSE)
p3_L2<-cSplit(p3_L2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L2$SNVID<-paste(p3_L2$V2,p3_L2$V4,sep=";")
p3_L2$SNVID<-paste(p3_L2$SNVID,p3_L2$V5,sep=">")
p3_L2<-p3_L2[,c(12,11)]
colnames(p3_L2)<-c("SNVID","freq")
p3_L2<-subset(p3_L2,p3_L2$freq>=0.003)

p13_L2<-merge(p12_L2, p3_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L2)<-c("SNVID","freq1","freq2","freq3")

p13_L2[is.na(p13_L2)]<-0


p1_L3<-read.table(file="/PATH/TO/rep1and2_refalignAP1_3.vcf", sep="\t", header=FALSE)
p1_L3<-cSplit(p1_L3,"V8",sep=";", type.convert=FALSE)
p1_L3<-cSplit(p1_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L3$SNVID<-paste(p1_L3$V2,p1_L3$V4,sep=";")
p1_L3$SNVID<-paste(p1_L3$SNVID,p1_L3$V5,sep=">")
p1_L3<-p1_L3[,c(12,11)]
colnames(p1_L3)<-c("SNVID","freq")
p1_L3<-subset(p1_L3,p1_L3$freq>=0.003)

p2_L3<-read.table(file="/PATH/TO/rep1and2_refalignAP2_3_5B.vcf", sep="\t", header=FALSE)
p2_L3<-cSplit(p2_L3,"V8",sep=";", type.convert=FALSE)
p2_L3<-cSplit(p2_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2_L3$SNVID<-paste(p2_L3$V2,p2_L3$V4,sep=";")
p2_L3$SNVID<-paste(p2_L3$SNVID,p2_L3$V5,sep=">")
p2_L3<-p2_L3[,c(12,11)]
colnames(p2_L3)<-c("SNVID","freq")
p2_L3<-subset(p2_L3,p2_L3$freq>=0.003)

p12_L3<-merge(p1_L3, p2_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12_L3)<-c("SNVID","freq1","freq2")

p3_L3<-read.table(file="/PATH/TO/rep1and2_refalignAP3_3.vcf", sep="\t", header=FALSE)
p3_L3<-cSplit(p3_L3,"V8",sep=";", type.convert=FALSE)
p3_L3<-cSplit(p3_L3, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L3$SNVID<-paste(p3_L3$V2,p3_L3$V4,sep=";")
p3_L3$SNVID<-paste(p3_L3$SNVID,p3_L3$V5,sep=">")
p3_L3<-p3_L3[,c(12,11)]
colnames(p3_L3)<-c("SNVID","freq")
p3_L3<-subset(p3_L3,p3_L3$freq>=0.003)

p13_L3<-merge(p12_L3, p3_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L3)<-c("SNVID","freq1","freq2","freq3")

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

p3_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP3_4.vcf", sep="\t", header=FALSE)
p3_L4<-cSplit(p3_L4,"V8",sep=";", type.convert=FALSE)
p3_L4<-cSplit(p3_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L4$SNVID<-paste(p3_L4$V2,p3_L4$V4,sep=";")
p3_L4$SNVID<-paste(p3_L4$SNVID,p3_L4$V5,sep=">")
p3_L4<-p3_L4[,c(12,11)]
colnames(p3_L4)<-c("SNVID","freq")
p3_L4<-subset(p3_L4,p3_L4$freq>=0.003)

p13_L4<-merge(p12_L4, p3_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L4)<-c("SNVID","freq1","freq2","freq3")

p4_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP4_4_B.vcf", sep="\t", header=FALSE)
p4_L4<-cSplit(p4_L4,"V8",sep=";", type.convert=FALSE)
p4_L4<-cSplit(p4_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L4$SNVID<-paste(p4_L4$V2,p4_L4$V4,sep=";")
p4_L4$SNVID<-paste(p4_L4$SNVID,p4_L4$V5,sep=">")
p4_L4<-p4_L4[,c(12,11)]
colnames(p4_L4)<-c("SNVID","freq")
p4_L4<-subset(p4_L4,p4_L4$freq>=0.003)

p14_L4<-merge(p13_L4, p4_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L4)<-c("SNVID","freq1","freq2","freq3","freq4")

p5_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP5_4.vcf", sep="\t", header=FALSE)
p5_L4<-cSplit(p5_L4,"V8",sep=";", type.convert=FALSE)
p5_L4<-cSplit(p5_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p5_L4$SNVID<-paste(p5_L4$V2,p5_L4$V4,sep=";")
p5_L4$SNVID<-paste(p5_L4$SNVID,p5_L4$V5,sep=">")
p5_L4<-p5_L4[,c(12,11)]
colnames(p5_L4)<-c("SNVID","freq")
p5_L4<-subset(p5_L4,p5_L4$freq>=0.003)

p15_L4<-merge(p14_L4, p5_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p15_L4)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p6_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP6_4_B.vcf", sep="\t", header=FALSE)
p6_L4<-cSplit(p6_L4,"V8",sep=";", type.convert=FALSE)
p6_L4<-cSplit(p6_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p6_L4$SNVID<-paste(p6_L4$V2,p6_L4$V4,sep=";")
p6_L4$SNVID<-paste(p6_L4$SNVID,p6_L4$V5,sep=">")
p6_L4<-p6_L4[,c(12,11)]
colnames(p6_L4)<-c("SNVID","freq")
p6_L4<-subset(p6_L4,p6_L4$freq>=0.003)

p16_L4<-merge(p15_L4, p6_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p16_L4)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6")

p7_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP7_4.vcf", sep="\t", header=FALSE)
p7_L4<-cSplit(p7_L4,"V8",sep=";", type.convert=FALSE)
p7_L4<-cSplit(p7_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p7_L4$SNVID<-paste(p7_L4$V2,p7_L4$V4,sep=";")
p7_L4$SNVID<-paste(p7_L4$SNVID,p7_L4$V5,sep=">")
p7_L4<-p7_L4[,c(12,11)]
colnames(p7_L4)<-c("SNVID","freq")
p7_L4<-subset(p7_L4,p7_L4$freq>=0.003)

p17_L4<-merge(p16_L4, p7_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p17_L4)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7")

p8_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP8_4_B.vcf", sep="\t", header=FALSE)
p8_L4<-cSplit(p8_L4,"V8",sep=";", type.convert=FALSE)
p8_L4<-cSplit(p8_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p8_L4$SNVID<-paste(p8_L4$V2,p8_L4$V4,sep=";")
p8_L4$SNVID<-paste(p8_L4$SNVID,p8_L4$V5,sep=">")
p8_L4<-p8_L4[,c(12,11)]
colnames(p8_L4)<-c("SNVID","freq")
p8_L4<-subset(p8_L4,p8_L4$freq>=0.003)

p18_L4<-merge(p17_L4, p8_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p18_L4)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8")

p9_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP9_4.vcf", sep="\t", header=FALSE)
p9_L4<-cSplit(p9_L4,"V8",sep=";", type.convert=FALSE)
p9_L4<-cSplit(p9_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p9_L4$SNVID<-paste(p9_L4$V2,p9_L4$V4,sep=";")
p9_L4$SNVID<-paste(p9_L4$SNVID,p9_L4$V5,sep=">")
p9_L4<-p9_L4[,c(12,11)]
colnames(p9_L4)<-c("SNVID","freq")
p9_L4<-subset(p9_L4,p9_L4$freq>=0.003)

p19_L4<-merge(p18_L4, p9_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p19_L4)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8","freq9")

p10_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP10_4_B.vcf", sep="\t", header=FALSE)
p10_L4<-cSplit(p10_L4,"V8",sep=";", type.convert=FALSE)
p10_L4<-cSplit(p10_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L4$SNVID<-paste(p10_L4$V2,p10_L4$V4,sep=";")
p10_L4$SNVID<-paste(p10_L4$SNVID,p10_L4$V5,sep=">")
p10_L4<-p10_L4[,c(12,11)]
colnames(p10_L4)<-c("SNVID","freq")
p10_L4<-subset(p10_L4,p10_L4$freq>=0.003)

p110_L4<-merge(p19_L4, p10_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L4)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10")

p11_L4<-read.table(file="/PATH/TO/rep1and2_refalignAP10_D.vcf", sep="\t", header=FALSE)
p11_L4<-cSplit(p11_L4,"V8",sep=";", type.convert=FALSE)
p11_L4<-cSplit(p11_L4, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L4$SNVID<-paste(p11_L4$V2,p11_L4$V4,sep=";")
p11_L4$SNVID<-paste(p11_L4$SNVID,p11_L4$V5,sep=">")
p11_L4<-p11_L4[,c(12,11)]
colnames(p11_L4)<-c("SNVID","freq")
p11_L4<-subset(p11_L4,p11_L4$freq>=0.003)

p111_L4<-merge(p110_L4, p11_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L4)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10","freqCC")

p111_L4[is.na(p111_L4)]<-0


p1_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP1_5.vcf", sep="\t", header=FALSE)
p1_L5<-cSplit(p1_L5,"V8",sep=";", type.convert=FALSE)
p1_L5<-cSplit(p1_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p1_L5$SNVID<-paste(p1_L5$V2,p1_L5$V4,sep=";")
p1_L5$SNVID<-paste(p1_L5$SNVID,p1_L5$V5,sep=">")
p1_L5<-p1_L5[,c(12,11)]
colnames(p1_L5)<-c("SNVID","freq")
p1_L5<-subset(p1_L5,p1_L5$freq>=0.003)

p2_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP2_5_3B.vcf", sep="\t", header=FALSE)
p2_L5<-cSplit(p2_L5,"V8",sep=";", type.convert=FALSE)
p2_L5<-cSplit(p2_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p2_L5$SNVID<-paste(p2_L5$V2,p2_L5$V4,sep=";")
p2_L5$SNVID<-paste(p2_L5$SNVID,p2_L5$V5,sep=">")
p2_L5<-p2_L5[,c(12,11)]
colnames(p2_L5)<-c("SNVID","freq")
p2_L5<-subset(p2_L5,p2_L5$freq>=0.003)

p12_L5<-merge(p1_L5, p2_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p12_L5)<-c("SNVID","freq1","freq2")

p3_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP3_5.vcf", sep="\t", header=FALSE)
p3_L5<-cSplit(p3_L5,"V8",sep=";", type.convert=FALSE)
p3_L5<-cSplit(p3_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p3_L5$SNVID<-paste(p3_L5$V2,p3_L5$V4,sep=";")
p3_L5$SNVID<-paste(p3_L5$SNVID,p3_L5$V5,sep=">")
p3_L5<-p3_L5[,c(12,11)]
colnames(p3_L5)<-c("SNVID","freq")
p3_L5<-subset(p3_L5,p3_L5$freq>=0.003)

p13_L5<-merge(p12_L5, p3_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p13_L5)<-c("SNVID","freq1","freq2","freq3")

p4_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP4_5_B.vcf", sep="\t", header=FALSE)
p4_L5<-cSplit(p4_L5,"V8",sep=";", type.convert=FALSE)
p4_L5<-cSplit(p4_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p4_L5$SNVID<-paste(p4_L5$V2,p4_L5$V4,sep=";")
p4_L5$SNVID<-paste(p4_L5$SNVID,p4_L5$V5,sep=">")
p4_L5<-p4_L5[,c(12,11)]
colnames(p4_L5)<-c("SNVID","freq")
p4_L5<-subset(p4_L5,p4_L5$freq>=0.003)

p14_L5<-merge(p13_L5, p4_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p14_L5)<-c("SNVID","freq1","freq2","freq3","freq4")

p5_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP5_5.vcf", sep="\t", header=FALSE)
p5_L5<-cSplit(p5_L5,"V8",sep=";", type.convert=FALSE)
p5_L5<-cSplit(p5_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p5_L5$SNVID<-paste(p5_L5$V2,p5_L5$V4,sep=";")
p5_L5$SNVID<-paste(p5_L5$SNVID,p5_L5$V5,sep=">")
p5_L5<-p5_L5[,c(12,11)]
colnames(p5_L5)<-c("SNVID","freq")
p5_L5<-subset(p5_L5,p5_L5$freq>=0.003)

p15_L5<-merge(p14_L5, p5_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p15_L5)<-c("SNVID","freq1","freq2","freq3","freq4","freq5")

p6_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP6_5_B.vcf", sep="\t", header=FALSE)
p6_L5<-cSplit(p6_L5,"V8",sep=";", type.convert=FALSE)
p6_L5<-cSplit(p6_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p6_L5$SNVID<-paste(p6_L5$V2,p6_L5$V4,sep=";")
p6_L5$SNVID<-paste(p6_L5$SNVID,p6_L5$V5,sep=">")
p6_L5<-p6_L5[,c(12,11)]
colnames(p6_L5)<-c("SNVID","freq")
p6_L5<-subset(p6_L5,p6_L5$freq>=0.003)

p16_L5<-merge(p15_L5, p6_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p16_L5)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6")

p7_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP7_5.vcf", sep="\t", header=FALSE)
p7_L5<-cSplit(p7_L5,"V8",sep=";", type.convert=FALSE)
p7_L5<-cSplit(p7_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p7_L5$SNVID<-paste(p7_L5$V2,p7_L5$V4,sep=";")
p7_L5$SNVID<-paste(p7_L5$SNVID,p7_L5$V5,sep=">")
p7_L5<-p7_L5[,c(12,11)]
colnames(p7_L5)<-c("SNVID","freq")
p7_L5<-subset(p7_L5,p7_L5$freq>=0.003)

p17_L5<-merge(p16_L5, p7_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p17_L5)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7")

p8_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP8_5_B.vcf", sep="\t", header=FALSE)
p8_L5<-cSplit(p8_L5,"V8",sep=";", type.convert=FALSE)
p8_L5<-cSplit(p8_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p8_L5$SNVID<-paste(p8_L5$V2,p8_L5$V4,sep=";")
p8_L5$SNVID<-paste(p8_L5$SNVID,p8_L5$V5,sep=">")
p8_L5<-p8_L5[,c(12,11)]
colnames(p8_L5)<-c("SNVID","freq")
p8_L5<-subset(p8_L5,p8_L5$freq>=0.003)

p18_L5<-merge(p17_L5, p8_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p18_L5)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8")

p9_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP9_5.vcf", sep="\t", header=FALSE)
p9_L5<-cSplit(p9_L5,"V8",sep=";", type.convert=FALSE)
p9_L5<-cSplit(p9_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p9_L5$SNVID<-paste(p9_L5$V2,p9_L5$V4,sep=";")
p9_L5$SNVID<-paste(p9_L5$SNVID,p9_L5$V5,sep=">")
p9_L5<-p9_L5[,c(12,11)]
colnames(p9_L5)<-c("SNVID","freq")
p9_L5<-subset(p9_L5,p9_L5$freq>=0.003)

p19_L5<-merge(p18_L5, p9_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p19_L5)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8","freq9")

p10_L5<-read.table(file="/PATH/TO/rep1and2_refalignAP10_5_B.vcf", sep="\t", header=FALSE)
p10_L5<-cSplit(p10_L5,"V8",sep=";", type.convert=FALSE)
p10_L5<-cSplit(p10_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p10_L5$SNVID<-paste(p10_L5$V2,p10_L5$V4,sep=";")
p10_L5$SNVID<-paste(p10_L5$SNVID,p10_L5$V5,sep=">")
p10_L5<-p10_L5[,c(12,11)]
colnames(p10_L5)<-c("SNVID","freq")
p10_L5<-subset(p10_L5,p10_L5$freq>=0.003)

p110_L5<-merge(p19_L5, p10_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p110_L5)<-c("SNVID","freq1","freq2","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10")

p110_L5[is.na(p110_L5)]<-0

#p11_L5<-read.table(file="/Volumes/HD1/R21_ZIKVConspecificPassage/Alternate_Passage/WGS/AP10_E_rep1_v3.1/reference_aligned/lofreq_reference_AP10_E_1.vcf", sep="\t", header=FALSE)
#p11_L5<-cSplit(p11_L5,"V8",sep=";", type.convert=FALSE)
#p11_L5<-cSplit(p11_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
#p11_L5$SNVID<-paste(p11_L5$V2,p11_L5$V4,sep=";")
#p11_L5$SNVID<-paste(p11_L5$SNVID,p11_L5$V5,sep=">")
#p11_L5<-p11_L5[,c(12,11)]
#colnames(p11_L5)<-c("SNVID","freq")
#p11_L5<-subset(p11_L5,p11_L5$freq>=0.003)

#p111_L5<-merge(p110_L5, p11_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
#colnames(p111_L5)<-c("SNVID","freq1","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10","freqCC")

#p111_L5[is.na(p111_L5)]<-0


L1andL2<-merge(p111_L1, p13_L2, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1andL2)<-c("SNVID","p1_L1","p2_L1","p3_L1","p4_L1","p5_L1","p6_L1","p7_L1","p8_L1","p9_L1","p10_L1","p11_L1","p1_L2","p2_L2","p3_L2")
L1thruL3<-merge(L1andL2,p13_L3, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1thruL3)<-c("SNVID","p1_L1","p2_L1","p3_L1","p4_L1","p5_L1","p6_L1","p7_L1","p8_L1","p9_L1","p10_L1","p11_L1","p1_L2","p2_L2","p3_L2","p1_L3","p2_L3","p3_L3")
L1thruL4<-merge(L1thruL3,p111_L4, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(L1thruL4)<-c("SNVID","p1_L1","p2_L1","p3_L1","p4_L1","p5_L1","p6_L1","p7_L1","p8_L1","p9_L1","p10_L1","p11_L1","p1_L2","p2_L2","p3_L2","p1_L3","p2_L3","p3_L3","p1_L4","p2_L4","p3_L4","p4_L4","p5_L4","p6_L4","p7_L4","p8_L4","p9_L4","p10_L4","p11_L4")
Alternate_L1thruL5<-merge(L1thruL4,p110_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(Alternate_L1thruL5)<-c("SNVID","p1_L1","p2_L1","p3_L1","p4_L1","p5_L1","p6_L1","p7_L1","p8_L1","p9_L1","p10_L1","p11_L1","p1_L2","p2_L2","p3_L2","p1_L3","p2_L3","p3_L3","p1_L4","p2_L4","p3_L4","p4_L4","p5_L4","p6_L4","p7_L4","p8_L4","p9_L4","p10_L4","p11_L4","p1_L5","p2_L5","p3_L5","p4_L5","p5_L5","p6_L5","p7_L5","p8_L5","p9_L5","p10_L5")

Alternate_L1thruL5[is.na(Alternate_L1thruL5)]<-0

Lineage1_p1<-Alternate_L1thruL5[,c(1,2)]
Lineage1_p1$passage<-1
colnames(Lineage1_p1)<-c("SNVID","freq","passage")
Lineage1_p2<-Alternate_L1thruL5[,c(1,3)]
Lineage1_p2$passage<-2
colnames(Lineage1_p2)<-c("SNVID","freq","passage")
Lineage1_p3<-Alternate_L1thruL5[,c(1,4)]
Lineage1_p3$passage<-3
colnames(Lineage1_p3)<-c("SNVID","freq","passage")
Lineage1_p4<-Alternate_L1thruL5[,c(1,5)]
Lineage1_p4$passage<-4
colnames(Lineage1_p4)<-c("SNVID","freq","passage")
Lineage1_p5<-Alternate_L1thruL5[,c(1,6)]
Lineage1_p5$passage<-5
colnames(Lineage1_p5)<-c("SNVID","freq","passage")
Lineage1_p6<-Alternate_L1thruL5[,c(1,7)]
Lineage1_p6$passage<-6
colnames(Lineage1_p6)<-c("SNVID","freq","passage")
Lineage1_p7<-Alternate_L1thruL5[,c(1,8)]
Lineage1_p7$passage<-7
colnames(Lineage1_p7)<-c("SNVID","freq","passage")
Lineage1_p8<-Alternate_L1thruL5[,c(1,9)]
Lineage1_p8$passage<-8
colnames(Lineage1_p8)<-c("SNVID","freq","passage")
Lineage1_p9<-Alternate_L1thruL5[,c(1,10)]
Lineage1_p9$passage<-9
colnames(Lineage1_p9)<-c("SNVID","freq","passage")
Lineage1_p10<-Alternate_L1thruL5[,c(1,11)]
Lineage1_p10$passage<-10
colnames(Lineage1_p10)<-c("SNVID","freq","passage")
Lineage1_p11<-Alternate_L1thruL5[,c(1,12)]
Lineage1_p11$passage<-11
colnames(Lineage1_p11)<-c("SNVID","freq","passage")

Alternate_Lineage1<-rbind(Lineage1_p1,Lineage1_p2)
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

Lineage2_p1<-Alternate_L1thruL5[,c(1,13)]
Lineage2_p1$passage<-1
colnames(Lineage2_p1)<-c("SNVID","freq","passage")
Lineage2_p2<-Alternate_L1thruL5[,c(1,14)]
Lineage2_p2$passage<-2
colnames(Lineage2_p2)<-c("SNVID","freq","passage")
Lineage2_p3<-Alternate_L1thruL5[,c(1,15)]
Lineage2_p3$passage<-3
colnames(Lineage2_p3)<-c("SNVID","freq","passage")

Alternate_Lineage2<-rbind(Lineage2_p1,Lineage2_p2)
Alternate_Lineage2<-rbind(Alternate_Lineage2,Lineage2_p3)

Alternate_Lineage2$freq<-as.numeric(Alternate_Lineage2$freq)
Alternate_Lineage2$passage<-as.numeric(Alternate_Lineage2$passage)

Lineage3_p1<-Alternate_L1thruL5[,c(1,16)]
Lineage3_p1$passage<-1
colnames(Lineage3_p1)<-c("SNVID","freq","passage")
Lineage3_p2<-Alternate_L1thruL5[,c(1,17)]
Lineage3_p2$passage<-2
colnames(Lineage3_p2)<-c("SNVID","freq","passage")
Lineage3_p3<-Alternate_L1thruL5[,c(1,18)]
Lineage3_p3$passage<-3
colnames(Lineage3_p3)<-c("SNVID","freq","passage")

Alternate_Lineage3<-rbind(Lineage3_p1,Lineage3_p2)
Alternate_Lineage3<-rbind(Alternate_Lineage3,Lineage3_p3)

Alternate_Lineage3$freq<-as.numeric(Alternate_Lineage3$freq)
Alternate_Lineage3$passage<-as.numeric(Alternate_Lineage3$passage)

Lineage4_p1<-Alternate_L1thruL5[,c(1,19)]
Lineage4_p1$passage<-1
colnames(Lineage4_p1)<-c("SNVID","freq","passage")
Lineage4_p2<-Alternate_L1thruL5[,c(1,20)]
Lineage4_p2$passage<-2
colnames(Lineage4_p2)<-c("SNVID","freq","passage")
Lineage4_p3<-Alternate_L1thruL5[,c(1,21)]
Lineage4_p3$passage<-3
colnames(Lineage4_p3)<-c("SNVID","freq","passage")
Lineage4_p4<-Alternate_L1thruL5[,c(1,22)]
Lineage4_p4$passage<-4
colnames(Lineage4_p4)<-c("SNVID","freq","passage")
Lineage4_p5<-Alternate_L1thruL5[,c(1,23)]
Lineage4_p5$passage<-5
colnames(Lineage4_p5)<-c("SNVID","freq","passage")
Lineage4_p6<-Alternate_L1thruL5[,c(1,24)]
Lineage4_p6$passage<-6
colnames(Lineage4_p6)<-c("SNVID","freq","passage")
Lineage4_p7<-Alternate_L1thruL5[,c(1,25)]
Lineage4_p7$passage<-7
colnames(Lineage4_p7)<-c("SNVID","freq","passage")
Lineage4_p8<-Alternate_L1thruL5[,c(1,26)]
Lineage4_p8$passage<-8
colnames(Lineage4_p8)<-c("SNVID","freq","passage")
Lineage4_p9<-Alternate_L1thruL5[,c(1,27)]
Lineage4_p9$passage<-9
colnames(Lineage4_p9)<-c("SNVID","freq","passage")
Lineage4_p10<-Alternate_L1thruL5[,c(1,28)]
Lineage4_p10$passage<-10
colnames(Lineage4_p10)<-c("SNVID","freq","passage")
Lineage4_p11<-Alternate_L1thruL5[,c(1,29)]
Lineage4_p11$passage<-11
colnames(Lineage4_p11)<-c("SNVID","freq","passage")

Lineage4<-rbind(Lineage4_p1,Lineage4_p2)
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

Lineage5_p1<-Alternate_L1thruL5[,c(1,30)]
Lineage5_p1$passage<-1
colnames(Lineage5_p1)<-c("SNVID","freq","passage")
Lineage5_p2<-Alternate_L1thruL5[,c(1,31)]
Lineage5_p2$passage<-2
colnames(Lineage5_p2)<-c("SNVID","freq","passage")
Lineage5_p3<-Alternate_L1thruL5[,c(1,32)]
Lineage5_p3$passage<-3
colnames(Lineage5_p3)<-c("SNVID","freq","passage")
Lineage5_p4<-Alternate_L1thruL5[,c(1,33)]
Lineage5_p4$passage<-4
colnames(Lineage5_p4)<-c("SNVID","freq","passage")
Lineage5_p5<-Alternate_L1thruL5[,c(1,34)]
Lineage5_p5$passage<-5
colnames(Lineage5_p5)<-c("SNVID","freq","passage")
Lineage5_p6<-Alternate_L1thruL5[,c(1,35)]
Lineage5_p6$passage<-6
colnames(Lineage5_p6)<-c("SNVID","freq","passage")
Lineage5_p7<-Alternate_L1thruL5[,c(1,36)]
Lineage5_p7$passage<-7
colnames(Lineage5_p7)<-c("SNVID","freq","passage")
Lineage5_p8<-Alternate_L1thruL5[,c(1,37)]
Lineage5_p8$passage<-8
colnames(Lineage5_p8)<-c("SNVID","freq","passage")
Lineage5_p9<-Alternate_L1thruL5[,c(1,38)]
Lineage5_p9$passage<-9
colnames(Lineage5_p9)<-c("SNVID","freq","passage")
Lineage5_p10<-Alternate_L1thruL5[,c(1,39)]
Lineage5_p10$passage<-10
colnames(Lineage5_p10)<-c("SNVID","freq","passage")
#Lineage5_p11<-Alternate_L1thruL5[,c(1,40)]
#Lineage5_p11$passage<-11
#colnames(Lineage5_p11)<-c("SNVID","freq","passage")

Lineage5<-rbind(Lineage5_p1,Lineage5_p2)
Lineage5<-rbind(Lineage5,Lineage5_p3)
Lineage5<-rbind(Lineage5,Lineage5_p4)
Lineage5<-rbind(Lineage5,Lineage5_p5)
Lineage5<-rbind(Lineage5,Lineage5_p6)
Lineage5<-rbind(Lineage5,Lineage5_p7)
Lineage5<-rbind(Lineage5,Lineage5_p8)
Lineage5<-rbind(Lineage5,Lineage5_p9)
Alternate_Lineage5<-rbind(Lineage5,Lineage5_p10)
#Alternate_Lineage5<-rbind(Lineage5,Lineage5_p11)

Alternate_Lineage5$freq<-as.numeric(Alternate_Lineage5$freq)
Alternate_Lineage5$passage<-as.numeric(Alternate_Lineage5$passage)

#### Shared SNVs across Mosquito lineages ####
Mosquito_L1thruL5_1<-Mosquito_L1thruL5[,c(1:5,7:10,12:15,17:20,22:25)]
Mouse_L1thruL5_1<-Mouse_L1thruL5[,c(1:5,7:10,12:15,17:20,22:25)]
Alternate_L1thruL5_1<-Alternate_L1thruL5[,c(1:11,13:28,30:39)]

Mosquito_SNVs_2<-Mosquito_L1thruL5_1 %>% filter_at(vars(starts_with("p")),any_vars(.>0.01))
Mouse_SNVs_2<-Mouse_L1thruL5_1 %>% filter_at(vars(starts_with("p")),any_vars(.>0.01))
Alternate_SNVs_2<-Alternate_L1thruL5_1 %>% filter_at(vars(starts_with("p")),any_vars(.>0.01))

Mosquito_SNVs_2$p1_L1<-as.numeric(Mosquito_SNVs_2$p1_L1)
Mosquito_SNVs_2$p3_L1<-as.numeric(Mosquito_SNVs_2$p3_L1)
Mosquito_SNVs_2$p4_L1<-as.numeric(Mosquito_SNVs_2$p4_L1)
Mosquito_SNVs_2$p10_L1<-as.numeric(Mosquito_SNVs_2$p10_L1)
Mosquito_SNVs_2$p1_L2<-as.numeric(Mosquito_SNVs_2$p1_L2)
Mosquito_SNVs_2$p3_L2<-as.numeric(Mosquito_SNVs_2$p3_L2)
Mosquito_SNVs_2$p4_L2<-as.numeric(Mosquito_SNVs_2$p4_L2)
Mosquito_SNVs_2$p10_L2<-as.numeric(Mosquito_SNVs_2$p10_L2)
Mosquito_SNVs_2$p1_L3<-as.numeric(Mosquito_SNVs_2$p1_L3)
Mosquito_SNVs_2$p3_L3<-as.numeric(Mosquito_SNVs_2$p3_L3)
Mosquito_SNVs_2$p4_L3<-as.numeric(Mosquito_SNVs_2$p4_L3)
Mosquito_SNVs_2$p10_L3<-as.numeric(Mosquito_SNVs_2$p10_L3)
Mosquito_SNVs_2$p1_L4<-as.numeric(Mosquito_SNVs_2$p1_L4)
Mosquito_SNVs_2$p3_L4<-as.numeric(Mosquito_SNVs_2$p3_L4)
Mosquito_SNVs_2$p4_L4<-as.numeric(Mosquito_SNVs_2$p4_L4)
Mosquito_SNVs_2$p10_L4<-as.numeric(Mosquito_SNVs_2$p10_L4)
Mosquito_SNVs_2$p1_L5<-as.numeric(Mosquito_SNVs_2$p1_L5)
Mosquito_SNVs_2$p3_L5<-as.numeric(Mosquito_SNVs_2$p3_L5)
Mosquito_SNVs_2$p4_L5<-as.numeric(Mosquito_SNVs_2$p4_L5)
Mosquito_SNVs_2$p10_L5<-as.numeric(Mosquito_SNVs_2$p10_L5)

## Mosquito Lineage 1_SNVs shared between mosquito lineages ##
Mosquito_Lineage1_only<-merge(Mosquito_Lineage1, Mosquito_SNVs_2, by = "SNVID")
Mosquito_Lineage1_only$lin1freqsum<-Mosquito_Lineage1_only$p1_L1 + Mosquito_Lineage1_only$p3_L1 + Mosquito_Lineage1_only$p4_L1 + Mosquito_Lineage1_only$p10_L1
Mosquito_Lineage1_only$alllinfreqsum<-Mosquito_Lineage1_only$p1_L1 + Mosquito_Lineage1_only$p3_L1 + Mosquito_Lineage1_only$p4_L1 + Mosquito_Lineage1_only$p10_L1 + 
  Mosquito_Lineage1_only$p1_L2 + Mosquito_Lineage1_only$p3_L2 + Mosquito_Lineage1_only$p4_L2 + Mosquito_Lineage1_only$p10_L2 + 
  Mosquito_Lineage1_only$p1_L3 + Mosquito_Lineage1_only$p3_L3 + Mosquito_Lineage1_only$p4_L3 + Mosquito_Lineage1_only$p10_L3 + 
  Mosquito_Lineage1_only$p1_L4 + Mosquito_Lineage1_only$p3_L4 + Mosquito_Lineage1_only$p4_L4 + Mosquito_Lineage1_only$p10_L4 + 
  Mosquito_Lineage1_only$p1_L5 + Mosquito_Lineage1_only$p3_L5 + Mosquito_Lineage1_only$p4_L5 + Mosquito_Lineage1_only$p10_L5
Mosquito_Lineage1_only$lin1only<-ifelse(Mosquito_Lineage1_only$lin1freqsum==Mosquito_Lineage1_only$alllinfreqsum,1,0)

Mosquito_Lineage1_only %>% count(lin1only)

ggplot(data=Mosquito_Lineage1_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin1only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/MP_linA_SNVs_overtime_line_MosquitoSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mosquito Lineage 2_SNVs shared between mosquito lineages ##
Mosquito_Lineage2_only<-merge(Mosquito_Lineage2, Mosquito_SNVs_2, by = "SNVID")
Mosquito_Lineage2_only$lin2freqsum<-Mosquito_Lineage2_only$p1_L2 + Mosquito_Lineage2_only$p3_L2 + Mosquito_Lineage2_only$p4_L2 + Mosquito_Lineage2_only$p10_L2
Mosquito_Lineage2_only$alllinfreqsum<-Mosquito_Lineage2_only$p1_L1 + Mosquito_Lineage2_only$p3_L1 + Mosquito_Lineage2_only$p4_L1 + Mosquito_Lineage2_only$p10_L1 + 
  Mosquito_Lineage2_only$p1_L2 + Mosquito_Lineage2_only$p3_L2 + Mosquito_Lineage2_only$p4_L2 + Mosquito_Lineage2_only$p10_L2 + 
  Mosquito_Lineage2_only$p1_L3 + Mosquito_Lineage2_only$p3_L3 + Mosquito_Lineage2_only$p4_L3 + Mosquito_Lineage2_only$p10_L3 + 
  Mosquito_Lineage2_only$p1_L4 + Mosquito_Lineage2_only$p3_L4 + Mosquito_Lineage2_only$p4_L4 + Mosquito_Lineage2_only$p10_L4 + 
  Mosquito_Lineage2_only$p1_L5 + Mosquito_Lineage2_only$p3_L5 + Mosquito_Lineage2_only$p4_L5 + Mosquito_Lineage2_only$p10_L5
Mosquito_Lineage2_only$lin2only<-ifelse(Mosquito_Lineage2_only$lin2freqsum==Mosquito_Lineage2_only$alllinfreqsum,1,0)

ggplot(data=Mosquito_Lineage2_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin2only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/MP_linB_SNVs_overtime_line_MosquitoSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mosquito Lineage 3_SNVs shared between mosquito lineages ##
Mosquito_Lineage3_only<-merge(Mosquito_Lineage3, Mosquito_SNVs_2, by = "SNVID")
Mosquito_Lineage3_only$lin3freqsum<-Mosquito_Lineage3_only$p1_L3 + Mosquito_Lineage3_only$p3_L3 + Mosquito_Lineage3_only$p4_L3 + Mosquito_Lineage3_only$p10_L3
Mosquito_Lineage3_only$alllinfreqsum<-Mosquito_Lineage3_only$p1_L1 + Mosquito_Lineage3_only$p3_L1 + Mosquito_Lineage3_only$p4_L1 + Mosquito_Lineage3_only$p10_L1 + 
  Mosquito_Lineage3_only$p1_L2 + Mosquito_Lineage3_only$p3_L2 + Mosquito_Lineage3_only$p4_L2 + Mosquito_Lineage3_only$p10_L2 + 
  Mosquito_Lineage3_only$p1_L3 + Mosquito_Lineage3_only$p3_L3 + Mosquito_Lineage3_only$p4_L3 + Mosquito_Lineage3_only$p10_L3 + 
  Mosquito_Lineage3_only$p1_L4 + Mosquito_Lineage3_only$p3_L4 + Mosquito_Lineage3_only$p4_L4 + Mosquito_Lineage3_only$p10_L4 + 
  Mosquito_Lineage3_only$p1_L5 + Mosquito_Lineage3_only$p3_L5 + Mosquito_Lineage3_only$p4_L5 + Mosquito_Lineage3_only$p10_L5
Mosquito_Lineage3_only$lin3only<-ifelse(Mosquito_Lineage3_only$lin3freqsum==Mosquito_Lineage3_only$alllinfreqsum,1,0)

ggplot(data=Mosquito_Lineage3_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin3only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/MP_linC_SNVs_overtime_line_MosquitoSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mosquito Lineage 4_SNVs shared between mosquito lineages ##
Mosquito_Lineage4_only<-merge(Mosquito_Lineage4, Mosquito_SNVs_2, by = "SNVID")
Mosquito_Lineage4_only$lin4freqsum<-Mosquito_Lineage4_only$p1_L4 + Mosquito_Lineage4_only$p3_L4 + Mosquito_Lineage4_only$p4_L4 + Mosquito_Lineage4_only$p10_L4
Mosquito_Lineage4_only$alllinfreqsum<-Mosquito_Lineage4_only$p1_L1 + Mosquito_Lineage4_only$p3_L1 + Mosquito_Lineage4_only$p4_L1 + Mosquito_Lineage4_only$p10_L1 + 
  Mosquito_Lineage4_only$p1_L2 + Mosquito_Lineage4_only$p3_L2 + Mosquito_Lineage4_only$p4_L2 + Mosquito_Lineage4_only$p10_L2 + 
  Mosquito_Lineage4_only$p1_L3 + Mosquito_Lineage4_only$p3_L3 + Mosquito_Lineage4_only$p4_L3 + Mosquito_Lineage4_only$p10_L3 + 
  Mosquito_Lineage4_only$p1_L4 + Mosquito_Lineage4_only$p3_L4 + Mosquito_Lineage4_only$p4_L4 + Mosquito_Lineage4_only$p10_L4 + 
  Mosquito_Lineage4_only$p1_L5 + Mosquito_Lineage4_only$p3_L5 + Mosquito_Lineage4_only$p4_L5 + Mosquito_Lineage4_only$p10_L5
Mosquito_Lineage4_only$lin4only<-ifelse(Mosquito_Lineage4_only$lin4freqsum==Mosquito_Lineage4_only$alllinfreqsum,1,0)

ggplot(data=Mosquito_Lineage4_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin4only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/MP_linD_SNVs_overtime_line_MosquitoSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mosquito Lineage 5_SNVs shared between mosquito lineages ##
Mosquito_Lineage5_only<-merge(Mosquito_Lineage5, Mosquito_SNVs_2, by = "SNVID")
Mosquito_Lineage5_only$lin5freqsum<-Mosquito_Lineage5_only$p1_L5 + Mosquito_Lineage5_only$p3_L5 + Mosquito_Lineage5_only$p4_L5 + Mosquito_Lineage5_only$p10_L5
Mosquito_Lineage5_only$alllinfreqsum<-Mosquito_Lineage5_only$p1_L1 + Mosquito_Lineage5_only$p3_L1 + Mosquito_Lineage5_only$p4_L1 + Mosquito_Lineage5_only$p10_L1 + 
  Mosquito_Lineage5_only$p1_L2 + Mosquito_Lineage5_only$p3_L2 + Mosquito_Lineage5_only$p4_L2 + Mosquito_Lineage5_only$p10_L2 + 
  Mosquito_Lineage5_only$p1_L3 + Mosquito_Lineage5_only$p3_L3 + Mosquito_Lineage5_only$p4_L3 + Mosquito_Lineage5_only$p10_L3 + 
  Mosquito_Lineage5_only$p1_L4 + Mosquito_Lineage5_only$p3_L4 + Mosquito_Lineage5_only$p4_L4 + Mosquito_Lineage5_only$p10_L4 + 
  Mosquito_Lineage5_only$p1_L5 + Mosquito_Lineage5_only$p3_L5 + Mosquito_Lineage5_only$p4_L5 + Mosquito_Lineage5_only$p10_L5
Mosquito_Lineage5_only$lin5only<-ifelse(Mosquito_Lineage5_only$lin5freqsum==Mosquito_Lineage5_only$alllinfreqsum,1,0)

ggplot(data=Mosquito_Lineage5_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin5only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/MP_linE_SNVs_overtime_line_MosquitoSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Concatenate all lineage only columns
linonly<-as.data.frame(cbind(Mosquito_Lineage1_only$lin1only,Mosquito_Lineage2_only$lin2only,Mosquito_Lineage3_only$lin3only,Mosquito_Lineage4_only$lin4only,Mosquito_Lineage5_only$lin5only))
linonly$combined<-(linonly$V1+linonly$V2+linonly$V3+linonly$V4+linonly$V5)
linonly$onelinsnv<-ifelse(linonly$combined>0,0,1)
linonly$morethanonelinsnv<-ifelse(linonly$combined>0,1,0)

onelinSNV<-sum(linonly$onelinsnv)/5
morethanonelinSNV<-sum(linonly$morethanonelinsnv)/5

onelinSNV
morethanonelinSNV

#### Shared SNVs across Mouse lineages at higher cut-off ####
Mouse_SNVs_2$p1_L1<-as.numeric(Mouse_SNVs_2$p1_L1)
Mouse_SNVs_2$p3_L1<-as.numeric(Mouse_SNVs_2$p3_L1)
Mouse_SNVs_2$p4_L1<-as.numeric(Mouse_SNVs_2$p4_L1)
Mouse_SNVs_2$p10_L1<-as.numeric(Mouse_SNVs_2$p10_L1)
Mouse_SNVs_2$p1_L2<-as.numeric(Mouse_SNVs_2$p1_L2)
Mouse_SNVs_2$p3_L2<-as.numeric(Mouse_SNVs_2$p3_L2)
Mouse_SNVs_2$p4_L2<-as.numeric(Mouse_SNVs_2$p4_L2)
Mouse_SNVs_2$p10_L2<-as.numeric(Mouse_SNVs_2$p10_L2)
Mouse_SNVs_2$p1_L3<-as.numeric(Mouse_SNVs_2$p1_L3)
Mouse_SNVs_2$p3_L3<-as.numeric(Mouse_SNVs_2$p3_L3)
Mouse_SNVs_2$p4_L3<-as.numeric(Mouse_SNVs_2$p4_L3)
Mouse_SNVs_2$p10_L3<-as.numeric(Mouse_SNVs_2$p10_L3)
Mouse_SNVs_2$p1_L4<-as.numeric(Mouse_SNVs_2$p1_L4)
Mouse_SNVs_2$p3_L4<-as.numeric(Mouse_SNVs_2$p3_L4)
Mouse_SNVs_2$p4_L4<-as.numeric(Mouse_SNVs_2$p4_L4)
Mouse_SNVs_2$p10_L4<-as.numeric(Mouse_SNVs_2$p10_L4)
Mouse_SNVs_2$p1_L5<-as.numeric(Mouse_SNVs_2$p1_L5)
Mouse_SNVs_2$p3_L5<-as.numeric(Mouse_SNVs_2$p3_L5)
Mouse_SNVs_2$p4_L5<-as.numeric(Mouse_SNVs_2$p4_L5)
Mouse_SNVs_2$p10_L5<-as.numeric(Mouse_SNVs_2$p10_L5)

## Mouse Lineage 1_SNVs shared between Mouse lineages ##
Mouse_Lineage1_only<-merge(Mouse_Lineage1, Mouse_SNVs_2, by = "SNVID")
Mouse_Lineage1_only$lin1freqsum<-Mouse_Lineage1_only$p1_L1 + Mouse_Lineage1_only$p3_L1 + Mouse_Lineage1_only$p4_L1 + Mouse_Lineage1_only$p10_L1
Mouse_Lineage1_only$alllinfreqsum<-Mouse_Lineage1_only$p1_L1 + Mouse_Lineage1_only$p3_L1 + Mouse_Lineage1_only$p4_L1 + Mouse_Lineage1_only$p10_L1 + 
  Mouse_Lineage1_only$p1_L2 + Mouse_Lineage1_only$p3_L2 + Mouse_Lineage1_only$p4_L2 + Mouse_Lineage1_only$p10_L2 + 
  Mouse_Lineage1_only$p1_L3 + Mouse_Lineage1_only$p3_L3 + Mouse_Lineage1_only$p4_L3 + Mouse_Lineage1_only$p10_L3 + 
  Mouse_Lineage1_only$p1_L4 + Mouse_Lineage1_only$p3_L4 + Mouse_Lineage1_only$p4_L4 + Mouse_Lineage1_only$p10_L4 + 
  Mouse_Lineage1_only$p1_L5 + Mouse_Lineage1_only$p3_L5 + Mouse_Lineage1_only$p4_L5 + Mouse_Lineage1_only$p10_L5
Mouse_Lineage1_only$lin1only<-ifelse(Mouse_Lineage1_only$lin1freqsum==Mouse_Lineage1_only$alllinfreqsum,1,0)

ggplot(data=Mouse_Lineage1_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin1only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/SCP_linA_SNVs_overtime_line_MouseSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mouse Lineage 2_SNVs shared between Mouse lineages ##
Mouse_Lineage2_only<-merge(Mouse_Lineage2, Mouse_SNVs_2, by = "SNVID")
Mouse_Lineage2_only$lin2freqsum<-Mouse_Lineage2_only$p1_L2 + Mouse_Lineage2_only$p3_L2 + Mouse_Lineage2_only$p4_L2 + Mouse_Lineage2_only$p10_L2
Mouse_Lineage2_only$alllinfreqsum<-Mouse_Lineage2_only$p1_L1 + Mouse_Lineage2_only$p3_L1 + Mouse_Lineage2_only$p4_L1 + Mouse_Lineage2_only$p10_L1 + 
  Mouse_Lineage2_only$p1_L2 + Mouse_Lineage2_only$p3_L2 + Mouse_Lineage2_only$p4_L2 + Mouse_Lineage2_only$p10_L2 + 
  Mouse_Lineage2_only$p1_L3 + Mouse_Lineage2_only$p3_L3 + Mouse_Lineage2_only$p4_L3 + Mouse_Lineage2_only$p10_L3 + 
  Mouse_Lineage2_only$p1_L4 + Mouse_Lineage2_only$p3_L4 + Mouse_Lineage2_only$p4_L4 + Mouse_Lineage2_only$p10_L4 + 
  Mouse_Lineage2_only$p1_L5 + Mouse_Lineage2_only$p3_L5 + Mouse_Lineage2_only$p4_L5 + Mouse_Lineage2_only$p10_L5
Mouse_Lineage2_only$lin2only<-ifelse(Mouse_Lineage2_only$lin2freqsum==Mouse_Lineage2_only$alllinfreqsum,1,0)

ggplot(data=Mouse_Lineage2_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin2only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/SCP_linB_SNVs_overtime_line_MouseSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mouse Lineage 3_SNVs shared between Mouse lineages ##
Mouse_Lineage3_only<-merge(Mouse_Lineage3, Mouse_SNVs_2, by = "SNVID")
Mouse_Lineage3_only$lin3freqsum<-Mouse_Lineage3_only$p1_L3 + Mouse_Lineage3_only$p3_L3 + Mouse_Lineage3_only$p4_L3 + Mouse_Lineage3_only$p10_L3
Mouse_Lineage3_only$alllinfreqsum<-Mouse_Lineage3_only$p1_L1 + Mouse_Lineage3_only$p3_L1 + Mouse_Lineage3_only$p4_L1 + Mouse_Lineage3_only$p10_L1 + 
  Mouse_Lineage3_only$p1_L2 + Mouse_Lineage3_only$p3_L2 + Mouse_Lineage3_only$p4_L2 + Mouse_Lineage3_only$p10_L2 + 
  Mouse_Lineage3_only$p1_L3 + Mouse_Lineage3_only$p3_L3 + Mouse_Lineage3_only$p4_L3 + Mouse_Lineage3_only$p10_L3 + 
  Mouse_Lineage3_only$p1_L4 + Mouse_Lineage3_only$p3_L4 + Mouse_Lineage3_only$p4_L4 + Mouse_Lineage3_only$p10_L4 + 
  Mouse_Lineage3_only$p1_L5 + Mouse_Lineage3_only$p3_L5 + Mouse_Lineage3_only$p4_L5 + Mouse_Lineage3_only$p10_L5
Mouse_Lineage3_only$lin3only<-ifelse(Mouse_Lineage3_only$lin3freqsum==Mouse_Lineage3_only$alllinfreqsum,1,0)

ggplot(data=Mouse_Lineage3_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin3only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/SCP_linC_SNVs_overtime_line_MouseSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mouse Lineage 4_SNVs shared between Mouse lineages ##
Mouse_Lineage4_only<-merge(Mouse_Lineage4, Mouse_SNVs_2, by = "SNVID")
Mouse_Lineage4_only$lin4freqsum<-Mouse_Lineage4_only$p1_L4 + Mouse_Lineage4_only$p3_L4 + Mouse_Lineage4_only$p4_L4 + Mouse_Lineage4_only$p10_L4
Mouse_Lineage4_only$alllinfreqsum<-Mouse_Lineage4_only$p1_L1 + Mouse_Lineage4_only$p3_L1 + Mouse_Lineage4_only$p4_L1 + Mouse_Lineage4_only$p10_L1 + 
  Mouse_Lineage4_only$p1_L2 + Mouse_Lineage4_only$p3_L2 + Mouse_Lineage4_only$p4_L2 + Mouse_Lineage4_only$p10_L2 + 
  Mouse_Lineage4_only$p1_L3 + Mouse_Lineage4_only$p3_L3 + Mouse_Lineage4_only$p4_L3 + Mouse_Lineage4_only$p10_L3 + 
  Mouse_Lineage4_only$p1_L4 + Mouse_Lineage4_only$p3_L4 + Mouse_Lineage4_only$p4_L4 + Mouse_Lineage4_only$p10_L4 + 
  Mouse_Lineage4_only$p1_L5 + Mouse_Lineage4_only$p3_L5 + Mouse_Lineage4_only$p4_L5 + Mouse_Lineage4_only$p10_L5
Mouse_Lineage4_only$lin4only<-ifelse(Mouse_Lineage4_only$lin4freqsum==Mouse_Lineage4_only$alllinfreqsum,1,0)

ggplot(data=Mouse_Lineage4_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin4only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/SCP_linD_SNVs_overtime_line_MouseSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mouse Lineage 5_SNVs shared between Mouse lineages ##
Mouse_Lineage5_only<-merge(Mouse_Lineage5, Mouse_SNVs_2, by = "SNVID")
Mouse_Lineage5_only$lin5freqsum<-Mouse_Lineage5_only$p1_L5 + Mouse_Lineage5_only$p3_L5 + Mouse_Lineage5_only$p4_L5 + Mouse_Lineage5_only$p10_L5
Mouse_Lineage5_only$alllinfreqsum<-Mouse_Lineage5_only$p1_L1 + Mouse_Lineage5_only$p3_L1 + Mouse_Lineage5_only$p4_L1 + Mouse_Lineage5_only$p10_L1 + 
  Mouse_Lineage5_only$p1_L2 + Mouse_Lineage5_only$p3_L2 + Mouse_Lineage5_only$p4_L2 + Mouse_Lineage5_only$p10_L2 + 
  Mouse_Lineage5_only$p1_L3 + Mouse_Lineage5_only$p3_L3 + Mouse_Lineage5_only$p4_L3 + Mouse_Lineage5_only$p10_L3 + 
  Mouse_Lineage5_only$p1_L4 + Mouse_Lineage5_only$p3_L4 + Mouse_Lineage5_only$p4_L4 + Mouse_Lineage5_only$p10_L4 + 
  Mouse_Lineage5_only$p1_L5 + Mouse_Lineage5_only$p3_L5 + Mouse_Lineage5_only$p4_L5 + Mouse_Lineage5_only$p10_L5
Mouse_Lineage5_only$lin5only<-ifelse(Mouse_Lineage5_only$lin5freqsum==Mouse_Lineage5_only$alllinfreqsum,1,0)

ggplot(data=Mouse_Lineage5_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin5only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/SCP_linE_SNVs_overtime_line_MouseSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Concatenate all lineage only columns
linonly<-as.data.frame(cbind(Mouse_Lineage1_only$lin1only,Mouse_Lineage2_only$lin2only,Mouse_Lineage3_only$lin3only,Mouse_Lineage4_only$lin4only,Mouse_Lineage5_only$lin5only))
linonly$combined<-(linonly$V1+linonly$V2+linonly$V3+linonly$V4+linonly$V5)
linonly$onelinsnv<-ifelse(linonly$combined>0,0,1)
linonly$morethanonelinsnv<-ifelse(linonly$combined>0,1,0)

onelinSNV<-sum(linonly$onelinsnv)/5
morethanonelinSNV<-sum(linonly$morethanonelinsnv)/5

onelinSNV
morethanonelinSNV

#### Shared SNVs across Alternate lineages with higher cut-off ####

Alternate_SNVs_2[,2:37] %<>% mutate_if(is.character,as.numeric)

## Alternate Lineage 1_SNVs shared between alternate lineages ##

Alternate_Lineage1_only<-merge(Alternate_Lineage1, Alternate_SNVs_2, by = "SNVID")
Alternate_Lineage1_only<-subset(Alternate_Lineage1_only,Alternate_Lineage1_only$passage<11)
Alternate_Lineage1_only<-Alternate_Lineage1_only %>% mutate(lin1freqsum=rowSums(na.rm=TRUE,select(.,p1_L1,p2_L1,p3_L1,p4_L1,p5_L1,p6_L1,p7_L1,p8_L1,p9_L1,p10_L1)))
Alternate_Lineage1_only<-Alternate_Lineage1_only %>% mutate(alllinfreqsum=rowSums(na.rm=TRUE,select(.,p1_L1,p2_L1,p3_L1,p4_L1,p5_L1,p6_L1,p7_L1,p8_L1,p9_L1,p10_L1,
                                                                                                    p1_L2,p2_L2,p3_L2,
                                                                                                    p1_L3,p2_L3,p3_L3,
                                                                                                    p1_L4,p2_L4,p3_L4,p4_L4,p5_L4,p6_L4,p7_L4,p8_L4,p9_L4,p10_L4,
                                                                                                    p1_L5,p2_L5,p3_L5,p4_L5,p5_L5,p6_L5,p7_L5,p8_L5,p9_L5,p10_L5)))
Alternate_Lineage1_only$lin1only<-ifelse(Alternate_Lineage1_only$lin1freqsum==Alternate_Lineage1_only$alllinfreqsum,1,0)

ggplot(data=Alternate_Lineage1_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin1only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/AP_linA_SNVs_overtime_line_AlternateSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Alternate Lineage 2_SNVs shared between alternate lineages  ##
Alternate_Lineage2_only<-merge(Alternate_Lineage2, Alternate_SNVs_2, by = "SNVID")
Alternate_Lineage2_only<-Alternate_Lineage2_only %>% mutate(lin2freqsum=rowSums(na.rm=TRUE,select(.,p1_L2,p2_L2,p3_L2)))
Alternate_Lineage2_only<-Alternate_Lineage2_only %>% mutate(alllinfreqsum=rowSums(na.rm=TRUE,select(.,p1_L1,p2_L1,p3_L1,p4_L1,p5_L1,p6_L1,p7_L1,p8_L1,p9_L1,p10_L1,
                                                                                                    p1_L2,p2_L2,p3_L2,
                                                                                                    p1_L3,p2_L3,p3_L3,
                                                                                                    p1_L4,p2_L4,p3_L4,p4_L4,p5_L4,p6_L4,p7_L4,p8_L4,p9_L4,p10_L4,
                                                                                                    p1_L5,p2_L5,p3_L5,p4_L5,p5_L5,p6_L5,p7_L5,p8_L5,p9_L5,p10_L5)))
Alternate_Lineage2_only$lin2only<-ifelse(Alternate_Lineage2_only$lin2freqsum==Alternate_Lineage2_only$alllinfreqsum,1,0)

ggplot(data=Alternate_Lineage2_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin2only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/AP_linB_SNVs_overtime_line_AlternateSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Alternate Lineage 3_SNVs shared between alternate lineages  ##
Alternate_Lineage3_only<-merge(Alternate_Lineage3, Alternate_SNVs_2, by = "SNVID")
Alternate_Lineage3_only<-Alternate_Lineage3_only %>% mutate(lin3freqsum=rowSums(na.rm=TRUE,select(.,p1_L3,p2_L3,p3_L3)))
Alternate_Lineage3_only<-Alternate_Lineage3_only %>% mutate(alllinfreqsum=rowSums(na.rm=TRUE,select(.,p1_L1,p2_L1,p3_L1,p4_L1,p5_L1,p6_L1,p7_L1,p8_L1,p9_L1,p10_L1,
                                                                                                    p1_L2,p2_L2,p3_L2,
                                                                                                    p1_L3,p2_L3,p3_L3,
                                                                                                    p1_L4,p2_L4,p3_L4,p4_L4,p5_L4,p6_L4,p7_L4,p8_L4,p9_L4,p10_L4,
                                                                                                    p1_L5,p2_L5,p3_L5,p4_L5,p5_L5,p6_L5,p7_L5,p8_L5,p9_L5,p10_L5)))
Alternate_Lineage3_only$lin3only<-ifelse(Alternate_Lineage3_only$lin3freqsum==Alternate_Lineage3_only$alllinfreqsum,1,0)

Alternate_Lineage3_only<-subset(Alternate_Lineage3_only,Alternate_Lineage3_only$passage!=2)


ggplot(data=Alternate_Lineage3_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin3only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/AP_linC_SNVs_overtime_line_AlternateSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


## Alternate Lineage 4_SNVs shared between alternate lineages  ##
Alternate_Lineage4_only<-merge(Alternate_Lineage4, Alternate_SNVs_2, by = "SNVID")
Alternate_Lineage4_only<-subset(Alternate_Lineage4_only,Alternate_Lineage4_only$passage<11)
Alternate_Lineage4_only<-Alternate_Lineage4_only %>% mutate(lin4freqsum=rowSums(na.rm=TRUE,select(.,p1_L4,p2_L4,p3_L4,p4_L4,p5_L4,p6_L4,p7_L4,p8_L4,p9_L4,p10_L4)))
Alternate_Lineage4_only<-Alternate_Lineage4_only %>% mutate(alllinfreqsum=rowSums(na.rm=TRUE,select(.,p1_L1,p2_L1,p3_L1,p4_L1,p5_L1,p6_L1,p7_L1,p8_L1,p9_L1,p10_L1,
                                                                                                    p1_L2,p2_L2,p3_L2,
                                                                                                    p1_L3,p2_L3,p3_L3,
                                                                                                    p1_L4,p2_L4,p3_L4,p4_L4,p5_L4,p6_L4,p7_L4,p8_L4,p9_L4,p10_L4,
                                                                                                    p1_L5,p2_L5,p3_L5,p4_L5,p5_L5,p6_L5,p7_L5,p8_L5,p9_L5,p10_L5)))
Alternate_Lineage4_only$lin4only<-ifelse(Alternate_Lineage4_only$lin4freqsum==Alternate_Lineage4_only$alllinfreqsum,1,0)

ggplot(data=Alternate_Lineage4_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin4only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/AP_linD_SNVs_overtime_line_AlternateSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Alternate Lineage 5_SNVs shared between alternate lineages  ##
Alternate_Lineage5_only<-merge(Alternate_Lineage5, Alternate_SNVs_2, by = "SNVID")
Alternate_Lineage5_only<-Alternate_Lineage5_only %>% mutate(lin5freqsum=rowSums(na.rm=TRUE,select(.,p1_L5,p2_L5,p3_L5,p4_L5,p5_L5,p6_L5,p7_L5,p8_L5,p9_L5,p10_L5)))
Alternate_Lineage5_only<-Alternate_Lineage5_only %>% mutate(alllinfreqsum=rowSums(na.rm=TRUE,select(.,p1_L1,p2_L1,p3_L1,p4_L1,p5_L1,p6_L1,p7_L1,p8_L1,p9_L1,p10_L1,
                                                                                                    p1_L2,p2_L2,p3_L2,
                                                                                                    p1_L3,p2_L3,p3_L3,
                                                                                                    p1_L4,p2_L4,p3_L4,p4_L4,p5_L4,p6_L4,p7_L4,p8_L4,p9_L4,p10_L4,
                                                                                                    p1_L5,p2_L5,p3_L5,p4_L5,p5_L5,p6_L5,p7_L5,p8_L5,p9_L5,p10_L5)))
Alternate_Lineage5_only$lin5only<-ifelse(Alternate_Lineage5_only$lin5freqsum==Alternate_Lineage5_only$alllinfreqsum,1,0)

Alternate_Lineage5_only<-subset(Alternate_Lineage5_only,Alternate_Lineage5_only$passage!=2)


ggplot(data=Alternate_Lineage5_only, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(lin5only)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits=c(1,10),breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"))

ggsave("PATH/TO/AP_linE_SNVs_overtime_line_AlternateSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Concatenate all lineage only columns
lin1only<-subset(Alternate_Lineage1_only,Alternate_Lineage1_only$passage==10)
lin2only<-subset(Alternate_Lineage2_only,Alternate_Lineage2_only$passage==3)
lin3only<-subset(Alternate_Lineage3_only,Alternate_Lineage3_only$passage==3)
lin4only<-subset(Alternate_Lineage4_only,Alternate_Lineage4_only$passage==10)
lin5only<-subset(Alternate_Lineage5_only,Alternate_Lineage5_only$passage==10)
linonly<-as.data.frame(cbind(lin1only$lin1only,lin2only$lin2only,lin3only$lin3only,lin4only$lin4only,lin5only$lin5only))
linonly$combined<-(linonly$V1+linonly$V2+linonly$V3+linonly$V4+linonly$V5)
linonly$onelinsnv<-ifelse(linonly$combined>0,0,1)
linonly$morethanonelinsnv<-ifelse(linonly$combined>0,1,0)

onelinSNV<-sum(linonly$onelinsnv)
morethanonelinSNV<-sum(linonly$morethanonelinsnv)

onelinSNV
morethanonelinSNV

#### COMPOSITE Shared SNVs across lineages ####
Mosquito_LineageOnly<-Mosquito_Lineage1_only
Mosquito_LineageOnly$lin2only<-Mosquito_Lineage2_only$lin2only
Mosquito_LineageOnly$lin3only<-Mosquito_Lineage3_only$lin3only
Mosquito_LineageOnly$lin4only<-Mosquito_Lineage4_only$lin4only
Mosquito_LineageOnly$lin5only<-Mosquito_Lineage5_only$lin5only
Mosquito_LineageOnly$Sumlinonly<-Mosquito_LineageOnly$lin1only+Mosquito_LineageOnly$lin2only+Mosquito_LineageOnly$lin3only+Mosquito_LineageOnly$lin4only+Mosquito_LineageOnly$lin5only
Mosquito_LineageOnly$Onelinonly<-ifelse(Mosquito_LineageOnly$Sumlinonly>0,1,0)

MosqL1<-subset(Mosquito_Lineage1_only,Mosquito_Lineage1_only$lin1only>0)
MosqL1$onelin<-MosqL1$lin1only
MosqL1<-MosqL1[,c(1,2,3,27)]
MosqL2<-subset(Mosquito_Lineage2_only,Mosquito_Lineage2_only$lin2only>0)
MosqL2$onelin<-MosqL2$lin2only
MosqL2<-MosqL2[,c(1,2,3,27)]
MosqL3<-subset(Mosquito_Lineage3_only,Mosquito_Lineage3_only$lin3only>0)
MosqL3$onelin<-MosqL3$lin3only
MosqL3<-MosqL3[,c(1,2,3,27)]
MosqL4<-subset(Mosquito_Lineage4_only,Mosquito_Lineage4_only$lin4only>0)
MosqL4$onelin<-MosqL4$lin4only
MosqL4<-MosqL4[,c(1,2,3,27)]
MosqL5<-subset(Mosquito_Lineage5_only,Mosquito_Lineage5_only$lin5only>0)
MosqL5$onelin<-MosqL5$lin5only
MosqL5<-MosqL5[,c(1,2,3,27)]

Mosq_1lin<-rbind(MosqL1,MosqL2)
Mosq_1lin<-rbind(Mosq_1lin,MosqL3)
Mosq_1lin<-rbind(Mosq_1lin,MosqL4)
Mosq_1lin<-rbind(Mosq_1lin,MosqL5)

MosqL1_shared<-subset(Mosquito_Lineage1_only,Mosquito_Lineage1_only$lin1only<1)
MosqL1_shared$onelin<-MosqL1_shared$lin1only
MosqL1_shared<-MosqL1_shared[,c(1,2,3,27)]
MosqL1_shared$SNVIDPASS<-paste(MosqL1_shared$SNVID,MosqL1_shared$passage,sep=",")
MosqL2_shared<-subset(Mosquito_Lineage2_only,Mosquito_Lineage2_only$lin2only<1)
MosqL2_shared$onelin<-MosqL2_shared$lin2only
MosqL2_shared<-MosqL2_shared[,c(1,2,3,27)]
MosqL2_shared$SNVIDPASS<-paste(MosqL2_shared$SNVID,MosqL2_shared$passage,sep=",")
MosqL3_shared<-subset(Mosquito_Lineage3_only,Mosquito_Lineage3_only$lin3only<1)
MosqL3_shared$onelin<-MosqL3_shared$lin3only
MosqL3_shared<-MosqL3_shared[,c(1,2,3,27)]
MosqL3_shared$SNVIDPASS<-paste(MosqL3_shared$SNVID,MosqL3_shared$passage,sep=",")
MosqL4_shared<-subset(Mosquito_Lineage4_only,Mosquito_Lineage4_only$lin4only<1)
MosqL4_shared$onelin<-MosqL4_shared$lin4only
MosqL4_shared<-MosqL4_shared[,c(1,2,3,27)]
MosqL4_shared$SNVIDPASS<-paste(MosqL4_shared$SNVID,MosqL4_shared$passage,sep=",")
MosqL5_shared<-subset(Mosquito_Lineage5_only,Mosquito_Lineage5_only$lin5only<1)
MosqL5_shared$onelin<-MosqL5_shared$lin5only
MosqL5_shared<-MosqL5_shared[,c(1,2,3,27)]
MosqL5_shared$SNVIDPASS<-paste(MosqL5_shared$SNVID,MosqL5_shared$passage,sep=",")

Mosq_shared<-merge(MosqL1_shared, MosqL2_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mosq_shared<-Mosq_shared[,c(1,3,7)]
colnames(Mosq_shared)<-c("SNVIDPASS","L1Freq","L2Freq")
Mosq_shared<-merge(Mosq_shared, MosqL3_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mosq_shared<-Mosq_shared[,c(1:3,5)]
colnames(Mosq_shared)<-c("SNVIDPASS","L1Freq","L2Freq","L3Freq")
Mosq_shared<-merge(Mosq_shared, MosqL4_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mosq_shared<-Mosq_shared[,c(1:4,6)]
colnames(Mosq_shared)<-c("SNVIDPASS","L1Freq","L2Freq","L3Freq","L4Freq")
Mosq_shared<-merge(Mosq_shared, MosqL5_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mosq_shared<-Mosq_shared[,c(1:5,7)]
colnames(Mosq_shared)<-c("SNVIDPASS","L1Freq","L2Freq","L3Freq","L4Freq","L5Freq")

Mosq_shared[Mosq_shared < 0.001]<-as.numeric("NA")

Mosq_shared<-Mosq_shared %>% mutate(freq=rowMeans(na.rm=TRUE,select(., L1Freq,L2Freq,L3Freq,L4Freq,L5Freq)))
Mosq_shared[is.na(Mosq_shared)]<-0
Mosq_shared<-cSplit(Mosq_shared, 'SNVIDPASS', sep=",", type.convert=FALSE)
Mosq_shared$SNVIDPASS_2<-as.numeric(Mosq_shared$SNVIDPASS_2)
Mosq_shared<-Mosq_shared[order(Mosq_shared$SNVIDPASS_1, Mosq_shared$SNVIDPASS_2),]

Mosq_shared$prezero<-ifelse(((dplyr::lead(Mosq_shared$freq)>0) | (Mosq_shared$freq==0 & dplyr::lag(Mosq_shared$SNVIDPASS_2)<10 & dplyr::lag(Mosq_shared$freq)>0)),1,0)

Mosq_shared<-subset(Mosq_shared,Mosq_shared$freq>0 | Mosq_shared$prezero>0)
Mosq_shared<-Mosq_shared[,c(7,6,8)]
colnames(Mosq_shared)<-c("SNVID","freq","passage")
Mosq_shared$onelin<-0

Mosq_ONELIN<-rbind(Mosq_1lin,Mosq_shared)
Mosq_ONELIN$passage<-as.numeric(Mosq_ONELIN$passage)
Mosq_ONELIN$SNVID<-as.factor(Mosq_ONELIN$SNVID)
Mosq_ONELIN$onelin<-as.factor(Mosq_ONELIN$onelin)
Mosq_ONELIN<-subset(Mosq_ONELIN,Mosq_ONELIN$passage!=11)

ggplot(data=Mosq_ONELIN, aes(x=passage, y=freq)) +
  geom_path(aes(group = SNVID, color = onelin)) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), legend.title = element_blank(), legend.key = element_blank()) +
  ggtitle("ZIKV Mosquito Composite") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"), labels = c("Convergent","1 Lineage"))

ggsave("PATH/TO/MP_COMPOSITE_SNVs_overtime_line_MosquitoSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


### Mice
Mouse_LineageOnly<-Mouse_Lineage1_only
Mouse_LineageOnly$lin2only<-Mouse_Lineage2_only$lin2only
Mouse_LineageOnly$lin3only<-Mouse_Lineage3_only$lin3only
Mouse_LineageOnly$lin4only<-Mouse_Lineage4_only$lin4only
Mouse_LineageOnly$lin5only<-Mouse_Lineage5_only$lin5only
Mouse_LineageOnly$Sumlinonly<-Mouse_LineageOnly$lin1only+Mouse_LineageOnly$lin2only+Mouse_LineageOnly$lin3only+Mouse_LineageOnly$lin4only+Mouse_LineageOnly$lin5only
Mouse_LineageOnly$Onelinonly<-ifelse(Mouse_LineageOnly$Sumlinonly>0,1,0)

MousL1<-subset(Mouse_Lineage1_only,Mouse_Lineage1_only$lin1only>0)
MousL1$onelin<-MousL1$lin1only
MousL1<-MousL1[,c(1,2,3,27)]
MousL2<-subset(Mouse_Lineage2_only,Mouse_Lineage2_only$lin2only>0)
MousL2$onelin<-MousL2$lin2only
MousL2<-MousL2[,c(1,2,3,27)]
MousL3<-subset(Mouse_Lineage3_only,Mouse_Lineage3_only$lin3only>0)
MousL3$onelin<-MousL3$lin3only
MousL3<-MousL3[,c(1,2,3,27)]
MousL4<-subset(Mouse_Lineage4_only,Mouse_Lineage4_only$lin4only>0)
MousL4$onelin<-MousL4$lin4only
MousL4<-MousL4[,c(1,2,3,27)]
MousL5<-subset(Mouse_Lineage5_only,Mouse_Lineage5_only$lin5only>0)
MousL5$onelin<-MousL5$lin5only
MousL5<-MousL5[,c(1,2,3,27)]

Mous_1lin<-rbind(MousL1,MousL2)
Mous_1lin<-rbind(Mous_1lin,MousL3)
Mous_1lin<-rbind(Mous_1lin,MousL4)
Mous_1lin<-rbind(Mous_1lin,MousL5)

MousL1_shared<-subset(Mouse_Lineage1_only,Mouse_Lineage1_only$lin1only<1)
MousL1_shared$onelin<-MousL1_shared$lin1only
MousL1_shared<-MousL1_shared[,c(1,2,3,27)]
MousL1_shared$SNVIDPASS<-paste(MousL1_shared$SNVID,MousL1_shared$passage,sep=",")
MousL2_shared<-subset(Mouse_Lineage2_only,Mouse_Lineage2_only$lin2only<1)
MousL2_shared$onelin<-MousL2_shared$lin2only
MousL2_shared<-MousL2_shared[,c(1,2,3,27)]
MousL2_shared$SNVIDPASS<-paste(MousL2_shared$SNVID,MousL2_shared$passage,sep=",")
MousL3_shared<-subset(Mouse_Lineage3_only,Mouse_Lineage3_only$lin3only<1)
MousL3_shared$onelin<-MousL3_shared$lin3only
MousL3_shared<-MousL3_shared[,c(1,2,3,27)]
MousL3_shared$SNVIDPASS<-paste(MousL3_shared$SNVID,MousL3_shared$passage,sep=",")
MousL4_shared<-subset(Mouse_Lineage4_only,Mouse_Lineage4_only$lin4only<1)
MousL4_shared$onelin<-MousL4_shared$lin4only
MousL4_shared<-MousL4_shared[,c(1,2,3,27)]
MousL4_shared$SNVIDPASS<-paste(MousL4_shared$SNVID,MousL4_shared$passage,sep=",")
MousL5_shared<-subset(Mouse_Lineage5_only,Mouse_Lineage5_only$lin5only<1)
MousL5_shared$onelin<-MousL5_shared$lin5only
MousL5_shared<-MousL5_shared[,c(1,2,3,27)]
MousL5_shared$SNVIDPASS<-paste(MousL5_shared$SNVID,MousL5_shared$passage,sep=",")

Mous_shared<-merge(MousL1_shared, MousL2_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mous_shared<-Mous_shared[,c(1,3,7)]
colnames(Mous_shared)<-c("SNVIDPASS","L1Freq","L2Freq")
Mous_shared<-merge(Mous_shared, MousL3_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mous_shared<-Mous_shared[,c(1:3,5)]
colnames(Mous_shared)<-c("SNVIDPASS","L1Freq","L2Freq","L3Freq")
Mous_shared<-merge(Mous_shared, MousL4_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mous_shared<-Mous_shared[,c(1:4,6)]
colnames(Mous_shared)<-c("SNVIDPASS","L1Freq","L2Freq","L3Freq","L4Freq")
Mous_shared<-merge(Mous_shared, MousL5_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mous_shared<-Mous_shared[,c(1:5,7)]
colnames(Mous_shared)<-c("SNVIDPASS","L1Freq","L2Freq","L3Freq","L4Freq","L5Freq")

Mous_shared[Mous_shared < 0.001]<-as.numeric("NA")

Mous_shared<-Mous_shared %>% mutate(freq=rowMeans(na.rm=TRUE,select(., L1Freq,L2Freq,L3Freq,L4Freq,L5Freq)))
Mous_shared[is.na(Mous_shared)]<-0
Mous_shared<-cSplit(Mous_shared, 'SNVIDPASS', sep=",", type.convert=FALSE)
Mous_shared$SNVIDPASS_2<-as.numeric(Mous_shared$SNVIDPASS_2)
Mous_shared<-Mous_shared[order(Mous_shared$SNVIDPASS_1, Mous_shared$SNVIDPASS_2),]

Mous_shared$prezero<-ifelse(((dplyr::lead(Mous_shared$freq)>0) | (Mous_shared$freq==0 & dplyr::lag(Mous_shared$SNVIDPASS_2)<10 & dplyr::lag(Mous_shared$freq)>0)),1,0)

Mous_shared<-subset(Mous_shared,Mous_shared$freq>0 | Mous_shared$prezero>0)
Mous_shared<-Mous_shared[,c(7,6,8)]
colnames(Mous_shared)<-c("SNVID","freq","passage")
Mous_shared$onelin<-0

Mous_ONELIN<-rbind(Mous_1lin,Mous_shared)
Mous_ONELIN$passage<-as.numeric(Mous_ONELIN$passage)
Mous_ONELIN$SNVID<-as.factor(Mous_ONELIN$SNVID)
Mous_ONELIN$onelin<-as.factor(Mous_ONELIN$onelin)
Mous_ONELIN<-subset(Mous_ONELIN,Mous_ONELIN$passage!=11)

ggplot(data=Mous_ONELIN, aes(x=passage, y=freq)) +
  geom_path(aes(group = SNVID, color = onelin)) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), legend.title = element_blank(), legend.key = element_blank()) +
  ggtitle("ZIKV Mouse Composite") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"), labels = c("Convergent","1 Lineage"))

ggsave("PATH/TO/SCP_COMPOSITE_SNVs_overtime_line_MouseSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

### Alternate

Alternate_LineageOnly<-merge(Alternate_Lineage1, Alternate_SNVs_2, by = "SNVID")
Alternate_LineageOnly<-subset(Alternate_LineageOnly,Alternate_LineageOnly$passage<11)
Alternate_LineageOnly<-Alternate_LineageOnly %>% mutate(lin1freqsum=rowSums(na.rm=TRUE,select(.,p1_L1,p2_L1,p3_L1,p4_L1,p5_L1,p6_L1,p7_L1,p8_L1,p9_L1,p10_L1)))
Alternate_LineageOnly<-Alternate_LineageOnly %>% mutate(alllinfreqsum=rowSums(na.rm=TRUE,select(.,p1_L1,p2_L1,p3_L1,p4_L1,p5_L1,p6_L1,p7_L1,p8_L1,p9_L1,p10_L1,
                                                                                                    p1_L2,p2_L2,p3_L2,
                                                                                                    p1_L3,p2_L3,p3_L3,
                                                                                                    p1_L4,p2_L4,p3_L4,p4_L4,p5_L4,p6_L4,p7_L4,p8_L4,p9_L4,p10_L4,
                                                                                                    p1_L5,p2_L5,p3_L5,p4_L5,p5_L5,p6_L5,p7_L5,p8_L5,p9_L5,p10_L5)))
Alternate_LineageOnly$lin1only<-ifelse(Alternate_LineageOnly$lin1freqsum==Alternate_LineageOnly$alllinfreqsum,1,0)

Alternate_LineageOnly<-Alternate_LineageOnly %>% mutate(lin2freqsum=rowSums(na.rm=TRUE,select(.,p1_L2,p2_L2,p3_L2)))
Alternate_LineageOnly$lin2only<-ifelse(Alternate_LineageOnly$lin2freqsum==Alternate_LineageOnly$alllinfreqsum,1,0)

Alternate_LineageOnly<-Alternate_LineageOnly %>% mutate(lin3freqsum=rowSums(na.rm=TRUE,select(.,p1_L3,p2_L3,p3_L3)))
Alternate_LineageOnly$lin3only<-ifelse(Alternate_LineageOnly$lin3freqsum==Alternate_LineageOnly$alllinfreqsum,1,0)

Alternate_LineageOnly<-Alternate_LineageOnly %>% mutate(lin4freqsum=rowSums(na.rm=TRUE,select(.,p1_L4,p2_L4,p3_L4,p4_L4,p5_L4,p6_L4,p7_L4,p8_L4,p9_L4,p10_L4)))
Alternate_LineageOnly$lin4only<-ifelse(Alternate_LineageOnly$lin4freqsum==Alternate_LineageOnly$alllinfreqsum,1,0)

Alternate_LineageOnly<-Alternate_LineageOnly %>% mutate(lin5freqsum=rowSums(na.rm=TRUE,select(.,p1_L5,p2_L5,p3_L5,p4_L5,p5_L5,p6_L5,p7_L5,p8_L5,p9_L5,p10_L5)))
Alternate_LineageOnly$lin5only<-ifelse(Alternate_LineageOnly$lin5freqsum==Alternate_LineageOnly$alllinfreqsum,1,0)

Alternate_LineageOnly[,2:50] %<>% mutate_if(is.character,as.numeric)
Alternate_LineageOnly<-Alternate_LineageOnly %>% mutate(Onelinonly=rowSums(na.rm=TRUE,select(.,lin1only,lin2only,lin3only,lin4only,lin5only)))

Alternate_LineageOnly_1lin<-subset(Alternate_LineageOnly,Alternate_LineageOnly$Onelinonly==1)
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin[,c(1,3:39,51)]
Alternate_LineageOnly_1lin[Alternate_LineageOnly_1lin < 0.001]<-as.numeric("NA")
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p1freq=rowMeans(na.rm=TRUE,select(.,p1_L1,p1_L2,p1_L3,p1_L4,p1_L5)))
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p2freq=rowMeans(na.rm=TRUE,select(.,p2_L1,p2_L2,p2_L3,p2_L4,p2_L5)))
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p3freq=rowMeans(na.rm=TRUE,select(.,p3_L1,p3_L2,p3_L3,p3_L4,p3_L5)))
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p4freq=rowMeans(na.rm=TRUE,select(.,p4_L1,p4_L4,p4_L5)))
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p5freq=rowMeans(na.rm=TRUE,select(.,p5_L1,p5_L4,p5_L5)))
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p6freq=rowMeans(na.rm=TRUE,select(.,p6_L1,p6_L4,p6_L5)))
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p7freq=rowMeans(na.rm=TRUE,select(.,p7_L1,p7_L4,p7_L5)))
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p8freq=rowMeans(na.rm=TRUE,select(.,p8_L1,p8_L4,p8_L5)))
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p9freq=rowMeans(na.rm=TRUE,select(.,p9_L1,p9_L4,p9_L5)))
Alternate_LineageOnly_1lin<-Alternate_LineageOnly_1lin %>% mutate(p10freq=rowMeans(na.rm=TRUE,select(.,p10_L1,p10_L4,p10_L5)))
Alternate_LineageOnly_1lin$freq<-ifelse(Alternate_LineageOnly_1lin$passage==1,Alternate_LineageOnly_1lin$p1freq,
                                          ifelse(Alternate_LineageOnly_1lin$passage==2,Alternate_LineageOnly_1lin$p2freq,
                                                 ifelse(Alternate_LineageOnly_1lin$passage==3,Alternate_LineageOnly_1lin$p3freq,
                                                        ifelse(Alternate_LineageOnly_1lin$passage==4,Alternate_LineageOnly_1lin$p4freq,
                                                               ifelse(Alternate_LineageOnly_1lin$passage==5,Alternate_LineageOnly_1lin$p5freq,
                                                                      ifelse(Alternate_LineageOnly_1lin$passage==6,Alternate_LineageOnly_1lin$p6freq,
                                                                             ifelse(Alternate_LineageOnly_1lin$passage==7,Alternate_LineageOnly_1lin$p7freq,
                                                                                    ifelse(Alternate_LineageOnly_1lin$passage==8,Alternate_LineageOnly_1lin$p8freq,
                                                                                           ifelse(Alternate_LineageOnly_1lin$passage==9,Alternate_LineageOnly_1lin$p9freq,
                                                                                                  ifelse(Alternate_LineageOnly_1lin$passage==10,Alternate_LineageOnly_1lin$p10freq,0))))))))))
Alternate_LineageOnly_1lin[is.na(Alternate_LineageOnly_1lin)]<-0 

Alternate_LineageOnly_shared<-subset(Alternate_LineageOnly,Alternate_LineageOnly$Onelinonly<1)
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared[,c(1,3:39,51)]
Alternate_LineageOnly_shared[Alternate_LineageOnly_shared < 0.001]<-as.numeric("NA")
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p1freq=rowMeans(na.rm=TRUE,select(.,p1_L1,p1_L2,p1_L3,p1_L4,p1_L5)))
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p2freq=rowMeans(na.rm=TRUE,select(.,p2_L1,p2_L2,p2_L3,p2_L4,p2_L5)))
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p3freq=rowMeans(na.rm=TRUE,select(.,p3_L1,p3_L2,p3_L3,p3_L4,p3_L5)))
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p4freq=rowMeans(na.rm=TRUE,select(.,p4_L1,p4_L4,p4_L5)))
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p5freq=rowMeans(na.rm=TRUE,select(.,p5_L1,p5_L4,p5_L5)))
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p6freq=rowMeans(na.rm=TRUE,select(.,p6_L1,p6_L4,p6_L5)))
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p7freq=rowMeans(na.rm=TRUE,select(.,p7_L1,p7_L4,p7_L5)))
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p8freq=rowMeans(na.rm=TRUE,select(.,p8_L1,p8_L4,p8_L5)))
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p9freq=rowMeans(na.rm=TRUE,select(.,p9_L1,p9_L4,p9_L5)))
Alternate_LineageOnly_shared<-Alternate_LineageOnly_shared %>% mutate(p10freq=rowMeans(na.rm=TRUE,select(.,p10_L1,p10_L4,p10_L5)))
Alternate_LineageOnly_shared$freq<-ifelse(Alternate_LineageOnly_shared$passage==1,Alternate_LineageOnly_shared$p1freq,
                                          ifelse(Alternate_LineageOnly_shared$passage==2,Alternate_LineageOnly_shared$p2freq,
                                                 ifelse(Alternate_LineageOnly_shared$passage==3,Alternate_LineageOnly_shared$p3freq,
                                                        ifelse(Alternate_LineageOnly_shared$passage==4,Alternate_LineageOnly_shared$p4freq,
                                                               ifelse(Alternate_LineageOnly_shared$passage==5,Alternate_LineageOnly_shared$p5freq,
                                                                      ifelse(Alternate_LineageOnly_shared$passage==6,Alternate_LineageOnly_shared$p6freq,
                                                                             ifelse(Alternate_LineageOnly_shared$passage==7,Alternate_LineageOnly_shared$p7freq,
                                                                                    ifelse(Alternate_LineageOnly_shared$passage==8,Alternate_LineageOnly_shared$p8freq,
                                                                                           ifelse(Alternate_LineageOnly_shared$passage==9,Alternate_LineageOnly_shared$p9freq,
                                                                                                  ifelse(Alternate_LineageOnly_shared$passage==10,Alternate_LineageOnly_shared$p10freq,0))))))))))
Alternate_LineageOnly_shared[is.na(Alternate_LineageOnly_shared)]<-0                                                                                                         

Alternate_ONELIN<-rbind(Alternate_LineageOnly_1lin,Alternate_LineageOnly_shared)
Alternate_ONELIN<-Alternate_ONELIN[,c(1,2,39,50)]
Alternate_ONELIN$passage<-as.numeric(Alternate_ONELIN$passage)
Alternate_ONELIN$SNVID<-as.factor(Alternate_ONELIN$SNVID)
Alternate_ONELIN$Onelinonly<-as.factor(Alternate_ONELIN$Onelinonly)


ggplot(data=Alternate_ONELIN, aes(x=passage, y=freq)) +
  geom_path(aes(group = SNVID, color = Onelinonly)) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), legend.title = element_blank(), legend.key = element_blank()) +
  ggtitle("ZIKV Alternate Composite") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","2","3","4","5","6","7","8","9","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red","grey1"), labels = c("Convergent","1 Lineage"))

ggsave("PATH/TO/AP_COMPOSITE_SNVs_overtime_line_AlternateSharedSNVs.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


Mouse_OneLin_Unique<-Mous_ONELIN %>% distinct(SNVID, .keep_all = TRUE)
Mouse_OneLin_Unique %>% count(onelin)

Mosq_OneLin_Unique<-Mosq_ONELIN %>% distinct(SNVID, .keep_all = TRUE)
Mosq_OneLin_Unique %>% count(onelin)

Alternate_OneLin_Unique<-Alternate_ONELIN %>% distinct(SNVID, .keep_all = TRUE)
Alternate_OneLin_Unique %>% count(Onelinonly)