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

p11_L5<-read.table(file="/Volumes/HD1/R21_ZIKVConspecificPassage/Alternate_Passage/WGS/AP10_E_rep1_v3.1/reference_aligned/lofreq_reference_AP10_E_1.vcf", sep="\t", header=FALSE)
p11_L5<-cSplit(p11_L5,"V8",sep=";", type.convert=FALSE)
p11_L5<-cSplit(p11_L5, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
p11_L5$SNVID<-paste(p11_L5$V2,p11_L5$V4,sep=";")
p11_L5$SNVID<-paste(p11_L5$SNVID,p11_L5$V5,sep=">")
p11_L5<-p11_L5[,c(12,11)]
colnames(p11_L5)<-c("SNVID","freq")
p11_L5<-subset(p11_L5,p11_L5$freq>=0.003)

p111_L5<-merge(p110_L5, p11_L5, by = "SNVID", all.x = TRUE, all.y = TRUE)
colnames(p111_L5)<-c("SNVID","freq1","freq3","freq4","freq5","freq6","freq7","freq8","freq9","freq10","freqCC")

p111_L5[is.na(p111_L5)]<-0


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
Lineage5_p11<-Alternate_L1thruL5[,c(1,40)]
Lineage5_p11$passage<-11
colnames(Lineage5_p11)<-c("SNVID","freq","passage")

Lineage5<-rbind(Lineage5_p1,Lineage5_p2)
Lineage5<-rbind(Lineage5,Lineage5_p3)
Lineage5<-rbind(Lineage5,Lineage5_p4)
Lineage5<-rbind(Lineage5,Lineage5_p5)
Lineage5<-rbind(Lineage5,Lineage5_p6)
Lineage5<-rbind(Lineage5,Lineage5_p7)
Lineage5<-rbind(Lineage5,Lineage5_p8)
Lineage5<-rbind(Lineage5,Lineage5_p9)
Alternate_Lineage5<-rbind(Lineage5,Lineage5_p10)
Alternate_Lineage5<-rbind(Lineage5,Lineage5_p11)

Alternate_Lineage5$freq<-as.numeric(Alternate_Lineage5$freq)
Alternate_Lineage5$passage<-as.numeric(Alternate_Lineage5$passage)

######  Define SNP species specificity           ######################################################

Mosquito_L1thruL5_1<-Mosquito_L1thruL5[,c(1:5,7:10,12:15,17:20,22:25)]
Mouse_L1thruL5_1<-Mouse_L1thruL5[,c(1:5,7:10,12:15,17:20,22:25)]
Alternate_L1thruL5_1<-Alternate_L1thruL5[,c(1:11,13:28,30:39)]

Mosquito_SNVs_2<-Mosquito_L1thruL5_1 %>% filter_at(vars(starts_with("p")),any_vars(.>0.01))
Mouse_SNVs_2<-Mouse_L1thruL5_1 %>% filter_at(vars(starts_with("p")),any_vars(.>0.01))
Alternate_SNVs_2<-Alternate_L1thruL5_1 %>% filter_at(vars(starts_with("p")),any_vars(.>0.01))

Mosquito_SNVs_2[,2:21] %<>% mutate_if(is.character,as.numeric)
Mouse_SNVs_2[,2:21] %<>% mutate_if(is.character,as.numeric)
Alternate_SNVs_2[,2:37] %<>% mutate_if(is.character,as.numeric)

Mosquito_SNVs_2$species<-1
Mouse_SNVs_2$species<-2
Mosquito_SNVs<-Mosquito_SNVs_2[,c(1,22)]
Mouse_SNVs<-Mouse_SNVs_2[,c(1,22)]
Alternate_SNVs_2$Species<-4
Alternate_SNVs<-Alternate_SNVs_2[,c(1,38)]

Mosq_mouse_SNVs<-merge(Mosquito_SNVs, Mouse_SNVs, by = "SNVID", all.x = TRUE, all.y = TRUE)
Mosq_mouse_SNVs[is.na(Mosq_mouse_SNVs)]<-0
Mosq_mouse_SNVs$species<-Mosq_mouse_SNVs$species.x + Mosq_mouse_SNVs$species.y
Mosq_mouse_SNVs %>% count(species)

Alt_mosq_mouse_SNVs<-merge(Mosq_mouse_SNVs,Alternate_SNVs, by = "SNVID", all.x = TRUE, all.y = TRUE)
Alt_mosq_mouse_SNVs[is.na(Alt_mosq_mouse_SNVs)]<-0
Alt_mosq_mouse_SNVs$speciesSUM<-Alt_mosq_mouse_SNVs$species.x + Alt_mosq_mouse_SNVs$species.y + Alt_mosq_mouse_SNVs$Species
Alt_mosq_mouse_SNVs %>% count(speciesSUM)

#### Species-specific SNVs in mosquito passage ####

## Mosquito Lineage 1 ##
Lineage1_p1<-Mosquito_SNVs_2[,c(1,2)]
Lineage1_p1$passage<-1
colnames(Lineage1_p1)<-c("SNVID","freq","passage")
Lineage1_p3<-Mosquito_SNVs_2[,c(1,3)]
Lineage1_p3$passage<-3
colnames(Lineage1_p3)<-c("SNVID","freq","passage")
Lineage1_p4<-Mosquito_SNVs_2[,c(1,4)]
Lineage1_p4$passage<-4
colnames(Lineage1_p4)<-c("SNVID","freq","passage")
Lineage1_p10<-Mosquito_SNVs_2[,c(1,5)]
Lineage1_p10$passage<-10
colnames(Lineage1_p10)<-c("SNVID","freq","passage")

Lineage1<-rbind(Lineage1_p1,Lineage1_p3)
Lineage1<-rbind(Lineage1,Lineage1_p4)
Mosquito_Lineage1_nop11<-rbind(Lineage1,Lineage1_p10)

Mosquito_Lineage1_species<-merge(Mosquito_Lineage1_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mosquito_Lineage1_species$species<-as.numeric(Mosquito_Lineage1_species$species)
Mosquito_Lineage1_species$freq<-as.numeric(as.character(Mosquito_Lineage1_species$freq))

Mosquito_Lineage1_species_present<-subset(Mosquito_Lineage1_species,Mosquito_Lineage1_species$freq>0)
Mosquito_Lineage1_species_present %>% count(species)

ggplot(data=Mosquito_Lineage1_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","red2"))

ggsave("/PATH/TO/MP_linA_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mosquito Lineage 2 ##
Lineage2_p1<-Mosquito_SNVs_2[,c(1,6)]
Lineage2_p1$passage<-1
colnames(Lineage2_p1)<-c("SNVID","freq","passage")
Lineage2_p3<-Mosquito_SNVs_2[,c(1,7)]
Lineage2_p3$passage<-3
colnames(Lineage2_p3)<-c("SNVID","freq","passage")
Lineage2_p4<-Mosquito_SNVs_2[,c(1,8)]
Lineage2_p4$passage<-4
colnames(Lineage2_p4)<-c("SNVID","freq","passage")
Lineage2_p10<-Mosquito_SNVs_2[,c(1,9)]
Lineage2_p10$passage<-10
colnames(Lineage2_p10)<-c("SNVID","freq","passage")

Lineage2<-rbind(Lineage2_p1,Lineage2_p3)
Lineage2<-rbind(Lineage2,Lineage2_p4)
Mosquito_Lineage2_nop11<-rbind(Lineage2,Lineage2_p10)

Mosquito_Lineage2_species<-merge(Mosquito_Lineage2_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mosquito_Lineage2_species$species<-as.numeric(Mosquito_Lineage2_species$species)
Mosquito_Lineage2_species$freq<-as.numeric(as.character(Mosquito_Lineage2_species$freq))

Mosquito_Lineage2_species_present<-subset(Mosquito_Lineage2_species,Mosquito_Lineage2_species$freq>0)
Mosquito_Lineage2_species_present %>% count(species)

ggplot(data=Mosquito_Lineage2_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","red2"))

ggsave("/PATH/TO/MP_linB_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mosquito Lineage 3 ##
Lineage3_p1<-Mosquito_SNVs_2[,c(1,10)]
Lineage3_p1$passage<-1
colnames(Lineage3_p1)<-c("SNVID","freq","passage")
Lineage3_p3<-Mosquito_SNVs_2[,c(1,11)]
Lineage3_p3$passage<-3
colnames(Lineage3_p3)<-c("SNVID","freq","passage")
Lineage3_p4<-Mosquito_SNVs_2[,c(1,12)]
Lineage3_p4$passage<-4
colnames(Lineage3_p4)<-c("SNVID","freq","passage")
Lineage3_p10<-Mosquito_SNVs_2[,c(1,13)]
Lineage3_p10$passage<-10
colnames(Lineage3_p10)<-c("SNVID","freq","passage")

Lineage3<-rbind(Lineage3_p1,Lineage3_p3)
Lineage3<-rbind(Lineage3,Lineage3_p4)
Mosquito_Lineage3_nop11<-rbind(Lineage3,Lineage3_p10)

Mosquito_Lineage3_species<-merge(Mosquito_Lineage3_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mosquito_Lineage3_species$species<-as.numeric(Mosquito_Lineage3_species$species)
Mosquito_Lineage3_species$freq<-as.numeric(as.character(Mosquito_Lineage3_species$freq))

Mosquito_Lineage3_species_present<-subset(Mosquito_Lineage3_species,Mosquito_Lineage3_species$freq>0)
Mosquito_Lineage3_species_present %>% count(species)

ggplot(data=Mosquito_Lineage3_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","red2"))

ggsave("/PATH/TO/MP_linC_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mosquito Lineage 4 ##
Lineage4_p1<-Mosquito_SNVs_2[,c(1,14)]
Lineage4_p1$passage<-1
colnames(Lineage4_p1)<-c("SNVID","freq","passage")
Lineage4_p3<-Mosquito_SNVs_2[,c(1,15)]
Lineage4_p3$passage<-3
colnames(Lineage4_p3)<-c("SNVID","freq","passage")
Lineage4_p4<-Mosquito_SNVs_2[,c(1,16)]
Lineage4_p4$passage<-4
colnames(Lineage4_p4)<-c("SNVID","freq","passage")
Lineage4_p10<-Mosquito_SNVs_2[,c(1,17)]
Lineage4_p10$passage<-10
colnames(Lineage4_p10)<-c("SNVID","freq","passage")

Lineage4<-rbind(Lineage4_p1,Lineage4_p3)
Lineage4<-rbind(Lineage4,Lineage4_p4)
Mosquito_Lineage4_nop11<-rbind(Lineage4,Lineage4_p10)

Mosquito_Lineage4_species<-merge(Mosquito_Lineage4_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mosquito_Lineage4_species$species<-as.numeric(Mosquito_Lineage4_species$species)
Mosquito_Lineage4_species$freq<-as.numeric(as.character(Mosquito_Lineage4_species$freq))

Mosquito_Lineage4_species_present<-subset(Mosquito_Lineage4_species,Mosquito_Lineage4_species$freq>0)
Mosquito_Lineage4_species_present %>% count(species)

ggplot(data=Mosquito_Lineage4_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","red2"))

ggsave("/PATH/TO/MP_linD_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mosquito Lineage 5 ##
Lineage5_p1<-Mosquito_SNVs_2[,c(1,18)]
Lineage5_p1$passage<-1
colnames(Lineage5_p1)<-c("SNVID","freq","passage")
Lineage5_p3<-Mosquito_SNVs_2[,c(1,19)]
Lineage5_p3$passage<-3
colnames(Lineage5_p3)<-c("SNVID","freq","passage")
Lineage5_p4<-Mosquito_SNVs_2[,c(1,20)]
Lineage5_p4$passage<-4
colnames(Lineage5_p4)<-c("SNVID","freq","passage")
Lineage5_p10<-Mosquito_SNVs_2[,c(1,21)]
Lineage5_p10$passage<-10
colnames(Lineage5_p10)<-c("SNVID","freq","passage")

Lineage5<-rbind(Lineage5_p1,Lineage5_p3)
Lineage5<-rbind(Lineage5,Lineage5_p4)
Mosquito_Lineage5_nop11<-rbind(Lineage5,Lineage5_p10)

Mosquito_Lineage5_species<-merge(Mosquito_Lineage5_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mosquito_Lineage5_species$species<-as.numeric(Mosquito_Lineage5_species$species)
Mosquito_Lineage5_species$freq<-as.numeric(as.character(Mosquito_Lineage5_species$freq))

Mosquito_Lineage5_species_present<-subset(Mosquito_Lineage5_species,Mosquito_Lineage5_species$freq>0)
Mosquito_Lineage5_species_present %>% count(species)

ggplot(data=Mosquito_Lineage5_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Mosquito") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","red2"))

ggsave("/PATH/TO/MP_linE_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

#### Species-specific SNVs in mouse passage ####

## Mouse Lineage 1 ##
Lineage1_p1<-Mouse_SNVs_2[,c(1,2)]
Lineage1_p1$passage<-1
colnames(Lineage1_p1)<-c("SNVID","freq","passage")
Lineage1_p3<-Mouse_SNVs_2[,c(1,3)]
Lineage1_p3$passage<-3
colnames(Lineage1_p3)<-c("SNVID","freq","passage")
Lineage1_p4<-Mouse_SNVs_2[,c(1,4)]
Lineage1_p4$passage<-4
colnames(Lineage1_p4)<-c("SNVID","freq","passage")
Lineage1_p10<-Mouse_SNVs_2[,c(1,5)]
Lineage1_p10$passage<-10
colnames(Lineage1_p10)<-c("SNVID","freq","passage")

Lineage1<-rbind(Lineage1_p1,Lineage1_p3)
Lineage1<-rbind(Lineage1,Lineage1_p4)
Mouse_Lineage1_nop11<-rbind(Lineage1,Lineage1_p10)

Mouse_Lineage1_species<-merge(Mouse_Lineage1_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mouse_Lineage1_species$species<-as.numeric(Mouse_Lineage1_species$species)
Mouse_Lineage1_species$freq<-as.numeric(as.character(Mouse_Lineage1_species$freq))

Mouse_Lineage1_species_present<-subset(Mouse_Lineage1_species,Mouse_Lineage1_species$freq>0)
Mouse_Lineage1_species_present %>% count(species)

ggplot(data=Mouse_Lineage1_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","red2"))

ggsave("/PATH/TO/SCP_linA_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mouse Lineage 2 ##
Lineage2_p1<-Mouse_SNVs_2[,c(1,6)]
Lineage2_p1$passage<-1
colnames(Lineage2_p1)<-c("SNVID","freq","passage")
Lineage2_p3<-Mouse_SNVs_2[,c(1,7)]
Lineage2_p3$passage<-3
colnames(Lineage2_p3)<-c("SNVID","freq","passage")
Lineage2_p4<-Mouse_SNVs_2[,c(1,8)]
Lineage2_p4$passage<-4
colnames(Lineage2_p4)<-c("SNVID","freq","passage")
Lineage2_p10<-Mouse_SNVs_2[,c(1,9)]
Lineage2_p10$passage<-10
colnames(Lineage2_p10)<-c("SNVID","freq","passage")

Lineage2<-rbind(Lineage2_p1,Lineage2_p3)
Lineage2<-rbind(Lineage2,Lineage2_p4)
Mouse_Lineage2_nop11<-rbind(Lineage2,Lineage2_p10)

Mouse_Lineage2_species<-merge(Mouse_Lineage2_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mouse_Lineage2_species$species<-as.numeric(Mouse_Lineage2_species$species)
Mouse_Lineage2_species$freq<-as.numeric(as.character(Mouse_Lineage2_species$freq))

Mouse_Lineage2_species_present<-subset(Mouse_Lineage2_species,Mouse_Lineage2_species$freq>0)
Mouse_Lineage2_species_present %>% count(species)

ggplot(data=Mouse_Lineage2_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","red2"))

ggsave("/PATH/TO/SCP_linB_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mouse Lineage 3 ##
Lineage3_p1<-Mouse_SNVs_2[,c(1,10)]
Lineage3_p1$passage<-1
colnames(Lineage3_p1)<-c("SNVID","freq","passage")
Lineage3_p3<-Mouse_SNVs_2[,c(1,11)]
Lineage3_p3$passage<-3
colnames(Lineage3_p3)<-c("SNVID","freq","passage")
Lineage3_p4<-Mouse_SNVs_2[,c(1,12)]
Lineage3_p4$passage<-4
colnames(Lineage3_p4)<-c("SNVID","freq","passage")
Lineage3_p10<-Mouse_SNVs_2[,c(1,13)]
Lineage3_p10$passage<-10
colnames(Lineage3_p10)<-c("SNVID","freq","passage")

Lineage3<-rbind(Lineage3_p1,Lineage3_p3)
Lineage3<-rbind(Lineage3,Lineage3_p4)
Mouse_Lineage3_nop11<-rbind(Lineage3,Lineage3_p10)

Mouse_Lineage3_species<-merge(Mouse_Lineage3_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mouse_Lineage3_species$species<-as.numeric(Mouse_Lineage3_species$species)
Mouse_Lineage3_species$freq<-as.numeric(as.character(Mouse_Lineage3_species$freq))

Mouse_Lineage3_species_present<-subset(Mouse_Lineage3_species,Mouse_Lineage3_species$freq>0)
Mouse_Lineage3_species_present %>% count(species)

ggplot(data=Mouse_Lineage3_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","red2"))

ggsave("/PATH/TO/SCP_linC_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mouse Lineage 4 ##
Lineage4_p1<-Mouse_SNVs_2[,c(1,14)]
Lineage4_p1$passage<-1
colnames(Lineage4_p1)<-c("SNVID","freq","passage")
Lineage4_p3<-Mouse_SNVs_2[,c(1,15)]
Lineage4_p3$passage<-3
colnames(Lineage4_p3)<-c("SNVID","freq","passage")
Lineage4_p4<-Mouse_SNVs_2[,c(1,16)]
Lineage4_p4$passage<-4
colnames(Lineage4_p4)<-c("SNVID","freq","passage")
Lineage4_p10<-Mouse_SNVs_2[,c(1,17)]
Lineage4_p10$passage<-10
colnames(Lineage4_p10)<-c("SNVID","freq","passage")

Lineage4<-rbind(Lineage4_p1,Lineage4_p3)
Lineage4<-rbind(Lineage4,Lineage4_p4)
Mouse_Lineage4_nop11<-rbind(Lineage4,Lineage4_p10)

Mouse_Lineage4_species<-merge(Mouse_Lineage4_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mouse_Lineage4_species$species<-as.numeric(Mouse_Lineage4_species$species)
Mouse_Lineage4_species$freq<-as.numeric(as.character(Mouse_Lineage4_species$freq))

Mouse_Lineage4_species_present<-subset(Mouse_Lineage4_species,Mouse_Lineage4_species$freq>0)
Mouse_Lineage4_species_present %>% count(species)

ggplot(data=Mouse_Lineage4_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","red2"))

ggsave("/PATH/TO/SCP_linD_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Mouse Lineage 5 ##
Lineage5_p1<-Mouse_SNVs_2[,c(1,18)]
Lineage5_p1$passage<-1
colnames(Lineage5_p1)<-c("SNVID","freq","passage")
Lineage5_p3<-Mouse_SNVs_2[,c(1,19)]
Lineage5_p3$passage<-3
colnames(Lineage5_p3)<-c("SNVID","freq","passage")
Lineage5_p4<-Mouse_SNVs_2[,c(1,20)]
Lineage5_p4$passage<-4
colnames(Lineage5_p4)<-c("SNVID","freq","passage")
Lineage5_p10<-Mouse_SNVs_2[,c(1,21)]
Lineage5_p10$passage<-10
colnames(Lineage5_p10)<-c("SNVID","freq","passage")

Lineage5<-rbind(Lineage5_p1,Lineage5_p3)
Lineage5<-rbind(Lineage5,Lineage5_p4)
Mouse_Lineage5_nop11<-rbind(Lineage5,Lineage5_p10)

Mouse_Lineage5_species<-merge(Mouse_Lineage5_nop11, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Mouse_Lineage5_species$species<-as.numeric(Mouse_Lineage5_species$species)
Mouse_Lineage5_species$freq<-as.numeric(as.character(Mouse_Lineage5_species$freq))

Mouse_Lineage5_species_present<-subset(Mouse_Lineage5_species,Mouse_Lineage5_species$freq>0)
Mouse_Lineage5_species_present %>% count(species)

ggplot(data=Mouse_Lineage5_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Mouse") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","red2"))

ggsave("/PATH/TO/SCP_linE_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

#### Species-specific SNVs in alternate passage ####

## Alternate Lineage 1 ##
Lineage1_p1<-Alternate_SNVs_2[,c(1,2)]
Lineage1_p1$passage<-1
colnames(Lineage1_p1)<-c("SNVID","freq","passage")
Lineage1_p2<-Alternate_SNVs_2[,c(1,3)]
Lineage1_p2$passage<-2
colnames(Lineage1_p2)<-c("SNVID","freq","passage")
Lineage1_p3<-Alternate_SNVs_2[,c(1,4)]
Lineage1_p3$passage<-3
colnames(Lineage1_p3)<-c("SNVID","freq","passage")
Lineage1_p4<-Alternate_SNVs_2[,c(1,5)]
Lineage1_p4$passage<-4
colnames(Lineage1_p4)<-c("SNVID","freq","passage")
Lineage1_p5<-Alternate_SNVs_2[,c(1,6)]
Lineage1_p5$passage<-5
colnames(Lineage1_p5)<-c("SNVID","freq","passage")
Lineage1_p6<-Alternate_SNVs_2[,c(1,7)]
Lineage1_p6$passage<-6
colnames(Lineage1_p6)<-c("SNVID","freq","passage")
Lineage1_p7<-Alternate_SNVs_2[,c(1,8)]
Lineage1_p7$passage<-7
colnames(Lineage1_p7)<-c("SNVID","freq","passage")
Lineage1_p8<-Alternate_SNVs_2[,c(1,9)]
Lineage1_p8$passage<-8
colnames(Lineage1_p8)<-c("SNVID","freq","passage")
Lineage1_p9<-Alternate_SNVs_2[,c(1,10)]
Lineage1_p9$passage<-9
colnames(Lineage1_p9)<-c("SNVID","freq","passage")
Lineage1_p10<-Alternate_SNVs_2[,c(1,11)]
Lineage1_p10$passage<-10
colnames(Lineage1_p10)<-c("SNVID","freq","passage")

Alternate_Lineage1<-do.call("rbind", list(Lineage1_p1,Lineage1_p2,Lineage1_p3,
                                          Lineage1_p4,Lineage1_p5,Lineage1_p6,
                                          Lineage1_p7,Lineage1_p8,Lineage1_p9,
                                          Lineage1_p10))

Alternate_Lineage1_species<-merge(Alternate_Lineage1, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Alternate_Lineage1_species[is.na(Alternate_Lineage1_species)]<-4
Alternate_Lineage1_species$freq<-as.numeric(Alternate_Lineage1_species$freq)

Alternate_Lineage1_species_present<-subset(Alternate_Lineage1_species,Alternate_Lineage1_species$freq>0)
Alternate_Lineage1_species_present %>% count(species)

ggplot(data=Alternate_Lineage1_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage A - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1,10), breaks=seq(1,10,by=1), labels = c("1","2","3","4","5","6","7","8","9","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","dodgerblue3","red2","forestgreen"))

ggsave("/PATH/TO/AP_linA_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Alternate Lineage 2 ##
Lineage2_p1<-Alternate_SNVs_2[,c(1,12)]
Lineage2_p1$passage<-1
colnames(Lineage2_p1)<-c("SNVID","freq","passage")
Lineage2_p2<-Alternate_SNVs_2[,c(1,13)]
Lineage2_p2$passage<-2
colnames(Lineage2_p2)<-c("SNVID","freq","passage")
Lineage2_p3<-Alternate_SNVs_2[,c(1,14)]
Lineage2_p3$passage<-3
colnames(Lineage2_p3)<-c("SNVID","freq","passage")

Alternate_Lineage2<-do.call("rbind", list(Lineage2_p1,Lineage2_p2,Lineage2_p3))

Alternate_Lineage2$freq<-as.numeric(Alternate_Lineage2$freq)
Alternate_Lineage2$passage<-as.numeric(Alternate_Lineage2$passage)

Alternate_Lineage2_species<-merge(Alternate_Lineage2, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Alternate_Lineage2_species[is.na(Alternate_Lineage2_species)]<-4

Alternate_Lineage2_species_present<-subset(Alternate_Lineage2_species,Alternate_Lineage2_species$freq>0)
Alternate_Lineage2_species_present %>% count(species)

ggplot(data=Alternate_Lineage2_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage B - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1,10), breaks=seq(1,10,by=1), labels = c("1","2","3","4","5","6","7","8","9","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","dodgerblue3","red2","forestgreen"))

ggsave("/PATH/TO/AP_linB_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Alternate Lineage 3 ##
Lineage3_p1<-Alternate_SNVs_2[,c(1,15)]
Lineage3_p1$passage<-1
colnames(Lineage3_p1)<-c("SNVID","freq","passage")
Lineage3_p2<-Alternate_SNVs_2[,c(1,16)]
Lineage3_p2$passage<-2
colnames(Lineage3_p2)<-c("SNVID","freq","passage")
Lineage3_p3<-Alternate_SNVs_2[,c(1,17)]
Lineage3_p3$passage<-3
colnames(Lineage3_p3)<-c("SNVID","freq","passage")

Alternate_Lineage3<-do.call("rbind", list(Lineage3_p1,Lineage3_p2,Lineage3_p3))

Alternate_Lineage3$freq<-as.numeric(Alternate_Lineage3$freq)
Alternate_Lineage3$passage<-as.numeric(Alternate_Lineage3$passage)

Alternate_Lineage3_species<-merge(Alternate_Lineage3, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Alternate_Lineage3_species[is.na(Alternate_Lineage3_species)]<-4

Alternate_Lineage3_species_present<-subset(Alternate_Lineage3_species,Alternate_Lineage3_species$freq>0)
Alternate_Lineage3_species_present %>% count(species)

Alternate_Lineage3_species<-subset(Alternate_Lineage3_species,Alternate_Lineage3_species$passage!=2)

ggplot(data=Alternate_Lineage3_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage C - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1,10), breaks=seq(1,10,by=1), labels = c("1","2","3","4","5","6","7","8","9","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","dodgerblue3","red2","forestgreen"))

ggsave("/PATH/TO/AP_linC_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Alternate Lineage 4 ##
Lineage4_p1<-Alternate_SNVs_2[,c(1,18)]
Lineage4_p1$passage<-1
colnames(Lineage4_p1)<-c("SNVID","freq","passage")
Lineage4_p2<-Alternate_SNVs_2[,c(1,19)]
Lineage4_p2$passage<-2
colnames(Lineage4_p2)<-c("SNVID","freq","passage")
Lineage4_p3<-Alternate_SNVs_2[,c(1,20)]
Lineage4_p3$passage<-3
colnames(Lineage4_p3)<-c("SNVID","freq","passage")
Lineage4_p4<-Alternate_SNVs_2[,c(1,21)]
Lineage4_p4$passage<-4
colnames(Lineage4_p4)<-c("SNVID","freq","passage")
Lineage4_p5<-Alternate_SNVs_2[,c(1,22)]
Lineage4_p5$passage<-5
colnames(Lineage4_p5)<-c("SNVID","freq","passage")
Lineage4_p6<-Alternate_SNVs_2[,c(1,23)]
Lineage4_p6$passage<-6
colnames(Lineage4_p6)<-c("SNVID","freq","passage")
Lineage4_p7<-Alternate_SNVs_2[,c(1,24)]
Lineage4_p7$passage<-7
colnames(Lineage4_p7)<-c("SNVID","freq","passage")
Lineage4_p8<-Alternate_SNVs_2[,c(1,25)]
Lineage4_p8$passage<-8
colnames(Lineage4_p8)<-c("SNVID","freq","passage")
Lineage4_p9<-Alternate_SNVs_2[,c(1,26)]
Lineage4_p9$passage<-9
colnames(Lineage4_p9)<-c("SNVID","freq","passage")
Lineage4_p10<-Alternate_SNVs_2[,c(1,27)]
Lineage4_p10$passage<-10
colnames(Lineage4_p10)<-c("SNVID","freq","passage")

Alternate_Lineage4<-do.call("rbind", list(Lineage4_p1,Lineage4_p2,Lineage4_p3,
                                          Lineage4_p4,Lineage4_p5,Lineage4_p6,
                                          Lineage4_p7,Lineage4_p8,Lineage4_p9,
                                          Lineage4_p10))

Alternate_Lineage4$freq<-as.numeric(Alternate_Lineage4$freq)
Alternate_Lineage4$passage<-as.numeric(Alternate_Lineage4$passage)

Alternate_Lineage4_species<-merge(Alternate_Lineage4, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Alternate_Lineage4_species[is.na(Alternate_Lineage4_species)]<-4

Alternate_Lineage4_species_present<-subset(Alternate_Lineage4_species,Alternate_Lineage4_species$freq>0)
Alternate_Lineage4_species_present %>% count(species)

ggplot(data=Alternate_Lineage4_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage D - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1,10), breaks=seq(1,10,by=1), labels = c("1","2","3","4","5","6","7","8","9","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","dodgerblue3","red2","forestgreen"))

ggsave("/PATH/TO/AP_linD_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)

## Alternate Lineage 5 ##
Lineage5_p1<-Alternate_SNVs_2[,c(1,28)]
Lineage5_p1$passage<-1
colnames(Lineage5_p1)<-c("SNVID","freq","passage")
Lineage5_p2<-Alternate_SNVs_2[,c(1,29)]
Lineage5_p2$passage<-2
colnames(Lineage5_p2)<-c("SNVID","freq","passage")
Lineage5_p3<-Alternate_SNVs_2[,c(1,30)]
Lineage5_p3$passage<-3
colnames(Lineage5_p3)<-c("SNVID","freq","passage")
Lineage5_p4<-Alternate_SNVs_2[,c(1,31)]
Lineage5_p4$passage<-4
colnames(Lineage5_p4)<-c("SNVID","freq","passage")
Lineage5_p5<-Alternate_SNVs_2[,c(1,32)]
Lineage5_p5$passage<-5
colnames(Lineage5_p5)<-c("SNVID","freq","passage")
Lineage5_p6<-Alternate_SNVs_2[,c(1,33)]
Lineage5_p6$passage<-6
colnames(Lineage5_p6)<-c("SNVID","freq","passage")
Lineage5_p7<-Alternate_SNVs_2[,c(1,34)]
Lineage5_p7$passage<-7
colnames(Lineage5_p7)<-c("SNVID","freq","passage")
Lineage5_p8<-Alternate_SNVs_2[,c(1,35)]
Lineage5_p8$passage<-8
colnames(Lineage5_p8)<-c("SNVID","freq","passage")
Lineage5_p9<-Alternate_SNVs_2[,c(1,36)]
Lineage5_p9$passage<-9
colnames(Lineage5_p9)<-c("SNVID","freq","passage")
Lineage5_p10<-Alternate_SNVs_2[,c(1,37)]
Lineage5_p10$passage<-10
colnames(Lineage5_p10)<-c("SNVID","freq","passage")

Alternate_Lineage5<-do.call("rbind", list(Lineage5_p1,Lineage5_p2,Lineage5_p3,
                                          Lineage5_p4,Lineage5_p5,Lineage5_p6,
                                          Lineage5_p7,Lineage5_p8,Lineage5_p9,
                                          Lineage5_p10))

Alternate_Lineage5$freq<-as.numeric(Alternate_Lineage5$freq)
Alternate_Lineage5$passage<-as.numeric(Alternate_Lineage5$passage)

Alternate_Lineage5_species<-merge(Alternate_Lineage5, Mosq_mouse_SNVs, by = "SNVID", all.x = TRUE)
Alternate_Lineage5_species[is.na(Alternate_Lineage5_species)]<-4

Alternate_Lineage5_species_present<-subset(Alternate_Lineage5_species,Alternate_Lineage5_species$freq>0)
Alternate_Lineage5_species_present %>% count(species)

Alternate_Lineage5_species<-subset(Alternate_Lineage5_species,Alternate_Lineage5_species$passage!=2)

ggplot(data=Alternate_Lineage5_species, aes(x=passage, y=freq, group=SNVID)) +
  geom_line(aes(color=as.factor(species)), show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("ZIKV lineage E - Alternate") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(limits = c(1,10), breaks=seq(1,10,by=1), labels = c("1","2","3","4","5","6","7","8","9","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","dodgerblue3","red2","forestgreen"))

ggsave("/PATH/TO/AP_linE_SNVs_overtime_line_SPECIFICSNVS_1percent.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


#### COMPOSITE Species-specific ####

## Mosquito ##
Mosquito_Species<-Mosquito_Lineage1_species
Mosquito_Species$lin2species<-Mosquito_Lineage2_species$species
Mosquito_Species$lin3species<-Mosquito_Lineage3_species$species
Mosquito_Species$lin4species<-Mosquito_Lineage4_species$species
Mosquito_Species$lin5species<-Mosquito_Lineage5_species$species
Mosquito_Species$Sumspecies<-Mosquito_Species$species+Mosquito_Species$lin2species+Mosquito_Species$lin3species+Mosquito_Species$lin4species+Mosquito_Species$lin5species
Mosquito_Species$Mosqspecific<-ifelse(Mosquito_Species$Sumspecies<6,1,0)

MosqL1<-subset(Mosquito_Lineage1_species,Mosquito_Lineage1_species$species==1)
MosqL2<-subset(Mosquito_Lineage2_species,Mosquito_Lineage2_species$species==1)
MosqL3<-subset(Mosquito_Lineage3_species,Mosquito_Lineage3_species$species==1)
MosqL4<-subset(Mosquito_Lineage4_species,Mosquito_Lineage4_species$species==1)
MosqL5<-subset(Mosquito_Lineage5_species,Mosquito_Lineage5_species$species==1)

Mosq_mosqspecific<-cbind(MosqL1,MosqL2)
Mosq_mosqspecific<-Mosq_mosqspecific[,c(1,3,2,8)]
Mosq_mosqspecific<-cbind(Mosq_mosqspecific,MosqL3)
Mosq_mosqspecific<-Mosq_mosqspecific[,c(1:4,6)]
Mosq_mosqspecific<-cbind(Mosq_mosqspecific,MosqL4)
Mosq_mosqspecific<-Mosq_mosqspecific[,c(1:5,7)]
Mosq_mosqspecific<-cbind(Mosq_mosqspecific,MosqL5)
Mosq_mosqspecific<-Mosq_mosqspecific[,c(1:6,8)]
Mosq_mosqspecific$SNVID<-as.factor(Mosq_mosqspecific$SNVID)
Mosq_mosqspecific$passage<-as.numeric(Mosq_mosqspecific$passage)
Mosq_mosqspecific<-Mosq_mosqspecific[order(Mosq_mosqspecific$SNVID, Mosq_mosqspecific$passage),]
colnames(Mosq_mosqspecific)<-c("SNVID","passage","freq1","freq2","freq3","freq4","freq5")

Mosq_mosqspecific[Mosq_mosqspecific<0.001]<-as.numeric("NA")

Mosq_mosqspecific<-Mosq_mosqspecific %>% mutate(meanfreq=rowMeans(na.rm=TRUE,select(.,freq1,freq2,freq3,freq4,freq5)))
Mosq_mosqspecific[is.na(Mosq_mosqspecific)]<-0
Mosq_mosqspecific$prezero<-ifelse(((dplyr::lead(Mosq_mosqspecific$meanfreq)>0) | (Mosq_mosqspecific$meanfreq==0 & dplyr::lag(Mosq_mosqspecific$passage)<10 & dplyr::lag(Mosq_mosqspecific$meanfreq)>0)),1,0)

Mosq_mosqspecific<-subset(Mosq_mosqspecific,Mosq_mosqspecific$meanfreq>0 | Mosq_mosqspecific$prezero>0)
Mosq_mosqspecific<-Mosq_mosqspecific[,c(1,8,2)]
colnames(Mosq_mosqspecific)<-c("SNVID","freq","passage")
Mosq_mosqspecific$species<-1

MosqL1_shared<-subset(Mosquito_Lineage1_species,Mosquito_Lineage1_species$species==3)
MosqL1_shared$SNVIDPASS<-paste(MosqL1_shared$SNVID,MosqL1_shared$passage,sep=",")
MosqL2_shared<-subset(Mosquito_Lineage2_species,Mosquito_Lineage2_species$species==3)
MosqL2_shared$SNVIDPASS<-paste(MosqL2_shared$SNVID,MosqL2_shared$passage,sep=",")
MosqL3_shared<-subset(Mosquito_Lineage3_species,Mosquito_Lineage3_species$species==3)
MosqL3_shared$SNVIDPASS<-paste(MosqL3_shared$SNVID,MosqL3_shared$passage,sep=",")
MosqL4_shared<-subset(Mosquito_Lineage4_species,Mosquito_Lineage4_species$species==3)
MosqL4_shared$SNVIDPASS<-paste(MosqL4_shared$SNVID,MosqL4_shared$passage,sep=",")
MosqL5_shared<-subset(Mosquito_Lineage5_species,Mosquito_Lineage5_species$species==3)
MosqL5_shared$SNVIDPASS<-paste(MosqL5_shared$SNVID,MosqL5_shared$passage,sep=",")

Mosq_shared<-merge(MosqL1_shared, MosqL2_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mosq_shared<-Mosq_shared[,c(1,3,9)]
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

nrow(Mosq_shared)/4

Mosq_shared[Mosq_shared < 0.001]<-as.numeric("NA")

Mosq_shared<-Mosq_shared %>% mutate(freq=rowMeans(na.rm=TRUE,select(., L1Freq,L2Freq,L3Freq,L4Freq,L5Freq)))
Mosq_shared[is.na(Mosq_shared)]<-0
Mosq_shared<-cSplit(Mosq_shared, 'SNVIDPASS', sep=",", type.convert=FALSE)
Mosq_shared$SNVIDPASS_2<-as.numeric(Mosq_shared$SNVIDPASS_2)
Mosq_shared<-Mosq_shared[order(Mosq_shared$SNVIDPASS_1, Mosq_shared$SNVIDPASS_2),]
Mosq_shared$prezero<-ifelse(dplyr::lead(Mosq_shared$freq)>0,1,0)
Mosq_shared<-subset(Mosq_shared,Mosq_shared$freq>0 | Mosq_shared$prezero>0)
Mosq_shared<-Mosq_shared[,c(7,6,8)]
colnames(Mosq_shared)<-c("SNVID","freq","passage")
Mosq_shared$species<-0

Mosq_SPECIFIC<-rbind(Mosq_mosqspecific,Mosq_shared)
Mosq_SPECIFIC$passage<-as.numeric(Mosq_SPECIFIC$passage)
Mosq_SPECIFIC$SNVID<-as.factor(Mosq_SPECIFIC$SNVID)
Mosq_SPECIFIC$species<-as.factor(Mosq_SPECIFIC$species)

ggplot(data=Mosq_SPECIFIC, aes(x=passage, y=freq)) +
  geom_path(aes(group = SNVID, color = species)) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), legend.title = element_blank(), legend.key = element_blank()) +
  ggtitle("ZIKV Mosquito Composite") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red2","grey0"), labels = c("Mosq + Mouse","Mosq-only"))

ggsave("/PATH/TO/MP_COMPOSITE_SNVs_overtime_line_MosquitoSpeciesSpecific.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


## Mice ##
Mouse_Species<-Mouse_Lineage1_species
Mouse_Species$lin2species<-Mouse_Lineage2_species$species
Mouse_Species$lin3species<-Mouse_Lineage3_species$species
Mouse_Species$lin4species<-Mouse_Lineage4_species$species
Mouse_Species$lin5species<-Mouse_Lineage5_species$species
Mouse_Species$Sumspecies<-Mouse_Species$species+Mouse_Species$lin2species+Mouse_Species$lin3species+Mouse_Species$lin4species+Mouse_Species$lin5species
Mouse_Species$Mousespecific<-ifelse(Mouse_Species$Sumspecies<6,1,0)

MouseL1<-subset(Mouse_Lineage1_species,Mouse_Lineage1_species$species==2)
MouseL2<-subset(Mouse_Lineage2_species,Mouse_Lineage2_species$species==2)
MouseL3<-subset(Mouse_Lineage3_species,Mouse_Lineage3_species$species==2)
MouseL4<-subset(Mouse_Lineage4_species,Mouse_Lineage4_species$species==2)
MouseL5<-subset(Mouse_Lineage5_species,Mouse_Lineage5_species$species==2)

Mouse_Mousespecific<-cbind(MouseL1,MouseL2)
Mouse_Mousespecific<-Mouse_Mousespecific[,c(1,3,2,8)]
Mouse_Mousespecific<-cbind(Mouse_Mousespecific,MouseL3)
Mouse_Mousespecific<-Mouse_Mousespecific[,c(1:4,6)]
Mouse_Mousespecific<-cbind(Mouse_Mousespecific,MouseL4)
Mouse_Mousespecific<-Mouse_Mousespecific[,c(1:5,7)]
Mouse_Mousespecific<-cbind(Mouse_Mousespecific,MouseL5)
Mouse_Mousespecific<-Mouse_Mousespecific[,c(1:6,8)]
Mouse_Mousespecific$SNVID<-as.factor(Mouse_Mousespecific$SNVID)
Mouse_Mousespecific$passage<-as.numeric(Mouse_Mousespecific$passage)
Mouse_Mousespecific<-Mouse_Mousespecific[order(Mouse_Mousespecific$SNVID, Mouse_Mousespecific$passage),]
colnames(Mouse_Mousespecific)<-c("SNVID","passage","freq1","freq2","freq3","freq4","freq5")

Mouse_Mousespecific[Mouse_Mousespecific < 0.001]<-as.numeric("NA")

Mouse_Mousespecific<-Mouse_Mousespecific %>% mutate(meanfreq=rowMeans(na.rm=TRUE,select(.,freq1,freq2,freq3,freq4,freq5)))
Mouse_Mousespecific[is.na(Mouse_Mousespecific)]<-0
Mouse_Mousespecific$prezero<-ifelse(((dplyr::lead(Mouse_Mousespecific$meanfreq)>0) | (Mouse_Mousespecific$meanfreq==0 & dplyr::lag(Mouse_Mousespecific$passage)<10 & dplyr::lag(Mouse_Mousespecific$meanfreq)>0)),1,0)

Mouse_Mousespecific<-subset(Mouse_Mousespecific,Mouse_Mousespecific$meanfreq>0 | Mouse_Mousespecific$prezero>0)
Mouse_Mousespecific<-Mouse_Mousespecific[,c(1,8,2)]
colnames(Mouse_Mousespecific)<-c("SNVID","freq","passage")
Mouse_Mousespecific$species<-1

MouseL1_shared<-subset(Mouse_Lineage1_species,Mouse_Lineage1_species$species==3)
MouseL1_shared$SNVIDPASS<-paste(MouseL1_shared$SNVID,MouseL1_shared$passage,sep=",")
MouseL2_shared<-subset(Mouse_Lineage2_species,Mouse_Lineage2_species$species==3)
MouseL2_shared$SNVIDPASS<-paste(MouseL2_shared$SNVID,MouseL2_shared$passage,sep=",")
MouseL3_shared<-subset(Mouse_Lineage3_species,Mouse_Lineage3_species$species==3)
MouseL3_shared$SNVIDPASS<-paste(MouseL3_shared$SNVID,MouseL3_shared$passage,sep=",")
MouseL4_shared<-subset(Mouse_Lineage4_species,Mouse_Lineage4_species$species==3)
MouseL4_shared$SNVIDPASS<-paste(MouseL4_shared$SNVID,MouseL4_shared$passage,sep=",")
MouseL5_shared<-subset(Mouse_Lineage5_species,Mouse_Lineage5_species$species==3)
MouseL5_shared$SNVIDPASS<-paste(MouseL5_shared$SNVID,MouseL5_shared$passage,sep=",")

Mouse_shared<-merge(MouseL1_shared, MouseL2_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mouse_shared<-Mouse_shared[,c(1,3,9)]
colnames(Mouse_shared)<-c("SNVIDPASS","L1Freq","L2Freq")
Mouse_shared<-merge(Mouse_shared, MouseL3_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mouse_shared<-Mouse_shared[,c(1:3,5)]
colnames(Mouse_shared)<-c("SNVIDPASS","L1Freq","L2Freq","L3Freq")
Mouse_shared<-merge(Mouse_shared, MouseL4_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mouse_shared<-Mouse_shared[,c(1:4,6)]
colnames(Mouse_shared)<-c("SNVIDPASS","L1Freq","L2Freq","L3Freq","L4Freq")
Mouse_shared<-merge(Mouse_shared, MouseL5_shared, by = "SNVIDPASS", all.x = TRUE, all.y = TRUE)
Mouse_shared<-Mouse_shared[,c(1:5,7)]
colnames(Mouse_shared)<-c("SNVIDPASS","L1Freq","L2Freq","L3Freq","L4Freq","L5Freq")

Mouse_shared[Mouse_shared < 0.001]<-as.numeric("NA")

Mouse_shared<-Mouse_shared %>% mutate(freq=rowMeans(na.rm=TRUE,select(., L1Freq,L2Freq,L3Freq,L4Freq,L5Freq)))
Mouse_shared[is.na(Mouse_shared)]<-0
Mouse_shared<-cSplit(Mouse_shared, 'SNVIDPASS', sep=",", type.convert=FALSE)
Mouse_shared$SNVIDPASS_2<-as.numeric(Mouse_shared$SNVIDPASS_2)
Mouse_shared<-Mouse_shared[order(Mouse_shared$SNVIDPASS_1, Mouse_shared$SNVIDPASS_2),]
Mouse_shared$prezero<-ifelse(dplyr::lead(Mouse_shared$freq)>0,1,0)
Mouse_shared<-subset(Mouse_shared,Mouse_shared$freq>0 | Mouse_shared$prezero>0)
Mouse_shared<-Mouse_shared[,c(7,6,8)]
colnames(Mouse_shared)<-c("SNVID","freq","passage")
Mouse_shared$species<-0

Mouse_SPECIFIC<-rbind(Mouse_Mousespecific,Mouse_shared)
Mouse_SPECIFIC$passage<-as.numeric(Mouse_SPECIFIC$passage)
Mouse_SPECIFIC$SNVID<-as.factor(Mouse_SPECIFIC$SNVID)
Mouse_SPECIFIC$species<-as.factor(Mouse_SPECIFIC$species)

ggplot(data=Mouse_SPECIFIC, aes(x=passage, y=freq)) +
  geom_path(aes(group = SNVID, color = species)) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),  legend.title = element_blank(), legend.key = element_blank()) +
  ggtitle("ZIKV Mouse Composite") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","","3","4","","","","","","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("red2","dodgerblue3"), labels = c("Mosq + Mouse","Mouse-only"))

ggsave("/PATH/TO/SCP_COMPOSITE_SNVs_overtime_line_MouseSpeciesSpecific.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


## Alternate ##
Alternate_Lineage1_species$SNVIDPASS<-paste(Alternate_Lineage1_species$SNVID,Alternate_Lineage1_species$passage,sep=",")
Alternate_Lineage2_species$SNVIDPASS<-paste(Alternate_Lineage2_species$SNVID,Alternate_Lineage2_species$passage,sep=",")
Alternate_Lineage3_species$SNVIDPASS<-paste(Alternate_Lineage3_species$SNVID,Alternate_Lineage3_species$passage,sep=",")
Alternate_Lineage4_species$SNVIDPASS<-paste(Alternate_Lineage4_species$SNVID,Alternate_Lineage4_species$passage,sep=",")
Alternate_Lineage5_species$SNVIDPASS<-paste(Alternate_Lineage5_species$SNVID,Alternate_Lineage5_species$passage,sep=",")

Alternate_Species<-Alternate_Lineage1_species
Alternate_Species<-merge(Alternate_Species,Alternate_Lineage2_species, by="SNVIDPASS", all.x=TRUE, all.y=TRUE)
Alternate_Species<-Alternate_Species[,c(1,3,9,7,13)]
colnames(Alternate_Species)<-c("SNVIDPASS","freq1","freq2","lin1species","lin2species")
Alternate_Species<-merge(Alternate_Species,Alternate_Lineage3_species, by="SNVIDPASS", all.x=TRUE, all.y=TRUE)
Alternate_Species<-Alternate_Species[,c(1:3,7,4:5,11)]
colnames(Alternate_Species)<-c("SNVIDPASS","freq1","freq2","freq3","lin1species","lin2species","lin3species")
Alternate_Species<-merge(Alternate_Species,Alternate_Lineage4_species, by="SNVIDPASS", all.x=TRUE, all.y=TRUE)
Alternate_Species<-Alternate_Species[,c(1:4,9,5:7,13)]
colnames(Alternate_Species)<-c("SNVIDPASS","freq1","freq2","freq3","freq4","lin1species","lin2species","lin3species","lin4species")
Alternate_Species<-merge(Alternate_Species,Alternate_Lineage5_species, by="SNVIDPASS", all.x=TRUE, all.y=TRUE)
Alternate_Species<-Alternate_Species[,c(1:5,11,6:9,15)]
colnames(Alternate_Species)<-c("SNVIDPASS","freq1","freq2","freq3","freq4","freq5","lin1species","lin2species","lin3species","lin4species","lin5species")

Alternate_Species[Alternate_Species < 0.001]<-as.numeric("NA")
Alternate_Species<-Alternate_Species %>% mutate(meanfreq=rowMeans(na.rm=TRUE,select(.,freq1,freq2,freq3,freq4,freq5)))
Alternate_Species[is.na(Alternate_Species)]<-0
Alternate_Species$Species<-apply(Alternate_Species[, 4:8], 1, max)

Alternate_Species<-cSplit(Alternate_Species, 'SNVIDPASS', sep=",", type.convert=FALSE)
colnames(Alternate_Species)<-c("freq1","freq2","freq3","freq4","freq5","lin1species","lin2species","lin3species","lin4species","lin5species","freq","species","SNVID","passage")
Alternate_Species$passage<-as.numeric(Alternate_Species$passage)
Alternate_Species<-Alternate_Species[with(Alternate_Species, order(SNVID, passage)),]
Alternate_Species$prezero<-ifelse(((dplyr::lead(Alternate_Species$freq)>0) | (Alternate_Species$freq==0 & dplyr::lag(Alternate_Species$passage)<10 & dplyr::lag(Alternate_Species$freq)>0)),1,0)

Alternate_Species<-subset(Alternate_Species,Alternate_Species$freq>0 | Alternate_Species$prezero>0)
Alternate_Species<-Alternate_Species[,c(13,14,11,12)]
colnames(Alternate_Species)<-c("SNVID","passage","freq","species")

Alternate_Species$passage<-as.numeric(Alternate_Species$passage)
Alternate_Species$SNVID<-as.factor(Alternate_Species$SNVID)
Alternate_Species$species<-as.factor(Alternate_Species$species)

ggplot(data=Alternate_Species, aes(x=passage, y=freq)) +
  geom_path(aes(group = SNVID, color = species)) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.50, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), legend.title = element_blank(), legend.key = element_blank()) +
  ggtitle("ZIKV Alternate Composite") +
  labs(x="Passage", y="SNV Frequency") +
  scale_x_continuous(breaks=seq(1,10,by=1), labels = c("1","2","3","4","5","6","7","8","9","10"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("grey0","dodgerblue3","red2","forestgreen"), labels = c("Mosq + Alt","Mouse + Alt","Mosq + Mouse + Alt","Alt-only"))
  

ggsave("/PATH/TO/AP_COMPOSITE_SNVs_overtime_line_AlternateSpeciesSpecific.png",plot=last_plot(),device = png(), scale = 1, dpi = 300)


Mouse_Species_Unique<-Mouse_Species %>% distinct(SNVID, .keep_all = TRUE)
Mouse_Species_Unique %>% count(species)

Alternate_LinOnly_Unique<-Alternate_LineageOnly %>% distinct(SNVID, .keep_all = TRUE)
Alternate_LinOnly_Unique %>% count(Onelinonly)

Mosquito_Species_Unique<-Mosquito_Species %>% distinct(SNVID, .keep_all = TRUE)
Mosquito_Species_Unique %>% count(species)
Alternate_Species_Unique<-Alternate_Species %>% distinct(SNVID, .keep_all = TRUE)
Alternate_Species_Unique %>% count(species)

