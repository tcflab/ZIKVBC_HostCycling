list.of.packages <- c("plyr","dplyr","splitstackshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("plyr")
library("dplyr")
library("splitstackshape")

#### Enter the sample name and output directory for the sample
args<-commandArgs()
homedir<-args[6]
name<-args[7]

outputsample<-paste(homedir,name,sep="/")

### Import vcf tables for each reference aligned vcf ###
rep1<-paste(name,"_1.vcf",sep="")
refvcf1_rep1<-paste(homedir,name,sep="/")
refvcf2_rep1<-paste(refvcf1_rep1,"_rep1_v3.1/reference_aligned/lofreq_reference_",sep="")
refvcf3_rep1<-paste(refvcf2_rep1,rep1,sep="")
vcf_ref_1<-read.table(file=refvcf3_rep1, sep="\t", header=FALSE)
vcf_ref_1<-cSplit(vcf_ref_1,"V8",sep=";", type.convert=FALSE)
vcf_ref_1<-cSplit(vcf_ref_1, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
vcf_ref_1<-vcf_ref_1[,c(1,2,3,4,5,6,7,11,13)]
colnames(vcf_ref_1)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFODP","INFOAF")
vcf_ref_1$INFODP<-as.numeric(vcf_ref_1$INFODP)
vcf_ref_1$INFOAF<-as.numeric(vcf_ref_1$INFOAF)
vcf_ref_1<-subset(vcf_ref_1,vcf_ref_1$INFODP>299)

rep2<-paste(name,"_2.vcf",sep="")
refvcf1_rep2<-paste(homedir,name,sep="/")
refvcf2_rep2<-paste(refvcf1_rep2,"_rep2_v3.1/reference_aligned/lofreq_reference_",sep="")
refvcf3_rep2<-paste(refvcf2_rep2,rep2,sep="")
vcf_ref_2<-read.table(file=refvcf3_rep2, sep="\t", header=FALSE)
vcf_ref_2<-cSplit(vcf_ref_2,"V8",sep=";", type.convert=FALSE)
vcf_ref_2<-cSplit(vcf_ref_2, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
vcf_ref_2<-vcf_ref_2[,c(1,2,3,4,5,6,7,11,13)]
colnames(vcf_ref_2)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFODP","INFOAF")
vcf_ref_2$INFODP<-as.numeric(vcf_ref_2$INFODP)
vcf_ref_2$INFOAF<-as.numeric(vcf_ref_2$INFOAF)
vcf_ref_2<-subset(vcf_ref_2,vcf_ref_2$INFODP>299)

## Merge VCF files by SNV ##
vcf_ref_1$POSREFALT<-paste(vcf_ref_1$POS,vcf_ref_1$REF,sep="")
vcf_ref_1$POSREFALT<-paste(vcf_ref_1$POSREFALT,vcf_ref_1$ALT,sep="")
vcf_ref_2$POSREFALT<-paste(vcf_ref_2$POS,vcf_ref_2$REF,sep="")
vcf_ref_2$POSREFALT<-paste(vcf_ref_2$POSREFALT,vcf_ref_2$ALT,sep="")

ref_1and2<-merge(vcf_ref_1, vcf_ref_2, by = "POSREFALT", all.x = TRUE, all.y = TRUE)
ref_1and2$INFODP.x[is.na(ref_1and2$INFODP.x)]<-0
ref_1and2$INFODP.y[is.na(ref_1and2$INFODP.y)]<-0
ref_1and2$INFOAF.x[is.na(ref_1and2$INFOAF.x)]<-0
ref_1and2$INFOAF.y[is.na(ref_1and2$INFOAF.y)]<-0
ref_1and2$POS.x[is.na(ref_1and2$POS.x)]<-0
ref_1and2$POS.y[is.na(ref_1and2$POS.y)]<-0
ref_1and2$QUAL.x[is.na(ref_1and2$QUAL.x)]<-0
ref_1and2$QUAL.y[is.na(ref_1and2$QUAL.y)]<-0

## Calculate coverage ##
ref_1and2$INFODP.x<-as.numeric(ref_1and2$INFODP.x)
ref_1and2$INFODP.y<-as.numeric(ref_1and2$INFODP.y)
ref_1and2$INFODP.xy<-ifelse(ref_1and2$INFODP.x>0 & ref_1and2$INFODP.y>0,ceiling((ref_1and2$INFODP.x+ref_1and2$INFODP.y)/2),pmax(ref_1and2$INFODP.x,ref_1and2$INFODP.y))

## Calculate allele frequency ##
ref_1and2$INFOAF.x<-as.numeric(ref_1and2$INFOAF.x)
ref_1and2$INFOAF.y<-as.numeric(ref_1and2$INFOAF.y)
ref_1and2$INFOAF.xy<-ifelse(ref_1and2$INFOAF.x>0 & ref_1and2$INFOAF.y>0, (ref_1and2$INFOAF.x+ref_1and2$INFOAF.y)/2, 0)
ref_1and2<-subset(ref_1and2, ref_1and2$INFOAF.xy>0)

## Fill in blanks in merged table
ref_1and2$FILTER.x<-"PASS"
ref_1and2$REF.x[is.na(ref_1and2$REF.x)]<-as.character(ref_1and2$REF.y[is.na(ref_1and2$REF.x)])
ref_1and2$ALT.x[is.na(ref_1and2$ALT.x)]<-as.character(ref_1and2$ALT.y[is.na(ref_1and2$ALT.x)])
ref_1and2$ID.x[is.na(ref_1and2$ID.x)]<-as.character(ref_1and2$ID.y[is.na(ref_1and2$ID.x)])
ref_1and2$CHROM.x[is.na(ref_1and2$CHROM.x)]<-as.character(ref_1and2$CHROM.y[is.na(ref_1and2$CHROM.x)])
ref_1and2$POS.x[ref_1and2$POS.x==0]<-as.numeric(ref_1and2$POS.y[ref_1and2$POS.x==0])
ref_1and2$QUAL.x[ref_1and2$QUAL.x==0]<-as.numeric(ref_1and2$QUAL.y[ref_1and2$QUAL.x==0])

## Reformat table in VCF format
con_1and2<-ref_1and2
con_1and2$INFOAF.xy<-ifelse(ref_1and2$INFOAF.xy>=0.5,1-ref_1and2$INFOAF.xy,ref_1and2$INFOAF.xy)
con_1and2$REF.x<-ifelse(ref_1and2$INFOAF.xy>=0.5,as.character(ref_1and2$ALT.x),as.character(ref_1and2$REF.x))
con_1and2$ALT.x<-ifelse(ref_1and2$INFOAF.xy>=0.5,as.character(ref_1and2$REF.x),as.character(ref_1and2$ALT.x))

ref_1and2$INFOAF.xy<-paste("AF=",ref_1and2$INFOAF.xy,sep="")
ref_1and2$INFODP.xy<-paste("DP=",ref_1and2$INFODP.xy,sep="")
ref_1and2$INFO<-paste(ref_1and2$INFODP.xy,ref_1and2$INFOAF.xy,sep=";")

con_1and2$INFOAF.xy<-paste("AF=",con_1and2$INFOAF.xy,sep="")
con_1and2$INFODP.xy<-paste("DP=",con_1and2$INFODP.xy,sep="")
con_1and2$INFO<-paste(con_1and2$INFODP.xy,con_1and2$INFOAF.xy,sep=";")

ref_1and2<-ref_1and2[,c(2,3,4,5,6,7,8,22)]
colnames(ref_1and2)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

con_1and2<-con_1and2[,c(2,3,4,5,6,7,8,22)]
colnames(con_1and2)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

ref_1and2<-ref_1and2[order(ref_1and2$POS),]
con_1and2<-con_1and2[order(con_1and2$POS),]

## Write VCF file for reference alignment ##
vcf1<-paste("rep1and2_refalign",name,sep="")
vcf2<-paste(vcf1,".vcf",sep="")
output_ref<-paste(homedir,vcf2,sep="/")
write.table(ref_1and2, file=output_ref, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

## Write VCF file for consensus alignment ##
vcf11<-paste("rep1and2_conalign",name,sep="")
vcf12<-paste(vcf11,".vcf",sep="")
output_ref<-paste(homedir,vcf12,sep="/")
write.table(con_1and2, file=output_ref, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")







