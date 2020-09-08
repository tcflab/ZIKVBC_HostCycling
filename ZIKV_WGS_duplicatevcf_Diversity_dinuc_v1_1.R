## load package libraries
list.of.packages <- c("diverse","plyr","dplyr","RcppRoll","splitstackshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("diverse")
library("plyr")
library("dplyr")
library("RcppRoll")
library("splitstackshape")

## import arguments
args<-commandArgs()
homedir<-args[6]
name<-args[7]

outputdir1<-paste(homedir,name,sep="/")
outputdir<-paste(outputdir1,"_combined",sep="")

#### Import duplicate combined vcf tables for each reference aligned vcf ####
ref1<-paste("rep1and2_refalign",name,sep="")
ref2<-paste(ref1,".vcf",sep="")
ref3<-paste(homedir,"Duplicate_refalign_vcf",sep="/")
ref4<-paste(ref3,ref2,sep="/")
refaligned<-read.table(file=ref4, sep="\t", header=FALSE)
refaligned<-cSplit(refaligned,"V8",sep=";", type.convert=FALSE)
refaligned<-cSplit(refaligned, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
refaligned<-refaligned[,c(1,2,3,4,5,6,7,9,11)]
colnames(refaligned)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFODP","INFOAF")

#### Import duplicate combined vcf tables for each consensus aligned vcf ####
con1<-paste("rep1and2_conalign",name,sep="")
con2<-paste(con1,".vcf",sep="")
con3<-paste(homedir,"Duplicate_conalign_vcf",sep="/")
con4<-paste(con3,con2,sep="/")
conaligned<-read.table(file=con4, sep="\t", header=FALSE)
conaligned<-cSplit(conaligned,"V8",sep=";", type.convert=FALSE)
conaligned<-cSplit(conaligned, c("V8_1","V8_2"), sep="=", type.convert=FALSE)
conaligned<-conaligned[,c(1,2,3,4,5,6,7,9,11)]
colnames(conaligned)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFODP","INFOAF")

rep1dir<-paste(name,"rep1",sep="_")
rep1dir2<-paste(rep1dir,"_v3.1",sep="")
inputdir<-paste(homedir,rep1dir2,sep="/")
rep1refdir<-paste(inputdir,"reference_aligned",sep="/")
rep1condir<-paste(inputdir,"consensus_aligned",sep="/")
rep1name<-paste(name,"1",sep="_")

rep2dir<-paste(name,"rep2",sep="_")
rep2dir2<-paste(rep2dir,"_v3.1",sep="")
inputdir<-paste(homedir,rep2dir2,sep="/")
rep2refdir<-paste(inputdir,"reference_aligned",sep="/")
rep2condir<-paste(inputdir,"consensus_aligned",sep="/")
rep2name<-paste(name,"2",sep="_")

## import reference_aligned ntcounts for rep 1
table_r_1<-paste("ntcounts_reference_",rep1name,sep="")
table_r_1<-paste(table_r_1,".csv",sep="")
table_r_1<-paste(rep1refdir,table_r_1,sep="/")
ntcounts_r_1<-read.csv(table_r_1, header=TRUE)

## import consensus_aligned ntcounts for rep1 
table_c_1<-paste("ntcounts_consensus_",rep1name,sep="")
table_c_1<-paste(table_c_1,".csv",sep="")
table_c_1<-paste(rep1condir,table_c_1,sep="/")
ntcounts_c_1<-read.csv(table_c_1, header=TRUE)

## import reference_aligned ntcounts for rep 2
table_r_2<-paste("ntcounts_reference_",rep2name,sep="")
table_r_2<-paste(table_r_2,".csv",sep="")
table_r_2<-paste(rep2refdir,table_r_2,sep="/")
ntcounts_r_2<-read.csv(table_r_2, header=TRUE)

## import consensus_aligned ntcounts for rep2
table_c_2<-paste("ntcounts_consensus_",rep2name,sep="")
table_c_2<-paste(table_c_2,".csv",sep="")
table_c_2<-paste(rep2condir,table_c_2,sep="/")
ntcounts_c_2<-read.csv(table_c_2, header=TRUE)

## rename column headers
columns<-c("pos","ref","reads_all","mismatches","deletions","insertions","A","C","T","G")
ntcounts_r_1<-ntcounts_r_1[columns]
ntcounts_c_1<-ntcounts_c_1[columns]
ntcounts_r_2<-ntcounts_r_2[columns]
ntcounts_c_2<-ntcounts_c_2[columns]
colnames(ntcounts_r_1)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G")
colnames(ntcounts_c_1)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G")
colnames(ntcounts_r_2)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G")
colnames(ntcounts_c_2)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G")

ntcounts_r<-merge(ntcounts_r_1,ntcounts_r_2,by="position")
ntcounts_c<-merge(ntcounts_c_1,ntcounts_c_2,by="position")

ntcounts_r$coverage.xy<-ntcounts_r$coverage.x+ntcounts_r$coverage.y
ntcounts_r$mismatches.xy<-ntcounts_r$mismatches.x+ntcounts_r$mismatches.y
ntcounts_r$deletions.xy<-ntcounts_r$deletions.x+ntcounts_r$deletions.y
ntcounts_r$insertions.xy<-ntcounts_r$insertions.x+ntcounts_r$insertions.y
ntcounts_r$A.xy<-ntcounts_r$A.x+ntcounts_r$A.y
ntcounts_r$C.xy<-ntcounts_r$C.x+ntcounts_r$C.y
ntcounts_r$T.xy<-ntcounts_r$T.x+ntcounts_r$T.y
ntcounts_r$G.xy<-ntcounts_r$G.x+ntcounts_r$G.y

ntcounts_c$coverage.xy<-ntcounts_c$coverage.x+ntcounts_c$coverage.y
ntcounts_c$mismatches.xy<-ntcounts_c$mismatches.x+ntcounts_c$mismatches.y
ntcounts_c$deletions.xy<-ntcounts_c$deletions.x+ntcounts_c$deletions.y
ntcounts_c$insertions.xy<-ntcounts_c$insertions.x+ntcounts_c$insertions.y
ntcounts_c$A.xy<-ntcounts_c$A.x+ntcounts_c$A.y
ntcounts_c$C.xy<-ntcounts_c$C.x+ntcounts_c$C.y
ntcounts_c$T.xy<-ntcounts_c$T.x+ntcounts_c$T.y
ntcounts_c$G.xy<-ntcounts_c$G.x+ntcounts_c$G.y

ntcounts_r<-ntcounts_r[,c(1,2,20,21,22,23,24,25,26,27)]
ntcounts_c<-ntcounts_c[,c(1,2,20,21,22,23,24,25,26,27)]
colnames(ntcounts_r)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G")
colnames(ntcounts_c)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G")


######   DIVERSITY METRICS     #######

## generate absolute SNV columns
nucs<-c("A","C","G","T")

for (i in nucs){
  ntcounts_c<-cbind(ntcounts_c,as.data.frame(ifelse(ntcounts_c$reference==i,0,
                                                ifelse(ntcounts_c$reference!=i,ntcounts_c[,i],
                                                       ifelse(ntcounts_c$reference=="",0,"na")))))
  
}

colnames(ntcounts_c)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G","A.SNV","C.SNV","G.SNV","T.SNV")
ntcounts_c$A.SNV<-as.numeric(ntcounts_c$A.SNV)
ntcounts_c$C.SNV<-as.numeric(ntcounts_c$C.SNV)
ntcounts_c$T.SNV<-as.numeric(ntcounts_c$T.SNV)
ntcounts_c$G.SNV<-as.numeric(ntcounts_c$G.SNV)

## generate SNV frequency columns
ntcounts_c$snv.freq<-(ntcounts_c$mismatches/ntcounts_c$coverage)
ntcounts_c$snvfreqA<-(ntcounts_c$A.SNV/ntcounts_c$coverage)
ntcounts_c$snvfreqT<-(ntcounts_c$T.SNV/ntcounts_c$coverage)
ntcounts_c$snvfreqC<-(ntcounts_c$C.SNV/ntcounts_c$coverage)
ntcounts_c$snvfreqG<-(ntcounts_c$G.SNV/ntcounts_c$coverage)

## generate Squared Deviation column
ntcounts_c$sq.dev<-((ntcounts_c$snv.freq)^2)

## generate unique SNV column
ntcounts_c$unique.A.SNV<-ifelse(ntcounts_c$A.SNV>0,1,0)
ntcounts_c$unique.C.SNV<-ifelse(ntcounts_c$C.SNV>0,1,0)
ntcounts_c$unique.G.SNV<-ifelse(ntcounts_c$G.SNV>0,1,0)
ntcounts_c$unique.T.SNV<-ifelse(ntcounts_c$T.SNV>0,1,0)

## subset data for coverage >= 300
highcov<-subset(ntcounts_c,ntcounts_c$coverage>299)
highcov$snv.freq<-((highcov$A.SNV + highcov$C.SNV + highcov$G.SNV + highcov$T.SNV)/highcov$coverage)
highcov$sq.dev<-((highcov$snv.freq)^2)

## subset highcov for positions called by lofreq
variant<-subset(highcov, highcov$position %in% conaligned$POS)

## calculate root mean squared deviation
rmsd<-sqrt((sum(variant$sq.dev)/nrow(highcov)))
rmsd

## calculate shannon entropy
nucvar<-subset(variant,select=c("A","C","G","T"))
nucvar<-as.matrix(nucvar)

entropy.nucvar<-diversity(nucvar, type="entropy")
entropy.mean<-sum(entropy.nucvar$entropy)/nrow(highcov)
#entropy.sd<-sd(entropy.nucvar$entropy)

## calculate Gini-Simpson index
ginisimpson<-diversity(nucvar, type="gini-simpson")
gs.mean<-sum(ginisimpson$gini.simpson.C)/nrow(highcov)
#ginisimpson$mean.gs.C<-gs.mean
#gs.sd<-sd(ginisimpson$gini.simpson.C)

## calculate mutation frequency (SNV/10,000 nt sequenced)
mut.freq<-(sum(variant$A.SNV,variant$C.SNV,variant$G.SNV,variant$T.SNV)/sum(highcov$coverage))*10000
mut.freq
unique.mut.freq<-(sum(variant$unique.A.SNV,variant$unique.C.SNV,variant$unique.G.SNV,variant$unique.T.SNV)/sum(highcov$coverage))*10000
unique.mut.freq

## merge diversity indices and write file
diversity<-c(mut.freq,unique.mut.freq,rmsd,entropy.mean,gs.mean)
diversity<-t(diversity)
colnames(diversity)<-c("Mut.Freq.per.10K","Unique.Mut.Freq","RMSD","Mean.Shannon.Entropy", "Mean.Gini.Simpson")
diversity<-as.data.frame(diversity)
diversityname<-paste(name,"_combined_diversity_metrics.csv",sep="")
outputdiversity<-paste(outputdir,diversityname,sep="/")
write.csv(diversity, file = outputdiversity, row.names = FALSE)


######   Dinucleotide frequencies   #######

####   ONLY CALLED SNPs   ####

refseq<-read.csv(file="~/Desktop/PRVABC59.csv",header=TRUE)

### Generate dinucleotide sites

dinucs<-c("ApA","ApC","ApG","ApT","CpA","CpC","CpG","CpT","GpA","GpC","GpG","GpT","TpA","TpC","TpG","TpT")

for (i in dinucs){
  assign(i,subset(refseq,refseq$reference==substring(i,1,1) & dplyr::lead(refseq$reference)==substring(i,3,3) | refseq$reference==substring(i,3,3) & dplyr::lag(refseq$reference)==substring(i,1,1)))
}

perbase1<-ntcounts_r
perbase2<-subset(perbase1, perbase1$position %in% refaligned$POS)
colnames(perbase2)<-c("pos","reference","reads_all","mismatches","deletions","insertions","A","C","T","G")
perbase3<-subset(perbase1, perbase1$coverage>299 & (perbase1$position<4002 | perbase1$position>4035))
colnames(perbase3)<-c("pos","reference","reads_all","mismatches","deletions","insertions","A","C","T","G")

### ApA
ApAlosses1<-subset(perbase3,(perbase3$reference=="A" & dplyr::lead(perbase3$reference=="A")) | (perbase3$reference=="A" & dplyr::lag(perbase3$reference=="A")))
ApAlosses1$losses<-ApAlosses1$reads_all-ApAlosses1$A
ApAgains1<-subset(perbase3,perbase3$reference!="A" & dplyr::lead(perbase3$reference=="A"))
ApAgains2<-subset(perbase3,perbase3$reference!="A" & dplyr::lag(perbase3$reference=="A"))

### ApC
ApClosses1<-subset(perbase3,perbase3$reference=="A" & dplyr::lead(perbase3$reference=="C"))
ApClosses1$losses<-ApClosses1$reads_all-ApClosses1$A
ApClosses2<-subset(perbase3,perbase3$reference=="C" & dplyr::lag(perbase3$reference=="A"))
ApClosses2$losses<-ApClosses2$reads_all-ApClosses2$C
ApCgains1<-subset(perbase3,perbase3$reference!="A" & dplyr::lead(perbase3$reference=="C"))
ApCgains2<-subset(perbase3,perbase3$reference!="C" & dplyr::lag(perbase3$reference=="A"))

### ApT
ApTlosses1<-subset(perbase3,perbase3$reference=="A" & dplyr::lead(perbase3$reference=="T"))
ApTlosses1$losses<-ApTlosses1$reads_all-ApTlosses1$A
ApTlosses2<-subset(perbase3,perbase3$reference=="T" & dplyr::lag(perbase3$reference=="A"))
ApTlosses2$losses<-ApTlosses2$reads_all-ApTlosses2$T
ApTgains1<-subset(perbase3,perbase3$reference!="A" & dplyr::lead(perbase3$reference=="T"))
ApTgains2<-subset(perbase3,perbase3$reference!="T" & dplyr::lag(perbase3$reference=="A"))

### ApG
ApGlosses1<-subset(perbase3,perbase3$reference=="A" & dplyr::lead(perbase3$reference=="G"))
ApGlosses1$losses<-ApGlosses1$reads_all-ApGlosses1$A
ApGlosses2<-subset(perbase3,perbase3$reference=="G" & dplyr::lag(perbase3$reference=="A"))
ApGlosses2$losses<-ApGlosses2$reads_all-ApGlosses2$G
ApGgains1<-subset(perbase3,perbase3$reference!="A" & dplyr::lead(perbase3$reference=="G"))
ApGgains2<-subset(perbase3,perbase3$reference!="G" & dplyr::lag(perbase3$reference=="A"))

### CpC
CpClosses1<-subset(perbase3,(perbase3$reference=="C" & dplyr::lead(perbase3$reference=="C")) | (perbase3$reference=="C" & dplyr::lag(perbase3$reference=="C")))
CpClosses1$losses<-CpClosses1$reads_all-CpClosses1$C
CpCgains1<-subset(perbase3,perbase3$reference!="C" & dplyr::lead(perbase3$reference=="C"))
CpCgains2<-subset(perbase3,perbase3$reference!="C" & dplyr::lag(perbase3$reference=="C"))

### CpA
CpAlosses1<-subset(perbase3,perbase3$reference=="C" & dplyr::lead(perbase3$reference=="A"))
CpAlosses1$losses<-CpAlosses1$reads_all-CpAlosses1$C
CpAlosses2<-subset(perbase3,perbase3$reference=="A" & dplyr::lag(perbase3$reference=="C"))
CpAlosses2$losses<-CpAlosses2$reads_all-CpAlosses2$A
CpAgains1<-subset(perbase3,perbase3$reference!="C" & dplyr::lead(perbase3$reference=="A"))
CpAgains2<-subset(perbase3,perbase3$reference!="A" & dplyr::lag(perbase3$reference=="C"))

### CpT
CpTlosses1<-subset(perbase3,perbase3$reference=="C" & dplyr::lead(perbase3$reference=="T"))
CpTlosses1$losses<-CpTlosses1$reads_all-CpTlosses1$C
CpTlosses2<-subset(perbase3,perbase3$reference=="T" & dplyr::lag(perbase3$reference=="C"))
CpTlosses2$losses<-CpTlosses2$reads_all-CpTlosses2$T
CpTgains1<-subset(perbase3,perbase3$reference!="C" & dplyr::lead(perbase3$reference=="T"))
CpTgains2<-subset(perbase3,perbase3$reference!="T" & dplyr::lag(perbase3$reference=="C"))

### CpG
CpGlosses1<-subset(perbase3,perbase3$reference=="C" & dplyr::lead(perbase3$reference=="G"))
CpGlosses1$losses<-CpGlosses1$reads_all-CpGlosses1$C
CpGlosses2<-subset(perbase3,perbase3$reference=="G" & dplyr::lag(perbase3$reference=="C"))
CpGlosses2$losses<-CpGlosses2$reads_all-CpGlosses2$G
CpGgains1<-subset(perbase3,perbase3$reference!="C" & dplyr::lead(perbase3$reference=="G"))
CpGgains2<-subset(perbase3,perbase3$reference!="G" & dplyr::lag(perbase3$reference=="C"))

### GpG
GpGlosses1<-subset(perbase3,(perbase3$reference=="G" & dplyr::lead(perbase3$reference=="G")) | (perbase3$reference=="G" & dplyr::lag(perbase3$reference=="G")))
GpGlosses1$losses<-GpGlosses1$reads_all-GpGlosses1$G
GpGgains1<-subset(perbase3,perbase3$reference!="G" & dplyr::lead(perbase3$reference=="G"))
GpGgains2<-subset(perbase3,perbase3$reference!="G" & dplyr::lag(perbase3$reference=="G"))

### GpC
GpClosses1<-subset(perbase3,perbase3$reference=="G" & dplyr::lead(perbase3$reference=="C"))
GpClosses1$losses<-GpClosses1$reads_all-GpClosses1$G
GpClosses2<-subset(perbase3,perbase3$reference=="C" & dplyr::lag(perbase3$reference=="G"))
GpClosses2$losses<-GpClosses2$reads_all-GpClosses2$C
GpCgains1<-subset(perbase3,perbase3$reference!="G" & dplyr::lead(perbase3$reference=="C"))
GpCgains2<-subset(perbase3,perbase3$reference!="C" & dplyr::lag(perbase3$reference=="G"))

### GpT
GpTlosses1<-subset(perbase3,perbase3$reference=="G" & dplyr::lead(perbase3$reference=="T"))
GpTlosses1$losses<-GpTlosses1$reads_all-GpTlosses1$G
GpTlosses2<-subset(perbase3,perbase3$reference=="T" & dplyr::lag(perbase3$reference=="G"))
GpTlosses2$losses<-GpTlosses2$reads_all-GpTlosses2$T
GpTgains1<-subset(perbase3,perbase3$reference!="G" & dplyr::lead(perbase3$reference=="T"))
GpTgains2<-subset(perbase3,perbase3$reference!="T" & dplyr::lag(perbase3$reference=="G"))

### GpA
GpAlosses1<-subset(perbase3,perbase3$reference=="G" & dplyr::lead(perbase3$reference=="A"))
GpAlosses1$losses<-GpAlosses1$reads_all-GpAlosses1$G
GpAlosses2<-subset(perbase3,perbase3$reference=="A" & dplyr::lag(perbase3$reference=="G"))
GpAlosses2$losses<-GpAlosses2$reads_all-GpAlosses2$A
GpAgains1<-subset(perbase3,perbase3$reference!="G" & dplyr::lead(perbase3$reference=="A"))
GpAgains2<-subset(perbase3,perbase3$reference!="A" & dplyr::lag(perbase3$reference=="G"))

### TpA
TpAlosses1<-subset(perbase3,perbase3$reference=="T" & dplyr::lead(perbase3$reference=="A"))
TpAlosses1$losses<-TpAlosses1$reads_all-TpAlosses1$T
TpAlosses2<-subset(perbase3,perbase3$reference=="A" & dplyr::lag(perbase3$reference=="T"))
TpAlosses2$losses<-TpAlosses2$reads_all-TpAlosses2$A
TpAgains1<-subset(perbase3,perbase3$reference!="T" & dplyr::lead(perbase3$reference=="A"))
TpAgains2<-subset(perbase3,perbase3$reference!="A" & dplyr::lag(perbase3$reference=="T"))

### TpC
TpClosses1<-subset(perbase3,perbase3$reference=="T" & dplyr::lead(perbase3$reference=="C"))
TpClosses1$losses<-TpClosses1$reads_all-TpClosses1$T
TpClosses2<-subset(perbase3,perbase3$reference=="C" & dplyr::lag(perbase3$reference=="T"))
TpClosses2$losses<-TpClosses2$reads_all-TpClosses2$C
TpCgains1<-subset(perbase3,perbase3$reference!="T" & dplyr::lead(perbase3$reference=="C"))
TpCgains2<-subset(perbase3,perbase3$reference!="C" & dplyr::lag(perbase3$reference=="T"))

### TpG
TpGlosses1<-subset(perbase3,perbase3$reference=="T" & dplyr::lead(perbase3$reference=="G"))
TpGlosses1$losses<-TpGlosses1$reads_all-TpGlosses1$T
TpGlosses2<-subset(perbase3,perbase3$reference=="G" & dplyr::lag(perbase3$reference=="T"))
TpGlosses2$losses<-TpGlosses2$reads_all-TpGlosses2$G
TpGgains1<-subset(perbase3,perbase3$reference!="T" & dplyr::lead(perbase3$reference=="G"))
TpGgains2<-subset(perbase3,perbase3$reference!="G" & dplyr::lag(perbase3$reference=="T"))

### TpT
TpTlosses1<-subset(perbase3,(perbase3$reference=="T" & dplyr::lead(perbase3$reference=="T")) | (perbase3$reference=="T" & dplyr::lag(perbase3$reference=="T")))
TpTlosses1$losses<-TpTlosses1$reads_all-TpTlosses1$T
TpTgains1<-subset(perbase3,perbase3$reference!="T" & dplyr::lead(perbase3$reference=="T"))
TpTgains2<-subset(perbase3,perbase3$reference!="T" & dplyr::lag(perbase3$reference=="T"))

#### Calculate gains and losses #####
ApC_gains_persite<-(sum(ApCgains1$A)+sum(ApCgains2$C))/(sum(ApCgains1$reads_all)+sum(ApCgains2$reads_all))*4000
ApC_losses_persite<-(sum(ApClosses1$losses)+sum(ApClosses2$losses))/(sum(ApClosses1$reads_all)+sum(ApClosses2$reads_all))*4000
ApC_gains_total<-ApC_gains_persite*(nrow(ApCgains1)+nrow(ApCgains2))
ApC_losses_total<-ApC_losses_persite*(nrow(ApClosses1)+nrow(ApClosses2))
ApC_ratio<-ApC_gains_persite/ApC_losses_persite

ApG_gains_persite<-(sum(ApGgains1$A)+sum(ApGgains2$G))/(sum(ApGgains1$reads_all)+sum(ApGgains2$reads_all))*4000
ApG_losses_persite<-(sum(ApGlosses1$losses)+sum(ApGlosses2$losses))/(sum(ApGlosses1$reads_all)+sum(ApGlosses2$reads_all))*4000
ApG_gains_total<-ApG_gains_persite*(nrow(ApGgains1)+nrow(ApGgains2))
ApG_losses_total<-ApG_losses_persite*(nrow(ApGlosses1)+nrow(ApGlosses2))
ApG_ratio<-ApG_gains_persite/ApG_losses_persite

ApT_gains_persite<-(sum(ApTgains1$A)+sum(ApTgains2$T))/(sum(ApTgains1$reads_all)+sum(ApTgains2$reads_all))*4000
ApT_losses_persite<-(sum(ApTlosses1$losses)+sum(ApTlosses2$losses))/(sum(ApTlosses1$reads_all)+sum(ApTlosses2$reads_all))*4000
ApT_gains_total<-ApT_gains_persite*(nrow(ApTgains1)+nrow(ApTgains2))
ApT_losses_total<-ApT_losses_persite*(nrow(ApTlosses1)+nrow(ApTlosses2))
ApT_ratio<-ApT_gains_persite/ApT_losses_persite

ApA_gains_persite<-(sum(ApAgains1$A)+sum(ApAgains2$A))/(sum(ApAgains1$reads_all)+sum(ApAgains2$reads_all))*4000
ApA_losses_persite<-sum(ApAlosses1$losses)/sum(ApAlosses1$reads_all)*4000
ApA_gains_total<-ApA_gains_persite*(nrow(ApAgains1)+nrow(ApAgains2))
ApA_losses_total<-ApA_losses_persite*nrow(ApAlosses1)
ApA_ratio<-ApA_gains_persite/ApA_losses_persite

CpA_gains_persite<-(sum(CpAgains1$C)+sum(CpAgains2$A))/(sum(CpAgains1$reads_all)+sum(CpAgains2$reads_all))*4000
CpA_losses_persite<-(sum(CpAlosses1$losses)+sum(CpAlosses2$losses))/(sum(CpAlosses1$reads_all)+sum(CpAlosses2$reads_all))*4000
CpA_gains_total<-CpA_gains_persite*(nrow(CpAgains1)+nrow(CpAgains2))
CpA_losses_total<-CpA_losses_persite*(nrow(CpAlosses1)+nrow(CpAlosses2))
CpA_ratio<-CpA_gains_persite/CpA_losses_persite

CpG_gains_persite<-(sum(CpGgains1$C)+sum(CpGgains2$G))/(sum(CpGgains1$reads_all)+sum(CpGgains2$reads_all))*4000
CpG_losses_persite<-(sum(CpGlosses1$losses)+sum(CpGlosses2$losses))/(sum(CpGlosses1$reads_all)+sum(CpGlosses2$reads_all))*4000
CpG_gains_total<-CpG_gains_persite*(nrow(CpGgains1)+nrow(CpGgains2))
CpG_losses_total<-CpG_losses_persite*(nrow(CpGlosses1)+nrow(CpGlosses2))
CpG_ratio<-CpG_gains_persite/CpG_losses_persite

CpT_gains_persite<-(sum(CpTgains1$C)+sum(CpTgains2$T))/(sum(CpTgains1$reads_all)+sum(CpTgains2$reads_all))*4000
CpT_losses_persite<-(sum(CpTlosses1$losses)+sum(CpTlosses2$losses))/(sum(CpTlosses1$reads_all)+sum(CpTlosses2$reads_all))*4000
CpT_gains_total<-CpT_gains_persite*(nrow(CpTgains1)+nrow(CpTgains2))
CpT_losses_total<-CpT_losses_persite*(nrow(CpTlosses1)+nrow(CpTlosses2))
CpT_ratio<-CpT_gains_persite/CpT_losses_persite

CpC_gains_persite<-(sum(CpCgains1$C)+sum(CpCgains2$C))/(sum(CpCgains1$reads_all)+sum(CpCgains2$reads_all))*4000
CpC_losses_persite<-sum(CpClosses1$losses)/sum(CpClosses1$reads_all)*4000
CpC_gains_total<-CpC_gains_persite*(nrow(CpCgains1)+nrow(CpCgains2))
CpC_losses_total<-CpC_losses_persite*nrow(CpClosses1)
CpC_ratio<-CpC_gains_persite/CpC_losses_persite

GpA_gains_persite<-(sum(GpAgains1$G)+sum(GpAgains2$A))/(sum(GpAgains1$reads_all)+sum(GpAgains2$reads_all))*4000
GpA_losses_persite<-(sum(GpAlosses1$losses)+sum(GpAlosses2$losses))/(sum(GpAlosses1$reads_all)+sum(GpAlosses2$reads_all))*4000
GpA_gains_total<-GpA_gains_persite*(nrow(GpAgains1)+nrow(GpAgains2))
GpA_losses_total<-GpA_losses_persite*(nrow(GpAlosses1)+nrow(GpAlosses2))
GpA_ratio<-GpA_gains_persite/GpA_losses_persite

GpC_gains_persite<-(sum(GpCgains1$G)+sum(GpCgains2$C))/(sum(GpCgains1$reads_all)+sum(GpCgains2$reads_all))*4000
GpC_losses_persite<-(sum(GpClosses1$losses)+sum(GpClosses2$losses))/(sum(GpClosses1$reads_all)+sum(GpClosses2$reads_all))*4000
GpC_gains_total<-GpC_gains_persite*(nrow(GpCgains1)+nrow(GpCgains2))
GpC_losses_total<-GpC_losses_persite*(nrow(GpClosses1)+nrow(GpClosses2))
GpC_ratio<-GpC_gains_persite/GpC_losses_persite

GpT_gains_persite<-(sum(GpTgains1$G)+sum(GpTgains2$T))/(sum(GpTgains1$reads_all)+sum(GpTgains2$reads_all))*4000
GpT_losses_persite<-(sum(GpTlosses1$losses)+sum(GpTlosses2$losses))/(sum(GpTlosses1$reads_all)+sum(GpTlosses2$reads_all))*4000
GpT_gains_total<-GpT_gains_persite*(nrow(GpTgains1)+nrow(GpTgains2))
GpT_losses_total<-GpT_losses_persite*(nrow(GpTlosses1)+nrow(GpTlosses2))
GpT_ratio<-GpT_gains_persite/GpT_losses_persite

GpG_gains_persite<-(sum(GpGgains1$G)+sum(GpGgains2$G))/(sum(GpGgains1$reads_all)+sum(GpGgains2$reads_all))*4000
GpG_losses_persite<-sum(GpGlosses1$losses)/sum(GpGlosses1$reads_all)*4000
GpG_gains_total<-GpG_gains_persite*(nrow(GpGgains1)+nrow(GpGgains2))
GpG_losses_total<-GpG_losses_persite*nrow(GpGlosses1)
GpG_ratio<-GpG_gains_persite/GpG_losses_persite

TpT_gains_persite<-(sum(TpTgains1$T)+sum(TpTgains2$T))/(sum(TpTgains1$reads_all)+sum(TpTgains2$reads_all))*4000
TpT_losses_persite<-sum(TpTlosses1$losses)/sum(TpTlosses1$reads_all)*4000
TpT_gains_total<-TpT_gains_persite*(nrow(TpTgains1)+nrow(TpTgains2))
TpT_losses_total<-TpT_losses_persite*nrow(TpTlosses1)
TpT_ratio<-TpT_gains_persite/TpT_losses_persite

TpA_gains_persite<-(sum(TpAgains1$T)+sum(TpAgains2$A))/(sum(TpAgains1$reads_all)+sum(TpAgains2$reads_all))*4000
TpA_losses_persite<-(sum(TpAlosses1$losses)+sum(TpAlosses2$losses))/(sum(TpAlosses1$reads_all)+sum(TpAlosses2$reads_all))*4000
TpA_gains_total<-TpA_gains_persite*(nrow(TpAgains1)+nrow(TpAgains2))
TpA_losses_total<-TpA_losses_persite*(nrow(TpAlosses1)+nrow(TpAlosses2))
TpA_ratio<-TpA_gains_persite/TpA_losses_persite

TpC_gains_persite<-(sum(TpCgains1$T)+sum(TpCgains2$C))/(sum(TpCgains1$reads_all)+sum(TpCgains2$reads_all))*4000
TpC_losses_persite<-(sum(TpClosses1$losses)+sum(TpClosses2$losses))/(sum(TpClosses1$reads_all)+sum(TpClosses2$reads_all))*4000
TpC_gains_total<-TpC_gains_persite*(nrow(TpCgains1)+nrow(TpCgains2))
TpC_losses_total<-TpC_losses_persite*(nrow(TpClosses1)+nrow(TpClosses2))
TpC_ratio<-TpC_gains_persite/TpC_losses_persite

TpG_gains_persite<-(sum(TpGgains1$T)+sum(TpGgains2$G))/(sum(TpGgains1$reads_all)+sum(TpGgains2$reads_all))*4000
TpG_losses_persite<-(sum(TpGlosses1$losses)+sum(TpGlosses2$losses))/(sum(TpGlosses1$reads_all)+sum(TpGlosses2$reads_all))*4000
TpG_gains_total<-TpG_gains_persite*(nrow(TpGgains1)+nrow(TpGgains2))
TpG_losses_total<-TpG_losses_persite*(nrow(TpGlosses1)+nrow(TpGlosses2))
TpG_ratio<-TpG_gains_persite/TpG_losses_persite

gainspersite<-as.data.frame(c(ApA_gains_persite,ApC_gains_persite,ApG_gains_persite,ApT_gains_persite,
                              CpA_gains_persite,CpC_gains_persite,CpG_gains_persite,CpT_gains_persite,
                              GpA_gains_persite,GpC_gains_persite,GpG_gains_persite,GpT_gains_persite,
                              TpA_gains_persite,TpC_gains_persite,TpG_gains_persite,TpT_gains_persite))
lossespersite<-as.data.frame(c(ApA_losses_persite,ApC_losses_persite,ApG_losses_persite,ApT_losses_persite,
                               CpA_losses_persite,CpC_losses_persite,CpG_losses_persite,CpT_losses_persite,
                               GpA_losses_persite,GpC_losses_persite,GpG_losses_persite,GpT_losses_persite,
                               TpA_losses_persite,TpC_losses_persite,TpG_losses_persite,TpT_losses_persite))
gainstotal<-as.data.frame(c(ApA_gains_total,ApC_gains_total,ApG_gains_total,ApT_gains_total,
                            CpA_gains_total,CpC_gains_total,CpG_gains_total,CpT_gains_total,
                            GpA_gains_total,GpC_gains_total,GpG_gains_total,GpT_gains_total,
                            TpA_gains_total,TpC_gains_total,TpG_gains_total,TpT_gains_total))
lossestotal<-as.data.frame(c(ApA_losses_total,ApC_losses_total,ApG_losses_total,ApT_losses_total,
                             CpA_losses_total,CpC_losses_total,CpG_losses_total,CpT_losses_total,
                             GpA_losses_total,GpC_losses_total,GpG_losses_total,GpT_losses_total,
                             TpA_losses_total,TpC_losses_total,TpG_losses_total,TpT_losses_total))
summary<-cbind(gainspersite,lossespersite)
summary<-cbind(summary,gainstotal)
summary<-cbind(summary,lossestotal)
colnames(summary)<-c("gains_per_site","losses_per_site","gains_total","losses_total")

#### Write dinuc_allreads table ####
dinucname<-paste(name,"_combined_dinuc_metrics_allreads.csv",sep="")
outputdinuc<-paste(outputdir,dinucname,sep="/")
write.csv(summary, file = outputdinuc, row.names = dinucs)


#### Mutational spectrum with all reads VS REFERENCE #####

for (i in nucs){
  ntcounts_r<-cbind(ntcounts_r,as.data.frame(ifelse(ntcounts_r$reference==i,0,
                                                    ifelse(ntcounts_r$reference!=i,ntcounts_r[,i],
                                                           ifelse(ntcounts_r$reference=="",0,"na")))))
  
}

colnames(ntcounts_r)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G","A.SNV","C.SNV","G.SNV","T.SNV")
ntcounts_r$A.SNV<-as.numeric(ntcounts_r$A.SNV)
ntcounts_r$C.SNV<-as.numeric(ntcounts_r$C.SNV)
ntcounts_r$T.SNV<-as.numeric(ntcounts_r$T.SNV)
ntcounts_r$G.SNV<-as.numeric(ntcounts_r$G.SNV)

highcov_r<-subset(ntcounts_r,ntcounts_r$coverage>299 & (ntcounts_r$pos<4007 | ntcounts_r$pos>4030))

## Calculate frequency of specific point mutations
highcov_r$AtoC<-as.numeric(ifelse(highcov_r$reference=="A",highcov_r$C.SNV,0))
highcov_r$AtoG<-as.numeric(ifelse(highcov_r$reference=="A",highcov_r$G.SNV,0))
highcov_r$AtoT<-as.numeric(ifelse(highcov_r$reference=="A",highcov_r$T.SNV,0))
highcov_r$CtoA<-as.numeric(ifelse(highcov_r$reference=="C",highcov_r$A.SNV,0))
highcov_r$CtoG<-as.numeric(ifelse(highcov_r$reference=="C",highcov_r$G.SNV,0))
highcov_r$CtoT<-as.numeric(ifelse(highcov_r$reference=="C",highcov_r$T.SNV,0))
highcov_r$GtoA<-as.numeric(ifelse(highcov_r$reference=="G",highcov_r$A.SNV,0))
highcov_r$GtoC<-as.numeric(ifelse(highcov_r$reference=="G",highcov_r$C.SNV,0))
highcov_r$GtoT<-as.numeric(ifelse(highcov_r$reference=="G",highcov_r$T.SNV,0))
highcov_r$TtoA<-as.numeric(ifelse(highcov_r$reference=="T",highcov_r$A.SNV,0))
highcov_r$TtoC<-as.numeric(ifelse(highcov_r$reference=="T",highcov_r$C.SNV,0))
highcov_r$TtoG<-as.numeric(ifelse(highcov_r$reference=="T",highcov_r$G.SNV,0))

AtoCfreq<-sum(highcov_r$AtoC)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
AtoGfreq<-sum(highcov_r$AtoG)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
AtoTfreq<-sum(highcov_r$AtoT)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
CtoAfreq<-sum(highcov_r$CtoA)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
CtoGfreq<-sum(highcov_r$CtoG)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
CtoTfreq<-sum(highcov_r$CtoT)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
GtoAfreq<-sum(highcov_r$GtoA)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
GtoCfreq<-sum(highcov_r$GtoC)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
GtoTfreq<-sum(highcov_r$GtoT)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
TtoAfreq<-sum(highcov_r$TtoA)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
TtoCfreq<-sum(highcov_r$TtoC)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))
TtoGfreq<-sum(highcov_r$TtoG)/(sum(highcov_r$A.SNV,highcov_r$C.SNV,highcov_r$G.SNV,highcov_r$T.SNV))

## Calculate mutation frequencies for specific subsititutions
highcov_r$Acoverage<-as.numeric(ifelse(highcov_r$reference=="A",highcov_r$coverage,0))
highcov_r$Ccoverage<-as.numeric(ifelse(highcov_r$reference=="C",highcov_r$coverage,0))
highcov_r$Gcoverage<-as.numeric(ifelse(highcov_r$reference=="G",highcov_r$coverage,0))
highcov_r$Tcoverage<-as.numeric(ifelse(highcov_r$reference=="T",highcov_r$coverage,0))

AtoCmutfreq<-sum(highcov_r$AtoC)/(sum(highcov_r$Acoverage))
AtoGmutfreq<-sum(highcov_r$AtoG)/(sum(highcov_r$Acoverage))
AtoTmutfreq<-sum(highcov_r$AtoT)/(sum(highcov_r$Acoverage))
CtoAmutfreq<-sum(highcov_r$CtoA)/(sum(highcov_r$Ccoverage))
CtoGmutfreq<-sum(highcov_r$CtoG)/(sum(highcov_r$Ccoverage))
CtoTmutfreq<-sum(highcov_r$CtoT)/(sum(highcov_r$Ccoverage))
GtoAmutfreq<-sum(highcov_r$GtoA)/(sum(highcov_r$Gcoverage))
GtoCmutfreq<-sum(highcov_r$GtoC)/(sum(highcov_r$Gcoverage))
GtoTmutfreq<-sum(highcov_r$GtoT)/(sum(highcov_r$Gcoverage))
TtoAmutfreq<-sum(highcov_r$TtoA)/(sum(highcov_r$Tcoverage))
TtoCmutfreq<-sum(highcov_r$TtoC)/(sum(highcov_r$Tcoverage))
TtoGmutfreq<-sum(highcov_r$TtoG)/(sum(highcov_r$Tcoverage))

## merge mutation spectrum and write file
spectrum<-c(AtoTfreq,AtoCfreq,AtoGfreq,CtoAfreq,CtoGfreq,CtoTfreq,GtoAfreq,GtoCfreq,GtoTfreq,TtoAfreq,TtoCfreq,TtoGfreq)
spectrum<-t(spectrum)
colnames(spectrum)<-c("AtoU","AtoC","AtoG","CtoA","CtoG","CtoU","GtoA","GtoC","GtoU","UtoA","UtoC","UtoG")
spectrum<-as.data.frame(spectrum)
spectrumname<-paste(name,"_combined_mutational_spectrum_allreads.csv",sep="")
csvfile<-paste(outputdir,spectrumname,sep="/")
write.csv(spectrum, file = csvfile, row.names = FALSE)

## merge mutation spectrum frequencies and write file
spectrum_freq<-c(AtoTmutfreq,AtoCmutfreq,AtoGmutfreq,CtoAmutfreq,CtoGmutfreq,CtoTmutfreq,GtoAmutfreq,GtoCmutfreq,GtoTmutfreq,TtoAmutfreq,TtoCmutfreq,TtoGmutfreq)
spectrum_freq<-t(spectrum_freq)
colnames(spectrum_freq)<-c("A>U","A>C","A>G","C>A","C>G","C>U","G>A","G>C","G>U","U>A","U>C","U>G")
spectrum_freq<-as.data.frame(spectrum_freq)
spectrumname<-paste(name,"_combined_mutational_spectrum_frequencies_allreads.csv",sep="")
csvfile<-paste(outputdir,spectrumname,sep="/")
write.csv(spectrum_freq, file = csvfile, row.names = FALSE)

#### Mutational spectrum with only called SNPs VS REFERENCE #####
Abovethresh<-subset(refaligned,refaligned$INFOAF>=0.003)
highcov_r_called<-subset(highcov_r, highcov_r$position %in% Abovethresh$POS)

highcov_r_called$A.SNVfreq<-highcov_r_called$A.SNV/highcov_r_called$coverage
highcov_r_called$C.SNVfreq<-highcov_r_called$C.SNV/highcov_r_called$coverage
highcov_r_called$T.SNVfreq<-highcov_r_called$T.SNV/highcov_r_called$coverage
highcov_r_called$G.SNVfreq<-highcov_r_called$G.SNV/highcov_r_called$coverage

highcov_r_called$A.SNVcalled<-ifelse(highcov_r_called$A.SNVfreq>=0.003,highcov_r_called$A.SNV,0)
highcov_r_called$C.SNVcalled<-ifelse(highcov_r_called$C.SNVfreq>=0.003,highcov_r_called$C.SNV,0)
highcov_r_called$T.SNVcalled<-ifelse(highcov_r_called$T.SNVfreq>=0.003,highcov_r_called$T.SNV,0)
highcov_r_called$G.SNVcalled<-ifelse(highcov_r_called$G.SNVfreq>=0.003,highcov_r_called$G.SNV,0)

## Calculate frequency of specific point mutations
highcov_r_called$AtoC<-as.numeric(ifelse(highcov_r_called$reference=="A",highcov_r_called$C.SNVcalled,0))
highcov_r_called$AtoG<-as.numeric(ifelse(highcov_r_called$reference=="A",highcov_r_called$G.SNVcalled,0))
highcov_r_called$AtoT<-as.numeric(ifelse(highcov_r_called$reference=="A",highcov_r_called$T.SNVcalled,0))
highcov_r_called$CtoA<-as.numeric(ifelse(highcov_r_called$reference=="C",highcov_r_called$A.SNVcalled,0))
highcov_r_called$CtoG<-as.numeric(ifelse(highcov_r_called$reference=="C",highcov_r_called$G.SNVcalled,0))
highcov_r_called$CtoT<-as.numeric(ifelse(highcov_r_called$reference=="C",highcov_r_called$T.SNVcalled,0))
highcov_r_called$GtoA<-as.numeric(ifelse(highcov_r_called$reference=="G",highcov_r_called$A.SNVcalled,0))
highcov_r_called$GtoC<-as.numeric(ifelse(highcov_r_called$reference=="G",highcov_r_called$C.SNVcalled,0))
highcov_r_called$GtoT<-as.numeric(ifelse(highcov_r_called$reference=="G",highcov_r_called$T.SNVcalled,0))
highcov_r_called$TtoA<-as.numeric(ifelse(highcov_r_called$reference=="T",highcov_r_called$A.SNVcalled,0))
highcov_r_called$TtoC<-as.numeric(ifelse(highcov_r_called$reference=="T",highcov_r_called$C.SNVcalled,0))
highcov_r_called$TtoG<-as.numeric(ifelse(highcov_r_called$reference=="T",highcov_r_called$G.SNVcalled,0))

AtoCfreq<-sum(highcov_r_called$AtoC)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
AtoGfreq<-sum(highcov_r_called$AtoG)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
AtoTfreq<-sum(highcov_r_called$AtoT)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
CtoAfreq<-sum(highcov_r_called$CtoA)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
CtoGfreq<-sum(highcov_r_called$CtoG)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
CtoTfreq<-sum(highcov_r_called$CtoT)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
GtoAfreq<-sum(highcov_r_called$GtoA)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
GtoCfreq<-sum(highcov_r_called$GtoC)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
GtoTfreq<-sum(highcov_r_called$GtoT)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
TtoAfreq<-sum(highcov_r_called$TtoA)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
TtoCfreq<-sum(highcov_r_called$TtoC)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))
TtoGfreq<-sum(highcov_r_called$TtoG)/(sum(highcov_r_called$A.SNVcalled,highcov_r_called$C.SNVcalled,highcov_r_called$G.SNVcalled,highcov_r_called$T.SNVcalled))

## Calculate mutation frequencies for specific subsititutions
highcov_r_called$Acoverage<-as.numeric(ifelse(highcov_r_called$reference=="A",highcov_r_called$coverage,0))
highcov_r_called$Ccoverage<-as.numeric(ifelse(highcov_r_called$reference=="C",highcov_r_called$coverage,0))
highcov_r_called$Gcoverage<-as.numeric(ifelse(highcov_r_called$reference=="G",highcov_r_called$coverage,0))
highcov_r_called$Tcoverage<-as.numeric(ifelse(highcov_r_called$reference=="T",highcov_r_called$coverage,0))

AtoCmutfreq<-sum(highcov_r_called$AtoC)/(sum(highcov_r$Acoverage))
AtoGmutfreq<-sum(highcov_r_called$AtoG)/(sum(highcov_r$Acoverage))
AtoTmutfreq<-sum(highcov_r_called$AtoT)/(sum(highcov_r$Acoverage))
CtoAmutfreq<-sum(highcov_r_called$CtoA)/(sum(highcov_r$Ccoverage))
CtoGmutfreq<-sum(highcov_r_called$CtoG)/(sum(highcov_r$Ccoverage))
CtoTmutfreq<-sum(highcov_r_called$CtoT)/(sum(highcov_r$Ccoverage))
GtoAmutfreq<-sum(highcov_r_called$GtoA)/(sum(highcov_r$Gcoverage))
GtoCmutfreq<-sum(highcov_r_called$GtoC)/(sum(highcov_r$Gcoverage))
GtoTmutfreq<-sum(highcov_r_called$GtoT)/(sum(highcov_r$Gcoverage))
TtoAmutfreq<-sum(highcov_r_called$TtoA)/(sum(highcov_r$Tcoverage))
TtoCmutfreq<-sum(highcov_r_called$TtoC)/(sum(highcov_r$Tcoverage))
TtoGmutfreq<-sum(highcov_r_called$TtoG)/(sum(highcov_r$Tcoverage))

## merge mutation spectrum and write file
spectrum<-c(AtoTfreq,AtoCfreq,AtoGfreq,CtoAfreq,CtoGfreq,CtoTfreq,GtoAfreq,GtoCfreq,GtoTfreq,TtoAfreq,TtoCfreq,TtoGfreq)
spectrum<-t(spectrum)
colnames(spectrum)<-c("AtoU","AtoC","AtoG","CtoA","CtoG","CtoU","GtoA","GtoC","GtoU","UtoA","UtoC","UtoG")
spectrum<-as.data.frame(spectrum)
spectrumname<-paste(name,"_combined_mutational_spectrum_calledSNPsonly.csv",sep="")
csvfile<-paste(outputdir,spectrumname,sep="/")
write.csv(spectrum, file = csvfile, row.names = FALSE)

## merge mutation spectrum frequencies and write file
spectrum_freq<-c(AtoTmutfreq,AtoCmutfreq,AtoGmutfreq,CtoAmutfreq,CtoGmutfreq,CtoTmutfreq,GtoAmutfreq,GtoCmutfreq,GtoTmutfreq,TtoAmutfreq,TtoCmutfreq,TtoGmutfreq)
spectrum_freq<-t(spectrum_freq)
colnames(spectrum_freq)<-c("A>U","A>C","A>G","C>A","C>G","C>U","G>A","G>C","G>U","U>A","U>C","U>G")
spectrum_freq<-as.data.frame(spectrum_freq)
spectrumname<-paste(name,"_combined_mutational_spectrum_frequencies_calledSNPsonly.csv",sep="")
csvfile<-paste(outputdir,spectrumname,sep="/")
write.csv(spectrum_freq, file = csvfile, row.names = FALSE)




