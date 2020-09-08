## load required packages
list.of.packages <- c("ggplot2","dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(dplyr)

#### Serial Mouse Samples ####
## Create list of samples ##
mouse<-c("SCP1_1","SCP1_2","SCP1_3","SCP1_4","SCP1_5",
         "SCP3_1","SCP3_2","SCP3_3","SCP3_4","SCP3_5",
         "SCP4_1","SCP4_2","SCP4_3","SCP4_4","SCP4_5",
         "SCP10_1","SCP10_2","SCP10_3","SCP10_4","SCP10_5",
         "SCP10_A","SCP10_B","SCP10_C","SCP10_D","SCP10_E")

## Read and merge ntcount files for both replicates ##
for (i in mouse){
  x1<-paste("/Volumes/HD1/R21_ZIKVConspecificPassage/Mouse_Passage/WGS/",i,sep="")
  y1<-paste(x1,"_rep1_v3.1/",sep="")
  z1<-paste(y1,"reference_aligned/ntcounts_reference_",sep="")
  a1<-paste(z1,i,sep="")
  b1<-paste(a1,"_1.csv",sep="")
  name1<-paste(i,"_1",sep="")
  df1<-read.csv(file=b1,header=TRUE)
  df2<-df1[,c(2,4)]
  df2$ID<-name1
  x2<-paste("/Volumes/HD1/R21_ZIKVConspecificPassage/Mouse_Passage/WGS/",i,sep="")
  y2<-paste(x2,"_rep2_v3.1/",sep="")
  z2<-paste(y2,"reference_aligned/ntcounts_reference_",sep="")
  a2<-paste(z2,i,sep="")
  b2<-paste(a2,"_2.csv",sep="")
  name2<-paste(i,"_2",sep="")
  df3<-read.csv(file=b2,header=TRUE)
  df4<-df3[,c(2,4)]
  df4$ID<-name2
  merge<-rbind(df2,df4)
  assign(i,merge)
}

#### Serial Mosquito Samples ####
## Create list of samples ##
mosquito<-c("MP1_1","MP1_2","MP1_3","MP1_4","MP1_5",
            "MP3_1","MP3_2","MP3_3","MP3_4","MP3_5",
            "MP4_1","MP4_2","MP4_3","MP4_4","MP4_5",
            "MP10_1","MP10_2","MP10_3","MP10_4","MP10_5",
            "MP10_A","MP10_B","MP10_C","MP10_D","MP10_E")

## Read and merge ntcount files for both replicates ##
for (i in mosquito){
  x1<-paste("/Volumes/HD1/R21_ZIKVConspecificPassage/Mosquito_Passage/WGS/",i,sep="")
  y1<-paste(x1,"_rep1_v3.1/",sep="")
  z1<-paste(y1,"reference_aligned/ntcounts_reference_",sep="")
  a1<-paste(z1,i,sep="")
  b1<-paste(a1,"_1.csv",sep="")
  name1<-paste(i,"_1",sep="")
  df1<-read.csv(file=b1,header=TRUE)
  df2<-df1[,c(2,4)]
  df2$ID<-name1
  x2<-paste("/Volumes/HD1/R21_ZIKVConspecificPassage/Mosquito_Passage/WGS/",i,sep="")
  y2<-paste(x2,"_rep2_v3.1/",sep="")
  z2<-paste(y2,"reference_aligned/ntcounts_reference_",sep="")
  a2<-paste(z2,i,sep="")
  b2<-paste(a2,"_2.csv",sep="")
  name2<-paste(i,"_2",sep="")
  df3<-read.csv(file=b2,header=TRUE)
  df4<-df3[,c(2,4)]
  df4$ID<-name2
  merge<-rbind(df2,df4)
  assign(i,merge)
}

#### Alternating Passage Samples ####
## Create list of samples ##
alternate<-c("AP1_1","AP1_2","AP1_3","AP1_4","AP1_5",
             "AP2_1_3B","AP2_1_3S","AP2_2_3B","AP2_3_3B","AP2_3_3S","AP2_3_5B","AP2_3_5S","AP2_3_6B","AP2_3_6S","AP2_4_1B","AP2_4_1S","AP2_4_2B","AP2_4_2S","AP2_4_4B","AP2_4_4S","AP2_4_5B","AP2_4_6B","AP2_4_6S","AP2_5_3B","AP2_5_3S",
             "AP3_1","AP3_2","AP3_3","AP3_4","AP3_5",
             "AP4_1_B","AP4_3_B","AP4_4_B","AP4_5_B",
             "AP5_1","AP5_4","AP5_5",
             "AP6_1_B","AP6_4_B","AP6_5_B",
             "AP7_1","AP7_4","AP7_5",
             "AP8_1_B","AP8_4_B","AP8_5_B",
             "AP9_1","AP9_4","AP9_5",
             "AP10_1_B","AP10_4_B","AP10_5_B")

## Read and merge ntcount files for both replicates ##
for (i in alternate){
  x1<-paste("/Volumes/HD1/R21_ZIKVConspecificPassage/Alternate_Passage/WGS/",i,sep="")
  y1<-paste(x1,"_rep1_v3.1/",sep="")
  z1<-paste(y1,"reference_aligned/ntcounts_reference_",sep="")
  a1<-paste(z1,i,sep="")
  b1<-paste(a1,"_1.csv",sep="")
  name1<-paste(i,"_1",sep="")
  df1<-read.csv(file=b1,header=TRUE)
  df2<-df1[,c(2,4)]
  df2$ID<-name1
  x2<-paste("/Volumes/HD1/R21_ZIKVConspecificPassage/Alternate_Passage/WGS/",i,sep="")
  y2<-paste(x2,"_rep2_v3.1/",sep="")
  z2<-paste(y2,"reference_aligned/ntcounts_reference_",sep="")
  a2<-paste(z2,i,sep="")
  b2<-paste(a2,"_2.csv",sep="")
  name2<-paste(i,"_2",sep="")
  df3<-read.csv(file=b2,header=TRUE)
  df4<-df3[,c(2,4)]
  df4$ID<-name2
  merge<-rbind(df2,df4)
  assign(i,merge)
}

### rbind the coverages
SCP_rbind<-do.call("rbind", list(SCP1_1,SCP1_2,SCP1_3,SCP1_4,SCP1_5,SCP3_1,SCP3_2,SCP3_3,SCP3_4,SCP3_5,
                                 SCP4_1,SCP4_2,SCP4_3,SCP4_4,SCP4_5,SCP10_1,SCP10_2,SCP10_3,SCP10_4,SCP10_5,
                                 SCP10_A,SCP10_B,SCP10_C,SCP10_D,SCP10_E))
MP_rbind<-do.call("rbind", list(MP1_1,MP1_2,MP1_3,MP1_4,MP1_5,MP3_1,MP3_2,MP3_3,MP3_4,MP3_5,
                                MP4_1,MP4_2,MP4_3,MP4_4,MP4_5,MP10_1,MP10_2,MP10_3,MP10_4,MP10_5,
                                MP10_A,MP10_B,MP10_C,MP10_D,MP10_E))
AP_rbind<-do.call("rbind", list(AP1_1,AP1_2,AP1_3,AP1_4,AP1_5,
                                AP2_1_3B,AP2_1_3S,AP2_2_3B,AP2_3_3B,AP2_3_3S,AP2_3_5B,AP2_3_5S,AP2_3_6B,AP2_3_6S,AP2_4_1B,AP2_4_1S,AP2_4_2B,AP2_4_2S,AP2_4_4B,AP2_4_4S,AP2_4_5B,AP2_4_6B,AP2_4_6S,AP2_5_3B,AP2_5_3S,
                                AP3_1,AP3_2,AP3_3,AP3_4,AP3_5,
                                AP4_1_B,AP4_3_B,AP4_4_B,AP4_5_B,
                                AP5_1,AP5_4,AP5_5,
                                AP6_1_B,AP6_4_B,AP6_5_B,
                                AP7_1,AP7_4,AP7_5,
                                AP8_1_B,AP8_4_B,AP8_5_B,
                                AP9_1,AP9_4,AP9_5,
                                AP10_1_B,AP10_4_B,AP10_5_B))


## Generate coverage plots
x<-(rep(c("red","green","blue","yellow","purple","grey","orange","brown","turquoise","pink","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","cyan","darkgrey"),113))

ggplot(data=SCP_rbind, aes(x=pos, y=reads_all, group=ID)) +
  geom_line(show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 300, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("Mouse Passage Samples") +
  labs(x="Genome Position", y="Coverage") +
  scale_x_continuous(limits = c(1, 10675), expand = c(0,0)) +
  scale_y_continuous(limits = c(1, 12000), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("PATH/TO/MouseSamples_downsampled_coverages.png",plot=last_plot(),device = png(), height = 3, width = 8, dpi = 300)

ggplot(data=MP_rbind, aes(x=pos, y=reads_all, group=ID)) +
  geom_line(show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 300, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("Mosquito Passage Samples") +
  labs(x="Genome Position", y="Coverage") +
  scale_x_continuous(limits = c(1, 10675), expand = c(0,0)) +
  scale_y_continuous(limits = c(1, 12000), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("PATH/TO/MosquitoSamples_downsampled_coverages.png",plot=last_plot(),device = png(), height = 3, width = 8, dpi = 300)

ggplot(data=AP_rbind, aes(x=pos, y=reads_all, group=ID)) +
  geom_line(show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 300, colour = "lightgrey", linetype = "dashed") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  ggtitle("Alternate Passage Samples") +
  labs(x="Genome Position", y="Coverage") +
  scale_x_continuous(limits = c(1, 10675), expand = c(0,0)) +
  scale_y_continuous(limits = c(1, 12000), expand = c(0,0)) +
  scale_color_manual(values=x)

ggsave("PATH/TO/AlternateSamples_downsampled_coverages.png",plot=last_plot(),device = png(), height = 3, width = 8, dpi = 300)


#Lessthan300<-subset(SCP_rbind,SCP_rbind$reads_all<300 & SCP_rbind$pos<10413)
#Lessthan300 %>% count(ID)

#Lessthan300<-subset(MP_rbind,MP_rbind$reads_all<300 & MP_rbind$pos<10413)
#Lessthan300 %>% count(ID)

#Lessthan300<-subset(AP_rbind,AP_rbind$reads_all<300 & AP_rbind$pos<10413)
#Lessthan300 %>% count(ID)





