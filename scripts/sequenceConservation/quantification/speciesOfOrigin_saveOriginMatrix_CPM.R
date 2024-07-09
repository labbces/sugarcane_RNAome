library(tximport)
library(reshape2)
library(ggplot2)
library(htmlwidgets)
library(plotly)
library(dplyr)
library(cowplot)

DIR="/Storage/data1/felipe.peres/Sugarcane_ncRNA/11_lncRNASequenceConservation/GenomicReads/PLOT"

setwd(DIR)

countsOrigin<-read.delim(file = "combined_matrix.sf",header=T,row.names=1)
head(countsOrigin)
colnames(countsOrigin)

tx2gene <- read.table(file.path(DIR, "updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv"), header = FALSE, sep = "\t")

# teste merge
countsOrigin <- merge(countsOrigin, tx2gene, by.x = "row.names", by.y = "V3")
colnames(countsOrigin) <- c("Transcript","S._barberi","S._officinarum","S._spontaneum","S._bicolor", "Category", "Gene", "Function", "Transcript Size", "Transcript Function", "Transcript Rfam family", "Transcript GO", "Gene GO")
#colnames(countsOrigin)<-c('S._barberi','S._officinarum','S._spontaneum')
origNumberGenes=dim(countsOrigin)[1]

dim(countsOrigin)
head(countsOrigin)

genomes <- c("S._barberi","S._officinarum","S._spontaneum","S._bicolor")
#table(rowSums(countsOrigin)==0)
table(rowSums(countsOrigin[, genomes])==0)
#countsOrigin<-countsOrigin[!rowSums(countsOrigin)==0,]               #remove transcript without reads mapped in both sps
countsOrigin<-countsOrigin[!rowSums(countsOrigin[, genomes])==0,]     #remove transcript without reads mapped in both sps

countsOrigin$CPM_SBAR<-((countsOrigin[,'S._barberi']+0.1)*10e6)/(sum(countsOrigin[,'S._barberi']))
countsOrigin$CPM_SOFF<-((countsOrigin[,'S._officinarum']+0.1)*10e6)/(sum(countsOrigin[,'S._officinarum']))
countsOrigin$CPM_SSPO<-((countsOrigin[,'S._spontaneum']+0.1)*10e6)/(sum(countsOrigin[,'S._spontaneum']))
#countsOrigin$CPM_SBIC<-((countsOrigin[,'S._bicolor']+0.1)*10e6)/(sum(countsOrigin[,'S._bicolor']))

#ratio spontanenum vs officinarum - AB
countsOrigin$ratioCPM_SSPO_vs_SOFF<-countsOrigin$CPM_SSPO/countsOrigin$CPM_SOFF
countsOrigin$log10ratioCPM_SSPO_vs_SOFF<-round(log10(countsOrigin$CPM_SSPO/countsOrigin$CPM_SOFF),2)

#ratio spontanenum vs barberi - AC
countsOrigin$ratioCPM_SSPO_vs_SBAR<-countsOrigin$CPM_SSPO/countsOrigin$CPM_SBAR
countsOrigin$log10ratioCPM_SSPO_vs_SBAR<-round(log10(countsOrigin$CPM_SSPO/countsOrigin$CPM_SBAR),2)

#ratio spontanenum vs sorghum - AD
#countsOrigin$ratioCPM_SSPO_vs_SBIC<-countsOrigin$CPM_SSPO/countsOrigin$CPM_SBIC
#countsOrigin$log10ratioCPM_SSPO_vs_SBIC<-round(log10(countsOrigin$CPM_SSPO/countsOrigin$CPM_SBIC),2)

# ratio officinarum vs barberi - BC
countsOrigin$ratioCPM_SOFF_vs_SBAR<-countsOrigin$CPM_SOFF/countsOrigin$CPM_SBAR
countsOrigin$log10ratioCPM_SOFF_vs_SBAR<-round(log10(countsOrigin$CPM_SOFF/countsOrigin$CPM_SBAR),2)

#ratio officinarum vs sorghum - BD
#countsOrigin$ratioCPM_SOFF_vs_SBIC<-countsOrigin$CPM_SOFF/countsOrigin$CPM_SBIC
#countsOrigin$log10ratioCPM_SOFF_vs_SBIC<-round(log10(countsOrigin$CPM_SOFF/countsOrigin$CPM_SBIC),2)

# ratio barberi vs sorghum - CD
#countsOrigin$ratioCPM_SBAR_vs_SBIC<-countsOrigin$CPM_SBAR/countsOrigin$CPM_SBIC
#countsOrigin$log10ratioCPM_SBAR_vs_SBIC<-round(log10(countsOrigin$CPM_SBAR/countsOrigin$CPM_SBIC),2)

countsOrigin$Fraction_SSPO<-round(countsOrigin$CPM_SSPO/(countsOrigin$CPM_SOFF+countsOrigin$CPM_SSPO+countsOrigin$CPM_SBAR),3)
countsOrigin$Fraction_SOFF<-round(countsOrigin$CPM_SOFF/(countsOrigin$CPM_SOFF+countsOrigin$CPM_SSPO+countsOrigin$CPM_SBAR),3)
countsOrigin$Fraction_SBAR<-round(countsOrigin$CPM_SBAR/(countsOrigin$CPM_SOFF+countsOrigin$CPM_SSPO+countsOrigin$CPM_SBAR),3)
#countsOrigin$Fraction_SBIC<-round(countsOrigin$CPM_SBIC/(countsOrigin$CPM_SOFF+countsOrigin$CPM_SSPO+countsOrigin$CPM_SBAR+countsOrigin$CPM_SBIC),3)

head(countsOrigin)

countsOrigin[which(countsOrigin$Fraction_SBAR >0.7),]
countsOrigin$Origin<-NA

countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & abs(countsOrigin$log10ratioCPM_SSPO_vs_SOFF) < 1 & abs(countsOrigin$log10ratioCPM_SOFF_vs_SBAR) < 1 & abs(countsOrigin$log10ratioCPM_SSPO_vs_SBAR) < 1),]),'Origin']<-'Common'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$log10ratioCPM_SSPO_vs_SOFF >= 1 & abs(countsOrigin$log10ratioCPM_SOFF_vs_SBAR) < 1 & countsOrigin$log10ratioCPM_SSPO_vs_SBAR >= 1),]),'Origin']<-'SSPO'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$log10ratioCPM_SSPO_vs_SOFF <= -1 & countsOrigin$log10ratioCPM_SOFF_vs_SBAR >= 1 & abs(countsOrigin$log10ratioCPM_SSPO_vs_SBAR) < 1),]),'Origin']<-'SOFF'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & abs(countsOrigin$log10ratioCPM_SSPO_vs_SOFF) < 1 & countsOrigin$log10ratioCPM_SOFF_vs_SBAR <= -1 & countsOrigin$log10ratioCPM_SSPO_vs_SBAR <= -1),]),'Origin']<-'SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._officinarum == 0 & countsOrigin$log10ratioCPM_SSPO_vs_SBAR <= -1),]),'Origin']<-'SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._officinarum == 0 & countsOrigin$log10ratioCPM_SSPO_vs_SBAR >= 1),]),'Origin']<-'SSPO'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._officinarum == 0 & abs(countsOrigin$log10ratioCPM_SSPO_vs_SBAR) < 1),]),'Origin']<-'CommonSSPO_SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._spontaneum == 0 & countsOrigin$log10ratioCPM_SOFF_vs_SBAR <= -1),]),'Origin']<-'SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._spontaneum == 0 & countsOrigin$log10ratioCPM_SOFF_vs_SBAR >= 1),]),'Origin']<-'SOFF'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._spontaneum == 0 & abs(countsOrigin$log10ratioCPM_SOFF_vs_SBAR) < 1),]),'Origin']<-'CommonSOFF_SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._barberi == 0 & countsOrigin$log10ratioCPM_SSPO_vs_SOFF <= -1),]),'Origin']<-'SOFF'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._barberi == 0 & countsOrigin$log10ratioCPM_SSPO_vs_SOFF >= 1),]),'Origin']<-'SSPO'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._barberi == 0 & abs(countsOrigin$log10ratioCPM_SSPO_vs_SOFF) <= 1),]),'Origin']<-'CommonSSPO_SOFF'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin)),]),'Origin']<-'UNK'

table(countsOrigin$Origin, useNA = 'always')*100/sum(table(countsOrigin$Origin, useNA = 'always'))
countsOrigin[which(is.na(countsOrigin$Origin)),]

write.table(countsOrigin, sep = "\t", file = "speciesOfOriginPanTranscriptome.tsv", quote = FALSE, row.names = TRUE)
