library(tximport)

DIR="/Storage/data1/felipe.peres/Sugarcane_ncRNA/11_lncRNASequenceConservation/GenomicReads/PLOT"

setwd(DIR)

countsOrigin<-read.delim(file = "combined_matrix.sf",header=T,row.names=1)
tx2gene <- read.table(file.path(DIR, "updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv"), header = FALSE, sep = "\t")

# merge countsOrigin with annotation table
countsOrigin <- merge(countsOrigin, tx2gene, by.x = "row.names", by.y = "V3")
colnames(countsOrigin) <- c("Transcript","S._barberi","S._officinarum","S._spontaneum","S._bicolor", "EffectiveLength", "Category", "Gene", "Function", "Transcript Size", "Transcript Function", "Transcript Rfam family", "Transcript GO", "Gene GO")
origNumberGenes=dim(countsOrigin)[1]

dim(countsOrigin)
head(countsOrigin)

genomes <- c("S._barberi","S._officinarum","S._spontaneum","S._bicolor")
table(rowSums(countsOrigin[, genomes])==0)
countsOrigin<-countsOrigin[!rowSums(countsOrigin[, genomes])==0,]     # remove transcript without reads mapped in both sps

# function from this blog post -> https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
countsToFPKM <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

counts_SBAR <- countsOrigin$S._barberi
counts_SOFF <- countsOrigin$S._officinarum
counts_SSPO <- countsOrigin$S._spontaneum
counts_SBIC <- countsOrigin$S._bicolor
effLen <- countsOrigin$EffectiveLength

# apply countsToFPKM
countsOrigin$FPKM_SBAR <- countsToFPKM(counts_SBAR, effLen)
countsOrigin$FPKM_SOFF <- countsToFPKM(counts_SOFF, effLen)
countsOrigin$FPKM_SSPO <- countsToFPKM(counts_SSPO, effLen)
countsOrigin$FPKM_SBIC <- countsToFPKM(counts_SBIC, effLen)

# ratio spontanenum vs officinarum
countsOrigin$ratioFPKM_SSPO_vs_SOFF <- countsOrigin$FPKM_SSPO/countsOrigin$FPKM_SOFF
countsOrigin$log10ratioFPKM_SSPO_vs_SOFF <- round(log10(countsOrigin$FPKM_SSPO/countsOrigin$FPKM_SOFF),2)

# ratio spontanenum vs barberi
countsOrigin$ratioFPKM_SSPO_vs_SBAR<-countsOrigin$FPKM_SSPO/countsOrigin$FPKM_SBAR
countsOrigin$log10ratioFPKM_SSPO_vs_SBAR<-round(log10(countsOrigin$FPKM_SSPO/countsOrigin$FPKM_SBAR),2)

# ratio officinarum vs barberi
countsOrigin$ratioFPKM_SOFF_vs_SBAR<-countsOrigin$FPKM_SOFF/countsOrigin$FPKM_SBAR
countsOrigin$log10ratioFPKM_SOFF_vs_SBAR<-round(log10(countsOrigin$FPKM_SOFF/countsOrigin$FPKM_SBAR),2)

countsOrigin$Fraction_SSPO<-round(countsOrigin$FPKM_SSPO/(countsOrigin$FPKM_SOFF+countsOrigin$FPKM_SSPO+countsOrigin$FPKM_SBAR),3)
countsOrigin$Fraction_SOFF<-round(countsOrigin$FPKM_SOFF/(countsOrigin$FPKM_SOFF+countsOrigin$FPKM_SSPO+countsOrigin$FPKM_SBAR),3)
countsOrigin$Fraction_SBAR<-round(countsOrigin$FPKM_SBAR/(countsOrigin$FPKM_SOFF+countsOrigin$FPKM_SSPO+countsOrigin$FPKM_SBAR),3)

head(countsOrigin)

countsOrigin[which(countsOrigin$Fraction_SBAR >0.7),]
countsOrigin$Origin<-NA

countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & abs(countsOrigin$log10ratioFPKM_SSPO_vs_SOFF) < 1 & abs(countsOrigin$log10ratioFPKM_SOFF_vs_SBAR) < 1 & abs(countsOrigin$log10ratioFPKM_SSPO_vs_SBAR) < 1),]),'Origin']<-'Common'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$log10ratioFPKM_SSPO_vs_SOFF >= 1 & abs(countsOrigin$log10ratioFPKM_SOFF_vs_SBAR) < 1 & countsOrigin$log10ratioFPKM_SSPO_vs_SBAR >= 1),]),'Origin']<-'SSPO'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$log10ratioFPKM_SSPO_vs_SOFF <= -1 & countsOrigin$log10ratioFPKM_SOFF_vs_SBAR >= 1 & abs(countsOrigin$log10ratioFPKM_SSPO_vs_SBAR) < 1),]),'Origin']<-'SOFF'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & abs(countsOrigin$log10ratioFPKM_SSPO_vs_SOFF) < 1 & countsOrigin$log10ratioFPKM_SOFF_vs_SBAR <= -1 & countsOrigin$log10ratioFPKM_SSPO_vs_SBAR <= -1),]),'Origin']<-'SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._officinarum == 0 & countsOrigin$log10ratioFPKM_SSPO_vs_SBAR <= -1),]),'Origin']<-'SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._officinarum == 0 & countsOrigin$log10ratioFPKM_SSPO_vs_SBAR >= 1),]),'Origin']<-'SSPO'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._officinarum == 0 & abs(countsOrigin$log10ratioFPKM_SSPO_vs_SBAR) < 1),]),'Origin']<-'CommonSSPO_SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._spontaneum == 0 & countsOrigin$log10ratioFPKM_SOFF_vs_SBAR <= -1),]),'Origin']<-'SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._spontaneum == 0 & countsOrigin$log10ratioFPKM_SOFF_vs_SBAR >= 1),]),'Origin']<-'SOFF'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._spontaneum == 0 & abs(countsOrigin$log10ratioFPKM_SOFF_vs_SBAR) < 1),]),'Origin']<-'CommonSOFF_SBAR'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._barberi == 0 & countsOrigin$log10ratioFPKM_SSPO_vs_SOFF <= -1),]),'Origin']<-'SOFF'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._barberi == 0 & countsOrigin$log10ratioFPKM_SSPO_vs_SOFF >= 1),]),'Origin']<-'SSPO'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin) & countsOrigin$S._barberi == 0 & abs(countsOrigin$log10ratioFPKM_SSPO_vs_SOFF) <= 1),]),'Origin']<-'CommonSSPO_SOFF'
countsOrigin[rownames(countsOrigin[which(is.na(countsOrigin$Origin)),]),'Origin']<-'UNK'

table(countsOrigin$Origin, useNA = 'always')*100/sum(table(countsOrigin$Origin, useNA = 'always'))
countsOrigin[which(is.na(countsOrigin$Origin)),]

# naming first column as "Name"
countsOrigin$Name <- rownames(countsOrigin)
countsOrigin <- countsOrigin[, c(ncol(countsOrigin), 1:(ncol(countsOrigin)-1))]

write.table(countsOrigin, sep = "\t", file = "speciesOfOriginPanTranscriptome_FPKM.tsv", quote = FALSE, row.names = FALSE)