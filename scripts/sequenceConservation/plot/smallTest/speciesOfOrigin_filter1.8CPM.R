library(tximport)
library(reshape2)
library(ggplot2)
library(htmlwidgets)
library(plotly)
library(dplyr)

#rm(list=ls())

#DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/sequenceConservation/plot/smallTest/"
DIR="/Storage/data1/felipe.peres/Sugarcane_ncRNA/11_lncRNASequenceConservation/GenomicReads/PLOT"

setwd(DIR)

countsOrigin<-read.delim(file = "combined_matrix_100k.sf",header=T,row.names=1)
countsOrigin<-read.delim(file = "combined_matrix.sf",header=T,row.names=1)
head(countsOrigin)
colnames(countsOrigin)

#colnames(countsOrigin)<-c('S._barberi','S._officinarum','S._spontaneum')
origNumberGenes=dim(countsOrigin)[1]

dim(countsOrigin)
head(countsOrigin)

table(rowSums(countsOrigin)==0)
countsOrigin<-countsOrigin[!rowSums(countsOrigin)==0,] #remove transcript without reads mapped in both sps

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
summary(countsOrigin$CPM_SBAR)

quantile(countsOrigin$Fraction_SBAR, c(.1,.2,.5,.8,.9,.99))
quantile(countsOrigin$CPM_SBAR, c(.1,.2,.5,.8,.9,.99))
quantile(countsOrigin$CPM_SSPO, c(.1,.2,.5,.8,.9,.99))
quantile(countsOrigin$CPM_SOFF, c(.1,.2,.5,.8,.9,.99))

colnames(countsOrigin)

# Log base 10 transformation
countsOrigin$log10_CPM_SOFF <- log10(countsOrigin$CPM_SOFF)
countsOrigin$log10_CPM_SSPO <- log10(countsOrigin$CPM_SSPO)
countsOrigin$log10_CPM_SBAR <- log10(countsOrigin$CPM_SBAR)

fig <- plot_ly(countsOrigin[which((countsOrigin$CPM_SOFF > 1.8 |
                                       countsOrigin$CPM_SSPO > 1.8 |
                                       countsOrigin$CPM_SBAR > 1.8)),], x = ~log10_CPM_SSPO,
               y = ~log10_CPM_SOFF,
               z = ~log10_CPM_SBAR,
               color = ~Origin)

fig <- fig %>% add_markers()

fig <- fig %>% layout(scene = list(xaxis = list(type = "log", title = 'S. barberi (log10 CPM)'),
                                   yaxis = list(type = "log", title = 'S. officinarum (log10 CPM)'),
                                   zaxis = list(type = "log", title = 'S. spontaneum (log10 CPM)')))

#saveWidget(fig, "speciesOfOriginPanRNAome_log10.html")

############# melhorando plot 3d
# Criar o gráfico 3D
fig <- plot_ly(filtered_data, 
               x = ~log10_CPM_SSPO, 
               y = ~log10_CPM_SOFF, 
               z = ~log10_CPM_SBAR, 
               color = ~Origin, 
               colors = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'),
               marker = list(size = 4, opacity = 0.7))

fig <- fig %>% add_markers()

# Melhorar o layout
fig <- fig %>% layout(
  title = "Conservação de Transcritos em Três Genomas Diferentes",
  scene = list(
    xaxis = list(type = "log", title = 'S. spontaneum (log10 CPM)'),
    yaxis = list(type = "log", title = 'S. officinarum (log10 CPM)'),
    zaxis = list(type = "log", title = 'S. barberi (log10 CPM)'),
    camera = list(
      eye = list(x = 1.5, y = 1.5, z = 1.5)
    )
  ),
  legend = list(title = list(text = 'Origin'), orientation = 'h', x = 0.5, y = -0.1, xanchor = 'center')
)

#fig
################### /melhorando plot 3d

ggplot(as.data.frame(countsOrigin),aes(x=Fraction_SBAR)) +
  theme_bw()+
  geom_histogram(bins=1000)+
  scale_x_log10()
  
  
#ggplot(as.data.frame(countsOrigin),aes(x=ratioCPM)) +
#  theme_bw()+
#  geom_histogram(bins=1000)+
#  scale_x_log10()+
#  geom_vline(xintercept = 1,color='red')

#ggplot(as.data.frame(countsOrigin),aes(x=log10ratioCPM)) +
#  theme_bw()+
#  geom_histogram(bins=1000)+
#  geom_vline(xintercept = 0,color='red')+
#  geom_vline(xintercept = -1,color='blue')+
#  geom_vline(xintercept = 1,color='blue')+
#  xlab('ratio of CPM, log10 scale')+
#  ylab('Number of transcripts')

countsOrigin$Origin

# interesting features
origin_SOFF_SSPO <- c("SOFF", "SSPO", "CommonSSPO_SOFF", "Common", "UNK")
origin_SOFF_SBAR <- c("SOFF", "SBAR", "CommonSOFF_SBAR", "Common", "UNK")
origin_SSPO_SBAR <- c("SSPO", "SBAR", "CommonSSPO_SBAR", "Common", "UNK")

#origin <- c("Common")

# filter for SOFF and SSPO
filtered_origin_SOFF_SSPO <- countsOrigin %>%
  filter(Origin %in% origin_SOFF_SSPO)
filtered_origin_SOFF_SSPO$Origin <- as.factor(filtered_origin_SOFF_SSPO$Origin)

# filter for SOFF and SBAR
filtered_origin_SOFF_SBAR <- countsOrigin %>%
  filter(Origin %in% origin_SOFF_SBAR)
filtered_origin_SOFF_SBAR$Origin <- as.factor(filtered_origin_SOFF_SBAR$Origin)

# filter for SSPO and SBAR
filtered_origin_SSPO_SBAR <- countsOrigin %>%
  filter(Origin %in% origin_SSPO_SBAR)
filtered_origin_SSPO_SBAR$Origin <- as.factor(filtered_origin_SSPO_SBAR$Origin)

ggplot(as.data.frame(filtered_origin_SOFF_SSPO),aes(x=CPM_SOFF, y=CPM_SSPO)) +
  theme_bw()+
  geom_jitter()+
  xlab('Log10 CPM S. spontaneum')+
  ylab('Log10 CPM S. officinarum')+
  scale_x_log10()+
  scale_y_log10()

#dim(countsOrigin[which(countsOrigin$log10ratioCPM>-1 & countsOrigin$log10ratioCPM<1),]) #possible recombinants, or common to both species
#dim(countsOrigin[which(countsOrigin$log10ratioCPM>=1),]) # Most likely from Spontaneum
#dim(countsOrigin[which(countsOrigin$log10ratioCPM<=-1),]) # Most likely from Officinarum

dim(countsOrigin)[1]
#Genes with mapped reads
dim(countsOrigin)[1]*100/origNumberGenes

#Proportion Genes with mapped reads assigned to common
#dim(countsOrigin[which(countsOrigin$log10ratioCPM>-1 & countsOrigin$log10ratioCPM<1),])[1]*100/dim(countsOrigin)[1]
#Proportion Genes with mapped reads assigned to spontaneum
#dim(countsOrigin[which(countsOrigin$log10ratioCPM>=1),])[1]*100/dim(countsOrigin)[1] 
#Proportion Genes with mapped reads assigned to offcinarum
#dim(countsOrigin[which(countsOrigin$log10ratioCPM<=-1),])[1]*100/dim(countsOrigin)[1] 

#countsOrigin$Origin<-NA
#countsOrigin[which(countsOrigin$log10ratioCPM>-1 & countsOrigin$log10ratioCPM<1),'Origin']<-'Common'
#countsOrigin[which(countsOrigin$log10ratioCPM>=1),'Origin']<-'S._spontaneum'
#countsOrigin[which(countsOrigin$log10ratioCPM<=-1),'Origin']<-'S._officinarum'
#head(countsOrigin)

ggplot(as.data.frame(filtered_origin_SOFF_SSPO),aes(x=CPM_SOFF, y=CPM_SSPO)) +
  theme(text=element_text(size=20))+
  theme_bw()+
  geom_jitter(aes(colour=Origin),alpha=0.1)+
  xlab('Log10 CPM S. officinarum')+
  ylab('Log10 CPM S. spontaneum')+
  scale_x_log10()+
  scale_y_log10()
  #scale_colour_manual(values = "blue")

# plot SOFF and SSPO
plot <- ggplot(as.data.frame(filtered_origin_SOFF_SSPO), aes(x=CPM_SOFF, y=CPM_SSPO)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "right") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. officinarum') +
  ylab('Log10 CPM S. spontaneum') +
  ggtitle("Common transcripts between SOFF and SSPO",
          subtitle = "Visualization of conserved/common transcripts") +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  annotate("text", x = 1e+04, y = 1e+00, label = "2D density lines show\nregions of high concentration", 
           size = 5, hjust = 0, color = "blue") 

# save as PNG
ggsave("common_SOFF_SSPO.png", plot = plot, width = 13, height = 11, dpi = 500, bg="white")
# save as SVG
#ggsave("common_SOFF_SSPO.svg", plot = plot, width = 13, height = 11)

# plot SOFF and SBAR
plot <- ggplot(as.data.frame(filtered_origin_SOFF_SBAR), aes(x=CPM_SOFF, y=CPM_SBAR)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "right") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. officinarum') +
  ylab('Log10 CPM S. barberi') +
  ggtitle("Common transcripts between SOFF and SBAR",
          subtitle = "Visualization of conserved/common transcripts") +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  annotate("text", x = 1e+04, y = 1e+00, label = "2D density lines show\nregions of high concentration", 
           size = 5, hjust = 0, color = "blue") 

# save as PNG
ggsave("common_SOFF_SBAR.png", plot = plot, width = 13, height = 11, dpi = 500, bg="white")
# save as SVG
#ggsave("common_SOFF_SBAR.svg", plot = plot, width = 13, height = 11)

# plot SSPO and SBAR
plot <- ggplot(as.data.frame(filtered_origin_SSPO_SBAR), aes(x=CPM_SSPO, y=CPM_SBAR)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "right") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. spontaneum') +
  ylab('Log10 CPM S. barberi') +
  ggtitle("Common transcripts between SSPO and SBAR",
          subtitle = "Visualization of conserved/common transcripts") +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  annotate("text", x = 1e+04, y = 1e+00, label = "2D density lines show\nregions of high concentration", 
           size = 5, hjust = 0, color = "blue") 

# save as PNG
ggsave("common_SSPO_SBAR.png", plot = plot, width = 13, height = 11, dpi = 500, bg="white")
# save as SVG
#ggsave("common_SSPO_SBAR.svg", plot = plot, width = 13, height = 11)

#write.table(countsOrigin, sep = "\t", file = "speciesOfOriginPanTranscriptome.csv")

#tx2gene<-read.delim("txp.group.tsv",sep=',',header=T)
#head(tx2gene)
#length(unique(tx2gene$GeneID))
#length(unique(tx2gene$TranscriptID))
