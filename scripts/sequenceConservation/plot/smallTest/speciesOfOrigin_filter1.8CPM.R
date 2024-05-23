library(tximport)
library(reshape2)
library(ggplot2)
library(htmlwidgets)
library(plotly)
library(dplyr)
library(cowplot)

#rm(list=ls())

#DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/sequenceConservation/plot/smallTest/"
DIR="/Storage/data1/felipe.peres/Sugarcane_ncRNA/11_lncRNASequenceConservation/GenomicReads/PLOT"

setwd(DIR)

#countsOrigin<-read.delim(file = "combined_matrix_100k.sf",header=T,row.names=1)
countsOrigin<-read.delim(file = "combined_matrix.sf",header=T,row.names=1)
head(countsOrigin)
colnames(countsOrigin)

tx2gene <- read.table(file.path(DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv"), header = FALSE, sep = "\t")

# teste merge
countsOrigin <- merge(countsOrigin, tx2gene, by.x = "row.names", by.y = "V3")
colnames(countsOrigin) <- c("Transcript","S._barberi","S._officinarum","S._spontaneum","S._bicolor", "Category", "Gene", "Function")
#colnames(countsOrigin)<-c('S._barberi','S._officinarum','S._spontaneum')
origNumberGenes=dim(countsOrigin)[1]

dim(countsOrigin)
head(countsOrigin)

genomes <- c("S._barberi","S._officinarum","S._spontaneum","S._bicolor")
#table(rowSums(countsOrigin)==0)
table(rowSums(countsOrigin[, genomes])==0)
#countsOrigin<-countsOrigin[!rowSums(countsOrigin)==0,] #remove transcript without reads mapped in both sps
countsOrigin<-countsOrigin[!rowSums(countsOrigin[, genomes])==0,] #remove transcript without reads mapped in both sps

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

#fig <- plot_ly(countsOrigin[which((countsOrigin$CPM_SOFF > 1.8 |
#                                       countsOrigin$CPM_SSPO > 1.8 |
#                                       countsOrigin$CPM_SBAR > 1.8)),], x = ~log10_CPM_SSPO,
#               y = ~log10_CPM_SOFF,
#               z = ~log10_CPM_SBAR,
#               color = ~Origin)
#
#fig <- fig %>% add_markers()

#fig <- fig %>% layout(scene = list(xaxis = list(type = "log", title = 'S. barberi (log10 CPM)'),
#                                   yaxis = list(type = "log", title = 'S. officinarum (log10 CPM)'),
#                                   zaxis = list(type = "log", title = 'S. spontaneum (log10 CPM)')))

#saveWidget(fig, "speciesOfOriginPanRNAome_log10.html")

############# melhorando plot 3d
# Criar o gráfico 3D
#fig <- plot_ly(filtered_data, 
#               x = ~log10_CPM_SSPO, 
#               y = ~log10_CPM_SOFF, 
#               z = ~log10_CPM_SBAR, 
#               color = ~Origin, 
#               colors = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'),
#               marker = list(size = 4, opacity = 0.7))

#fig <- fig %>% add_markers()

# Melhorar o layout
#fig <- fig %>% layout(
#  title = "Conservação de Transcritos em Três Genomas Diferentes",
#  scene = list(
#    xaxis = list(type = "log", title = 'S. spontaneum (log10 CPM)'),
#    yaxis = list(type = "log", title = 'S. officinarum (log10 CPM)'),
#    zaxis = list(type = "log", title = 'S. barberi (log10 CPM)'),
#    camera = list(
#      eye = list(x = 1.5, y = 1.5, z = 1.5)
#    )
#  ),
#  legend = list(title = list(text = 'Origin'), orientation = 'h', x = 0.5, y = -0.1, xanchor = 'center')
#)

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



#################################### SOFF and SSPO - all transcripts

filtered_origin_SOFF_SSPO_coding <- subset(filtered_origin_SOFF_SSPO, Function == "protein-coding")
filtered_origin_SOFF_SSPO_noncoding <- subset(filtered_origin_SOFF_SSPO, Function == "non-coding")

# Contagem de transcritos por função
function_counts <- filtered_origin_SOFF_SSPO %>%
  group_by(Origin) %>%
  summarise(Count = n())

# Definir a paleta de cores comum
color_palette <- scale_fill_brewer(palette = "Set1")

# Gráfico de barras estilizado para as contagens de funções
bar_plot <- ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
  geom_bar(stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=Count), vjust=-0.5, size=5) +  # Adiciona as contagens em cada barra
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey")
  ) +
  ylab('Transcripts') +
  xlab(NULL) +  # Remover o rótulo do eixo X
  color_palette +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Gráfico de dispersão
scatter_plot <- ggplot(as.data.frame(filtered_origin_SOFF_SSPO), aes(x=CPM_SOFF, y=CPM_SSPO)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. officinarum') +
  ylab('Log10 CPM S. spontaneum') +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Extrair a legenda do gráfico de dispersão
legend <- get_legend(
  scatter_plot + theme(legend.position = "right")
)

# Remover a legenda do gráfico de dispersão principal
scatter_plot <- scatter_plot + theme(legend.position = "none")

# Combinar os gráficos lado a lado com a legenda
combined_plots <- plot_grid(
  scatter_plot, bar_plot, legend,
  ncol = 3, rel_widths = c(3, 2, 0.8)  # Ajuste das larguras relativas para evitar sobreposição
#  align = "h"
)

# Adicionar o título e subtítulo principal
title <- ggdraw() + 
  draw_label("Common transcripts between SOFF and SSPO", fontface = 'bold', size = 22, hjust = 0.5) +
  draw_label("Visualization of conserved/common transcripts", fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
  theme(plot.margin = margin(t = 10, b = 20))

# Combinar o título e os gráficos
final_plot <- plot_grid(
  title, combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Salvar a combinação como PNG
ggsave("common_SOFF_SSPO.png", plot = final_plot, width = 18, height = 11, dpi = 500, bg="white")






#################################### SOFF and SSPO - coding

# Contagem de transcritos por função
function_counts <- filtered_origin_SOFF_SSPO_coding %>%
  group_by(Origin) %>%
  summarise(Count = n())

# Definir a paleta de cores comum
color_palette <- scale_fill_brewer(palette = "Set1")

# Gráfico de barras estilizado para as contagens de funções
bar_plot <- ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
  geom_bar(stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=Count), vjust=-0.5, size=5) +  # Adiciona as contagens em cada barra
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey")
  ) +
  ylab('Transcripts') +
  xlab(NULL) +  # Remover o rótulo do eixo X
  color_palette +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Gráfico de dispersão
scatter_plot <- ggplot(as.data.frame(filtered_origin_SOFF_SSPO_coding), aes(x=CPM_SOFF, y=CPM_SSPO)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. officinarum') +
  ylab('Log10 CPM S. spontaneum') +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Extrair a legenda do gráfico de dispersão
legend <- get_legend(
  scatter_plot + theme(legend.position = "right")
)

# Remover a legenda do gráfico de dispersão principal
scatter_plot <- scatter_plot + theme(legend.position = "none")

# Combinar os gráficos lado a lado com a legenda
combined_plots <- plot_grid(
  scatter_plot, bar_plot, legend,
  ncol = 3, rel_widths = c(3, 2, 0.8)  # Ajuste das larguras relativas para evitar sobreposição
  #  align = "h"
)

# Adicionar o título e subtítulo principal
title <- ggdraw() + 
  draw_label("Common coding transcripts between SOFF and SSPO", fontface = 'bold', size = 22, hjust = 0.5) +
  draw_label("Visualization of conserved/common transcripts", fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
  theme(plot.margin = margin(t = 10, b = 20))

# Combinar o título e os gráficos
final_plot <- plot_grid(
  title, combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Salvar a combinação como PNG
ggsave("common_SOFF_SSPO_coding.png", plot = final_plot, width = 18, height = 11, dpi = 500, bg="white")







#################################### SOFF and SSPO - non-coding

# Contagem de transcritos por função
function_counts <- filtered_origin_SOFF_SSPO_noncoding %>%
  group_by(Origin) %>%
  summarise(Count = n())

# Definir a paleta de cores comum
color_palette <- scale_fill_brewer(palette = "Set1")

# Gráfico de barras estilizado para as contagens de funções
bar_plot <- ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
  geom_bar(stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=Count), vjust=-0.5, size=5) +  # Adiciona as contagens em cada barra
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey")
  ) +
  ylab('Transcripts') +
  xlab(NULL) +  # Remover o rótulo do eixo X
  color_palette +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Gráfico de dispersão
scatter_plot <- ggplot(as.data.frame(filtered_origin_SOFF_SSPO_noncoding), aes(x=CPM_SOFF, y=CPM_SSPO)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. officinarum') +
  ylab('Log10 CPM S. spontaneum') +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Extrair a legenda do gráfico de dispersão
legend <- get_legend(
  scatter_plot + theme(legend.position = "right")
)

# Remover a legenda do gráfico de dispersão principal
scatter_plot <- scatter_plot + theme(legend.position = "none")

# Combinar os gráficos lado a lado com a legenda
combined_plots <- plot_grid(
  scatter_plot, bar_plot, legend,
  ncol = 3, rel_widths = c(3, 2, 0.8)  # Ajuste das larguras relativas para evitar sobreposição
  #  align = "h"
)

# Adicionar o título e subtítulo principal
title <- ggdraw() + 
  draw_label("Common non-coding transcripts between SOFF and SSPO", fontface = 'bold', size = 22, hjust = 0.5) +
  draw_label("Visualization of conserved/common transcripts", fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
  theme(plot.margin = margin(t = 10, b = 20))

# Combinar o título e os gráficos
final_plot <- plot_grid(
  title, combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Salvar a combinação como PNG
ggsave("common_SOFF_SSPO_noncoding.png", plot = final_plot, width = 18, height = 11, dpi = 500, bg="white")










############################### plot SOFF and SBAR - all transcripts


filtered_origin_SOFF_SBAR_coding <- subset(filtered_origin_SOFF_SBAR, Function == "protein-coding")
filtered_origin_SOFF_SBAR_noncoding <- subset(filtered_origin_SOFF_SBAR, Function == "non-coding")

# Contagem de transcritos por função
function_counts <- filtered_origin_SOFF_SBAR %>%
  group_by(Origin) %>%
  summarise(Count = n())

# Definir a paleta de cores comum
color_palette <- scale_fill_brewer(palette = "Set1")

# Gráfico de barras estilizado para as contagens de funções
bar_plot <- ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
  geom_bar(stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=Count), vjust=-0.5, size=5) +  # Adiciona as contagens em cada barra
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey")
  ) +
  ylab('Transcripts') +
  xlab(NULL) +  # Remover o rótulo do eixo X
  color_palette +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Gráfico de dispersão
scatter_plot <- ggplot(as.data.frame(filtered_origin_SOFF_SBAR), aes(x=CPM_SOFF, y=CPM_SBAR)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. officinarum') +
  ylab('Log10 CPM S. barberi') +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Extrair a legenda do gráfico de dispersão
legend <- get_legend(
  scatter_plot + theme(legend.position = "right")
)

# Remover a legenda do gráfico de dispersão principal
scatter_plot <- scatter_plot + theme(legend.position = "none")

# Combinar os gráficos lado a lado com a legenda
combined_plots <- plot_grid(
  scatter_plot, bar_plot, legend,
  ncol = 3, rel_widths = c(3, 2, 0.8)  # Ajuste das larguras relativas para evitar sobreposição
  #  align = "h"
)

# Adicionar o título e subtítulo principal
title <- ggdraw() + 
  draw_label("Common transcripts between SOFF and SBAR", fontface = 'bold', size = 22, hjust = 0.5) +
  draw_label("Visualization of conserved/common transcripts", fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
  theme(plot.margin = margin(t = 10, b = 20))

# Combinar o título e os gráficos
final_plot <- plot_grid(
  title, combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Salvar a combinação como PNG
ggsave("common_SOFF_SBAR.png", plot = final_plot, width = 18, height = 11, dpi = 500, bg="white")






#################################### SOFF and SBAR - coding

# Contagem de transcritos por função
function_counts <- filtered_origin_SOFF_SBAR_coding %>%
  group_by(Origin) %>%
  summarise(Count = n())

# Definir a paleta de cores comum
color_palette <- scale_fill_brewer(palette = "Set1")

# Gráfico de barras estilizado para as contagens de funções
bar_plot <- ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
  geom_bar(stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=Count), vjust=-0.5, size=5) +  # Adiciona as contagens em cada barra
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey")
  ) +
  ylab('Transcripts') +
  xlab(NULL) +  # Remover o rótulo do eixo X
  color_palette +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Gráfico de dispersão
scatter_plot <- ggplot(as.data.frame(filtered_origin_SOFF_SBAR_coding), aes(x=CPM_SOFF, y=CPM_SBAR)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. officinarum') +
  ylab('Log10 CPM S. barberi') +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Extrair a legenda do gráfico de dispersão
legend <- get_legend(
  scatter_plot + theme(legend.position = "right")
)

# Remover a legenda do gráfico de dispersão principal
scatter_plot <- scatter_plot + theme(legend.position = "none")

# Combinar os gráficos lado a lado com a legenda
combined_plots <- plot_grid(
  scatter_plot, bar_plot, legend,
  ncol = 3, rel_widths = c(3, 2, 0.8)  # Ajuste das larguras relativas para evitar sobreposição
  #  align = "h"
)

# Adicionar o título e subtítulo principal
title <- ggdraw() + 
  draw_label("Common coding transcripts between SOFF and SBAR", fontface = 'bold', size = 22, hjust = 0.5) +
  draw_label("Visualization of conserved/common transcripts", fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
  theme(plot.margin = margin(t = 10, b = 20))

# Combinar o título e os gráficos
final_plot <- plot_grid(
  title, combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Salvar a combinação como PNG
ggsave("common_SOFF_SBAR_coding.png", plot = final_plot, width = 18, height = 11, dpi = 500, bg="white")







#################################### SOFF and SBAR - non-coding

# Contagem de transcritos por função
function_counts <- filtered_origin_SOFF_SBAR_noncoding %>%
  group_by(Origin) %>%
  summarise(Count = n())

# Definir a paleta de cores comum
color_palette <- scale_fill_brewer(palette = "Set1")

# Gráfico de barras estilizado para as contagens de funções
bar_plot <- ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
  geom_bar(stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=Count), vjust=-0.5, size=5) +  # Adiciona as contagens em cada barra
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey")
  ) +
  ylab('Transcripts') +
  xlab(NULL) +  # Remover o rótulo do eixo X
  color_palette +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Gráfico de dispersão
scatter_plot <- ggplot(as.data.frame(filtered_origin_SOFF_SBAR_noncoding), aes(x=CPM_SOFF, y=CPM_SBAR)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. officinarum') +
  ylab('Log10 CPM S. barberi') +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Extrair a legenda do gráfico de dispersão
legend <- get_legend(
  scatter_plot + theme(legend.position = "right")
)

# Remover a legenda do gráfico de dispersão principal
scatter_plot <- scatter_plot + theme(legend.position = "none")

# Combinar os gráficos lado a lado com a legenda
combined_plots <- plot_grid(
  scatter_plot, bar_plot, legend,
  ncol = 3, rel_widths = c(3, 2, 0.8)  # Ajuste das larguras relativas para evitar sobreposição
  #  align = "h"
)

# Adicionar o título e subtítulo principal
title <- ggdraw() + 
  draw_label("Common non-coding transcripts between SOFF and SBAR", fontface = 'bold', size = 22, hjust = 0.5) +
  draw_label("Visualization of conserved/common transcripts", fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
  theme(plot.margin = margin(t = 10, b = 20))

# Combinar o título e os gráficos
final_plot <- plot_grid(
  title, combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Salvar a combinação como PNG
ggsave("common_SOFF_SBAR_noncoding.png", plot = final_plot, width = 18, height = 11, dpi = 500, bg="white")








############################### plot SSPO and SBAR - all transcripts

filtered_origin_SSPO_SBAR_coding <- subset(filtered_origin_SSPO_SBAR, Function == "protein-coding")
filtered_origin_SSPO_SBAR_noncoding <- subset(filtered_origin_SSPO_SBAR, Function == "non-coding")

# Contagem de transcritos por função
function_counts <- filtered_origin_SSPO_SBAR %>%
  group_by(Origin) %>%
  summarise(Count = n())

# Definir a paleta de cores comum
color_palette <- scale_fill_brewer(palette = "Set1")

# Gráfico de barras estilizado para as contagens de funções
bar_plot <- ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
  geom_bar(stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=Count), vjust=-0.5, size=5) +  # Adiciona as contagens em cada barra
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey")
  ) +
  ylab('Transcripts') +
  xlab(NULL) +  # Remover o rótulo do eixo X
  color_palette +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Gráfico de dispersão
scatter_plot <- ggplot(as.data.frame(filtered_origin_SSPO_SBAR), aes(x=CPM_SSPO, y=CPM_SBAR)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. spontaneum') +
  ylab('Log10 CPM S. barberi') +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Extrair a legenda do gráfico de dispersão
legend <- get_legend(
  scatter_plot + theme(legend.position = "right")
)

# Remover a legenda do gráfico de dispersão principal
scatter_plot <- scatter_plot + theme(legend.position = "none")

# Combinar os gráficos lado a lado com a legenda
combined_plots <- plot_grid(
  scatter_plot, bar_plot, legend,
  ncol = 3, rel_widths = c(3, 2, 0.8)  # Ajuste das larguras relativas para evitar sobreposição
  #  align = "h"
)

# Adicionar o título e subtítulo principal
title <- ggdraw() + 
  draw_label("Common transcripts between SSPO and SBAR", fontface = 'bold', size = 22, hjust = 0.5) +
  draw_label("Visualization of conserved/common transcripts", fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
  theme(plot.margin = margin(t = 10, b = 20))

# Combinar o título e os gráficos
final_plot <- plot_grid(
  title, combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Salvar a combinação como PNG
ggsave("common_SSPO_SBAR.png", plot = final_plot, width = 18, height = 11, dpi = 500, bg="white")






#################################### SSPO and SBAR - coding

# Contagem de transcritos por função
function_counts <- filtered_origin_SSPO_SBAR_coding %>%
  group_by(Origin) %>%
  summarise(Count = n())

# Definir a paleta de cores comum
color_palette <- scale_fill_brewer(palette = "Set1")

# Gráfico de barras estilizado para as contagens de funções
bar_plot <- ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
  geom_bar(stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=Count), vjust=-0.5, size=5) +  # Adiciona as contagens em cada barra
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey")
  ) +
  ylab('Transcripts') +
  xlab(NULL) +  # Remover o rótulo do eixo X
  color_palette +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Gráfico de dispersão
scatter_plot <- ggplot(as.data.frame(filtered_origin_SSPO_SBAR_coding), aes(x=CPM_SSPO, y=CPM_SBAR)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. spontaneum') +
  ylab('Log10 CPM S. barberi') +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Extrair a legenda do gráfico de dispersão
legend <- get_legend(
  scatter_plot + theme(legend.position = "right")
)

# Remover a legenda do gráfico de dispersão principal
scatter_plot <- scatter_plot + theme(legend.position = "none")

# Combinar os gráficos lado a lado com a legenda
combined_plots <- plot_grid(
  scatter_plot, bar_plot, legend,
  ncol = 3, rel_widths = c(3, 2, 0.8)  # Ajuste das larguras relativas para evitar sobreposição
  #  align = "h"
)

# Adicionar o título e subtítulo principal
title <- ggdraw() + 
  draw_label("Common coding transcripts between SSPO and SBAR", fontface = 'bold', size = 22, hjust = 0.5) +
  draw_label("Visualization of conserved/common transcripts", fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
  theme(plot.margin = margin(t = 10, b = 20))

# Combinar o título e os gráficos
final_plot <- plot_grid(
  title, combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Salvar a combinação como PNG
ggsave("common_SSPO_SBAR_coding.png", plot = final_plot, width = 18, height = 11, dpi = 500, bg="white")







#################################### SSPO and SBAR - non-coding

# Contagem de transcritos por função
function_counts <- filtered_origin_SSPO_SBAR_noncoding %>%
  group_by(Origin) %>%
  summarise(Count = n())

# Definir a paleta de cores comum
color_palette <- scale_fill_brewer(palette = "Set1")

# Gráfico de barras estilizado para as contagens de funções
bar_plot <- ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
  geom_bar(stat="identity", position="dodge", width=0.5) +
  geom_text(aes(label=Count), vjust=-0.5, size=5) +  # Adiciona as contagens em cada barra
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.2, color = "grey")
  ) +
  ylab('Transcripts') +
  xlab(NULL) +  # Remover o rótulo do eixo X
  color_palette +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Gráfico de dispersão
scatter_plot <- ggplot(as.data.frame(filtered_origin_SSPO_SBAR_noncoding), aes(x=CPM_SSPO, y=CPM_SBAR)) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
  scale_color_brewer(palette="Set1") +
  xlab('Log10 CPM S. spontaneum') +
  ylab('Log10 CPM S. barberi') +
  scale_x_log10() +
  scale_y_log10() +
  geom_density2d(size=0.3) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

# Extrair a legenda do gráfico de dispersão
legend <- get_legend(
  scatter_plot + theme(legend.position = "right")
)

# Remover a legenda do gráfico de dispersão principal
scatter_plot <- scatter_plot + theme(legend.position = "none")

# Combinar os gráficos lado a lado com a legenda
combined_plots <- plot_grid(
  scatter_plot, bar_plot, legend,
  ncol = 3, rel_widths = c(3, 2, 0.8)  # Ajuste das larguras relativas para evitar sobreposição
  #  align = "h"
)

# Adicionar o título e subtítulo principal
title <- ggdraw() + 
  draw_label("Common non-coding transcripts between SSPO and SBAR", fontface = 'bold', size = 22, hjust = 0.5) +
  draw_label("Visualization of conserved/common transcripts", fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
  theme(plot.margin = margin(t = 10, b = 20))

# Combinar o título e os gráficos
final_plot <- plot_grid(
  title, combined_plots,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Salvar a combinação como PNG
ggsave("common_SSPO_SBAR_noncoding.png", plot = final_plot, width = 18, height = 11, dpi = 500, bg="white")










#write.table(countsOrigin, sep = "\t", file = "speciesOfOriginPanTranscriptome.csv")

#tx2gene<-read.delim("txp.group.tsv",sep=',',header=T)
#head(tx2gene)
#length(unique(tx2gene$GeneID))
#length(unique(tx2gene$TranscriptID))