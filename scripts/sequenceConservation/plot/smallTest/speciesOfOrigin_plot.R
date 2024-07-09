library(ggplot2)
library(dplyr)
library(cowplot)

# install ggplot2 additional svglite to save svg files!
# install.packages('svglite')

rm(list=ls())

DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/sequenceConservation/plot/smallTest"
DIR="/Storage/data1/felipe.peres/Sugarcane_ncRNA/11_lncRNASequenceConservation/GenomicReads/PLOT"

setwd(DIR)

countsOrigin<-read.delim(file = "5000_speciesOfOriginPanTranscriptome.tsv",header=T,row.names=1)
head(countsOrigin)
colnames(countsOrigin)

origNumberGenes=dim(countsOrigin)[1]
dim(countsOrigin)
head(countsOrigin)

table(countsOrigin$Origin, useNA = 'always')*100/sum(table(countsOrigin$Origin, useNA = 'always'))

quantile(countsOrigin$CPM_SBAR, c(.1,.2,.5,.8,.9,.99))
quantile(countsOrigin$CPM_SSPO, c(.1,.2,.5,.8,.9,.99))
quantile(countsOrigin$CPM_SOFF, c(.1,.2,.5,.8,.9,.99))

colnames(countsOrigin)

colors <- c("SOFF" = "red", "SSPO" = "blue", "SBAR" = "green", "Common" = "yellow", "UNK" = "grey",
            "CommonSOFF_SBAR" = "brown", "CommonSSPO_SBAR" = "orange", "CommonSSPO_SOFF" = "purple")

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


# filter transcript by function
filter_transcripts <- function(data, function_type) {
  subset(data, Function == function_type)
}

# filter transcript by origin
count_transcripts <- function(data, function_type) {
  data %>%
    filter(Function == function_type) %>%
    group_by(Origin) %>%
    summarise(Count = n())
}

# create histogram by function
create_bar_plot <- function(data, function_types) {
  function_counts <- count_transcripts(data, function_types)
  ggplot(function_counts, aes(x=Origin, y=Count, fill=Origin)) +
    geom_bar(stat="identity", position="dodge", width=0.5) +
    geom_text(aes(label=Count), vjust=-0.5, size=5) +
    scale_fill_manual(values = colors) + 
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(size = 0.2, color = "grey")
    ) +
    ylab('Transcripts') +
    xlab(NULL) +
    #scale_fill_brewer(palette = "Set1") +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14)
    )
}

# create scatterplot by function
create_scatter_plot <- function(data, xy_columns) {
  ggplot(as.data.frame(data), aes_string(x = xy_columns[1], y = xy_columns[2])) +
    theme_minimal(base_size = 20) +
    theme(legend.position = "none") +
    geom_jitter(aes(colour=Origin), alpha=0.2, size=1.5) +
    scale_color_manual(values = colors) + 
    #scale_color_brewer(palette="Set1") +
    xlab(paste('Log10 CPM', gsub('CPM_', '', xy_columns[1]))) +
    ylab(paste('Log10 CPM', gsub('CPM_', '', xy_columns[2]))) +
    scale_x_log10() +
    scale_y_log10() +
    geom_density2d(size=0.3) +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14)
    )
}

# combine histogram with scatterplot 
combine_plots <- function(bar_plot, scatter_plot, title_text, subtitle_text) {
  legend <- get_legend(scatter_plot + theme(legend.position = "right"))
  scatter_plot <- scatter_plot + theme(legend.position = "none")
  
  combined_plots <- plot_grid(
    scatter_plot, bar_plot, legend,
    ncol = 3, rel_widths = c(3, 2, 0.8)
  )
  
  title <- ggdraw() + 
    draw_label(title_text, fontface = 'bold', size = 22, hjust = 0.5) +
    draw_label(subtitle_text, fontface = 'plain', size = 18, hjust = 0.5, vjust = 2) +
    theme(plot.margin = margin(t = 10, b = 20))
  
  plot_grid(title, combined_plots, ncol = 1, rel_heights = c(0.1, 1))
}

# save plots
save_plot <- function(plot, filename) {
  #ggsave(filename, plot = plot, width = 18, height = 11, dpi = 500, bg="white")
  ggsave(filename, plot = plot, device = "svg", width = 18, height = 11, units = "in", dpi = 300)
  
}

# iterate over datasets and save plots
generate_plots <- function(datasets, function_types, xy_mappings) {
  for (i in seq_along(datasets)) {
    filtered_data <- datasets[[i]]
    xy_columns <- xy_mappings[[i]]
    
    for (function_type in names(function_types)) {
      filtered_subset <- if (is.null(function_types[[function_type]])) {
        filtered_data
      } else {
        filter_transcripts(filtered_data, function_types[[function_type]])
      }
      
      #bar_plot <- create_bar_plot(filtered_subset, function_types)
      bar_plot <- create_bar_plot(countsOrigin, function_types)
      scatter_plot <- create_scatter_plot(filtered_subset, xy_columns)
      
      title_text <- paste("Common", function_type, "transcripts between", names(xy_mappings)[i])
      subtitle_text <- "Visualization of conserved/common transcripts"
      
      final_plot <- combine_plots(bar_plot, scatter_plot, title_text, subtitle_text)
      #filename <- paste0("common_", names(xy_mappings)[i], "_", function_type, ".png")   # save as png
      filename <- paste0("common_", names(xy_mappings)[i], "_", function_type, ".svg")    # save as svg
      
      #print(final_plot)
      #readline(prompt = "")                                                              # press enter in the console to show more plots
      
      save_plot(final_plot, filename)
    }
  }
}

datasets <- list(filtered_origin_SOFF_SSPO, filtered_origin_SOFF_SBAR, filtered_origin_SSPO_SBAR)
function_types <- list("all" = NULL, "coding" = "protein-coding", "noncoding" = "ncRNA", "lncRNA" = "lncRNA")
#function_types <- list("lncRNA" = "lncRNA")

xy_mappings <- list(
  filtered_origin_SOFF_SSPO = c("CPM_SOFF", "CPM_SSPO"),
  filtered_origin_SOFF_SBAR = c("CPM_SOFF", "CPM_SBAR"),
  filtered_origin_SSPO_SBAR = c("CPM_SSPO", "CPM_SBAR")
)

generate_plots(datasets, function_types, xy_mappings)
