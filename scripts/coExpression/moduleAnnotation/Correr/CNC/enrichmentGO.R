library(topGO)
library(ggplot2)

rm(list=ls())
setwd(getwd())

args <- commandArgs(trailingOnly = TRUE)                                                 # args from command line

if (length(args) != 3) {                                                                 # verify number of args
  stop("Invalid number of arguments. Please provide the file names.")
}

results_path <- "results/"                                                               # output directory 
dir.create(results_path, recursive = TRUE)                                               # create recursive dir

load_files <- function(clusterSize, GO_universe, Cliques){                               # function to load files
  number_clusters <<- read.table(clusterSize, header = F)                                # read clusters size
  geneID2GO <<- readMappings(file = GO_universe)                                         # read list of annotated genes for all the network
  genes_in_graph <<- as.data.frame(read.table(Cliques, header = F))                      # read genes in graph
  colnames(genes_in_graph) <<- c("gene", "category", "module_No")
  geneID2GO_filtered <<- geneID2GO[genes_in_graph[,1]]
  geneNames <<- genes_in_graph$gene
}

load_files(args[1], args[2], args[3])                                                    # call the load_files function with the provided arguments

anot_modules <- function(Module_no, results_path){                                       # function make annot fun
  myInterestingGenes <- genes_in_graph$gene[genes_in_graph$module_No==Module_no]         # get interesting genes from each module
  geneList <- as.factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  
  GOdata <- new("topGOdata",                                                             # create GO dataset 
                 ontology = "BP",
                 allGenes = geneList,
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GO_filtered)
  
  allGO=usedGO(GOdata)
  
  Classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")                                  # run fisher classic
  #resultsWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic')              # make results  table
  #table <- GenTable(GOdata, Classic = resultClassic, Weight01 = resultsWeight01, topNodes = length(allGO), orderBy = 'Classic')
  
  table$Classic[table$Classic == "< 1e-30"] <- 1e-30                                                       # replace the values "< 1e-30" with a very small number ( < 1e-30 cannot be corrected with BH)
  table$Classic <- as.numeric(table$Classic)                                                               # convert the Classic column to numeric
  
  #table1 <- filter(table, Classic < 0.05)                                                                 # filter not significant values for classic algorithm
  
  #p.adj <- round(p.adjust(table1$Classic,method="BH"),digits = 4)                                         # performing BH correction on p values FDR
  p.adj <- round(p.adjust(table$Classic,method="BH"),digits = 4)
  
  #all_res_final <<- cbind(table1,p.adj)                                                                   # create the file with all the statistics from GO analysis
  all_res_final <- cbind(table,p.adj)
  all_res_final <- all_res_final[order(all_res_final$p.adj),]                                              # order results
  
  results.table.p = all_res_final[which(all_res_final$Classic <=0.05),]                                    # get list of significant GO before multiple testing correction
  results.table.bh = all_res_final[which(all_res_final$p.adj  <=0.05),]
  
  # save top 50 ontologies sorted by adjusted pvalues
  write.table(all_res_final[1:50,], file = paste0(results_path, "module_", Module_no, ".csv"), quote=FALSE, row.names=FALSE, sep = ",")

  module <- paste0("module_", Module_no)                                                                   # create plots for each module
  topGO_all_table <- all_res_final
  colnames(topGO_all_table) <- c("GO.ID","Term","Annotated","Significant","Expected","Classic","p.adj")   
  topGO_all_table <- topGO_all_table[order(topGO_all_table$Classic),]
  
  ntop <- 30                                                                                               # only plot top 30 out 50 terms
  ggdata <- topGO_all_table[1:ntop,]
  #ggdata <- ggdata[!duplicated(ggdata$Term), ]                                                            # remove duplicated Terms                               
  #ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))                                           # fixes order
  #ggdata$Classic <- as.numeric(ggdata$Classic) + 0.000001                                                 # Add small number for log operations (log10 must be > 0) 
  ggdata
  
  max_y <- max(-log10(ggdata$Classic))                                                                     # calculate the maximum value on the y-axis
  div_points <- quantile(-log10(ggdata$Classic), probs = c(0.95, 0.5, 0.05))                               # calculate quantiles for dividing the plot
  
  ggplot(ggdata,
         aes(x = ggdata$Term, y = -log10(Classic), size = -log10(Classic), fill = -log10(Classic))) +
    
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5)) +
    scale_fill_continuous(low = 'royalblue', high = 'red4') +
    
    xlab('') + ylab('Enrichment score') +
    labs(
      title = module,
      subtitle = 'Top 30 terms ordered by Fisher Exact p-value',
      caption = paste('Cut-off lines drawn at equivalents of p =', round(div_points[3], 1), ',', round(div_points[2], 1), ',', round(div_points[1], 1))) +

    geom_hline(yintercept = div_points,                                                                      # draw vertical lines
               linetype = c("dotted", "longdash", "solid"),
               colour = c("black", "black", "black"),
               size = c(0.5, 1.5, 3)) +
    
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      
      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
      axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
      axis.title = element_text(size = 12, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.line = element_line(colour = 'black'),

      legend.key = element_blank(),                                                                          # removes legend the border
      legend.key.size = unit(1, "cm"),                                                                       # sets overall area/size of the legend
      legend.text = element_text(size = 14, face = "bold"),                                                  # text size
      title = element_text(size = 14, face = "bold")) +
    
    coord_flip()
  
  plot_path <- paste0(results_path, "module_", Module_no, "_GO_30Terms_Fisher", ".png")
  
  ggplot2::ggsave(plot_path,
                  device = NULL,
                  height = 8.5,
                  width = 12)
}

# iterate over all modules that have more than 5 genes
for (i in 1:nrow(number_clusters)){                                                                          # for testing: for (i in 1:2){
  if(number_clusters[i,1]>=5){
    print(i)
    anot_modules(i, results_path)
  }
}