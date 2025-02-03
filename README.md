# Multi-genotype analyses of long-ncRNAs in Sugarcane

This work addresses gaps in the understanding of non-coding RNAs (ncRNAs) and long non-coding RNAs (lncRNAs) in sugarcane, providing a comprehensive exploration of the variability and functional roles of these transcripts. 

Developed as part of a master's thesis at [Center for Nuclear Energy in Agriculture - University of São Paulo](http://www.cena.usp.br/) in 2024, conducted in the [Computational, Evolutionary and Systems Biology Laboratory](https://labbces.cena.usp.br/), this research leveraged publicly available RNA-Seq data to characterize nc/lncRNAs in sugarcane, uncovering their classification, conservation, co-expression patterns, and potential functions.

>[!NOTE] 
>This repository serves as a resource for accessing the code, technical decisions, and results generated during the study.

## Approach and Methods

* *How can we explore publicly available RNA-Seq data for sugarcane to identify nc/lncRNAs?*
  * [Predicting putative lncRNAs with CPC2](https://github.com/labbces/sugarcane_RNAome/wiki/Predicting-putative-lncRNAs-with-CPC2)
  * [Predicting putative lncRNAs with PLncPRO](https://github.com/labbces/sugarcane_RNAome/wiki/Predicting-putative-lncRNAs-with-PLncPRO)
  * [Predicting putative lncRNAs with RNAplonc](https://github.com/labbces/sugarcane_RNAome/wiki/Predicting-putative-lncRNAs-with-RNAplonc)
  * [Aggregating the consensus of putative ncRNAs predicted by the 3 tools](https://github.com/labbces/sugarcane_RNAome/wiki/Aggregating-the-consensus-of-putative-ncRNAs-predicted-by-the-3-tools)

* *What are the main families of nc/lncRNAs present in sugarcane, and how can they be classified?*
  * [Annotation of ncRNAs families](https://github.com/labbces/sugarcane_RNAome/wiki/Annotation-of-ncRNAs-families)

* *pan-RNAome inference: Are there conserved nc/lncRNAs shared among different genotypes?*
  * [Clustering sequences by identity](https://github.com/labbces/sugarcane_RNAome/wiki/Clustering-sequences-by-identity)
  * [Inference of the pan-RNAome](https://github.com/labbces/sugarcane_RNAome/wiki/Inference-of-the-pan%E2%80%90ncRNAome)

* *What is the probable origin of the lncRNAs in modern sugarcane hybrids, and how do they trace back to the ancestral species the Saccharum complex?*
  * [Estimating lncRNA origin](https://github.com/labbces/sugarcane_RNAome/wiki/Estimating-lncRNA-origin)

* *What are the expression and co-expression patterns of nc/lncRNAs in sugarcane across different conditions, and how can co-expression networks be utilized to annotate and infer the functional roles of these nc/lncRNAs?*
  * [Contrasting genotypes selection](https://github.com/labbces/sugarcane_RNAome/wiki/Contrasting-genotypes-selection)
  * [Calculating the RNA expression for each sample](https://github.com/labbces/sugarcane_RNAome/wiki/Calculating-the-RNA-expression-for-each-sample)
  * [Filtering the RNAs expression matrix of the contrasting genotypes](https://github.com/labbces/sugarcane_RNAome/wiki/Filtering-the-RNAs-expression-matrix-of-the-contrasting-genotypes)
  * [Filtering the RNAs expression matrix by function](https://github.com/labbces/sugarcane_RNAome/wiki/Filtering-the-RNAs-expression-matrix-by-function)
  * [Analyzing RNA expression patterns after filtering](https://github.com/labbces/sugarcane_RNAome/wiki/Analyzing-RNA-expression-patterns-after-filtering)
  * [Principal Component Analysis of lncRNAs expression in contrasting genotypes](https://github.com/labbces/sugarcane_RNAome/wiki/Principal-Component-Analysis-of-lncRNAs-expression-in-contrasting-genotypes)
  * [Computing gene pair relationships with Pearson correlation](https://github.com/labbces/sugarcane_RNAome/wiki/Computing-gene-pair-relationships-with-Pearson-correlation)
  * [Analyzing Network Metrics](https://github.com/labbces/sugarcane_RNAome/wiki/Analyzing-Network-Metrics)
  * [Analyzing MCL clusterings (effects of Inflation value)](https://github.com/labbces/sugarcane_RNAome/wiki/Analyzing-MCL-clusterings-(effects-of-Inflation-value))
  * [Finding conserved modules of co-expression (complete subgraphs - cliques)](https://github.com/labbces/sugarcane_RNAome/wiki/Finding-conserved-modules-of-co%E2%80%90expression-(complete-subgraphs-%E2%80%90-cliques))
  * [Visualizing co-expressed gene modules in contrasting genotypes](https://github.com/labbces/sugarcane_RNAome/wiki/Visualizing-co%E2%80%90expressed-gene-modules-in-contrasting-genotypes)
  * [Gene Ontology enrichment analysis: exploring biological function of modules](https://github.com/labbces/sugarcane_RNAome/wiki/Gene-Ontology-enrichment-analysis:-exploring-biological-function-of-modules)

## Financial Support

This research was supported by Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq), process number: 130748/2022-6

## Citation
```bash
@article {,
      title={}, 
      author={},
      year={},
      eprint={},
      archivePrefix={},
      primaryClass={},
      url={}, 
}
```
