# Gene ontology analysis with clusterProfiler
library(clusterProfiler)

setwd('/Users/juanjovel/jj/data_analysis/walterLopez/Carolina_Tomato_Nematode/new_results_paired-end')

# List all files ending with '_sorted_genes.txt'
files <- list.files(pattern = "_sorted_genes\\.txt$")

# Read universe and gene-to-GO data just once as they are common for all analyses
universe     <- read.table("universe_sorted.txt", header = T)
universe     <- as.character(universe$genes)
gene_to_go   <- read.table("gene-to-GO-ALL_sorted.tsv", header = T, sep = '\t')

# Function for GO enrichment

GOEnrichment <- function(genes_of_interest, universe, gene_to_GO, ontology){
  
  if (ontology == 'ALL') {
    
    TERM2GENE = gene_to_GO[,c('GOterm', 'LocusID')]
    TERM2GENE$LocusID = substr(TERM2GENE$LocusID, start = 1, stop = 16)
    TERM2NAME = gene_to_GO[,c('GOterm', 'GOdesc')]
    
  } else {
    
    TERM2GENE = gene_to_GO[gene_to_GO$Ontology == ontology, c('GOterm', 'LocusID')]
    TERM2GENE$LocusID = substr(TERM2GENE$LocusID, start = 1, stop = 16)
    TERM2NAME = gene_to_GO[gene_to_GO$Ontology == ontology, c('GOterm', 'GOdesc')]
    
  }
  
  results <- enricher(gene = genes_of_interest,
                                       universe = universe,
                                       TERM2GENE = TERM2GENE,
                                       TERM2NAME = TERM2NAME,
                                       pAdjustMethod = 'BH',
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.1,
                                       minGSSize = 7,
                                       maxGSSize = 500)
  
  return(results)
}

# Iterate over each file
for(file in files) {
  genes <- read.table(file, header = T)
  genes <- gsub("\\.\\d\\.\\d+", "", genes$genes)
  
  results <- GOEnrichment(genes_of_interest = genes, universe = universe, 
                          gene_to_GO = gene_to_go, ontology = 'ALL')
  
  # Generate dynamic file names for output based on input file name
  output_base <- gsub("_sorted_genes\\.txt", "", file)
  results_file <- paste0(output_base, "_q0.05_GO_results.txt")
  dotplot_file <- paste0(output_base, "_q0.05_dotplot.png")
  
  # Write results and generate dotplot
  write.table(results, results_file, sep = '\t', quote = F)
  
  # Create the ggplot object for the dotplot
  plot_object <- dotplot(results, showCategory = 10, font.size = 8)
  
  # Save the plot
  png(dotplot_file)
  print(plot_object)  # Explicitly print the plot object
  dev.off()
}