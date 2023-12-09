library (dplyr)
library (ggplot2)
library (wesanderson)

setwd('/Users/juanjovel/jj/data_analysis/walterLopez/Carolina_Tomato_Nematode/new_results_paired-end/boxplots')
ids_file <- 'all.ids'
file <- file.path(getwd(), 'all_samples_tpms.tsv')

metafile <- 'metadata.txt'
metadata <- read.table(metafile, header = T, sep = '\t', row.names = 1)
counts_data_df <- read.table(file, header = T, sep = '\t', row.names = 1)

# Import a list of transcripts IDs that were found differentially expressed as above
transc_ids <- read.table (ids_file, header = F)

counter = 0

for (i in 1:nrow(transc_ids)) {
  
  counter = counter + 1
  trID <- transc_ids [i,]
  match <- grep (trID, readLines(file), value = T)
  sel_row   <- counts_data_df[trID,]
  row_names <- colnames(counts_data_df)
  group     <- metadata$group
  counts <- as.numeric(sel_row[1:ncol(sel_row)])

  df = data.frame(row.names = row_names, 
                  condition = group,
                  counts = counts)
  
  gene <- row.names(sel_row)
  gene.name <- paste(gene, "Expression", sep = " ")
  fileName <- paste (gene, "pdf", sep = ".")
  message <- paste("Plot#:", counter, "Current plot:", fileName, sep = ' ')
  print (message )
  pdf(file = fileName)
  p <- ggplot(df, aes(condition, counts))  
    
  print (p + geom_boxplot(notch=F,  aes(fill = factor(group))) + theme_bw() +
    #scale_fill_brewer(palette="Dark2") +
      scale_fill_manual(values=wes_palette(n=4, name="Darjeeling2")) +
      ggtitle (gene.name) +
      xlab("Condition") +
      ylab("Abundance (TPM)") +
      #scale_y_log10 () +
      # Control the text features
      theme( axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
           axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
           axis.title.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
           axis.title.y = element_text(colour="black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))
  )
  dev.off()
}
