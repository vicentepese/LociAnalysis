# Import
library(readr)
library(rjson)


# Set working directory 
setwd("~/Documents/LociAnalysis")

# Read options 
options <- fromJSON(file = 'options.json')

# Import protein data
genes <- read_csv(options$file$geneData)
genes <- colnames(genes)
genes <- genes[3:length(genes)]

# Import SNPS
snps <- read_csv('tst.csv', col_names = FALSE)
colnames(snps) <- c('SNP','Gene','Statistic','pvalue','FDR','Beta')
snps$log10pval <- -log10(snps$pvalue)

# Position range 
pos_range <- c(105863196, 105863260)
pos_range <- pos_range + c(-500e3, 500e3)
pos_array <- c(pos_range[1]:pos_range[2])
pval_array <- rep(NA, length(pos_array))

i <- 1
for (snp in snps$SNP){
  
  if (grepl('rs', snp)){
    # Get position
    pos <- strsplit(snp, ':')[[1]][2]
    
    # Find idx
    idx <- which(pos_array == as.numeric(pos))
    
    # Add pvalue to pval array 
    pval_array[idx] <- snps$log10pval[i] 
  }
  
  i <- i+1
  
}


# GGPLOT
