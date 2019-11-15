# Import libraries 
library(rjson)
library(tidyverse)
library(MatrixEQTL)
library(readr)
library(data.table)
library(ggrepel)
library(viridis)
library(hrbrthemes)

# Set working directory 
setwd("~/Documents/LociAnalysis")

################ IMPORT ###################

# Import options 
options <- fromJSON(file = 'options.json')

# Import CHR14
chr14.data <- fread(options$Outputs$CHR14, sep = ',')

# Reorder
chr14.data.ordered <- chr14.data %>% arrange(SNP)

################ PLOT ####################

# Create array 
x <- c(1:length(chr14.data$SNP))
y <- -log10(chr14.data.ordered$pvalue)
xlab <- chr14.data.ordered$SNP

# Create labels
labels <- rep(x = NA, length(chr14.data$SNP))
genes <- unique(chr14.data$gene)
for (gene in genes){
  chr14.data.gene <- chr14.data.ordered[chr14.data.ordered$gene == gene,]
  minpval <- chr14.data.gene[chr14.data.gene$pvalue == min(chr14.data.gene$pvalue),]
  idx <- which(chr14.data.ordered$SNP == minpval$SNP[1] & chr14.data.ordered$gene == gene)
  labels[idx[1]] <- gene
}
chr14.data.ordered$labels <- labels

# Plot CHR14
data.plot <- data.frame(x = x, y = y)
ggplot(data = data.plot, aes(x, y, label = labels)) + geom_point(col = 'blue') + ylim(-0,120) +
  geom_text_repel( nudge_y = 300, direction = "x", angle = 90,
                     vjust = 0, segment.size = 0.2, size=4) +
  labs(x='Position',y='-log10(pval)', title='Chromosome 14') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

############## COMPARE TO NMDA ################

NMDA.data <- fread(options$file$NMDA)
NMDA.chr14 <- NMDA.data[NMDA.data$chr == 14,]
chr14.data.ordered$NMDA <- rep("Not Match", length(chr14.data.ordered$SNP))
for (pos in NMDA.chr14$pos){
  idx <- which(chr14.data.ordered$SNP == pos)
  if (length(idx) > 1){
    i <- which(chr14.data.ordered$pvalue[idx] == min(chr14.data.ordered$pvalue[idx]))
    chr14.data.ordered$NMDA[idx[i]] <- "Match"
  }else{
    chr14.data.ordered$NMDA[idx] <- "Match"
  }
}
tst <- chr14.data.ordered[chr14.data.ordered$NMDA == "Match",]
data.plot <- data.frame(x = x, y = y, NMDA =  chr14.data.ordered$NMDA, labels = labels)

# Filter data
chr14.data.filtered <- data.plot %>% 
  group_by(NMDA) %>%
  filter(NMDA == 'Match') %>%
  ungroup()

# Plot
ggplot() + geom_point(data = data.plot, aes(x, y), color = 'blue', alpha = '0.1') +
  geom_point(data = chr14.data.filtered, aes(x,y), color = 'red') +
  ylim(-0,120) +
  geom_text_repel(aes(x,y,label=labels), nudge_y = 300, direction = "x", angle = 90,
                   vjust = 0, segment.size = 0.2, size=4) +
  labs(x='Position',y='-log10(pval)', title='Chromosome 14 matches with NMDA (red)') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position = "top")

################## PLOT ZOOM IN ####################

lowPos <- min(chr14.data.ordered[chr14.data.ordered$NMDA == 'Match',]$SNP)
highPos <- max(chr14.data.ordered[chr14.data.ordered$NMDA == 'Match',]$SNP)

data.plot.zoom <- chr14.data.ordered[chr14.data.ordered$SNP > (lowPos -1000) &
                                      chr14.data.ordered$SNP < (highPos + 1000),]
x <- c(1:length(data.plot.zoom$SNP))
y <- -log10(data.plot.zoom$pvalue)
data.plot <- data.frame(x = x, y = y, NMDA =  data.plot.zoom$NMDA, labels = data.plot.zoom$labels)
chr14.data.filtered.zoom <- data.plot %>% 
  group_by(NMDA) %>%
  filter(NMDA == 'Match') %>%
  ungroup()

# Plot
ggplot() + geom_point(data = data.plot, aes(x, y), color = 'blue', alpha = '0.1') +
  geom_point(data = chr14.data.filtered.zoom, aes(x,y), color = 'red') +
  geom_text_repel(aes(x,y,label=data.plot.zoom$labels), nudge_y = 300, direction = "x", angle = 90,
                  vjust = 0, segment.size = 0.2, size=4) +
  labs(x='Position',y='-log10(pval)', title='Chromosome 14 zoomed in NMDA matches') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position = "top")
