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

############# IMPORT ###############

# Import options 
options <- jsonlite::fromJSON('options.json')

# Get list of files 
plot.files <- list.files(options$folder$filtChrPlot)

# Sort
plot.files <- plot.files[order(nchar(plot.files), plot.files)]

############ PLOT ################

pvals <- c()
colors <- c()
chr <- c()
i <- 1
for (plot.file in plot.files){
  
  # Import file 
  plot.data <- read_csv(file = paste(options$folder$filtChrPlot, plot.file, sep = ''), col_names = FALSE)
  colnames(plot.data) <- c("pos","gene","pval")
  
  # Order positions chronologically
  plot.data <- plot.data[order(plot.data$pos),]
  
  # Keep pvals
  pvals <- c(pvals, plot.data$pval)
  
  # Keep chromosome number 
  chr <- c(chr, rep(strsplit(strsplit(plot.file,'_')[[1]][2], '.', fixed = TRUE)[[1]][1], length(plot.data$pos)))
  
  # Colors
  if (i == 1){
    i <- 0
    colors <- c(colors, rep('#F8766D', length(plot.data$pos)))
  } else{
    i <- 1 
    colors <- c(colors, rep('#00BFC4', length(plot.data$pos)))
  }
  
}

# Create arrays 
x <- c(1:length(pvals))

# Creat datatable
plot.data <- data.table('pval' = -log10(pvals), 'pos' = x, 'colors' = colors, 'chr' = chr)
axisdf <- plot.data %>% group_by(chr) %>% summarize(center=( max(pos) + min(pos)) / 2 )

# Plot
png("GWAS.png", width = 1920, height = 1080)
ggplot() + geom_point(data = plot.data, aes(pos, pval, color = as.factor(colors)), alpha = 0.3) +
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center ) + 
  labs(x='Chromosome',y='-log10(pval)', title='GWAS') +
  geom_hline(yintercept = -log10(5e-8)) + 
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
dev.off()

############## MANHATTAN PLOT TEMPLATE ################

ggManh=function(metaout){
  df.temp= metaout%>% group_by(chr) %>% 
    summarise(chr.len= max(pos)) %>% 
    mutate(tot=cumsum(as.numeric(chr.len))-chr.len) %>% 
    select(-chr.len) %>% left_join(metaout, ., by=c("chr"="chr")) %>% arrange(chr, pos) %>% mutate(BPcum=pos+tot)
  axisdf <- df.temp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  setDT(df.temp)
  a=ggplot(df.temp[pvals > 1e-20], aes(x=BPcum, y=-log10(pval))) +
    geom_point(aes(color=as.factor(chr)), alpha=0.3, size=1.4)+geom_hline(yintercept = -log10(5e-8))+
    scale_x_continuous(label = axisdf$chr, breaks= axisdf$center )+ labs(x = "Chromosome")+theme_bw(base_size = 22) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) 
  return(a)
}
