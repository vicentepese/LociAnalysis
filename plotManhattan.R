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
options <- fromJson(file = 'options.json')