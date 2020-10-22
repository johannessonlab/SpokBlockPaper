#!/usr/bin/env Rscript

### TSD histogram in the *Podospora* species
#############################################################################
# 
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020-09-22
# Version 1
# =======================================

library(ggplot2)
library(dplyr)

# ============================
# Data
# ============================
## Snakemake
countsfile <- snakemake@input$counts
kmers <- snakemake@params$kmers

## Read and clean data
counts <- read.table(countsfile)
names(counts) <- c("Kmer", "Count")

## Get the relevant kmers
marked <- counts %>% filter(Kmer %in% kmers)

distkmer <- ggplot(counts, aes(x = Count)) + geom_histogram() +
  theme_bw() +
  xlab("Kmer frequency") + 
  geom_vline(data = marked, aes(xintercept = Count, colour = Kmer), size=1)

ggsave(plot = distkmer, snakemake@output$hist, width = 5, height = 3)
