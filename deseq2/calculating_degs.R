

rm(list = ls())
library(tidyverse)
library(ggpubr)
library(reshape2)
library(ggbeeswarm)
library(DESeq2)

setwd('C:\\Users\\Ben\\Dropbox\\bilbo_lab_spr2020\\microglia-seq_website\\microglia-seq\\deseq2')

sem <- function(x) sqrt(var(x)/length(x))

coldata <- read.csv('col_data_combined.csv')
df <- read.csv('GSE99622_hanamsagar2017_raw_reads_p60.csv', row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = coldata,
                              design = ~ tx)

dds <- DESeq(dds)
res <- results(dds)

write.csv(res, 'combined_degs.csv')
  
#################################################################################################
#################################################################################################
#################################################################################################


coldata <- read.csv('col_data_males.csv')
df <- read.csv('GSE99622_hanamsagar2017_raw_reads_p60_male.csv', row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = coldata,
                              design = ~ tx)

dds <- DESeq(dds)
res <- results(dds)

write.csv(res, 'male_degs.csv')



#################################################################################################
#################################################################################################
#################################################################################################


coldata <- read.csv('col_data_females.csv')
df <- read.csv('GSE99622_hanamsagar2017_raw_reads_p60_female.csv', row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = coldata,
                              design = ~ tx)

dds <- DESeq(dds)
res <- results(dds)

write.csv(res, 'female_degs.csv')




