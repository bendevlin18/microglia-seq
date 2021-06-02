

######################################
###                                ###
###     global.R file for          ###
###       microglia-seq            ###
###        application             ###
###                                ###
######################################

##### Load in necessary packages #####

library(shiny)
library(ggplot2)
library(dplyr)
library(shinydashboard)
library(plotly)
library(ggpubr)
library(DT)
library(tidyverse)
library(shinycssloaders)



### Load in the data and generate objects ###
df <- read.csv('GSE99622_hanamsagar2017_tpm_v2.csv')
df2 <- readRDS('gene_ensembl_ids_opentarget.rds')
df3 <- read.csv('unique_data_index.csv')
df4 <- read.csv('unique_data_index_gene_list.csv')


genes <- unique(df$gene)
std_err <- function(x) sd(x) / sqrt(length(x))
options(spinner.color='#222d32', spinner.color.background = '#ecf0f5', spinner.type = 2)

graph_theme_settings <- list(
    xlab('\r\nAge/Treatment'),
    ylab('TPM (Transcripts per million) \r\n\r\n'),
    theme_classic(),
    theme(rect = element_rect(fill = 'transparent'), 
          text = element_text(size = 13),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = .6),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(color = NA, col = 0),
          axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
          strip.background = element_rect(color="transparent", fill="transparent", linetype="solid")))

dotplot_theme <- list(
  ylab('\r\n *Note: Size of Dots = Standard Error'),
  xlab('Gene'),
  theme_bw(),
  theme(rect = element_rect(fill = 'transparent'), 
        text = element_text(size = 13), 
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = .6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.ticks.length = unit(.2, 'cm'),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(color = NA, col = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 75, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 0, r = 75, b = 0, l = 0)),
        panel.spacing = unit(1, 'cm'))
  )


