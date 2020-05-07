

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
library(shinythemes)


### Load in the data and generate objects ###
df <- as_tibble(readRDS('GSE99622_hanamsagar2017_cleaned_melted.rds'))
df2 <- as_tibble(readRDS('gene_ensembl_ids_opentarget.rds'))

genes <- unique(df$gene)
std_err <- function(x) sd(x) / sqrt(length(x))