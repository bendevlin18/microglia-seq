

## Global.R file to load in all data upfront and perform preprocessing ##

### Load in the data and generate objects ###
df <- as_tibble(read.csv('GSE99622_hanamsagar2017_cleaned_melted.csv'))
df2 <- as_tibble(read.csv('gene_ensembl_ids_opentarget.csv'))

genes <- unique(df$gene)
std_err <- function(x) sd(x) / sqrt(length(x))