# The goal of this code is for providing generalizable, useful, working functions intente
# the creation of developmental indices on any gene expression data (preferably a dataset that
#  has already been adjusted to appropriate values i.e. RPKM/FPKM/TPM)


### (expression value of the gene in a sample—minimum expression for the gene across all samples)/(maximum expression for the gene across all
### samples—minimum expression for the gene across all samples),
### scale all values so that they add equal weight to the index
def scale_expression(df):  
    """Scales the expression of each gene based on minimum expression of gene across all samples present.

    Args: Takes a dataframe as sole argument, where each rows corresponds to a gene, and each column a sample/cell.
    
    """
    import numpy as np 
    import pandas as pd

    scaled_expression_df = df
    
    ## iterating through the rows of the dataframe to scale expression of each gene (by row)
    for row in range(len(scaled_expression_df)):
        scaled_expression_df.iloc[row] = (scaled_expression_df.iloc[row] - np.min(scaled_expression_df.iloc[row])) / np.max(scaled_expression_df.iloc[row] - np.min(scaled_expression_df.iloc[row]))        
    
    return scaled_expression_df



### just a quick function for cleaning data and making sure you don't have and unexpressed genes
def drop_unexpressed_genes(df):

    """
    Just a quick function for cleaning data and making sure you don't have and unexpressed genes

    Args: Takes a dataframe as sole argument, where each rows corresponds to a gene, and each column a sample/cell.

    """
    import numpy as np 
    import pandas as pd

    df.dropna(inplace = True)
    
    return df


### master function that takes in the columns containing young and old data, and calculates genes that are significantly regulated
### by development
def identify_significant_genes(df, young, old):

    """
    Master function that takes in the columns containing young and old data, and calculates genes that are significantly regulated by development

    Args: Takes a dataframe, young samples, and old samples (which are slices of the original dataframe) as the three required arguments

    """
    import numpy as np 
    import pandas as pd
    import scipy.stats as stats

    for gene, series in young.iterrows():
        if np.mean(young.loc[gene]) == 0:
            young.drop(gene, inplace = True)

    for gene, series in old.iterrows():
        if np.mean(old.loc[gene]) == 0:
            old.drop(gene, inplace = True)

    pvals = np.zeros(shape = len(df))
    sig = np.zeros(shape = len(df))
    logdiff = np.zeros(shape = len(df))
    row = -1

    for gene, series in df.iterrows():
        row = row + 1
        if gene in young.index:
            if gene in old.index:
                pvals[row] = stats.ttest_ind(young.loc[gene], old.loc[gene])[1]
                sig[row] = stats.ttest_ind(young.loc[gene], old.loc[gene])[0]
                logdiff[row] = np.log2(np.mean(old.loc[gene])/np.mean(young.loc[gene]))
            else:
                pvals[row] = np.NAN
                sig[row] = np.NAN
                logdiff[row] = np.NAN
        else:
            pvals[row] = np.NAN
            sig[row] = np.NAN
            logdiff[row] = np.NAN


    df['pvals'] = pvals
    df['sig'] = sig
    df['logdiff'] = logdiff

    direction = [0] * len(df)
    for row in range(len(df)):
        if df['pvals'][row] < 0.05:
            if df['logdiff'][row] > 0:
                direction[row] = 'UP'
            else:
                direction[row] = 'DOWN'
        else:
            direction[row] = 'N/A'
    df['direction'] = direction


    return df



### housekeeping function for removing rows (genes) that are not developmentally regulated
def remove_insignificant_rows(df):

    """
    Housekeeping function for removing rows (genes) that are not developmentally regulated

    Args: Takes a dataframe as its sole argument. The dataframe must have 'identify significant genes' computed

    """
    import numpy as np 
    import pandas as pd

    for gene, series in df.iterrows():
        if df['direction'][gene] == 'UP':
            pass
        if df['direction'][gene] == 'DOWN':
            pass
        else:
            df.drop(gene, inplace = True)  
    df.reset_index(inplace=True)
    df.rename(columns = {'index' : 'gene'}, inplace = True)
    return df


### function for extracting the raw info about genes that were developmentally regulated
### i.e. the genes that make up the developmental index
def extract_regulated_genes(df):

    """
    Function for extracting the raw info about genes that were developmentally regulated
    i.e. the genes that make up the developmental index

    Args: Takes a dataframe as its sole argument. The dataframe must have 'identify significant genes' computed

    """

    import numpy as np 
    import pandas as pd
    
    down_genes = df['gene'][df['direction'] == 'DOWN'].to_list()
    up_genes = df['gene'][df['direction'] == 'UP'].to_list()
    regulated_genes = pd.DataFrame((dict([ (k,pd.Series(v)) for k,v in {'gene' : df['gene'], 
                                                                        'direction' : df['direction'], 
                                                                        'valence' : df['logdiff']}.items() ])))
    
    regulated_genes = regulated_genes.sort_values(by = 'valence', ascending = False).set_index('gene')
    
    return regulated_genes


### takes sample columns as an argument and calculates the index value for each sample
def generate_index(df, sample_cols):

    """
    Takes sample columns as an argument and calculates the index value for each sample

    Args: Takes a dataframe and the iloc of columns containing the samples as the two arguments.
    This iloc is usually df.columns[1:-4] 
    The dataframe must have 'identify significant genes' computed

    """
    
    import numpy as np 
    import pandas as pd
    
    ## initialize the dataframe with the sample columns and create an index per sample array
    ## also, create an int counter
    samples = sample_cols
    index_per_sample = [0] * len(samples)
    i = -1

    ## iterate through the sampels and calculate the index given the index genes
    for sample in samples:
        i = i + 1
        if np.mean(df[sample][df['direction'] == 'UP']) > 0:
            if np.mean(df[sample][df['direction'] == 'DOWN']) > 0:
                index_per_sample[i] = np.mean(df[sample][df['direction'] == 'UP']) / np.mean(df[sample][df['direction'] == 'DOWN'])
            else:
                pass
        else:
            pass
    ## create the output dataframe by stitching this index per sample array to the corresponding sample names    
    final_df = pd.DataFrame([index_per_sample], columns = samples)

    return final_df


### function that takes in the index output and normalizes all values to between 0 and 1
### so that the peak development == 1, and 'immature' == 0
def scale_index(df):

    """
    Function that takes in the index output and normalizes all values to between 0 and 1
    So that the peak development == 1, and 'immature' == 0

    Args: Takes a dataframe as its sole argument. This should be done on the dataframe output from 'generate index' function

    """
    
    import numpy as np 
    import pandas as pd

    ## scale data to between 0 and 1
    final_df_scaled = pd.DataFrame(df.iloc[0] - np.min(df.iloc[0])) / (np.max(df.iloc[0] - np.min(df.iloc[0])))
    final_df_scaled.reset_index(inplace = True)
    
    return final_df_scaled


### just a quick merge, make sure that the gene column (join_col) from the two input dataframes are labelled correctly

def import_index(df, gene_list_df):
    
    """
    Just a quick merge of the main dataframe and an imported index gene list. Make sure that both dataframes have original gene IDs as the 
    row indexes.

    Args: Takes main dataframe and dataframe with index genes as the arguments.

    """
    
    import numpy as np 
    import pandas as pd

    output_df = df.join(gene_list_df)

    return output_df