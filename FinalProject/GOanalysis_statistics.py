import pandas as pd
from scipy.stats import mannwhitneyu
import numpy as np


def run_mannwhitneyu(exist, nonexist):
    U1, p_val = mannwhitneyu(exist, nonexist, alternative='two-sided')
    return p_val

def cal_foldchange(exist, nonexist):
    exist_mean = np.mean(exist)
    nonexist_mean = np.mean(nonexist)
    if (nonexist_mean == 0) or (exist_mean == 0):
        return np.nan
    return np.log2(exist_mean / nonexist_mean)


def run_statistics_per_GO(GOterm,existDF, nonexistDF):    
    # extract clip and ribosome density value from df -> calculate mannwhitney
    clip_exist = existDF['clip_enrichment'].to_list()
    clip_nonexist = nonexistDF['clip_enrichment'].to_list()
    RNAden_exist = existDF['ribosome_density_change'].to_list()
    RNAden_nonexist = nonexistDF['ribosome_density_change'].to_list()
    
    #Run mannwhitney U
    clip_pval = run_mannwhitneyu(clip_exist, clip_nonexist)
    RNAden_pval = run_mannwhitneyu(RNAden_exist, RNAden_nonexist)
    
    #Run foldchange
    clip_foldchange = cal_foldchange(clip_exist, clip_nonexist)
    RNAden_foldchange = cal_foldchange(RNAden_exist,RNAden_nonexist)
    
    #Gene number
    gene_number = len(existDF)
    
    return [GOterm, gene_number, clip_pval, clip_foldchange, RNAden_pval, RNAden_foldchange]
