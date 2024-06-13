import pandas as pd
import numpy as np
import os
from UniProtMapper import ProtMapper # Get uniprot entry
from mapGOA import map_uniprot_GO, map_uniprot_GO_annotation
from GOanalysis_statistics import run_statistics_per_GO
from statsmodels.stats.multitest import multipletests
import numpy as np 


mapper = ProtMapper()
# from IPython.core.interactiveshell import InteractiveShell 
READ_COUNTS_PATH = '/home/imgyuri/BIOINFO1/binfo1-work/read-counts.txt'
OUPUT_PATH = '/home/imgyuri/BIOINFO1/Bioinfo1_Labsession/FinalProject'

def cnt_preprocess(cnts_ori):
    ## CNT preprocessing
    cnts_ori['Geneid'] = cnts_ori['Geneid'].str.split('.').str[0]
    ## Filter by read counts : RNA-seq (<30 raw reads) or low ribosome footprints (<80 raw footprint tags in siLuc library)
    # cnts = cnts_ori[(cnts_ori['RNA-control.bam'] >= 30)& (cnts_ori['CLIP-35L33G.bam'] >= 30) & (cnts_ori['RNA-siLin28a.bam'] >= 30) & (cnts_ori['RNA-siLuc.bam'] >= 30)&(cnts_ori['RPF-siLuc.bam'] >= 80)].copy()
    ##Trial
    cnts = cnts_ori[(cnts_ori['RNA-control.bam'] >= 30) & (cnts_ori['RNA-siLuc.bam'] >= 30) &(cnts_ori['RPF-siLuc.bam'] >= 80)].copy()

    # cnts['RNA-control.bam'] = (cnts['RNA-control.bam'] / (cnts['Length'] / 1000 * cnts['RNA-control.bam'].sum() / 1e6))
    # cnts['RNA-siLin28a.bam'] = (cnts['RNA-siLin28a.bam'] / (cnts['Length'] / 1000 * cnts['RNA-siLin28a.bam'].sum() / 1e6))

    ## normalize CLIP library
    cnts['clip_enrichment'] = (cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam'])
    cnts['ribosome_density_change'] = (cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam'] / (cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam']))
    cnts.dropna(axis=0, inplace=True)
    print(len(cnts))
    return cnts

def map_ensembl_uniprot(cnts):
    gene_list = cnts['Geneid'].tolist()
    uniprot_scrap_result, failed = mapper.get(
        ids=gene_list, from_db="Ensembl", to_db="UniProtKB"
        )
    print(f"Success : {len(uniprot_scrap_result)}, Failed : {len(failed)}")
    uniprot_scrap_result.dropna(axis=0, inplace=True)
    uniprot_scrap_result.to_parquet(os.path.join(OUPUT_PATH, 'uniprot_mapping_result.parquet'), index=False)
    return uniprot_scrap_result

def extract_uniprot_entryNgenename(uniprot_scrap_result):
    ensembl_uniprot = uniprot_scrap_result.groupby('From').agg({
            'Entry': lambda x: list(set(x)),
            'Gene Names': lambda x: list(set(gene for sublist in x for gene in sublist.split()))
        }).reset_index()
    ensembl_uniprot.columns = ['Geneid', 'UniprotEntry', 'GeneNames']
    ensembl_uniprot.to_parquet(os.path.join(OUPUT_PATH, 'Uniprot_mapping_result_grouped.parquet'))
    # Extract unique genes
    Uniprot_Genes_all = []
    for entry in ensembl_uniprot['UniprotEntry']:
        Uniprot_Genes_all.extend(entry)
    Uniprot_Genes_unique = list(set(Uniprot_Genes_all))
    return ensembl_uniprot, Uniprot_Genes_unique

def get_GO_terms(uniprot_ids, uniprot_GO_dict):
    GO_terms = []
    for id in uniprot_ids:
        # print(id)
        if id in uniprot_GO_dict:
            GO_terms.extend(uniprot_GO_dict[id])  
    if len(GO_terms) == 0:
        return np.nan
    return list(set(GO_terms))

def map_ensembl_GO(uniprot_GO, ensembl_uniprot):
    uniprot_GO_dict = dict(zip(uniprot_GO['UniprotID'], uniprot_GO['GO']))

    ensembl_uniprot['GO'] = ensembl_uniprot['UniprotEntry'].apply(lambda x: get_GO_terms(x, uniprot_GO_dict))
    ensembl_uniprot.dropna(axis=0, inplace=True)
    ensembl_uniprot.to_parquet(os.path.join(OUPUT_PATH, 'Ensembl_GO_mapping_result.parquet'))
    return ensembl_uniprot


def run_stats(Ensembl_GO_CNT_DB):
    Ensembl_GO_CNT_DB = Ensembl_GO_CNT_DB[~Ensembl_GO_CNT_DB['GO'].apply(lambda x: 'GO:0031966' in x)]
    # allGenes = Ensembl_GO_CNT_DB['Geneid'].tolist()
    uniqueGOterms = set(go for sublist in Ensembl_GO_CNT_DB['GO'] for go in sublist)

    rows = []
    for go in uniqueGOterms:
        protein_exist = Ensembl_GO_CNT_DB[Ensembl_GO_CNT_DB['GO'].apply(lambda x: go in x)]
        protein_not_exist = Ensembl_GO_CNT_DB[~Ensembl_GO_CNT_DB['GO'].apply(lambda x: go in x)]
    
        rows.append(run_statistics_per_GO(go, protein_exist, protein_not_exist))
    result_stat = pd.DataFrame(rows, columns=['GOterm', 'GeneNumber', 'ClipPval', 'ClipFoldChange', 'RNAdenPval', 'RNAdenFoldChange'])
    
    ###BH correction
     # Collect all p-values for BH correction
    all_pvals = result_stat[['ClipPval', 'RNAdenPval']].values.flatten()
    
    # Apply BH correction
    _, corrected_pvals, _, _ = multipletests(all_pvals, method='fdr_bh')
    
    # Reshape corrected p-values and assign them back to the DataFrame
    print("Start BH correction")
    corrected_pvals = np.reshape(corrected_pvals, result_stat[['ClipPval', 'RNAdenPval']].shape)
    result_stat[['ClipPval', 'RNAdenPval']] = corrected_pvals
    
    # Filter
    result_stat.dropna(axis=0, inplace=True)
    result_stat = result_stat.sort_values(by='GeneNumber', ascending=False)
    result_stat = result_stat[(result_stat['ClipPval'] <= 0.05) & (result_stat['RNAdenPval'] <= 0.05)] #& (result_stat['GeneNumber']>=10)]
    result_stat.to_parquet(os.path.join(OUPUT_PATH, 'Statistics_result.parquet'))
    return result_stat

def main():
    cnts_ori = pd.read_csv(READ_COUNTS_PATH, sep='\t', comment='#')
    cnts = cnt_preprocess(cnts_ori)
    print("Start ENSEMBL to UNIPROT mapping")
    # ensemblNuniprot = map_ensembl_uniprot(cnts)
    # Load ensemble to uniprot mapping result
    ensemblNuniprot = pd.read_parquet(os.path.join(OUPUT_PATH, 'uniprot_mapping_result.parquet'), engine='pyarrow')
    ensembl_uniprot_filtered, unique_uniprot_entries = extract_uniprot_entryNgenename(ensemblNuniprot)
    print("Start UNIPROT to GO mapping")
    # uniprot_GO_DB = map_uniprot_GO(unique_uniprot_entries)
    uniprot_GO_DB = map_uniprot_GO_annotation(unique_uniprot_entries)
    ensembl_GO_DB = map_ensembl_GO(uniprot_GO_DB, ensembl_uniprot_filtered)
    Ensembl_GO_CNT_DB = pd.merge(cnts, ensembl_GO_DB, left_on='Geneid', right_on='Geneid', how='inner')
    ## Run Statistics
    print("Start Statistics")
    Ensembl_GO_CNT_DB.to_parquet(os.path.join(OUPUT_PATH, 'Ensembl_GO_CNT_DB_main.parquet'))
    run_stats(Ensembl_GO_CNT_DB)
    print("Done. Now move to notebook for plotting.")

if __name__ == '__main__':
    main()
    
    ## Remove GO:0031966