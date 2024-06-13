import pandas as pd
from Bio.UniProt.GOA import gafiterator, gpa_iterator
import gzip
import pickle

def map_uniprot_GO(Uniprot_Genes_unique):
    goa_path="/home/imgyuri/BIOINFO1/goa_mouse.gaf"
    GOmapping_result = {unique_id:[] for unique_id in Uniprot_Genes_unique}

    idlist=GOmapping_result.keys()
    with open(goa_path, 'rt') as goa_file:
        for annotation in gafiterator(goa_file):
            id = annotation['DB_Object_ID']
            if id in idlist:
                GOmapping_result[id].append(annotation['GO_ID'])
    
    GOmapping_filtered = {} # remove uniprot id that does not have GO term
    for key, go in GOmapping_result.items():
        if len(go) > 0:
            GOmapping_filtered[key] = go
    GO_mapfiltered_db = pd.DataFrame(GOmapping_filtered.items(), columns=['UniprotID', 'GO'])
    GO_mapfiltered_db.to_parquet('/home/imgyuri/BIOINFO1/Bioinfo1_Labsession/FinalProject/GO_mapping_result.parquet')
    return GO_mapfiltered_db

def map_uniprot_GO_annotation(Uniprot_Genes_unique):
    GOmapping_result = {unique_id:[] for unique_id in Uniprot_Genes_unique}
    idlist=GOmapping_result.keys()
    gpa_mouse_path = "/home/imgyuri/BIOINFO1/gene_association.goa_mouse.gaf"
    with open(gpa_mouse_path, 'rt') as gpa_file:
        for annotation in gafiterator(gpa_file):
            id = annotation['DB_Object_ID']
            if id in idlist:
                GOmapping_result[id].append(annotation['GO_ID'])
    GOmapping_filtered = {} # remove uniprot id that does not have GO term
    for key, go in GOmapping_result.items():
        if len(go) > 0:
            GOmapping_filtered[key] = go
    GO_mapfiltered_db = pd.DataFrame(GOmapping_filtered.items(), columns=['UniprotID', 'GO'])
    GO_mapfiltered_db.to_parquet('/home/imgyuri/BIOINFO1/Bioinfo1_Labsession/FinalProject/GAF_GO_mapping_result.parquet')
    return GO_mapfiltered_db
