import argparse
from pathlib import Path
import pandas as pd
import os
from Bio import SeqIO
from check_mapping_accuracy import check_mapping_accuracy, preprocess_mappings_for_tests 

def check_file(path):
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")

def get_aro_for_hits(fa, rgi_output, database):
    """
    Generates ARO mapping tables by copying Best_Hit_ARO, ARO and ORF_ID/Contig columns from RGI output
    """
    database_entries = []
    for record in SeqIO.parse(str(fa), 'fasta'):
        if 'mutation' not in record.id:
            if database != 'resfinder_fg':
                database_entries.append(record.id)
            else:
                database_entries.append(record.description)

    rgi_hits = pd.read_csv(rgi_output, sep='\t')

    if database == 'sarg':
        rgi_hits['Original ID'] = rgi_hits['ORF_ID'].apply(lambda x: x.split()[0])
    else:
        rgi_hits['Original ID'] = rgi_hits['ORF_ID']

    mapping = rgi_hits[['Original ID', "Best_Hit_ARO", 'ARO', 'Cut_Off']]
    mapping = mapping.astype({'ARO': 'str'})
    mapping = mapping.rename(columns={'Best_Hit_ARO': 'Gene Name in CARD'})
    
    checked_mappings = check_mapping_accuracy(preprocess_mappings_for_tests(mapping.copy(), database))
    metal_biocide_and_virulence_genes = checked_mappings['metal_biocide_virulence_genes']
    mismatched_genes = checked_mappings['mismatched_genes']

    os.makedirs('manual_curation', exist_ok=True)
    manual_curation = set(database_entries) - set(mapping['Original ID'].unique())
    manual_curation = list(manual_curation) + metal_biocide_and_virulence_genes
    manual_curation = pd.Series(manual_curation)
    manual_curation.name = "Original ID"
    manual_curation.to_csv(os.path.join('./manual_curation/', f'./{database}_manual_curation_raw.tsv'), sep='\t', index=False)
    mismatched_genes.to_csv(os.path.join('./manual_curation/', f'./{database}_mismatched_genes.tsv'), sep='\t', index=False)

    database_entries = pd.Series(list(database_entries))
    database_entries.name = "Original ID"

    mapping['Database'] = database
    mapping = mapping.sort_values(by=['Original ID'])

    return mapping