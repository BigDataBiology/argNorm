import argparse
from pathlib import Path
import pandas as pd
import os
from Bio import SeqIO

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
    database_entries = []
    for record in SeqIO.parse(str(fa), 'fasta'):
        if 'mutation' not in record.id:
            if database != 'resfinder_fg':
                database_entries.append(record.id)
            else:
                database_entries.append(record.description)

    rgi_hits = pd.read_csv(rgi_output, sep='\t')

    if database == 'resfinder':
        rgi_hits['Original ID'] = rgi_hits['Contig'].apply(lambda x: "_".join(x.split('_')[:-1]))
    elif database == 'ncbi':
        rgi_hits['Original ID'] = rgi_hits['ORF_ID']
    elif database == 'sarg':
        rgi_hits['Original ID'] = rgi_hits['ORF_ID'].apply(lambda x: x.split()[0])
    elif database == 'deeparg':
        rgi_hits['Original ID'] = rgi_hits['ORF_ID']
    elif database == 'resfinder_fg':
        rgi_hits['Original ID'] = rgi_hits['ORF_ID']
    elif database == 'argannot':
        rgi_hits['Original ID'] = rgi_hits['Contig'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    elif database == 'megares':
        rgi_hits['Original ID'] = rgi_hits['Contig'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    else:
        raise ValueError(f"Unknown database {database}")


    mapping = rgi_hits[['Original ID', "Best_Hit_ARO", 'ARO']]
    mapping = mapping.astype({'ARO': 'str'})
    mapping = mapping.rename(columns={'Best_Hit_ARO': 'Gene Name in CARD'})

    os.makedirs('manual_curation', exist_ok=True)
    missing_hits_from_original = set(database_entries) - set(mapping['Original ID'].unique())
    missing_hits_from_original = pd.Series(list(missing_hits_from_original))
    missing_hits_from_original.name = "Original ID"
    missing_hits_from_original.to_csv(os.path.join('./manual_curation/', f'./{database}_manual_curation_raw.tsv'), sep='\t', index=False)

    database_entries = pd.Series(list(database_entries))
    database_entries.name = "Original ID"

    mapping['Database'] = database
    mapping = mapping.sort_values(by=['Original ID'])

    return mapping