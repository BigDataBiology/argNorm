import argparse
from pathlib import Path
import pandas as pd
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
    print(mapping)

    database_entries = pd.Series(list(database_entries))
    database_entries.name = "Original ID"
    print(database_entries)
    mapping = pd.merge(database_entries, mapping, how='outer')
    mapping['Database'] = database

    return mapping

get_aro_for_hits('./dbs/resfinder_fg.faa', './mapping/resfinder_fg_rgi.txt', 'resfinder_fg').to_csv('./mapping/resfinder_fg_ARO_mapping.tsv', sep='\t')