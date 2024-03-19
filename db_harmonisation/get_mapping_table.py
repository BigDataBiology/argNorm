import argparse
from pathlib import Path
import pandas as pd

def check_file(path):
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")
    
def get_aro_for_hits(rgi_output, database):
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
    mapping['Database'] = database

    return mapping