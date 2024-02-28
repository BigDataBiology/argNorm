#!/usr/bin/env python

import argparse
from pathlib import Path
from Bio import SeqIO
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


    mapping = rgi_hits[['Original ID', "Best_Hit_ARO", 'ARO']]
    mapping = mapping.astype({'ARO': 'str'})
    mapping = mapping.rename(columns={'Best_Hit_ARO': 'Gene Name in CARD'})
    mapping['Database'] = database

    return mapping

if __name__ == '__main__':

    parser = argparse.ArgumentParser("Recover cure ARO mapping for an AMR "
                                     "database using the fasta and an RGI "
                                     "tsv (homolog models only)")
    parser.add_argument("-f", "--fasta", required=False, type=check_file,
						 help="Fasta file containing sequences run through "
							  "RGI")

    parser.add_argument("-r", "--rgi", required=True, type=check_file,
						help="Corresponding rgi output tsv for the fasta file")
    parser.add_argument("-d", "--database", required=True, type=str,
						help="Name of the database",
                        choices=['resfinder', 'ncbi', 'sarg', 'deeparg', 'resfinder_fg', 'megares', 'argannot'])


    args = parser.parse_args()

    mapping = get_aro_for_hits(args.rgi, args.database)

    output_file = f"{args.database}_ARO_mapping.tsv"
    print(f"Writing mapping to {output_file}")
    mapping.to_csv(output_file, sep='\t', index=False)