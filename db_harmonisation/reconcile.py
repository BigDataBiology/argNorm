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


def get_aro_for_hits(fasta, rgi_output, database):
    """
    Extract the mappings and any missed mappings for other AMR databases
    vs NCBI
    """
    database_entries = []
    for record in SeqIO.parse(str(fasta), 'fasta'):
        if "mutation"  not in record.id:
            database_entries.append(record.id)
    database_entries = set(database_entries)

    rgi_hits = pd.read_csv(rgi_output, sep='\t')

    if database == 'resfinder':
        rgi_hits['Original ID'] = rgi_hits['Contig'].apply(lambda x: "_".join(x.split('_')[:-1]))
    elif database == 'ncbi':
        rgi_hits['Original ID'] = rgi_hits['ORF_ID']

    # homolog models only for now
    rgi_hits = rgi_hits[rgi_hits['Model_type'] == "protein homolog model"]


    # tidy up "ORF ID"
    mapping = rgi_hits[['Original ID', "Best_Hit_ARO", 'ARO']]
    mapping = mapping.astype({'ARO': 'str'})
    mapping = mapping.rename(columns={'Best_Hit_ARO': 'Gene Name in CARD'})

    # fill in any missing mappings with NAs
    missing_hits_from_original = database_entries - set(mapping['Original ID'].unique())
    print(f"Warning {len(missing_hits_from_original)} in {database} without "
           "trivial mapping to ARO. Listed as 'nan' in output mapping tsv")

    # merge with full set of database entries to add appropriate NA values
    database_entries = pd.Series(list(database_entries))
    database_entries.name = "Original ID"
    mapping = pd.merge(database_entries, mapping, how='outer')
    mapping['Database'] = database

    return mapping

if __name__ == '__main__':

    parser = argparse.ArgumentParser("Recover cure ARO mapping for an AMR "
                                     "database using the fasta and an RGI "
                                     "tsv (homolog models only)")
    parser.add_argument("-f", "--fasta", required=True, type=check_file,
						 help="Fasta file containing sequences run through "
							  "RGI")

    parser.add_argument("-r", "--rgi", required=True, type=check_file,
						help="Corresponding rgi output tsv for the fasta file")
    parser.add_argument("-d", "--database", required=True, type=str,
						help="Name of the database", choices=['resfinder', 'ncbi'])


    args = parser.parse_args()

    mapping = get_aro_for_hits(args.fasta, args.rgi, args.database)

    output_file = f"{args.database}_ARO_mapping.tsv"
    print(f"Writing mapping to {output_file}")
    mapping.to_csv(output_file, sep='\t')
