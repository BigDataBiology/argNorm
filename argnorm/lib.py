import os
import pandas as pd
import pronto

ORIGINAL_ID_COL = 'Original ID'
MAPPING_TABLE_ARO_COL = 'ARO'
TARGET_ARO_COL = 'ARO'

_ROOT = os.path.abspath(os.path.dirname(__file__))

def is_number(num):
    try:
        int(num)
    except ValueError:
        return False

    return True

def get_aro_mapping_table(database):
    df = pd.read_csv(os.path.join(_ROOT, 'data', f'{database}_ARO_mapping.tsv'), sep='\t')

    missing_from_rgi = pd.read_csv(os.path.join(_ROOT, 'data/manual_curation', f'{database}_curation.tsv'), sep='\t')
    missing_from_rgi['Database'] = df['Database']

    aro_mapping_table = pd.concat([df, missing_from_rgi])

    # Handle gene clusters and reverse complements
    if database == 'resfinder':
        aro_mapping_table = aro_mapping_table.drop_duplicates(subset=['Original ID']).set_index('Original ID')
        cluster_rc_correction = pd.read_csv(os.path.join(_ROOT, 'data/cluster_rc_correction', f'{database}_cluster_rc_correction.tsv'), sep='\t')

        for i in cluster_rc_correction['Original ID']:
            aro_mapping_table.loc[i, 'ARO'] = cluster_rc_correction.set_index('Original ID').loc[i, 'ARO']
            aro_mapping_table.loc[i, 'Gene Name in CARD'] = cluster_rc_correction.set_index('Original ID').loc[i, 'Gene Name in CARD']

    aro_mapping_table[TARGET_ARO_COL] = aro_mapping_table[TARGET_ARO_COL].map(lambda a: f'ARO:{int(a)}' if is_number(a) else a)
    return aro_mapping_table.reset_index()

def map_to_aro(gene, database):
    if database not in ['ncbi', 'deeparg', 'resfinder', 'sarg', 'megares', 'argannot']:
        raise Exception(f'{database} is not a supported database.')

    mapping_table = get_aro_mapping_table(database).set_index('Original ID')

    try:
        result = mapping_table.loc[gene, 'ARO']
    except KeyError:
        raise Exception(f'{gene} is not in {database} database')
    else:
        # Dealing with duplicated genes in ARO mapping table.
        # Getting only one ARO number
        ARO = pronto.Ontology.from_obo_library('aro.obo')
        if type(result) != str:
            return ARO[list(set(result))[0]]
        else:
            return ARO[result]