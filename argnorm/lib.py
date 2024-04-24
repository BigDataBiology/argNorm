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

    manual_curation = pd.read_csv(os.path.join(_ROOT, 'data/manual_curation', f'{database}_curation.tsv'), sep='\t')
    manual_curation['Database'] = df['Database']
    aro_mapping_table = df

    if database != 'resfinder':
        aro_mapping_table = pd.concat([df, manual_curation])
    else:
        # Handle gene clusters and reverse complements
        aro_mapping_table = aro_mapping_table.drop_duplicates(subset=['Original ID'], ignore_index=True).set_index('Original ID')

        for i in manual_curation['Original ID']:
            if i in aro_mapping_table.index:
                aro_mapping_table.loc[i, 'ARO'] = manual_curation.set_index('Original ID').loc[i, 'ARO']
                aro_mapping_table.loc[i, 'Gene Name in CARD'] = manual_curation.set_index('Original ID').loc[i, 'Gene Name in CARD']
            else:
                aro_mapping_table.loc[i] = manual_curation.set_index('Original ID').loc[i]

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