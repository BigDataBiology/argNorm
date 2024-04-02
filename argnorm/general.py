import os
import pandas as pd

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

def get_data_path(path, getting_manual_curation):
    if getting_manual_curation:
        return os.path.join(_ROOT, 'data/manual_curation', path)

    return os.path.join(_ROOT, 'data', path)

def get_aro_mapping_table(database):
    df = pd.read_csv(get_data_path(f'{database}_ARO_mapping.tsv', False), sep='\t')

    manual_curation = pd.read_csv(get_data_path(f'{database}_curation.tsv', True), sep='\t')
    manual_curation['Database'] = df['Database']

    aro_mapping_table = pd.concat([df, manual_curation])
    aro_mapping_table[TARGET_ARO_COL] = aro_mapping_table[TARGET_ARO_COL].map(lambda a: f'ARO:{int(a)}' if is_number(a) else a)
    
    return aro_mapping_table

def map_to_aro(gene, database):
    mapping_table = get_aro_mapping_table(database).set_index('Original ID')
    result = mapping_table.loc[gene, 'ARO']

    # Dealing with duplicated genes in ARO mapping table.
    # Getting only one ARO number
    if type(result) != str:
        return list(set(result))[0]
    else:
        return result