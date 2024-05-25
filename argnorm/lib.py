import os
import pandas as pd
# Pandas <1.5 does not have copy_on_write option
try:
    pd.options.mode.copy_on_write = True
except pd.errors.OptionError:
    pass

ORIGINAL_ID_COL = 'Original ID'
MAPPING_TABLE_ARO_COL = 'ARO'
TARGET_ARO_COL = 'ARO'

_ROOT = os.path.abspath(os.path.dirname(__file__))
_ARO = None

def get_aro_ontology():
    """
    The ARO ontology used by (bundled with) argNorm

    Returns:
        ARO (pronto.Ontology): A pronto ontology object with ARO terms.
    """
    import pronto
    import importlib.resources
    global _ARO
    if _ARO is None:
        _ARO = pronto.Ontology(os.path.join(_ROOT, 'data/aro.obo'))
    return _ARO

def get_aro_mapping_table(database):
    """
    Description: Returns the ARO mapping table for a specific supported databases.

    Parameters:
        database (str): name of database. Can be: argannot, deeparg, megares, ncbi, resfinderfg and sarg

    Returns:
        aro_mapping_table (DataFrame): A pandas dataframe with ARGs mapped to AROs.
    """

    aro_mapping_table = pd.read_csv(
            os.path.join(_ROOT, 'data', f'{database}_ARO_mapping.tsv'),
            sep='\t')
    aro_mapping_table.drop_duplicates(subset=['Original ID'], inplace=True)
    aro_mapping_table.set_index('Original ID', inplace=True)

    manual_curation = pd.read_csv(
                    os.path.join(_ROOT, 'data/manual_curation', f'{database}_curation.tsv'),
                    sep='\t', index_col=0)
    manual_curation['Database'] = aro_mapping_table['Database'].iloc[0]
    aro_mapping_table.drop(index=set(manual_curation.index) & set(aro_mapping_table.index), inplace=True)
    aro_mapping_table = pd.concat([aro_mapping_table, manual_curation])

    aro_mapping_table['ARO'] = aro_mapping_table['ARO'].map(lambda a: f'ARO:{int(a)}', na_action='ignore')
    return aro_mapping_table

def map_to_aro(gene, database):
    """
    Description: Gets ARO mapping for a specific gene in a database.

    Parameters:
        gene (str): The original ID of the gene as mentioned in source database.
        database (str): name of database. Can be: argannot, deeparg, megares, ncbi, resfinderfg and sarg

    Returns:
        ARO[result] (pronto.term.Term): A pronto term with the ARO number of input gene. ARO number can be accessed using 'id' attribute and gene name can be accessed using 'name' attribute.

        If ARO mapping is doesn't exist, None is returned.
    """

    if database not in ['ncbi', 'deeparg', 'resfinder', 'sarg', 'megares', 'argannot']:
        raise Exception(f'{database} is not a supported database.')

    mapping_table = get_aro_mapping_table(database)

    try:
        result = mapping_table.loc[gene, 'ARO']
    except KeyError:
        raise Exception(f'{gene} is not in {database} database')
    else:
        if pd.isna(result):
            return None
        ARO = get_aro_ontology()
        return ARO.get(result)
