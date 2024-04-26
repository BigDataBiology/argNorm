import os
import pandas as pd
import pronto

ORIGINAL_ID_COL = 'Original ID'
MAPPING_TABLE_ARO_COL = 'ARO'
TARGET_ARO_COL = 'ARO'

_ROOT = os.path.abspath(os.path.dirname(__file__))

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
    """

    if database not in ['ncbi', 'deeparg', 'resfinder', 'sarg', 'megares', 'argannot']:
        raise Exception(f'{database} is not a supported database.')

    mapping_table = get_aro_mapping_table(database)

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
