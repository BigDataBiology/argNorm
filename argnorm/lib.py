import os
# do not import pandas here because it is slow to import

# import here so that __version__ is available in the package
from .argnorm_version import __version__

ORIGINAL_ID_COL = 'Original ID'
MAPPING_TABLE_ARO_COL = 'ARO'
TARGET_ARO_COL = 'ARO'
CUT_OFF_COL = "Cut_Off"

DATABASES = [
    'argannot',
    'deeparg',
    'megares',
    'ncbi',
    'resfinder',
    'resfinderfg',
    'sarg',
    'groot-db',
    'groot-core-db',
    'groot-argannot',
    'groot-resfinder',
    'groot-card',
]

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
    import pandas as pd

    if 'groot' in database:
        database = 'groot'

    aro_mapping_table = pd.read_csv(
            os.path.join(_ROOT, 'data', f'{database}_ARO_mapping.tsv'),
            sep='\t', dtype={'ARO': str})
    aro_mapping_table.drop_duplicates(subset=['Original ID'], inplace=True)
    aro_mapping_table.set_index('Original ID', inplace=True)

    manual_curation = pd.read_csv(
                    os.path.join(_ROOT, 'data/manual_curation', f'{database}_curation.tsv'),
                    sep='\t', index_col=0, dtype={'ARO': str})
    if database != 'groot':
        manual_curation['Database'] = aro_mapping_table['Database'].iloc[0]
    manual_curation[CUT_OFF_COL] = 'Manual'
    aro_mapping_table.drop(index=set(manual_curation.index) & set(aro_mapping_table.index), inplace=True)
    aro_mapping_table = pd.concat([aro_mapping_table, manual_curation])

    aro_mapping_table['ARO'] = aro_mapping_table['ARO'].map(lambda a: f'ARO:{a}', na_action='ignore')
    return aro_mapping_table

def map_to_aro(gene, database):
    """
    Description: Gets ARO mapping for a specific gene in a database.

    Parameters:
        gene (str): The original ID of the gene as mentioned in source database.
        database (str): name of database. Can be: argannot, deeparg, megares, ncbi, resfinderfg, sarg, groot-db, groot-core-db, groot-argannot, groot-resfinder, groot-card

    Returns:
        ARO[result] (pronto.term.Term): A pronto term with the ARO number of input gene. ARO number can be accessed using 'id' attribute and gene name can be accessed using 'name' attribute.

        If ARO mapping is doesn't exist, None is returned.
    """
    import pandas as pd

    if database not in DATABASES:
        raise Exception(f'{database} is not a supported database.')

    mapping_table = get_aro_mapping_table(database)

    # Preprocess input gene & mapping table original ids if groot is being used
    if database == 'groot-argannot':
        gene = gene.split('~~~')[-1]
        mapping_table.index = mapping_table.index.map(lambda x: ':'.join(str(x).split(':')[1:3]))
    if database == 'groot-card':
        gene = gene.split('.')[0]
    if database in ['groot-db', 'groot-core-db']:
        if 'card' in gene.lower():
            gene = gene.split('|')[-1]
        else:
            gene = gene.split('__')[1]

    try:
        result = mapping_table.loc[gene, 'ARO']
    except KeyError:
        raise Exception(f'{gene} is not in {database} database')
    else:
        if pd.isna(result):
            return None
        ARO = get_aro_ontology()
        return ARO.get(result)
