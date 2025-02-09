import pytest
from argnorm import lib
from argnorm.lib import map_to_aro, get_aro_mapping_table
import pandas as pd

ARO = lib.get_aro_ontology()

def test_map_to_aro():
    test_cases = [
        ["(AGly)AAC(6')-Isa:NG_047311:101-574:474", 'argannot'],
        ["MEG_21|Drugs|Aminoglycosides|Aminoglycoside_N-acetyltransferases|AAC3", 'megares'],
        ["1028085756|WP_063844287.1|1|1|cpt|cpt|phosphotransferase|2|CHLORAMPHENICOL|PHENICOL|chloramphenicol_phosphotransferase_CPT", 'ncbi'],
        ["gb|AAG57600.1|ARO:3000318|mphB", "sarg"],
        ["(Phe)cpt_strepv:U09991:AAB36569:1412-1948:537", "argannot"],
        ["MEG_4060|Metals|Multi-metal_resistance|Multi-metal_resistance_protein|MREA", "megares"],
        ["gi:447201629:ref:WP_001278885.1:|FEATURES|cob(I)alamin_adenolsyltransferase|unclassified|cob(I)alamin_adenolsyltransferase", "deeparg"],
        ["argannot~~~(Bla)cfxA4~~~AY769933:1-966", 'groot-argannot'],
        ["ErmF.3000498.M17124.1181-1982.593", 'groot-card'],
        ["groot-db_RESFINDER__tet(W)_1_DQ060146", 'groot-db']
    ]

    expected_output = [
        ARO.get_term('ARO:3002563'),
        ARO.get_term('ARO:3004623'),
        ARO.get_term('ARO:3000249'),
        ARO.get_term('ARO:3000318'),
        ARO.get_term('ARO:3000249'),
        None,
        ARO.get_term('ARO:0010004'),
        ARO.get_term('ARO:3003005'),
        ARO.get_term('ARO:3000498'),
        ARO.get_term('ARO:3000194')
    ]

    for t, e in zip(test_cases, expected_output):
        if t[1] == 'groot':
            assert map_to_aro(t[0], t[1], t[2]) == e
        else:
            assert map_to_aro(t[0], t[1]) == e

@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi', 'resfinder', 'resfinderfg', 'groot', 'groot-argannot'])
def test_get_aro_mapping_table_smoke(database):
    df = get_aro_mapping_table(database)
    assert len(df) > 0

def _search_argnorm_mappings(mappings, db):
    mappings.drop(columns=['UpdatedHeader', 'Database'], inplace=True)
    aros = []
    mapping_table = get_aro_mapping_table(db)
    mapping_table.index = list(mapping_table.index.map(lambda x: str(x).lower()))

    for header in mappings['Source_header']:
        try:
            aros.append(str(mapping_table.loc[header.strip().lower(), 'ARO'])[4:])
        except KeyError:
            aros.append(None)

    mappings['ARO'] = aros
    mappings['Database'] = db
    mappings.rename(columns={'MEGARes_header': 'Original ID'}, inplace=True)
    return mappings.dropna()

def _get_resfinder_mappings(megares_headers):
    mappings = _search_argnorm_mappings(megares_headers[
        ((megares_headers.Database == 'ResFinder') | (megares_headers.Database == 'Resfinder'))\
        & (~ megares_headers.MEGARes_header.str.contains('Multi-compound'))\
        & (~ megares_headers.MEGARes_header.str.contains('Metals'))\
        & (~ megares_headers.MEGARes_header.str.contains('Biocides'))], 'resfinder')
    return mappings

def _get_argannot_mappings(megares_headers):
    mappings = _search_argnorm_mappings(megares_headers[megares_headers.Database == 'ARG-ANNOT'], 'argannot')
    return mappings

def test_megares_mappings():
    megares_headers = pd.read_csv('./tests/megares_mappings/megares_headers.csv')   
    argannot_and_resfinder_mappings = pd.concat([_get_argannot_mappings(megares_headers), _get_resfinder_mappings(megares_headers)])
    
    for i in range(argannot_and_resfinder_mappings.shape[0]):
        gene = argannot_and_resfinder_mappings.iloc[i]['Original ID']
        aro = argannot_and_resfinder_mappings.iloc[i]['ARO']
        assert map_to_aro(gene, 'megares') == ARO.get_term(f'ARO:{aro}')
