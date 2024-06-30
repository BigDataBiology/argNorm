import pytest
from argnorm import lib
from argnorm.lib import map_to_aro, get_aro_mapping_table

def test_map_to_aro():
    test_cases = [
        ["(AGly)AAC(6')-Isa:NG_047311:101-574:474", 'argannot'],
        ["MEG_21|Drugs|Aminoglycosides|Aminoglycoside_N-acetyltransferases|AAC3", 'megares'],
        ["1028085756|WP_063844287.1|1|1|cpt|cpt|phosphotransferase|2|CHLORAMPHENICOL|PHENICOL|chloramphenicol_phosphotransferase_CPT", 'ncbi'],
        ["gb|AAG57600.1|ARO:3000318|mphB", "sarg"],
        ["(Phe)cpt_strepv:U09991:AAB36569:1412-1948:537", "argannot"],
        ["MEG_4060|Metals|Multi-metal_resistance|Multi-metal_resistance_protein|MREA", "megares"],
        ["gi:447201629:ref:WP_001278885.1:|FEATURES|cob(I)alamin_adenolsyltransferase|unclassified|cob(I)alamin_adenolsyltransferase", "deeparg"],
        ['cat(pC194)_1_NC_002013', 'resfinder']
    ]

    ARO = lib.get_aro_ontology()
    expected_output = [
        ARO.get_term('ARO:3002563'),
        ARO.get_term('ARO:3004623'),
        ARO.get_term('ARO:3000249'),
        ARO.get_term('ARO:3000318'),
        ARO.get_term('ARO:3000249'),
        None,
        ARO.get_term('ARO:0010004'),
        ARO.get_term('ARO:3004455')
    ]

    for t, e in zip(test_cases, expected_output):
        assert map_to_aro(t[0], t[1]) == e

@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi', 'resfinder', 'resfinderfg'])
def test_get_aro_mapping_table_smoke(database):
    df = get_aro_mapping_table(database)
    assert len(df) > 0

