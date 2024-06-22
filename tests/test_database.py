import pytest
from argnorm import lib
from argnorm import drug_categorization

ARO = lib.get_aro_ontology()

@pytest.mark.parametrize('database', lib.DATABASES)
def test_database(database):
    determinant_of_ab_resistance = ARO['ARO:3000000']
    antibiotic_molecule = ARO['ARO:1000003']
    assert determinant_of_ab_resistance.name == 'determinant of antibiotic resistance'
    assert antibiotic_molecule.name == 'antibiotic molecule'

    db = lib.get_aro_mapping_table(database)
    assert len(db.index) == len(set(db.index)), f'Duplicate gene names in {database}'
    for ar in db['ARO'].dropna():
        assert ar in ARO, f'ARO not found in the ontology: {ar}'
        assert determinant_of_ab_resistance in ARO[ar].superclasses(), f'ARO not a determinant of antibiotic resistance: {ar}'
        for d in drug_categorization.confers_resistance_to(ar):
            assert d in ARO, f'ARO not found in the ontology: {d}'
            assert antibiotic_molecule in ARO[d].superclasses(), f'ARO not an antibiotic molecule: {d}'
            for dc in drug_categorization.drugs_to_drug_classes([d]):
                assert dc in ARO, f'ARO not found in the ontology: {dc}'
                assert antibiotic_molecule in  ARO[dc].superclasses(1), f'ARO term not an immediate child of antibiotic molecule: {dc}'

def test_misannotated_genes():
    assert lib.map_to_aro('(Tet)tetH:EF460464:6286-7839:1554', 'argannot') == ARO['ARO:3000175']