from argnorm.drug_categorization import (
    confers_resistance_to,
    drugs_to_drug_classes
)
from argnorm.lib import get_aro_ontology

def _assert_aro_name(aro_id, name):
    ARO = get_aro_ontology()
    assert ARO[aro_id].name == name
    return aro_id

def test_confers_resistance_to():
    test_cases = [
        "ARO:3003938",
        "ARO:3000014",
        "ARO:3003725",
        "ARO:3000230",
        "ARO:3005046",
        "ARO:3001414"
    ]

    expected_output = [
        ["ARO:3000008", "ARO:0000032", "ARO:0000044"],
        ["ARO:3000008", "ARO:0000004", "ARO:3003706", "ARO:0000032", "ARO:3000637"],
        ["ARO:3000081"],
        ["ARO:3007382", "ARO:0000052", "ARO:0000007", "ARO:0000035", "ARO:0000049"],
        ["ARO:0000015", "ARO:3000008"],
        ['ARO:0000032', 'ARO:0000056', 'ARO:3000008']
    ]

    for t,e in zip(test_cases, expected_output):
        assert sorted(confers_resistance_to(t)) == sorted(e)

def test_oxa19_drugs():
    oxa_19 = _assert_aro_name('ARO:3001414', 'OXA-19')
    cephalosporin = _assert_aro_name('ARO:0000032', 'cephalosporin')
    penam = _assert_aro_name('ARO:3000008', 'penam')
    oxacillin = _assert_aro_name('ARO:0000056', 'oxacillin')
    assert sorted(confers_resistance_to(oxa_19)) == sorted([cephalosporin, penam, oxacillin])

def test_aadb_drugs():
    # aadb is a synonym for ANT(2'')-Ia
    aadb = _assert_aro_name('ARO:3000230', "ANT(2'')-Ia")
    gentamicin = _assert_aro_name('ARO:3007382', 'gentamicin')
    tobramycin = _assert_aro_name('ARO:0000052', 'tobramycin')
    dibekacin = _assert_aro_name('ARO:0000007', 'dibekacin')
    sisomicin = _assert_aro_name('ARO:0000035', 'sisomicin')
    kanamycin_a = _assert_aro_name('ARO:0000049', 'kanamycin A')
    assert sorted(confers_resistance_to(aadb)) == sorted([
        gentamicin,
        tobramycin,
        dibekacin,
        sisomicin,
        kanamycin_a
    ])

def test_drug_to_drug_classes():
    test_cases = [
        ["ARO:0000004"],
        ["ARO:3000157"],
        ["ARO:0000030"],
        ["ARO:0000036"],
        ["ARO:3000081"]
    ]
    expected_output = [
        ["ARO:3000007"],
        ["ARO:3000157"],
        ["ARO:3000050"],
        ["ARO:0000001"],
        ["ARO:3000081"]
    ]
    for t, e in zip(test_cases, expected_output):
        assert drugs_to_drug_classes(t) == e
