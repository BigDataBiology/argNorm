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
        "ARO:3001414",
        "ARO:3003548",
        "ARO:3000826",
        "ARO:3003066",
        "ARO:3003745"
    ]

    expected_output = [
        ["ARO:3000008", "ARO:0000032"],
        ["ARO:3000008", "ARO:0000004", "ARO:0000032", "ARO:3000637"],
        ['ARO:0000028', 'ARO:0000029'],
        ["ARO:3007382", "ARO:0000052", "ARO:0000007", "ARO:0000035", "ARO:0000049"],
        ["ARO:0000015", "ARO:3000008"],
        ['ARO:0000032', 'ARO:0000056', 'ARO:3000008'],
        ['ARO:0000045', 'ARO:0000047'],
        ['ARO:0000001', 'ARO:0000030', 'ARO:0000051', 'ARO:3000169', 'ARO:3000385', 'ARO:3000637', 'ARO:3000704', 'ARO:3000870'],
        ['ARO:0000016', 'ARO:0000032', 'ARO:3000008'],
        ['ARO:0000000', 'ARO:3000668']
    ]

    for t,e in zip(test_cases, expected_output):
        assert sorted(confers_resistance_to(t)) == sorted(e)

def test_oxa19_drugs():
    oxa_19 = _assert_aro_name('ARO:3001414', 'OXA-19')
    cephalosporin = _assert_aro_name('ARO:0000032', 'cephalosporin')
    penam = _assert_aro_name('ARO:3000008', 'penicillin beta-lactam')
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
        ["ARO:3000081"],
        ["ARO:3007382"],
        ["ARO:3004022"],
        ["ARO:3000199"],
        ["ARO:3003998"],
        ["ARO:3007382", "ARO:0000052", "ARO:0000007", "ARO:0000035", "ARO:0000049"],
        ['ARO:3007145']
    ]

    expected_output = [
        ["ARO:3000007"],
        ["ARO:3000157"],
        ["ARO:3000050"],
        ["ARO:0000001"],
        ["ARO:3000081"],
        ["ARO:0000016"],
        ["ARO:0000026"],
        ["ARO:3000053"],
        ["ARO:3000007"],
        ["ARO:0000016", "ARO:0000016", "ARO:0000016", "ARO:0000016", "ARO:0000016"],
        ["ARO:3000007"]
    ]

    for t, e in zip(test_cases, expected_output):
        assert drugs_to_drug_classes(t) == e

def test_betalactams():
    test_cases = [
        _assert_aro_name('ARO:0000004', 'monobactam'),
        _assert_aro_name('ARO:0000032', 'cephalosporin'),
        _assert_aro_name('ARO:3007145', 'imipenem-cilastatin-relebactam'),
        _assert_aro_name('ARO:3003998', 'ampicillin-sulbactam'),
        _assert_aro_name('ARO:3004705', 'ceftazidime-clavulanic acid'),
        _assert_aro_name('ARO:3007040', 'cefotaxime-ceftiofur-tazobactam-clavulanate')
    ]

    beta_lactam_antibiotic = _assert_aro_name('ARO:3000007', 'beta-lactam antibiotic')
    assert set(drugs_to_drug_classes(test_cases)) == set([beta_lactam_antibiotic])

def test_peptide_antibiotics():
    test_cases = [
        _assert_aro_name('ARO:3000123', 'gramicidin'),
        _assert_aro_name('ARO:3000119', 'edeine'),
        _assert_aro_name('ARO:0000041', 'bacitracin'),
        _assert_aro_name('ARO:3000199', 'gramicidin D')
    ]

    peptide_antibiotic = _assert_aro_name('ARO:3000053', 'peptide antibiotic')
    assert set(drugs_to_drug_classes(test_cases)) == set([peptide_antibiotic])

def test_aminoglycosides():
    test_cases = [
        _assert_aro_name('ARO:3007382', 'gentamicin'),
        _assert_aro_name('ARO:0000052', 'tobramycin'),
        _assert_aro_name('ARO:0000007', 'dibekacin'),
        _assert_aro_name('ARO:0000035', 'sisomicin'),
        _assert_aro_name('ARO:0000049', 'kanamycin A')
    ]

    aminoglycoside_antibiotic = _assert_aro_name('ARO:0000016', 'aminoglycoside antibiotic')
    assert set(drugs_to_drug_classes(test_cases)) == set([aminoglycoside_antibiotic])

def test_streptogramin_mixture():
    quinupristin_dalfopristin = _assert_aro_name('ARO:3004022', 'quinupristin-dalfopristin')

    streptogramin_antibiotic = _assert_aro_name('ARO:0000026', 'streptogramin antibiotic')
    assert drugs_to_drug_classes([quinupristin_dalfopristin]) == [streptogramin_antibiotic]


def test_multiple_parents():
    # ARO:3000130 has multiple parents
    edeine_A = _assert_aro_name('ARO:3000130', 'edeine A')
    assert drugs_to_drug_classes([edeine_A]) == [_assert_aro_name('ARO:3000053', 'peptide antibiotic')]


def test_categories_all():
    import pandas as pd
    golden = pd.read_csv('tests/all_aros_category.tsv.gz',
                         sep='\t',
                         index_col=0)
    golden = golden.T.to_dict(orient='list')

    ARO = get_aro_ontology()

    for ar in ARO.terms():
        ar = ar.id
        drugs = list(confers_resistance_to(ar))
        # the walrus operator is only available in python 3.8+, so we have to
        # do this in two steps
        expected = golden.get(ar)
        if expected:
            cats = drugs_to_drug_classes(drugs)
            drugs = ';'.join(drugs)
            cats = ';'.join(cats)
            assert [drugs, cats] == expected
        else:
            assert not drugs

