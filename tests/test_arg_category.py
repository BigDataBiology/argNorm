import argnorm.drug_categorization as drug_categorization

def test_get_immediate_drug_classes():
    test_cases = ['ARO:3003938', 'ARO:3000014']

    expected_output = [
        [['ARO:3000008', 'penam'], ['ARO:0000032', 'cephalosporin'], ['ARO:0000044', 'cephamycin']],
        [['ARO:3000008', 'penam'], ['ARO:0000004', 'monobactam'], ['ARO:3003706', 'penem'], ['ARO:0000032', 'cephalosporin'], ['ARO:3000637', 'ampicillin']]
    ]

    for i in range(len(test_cases)):
        assert sorted(drug_categorization.get_immediate_drug_classes(test_cases[i])) == sorted(expected_output[i]), drug_categorization.get_immediate_drug_classes(test_cases[i])

def test_get_drug_class_category():
    test_cases = [
        [['ARO:0000004', 'monobactam']],
        [['ARO:3000157', 'rifamycin antibiotic']],
        [['ARO:0000030', 'tigecycline']],
        [['ARO:0000036', 'ciprofloxacin']]
    ]

    expected_output = [['beta-lactam antibiotic'], ['rifamycin antibiotic'], ['tetracycline antibiotic'], ['fluoroquinolone antibiotic']]

    for i in range(len(test_cases)):
        assert drug_categorization.get_drug_class_category(test_cases[i]) == expected_output[i], i
