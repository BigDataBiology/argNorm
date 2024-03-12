from argnorm.drug_categorization import (
    confers_resistance_to,
    drugs_to_drug_classes
)

def test_confers_resistance_to():
    test_cases = ["ARO:3003938", "ARO:3000014"]

    expected_output = [
        [
            ("ARO:3000008", "penam"),
            ("ARO:0000032", "cephalosporin"),
            ("ARO:0000044", "cephamycin"),
        ],
        [
            ("ARO:3000008", "penam"),
            ("ARO:0000004", "monobactam"),
            ("ARO:3003706", "penem"),
            ("ARO:0000032", "cephalosporin"),
            ("ARO:3000637", "ampicillin"),
        ],
    ]

    zipped_test_cases_and_output = zip(test_cases, expected_output)

    for i in zipped_test_cases_and_output:
        assert sorted(confers_resistance_to(i[0])) == sorted(i[1])


def test_get_inaffective_drug_classes():
    test_cases = [
        [("ARO:0000004", "monobactam")],
        [("ARO:3000157", "rifamycin antibiotic")],
        [("ARO:0000030", "tigecycline")],
        [("ARO:0000036", "ciprofloxacin")],
    ]

    expected_output = [
        [("ARO:3000007", "beta-lactam antibiotic")],
        [("ARO:3000157", "rifamycin antibiotic")],
        [("ARO:3000050", "tetracycline antibiotic")],
        [("ARO:0000001", "fluoroquinolone antibiotic")],
    ]

    zipped_test_cases_and_output = zip(test_cases, expected_output)

    for i in zipped_test_cases_and_output:
        assert drugs_to_drug_classes(i[0]) == i[1]
