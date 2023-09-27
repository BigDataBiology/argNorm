from argnorm.drug_categorization import (
    get_immediate_drug_classes,
    get_drug_class_category,
)
import pronto

def test_get_immediate_drug_classes():
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
        assert sorted(get_immediate_drug_classes(i[0])) == sorted(i[1])


def test_get_drug_class_category():
    test_cases = [
        [("ARO:0000004", "monobactam")],
        [("ARO:3000157", "rifamycin antibiotic")],
        [("ARO:0000030", "tigecycline")],
        [("ARO:0000036", "ciprofloxacin")],
    ]

    expected_output = [
        ["beta-lactam antibiotic"],
        ["rifamycin antibiotic"],
        ["tetracycline antibiotic"],
        ["fluoroquinolone antibiotic"],
    ]

    zipped_test_cases_and_output = zip(test_cases, expected_output)

    for i in zipped_test_cases_and_output:
        assert get_drug_class_category(i[0]) == i[1]
