from argnorm.drug_categorization import (
    confers_resistance_to,
    drugs_to_drug_classes
)

def test_confers_resistance_to():
    test_cases = [
        "ARO:3003938",
        "ARO:3000014",
        "ARO:3003725"
    ]

    expected_output = [
        ["ARO:3000008", "ARO:0000032", "ARO:0000044"],
        ["ARO:3000008", "ARO:0000004", "ARO:3003706", "ARO:0000032", "ARO:3000637"],
        ["ARO:3000081"]
    ]

    for t,e in zip(test_cases, expected_output):
        assert sorted(confers_resistance_to(t)) == sorted(e)


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
